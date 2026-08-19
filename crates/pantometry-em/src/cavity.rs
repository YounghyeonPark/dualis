//! A rectangular box with conducting walls, marched by the Yee update.

use glam::DVec3;
use pantometry_core::conserved::quantity;
use pantometry_core::{Domain, Exchange, Kind, Ledger, Reading, Violation};
use pantometry_units::{Energy, Frequency, Length, LengthVec, Time};

use crate::{cavity_frequency, Medium};

/// The Courant number a three-dimensional Yee grid is stable up to.
///
/// `dt ≤ dx / (c√3)`, so the number is `1/√3`. The same shape as the acoustic wave equation's, and
/// for the same reason: three axes each contribute a term to the amplification factor and the
/// worst case is the mode that is sharpest along all three at once.
pub const COURANT_3D: f64 = 0.577_350_269_189_625_8;

/// One of the six faces of the box.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Wall {
    /// `x = 0`.
    XLow,
    /// `x = L`.
    XHigh,
    /// `y = 0`.
    YLow,
    /// `y = L`.
    YHigh,
    /// `z = 0`.
    ZLow,
    /// `z = L`.
    ZHigh,
}

impl Wall {
    /// All six, in index order.
    pub const ALL: [Wall; 6] = [
        Wall::XLow,
        Wall::XHigh,
        Wall::YLow,
        Wall::YHigh,
        Wall::ZLow,
        Wall::ZHigh,
    ];

    fn index(self) -> usize {
        match self {
            Wall::XLow => 0,
            Wall::XHigh => 1,
            Wall::YLow => 2,
            Wall::YHigh => 3,
            Wall::ZLow => 4,
            Wall::ZHigh => 5,
        }
    }
}

/// What a face does to a wave that reaches it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Boundary {
    /// A perfect conductor: the tangential electric field is zero, and everything comes back.
    Conducting,
    /// A perfect *magnetic* conductor: the tangential magnetic field is zero, so the tangential
    /// electric field has no normal gradient and is mirrored across the face.
    ///
    /// # What it is for
    ///
    /// Not a material. It is the boundary a **symmetry plane** has, and it is the only way to put
    /// a plane wave in a box without periodicity: a wave travelling along `z` with `E ∥ ŷ` has `Ey`
    /// tangential to the `x` faces, so a conductor there would force it to zero and turn the plane
    /// wave into a waveguide mode with a cutoff and the wrong phase velocity.
    ///
    /// A conductor on the `y` faces and a magnetic conductor on the `x` ones leaves that wave
    /// exactly uniform across the box, which is what makes a one-dimensional problem statable in a
    /// three-dimensional grid.
    Magnetic,
    /// Mur's first-order absorbing condition: the wave leaves, mostly.
    ///
    /// # How open "open" is, and why the number is not zero
    ///
    /// Mur's condition is the discrete form of the one-way wave equation `∂E/∂t = −c ∂E/∂n`, which
    /// a wave arriving **along the normal** satisfies exactly. Such a wave leaves with no
    /// reflection at all in the continuum, and with `O((kΔ)²)` of one on a grid.
    ///
    /// A wave arriving at an angle does not satisfy it. The reflection coefficient goes as
    /// `(1−cos θ)/(1+cos θ)`, so grazing incidence reflects almost everything — which is why a
    /// point source in a box, radiating at every angle at once, leaves a few percent of its energy
    /// behind rather than none. That figure is measured rather than claimed, and it is the honest
    /// alternative to a perfectly matched layer, which would be right to four digits and is a
    /// larger piece of machinery than this crate has earned.
    Open,
}

/// A rectangular cavity with perfectly conducting walls.
///
/// # The grid
///
/// `counts` is cells. `E` lives on the cell edges and `H` on the cell faces, half a cell apart —
/// see the crate docs for the layout and for why it is the only one on which the divergence of the
/// curl vanishes identically.
///
/// # The walls
///
/// Perfect conductor: the **tangential** electric field is zero on every wall. That is a real
/// boundary — a copper box, near enough — and not an approximation to open space. A pulse released
/// here reflects and comes back, on purpose.
#[derive(Clone, Debug)]
pub struct Cavity {
    name: String,
    counts: (usize, usize, usize),
    dx: f64,
    medium: Medium,
    /// `Ex` on edges along x: `nx · (ny+1) · (nz+1)`.
    ex: Vec<f64>,
    /// `Ey`: `(nx+1) · ny · (nz+1)`.
    ey: Vec<f64>,
    /// `Ez`: `(nx+1) · (ny+1) · nz`.
    ez: Vec<f64>,
    /// `Hx` on faces normal to x: `(nx+1) · ny · nz`.
    hx: Vec<f64>,
    /// `Hy`: `nx · (ny+1) · nz`.
    hy: Vec<f64>,
    /// `Hz`: `nx · ny · (nz+1)`.
    hz: Vec<f64>,
    /// What each face does, in [`Wall::ALL`] order.
    boundary: [Boundary; 6],
    /// What fills each cell. Uniform unless [`Cavity::fill`] says otherwise.
    cells: Vec<Medium>,
    /// Cells that are perfect conductors: every electric edge touching one is held at zero.
    solid: Vec<bool>,
    /// Whether any cell is, so the common case costs nothing.
    any_solid: bool,
    /// Permittivity, conductivity and permeability **at each field component's own place**, which
    /// is not a cell centre — see [`Cavity::refresh_media`].
    eps: [Vec<f64>; 3],
    sigma: [Vec<f64>; 3],
    mu: [Vec<f64>; 3],
    /// The tangential electric field one step ago, for the faces that need it.
    ///
    /// Mur's condition relates the boundary at `n+1` to the boundary and its neighbour at `n`, so
    /// one step of history is the whole of what an absorbing face costs.
    previous: Option<Box<Tangential>>,
    /// Whether `H` has been advanced the half step the leapfrog opens with.
    started: bool,
    /// What leapfrog actually conserves, evaluated at the one instant it exists.
    ///
    /// **`½εE² + ½μH²` is not an invariant of this scheme**, and the crate's own closed-form test
    /// says so with the arithmetic: `E` and `H` are half a step apart, so the naive sum oscillates
    /// about the true one by `2 sin(ωΔt/2)` peak to peak — 10.7% for a 1.77 GHz mode at the Courant
    /// limit. A ledger reporting the naive sum makes the audit stop every correct cavity run, and
    /// loosening the scene's budget to let it through would be turning the audit off for a reason
    /// that has a closed form.
    ///
    /// `½ε|Eⁿ|² + ½μ Hⁿ⁻¹ᐟ²·Hⁿ⁺¹ᐟ²`, and **the timing is the whole difficulty**. After a step the
    /// state is `Eⁿ⁺¹` and `Hⁿ⁺¹ᐟ²`, which do not pair: the invariant at `n+1` would need `Hⁿ⁺³ᐟ²`,
    /// which has not happened. The one moment all three terms exist is *between* the magnetic and
    /// electric updates, so it is computed there and kept. The first version paired the post-step
    /// electric energy with them and made the drift **worse** — 9.1% against the naive sum's 4.5%
    /// — which is what a half-step error looks like when it is mistaken for a fix.
    ///
    /// `None` before the first step, where the naive sum is the honest answer.
    invariant: Option<f64>,
    dissipated: f64,
    saved: Option<Box<Saved>>,
}

/// The three electric components one step ago, which is all Mur's condition needs of the past.
#[derive(Clone, Debug)]
struct Tangential {
    ex: Vec<f64>,
    ey: Vec<f64>,
    ez: Vec<f64>,
}

#[derive(Clone, Debug)]
struct Saved {
    ex: Vec<f64>,
    ey: Vec<f64>,
    ez: Vec<f64>,
    hx: Vec<f64>,
    hy: Vec<f64>,
    hz: Vec<f64>,
    started: bool,
    dissipated: f64,
}

impl pantometry_core::ScalarField for Cavity {
    /// Volts per metre — a **magnitude**, not a vector. See [`Cavity::as_field`].
    fn unit(&self) -> &'static str {
        "V/m"
    }

    /// `|E|` in the cell `p` falls in.
    ///
    /// Nearest cell rather than interpolated: the components already live on three different edges
    /// and `electric_at` averages each pair to centre them, so interpolating that average again
    /// would be a blur of a blur — the same reason `Channel` samples its speed this way.
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let (nx, ny, nz) = self.counts;
        let q = p.to_si() / self.dx;
        if q.is_nan() {
            return f64::NAN;
        }
        let pick = |v: f64, n: usize| (v.floor().max(0.0) as usize).min(n.saturating_sub(1));
        self.electric_at(pick(q.x, nx), pick(q.y, ny), pick(q.z, nz))
            .length()
    }
}

impl Cavity {
    /// An empty box of `counts` cubic cells of side `cell`.
    pub fn new(
        name: impl Into<String>,
        counts: (usize, usize, usize),
        cell: Length,
        medium: Medium,
    ) -> Cavity {
        let counts = (counts.0.max(1), counts.1.max(1), counts.2.max(1));
        let (nx, ny, nz) = counts;
        let mut built = Cavity {
            name: name.into(),
            counts,
            dx: cell.to_si(),
            medium,
            ex: vec![0.0; nx * (ny + 1) * (nz + 1)],
            ey: vec![0.0; (nx + 1) * ny * (nz + 1)],
            ez: vec![0.0; (nx + 1) * (ny + 1) * nz],
            hx: vec![0.0; (nx + 1) * ny * nz],
            hy: vec![0.0; nx * (ny + 1) * nz],
            hz: vec![0.0; nx * ny * (nz + 1)],
            boundary: [Boundary::Conducting; 6],
            cells: vec![medium; nx * ny * nz],
            solid: vec![false; nx * ny * nz],
            any_solid: false,
            eps: [Vec::new(), Vec::new(), Vec::new()],
            sigma: [Vec::new(), Vec::new(), Vec::new()],
            mu: [Vec::new(), Vec::new(), Vec::new()],
            previous: None,
            started: false,
            invariant: None,
            dissipated: 0.0,
            saved: None,
        };
        built.refresh_media();
        built
    }

    /// Fill the cells a predicate selects with a different medium.
    ///
    /// How a slab, a coating or a lens goes into the box. The predicate is handed cell indices.
    pub fn fill(
        &mut self,
        medium: Medium,
        which: impl Fn(usize, usize, usize) -> bool,
    ) -> &mut Cavity {
        let (nx, ny, nz) = self.counts;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if which(i, j, k) {
                        self.cells[i + nx * (j + ny * k)] = medium;
                    }
                }
            }
        }
        self.refresh_media();
        self
    }

    /// Make the cells a predicate selects perfect conductors.
    ///
    /// # A staircase, and what that costs
    ///
    /// Every electric edge touching a solid cell is held at zero after each update, which is the
    /// boundary condition a perfect conductor has. It is exact on a surface that lies along the
    /// grid and a **staircase** on one that does not, so a flat screen normal to an axis — a slit, a
    /// mask, a ground plane — is represented exactly and a curved one is not.
    ///
    /// The alternative is a very large conductivity through [`Cavity::fill`], and it is worse: the
    /// semi-implicit update's decay factor is `(1−σΔt/2ε)/(1+σΔt/2ε)`, which tends to **−1** rather
    /// than 0, so a good conductor made that way flips sign every step instead of holding zero.
    pub fn obstruct(&mut self, which: impl Fn(usize, usize, usize) -> bool) -> &mut Cavity {
        let (nx, ny, nz) = self.counts;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if which(i, j, k) {
                        self.solid[i + nx * (j + ny * k)] = true;
                        self.any_solid = true;
                    }
                }
            }
        }
        self.enforce_solid();
        self
    }

    /// Whether one cell is a perfect conductor.
    pub fn is_solid(&self, i: usize, j: usize, k: usize) -> bool {
        let (nx, ny, nz) = self.counts;
        self.solid[i.min(nx - 1) + nx * (j.min(ny - 1) + ny * k.min(nz - 1))]
    }

    /// Zero every electric edge that touches a solid cell.
    fn enforce_solid(&mut self) {
        if !self.any_solid {
            return;
        }
        let (nx, ny, nz) = self.counts;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if !self.solid[i + nx * (j + ny * k)] {
                        continue;
                    }
                    // The twelve edges of the cell: four along each axis.
                    for (dj, dk) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                        let at = self.iex(i, j + dj, k + dk);
                        self.ex[at] = 0.0;
                    }
                    for (di, dk) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                        let at = self.iey(i + di, j, k + dk);
                        self.ey[at] = 0.0;
                    }
                    for (di, dj) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                        let at = self.iez(i + di, j + dj, k);
                        self.ez[at] = 0.0;
                    }
                }
            }
        }
    }

    /// What fills one cell.
    pub fn medium_at(&self, i: usize, j: usize, k: usize) -> Medium {
        let (nx, ny, nz) = self.counts;
        self.cells[i.min(nx - 1) + nx * (j.min(ny - 1) + ny * k.min(nz - 1))]
    }

    /// Recompute each field component's material from the cells around it.
    ///
    /// # A field component does not sit in a cell
    ///
    /// `Ey` lives on an edge shared by four cells and `Hx` on a face shared by two, so neither has
    /// *a* permittivity — it has the average of the ones it borders. That averaging is what puts a
    /// sharp interface exactly on a grid plane and keeps it second order; taking the material of
    /// the nearest cell instead moves the interface half a cell and is first order, which shows up
    /// as a reflectance that is right to one digit and converges at the wrong rate.
    fn refresh_media(&mut self) {
        let (nx, ny, nz) = self.counts;
        let at = |i: isize, j: isize, k: isize| -> usize {
            let c = |v: isize, n: usize| v.clamp(0, n as isize - 1) as usize;
            c(i, nx) + nx * (c(j, ny) + ny * c(k, nz))
        };
        let mean = |cells: &[Medium], of: &[usize], f: fn(&Medium) -> f64| -> f64 {
            of.iter().map(|c| f(&cells[*c])).sum::<f64>() / of.len() as f64
        };
        let permittivity = |m: &Medium| m.permittivity();
        let conductivity = |m: &Medium| m.conductivity.max(0.0);
        let permeability = |m: &Medium| m.permeability();

        let mut eps = [
            vec![0.0; self.ex.len()],
            vec![0.0; self.ey.len()],
            vec![0.0; self.ez.len()],
        ];
        let mut sigma = eps.clone();
        for k in 0..=nz {
            for j in 0..=ny {
                for i in 0..=nx {
                    if i < nx {
                        let of = [
                            at(i as isize, j as isize - 1, k as isize - 1),
                            at(i as isize, j as isize, k as isize - 1),
                            at(i as isize, j as isize - 1, k as isize),
                            at(i as isize, j as isize, k as isize),
                        ];
                        let idx = self.iex(i, j, k);
                        eps[0][idx] = mean(&self.cells, &of, permittivity);
                        sigma[0][idx] = mean(&self.cells, &of, conductivity);
                    }
                    if j < ny {
                        let of = [
                            at(i as isize - 1, j as isize, k as isize - 1),
                            at(i as isize, j as isize, k as isize - 1),
                            at(i as isize - 1, j as isize, k as isize),
                            at(i as isize, j as isize, k as isize),
                        ];
                        let idx = self.iey(i, j, k);
                        eps[1][idx] = mean(&self.cells, &of, permittivity);
                        sigma[1][idx] = mean(&self.cells, &of, conductivity);
                    }
                    if k < nz {
                        let of = [
                            at(i as isize - 1, j as isize - 1, k as isize),
                            at(i as isize, j as isize - 1, k as isize),
                            at(i as isize - 1, j as isize, k as isize),
                            at(i as isize, j as isize, k as isize),
                        ];
                        let idx = self.iez(i, j, k);
                        eps[2][idx] = mean(&self.cells, &of, permittivity);
                        sigma[2][idx] = mean(&self.cells, &of, conductivity);
                    }
                }
            }
        }

        let mut mu = [
            vec![0.0; self.hx.len()],
            vec![0.0; self.hy.len()],
            vec![0.0; self.hz.len()],
        ];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = self.ihx(i, j, k);
                    mu[0][idx] = mean(
                        &self.cells,
                        &[
                            at(i as isize - 1, j as isize, k as isize),
                            at(i as isize, j as isize, k as isize),
                        ],
                        permeability,
                    );
                    let idx = self.ihy(i, j, k);
                    mu[1][idx] = mean(
                        &self.cells,
                        &[
                            at(i as isize, j as isize - 1, k as isize),
                            at(i as isize, j as isize, k as isize),
                        ],
                        permeability,
                    );
                    let idx = self.ihz(i, j, k);
                    mu[2][idx] = mean(
                        &self.cells,
                        &[
                            at(i as isize, j as isize, k as isize - 1),
                            at(i as isize, j as isize, k as isize),
                        ],
                        permeability,
                    );
                }
            }
        }
        // The `H` faces the loop above does not reach lie outside the box in one index. Nothing
        // reads them, but a zero permeability would divide by zero if anything ever did.
        let fallback = self.medium.permeability();
        for component in mu.iter_mut() {
            for v in component.iter_mut() {
                if *v <= 0.0 {
                    *v = fallback;
                }
            }
        }
        self.eps = eps;
        self.sigma = sigma;
        self.mu = mu;
    }

    /// Launch a plane wave travelling along `+z`, polarised along `y`.
    ///
    /// The envelope is given as a function of `z` in metres. `H` is set from it exactly:
    ///
    /// ```text
    ///   E = f(z − ct) ŷ,      H = (1/η) ẑ × E = −f(z − ct)/η · x̂
    /// ```
    ///
    /// evaluated half a cell along and half a step back, which is where the staggering puts it.
    /// Getting either offset wrong leaves a backward-travelling remainder that a reflectance
    /// measurement counts as a reflection.
    ///
    /// The envelope is a function of **position**, not only of `z`, so a guided mode's transverse
    /// profile can be given with it: `|x, _, z| (PI * x / a).sin() * pulse(z)` launches TE₁₀ down a
    /// waveguide, and `|_, _, z| pulse(z)` launches a plane wave.
    ///
    /// A plane wave wants [`Boundary::Magnetic`] on the `x` faces — see that variant for why. A
    /// guided mode does not: its profile already vanishes on the conducting walls, which is what
    /// makes it a mode.
    ///
    /// # `H` is the plane wave's, even for a mode
    ///
    /// A mode's wave impedance is `η/√(1−(f_c/f)²)` and it carries a longitudinal `H` besides, so
    /// this launches a superposition rather than a pure mode. That is fine and is *why* the
    /// measurements downstream are of the field far from the source: everything that is not a
    /// propagating mode is either below its own cutoff and gone within a few wavelengths, or
    /// orthogonal to the profile and never excited. What survives is the mode, whatever was
    /// launched.
    pub fn launch_along_z(
        &mut self,
        dt: Time,
        envelope: impl Fn(f64, f64, f64) -> f64,
    ) -> &mut Cavity {
        let (nx, ny, nz) = self.counts;
        let c = self.medium.wave_speed().to_si();
        let eta = (self.medium.permeability() / self.medium.permittivity()).sqrt();
        for k in 0..=nz {
            for j in 0..ny {
                for i in 0..=nx {
                    let idx = self.iey(i, j, k);
                    self.ey[idx] = envelope(
                        i as f64 * self.dx,
                        (j as f64 + 0.5) * self.dx,
                        k as f64 * self.dx,
                    );
                }
            }
        }
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..=nx {
                    let z = (k as f64 + 0.5) * self.dx + c * dt.to_si() / 2.0;
                    let idx = self.ihx(i, j, k);
                    self.hx[idx] =
                        -envelope(i as f64 * self.dx, (j as f64 + 0.5) * self.dx, z) / eta;
                }
            }
        }
        self.enforce_walls_where(|b| b == Boundary::Conducting);
        self.started = true;
        self
    }

    /// Add to `Ey` on one plane of constant `z`, as a driven source.
    ///
    /// A **soft** source: it adds rather than assigns, so a wave arriving at the plane passes
    /// through instead of being erased. It radiates in both directions, which an absorbing face
    /// behind it takes care of.
    ///
    /// # What an initial condition cannot do
    ///
    /// [`Cavity::launch_along_z`] sets a field and lets it go, which is everything for a
    /// propagating wave and nothing for an **evanescent** one: a field below a waveguide's cutoff
    /// is a near field, it does not travel, and with no source it simply decays to zero. There is no
    /// steady state for it to settle into and therefore no spatial profile to measure.
    ///
    /// So a driven source is not a convenience here; it is the only way to ask the question.
    pub fn impress(&mut self, plane: usize, add: impl Fn(f64, f64) -> f64) -> &mut Cavity {
        let (nx, ny, nz) = self.counts;
        let k = plane.min(nz);
        for j in 0..ny {
            for i in 0..=nx {
                let idx = self.iey(i, j, k);
                self.ey[idx] += add(i as f64 * self.dx, (j as f64 + 0.5) * self.dx);
            }
        }
        self.enforce_walls_where(|b| b == Boundary::Conducting);
        self
    }

    /// The cutoff frequency of a rectangular guide's `(m, n)` mode.
    ///
    /// ```text
    ///   f_c = (c/2) √((m/a)² + (n/b)²)
    /// ```
    ///
    /// The cavity resonance with the propagating index removed, which is not a coincidence: a
    /// cavity is a guide with two more walls, and a mode of one is a mode of the other standing
    /// still.
    pub fn cutoff_frequency(&self, mode: (u32, u32)) -> Frequency {
        let s = self.size().to_si();
        let sum = (mode.0 as f64 / s.x).powi(2) + (mode.1 as f64 / s.y).powi(2);
        Frequency::from_si(0.5 * self.medium.wave_speed().to_si() * sum.sqrt())
    }

    /// Make a face absorb rather than reflect, or the other way round.
    pub fn set_boundary(&mut self, wall: Wall, boundary: Boundary) -> &mut Cavity {
        self.boundary[wall.index()] = boundary;
        self
    }

    /// Make every face absorb. What a piece of open space is, as near as Mur gets.
    pub fn open(&mut self) -> &mut Cavity {
        for wall in Wall::ALL {
            self.set_boundary(wall, Boundary::Open);
        }
        self
    }

    /// What a face currently does.
    pub fn boundary(&self, wall: Wall) -> Boundary {
        self.boundary[wall.index()]
    }

    /// Mur's coefficient, `(cΔt − Δ)/(cΔt + Δ)`.
    ///
    /// A function of the Courant number and nothing else. At the 3D stability limit it is −0.268;
    /// at half of it, −0.552.
    pub fn mur_coefficient(&self, dt: Time) -> f64 {
        let c_dt = self.medium.wave_speed().to_si() * dt.to_si();
        (c_dt - self.dx) / (c_dt + self.dx)
    }

    /// Release a compact, **divergence-free** pulse to watch leave.
    ///
    /// A blob of `Hz`, uniform along `z` and Gaussian in `x` and `y`. That combination has
    /// `∇·H = 0` exactly on the discrete grid — the only non-zero term would be `∂Hz/∂z`, and there
    /// is none — so it is a real radiating field and not a charge.
    ///
    /// # Why it cannot be a blob of `E`
    ///
    /// That was the first version, and it left **34%** of its energy in an open box where the
    /// radiating part leaves in one crossing. A localised `Ey` has `∂Ey/∂y ≠ 0`, so `∇·D ≠ 0`, so
    /// it *is* a charge distribution — and the near field of a charge does not propagate. No
    /// absorbing boundary can remove what was never going anywhere, and a measurement of one that
    /// used such a source would be measuring the electrostatics instead.
    ///
    /// Uniform along `z` makes it a line source, radiating cylindrically. A genuinely compact
    /// three-dimensional source needs a current loop, which is more machinery than a boundary test
    /// needs.
    pub fn pulse(&mut self, at: (usize, usize), amplitude: f64, width_cells: f64) -> &mut Cavity {
        let (nx, ny, nz) = self.counts;
        let w = width_cells.max(0.5);
        for k in 0..=nz {
            for j in 0..ny {
                for i in 0..nx {
                    let r2 = ((i as f64 + 0.5 - at.0 as f64).powi(2)
                        + (j as f64 + 0.5 - at.1 as f64).powi(2))
                        / (w * w);
                    let idx = self.ihz(i, j, k);
                    self.hz[idx] += amplitude * (-r2).exp();
                }
            }
        }
        self.started = true;
        self
    }

    /// Cells along each axis.
    pub fn counts(&self) -> (usize, usize, usize) {
        self.counts
    }

    /// The cell side.
    pub fn cell(&self) -> Length {
        Length::from_si(self.dx)
    }

    /// The box's outside dimensions.
    pub fn size(&self) -> LengthVec {
        LengthVec::from_si(
            DVec3::new(
                self.counts.0 as f64,
                self.counts.1 as f64,
                self.counts.2 as f64,
            ) * self.dx,
        )
    }

    /// What fills it.
    pub fn medium(&self) -> Medium {
        self.medium
    }

    /// The closed-form frequency of one of the box's modes.
    pub fn mode_frequency(&self, mode: [u32; 3]) -> Frequency {
        let s = self.size().to_si();
        cavity_frequency([s.x, s.y, s.z], mode, &self.medium)
    }

    /// Fill the cavity with the `(m, 0, p)` standing mode, at an amplitude.
    ///
    /// # Why this family and not a general `(m, n, p)`
    ///
    /// `E = ŷ E₀ sin(mπx/a) sin(pπz/d)` is a mode of the box and is **exactly** representable on
    /// the grid: it is tangential to the four walls it has to vanish on and normal to the two it
    /// does not, so the boundary condition is satisfied by the function rather than by clamping it
    /// afterwards. A general `(m, n, p)` needs the full TE/TM construction, and releasing an
    /// approximation to a mode would make every frequency measurement a measurement of the
    /// approximation.
    ///
    /// `H` is set a **half step behind**, which is where the leapfrog wants it. Getting that wrong
    /// costs first-order accuracy and shows up as a frequency that converges at the wrong rate
    /// rather than as anything visibly broken — the defect this workspace has already found once,
    /// in `Room`.
    pub fn release_mode(&mut self, mode: (u32, u32), amplitude: f64, dt: Time) -> &mut Cavity {
        let (m, p) = (mode.0.max(1) as f64, mode.1.max(1) as f64);
        let (nx, ny, nz) = self.counts;
        let s = self.size().to_si();
        let (kx, kz) = (
            std::f64::consts::PI * m / s.x,
            std::f64::consts::PI * p / s.z,
        );
        self.ex.fill(0.0);
        self.ey.fill(0.0);
        self.ez.fill(0.0);
        self.hx.fill(0.0);
        self.hy.fill(0.0);
        self.hz.fill(0.0);

        // Ey at (i, j+½, k): sin(kx·x) sin(kz·z), at t = 0.
        for k in 0..=nz {
            for j in 0..ny {
                for i in 0..=nx {
                    let (x, z) = (i as f64 * self.dx, k as f64 * self.dx);
                    let at = self.iey(i, j, k);
                    self.ey[at] = amplitude * (kx * x).sin() * (kz * z).sin();
                }
            }
        }

        self.enforce_walls();

        // **`H` is the *discrete* curl of `E`, not a sampling of the continuum one**, and that
        // distinction is the whole of why `∇·B` is zero rather than small.
        //
        // `∇·(∇×F) = 0` is an identity of the *discrete* operators on this grid — every term
        // appears twice and cancels — but sampling a continuum field that happens to be
        // divergence-free does not inherit it. The first version set `H` from the analytic curl
        // and measured `3.3e-4` where the identity says `0`.
        //
        // Applying `advance_magnetic(−dt/2)` to a zero field gives `H = +(dt/2μ)∇×E` with the
        // scheme's own stencil, which is both divergence-free by construction and exactly the
        // half step the leapfrog wants: for a standing mode `E ∝ cos ωt` and `H ∝ sin ωt`, so
        // `H(+dt/2) = −H(−dt/2)` and the single update between them fixes both.
        self.advance_magnetic(-0.5 * dt.to_si());
        self.enforce_walls();
        self.started = true;
        // **And the invariant at `n = 0`, which needs a field that has not happened yet.**
        // `½ε|E⁰|² + ½μ Hⁿ⁻¹ᐟ²·Hⁿ⁺¹ᐟ²` pairs the `H` just set with the one the first step will
        // produce, so it cannot be read off the state — it is produced here by taking that step
        // and putting it back. The alternative is an opening ledger holding the *naive* energy
        // while every later one holds the invariant, and the audit then sees the difference
        // between the two as a loss: measured 5.369e-3 on scene 27, which is half the swing and
        // has nothing to do with the run.
        let keep = (self.hx.clone(), self.hy.clone(), self.hz.clone());
        self.advance_magnetic(dt.to_si());
        self.invariant = Some(self.electric_energy().to_si() + self.paired_magnetic_energy(&keep));
        self.hx = keep.0;
        self.hy = keep.1;
        self.hz = keep.2;
        self
    }

    /// The electric field at a point on the grid, by nearest edge.
    pub fn electric_at(&self, i: usize, j: usize, k: usize) -> DVec3 {
        let (nx, ny, nz) = self.counts;
        DVec3::new(
            self.ex[self.iex(i.min(nx - 1), j.min(ny), k.min(nz))],
            self.ey[self.iey(i.min(nx), j.min(ny - 1), k.min(nz))],
            self.ez[self.iez(i.min(nx), j.min(ny), k.min(nz - 1))],
        )
    }

    /// The magnetic field strength at a point on the grid, by nearest face.
    pub fn magnetic_at(&self, i: usize, j: usize, k: usize) -> DVec3 {
        let (nx, ny, nz) = self.counts;
        DVec3::new(
            self.hx[self.ihx(i.min(nx), j.min(ny - 1), k.min(nz - 1))],
            self.hy[self.ihy(i.min(nx - 1), j.min(ny), k.min(nz - 1))],
            self.hz[self.ihz(i.min(nx - 1), j.min(ny - 1), k.min(nz))],
        )
    }

    /// The trapezoidal weight of a sample: half at each end of an axis it is *node*-centred on.
    ///
    /// # The half cells at the wall are half cells
    ///
    /// A Yee component is cell-centred along the axes it points across and node-centred along the
    /// one it points along. A node-centred sample at an end of the box represents **half** a
    /// cell's worth of field, and counting it whole over-states the energy by the boundary layer's
    /// share — measured, 12.8% of a swing where the leapfrog's own `O((ωdt)²)` predicts 0.6%.
    ///
    /// It shows up as an energy that appears to oscillate wildly rather than as anything visibly
    /// wrong, because the electric field of a cavity mode vanishes at the wall and the magnetic
    /// one does not: only one of the two was over-counted, so the pair no longer balanced.
    fn weight(index: usize, last: usize, node_centred: bool) -> f64 {
        if node_centred && (index == 0 || index == last) {
            0.5
        } else {
            1.0
        }
    }

    /// Sum `material · f²` over one component array, with the boundary half-cells weighted as
    /// such.
    fn weighted_square_with(
        &self,
        field: &[f64],
        material: &[f64],
        extent: (usize, usize, usize),
        node: (bool, bool, bool),
    ) -> f64 {
        let (ex, ey, ez) = extent;
        let mut total = 0.0;
        for k in 0..ez {
            for j in 0..ey {
                for i in 0..ex {
                    let w = Cavity::weight(i, ex - 1, node.0)
                        * Cavity::weight(j, ey - 1, node.1)
                        * Cavity::weight(k, ez - 1, node.2);
                    let idx = i + ex * (j + ey * k);
                    total += w * material[idx] * field[idx] * field[idx];
                }
            }
        }
        total
    }

    /// Sum `f²` over one component array, with the boundary half-cells weighted as such.
    #[allow(dead_code)]
    fn weighted_square(
        &self,
        field: &[f64],
        extent: (usize, usize, usize),
        node: (bool, bool, bool),
    ) -> f64 {
        let (ex, ey, ez) = extent;
        let mut total = 0.0;
        for k in 0..ez {
            for j in 0..ey {
                for i in 0..ex {
                    let w = Cavity::weight(i, ex - 1, node.0)
                        * Cavity::weight(j, ey - 1, node.1)
                        * Cavity::weight(k, ez - 1, node.2);
                    let v = field[i + ex * (j + ey * k)];
                    total += w * v * v;
                }
            }
        }
        total
    }

    /// Electric energy, `½∫εE² dV`.
    pub fn electric_energy(&self) -> Energy {
        let (nx, ny, nz) = self.counts;
        let v = self.dx.powi(3);
        // Ex is cell-centred in x and node-centred in y and z; the others rotate that. The
        // permittivity is per component, so a box with a slab in it reports the energy the slab
        // actually holds rather than the box's nominal one.
        let sum = self.weighted_square_with(
            &self.ex,
            &self.eps[0],
            (nx, ny + 1, nz + 1),
            (false, true, true),
        ) + self.weighted_square_with(
            &self.ey,
            &self.eps[1],
            (nx + 1, ny, nz + 1),
            (true, false, true),
        ) + self.weighted_square_with(
            &self.ez,
            &self.eps[2],
            (nx + 1, ny + 1, nz),
            (true, true, false),
        );
        Energy::from_si(0.5 * sum * v)
    }

    /// `½μ H_before·H_now`, cell by cell — the magnetic half of the leapfrog invariant.
    ///
    /// The same weighting `magnetic_energy` uses, with the square replaced by the product of the
    /// two half-step-separated fields. Written beside it rather than in terms of it because a
    /// product is not a square and sharing the code would have meant a flag saying which.
    fn paired_magnetic_energy(&self, before: &(Vec<f64>, Vec<f64>, Vec<f64>)) -> f64 {
        let (nx, ny, nz) = self.counts;
        let v = self.dx.powi(3);
        let pair = |now: &[f64],
                    was: &[f64],
                    mu: &[f64],
                    counts: (usize, usize, usize),
                    edge: (bool, bool, bool)| {
            let (cx, cy, cz) = counts;
            let mut sum = 0.0;
            for k in 0..cz {
                for j in 0..cy {
                    for i in 0..cx {
                        let at = i + cx * (j + cy * k);
                        // Half weight on a face that is shared, exactly as the square does.
                        let mut w = 1.0;
                        if edge.0 && (i == 0 || i + 1 == cx) {
                            w *= 0.5;
                        }
                        if edge.1 && (j == 0 || j + 1 == cy) {
                            w *= 0.5;
                        }
                        if edge.2 && (k == 0 || k + 1 == cz) {
                            w *= 0.5;
                        }
                        sum += w * mu[at] * now[at] * was[at];
                    }
                }
            }
            sum
        };
        0.5 * v
            * (pair(
                &self.hx,
                &before.0,
                &self.mu[0],
                (nx + 1, ny, nz),
                (true, false, false),
            ) + pair(
                &self.hy,
                &before.1,
                &self.mu[1],
                (nx, ny + 1, nz),
                (false, true, false),
            ) + pair(
                &self.hz,
                &before.2,
                &self.mu[2],
                (nx, ny, nz + 1),
                (false, false, true),
            ))
    }

    /// What this scheme actually conserves: `½εE(t)² + ½μ H(t−Δt/2)·H(t+Δt/2)`.
    ///
    /// [`energy`](Cavity::energy) is the physical field energy and is what a reader wants; this is
    /// what the *audit* has to watch, because the difference between them is a real oscillation
    /// with a closed form and not a drift.
    pub fn invariant_energy(&self) -> Energy {
        match self.invariant {
            Some(u) => Energy::from_si(u),
            None => self.energy(),
        }
    }

    /// Magnetic energy, `½∫μH² dV`.
    pub fn magnetic_energy(&self) -> Energy {
        let (nx, ny, nz) = self.counts;
        let v = self.dx.powi(3);
        let sum = self.weighted_square_with(
            &self.hx,
            &self.mu[0],
            (nx + 1, ny, nz),
            (true, false, false),
        ) + self.weighted_square_with(
            &self.hy,
            &self.mu[1],
            (nx, ny + 1, nz),
            (false, true, false),
        ) + self.weighted_square_with(
            &self.hz,
            &self.mu[2],
            (nx, ny, nz + 1),
            (false, false, true),
        );
        Energy::from_si(0.5 * sum * v)
    }

    /// The field energy in the box.
    pub fn energy(&self) -> Energy {
        Energy::from_si(self.electric_energy().to_si() + self.magnetic_energy().to_si())
    }

    /// What a lossy medium has turned into heat over the run.
    pub fn dissipated(&self) -> Energy {
        Energy::from_si(self.dissipated)
    }

    /// The largest `|∇·B|` anywhere, relative to what a single cell's field would give.
    ///
    /// **Zero, and not to a tolerance.** Every term in the discrete divergence of the discrete curl
    /// appears twice with opposite signs, so the update cannot change this quantity at all —
    /// whatever it was when the field was set, it stays. A number that drifts here is a scheme that
    /// is not Yee's, however much it looks like it.
    ///
    /// Reported relative to `|H|·dx²`, the scale a divergence of one cell's worth would have, so it
    /// is a number and not a unit.
    pub fn magnetic_divergence(&self) -> f64 {
        let scale = self
            .hx
            .iter()
            .chain(&self.hy)
            .chain(&self.hz)
            .fold(0.0f64, |a, b| a.max(b.abs()))
            .max(f64::MIN_POSITIVE);
        self.peak_magnetic_divergence() / scale
    }

    /// The largest `|∇·B|` anywhere, **unnormalised**, in A/m.
    ///
    /// The absolute figure, because the relative one divides by `max |H|` and that oscillates
    /// through the cycle: comparing a normalised divergence at two instants compares two different
    /// denominators, and reports a factor of two where the numerator did not move at all. Measured
    /// — the first version of the injection test read `1.099` and then `0.5387` for a quantity the
    /// update cannot touch.
    pub fn peak_magnetic_divergence(&self) -> f64 {
        let (nx, ny, nz) = self.counts;
        let mut worst: f64 = 0.0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let d = (self.hx[self.ihx(i + 1, j, k)] - self.hx[self.ihx(i, j, k)])
                        + (self.hy[self.ihy(i, j + 1, k)] - self.hy[self.ihy(i, j, k)])
                        + (self.hz[self.ihz(i, j, k + 1)] - self.hz[self.ihz(i, j, k)]);
                    worst = worst.max(d.abs());
                }
            }
        }
        worst
    }

    /// Add a magnetic monopole's worth of divergence at one cell, for a test to watch.
    ///
    /// There is no physical reason to want this. It exists so that the claim "the update cannot
    /// change `∇·B`" can be checked from the other side: put some in, and it neither grows nor
    /// heals. A scheme without the identity does one or the other.
    pub fn inject_divergence(&mut self, i: usize, j: usize, k: usize, amount: f64) {
        let (nx, ny, nz) = self.counts;
        let at = self.ihx(i.min(nx), j.min(ny - 1), k.min(nz - 1));
        self.hx[at] += amount;
    }

    /// The largest `|∇·D|` at an interior node, relative to `|E|`.
    ///
    /// Preserved by the same identity as the magnetic one, and zero for the same reason, but only
    /// where the node has all six neighbours: at a wall the constraint involves the surface charge
    /// the conductor carries, which this does not track.
    pub fn electric_divergence(&self) -> f64 {
        let (nx, ny, nz) = self.counts;
        let mut worst: f64 = 0.0;
        let scale = self
            .ex
            .iter()
            .chain(&self.ey)
            .chain(&self.ez)
            .fold(0.0f64, |a, b| a.max(b.abs()))
            .max(f64::MIN_POSITIVE);
        for k in 1..nz {
            for j in 1..ny {
                for i in 1..nx {
                    let d = (self.ex[self.iex(i, j, k)] - self.ex[self.iex(i - 1, j, k)])
                        + (self.ey[self.iey(i, j, k)] - self.ey[self.iey(i, j - 1, k)])
                        + (self.ez[self.iez(i, j, k)] - self.ez[self.iez(i, j, k - 1)]);
                    worst = worst.max(d.abs());
                }
            }
        }
        worst / scale
    }

    /// The largest step this grid is stable at, `dx/(c√3)`.
    pub fn courant_limit(&self) -> Time {
        Time::from_si(COURANT_3D * self.dx / self.medium.wave_speed().to_si())
    }

    // --- indexing ----------------------------------------------------------
    //
    // Written out rather than generated, because every one of the six has a different extent and
    // a transposition here is a scheme that still runs.

    fn iex(&self, i: usize, j: usize, k: usize) -> usize {
        let (nx, ny, _) = self.counts;
        i + nx * (j + (ny + 1) * k)
    }
    fn iey(&self, i: usize, j: usize, k: usize) -> usize {
        let (nx, ny, _) = self.counts;
        i + (nx + 1) * (j + ny * k)
    }
    fn iez(&self, i: usize, j: usize, k: usize) -> usize {
        let (nx, ny, _) = self.counts;
        i + (nx + 1) * (j + (ny + 1) * k)
    }
    fn ihx(&self, i: usize, j: usize, k: usize) -> usize {
        let (nx, ny, _) = self.counts;
        i + (nx + 1) * (j + ny * k)
    }
    fn ihy(&self, i: usize, j: usize, k: usize) -> usize {
        let (nx, ny, _) = self.counts;
        i + nx * (j + (ny + 1) * k)
    }
    fn ihz(&self, i: usize, j: usize, k: usize) -> usize {
        let (nx, ny, _) = self.counts;
        i + nx * (j + ny * k)
    }

    /// Apply each face's condition to the tangential electric field.
    ///
    /// A conducting face is held at zero. An absorbing one takes Mur's update, which needs the
    /// previous step's values at the face and one cell in — so this runs after the interior update
    /// and before anything reads the field.
    fn apply_boundaries(&mut self, dt: f64) {
        let k = {
            let c_dt = self.medium.wave_speed().to_si() * dt;
            (c_dt - self.dx) / (c_dt + self.dx)
        };
        let previous = self.previous.take();
        let (nx, ny, nz) = self.counts;

        // Every face that is conducting keeps the old behaviour, and doing that first means an
        // edge shared by a conducting face and an open one is held — which is the right answer,
        // because a perfect conductor's condition is not negotiable and Mur's is an approximation.
        self.enforce_walls_where(|b| b == Boundary::Conducting);
        self.mirror_magnetic_walls();
        if previous.is_none() {
            self.enforce_walls_where(|b| b == Boundary::Open);
        }

        if let Some(prev) = previous.as_deref() {
            let (pex, pey, pez) = (&prev.ex, &prev.ey, &prev.ez);
            // For each open face, the two tangential components, from the value one cell in.
            for wall in Wall::ALL {
                if self.boundary[wall.index()] != Boundary::Open {
                    continue;
                }
                match wall {
                    Wall::XLow | Wall::XHigh => {
                        let (edge, inner) = if wall == Wall::XLow {
                            (0, 1)
                        } else {
                            (nx, nx - 1)
                        };
                        for kk in 0..=nz {
                            for j in 0..ny {
                                let (a, b) = (self.iey(edge, j, kk), self.iey(inner, j, kk));
                                self.ey[a] = pey[b] + k * (self.ey[b] - pey[a]);
                            }
                        }
                        for kk in 0..nz {
                            for j in 0..=ny {
                                let (a, b) = (self.iez(edge, j, kk), self.iez(inner, j, kk));
                                self.ez[a] = pez[b] + k * (self.ez[b] - pez[a]);
                            }
                        }
                    }
                    Wall::YLow | Wall::YHigh => {
                        let (edge, inner) = if wall == Wall::YLow {
                            (0, 1)
                        } else {
                            (ny, ny - 1)
                        };
                        for kk in 0..=nz {
                            for i in 0..nx {
                                let (a, b) = (self.iex(i, edge, kk), self.iex(i, inner, kk));
                                self.ex[a] = pex[b] + k * (self.ex[b] - pex[a]);
                            }
                        }
                        for kk in 0..nz {
                            for i in 0..=nx {
                                let (a, b) = (self.iez(i, edge, kk), self.iez(i, inner, kk));
                                self.ez[a] = pez[b] + k * (self.ez[b] - pez[a]);
                            }
                        }
                    }
                    Wall::ZLow | Wall::ZHigh => {
                        let (edge, inner) = if wall == Wall::ZLow {
                            (0, 1)
                        } else {
                            (nz, nz - 1)
                        };
                        for j in 0..=ny {
                            for i in 0..nx {
                                let (a, b) = (self.iex(i, j, edge), self.iex(i, j, inner));
                                self.ex[a] = pex[b] + k * (self.ex[b] - pex[a]);
                            }
                        }
                        for j in 0..ny {
                            for i in 0..=nx {
                                let (a, b) = (self.iey(i, j, edge), self.iey(i, j, inner));
                                self.ey[a] = pey[b] + k * (self.ey[b] - pey[a]);
                            }
                        }
                    }
                }
            }
        }
        self.previous = Some(Box::new(Tangential {
            ex: self.ex.clone(),
            ey: self.ey.clone(),
            ez: self.ez.clone(),
        }));
    }

    /// Mirror the tangential electric field across every magnetic-conductor face.
    ///
    /// `∂E_t/∂n = 0`, which is what "no tangential H" means for the electric field. A plane wave
    /// uniform across such a face stays uniform, which is the whole point of having one.
    fn mirror_magnetic_walls(&mut self) {
        let (nx, ny, nz) = self.counts;
        for wall in Wall::ALL {
            if self.boundary[wall.index()] != Boundary::Magnetic {
                continue;
            }
            match wall {
                Wall::XLow | Wall::XHigh => {
                    let (edge, inner) = if wall == Wall::XLow {
                        (0, 1)
                    } else {
                        (nx, nx - 1)
                    };
                    for k in 0..=nz {
                        for j in 0..ny {
                            let (a, b) = (self.iey(edge, j, k), self.iey(inner, j, k));
                            self.ey[a] = self.ey[b];
                        }
                    }
                    for k in 0..nz {
                        for j in 0..=ny {
                            let (a, b) = (self.iez(edge, j, k), self.iez(inner, j, k));
                            self.ez[a] = self.ez[b];
                        }
                    }
                }
                Wall::YLow | Wall::YHigh => {
                    let (edge, inner) = if wall == Wall::YLow {
                        (0, 1)
                    } else {
                        (ny, ny - 1)
                    };
                    for k in 0..=nz {
                        for i in 0..nx {
                            let (a, b) = (self.iex(i, edge, k), self.iex(i, inner, k));
                            self.ex[a] = self.ex[b];
                        }
                    }
                    for k in 0..nz {
                        for i in 0..=nx {
                            let (a, b) = (self.iez(i, edge, k), self.iez(i, inner, k));
                            self.ez[a] = self.ez[b];
                        }
                    }
                }
                Wall::ZLow | Wall::ZHigh => {
                    let (edge, inner) = if wall == Wall::ZLow {
                        (0, 1)
                    } else {
                        (nz, nz - 1)
                    };
                    for j in 0..=ny {
                        for i in 0..nx {
                            let (a, b) = (self.iex(i, j, edge), self.iex(i, j, inner));
                            self.ex[a] = self.ex[b];
                        }
                    }
                    for j in 0..ny {
                        for i in 0..=nx {
                            let (a, b) = (self.iey(i, j, edge), self.iey(i, j, inner));
                            self.ey[a] = self.ey[b];
                        }
                    }
                }
            }
        }
    }

    /// Zero the tangential electric field on the faces a predicate selects.
    fn enforce_walls_where(&mut self, which: impl Fn(Boundary) -> bool) {
        let (nx, ny, nz) = self.counts;
        for wall in Wall::ALL {
            if !which(self.boundary[wall.index()]) {
                continue;
            }
            match wall {
                Wall::XLow | Wall::XHigh => {
                    let edge = if wall == Wall::XLow { 0 } else { nx };
                    for k in 0..=nz {
                        for j in 0..ny {
                            let a = self.iey(edge, j, k);
                            self.ey[a] = 0.0;
                        }
                    }
                    for k in 0..nz {
                        for j in 0..=ny {
                            let a = self.iez(edge, j, k);
                            self.ez[a] = 0.0;
                        }
                    }
                }
                Wall::YLow | Wall::YHigh => {
                    let edge = if wall == Wall::YLow { 0 } else { ny };
                    for k in 0..=nz {
                        for i in 0..nx {
                            let a = self.iex(i, edge, k);
                            self.ex[a] = 0.0;
                        }
                    }
                    for k in 0..nz {
                        for i in 0..=nx {
                            let a = self.iez(i, edge, k);
                            self.ez[a] = 0.0;
                        }
                    }
                }
                Wall::ZLow | Wall::ZHigh => {
                    let edge = if wall == Wall::ZLow { 0 } else { nz };
                    for j in 0..=ny {
                        for i in 0..nx {
                            let a = self.iex(i, j, edge);
                            self.ex[a] = 0.0;
                        }
                    }
                    for j in 0..ny {
                        for i in 0..=nx {
                            let a = self.iey(i, j, edge);
                            self.ey[a] = 0.0;
                        }
                    }
                }
            }
        }
    }

    /// Zero the tangential electric field on every wall.
    fn enforce_walls(&mut self) {
        let (nx, ny, nz) = self.counts;
        // Ex is tangential to the y and z walls.
        for k in 0..=nz {
            for i in 0..nx {
                let (a, b) = (self.iex(i, 0, k), self.iex(i, ny, k));
                self.ex[a] = 0.0;
                self.ex[b] = 0.0;
            }
        }
        for j in 0..=ny {
            for i in 0..nx {
                let (a, b) = (self.iex(i, j, 0), self.iex(i, j, nz));
                self.ex[a] = 0.0;
                self.ex[b] = 0.0;
            }
        }
        // Ey is tangential to the x and z walls.
        for k in 0..=nz {
            for j in 0..ny {
                let (a, b) = (self.iey(0, j, k), self.iey(nx, j, k));
                self.ey[a] = 0.0;
                self.ey[b] = 0.0;
            }
        }
        for j in 0..ny {
            for i in 0..=nx {
                let (a, b) = (self.iey(i, j, 0), self.iey(i, j, nz));
                self.ey[a] = 0.0;
                self.ey[b] = 0.0;
            }
        }
        // Ez is tangential to the x and y walls.
        for k in 0..nz {
            for j in 0..=ny {
                let (a, b) = (self.iez(0, j, k), self.iez(nx, j, k));
                self.ez[a] = 0.0;
                self.ez[b] = 0.0;
            }
        }
        for k in 0..nz {
            for i in 0..=nx {
                let (a, b) = (self.iez(i, 0, k), self.iez(i, ny, k));
                self.ez[a] = 0.0;
                self.ez[b] = 0.0;
            }
        }
    }

    /// Advance `H` by `dt`, from `E`.
    fn advance_magnetic(&mut self, dt: f64) {
        let (nx, ny, nz) = self.counts;
        let g = dt / self.dx;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..=nx {
                    let curl = (self.ez[self.iez(i, j + 1, k)] - self.ez[self.iez(i, j, k)])
                        - (self.ey[self.iey(i, j, k + 1)] - self.ey[self.iey(i, j, k)]);
                    let idx = self.ihx(i, j, k);
                    self.hx[idx] -= g / self.mu[0][idx] * curl;
                }
            }
        }
        for k in 0..nz {
            for j in 0..=ny {
                for i in 0..nx {
                    let curl = (self.ex[self.iex(i, j, k + 1)] - self.ex[self.iex(i, j, k)])
                        - (self.ez[self.iez(i + 1, j, k)] - self.ez[self.iez(i, j, k)]);
                    let idx = self.ihy(i, j, k);
                    self.hy[idx] -= g / self.mu[1][idx] * curl;
                }
            }
        }
        for k in 0..=nz {
            for j in 0..ny {
                for i in 0..nx {
                    let curl = (self.ey[self.iey(i + 1, j, k)] - self.ey[self.iey(i, j, k)])
                        - (self.ex[self.iex(i, j + 1, k)] - self.ex[self.iex(i, j, k)]);
                    let idx = self.ihz(i, j, k);
                    self.hz[idx] -= g / self.mu[2][idx] * curl;
                }
            }
        }
    }

    /// Advance `E` by `dt`, from `H`. Interior edges only; the walls are held by
    /// [`Cavity::enforce_walls`].
    fn advance_electric(&mut self, dt: f64) -> f64 {
        let (nx, ny, nz) = self.counts;
        // The semi-implicit form, which is stable for any conductivity rather than only for
        // `dt < 2ε/σ`. An explicit `E -= (σ/ε)E dt` is the obvious version and it explodes in a
        // good conductor at any step this grid would use.
        //
        // Per component rather than per cavity, because a field on the edge between two materials
        // has the average of both and not the box's own.
        let coefficients = |eps: f64, sigma: f64| -> (f64, f64) {
            let half = sigma * dt / (2.0 * eps);
            (
                (1.0 - half) / (1.0 + half),
                (dt / (eps * self.dx)) / (1.0 + half),
            )
        };
        // **Joule heating, `∫σE²dt`, and not the change in the field's own energy.** The first
        // version took `½ε(E² − E'²)`, which is the loss *plus* whatever flowed in from the
        // magnetic field on the same step — and reported dissipating thirty-six times the energy
        // the cavity ever held. `E` at the midpoint is the semi-implicit scheme's own effective
        // value, which is what makes this the second-order figure rather than a first-order one.
        let mut lost = 0.0;
        let volume = self.dx.powi(3);

        for k in 1..nz {
            for j in 1..ny {
                for i in 0..nx {
                    let curl = (self.hz[self.ihz(i, j, k)] - self.hz[self.ihz(i, j - 1, k)])
                        - (self.hy[self.ihy(i, j, k)] - self.hy[self.ihy(i, j, k - 1)]);
                    let idx = self.iex(i, j, k);
                    let (eps, sigma) = (self.eps[0][idx], self.sigma[0][idx]);
                    let (decay, gain) = coefficients(eps, sigma);
                    let before = self.ex[idx];
                    self.ex[idx] = decay * before + gain * curl;
                    let mid = 0.5 * (before + self.ex[idx]);
                    lost += sigma * mid * mid * dt * volume;
                }
            }
        }
        for k in 1..nz {
            for j in 0..ny {
                for i in 1..nx {
                    let curl = (self.hx[self.ihx(i, j, k)] - self.hx[self.ihx(i, j, k - 1)])
                        - (self.hz[self.ihz(i, j, k)] - self.hz[self.ihz(i - 1, j, k)]);
                    let idx = self.iey(i, j, k);
                    let (eps, sigma) = (self.eps[1][idx], self.sigma[1][idx]);
                    let (decay, gain) = coefficients(eps, sigma);
                    let before = self.ey[idx];
                    self.ey[idx] = decay * before + gain * curl;
                    let mid = 0.5 * (before + self.ey[idx]);
                    lost += sigma * mid * mid * dt * volume;
                }
            }
        }
        for k in 0..nz {
            for j in 1..ny {
                for i in 1..nx {
                    let curl = (self.hy[self.ihy(i, j, k)] - self.hy[self.ihy(i - 1, j, k)])
                        - (self.hx[self.ihx(i, j, k)] - self.hx[self.ihx(i, j - 1, k)]);
                    let idx = self.iez(i, j, k);
                    let (eps, sigma) = (self.eps[2][idx], self.sigma[2][idx]);
                    let (decay, gain) = coefficients(eps, sigma);
                    let before = self.ez[idx];
                    self.ez[idx] = decay * before + gain * curl;
                    let mid = 0.5 * (before + self.ez[idx]);
                    lost += sigma * mid * mid * dt * volume;
                }
            }
        }
        lost.max(0.0)
    }
}

impl Domain for Cavity {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    fn max_stable_dt(&self, _now: Time) -> Time {
        self.courant_limit()
    }

    fn step(&mut self, _t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let h = dt.to_si();
        let limit = self.courant_limit().to_si();
        if h > limit {
            return Err(Violation {
                quantity: "Courant number".into(),
                site: self.name.clone(),
                before: limit,
                after: h,
                scale: limit,
                tolerance: 0.0,
            });
        }
        // The leapfrog's opening half step, which is the defect this workspace already paid for
        // once: `H` is defined half a step behind `E`, and a first update of a whole step puts it
        // a half step wrong for the rest of the run. `release_mode` sets `H` at `−dt/2` already
        // and says so; a cavity that was driven rather than released starts from rest, where the
        // half step costs nothing because there is nothing to be wrong about.
        if !self.started {
            self.advance_magnetic(0.5 * h);
            self.started = true;
        }
        // The two half-step-separated `H` fields, paired in the same statement that produces the
        // second. Nothing else can pair them: after this line the earlier one is gone.
        let before = (self.hx.clone(), self.hy.clone(), self.hz.clone());
        self.advance_magnetic(h);
        // Here, and only here: `E` is still `Eⁿ`, `H` is `Hⁿ⁺¹ᐟ²` and `before` is `Hⁿ⁻¹ᐟ²`.
        self.invariant =
            Some(self.electric_energy().to_si() + self.paired_magnetic_energy(&before));
        self.dissipated += self.advance_electric(h);
        self.apply_boundaries(h);
        self.enforce_solid();
        Ok(())
    }

    /// The field energy, plus whatever a lossy medium has turned into heat.
    ///
    /// Arranged so the total is constant: in a lossless box nothing leaves, and in a lossy one what
    /// the field lost is on the books beside it rather than gone.
    /// The **leapfrog invariant** plus what a lossy medium has turned into heat.
    ///
    /// Not `energy()`, and the difference is the point. `½εE² + ½μH²` is the physical field energy
    /// and is what a reader wants to see; it is *not* what this scheme conserves, because `E` and
    /// `H` are half a step apart. The naive sum oscillates about the true invariant by
    /// `2 sin(ωΔt/2)` peak to peak — a closed form this crate's own test checks, and 10.7% for a
    /// 1.77 GHz mode at the Courant limit. An audit watching it would stop every correct cavity
    /// run, and widening the scene's budget past 10% to let it through would be switching the
    /// audit off over something that has an exact expression.
    ///
    /// Two contributions rather than their sum, so `Ledger`'s scale is the larger of the two: a
    /// cavity that has dissipated most of what it held would otherwise have a near-zero net and no
    /// scale for the audit to judge against.
    fn ledger(&self) -> Ledger {
        Ledger::new()
            .with(quantity::ENERGY, self.invariant_energy().to_si())
            .with(quantity::ENERGY, self.dissipated)
    }

    fn readings(&self) -> Vec<Reading> {
        vec![
            Reading::new(&self.name, "field energy", self.energy().to_si(), "J"),
            // Beside it, because a reader comparing the two is the shortest way to see that the
            // wobble in the first is the scheme's staggering and not a leak.
            Reading::new(
                &self.name,
                "invariant",
                self.invariant_energy().to_si(),
                "J",
            ),
            Reading::new(&self.name, "electric", self.electric_energy().to_si(), "J"),
            Reading::new(&self.name, "magnetic", self.magnetic_energy().to_si(), "J"),
            Reading::new(&self.name, "div B", self.magnetic_divergence(), ""),
            Reading::new(&self.name, "dissipated", self.dissipated, "J"),
        ]
    }

    /// **How strong the electric field is**, so a resonance can be looked at.
    ///
    /// The caveat is the one a magnitude always carries: `E` is a vector on cell edges and this is
    /// its length at the cell centre, so a picture shows the standing pattern — which is what a
    /// mode *is* — and not the polarisation. `electric_at` gives the vector.
    fn as_field(&self) -> Option<&dyn pantometry_core::ScalarField> {
        Some(self)
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    /// **And mutably**, because a domain that can be read and not written is one a coupling can
    /// only fail at silently. `Simulation::domain_as_mut` returns `None` when this is not
    /// implemented, and a caller that wrote through it would do nothing and report nothing —
    /// which is how a whole coupled body came back reporting zero strain that read as *no stress*
    /// rather than as *not connected*.
    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }

    fn checkpoint(&mut self) {
        self.saved = Some(Box::new(Saved {
            ex: self.ex.clone(),
            ey: self.ey.clone(),
            ez: self.ez.clone(),
            hx: self.hx.clone(),
            hy: self.hy.clone(),
            hz: self.hz.clone(),
            started: self.started,
            dissipated: self.dissipated,
        }));
    }

    fn restore(&mut self) {
        if let Some(s) = self.saved.take() {
            self.ex = s.ex;
            self.ey = s.ey;
            self.ez = s.ez;
            self.hx = s.hx;
            self.hy = s.hy;
            self.hz = s.hz;
            self.started = s.started;
            self.dissipated = s.dissipated;
        }
    }
}
