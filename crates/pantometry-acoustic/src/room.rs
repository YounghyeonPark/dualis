//! Sound in two dimensions: a room rather than a tube.
//!
//! [`Tube`](crate::Tube) is one-dimensional, which covers a duct and an organ pipe and
//! nothing with a shape. A room has modes in every direction at once, and they interact
//! in a way a one-dimensional model cannot show: the frequencies are not a harmonic
//! series, they crowd together as they rise, and the low ones are far enough apart to be
//! heard individually as the boom a small room has at one particular note.
//!
//! # The closed form, which is exact and not obvious
//!
//! A rectangular room with rigid walls resonates at
//!
//! ```text
//! f(nx, ny) = (c/2) √((nx/Lx)² + (ny/Ly)²)
//! ```
//!
//! for every pair of non-negative integers. Two things follow that a tube does not show.
//! The series is not harmonic — `f(1,1)` is not a whole multiple of `f(1,0)` — so a room
//! does not ring on a note, it rings on a chord that is not one. And the modes get denser
//! as the frequency rises, going as `f²` in two dimensions rather than staying evenly
//! spaced, which is why a room's colouration is audible at the bottom of the spectrum and
//! not at the top.
//!
//! # What this shares with the tube, and what it does not
//!
//! The same staggered scheme, extended: pressures at cell centres, velocities on the
//! faces between them in each direction. The CFL condition tightens, and by a factor that
//! is worth knowing — `dt ≤ dx/(c√2)` in two dimensions rather than `dx/c`, because a
//! wave travelling diagonally covers `√2` cells while crossing one. In three dimensions
//! it would be `√3`. Every explicit wave solver pays that, and it is the reason
//! three-dimensional acoustics is expensive rather than merely large.

use glam::DVec3;
use pantometry_core::conserved::quantity;
use pantometry_core::{Domain, Exchange, Kind, Ledger, Reading, ScalarField, Violation};
use pantometry_units::{
    Area, Density, Energy, Frequency, Length, LengthVec, Pressure, Time, Velocity,
};

/// Everything [`Room::checkpoint`] has to put back: the pressure, both velocity components,
/// and whether the velocity has been staggered yet.
///
/// A named type because clippy is right that four things in a tuple is where a reader starts
/// counting commas.
type Snapshot = (Vec<f64>, Vec<f64>, Vec<f64>, bool);

/// A rectangular room, discretised on a uniform grid with rigid walls.
pub struct Room {
    name: String,
    /// Pressure at cell centres, row-major.
    pressure: Vec<f64>,
    /// Pressure one whole step earlier, kept only so that [`Room::energy`] can be the
    /// quantity the scheme actually conserves rather than one that wobbles.
    pressure_prev: Vec<f64>,
    /// Velocity on the faces between horizontally adjacent cells: `(nx - 1)` per row.
    vx: Vec<f64>,
    /// Velocity on the faces between vertically adjacent cells: `(ny - 1)` rows.
    vy: Vec<f64>,
    nx: usize,
    ny: usize,
    dx: f64,
    speed: f64,
    density: f64,
    /// Depth, so a two-dimensional room still has a volume and an energy in joules.
    depth: f64,
    /// Whether `vx`/`vy` have been offset to the half step the scheme carries them at.
    ///
    /// A staggered leapfrog wants velocity half a step *behind* the pressure it is about to
    /// update: `v^(n+1/2) = v^(n-1/2) - (h/rho) grad p^n`. An initial condition gives
    /// velocity at `t = 0`, not at `t = -h/2`, so the very first velocity update has only
    /// half a step to travel. Using a whole one leaves an `O(h)` error that never decays —
    /// and since `h` follows `dx` through the CFL condition, it drags a second-order scheme
    /// to first order. See the note on [`Room::released_from`].
    velocity_staggered: bool,
    /// What the released state's energy was, minus what the scheme's invariant came out as
    /// once it had been staggered. See [`Room::startup_adjustment`].
    energy_offset: f64,
    /// The physical energy at release, which is the datum the reported energy is measured
    /// against.
    energy_datum: f64,
    saved: Option<Snapshot>,
}

impl Room {
    /// A room of the given size, at rest, with rigid walls.
    ///
    /// The grid spacing is taken from the width; the height is quantised to the nearest
    /// whole number of the same cells, so the cells stay square and the CFL limit stays
    /// isotropic.
    pub fn new(
        name: impl Into<String>,
        width: Length,
        height: Length,
        cells_across: usize,
        depth: Length,
        density: Density,
        sound_speed: Velocity,
    ) -> Room {
        let nx = cells_across.max(3);
        let dx = width.to_si() / (nx - 1) as f64;
        let ny = ((height.to_si() / dx).round() as usize + 1).max(3);
        Room {
            name: name.into(),
            pressure: vec![0.0; nx * ny],
            pressure_prev: vec![0.0; nx * ny],
            vx: vec![0.0; (nx - 1) * ny],
            vy: vec![0.0; nx * (ny - 1)],
            nx,
            ny,
            dx,
            speed: sound_speed.to_si(),
            density: density.to_si(),
            depth: depth.to_si(),
            velocity_staggered: false,
            energy_offset: 0.0,
            energy_datum: 0.0,
            saved: None,
        }
    }

    /// A room of air at 20 °C, one metre deep.
    pub fn of_air(
        name: impl Into<String>,
        width: Length,
        height: Length,
        cells_across: usize,
    ) -> Room {
        Room::new(
            name,
            width,
            height,
            cells_across,
            Length::m(1.0),
            Density::kg_per_m3(1.204),
            Velocity::m_per_s(343.0),
        )
    }

    /// Start with a pressure field and no motion.
    pub fn released_from(mut self, profile: impl Fn(Length, Length) -> Pressure) -> Room {
        for j in 0..self.ny {
            for i in 0..self.nx {
                let x = Length::from_si(i as f64 * self.dx);
                let y = Length::from_si(j as f64 * self.dx);
                self.pressure[j * self.nx + i] = profile(x, y).to_si();
            }
        }
        self.pressure_prev = self.pressure.clone();
        self.vx.iter_mut().for_each(|v| *v = 0.0);
        self.vy.iter_mut().for_each(|v| *v = 0.0);
        // At rest at `t = 0`, which is not where the scheme carries velocity. The first step
        // makes up the difference; see `velocity_staggered`.
        self.velocity_staggered = false;
        // The datum: the energy this state physically has, with the velocity where the
        // initial condition put it. What the scheme will conserve is a slightly different
        // number, and `startup_adjustment` is the difference.
        self.energy_offset = 0.0;
        self.energy_datum = self.invariant_si();
        self
    }

    /// Excite one rigid-wall mode exactly: `cos(nx π x/Lx) cos(ny π y/Ly)`.
    ///
    /// The shape whose oscillation frequency the closed form predicts, so that the
    /// prediction can be checked against the integration rather than against itself.
    pub fn released_in_mode(self, nx: u32, ny: u32, amplitude: Pressure) -> Room {
        let (lx, ly) = (self.width().to_si(), self.height().to_si());
        let a = amplitude.to_si();
        self.released_from(|x, y| {
            let px = if lx > 0.0 {
                (nx as f64 * std::f64::consts::PI * x.to_si() / lx).cos()
            } else {
                1.0
            };
            let py = if ly > 0.0 {
                (ny as f64 * std::f64::consts::PI * y.to_si() / ly).cos()
            } else {
                1.0
            };
            Pressure::from_si(a * px * py)
        })
    }

    /// Width of the room. `(nx − 1)·dx`, because the samples sit on the walls.
    pub fn width(&self) -> Length {
        Length::from_si((self.nx - 1) as f64 * self.dx)
    }

    /// Height of the room.
    pub fn height(&self) -> Length {
        Length::from_si((self.ny - 1) as f64 * self.dx)
    }

    /// Samples across and up, as `(nx, ny)`.
    pub fn cells(&self) -> (usize, usize) {
        (self.nx, self.ny)
    }

    /// Acoustic pressure at one node, clamped to the walls. Signed.
    pub fn pressure_at(&self, i: usize, j: usize) -> Pressure {
        let i = i.min(self.nx - 1);
        let j = j.min(self.ny - 1);
        Pressure::from_si(self.pressure[j * self.nx + i])
    }

    /// Largest pressure magnitude anywhere in the room.
    pub fn peak_pressure(&self) -> Pressure {
        Pressure::from_si(self.pressure.iter().fold(0.0f64, |a, p| a.max(p.abs())))
    }

    /// Acoustic energy, in the form the scheme actually conserves.
    ///
    /// The obvious expression — `p²` at one time level plus `u²` at another — is not it.
    /// On a staggered grid the pressures sit at whole steps and the velocities half a
    /// step later, so reading both at once carries an `O(ω dt)` wobble: measured at 4% for
    /// a pulse at five substeps a millisecond, which does not drift but is more than
    /// enough to make a conservation audit meaningless.
    ///
    /// The time-centred form uses the *product* of the pressure either side of the step,
    /// `pⁿ pⁿ⁺¹`, against the velocity at the half step between them. That combination is
    /// exactly conserved by the lossless leapfrog rather than approximately, which is what
    /// lets the audit be tight. It costs one extra array.
    ///
    /// A consequence worth knowing: a single cell's contribution can be negative, where
    /// the pressure changed sign across the step. The total cannot.
    pub fn energy(&self) -> Energy {
        Energy::from_si(self.invariant_si() + self.energy_offset)
    }

    /// How much the discrete invariant differs from the released state's physical energy.
    ///
    /// Zero until the first step, because the invariant is a function of the step size and
    /// there is no step size before then. `O(h²)` and therefore converging: 0.42% of the
    /// total on a 31-cell room, 0.10% at 61, 0.026% at 121, 0.0065% at 241.
    ///
    /// Reported rather than hidden, because it is the one place where what this domain tells
    /// the audit is not simply what its state says.
    pub fn startup_adjustment(&self) -> Energy {
        Energy::from_si(self.energy_offset)
    }

    /// The scheme's invariant itself, before the release datum is applied.
    fn invariant_si(&self) -> f64 {
        let volume = self.dx * self.dx * self.depth;
        let rc2 = self.density * self.speed * self.speed;
        if rc2 <= 0.0 {
            return 0.0;
        }
        // Every term carries the volume it actually occupies. A node on a wall owns half a
        // cell, one in a corner a quarter, and a face between two boundary rows likewise —
        // the same weights the update uses, which is what keeps the conservation exact
        // rather than approximate.
        let (nx, ny) = (self.nx, self.ny);
        let potential: f64 = self
            .pressure
            .iter()
            .zip(self.pressure_prev.iter())
            .enumerate()
            .map(|(k, (p, prev))| {
                let (i, j) = (k % nx, k / nx);
                let share = 1.0 / (wall_weight(i, nx) * wall_weight(j, ny));
                share * p * prev / (2.0 * rc2) * volume
            })
            .sum();
        // A vx face spans a whole cell in x and sits in row `j`, so it takes that row's
        // share; a vy face is the other way round.
        let kinetic_x: f64 = self
            .vx
            .iter()
            .enumerate()
            .map(|(k, u)| self.density * u * u / 2.0 * volume / wall_weight(k / (nx - 1), ny))
            .sum();
        let kinetic_y: f64 = self
            .vy
            .iter()
            .enumerate()
            .map(|(k, u)| self.density * u * u / 2.0 * volume / wall_weight(k % nx, nx))
            .sum();
        potential + kinetic_x + kinetic_y
    }

    /// Cross-sectional area of one cell, for reporting.
    pub fn cell_area(&self) -> Area {
        Area::from_si(self.dx * self.depth)
    }

    /// The `(nx, ny)` rigid-wall mode frequency: `(c/2)√((nx/Lx)² + (ny/Ly)²)`.
    ///
    /// Exact, and not a harmonic series: `f(1,1)` of a square room is `√2` times
    /// `f(1,0)`, which is a tritone and not a note anyone would call in tune.
    pub fn mode_frequency(&self, nx: u32, ny: u32) -> Frequency {
        let (lx, ly) = (self.width().to_si(), self.height().to_si());
        if lx <= 0.0 || ly <= 0.0 {
            return Frequency::from_si(0.0);
        }
        let a = nx as f64 / lx;
        let b = ny as f64 / ly;
        Frequency::from_si(self.speed / 2.0 * (a * a + b * b).sqrt())
    }

    /// Every mode below a frequency, in ascending order.
    ///
    /// The count grows as `f²` in two dimensions — the modes fill an area of the `(nx,
    /// ny)` lattice — which is why a room's individual resonances are audible at the
    /// bottom of the range and merge into a statistical hiss further up.
    pub fn modes_below(&self, limit: Frequency) -> Vec<(u32, u32, Frequency)> {
        let mut found = Vec::new();
        // Bound the search from the largest index either axis could contribute.
        let max_index = |l: f64| {
            if l <= 0.0 {
                0
            } else {
                (2.0 * limit.to_si() * l / self.speed).floor() as u32
            }
        };
        let (mx, my) = (
            max_index(self.width().to_si()),
            max_index(self.height().to_si()),
        );
        for nx in 0..=mx {
            for ny in 0..=my {
                if nx == 0 && ny == 0 {
                    continue;
                }
                let f = self.mode_frequency(nx, ny);
                if f <= limit {
                    found.push((nx, ny, f));
                }
            }
        }
        found.sort_by(|a, b| a.2.to_si().total_cmp(&b.2.to_si()));
        found
    }

    /// Courant number, `c dt √2/dx` in two dimensions.
    pub fn courant(&self, dt: Time) -> f64 {
        self.speed * dt.to_si() * 2f64.sqrt() / self.dx
    }
}

impl Domain for Room {
    fn books_balance(&self) -> bool {
        true
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// `dx/(c√2)`, tighter than the tube's `dx/c` by exactly the diagonal.
    ///
    /// A wave going corner to corner covers `√2` cells while crossing one, and the
    /// three-point stencil in each direction cannot know about it. In three dimensions the
    /// factor is `√3`, so the same room at the same resolution costs 22% more steps again
    /// on top of the extra dimension — which is the real reason three-dimensional
    /// acoustics is dear.
    fn max_stable_dt(&self, _now: Time) -> Time {
        if self.speed <= 0.0 {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(self.dx / (self.speed * 2f64.sqrt()))
    }

    fn step(&mut self, _t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let courant = self.courant(dt);
        if courant > 1.0 + 1e-12 {
            return Err(Violation {
                quantity: "Courant number".to_string(),
                site: format!("{} (two-dimensional wave equation)", self.name),
                before: 1.0,
                after: courant,
                scale: 1.0,
                tolerance: 1e-12,
            });
        }
        let h = dt.to_si();
        if h <= 0.0 {
            return Ok(());
        }
        let rc2 = self.density * self.speed * self.speed;
        let (nx, ny) = (self.nx, self.ny);

        // Half a step the first time, a whole one after that. The stored velocity starts at
        // `t = 0` where the initial condition put it, and the scheme wants it at `t = -h/2`,
        // so the first update has half the distance to cover. Kicking it a whole step
        // instead is the classic leapfrog startup error: `O(h)`, permanent, and enough to
        // make a second-order scheme converge at first order.
        let vh = if self.velocity_staggered { h } else { 0.5 * h };
        let starting = !self.velocity_staggered;
        self.velocity_staggered = true;

        // Faces: rho du/dt = -grad p.
        for j in 0..ny {
            for i in 0..nx - 1 {
                let face = j * (nx - 1) + i;
                self.vx[face] -= vh / (self.density * self.dx)
                    * (self.pressure[j * nx + i + 1] - self.pressure[j * nx + i]);
            }
        }
        for j in 0..ny - 1 {
            for i in 0..nx {
                let face = j * nx + i;
                self.vy[face] -= vh / (self.density * self.dx)
                    * (self.pressure[(j + 1) * nx + i] - self.pressure[j * nx + i]);
            }
        }

        // Cells: dp/dt = -rho c^2 div u. Rigid walls mean no velocity through them, so
        // the flux outside the domain is zero — which is the boundary condition and needs
        // no ghost cells.
        self.pressure_prev.copy_from_slice(&self.pressure);
        for j in 0..ny {
            for i in 0..nx {
                let left = if i == 0 {
                    0.0
                } else {
                    self.vx[j * (nx - 1) + i - 1]
                };
                let right = if i == nx - 1 {
                    0.0
                } else {
                    self.vx[j * (nx - 1) + i]
                };
                let below = if j == 0 {
                    0.0
                } else {
                    self.vy[(j - 1) * nx + i]
                };
                let above = if j == ny - 1 {
                    0.0
                } else {
                    self.vy[j * nx + i]
                };
                // A wall node's control volume is only half a cell wide in the direction
                // normal to that wall, so its divergence is divided by half the spacing.
                // Getting this wrong makes the walls twice as heavy as they are and drops
                // the whole scheme to first order — see `wall_weight`.
                self.pressure[j * nx + i] -= h * rc2 / self.dx
                    * (wall_weight(i, nx) * (right - left) + wall_weight(j, ny) * (above - below));
            }
        }

        // Startup, once: the released state has its velocity at t = 0 and the scheme's
        // invariant is defined on a staggered state, so converting between them moves the
        // number by O(h²). Measured on a 61-cell room: 0.10%, and it quarters on refinement
        // — 0.42% at 31 cells, 0.026% at 121, 0.0065% at 241.
        //
        // That is discretisation, not a leak: from this step onward the invariant holds to
        // about 1e-15. So the energy is reported against the released state as its datum and
        // the difference is kept here, where `startup_adjustment` can be asked for it, rather
        // than being handed to the audit as a violation it would be right to refuse.
        //
        // The bound is what stops this from being a place for a real bug to hide. A genuine
        // error in the first step — a sign, a factor of two, a boundary — is not O(h²), and
        // a quarter of the energy is far outside anything the Courant condition permits.
        if starting {
            let invariant = self.invariant_si();
            let offset = self.energy_datum - invariant;
            let scale = self.energy_datum.abs().max(invariant.abs());
            if scale > 0.0 && offset.abs() / scale > 0.25 {
                return Err(Violation {
                    quantity: "energy".to_string(),
                    site: format!("{} (starting the leapfrog)", self.name),
                    before: self.energy_datum,
                    after: invariant,
                    scale,
                    tolerance: 0.25,
                });
            }
            self.energy_offset = offset;
        }

        Ok(())
    }

    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.energy().to_si())
    }

    fn checkpoint(&mut self) {
        self.saved = Some((
            self.pressure.clone(),
            self.vx.clone(),
            self.vy.clone(),
            self.velocity_staggered,
        ));
    }

    fn restore(&mut self) {
        if let Some((p, vx, vy, staggered)) = self.saved.clone() {
            self.pressure_prev.copy_from_slice(&p);
            self.pressure = p;
            self.vx = vx;
            self.vy = vy;
            self.velocity_staggered = staggered;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    /// The peak pressure, which is what a mode's decay shows.
    ///
    /// Not the mean: a pressure field's mean is zero by symmetry, and a column of it would be a
    /// column of rounding.
    fn readings(&self) -> Vec<Reading> {
        let peak = self.pressure.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        vec![Reading::new(&self.name, "peak", peak, "Pa")]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }

    /// The room reads as a pressure field, so a renderer never has to know it is a room.
    fn as_field(&self) -> Option<&dyn pantometry_core::ScalarField> {
        Some(self)
    }
}

/// The room as a pressure field, in pascals.
///
/// The second implementation of [`ScalarField`] in the workspace, and it exists to answer
/// the question the first one could not: whether the trait fits anything but a
/// one-dimensional diffusion. `Bar1D` is a line of cells governed by an equation that is
/// first order in time; this is a plane of nodes governed by one that is second order. They
/// disagree in three ways that are worth writing down, because a visualiser will meet all
/// three.
///
/// # `rate` needs the velocity, not the Laplacian
///
/// For diffusion the governing equation gives `∂T/∂t = α∇²T`, so the field's own curvature
/// *is* its rate of change and one array answers everything. A wave is second order:
/// `∂²p/∂t² = c²∇²p`. The Laplacian gives the acceleration, and the *rate* comes from the
/// companion velocity field instead — `∂p/∂t = −ρc²∇·u`, which this scheme stores on the
/// faces and which [`Domain::step`] uses verbatim.
///
/// So both domains can report an exact rate from their present state and neither needs
/// history, but by different routes. The general statement is the useful one: a marched
/// domain can answer `rate` whenever its stored state is enough to evaluate the governing
/// equation, which is the same condition as being well posed as a first-order system. It is
/// not a property of diffusion.
///
/// # The gradient really is zero at these walls
///
/// The thermal crate's `Bar1D` — named in prose rather than linked, because these two
/// domains do not depend on each other and rustdoc rightly refuses to pretend otherwise —
/// uses a cell-centred grid, which puts its first sample half a cell
/// inside the wall, where the temperature is still changing, so a gradient reported there is
/// not zero. This grid is node-centred: a node sits *on* the wall, and a rigid wall's
/// boundary condition is exactly `∂p/∂n = 0`. Mirroring gives `(p₁ − p₁)/2dx`, which is zero
/// to the last bit.
///
/// Neither is a mistake. The difference is where the samples are, and that is the kind of
/// thing only a second implementor could have shown.
///
/// # And the values are signed
///
/// A temperature is positive and climbs; a pressure swings either side of zero and averages
/// out. Any colour map worth having for this one is diverging with a fixed midpoint, and one
/// built for temperature will draw a room as though half of it were cold. That is a
/// requirement on a renderer that a single monotonic field would never have surfaced.
///
/// # One sharp edge to know about
///
/// `at` interpolates between nodes; `gradient`, `laplacian` and `rate` snap to the nearest
/// one, because they are the scheme's own stencils and the scheme only has values at nodes.
/// So combining them at an arbitrary point — `laplacian(p) / at(p)`, say — divides a number
/// computed at one place by one computed at another. Sample on a node when the ratio
/// matters. The same is true of any field backed by a grid.
///
/// Depth is ignored: this is a two-dimensional model and says so.
impl ScalarField for Room {
    fn unit(&self) -> &'static str {
        "Pa"
    }

    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let v = p.to_si();
        let (u, w) = (v.x / self.dx, v.y / self.dx);
        // Bilinear between the four surrounding nodes, clamped to the walls.
        let (i, fx) = clamp_index(u, self.nx);
        let (j, fy) = clamp_index(w, self.ny);
        let (i1, j1) = ((i + 1).min(self.nx - 1), (j + 1).min(self.ny - 1));
        let p00 = self.pressure[j * self.nx + i];
        let p10 = self.pressure[j * self.nx + i1];
        let p01 = self.pressure[j1 * self.nx + i];
        let p11 = self.pressure[j1 * self.nx + i1];
        let bottom = p00 * (1.0 - fx) + p10 * fx;
        let top = p01 * (1.0 - fx) + p11 * fx;
        bottom * (1.0 - fy) + top * fy
    }

    /// `∇p`, by the mirrored central difference the walls impose. Exactly zero on a wall,
    /// in the direction normal to it.
    fn gradient(&self, p: LengthVec, _t: Time, _h: Length) -> DVec3 {
        let (i, j) = self.node_at(p);
        let sample = |i: usize, j: usize| self.pressure[j * self.nx + i];
        let mirror = |k: usize, n: usize| {
            (
                if k == 0 { 1 } else { k - 1 },
                if k + 1 >= n {
                    n.saturating_sub(2)
                } else {
                    k + 1
                },
            )
        };
        let (il, ir) = mirror(i, self.nx);
        let (jl, jr) = mirror(j, self.ny);
        DVec3::new(
            (sample(ir, j) - sample(il, j)) / (2.0 * self.dx),
            (sample(i, jr) - sample(i, jl)) / (2.0 * self.dx),
            0.0,
        )
    }

    /// `∇²p`, five-point with rigid walls. This is `∂²p/∂t²` over `c²` — the acceleration,
    /// not the rate; see [`Room::rate`](ScalarField::rate).
    fn laplacian(&self, p: LengthVec, _t: Time, _h: Length) -> f64 {
        let (i, j) = self.node_at(p);
        let sample = |i: usize, j: usize| self.pressure[j * self.nx + i];
        let mirror = |k: usize, n: usize| {
            (
                if k == 0 { 1 } else { k - 1 },
                if k + 1 >= n {
                    n.saturating_sub(2)
                } else {
                    k + 1
                },
            )
        };
        let (il, ir) = mirror(i, self.nx);
        let (jl, jr) = mirror(j, self.ny);
        let centre = sample(i, j);
        (sample(il, j) + sample(ir, j) + sample(i, jl) + sample(i, jr) - 4.0 * centre)
            / (self.dx * self.dx)
    }

    /// `∂p/∂t = −ρc²∇·u`, read off the stored face velocities.
    ///
    /// The same expression [`Domain::step`] integrates, on the same stencil, with the same
    /// zero flux through the walls — so it is exact rather than approximate. But it is exact
    /// about the step just *taken*, not the one coming, and that is worth being precise
    /// about because the bar is the other way round.
    ///
    /// A leapfrog stores its velocities half a step behind its pressures. The update is
    /// `pⁿ⁺¹ = pⁿ − h·ρc²∇·uⁿ⁺¹⁄²`, and `uⁿ⁺¹⁄²` does not exist until the next step computes
    /// it — what is in the arrays now is `uⁿ⁻¹⁄²`. So this returns `(pⁿ − pⁿ⁻¹)/h` exactly:
    /// the centred derivative at `t − dt/2`.
    ///
    /// Being half a step stale is the better trade, not a compromise. A centred difference
    /// at the midpoint is second-order accurate where a forward one at the present instant
    /// is first-order, and half a step is nothing next to a period — twenty microseconds
    /// against a fiftieth of a second for a room mode.
    ///
    /// Nearest node rather than interpolated, because the divergence is a cell quantity and
    /// inventing values between cells would report a rate the scheme is not going to
    /// produce.
    fn rate(&self, p: LengthVec, _t: Time, _dt: Time) -> f64 {
        let (i, j) = self.node_at(p);
        let (nx, ny) = (self.nx, self.ny);
        let left = if i == 0 {
            0.0
        } else {
            self.vx[j * (nx - 1) + i - 1]
        };
        let right = if i == nx - 1 {
            0.0
        } else {
            self.vx[j * (nx - 1) + i]
        };
        let below = if j == 0 {
            0.0
        } else {
            self.vy[(j - 1) * nx + i]
        };
        let above = if j == ny - 1 {
            0.0
        } else {
            self.vy[j * nx + i]
        };
        // The same wall weighting the update uses, or this would report a rate the scheme
        // is not going to produce — and only at the boundary, where it is hardest to see.
        -self.density * self.speed * self.speed / self.dx
            * (wall_weight(i, nx) * (right - left) + wall_weight(j, ny) * (above - below))
    }
}

impl Room {
    /// Nearest node to a point, clamped into the room.
    fn node_at(&self, p: LengthVec) -> (usize, usize) {
        let v = p.to_si();
        let round = |q: f64, n: usize| {
            if q.is_nan() || q < 0.0 {
                0
            } else {
                (q / self.dx).round().min((n - 1) as f64) as usize
            }
        };
        (round(v.x, self.nx), round(v.y, self.ny))
    }
}

/// How much of a whole cell's spacing a node's control volume spans, inverted.
///
/// `1` inside, `2` on a wall — because a wall node owns half a cell in that direction, so
/// its divergence is divided by `dx/2` rather than `dx`. The same number appears reciprocally
/// in the energy, where that node holds half as much.
///
/// Leaving it out is what made every mode of this room read low. The interior scheme is
/// second order and the boundary was first, so the whole thing converged at first order:
/// 5.4% at 23 cells, 2.8% at 45, 1.4% at 89, halving rather than quartering. `room_modes`
/// measures that convergence, and it is now second order.
///
/// A cross-check that the factor is the one the physics wants: with it, the acceleration the
/// scheme produces at a wall node matches the mirrored five-point Laplacian,
/// `(p₁ − 2p₀ + p₁)/dx² = 2(p₁ − p₀)/dx²`, which is the standard second-order treatment of a
/// zero-gradient boundary. Without it the scheme produced exactly half of that, and disagreed
/// with the Laplacian [`ScalarField::laplacian`] reports for the same field.
fn wall_weight(k: usize, n: usize) -> f64 {
    if k == 0 || k + 1 == n {
        2.0
    } else {
        1.0
    }
}

/// Index of the node below a position in grid units, and the fraction beyond it.
fn clamp_index(q: f64, n: usize) -> (usize, f64) {
    if q.is_nan() || q <= 0.0 {
        return (0, 0.0);
    }
    let last = (n - 1) as f64;
    if q >= last {
        return (n - 1, 0.0);
    }
    let i = q.floor();
    (i as usize, q - i)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square(cells: usize) -> Room {
        Room::of_air("room", Length::m(4.0), Length::m(4.0), cells)
    }

    fn at(x: f64, y: f64) -> LengthVec {
        LengthVec::m(x, y, 0.0)
    }

    /// The field reads the room where the room is, interpolates between nodes, and holds
    /// its value outside the walls.
    #[test]
    fn the_field_samples_the_room_and_stops_at_its_walls() {
        let room = square(33).released_in_mode(1, 0, Pressure::from_si(1.0));
        let (lx, _) = (room.width().to_si(), room.height().to_si());

        // The mode is cos(pi x / Lx), so it is +1 at one wall and -1 at the other.
        assert!((room.at(at(0.0, 2.0), Time::ZERO) - 1.0).abs() < 1e-12);
        assert!((room.at(at(lx, 2.0), Time::ZERO) + 1.0).abs() < 1e-12);
        assert!(room.at(at(lx / 2.0, 2.0), Time::ZERO).abs() < 1e-12);

        // Between two nodes, bilinear rather than nearest — a visualiser sampling finer
        // than the grid must not get stairs.
        let dx = lx / 32.0;
        let midpoint = room.at(at(dx * 0.5, 2.0), Time::ZERO);
        let ends = (room.at(at(0.0, 2.0), Time::ZERO) + room.at(at(dx, 2.0), Time::ZERO)) / 2.0;
        assert!((midpoint - ends).abs() < 1e-12, "{midpoint} against {ends}");

        // Outside the walls the value is held rather than extrapolated or indexed out of
        // range, and a NaN reads as the corner rather than casting to nonsense.
        assert!((room.at(at(-9.0, 2.0), Time::ZERO) - 1.0).abs() < 1e-12);
        assert!((room.at(at(lx + 9.0, 2.0), Time::ZERO) + 1.0).abs() < 1e-12);
        assert!(room.at(at(f64::NAN, f64::NAN), Time::ZERO).is_finite());

        // And it is uniform in depth, which is what a two-dimensional model claims.
        assert_eq!(
            room.at(LengthVec::m(1.0, 2.0, 0.0), Time::ZERO),
            room.at(LengthVec::m(1.0, 2.0, 55.0), Time::ZERO)
        );
    }

    /// **Where this differs from the bar, and why neither is wrong.** A rigid wall's
    /// boundary condition is `∂p/∂n = 0`, and on a node-centred grid a node sits *on* the
    /// wall — so the reported gradient there is zero to the last bit.
    ///
    /// `Bar1D` reports a nonzero gradient at its insulated end for the opposite reason: its
    /// grid is cell-centred, so the nearest sample is half a cell inside, where the field is
    /// still changing. The physics agrees in both cases; the grids sample it differently,
    /// and only a second implementor could have shown that.
    #[test]
    fn the_gradient_vanishes_on_a_rigid_wall() {
        let room = square(33).released_in_mode(1, 1, Pressure::from_si(1.0));
        let (lx, ly) = (room.width().to_si(), room.height().to_si());

        for (x, y) in [(0.0, 1.5), (lx, 2.5), (1.5, 0.0), (2.5, ly)] {
            let g = room.gradient(at(x, y), Time::ZERO, Length::from_si(room.dx));
            // Only the component normal to that wall has to vanish.
            let normal = if x == 0.0 || x == lx { g.x } else { g.y };
            assert!(
                normal == 0.0,
                "at ({x}, {y}) the normal gradient was {normal}"
            );
            assert_eq!(g.z, 0.0, "a two-dimensional room has no gradient in z");
        }

        // Away from the walls it is emphatically not zero, so the test above is not passing
        // on a field that is flat everywhere.
        let inside = room.gradient(at(1.0, 1.0), Time::ZERO, Length::from_si(room.dx));
        assert!(inside.length() > 0.1, "got {inside}");
    }

    /// The Laplacian against the closed form the mode shape provides.
    ///
    /// For `p = cos(nπx/Lx)cos(mπy/Ly)` the exact Laplacian is `−[(nπ/Lx)² + (mπ/Ly)²]p`,
    /// so the ratio of the two is a constant the grid should reproduce. It reproduces it to
    /// second order in `dx`, which is checked by refining rather than by asserting a
    /// tolerance that happened to pass.
    #[test]
    fn the_laplacian_matches_the_mode_it_was_given() {
        let error_at = |cells: usize| {
            let room = square(cells).released_in_mode(1, 1, Pressure::from_si(1.0));
            let (lx, ly) = (room.width().to_si(), room.height().to_si());
            let k2 = (std::f64::consts::PI / lx).powi(2) + (std::f64::consts::PI / ly).powi(2);
            // A quarter of the way in, where the mode is neither at a peak nor a node.
            let p = at(lx * 0.25, ly * 0.25);
            let expected = -k2 * room.at(p, Time::ZERO);
            let got = room.laplacian(p, Time::ZERO, Length::from_si(room.dx));
            (got / expected - 1.0).abs()
        };
        let coarse = error_at(17);
        let fine = error_at(33);
        assert!(coarse < 0.02, "17 cells was already off by {coarse}");
        // Halving dx should quarter the error. Allowing 3x rather than 4x leaves room for
        // the height being requantised to a whole number of cells.
        assert!(
            fine < coarse / 3.0,
            "refining must converge second order: {coarse} then {fine}"
        );
    }

    /// **The check that makes the field worth having, and it is not the Laplacian.**
    ///
    /// A wave equation is second order in time, so `α∇²p` is not the rate — it is the
    /// acceleration. The rate comes from the velocity divergence the scheme stores, and
    /// because it is the same expression on the same stencil it is exact.
    ///
    /// Exact about the step just *taken*, though, where the bar's was exact about the one
    /// coming. A leapfrog holds its velocities half a step behind its pressures, so the
    /// velocity for the next update does not exist yet. That is the sign convention this
    /// test exists to pin down, and getting it backwards is invisible in a picture.
    #[test]
    fn the_reported_rate_is_exactly_the_step_just_taken() {
        let mut room = square(21).released_in_mode(2, 1, Pressure::from_si(1.0));
        let dt = Time::from_si(room.dx / (343.0 * std::f64::consts::SQRT_2) * 0.9);

        let (nx, ny) = room.cells();
        let nodes: Vec<LengthVec> = (0..ny)
            .flat_map(|j| (0..nx).map(move |i| (i, j)))
            .map(|(i, j)| LengthVec::from_si(DVec3::new(i as f64, j as f64, 0.0) * room.dx))
            .collect();

        // Twice, so this cannot be an accident of starting from rest.
        for round in 0..2 {
            let before: Vec<f64> = room.pressure.clone();
            room.step(Time::ZERO, dt, &mut Exchange::new()).unwrap();
            let reported: Vec<f64> = nodes
                .iter()
                .map(|p| room.rate(*p, Time::ZERO, dt))
                .collect();

            let scale = reported.iter().fold(0.0f64, |a, v| a.max(v.abs()));
            assert!(
                scale > 1.0,
                "round {round}: nothing moved, so nothing was tested"
            );
            for (k, p) in nodes.iter().enumerate() {
                let observed = (room.pressure[k] - before[k]) / dt.to_si();
                let _ = p;
                assert!(
                    (observed - reported[k]).abs() < scale * 1e-12,
                    "round {round}, node {k}: reported {} but the step did {observed}",
                    reported[k]
                );
            }
        }

        // And the Laplacian is not that rate but the *acceleration*, which has a closed form
        // of its own: `c²∇²p = ∂²p/∂t² = −ω²p` for a mode. So dividing the field's curvature
        // by its value has to give back the mode frequency that `mode_frequency` computes
        // from the room's dimensions alone — two numbers that share no code arriving at the
        // same answer.
        //
        // Comparing the Laplacian against the rate directly would have been meaningless:
        // one is Pa/m² and the other Pa/s, so their ratio says more about the choice of
        // units than about the physics.
        let probe = at(room.width().to_si() * 0.15, room.height().to_si() * 0.3);
        let value = room.at(probe, Time::ZERO);
        let curvature = room.laplacian(probe, Time::ZERO, Length::from_si(room.dx));
        assert!(
            value.abs() > 0.1,
            "pick a probe where the mode is alive: {value}"
        );
        let omega = -room.speed * room.speed * curvature / value;
        let expected = (std::f64::consts::TAU * room.mode_frequency(2, 1).to_si()).powi(2);
        assert!(
            (omega / expected - 1.0).abs() < 0.02,
            "omega squared {omega} against the closed form {expected}"
        );
    }

    /// **What the wall weighting is, stated exactly.** One step from rest moves every
    /// pressure by `½h²c²∇²p`, walls included, where `∇²` is the mirrored five-point stencil
    /// [`ScalarField::laplacian`] reports.
    ///
    /// **The ½ is Taylor's.** From rest `ṗ(0) = −ρc²∇·v(0) = 0`, so `p(h) − p(0) =
    /// ½h²p̈(0) = ½h²c²∇²p`. This test used to assert `h²c²∇²p` with no half, and passed,
    /// because the scheme kicked the velocity a whole step on its first update where the
    /// initial condition only entitled it to half — the startup error that made a
    /// second-order scheme converge at first order. The test had turned the bug into the
    /// specification, which is what a test written from the implementation does.
    ///
    /// After the first half-step the velocities are `−h∇p/2ρdx`, so the pressure change is
    /// `½h²c²` times whatever Laplacian the scheme implies — which makes this a direct read
    /// of that operator, boundary treatment and all. Without the wall weighting a boundary
    /// node moves by exactly half as much again, and the two would disagree there and only
    /// there.
    ///
    /// Cheaper and sharper than measuring a frequency: it is exact, it needs one step, and it
    /// names the wrong node rather than reporting that a number came out low.
    #[test]
    fn one_step_from_rest_is_the_laplacian_the_field_reports() {
        let mut room = square(17).released_in_mode(2, 1, Pressure::from_si(1.0));
        let h = room.max_stable_dt(Time::ZERO).to_si() * 0.7;
        let before = room.pressure.clone();
        // Read the operator *before* stepping: the step is what is being measured, so
        // comparing it against the Laplacian of the field it produced would be circular and
        // off by a few percent besides.
        let (nx, ny) = room.cells();
        let reported: Vec<f64> = (0..ny)
            .flat_map(|j| (0..nx).map(move |i| (i, j)))
            .map(|(i, j)| {
                let p = LengthVec::from_si(DVec3::new(i as f64, j as f64, 0.0) * room.dx);
                room.laplacian(p, Time::ZERO, Length::from_si(room.dx))
            })
            .collect();
        room.step(Time::ZERO, Time::from_si(h), &mut Exchange::new())
            .unwrap();

        let mut walls_checked = 0;
        let scale = before.iter().fold(0.0f64, |a, v| a.max(v.abs())) / (room.dx * room.dx);
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let implied =
                    (room.pressure[k] - before[k]) / (0.5 * h * h * room.speed * room.speed);
                assert!(
                    (implied - reported[k]).abs() < scale * 1e-12,
                    "node ({i},{j}): the step implies {implied} but the field says {}",
                    reported[k]
                );
                if i == 0 || j == 0 || i == nx - 1 || j == ny - 1 {
                    walls_checked += 1;
                }
            }
        }
        assert_eq!(
            walls_checked,
            2 * (nx + ny) - 4,
            "every wall node was covered"
        );
    }

    /// A pressure field swings either side of zero, which a temperature never does. Worth
    /// pinning because it is a requirement on anything that draws one.
    #[test]
    fn the_field_is_signed_about_zero() {
        let room = square(33).released_in_mode(1, 1, Pressure::from_si(2.0));
        let (lx, ly) = (room.width().to_si(), room.height().to_si());
        let corners = [
            room.at(at(0.0, 0.0), Time::ZERO),
            room.at(at(lx, 0.0), Time::ZERO),
            room.at(at(0.0, ly), Time::ZERO),
            room.at(at(lx, ly), Time::ZERO),
        ];
        assert!(corners.iter().any(|v| *v > 1.9), "{corners:?}");
        assert!(corners.iter().any(|v| *v < -1.9), "{corners:?}");
        // And it averages to nothing, so a colour map anchored at the minimum would put the
        // midpoint in the wrong place.
        let mean: f64 = room.pressure.iter().sum::<f64>() / room.pressure.len() as f64;
        assert!(mean.abs() < 1e-12, "mean {mean}");
    }

    /// The closed form, and the two things about it a tube cannot show.
    #[test]
    fn room_modes_are_not_a_harmonic_series() {
        let room = Room::of_air("room", Length::m(5.0), Length::m(3.0), 51);
        assert!((room.width().to_si() - 5.0).abs() < 1e-12);
        // The height is quantised to whole cells, so it lands near 3 m rather than on it.
        assert!(
            (room.height().to_si() - 3.0).abs() < 0.11,
            "{}",
            room.height().to_si()
        );

        // (c/2) sqrt((1/5)^2) = 34.3 Hz along the length.
        assert!((room.mode_frequency(1, 0).to_si() - 34.3).abs() < 0.1);
        // And (c/2)(1/3) = 57.2 Hz across the width.
        assert!(
            (room.mode_frequency(0, 1).to_si() - 57.2).abs() < 1.5,
            "{}",
            room.mode_frequency(0, 1).to_si()
        );

        // The oblique mode is the quadrature sum, not the arithmetic one — which is what
        // makes the series inharmonic.
        let f10 = room.mode_frequency(1, 0).to_si();
        let f01 = room.mode_frequency(0, 1).to_si();
        let f11 = room.mode_frequency(1, 1).to_si();
        assert!(
            (f11 - (f10 * f10 + f01 * f01).sqrt()).abs() < 1e-9,
            "the oblique mode is the quadrature sum"
        );
        assert!(f11 < f10 + f01, "and so it is below the arithmetic sum");

        // In a square room the diagonal mode is exactly root two times the axial one: a
        // tritone, which is why a square room sounds wrong rather than merely loud.
        let sq = square(41);
        assert!((sq.mode_frequency(1, 1) / sq.mode_frequency(1, 0) - 2f64.sqrt()).abs() < 1e-12);
        assert_eq!(sq.mode_frequency(0, 0).to_si(), 0.0);
    }

    /// The modes crowd together as they rise, going as `f²` in two dimensions. That is
    /// why a small room's colouration is heard at the bottom of the range and turns into
    /// an even hiss further up.
    #[test]
    fn the_modes_get_denser_as_the_frequency_rises() {
        let room = Room::of_air("room", Length::m(5.0), Length::m(4.0), 51);
        let count = |f: f64| room.modes_below(Frequency::hz(f)).len();

        // Doubling the frequency should roughly quadruple the count.
        let (low, high) = (count(200.0), count(400.0));
        let ratio = high as f64 / low as f64;
        assert!(
            (ratio - 4.0).abs() < 1.0,
            "the count should go as f squared: {low} below 200 Hz and {high} below 400, \
             a ratio of {ratio:.2}"
        );

        // They come out in ascending order and the first is the longest dimension's.
        let modes = room.modes_below(Frequency::hz(120.0));
        assert!(modes.len() > 4);
        for pair in modes.windows(2) {
            assert!(pair[0].2 <= pair[1].2, "modes should be sorted");
        }
        assert_eq!((modes[0].0, modes[0].1), (1, 0));
        // Nothing below the fundamental, and the zero mode is not a mode.
        assert!(room.modes_below(Frequency::hz(30.0)).is_empty());
    }

    /// The CFL condition tightens by exactly the diagonal, and going past it is refused.
    #[test]
    fn two_dimensions_tighten_the_courant_limit_by_root_two() {
        let room = square(41);
        // 4 m over 40 intervals is 100 mm; 343 m/s crosses that in 292 us, and the
        // two-dimensional limit is that over root two.
        let limit = room.max_stable_dt(Time::ZERO);
        assert!(
            (limit.in_us() - 206.2).abs() < 0.5,
            "limit {} us",
            limit.in_us()
        );
        let tube_limit = 0.1 / 343.0;
        assert!(
            (tube_limit / limit.to_si() - 2f64.sqrt()).abs() < 1e-9,
            "exactly root two tighter than one dimension"
        );
        assert!((room.courant(limit) - 1.0).abs() < 1e-12);

        let mut room = square(41);
        let mut bus = Exchange::new();
        assert!(room.step(Time::ZERO, limit, &mut bus).is_ok());
        let err = room
            .step(Time::ZERO, limit * 1.01, &mut bus)
            .expect_err("past the limit must be refused");
        assert_eq!(err.quantity, "Courant number");
    }

    /// A mode excited exactly should oscillate at exactly its predicted frequency: half a
    /// period inverts it, a whole one restores it.
    ///
    /// This is what ties the closed form to the code actually stepping, and it is the test
    /// that would fail if the stencil, the staggering or the wall condition were wrong.
    #[test]
    fn an_excited_mode_oscillates_at_its_predicted_frequency() {
        for (nx, ny) in [(1u32, 0u32), (1, 1), (2, 1)] {
            let mut room = square(81).released_in_mode(nx, ny, Pressure::from_si(100.0));
            let f = room.mode_frequency(nx, ny).to_si();
            let period = 1.0 / f;
            // Well inside the limit, so the numerical dispersion stays small.
            let dt = room.max_stable_dt(Time::ZERO) * 0.5;
            let steps = (period / dt.to_si()).round() as u32;
            let start = room.pressure_at(0, 0).to_si();
            assert!(
                start.abs() > 50.0,
                "the corner is an antinode of every mode"
            );

            let mut bus = Exchange::new();
            for _ in 0..steps / 2 {
                room.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            let half = room.pressure_at(0, 0).to_si();
            assert!(
                half < -0.85 * start,
                "mode ({nx},{ny}): half a period should invert it, {start:.1} to {half:.1}"
            );

            for _ in 0..steps - steps / 2 {
                room.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            let whole = room.pressure_at(0, 0).to_si();
            assert!(
                (whole / start - 1.0).abs() < 0.1,
                "mode ({nx},{ny}): a whole period should restore it, {start:.1} to \
                 {whole:.1}"
            );
        }
    }

    /// Rigid walls keep the energy in, so it holds rather than draining.
    #[test]
    fn a_closed_room_keeps_its_energy() {
        let mut room = square(61).released_from(|x, y| {
            let (dx, dy) = (x.to_si() - 2.0, y.to_si() - 2.0);
            let r2 = dx * dx + dy * dy;
            Pressure::from_si(200.0 * (-r2 / 0.09).exp())
        });
        let start = room.energy().to_si();
        assert!(start > 0.0);

        let dt = room.max_stable_dt(Time::ZERO) * 0.5;
        let mut bus = Exchange::new();
        for _ in 0..3000 {
            room.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        let now = room.energy().to_si();
        assert!(
            (now / start - 1.0).abs() < 0.02,
            "energy went from {start:e} to {now:e}"
        );
        // Rigid walls publish nothing.
        assert!(bus.peek(quantity::ENERGY).abs() < 1e-30);
    }

    /// A pulse in the middle spreads outwards in a circle rather than along the axes,
    /// which is the thing a pair of one-dimensional tubes could not produce.
    #[test]
    fn a_pulse_spreads_isotropically() {
        let mut room = square(81).released_from(|x, y| {
            let (dx, dy) = (x.to_si() - 2.0, y.to_si() - 2.0);
            let r2 = dx * dx + dy * dy;
            Pressure::from_si(200.0 * (-r2 / 0.02).exp())
        });
        let dt = room.max_stable_dt(Time::ZERO) * 0.5;
        let mut bus = Exchange::new();
        // Long enough for the front to travel 0.8 m, well short of the walls.
        let steps = (0.8 / 343.0 / dt.to_si()).round() as u32;
        for _ in 0..steps {
            room.step(Time::ZERO, dt, &mut bus).unwrap();
        }

        // Sample the same distance from the centre along an axis and along the diagonal.
        let (cx, cy) = (40usize, 40usize);
        let offset = 16usize; // 0.8 m at 50 mm a cell
        let diag = (offset as f64 / 2f64.sqrt()).round() as usize;
        let along_axis = room.pressure_at(cx + offset, cy).to_si().abs();
        let along_diagonal = room.pressure_at(cx + diag, cy + diag).to_si().abs();
        assert!(
            along_axis > 1.0,
            "the front should have arrived: {along_axis}"
        );
        assert!(
            (along_axis - along_diagonal).abs() / along_axis.max(along_diagonal) < 0.25,
            "a circular front should reach both at once: {along_axis:.3} along the axis \
             and {along_diagonal:.3} along the diagonal"
        );
    }

    /// Under the scheduler, subcycling to the two-dimensional limit.
    #[test]
    fn the_domain_runs_under_the_scheduler() {
        use pantometry_core::{Schedule, Simulation};

        let room = square(41).released_from(|x, y| {
            let (dx, dy) = (x.to_si() - 2.0, y.to_si() - 2.0);
            Pressure::from_si(200.0 * (-(dx * dx + dy * dy) / 0.5).exp())
        });
        // Tight, because `Room::energy` reports the time-centred quantity the leapfrog
        // conserves exactly rather than the one that wobbles by a few percent. The first
        // version of that function needed 5e-2 here and was still failing.
        let mut sim = Simulation::new(Schedule::Multirate)
            .conservation_tolerance(1e-9)
            .with(room);
        let report = sim.advance(Time::ms(1.0)).unwrap();
        assert_eq!(report.substeps[0].0, "room");
        // 1 ms at a 206 us limit is 5 substeps.
        assert_eq!(report.substeps[0].1, 5);
        for _ in 0..20 {
            sim.advance(Time::ms(1.0)).expect("energy should hold");
        }
    }

    /// The leapfrog is started correctly, and the proof is the rate.
    ///
    /// A mode released from rest follows `|cos(2 pi f t)|`; the departure from it must
    /// *quarter* on refinement, not halve. It halved until the velocity was given its missing
    /// half step, and this pins the fix — because losing it again produces a scheme that
    /// still runs, still conserves energy to 1e-15, and is quietly first order. That
    /// combination is why it survived: nothing here was checking a rate, and everything that
    /// was being checked passed.
    ///
    /// Measured, at a fixed 20 ms rather than a fixed step count: 7.83e-4 at 31 cells,
    /// 1.76e-4 at 61, 3.83e-5 at 121, 9.95e-6 at 241 — ratios of 4.5, 4.6, 3.9.
    ///
    /// **Two ways this test could have measured nothing**, both met while writing it. A fixed
    /// number of steps covers less physical time as the grid refines, which flattered the
    /// ratio to 8.7 by comparing eight periods against one. And a *square* room stepped at
    /// exactly its stability limit solves its symmetric modes to 1e-13 — the two-dimensional
    /// magic time step — so the room here is 4.4 by 3.1 and not 4 by 4.
    #[test]
    fn starting_the_leapfrog_correctly_makes_the_scheme_second_order() {
        let worst = |cells: usize| {
            let mut room = Room::of_air("room", Length::m(4.4), Length::m(3.1), cells)
                .released_in_mode(1, 1, Pressure::from_si(1.0));
            let f = room.mode_frequency(1, 1).to_si();
            let steps = (0.02 / room.max_stable_dt(Time::ZERO).to_si()).ceil() as usize;
            let h = Time::from_si(0.02 / steps as f64);
            let mut bus = Exchange::new();
            let (mut t, mut worst) = (0.0f64, 0.0f64);
            for _ in 0..steps {
                room.step(Time::from_si(t), h, &mut bus).unwrap();
                t += h.to_si();
                let want = (2.0 * std::f64::consts::PI * f * t).cos().abs();
                worst = worst.max((room.peak_pressure().to_si() - want).abs());
            }
            worst
        };
        let (coarse, fine) = (worst(31), worst(241));
        let fall = coarse / fine;
        // Three doublings: second order is 64x, first order 8x. Measured 78.7.
        assert!(
            fall > 30.0,
            "31 -> 241 cells fell only {fall:.1}x ({coarse:.4e} -> {fine:.4e});              a first-order startup would give about 8"
        );
    }

    /// The startup adjustment is real, small, and second order.
    ///
    /// It is the one number this domain reports that its state does not directly say, so it
    /// is worth pinning rather than trusting. Converting the released velocity from `t = 0`
    /// to the half step moves the invariant by `O(h²)`: 0.42% of the total at 31 cells,
    /// 0.10% at 61, 0.026% at 121.
    #[test]
    fn the_startup_adjustment_is_second_order_and_declared() {
        let relative = |cells: usize| {
            let mut room = square(cells).released_in_mode(1, 1, Pressure::from_si(1.0));
            let total = room.energy().to_si();
            room.step(
                Time::ZERO,
                room.max_stable_dt(Time::ZERO),
                &mut Exchange::new(),
            )
            .unwrap();
            room.startup_adjustment().to_si().abs() / total
        };
        let (a, b) = (relative(31), relative(121));
        assert!(
            a > 1e-4,
            "the adjustment should not be nothing, got {a:.2e}"
        );
        // Two doublings of the grid, so sixteenfold if it is second order.
        let fall = a / b;
        assert!(
            (10.0..24.0).contains(&fall),
            "31 -> 121 cells should quarter it twice: {a:.3e} -> {b:.3e}, a factor of {fall:.1}"
        );
    }

    /// Degenerate rooms are handled.
    #[test]
    fn degenerate_rooms_are_handled() {
        let tiny = Room::of_air("tiny", Length::m(0.1), Length::m(0.1), 1);
        let (nx, ny) = tiny.cells();
        assert!(nx >= 3 && ny >= 3);

        let mut silent = Room::new(
            "silent",
            Length::m(1.0),
            Length::m(1.0),
            11,
            Length::m(1.0),
            Density::kg_per_m3(1.0),
            Velocity::m_per_s(0.0),
        );
        assert!(!silent.max_stable_dt(Time::ZERO).to_si().is_finite());
        let mut bus = Exchange::new();
        assert!(silent.step(Time::ZERO, Time::s(1.0), &mut bus).is_ok());
        assert_eq!(silent.peak_pressure().to_si(), 0.0);

        // A zero step changes nothing.
        let mut room = square(21).released_from(|_, _| Pressure::from_si(7.0));
        let before = room.pressure_at(5, 5);
        room.step(Time::ZERO, Time::ZERO, &mut bus).unwrap();
        assert_eq!(room.pressure_at(5, 5), before);
    }
}
