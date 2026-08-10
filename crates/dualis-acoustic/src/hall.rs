//! The wave equation with a ceiling.
//!
//! [`Room`](crate::Room) is two-dimensional, and says so in four places, because a third
//! dimension is genuinely expensive: `√3` in the Courant limit as well as the obvious factor in
//! cells. That reduction is right for a floor plan and wrong for the two things a room's acoustics
//! are usually *about* — the vertical modes between floor and ceiling, and the oblique modes that
//! need all three axes at once.
//!
//! A 4.4 × 3.1 m room has its first axial mode at 39 Hz and its first vertical one at 71 Hz for a
//! 2.4 m ceiling. A two-dimensional model does not have the second one at all. It is not less
//! accurate about it; it does not have it.
//!
//! # What the third dimension costs, exactly
//!
//! ```text
//!   1D   dt ≤ dx/c            Tube
//!   2D   dt ≤ dx/(c√2)        Room
//!   3D   dt ≤ dx/(c√3)        Hall
//! ```
//!
//! A wave crossing a cell diagonally covers `√d` cells while the stencil sees one. So the same
//! room at the same spacing is `n` times the cells **and** 22% more steps than the 2D version —
//! and a room resolved to 61 cells across is 61 × 43 = 2623 nodes in two dimensions and 61 × 43 ×
//! 34 = 89,182 in three.
//!
//! The mode count grows the same way. Below a frequency `f`, the number of modes goes as `f²` in
//! two dimensions and as `f³` in three, which is *why* a real room's resonances merge into a hiss
//! above the Schroeder frequency while a two-dimensional model keeps them separable much further
//! up. That is the qualitative thing the reduction changes, not merely a number.
//!
//! # What it is checked against
//!
//! The rigid-wall mode frequency, which is exact:
//!
//! ```text
//!   f(a,b,c) = (c_sound/2)·√( (a/Lx)² + (b/Ly)² + (c/Lz)² )
//! ```
//!
//! A mode released at its own shape stays that shape and every point rides `cos(2πft)`, so the
//! peak of the field follows `|cos(2πft)|` — a closed form the integration never sees. The check
//! is on the **rate** at which the integration approaches it, because a scheme that is first order
//! where it claims second is the defect this crate has already shipped once: `Room` and `Tube`
//! both started their leapfrog with the velocity at the wrong time level, `O(h)` and permanent,
//! and it survived 345 passing tests.
//!
//! This domain carries the fix from its first line — see [`Hall::released_from`] — rather than
//! inheriting the bug and being corrected later.

use dualis_core::conserved::quantity;
use dualis_core::{Domain, Exchange, Kind, Ledger, Reading, ScalarField, Violation};
use dualis_units::{Density, Energy, Frequency, Length, LengthVec, Pressure, Time, Velocity};
use glam::DVec3;

/// A rectangular room with rigid walls, floor and ceiling, on a node-centred staggered grid.
///
/// Pressure sits on the nodes; velocity sits on the faces between them, one array per axis. The
/// walls are perfectly rigid, so no velocity passes through them and the domain conserves energy
/// exactly — which is what makes the audit worth running rather than a formality.
///
/// Cells are **cubes**. The height and depth are quantised to a whole number of the spacing the
/// width sets, so the Courant limit stays isotropic; a 3.1 m ceiling on a 4.4 m room at 45 nodes
/// across becomes 3.1 m to the nearest 0.1 m. [`Hall::height`] and [`Hall::depth`] report what it
/// actually is, and [`Hall::mode_frequency`] uses those rather than what was asked for.
#[derive(Clone, Debug)]
pub struct Hall {
    name: String,
    /// Pressure at nodes, indexed `i + nx*(j + ny*k)`.
    pressure: Vec<f64>,
    /// Pressure one whole step earlier, so [`Hall::energy`] can report the quantity the scheme
    /// actually conserves rather than one that wobbles at the step frequency.
    pressure_prev: Vec<f64>,
    /// Velocity on the faces between neighbours along x: `(nx - 1) * ny * nz`.
    vx: Vec<f64>,
    /// Along y: `nx * (ny - 1) * nz`.
    vy: Vec<f64>,
    /// Along z: `nx * ny * (nz - 1)`.
    vz: Vec<f64>,
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    speed: f64,
    density: f64,
    /// Whether the velocities have been offset to the half step the scheme carries them at.
    velocity_staggered: bool,
    /// The released state's physical energy minus the scheme's invariant once staggered.
    energy_offset: f64,
    energy_datum: f64,
    saved: Option<Snapshot>,
}

type Snapshot = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, bool);

impl Hall {
    /// A room of the given size, at rest, with rigid surfaces.
    ///
    /// `across` is the node count along the width and sets the spacing; height and depth are
    /// quantised to whole cells of that spacing so the cells stay cubic.
    pub fn new(
        name: impl Into<String>,
        width: Length,
        height: Length,
        depth: Length,
        across: usize,
        speed: Velocity,
        density: Density,
    ) -> Hall {
        let nx = across.max(2);
        let dx = (width.to_si() / (nx - 1) as f64).max(f64::MIN_POSITIVE);
        // Cells, then nodes. A zero extent is **one** node and no faces along that axis, which
        // is what makes a flat hall reduce exactly to a room rather than to a two-node slab with
        // a spurious mode across it.
        let quantise = |l: Length| (l.to_si() / dx).round().max(0.0) as usize + 1;
        let (ny, nz) = (quantise(height), quantise(depth));
        let nodes = nx * ny * nz;
        Hall {
            name: name.into(),
            pressure: vec![0.0; nodes],
            pressure_prev: vec![0.0; nodes],
            vx: vec![0.0; (nx - 1) * ny * nz],
            vy: vec![0.0; nx * (ny - 1) * nz],
            vz: vec![0.0; nx * ny * (nz - 1)],
            nx,
            ny,
            nz,
            dx,
            speed: speed.to_si(),
            density: density.to_si(),
            velocity_staggered: false,
            energy_offset: 0.0,
            energy_datum: 0.0,
            saved: None,
        }
    }

    /// A room of air at 20 °C: 343 m/s, 1.204 kg/m³.
    pub fn of_air(
        name: impl Into<String>,
        width: Length,
        height: Length,
        depth: Length,
        across: usize,
    ) -> Hall {
        Hall::new(
            name,
            width,
            height,
            depth,
            across,
            Velocity::from_si(343.0),
            Density::from_si(1.204),
        )
    }

    /// Set the pressure everywhere from a profile of position, and release it from rest.
    ///
    /// **The velocity is left at `t = 0`, where an initial condition puts it, and the scheme is
    /// told so.** A staggered leapfrog carries velocity at `t = −h/2`, so the very first velocity
    /// update has only half a step to travel; taking a whole one leaves an `O(h)` error that
    /// never decays, and since `h` follows `dx` through the Courant condition it drags the whole
    /// second-order scheme to first order.
    ///
    /// That is not hypothetical here. `Room` and `Tube` both shipped with it, it survived 345
    /// passing tests, and two of those tests had turned the bug into the specification. This
    /// domain starts correct.
    pub fn released_from(mut self, profile: impl Fn(LengthVec) -> Pressure) -> Hall {
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let p = LengthVec::from_si(DVec3::new(i as f64, j as f64, k as f64) * self.dx);
                    let idx = self.at(i, j, k);
                    self.pressure[idx] = profile(p).to_si();
                }
            }
        }
        self.pressure_prev.clone_from(&self.pressure);
        self.vx.iter_mut().for_each(|v| *v = 0.0);
        self.vy.iter_mut().for_each(|v| *v = 0.0);
        self.vz.iter_mut().for_each(|v| *v = 0.0);
        self.velocity_staggered = false;
        self.energy_offset = 0.0;
        self.energy_datum = self.invariant_si();
        self
    }

    /// Excite one rigid-wall mode exactly: `cos(aπx/Lx)·cos(bπy/Ly)·cos(cπz/Lz)`.
    ///
    /// The shape whose frequency [`mode_frequency`](Hall::mode_frequency) predicts in closed
    /// form, so the prediction can be checked against the integration rather than against itself.
    pub fn released_in_mode(self, mode: (u32, u32, u32), amplitude: Pressure) -> Hall {
        let (lx, ly, lz) = (
            self.width().to_si(),
            self.height().to_si(),
            self.depth().to_si(),
        );
        let a = amplitude.to_si();
        self.released_from(|p| {
            let v = p.to_si();
            let phase = |n: u32, x: f64, l: f64| {
                if l > 0.0 {
                    (n as f64 * std::f64::consts::PI * x / l).cos()
                } else {
                    1.0
                }
            };
            Pressure::from_si(
                a * phase(mode.0, v.x, lx) * phase(mode.1, v.y, ly) * phase(mode.2, v.z, lz),
            )
        })
    }

    /// Node counts along x, y and z.
    pub fn nodes(&self) -> (usize, usize, usize) {
        (self.nx, self.ny, self.nz)
    }

    /// The node spacing.
    pub fn spacing(&self) -> Length {
        Length::from_si(self.dx)
    }

    /// The width, which is what was asked for.
    pub fn width(&self) -> Length {
        Length::from_si((self.nx - 1) as f64 * self.dx)
    }

    /// The height, **quantised** to a whole number of cells.
    pub fn height(&self) -> Length {
        Length::from_si((self.ny - 1) as f64 * self.dx)
    }

    /// The depth, quantised the same way.
    pub fn depth(&self) -> Length {
        Length::from_si((self.nz - 1) as f64 * self.dx)
    }

    /// The pressure at one node.
    pub fn pressure_at(&self, i: usize, j: usize, k: usize) -> Pressure {
        let idx = self.at(i.min(self.nx - 1), j.min(self.ny - 1), k.min(self.nz - 1));
        Pressure::from_si(self.pressure[idx])
    }

    /// The largest pressure anywhere, which is what a mode's oscillation shows.
    pub fn peak_pressure(&self) -> Pressure {
        Pressure::from_si(self.pressure.iter().fold(0.0f64, |m, v| m.max(v.abs())))
    }

    /// The energy the scheme conserves, measured from the released state.
    pub fn energy(&self) -> Energy {
        Energy::from_si(self.invariant_si() + self.energy_offset)
    }

    /// How much the discrete invariant differs from the released state's physical energy.
    ///
    /// Zero until the first step, because the invariant is a function of the step size and there
    /// is no step size before then. `O(h²)`, and therefore converging — the same quantity
    /// `Room::startup_adjustment` reports, for the same reason.
    pub fn startup_adjustment(&self) -> Energy {
        Energy::from_si(self.energy_offset)
    }

    /// The rigid-wall mode frequency, in closed form.
    ///
    /// `(c/2)·√((a/Lx)² + (b/Ly)² + (c/Lz)²)`, using the **quantised** dimensions rather than
    /// the ones asked for — otherwise the closed form would be about a room the grid is not.
    pub fn mode_frequency(&self, mode: (u32, u32, u32)) -> Frequency {
        let lengths = [
            self.width().to_si(),
            self.height().to_si(),
            self.depth().to_si(),
        ];
        let counts = [mode.0 as f64, mode.1 as f64, mode.2 as f64];
        let sum: f64 = counts
            .iter()
            .zip(lengths)
            .map(|(n, l)| if l > 0.0 { (n / l).powi(2) } else { 0.0 })
            .sum();
        Frequency::from_si(self.speed / 2.0 * sum.sqrt())
    }

    /// Every mode below a frequency, ascending.
    ///
    /// **The count goes as `f³`**, because the modes fill a *volume* of the `(a,b,c)` lattice
    /// rather than an area. That is the qualitative thing a two-dimensional model gets wrong: it
    /// keeps resonances separable far above the frequency at which a real room's have merged.
    pub fn modes_below(&self, limit: Frequency) -> Vec<((u32, u32, u32), Frequency)> {
        let lengths = [
            self.width().to_si(),
            self.height().to_si(),
            self.depth().to_si(),
        ];
        let bound = |l: f64| {
            if l <= 0.0 || self.speed <= 0.0 {
                0
            } else {
                (2.0 * limit.to_si() * l / self.speed).floor() as u32
            }
        };
        let [mx, my, mz] = [bound(lengths[0]), bound(lengths[1]), bound(lengths[2])];
        let mut found = Vec::new();
        for a in 0..=mx {
            for b in 0..=my {
                for c in 0..=mz {
                    if a == 0 && b == 0 && c == 0 {
                        continue;
                    }
                    let f = self.mode_frequency((a, b, c));
                    if f.to_si() <= limit.to_si() {
                        found.push(((a, b, c), f));
                    }
                }
            }
        }
        found.sort_by(|x, y| x.1.to_si().total_cmp(&y.1.to_si()));
        found
    }

    /// `c·dt·√3/dx`. Must stay at or under 1.
    pub fn courant(&self, dt: Time) -> f64 {
        self.speed * dt.to_si() * 3f64.sqrt() / self.dx
    }

    fn at(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.nx * (j + self.ny * k)
    }

    /// The scheme's invariant, before the release datum is applied.
    ///
    /// Every term carries the volume it actually occupies. A node on a wall owns half a cell, one
    /// on an edge a quarter and one in a corner an eighth — the same weights the update uses,
    /// which is what makes the conservation exact rather than approximate. In two dimensions
    /// getting this wrong dropped the whole scheme to first order; in three there are three
    /// weights to get wrong instead of two.
    fn invariant_si(&self) -> f64 {
        let volume = self.dx * self.dx * self.dx;
        let rc2 = self.density * self.speed * self.speed;
        if rc2 <= 0.0 {
            return 0.0;
        }
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let mut total = 0.0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let share =
                        1.0 / (wall_weight(i, nx) * wall_weight(j, ny) * wall_weight(k, nz));
                    let idx = self.at(i, j, k);
                    total +=
                        share * self.pressure[idx] * self.pressure_prev[idx] / (2.0 * rc2) * volume;
                }
            }
        }
        // A vx face spans a whole cell along x and sits at a node in y and z, so it takes those
        // two shares and not the x one. Each axis is the same statement rotated.
        let half_rho = 0.5 * self.density;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx - 1 {
                    let v = self.vx[i + (nx - 1) * (j + ny * k)];
                    let share = 1.0 / (wall_weight(j, ny) * wall_weight(k, nz));
                    total += share * half_rho * v * v * volume;
                }
            }
        }
        for k in 0..nz {
            for j in 0..ny - 1 {
                for i in 0..nx {
                    let v = self.vy[i + nx * (j + (ny - 1) * k)];
                    let share = 1.0 / (wall_weight(i, nx) * wall_weight(k, nz));
                    total += share * half_rho * v * v * volume;
                }
            }
        }
        for k in 0..nz - 1 {
            for j in 0..ny {
                for i in 0..nx {
                    let v = self.vz[i + nx * (j + ny * k)];
                    let share = 1.0 / (wall_weight(i, nx) * wall_weight(j, ny));
                    total += share * half_rho * v * v * volume;
                }
            }
        }
        total
    }
}

impl Domain for Hall {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// `dx/(c√3)` — the diagonal of a cube rather than of a square.
    fn max_stable_dt(&self, _now: Time) -> Time {
        if self.speed <= 0.0 {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(self.dx / (self.speed * 3f64.sqrt()))
    }

    fn step(&mut self, _t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let courant = self.courant(dt);
        if courant > 1.0 + 1e-12 {
            return Err(Violation {
                quantity: "Courant number".to_string(),
                site: format!("{} (three-dimensional wave equation)", self.name),
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
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);

        // Half a step the first time; see `released_from`.
        let vh = if self.velocity_staggered { h } else { 0.5 * h };
        let starting = !self.velocity_staggered;
        self.velocity_staggered = true;

        // Faces: rho du/dt = -grad p.
        let gain = vh / (self.density * self.dx);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx - 1 {
                    let d = self.pressure[self.at(i + 1, j, k)] - self.pressure[self.at(i, j, k)];
                    self.vx[i + (nx - 1) * (j + ny * k)] -= gain * d;
                }
            }
        }
        for k in 0..nz {
            for j in 0..ny - 1 {
                for i in 0..nx {
                    let d = self.pressure[self.at(i, j + 1, k)] - self.pressure[self.at(i, j, k)];
                    self.vy[i + nx * (j + (ny - 1) * k)] -= gain * d;
                }
            }
        }
        for k in 0..nz - 1 {
            for j in 0..ny {
                for i in 0..nx {
                    let d = self.pressure[self.at(i, j, k + 1)] - self.pressure[self.at(i, j, k)];
                    self.vz[i + nx * (j + ny * k)] -= gain * d;
                }
            }
        }

        // Nodes: dp/dt = -rho c^2 div u. A rigid surface passes no velocity, so the flux outside
        // the domain is zero and no ghost nodes are needed. A node on a wall has a control volume
        // half a cell thick normal to that wall, so its divergence is divided by half the spacing
        // — which is what `wall_weight` is, per axis.
        self.pressure_prev.copy_from_slice(&self.pressure);
        let scale = h * rc2 / self.dx;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let lo_x = if i > 0 {
                        self.vx[i - 1 + (nx - 1) * (j + ny * k)]
                    } else {
                        0.0
                    };
                    let hi_x = if i + 1 < nx {
                        self.vx[i + (nx - 1) * (j + ny * k)]
                    } else {
                        0.0
                    };
                    let lo_y = if j > 0 {
                        self.vy[i + nx * (j - 1 + (ny - 1) * k)]
                    } else {
                        0.0
                    };
                    let hi_y = if j + 1 < ny {
                        self.vy[i + nx * (j + (ny - 1) * k)]
                    } else {
                        0.0
                    };
                    let lo_z = if k > 0 {
                        self.vz[i + nx * (j + ny * (k - 1))]
                    } else {
                        0.0
                    };
                    let hi_z = if k + 1 < nz {
                        self.vz[i + nx * (j + ny * k)]
                    } else {
                        0.0
                    };
                    let idx = self.at(i, j, k);
                    self.pressure[idx] -= scale
                        * (wall_weight(i, nx) * (hi_x - lo_x)
                            + wall_weight(j, ny) * (hi_y - lo_y)
                            + wall_weight(k, nz) * (hi_z - lo_z));
                }
            }
        }

        // Startup, once. The released state carries velocity at `t = 0` and the invariant is
        // defined on a staggered state, so converting between them moves the number by `O(h²)`.
        // That is discretisation and not a leak, so it is reported by `startup_adjustment` rather
        // than handed to the audit — but it is bounded, because a real first-step error (a sign,
        // a factor of two, a boundary) is not `O(h²)` and a quarter of the energy is far outside
        // anything the Courant condition permits.
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
            self.vz.clone(),
            self.velocity_staggered,
        ));
    }

    fn restore(&mut self) {
        if let Some((p, vx, vy, vz, staggered)) = self.saved.clone() {
            self.pressure_prev.copy_from_slice(&p);
            self.pressure = p;
            self.vx = vx;
            self.vy = vy;
            self.vz = vz;
            self.velocity_staggered = staggered;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    /// The peak pressure. Not the mean, which is zero by symmetry and would be a column of
    /// rounding.
    fn readings(&self) -> Vec<Reading> {
        vec![Reading::new(
            &self.name,
            "peak",
            self.peak_pressure().to_si(),
            "Pa",
        )]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }

    fn as_field(&self) -> Option<&dyn ScalarField> {
        Some(self)
    }
}

impl ScalarField for Hall {
    fn unit(&self) -> &'static str {
        "Pa"
    }

    /// Trilinear between nodes, clamped at the surfaces.
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let q = p.to_si() / self.dx;
        if q.is_nan() {
            return self.pressure[0];
        }
        let (i, fx) = clamp_index(q.x, self.nx);
        let (j, fy) = clamp_index(q.y, self.ny);
        let (k, fz) = clamp_index(q.z, self.nz);
        let (i1, j1, k1) = (
            (i + 1).min(self.nx - 1),
            (j + 1).min(self.ny - 1),
            (k + 1).min(self.nz - 1),
        );
        let g = |a: usize, b: usize, c: usize| self.pressure[self.at(a, b, c)];
        let lerp = |lo: f64, hi: f64, t: f64| lo * (1.0 - t) + hi * t;
        let z0 = lerp(
            lerp(g(i, j, k), g(i1, j, k), fx),
            lerp(g(i, j1, k), g(i1, j1, k), fx),
            fy,
        );
        let z1 = lerp(
            lerp(g(i, j, k1), g(i1, j, k1), fx),
            lerp(g(i, j1, k1), g(i1, j1, k1), fx),
            fy,
        );
        lerp(z0, z1, fz)
    }

    /// `∂p/∂t = −ρc²∇·u`, from the velocities the scheme stores.
    ///
    /// Not from the Laplacian. A wave is second order in time, so the curvature of the pressure
    /// gives its *acceleration*; the rate comes from the companion velocity field. That is the
    /// same distinction `Room` documents, and it is the reason `rate` is answerable at all
    /// without keeping history.
    fn rate(&self, p: LengthVec, _t: Time, _dt: Time) -> f64 {
        let q = p.to_si() / self.dx;
        if q.is_nan() {
            return 0.0;
        }
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let (i, _) = clamp_index(q.x, nx);
        let (j, _) = clamp_index(q.y, ny);
        let (k, _) = clamp_index(q.z, nz);
        let lo_x = if i > 0 {
            self.vx[i - 1 + (nx - 1) * (j + ny * k)]
        } else {
            0.0
        };
        let hi_x = if i + 1 < nx {
            self.vx[i + (nx - 1) * (j + ny * k)]
        } else {
            0.0
        };
        let lo_y = if j > 0 {
            self.vy[i + nx * (j - 1 + (ny - 1) * k)]
        } else {
            0.0
        };
        let hi_y = if j + 1 < ny {
            self.vy[i + nx * (j + (ny - 1) * k)]
        } else {
            0.0
        };
        let lo_z = if k > 0 {
            self.vz[i + nx * (j + ny * (k - 1))]
        } else {
            0.0
        };
        let hi_z = if k + 1 < nz {
            self.vz[i + nx * (j + ny * k)]
        } else {
            0.0
        };
        let rc2 = self.density * self.speed * self.speed;
        -rc2 / self.dx
            * (wall_weight(i, nx) * (hi_x - lo_x)
                + wall_weight(j, ny) * (hi_y - lo_y)
                + wall_weight(k, nz) * (hi_z - lo_z))
    }
}

/// Two on a wall, one inside. A node on a boundary owns half the control volume of one inside,
/// so its divergence is divided by half the spacing — per axis, so a corner node picks up the
/// factor three times.
fn wall_weight(k: usize, n: usize) -> f64 {
    if k == 0 || k + 1 == n {
        2.0
    } else {
        1.0
    }
}

/// The node below a position in grid units, and the fraction beyond it.
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
