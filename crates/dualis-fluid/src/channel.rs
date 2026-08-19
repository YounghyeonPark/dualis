//! A box of fluid on a MAC grid, marched by projection.

use dualis_core::conserved::quantity;
use dualis_core::{Domain, Exchange, Kind, Ledger, Reading, Violation};
use dualis_units::{Energy, Length, LengthVec, Time, Velocity};
use glam::DVec3;

use crate::Fluid;

/// How many conjugate-gradient iterations per cell the pressure solve gets before it gives up.
const ITERATION_BUDGET: usize = 4;

/// The largest cell Reynolds number central differences stay stable at.
///
/// Two, and it is a property of the **mesh** rather than of the step: advection sharpens what
/// viscosity smooths, and past this the cell is too coarse for the smoothing to keep up. No amount
/// of shortening the time step helps. The symptom is a sawtooth that a reader takes for turbulence.
pub const CELL_REYNOLDS_LIMIT: f64 = 2.0;

/// What the `y` faces of the box are.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Walls {
    /// Periodic in every direction. What Taylor–Green lives in.
    None,
    /// No-slip walls at `y = 0` and `y = h`, moving at the given speeds along `x`.
    ///
    /// Both zero is a channel; one moving is Couette flow.
    Sliding {
        /// Speed of the `y = 0` wall along `x`.
        low: f64,
        /// Speed of the `y = h` wall along `x`.
        high: f64,
    },
}

/// A rectangular box of incompressible fluid.
///
/// # The grid
///
/// Velocities on cell faces and pressure at cell centres. `u` sits on the `x` faces, `v` on the
/// `y` faces and `w` on the `z` faces, so a divergence lands at a cell centre and a pressure
/// gradient lands on a face, with no interpolation in either.
///
/// `x` and `z` are always periodic. `y` is periodic or walled, by [`Walls`].
#[derive(Clone, Debug)]
pub struct Channel {
    name: String,
    counts: (usize, usize, usize),
    dx: f64,
    fluid: Fluid,
    walls: Walls,
    /// Body force per unit mass, m/s².
    force: DVec3,
    /// Work the body force has done on the fluid, in joules.
    ///
    /// **A driven channel is not a closed system.** `drive` is a pressure gradient written as what
    /// it does, and the pump behind it is outside this domain — so the kinetic energy it puts in
    /// arrived from nowhere the bus can see, and the audit is right to notice. Counted here so the
    /// books close, the way `Solid3D` counts what it generates and `Conductor` what it gives away.
    ///
    /// Viscosity takes it straight back out again as heat, which this domain does not model and
    /// does not pretend to: at steady state the drive's power and the dissipation are equal, the
    /// kinetic energy stops moving, and this counter keeps climbing at exactly the rate the fluid
    /// is warming somewhere that is not here.
    driven: f64,
    /// The saved counterpart of [`Channel::driven`].
    saved_driven: f64,
    /// `u` on x faces: `nx · ny · nz`, periodic in x.
    u: Vec<f64>,
    /// `v` on y faces: `nx · (ny+1) · nz`.
    v: Vec<f64>,
    /// `w` on z faces: `nx · ny · nz`, periodic in z.
    w: Vec<f64>,
    /// Pressure at centres: `nx · ny · nz`.
    p: Vec<f64>,
    tolerance: f64,
    residual: f64,
    converged: bool,
    saved: Option<Box<Saved>>,
}

#[derive(Clone, Debug)]
struct Saved {
    u: Vec<f64>,
    v: Vec<f64>,
    w: Vec<f64>,
    p: Vec<f64>,
}

impl dualis_core::ScalarField for Channel {
    /// Metres per second — a **speed**, not a velocity. See [`Channel::as_field`].
    fn unit(&self) -> &'static str {
        "m/s"
    }

    /// The speed at `p`, from the cell it falls in.
    ///
    /// Nearest cell rather than trilinear, and the reason is the grid: velocity lives on the faces
    /// and [`velocity_at`](Channel::velocity_at) already averages the pair across each cell to get
    /// a centred vector. Interpolating that average again would smooth a field that has already
    /// been smoothed once, and a viewer would be looking at a blur of a blur.
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let (nx, ny, nz) = self.counts;
        let q = p.to_si() / self.dx;
        if q.is_nan() {
            return f64::NAN;
        }
        let pick = |v: f64, n: usize| (v.floor().max(0.0) as usize).min(n.saturating_sub(1));
        self.velocity_at(pick(q.x, nx), pick(q.y, ny), pick(q.z, nz))
            .length()
    }
}

impl Channel {
    /// A box of `counts` cubic cells of side `cell`, at rest.
    pub fn new(
        name: impl Into<String>,
        counts: (usize, usize, usize),
        cell: Length,
        fluid: Fluid,
        walls: Walls,
    ) -> Channel {
        let counts = (counts.0.max(1), counts.1.max(1), counts.2.max(1));
        let (nx, ny, nz) = counts;
        Channel {
            name: name.into(),
            counts,
            dx: cell.to_si(),
            fluid,
            walls,
            force: DVec3::ZERO,
            driven: 0.0,
            saved_driven: 0.0,
            u: vec![0.0; nx * ny * nz],
            v: vec![0.0; nx * (ny + 1) * nz],
            w: vec![0.0; nx * ny * nz],
            p: vec![0.0; nx * ny * nz],
            tolerance: 1e-12,
            residual: f64::INFINITY,
            converged: true,
            saved: None,
        }
    }

    /// Cells along each axis.
    pub fn counts(&self) -> (usize, usize, usize) {
        self.counts
    }

    /// The cell side.
    pub fn cell(&self) -> Length {
        Length::from_si(self.dx)
    }

    /// The box's dimensions.
    pub fn size(&self) -> LengthVec {
        LengthVec::from_si(
            DVec3::new(
                self.counts.0 as f64,
                self.counts.1 as f64,
                self.counts.2 as f64,
            ) * self.dx,
        )
    }

    /// What it is full of.
    pub fn fluid(&self) -> Fluid {
        self.fluid
    }

    /// The gap between the walls, for a walled box.
    pub fn gap(&self) -> Length {
        Length::from_si(self.counts.1 as f64 * self.dx)
    }

    /// Drive the flow with a uniform body force per unit mass, in m/s².
    ///
    /// A pressure gradient in disguise, and the form every closed form for channel flow is written
    /// in: a periodic box cannot carry a mean pressure gradient, so the drive has to be a force.
    pub fn drive(&mut self, force: DVec3) -> &mut Channel {
        self.force = force;
        self
    }

    /// Set every `u` face from a function of position, for releasing an exact solution.
    ///
    /// The callback is handed the face's own centre. `v` and `w` are set the same way by
    /// [`Channel::set_velocity`], which takes all three at once and samples each component where
    /// that component lives — the staggering is the whole point and a function evaluated at one
    /// place for all three would be a different field.
    pub fn set_velocity(&mut self, field: impl Fn(DVec3) -> DVec3) -> &mut Channel {
        let (nx, ny, nz) = self.counts;
        let h = self.dx;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let at = DVec3::new(i as f64, j as f64 + 0.5, k as f64 + 0.5) * h;
                    let at_i = self.iu(i, j, k);
                    self.u[at_i] = field(at).x;
                    let at = DVec3::new(i as f64 + 0.5, j as f64 + 0.5, k as f64) * h;
                    let at_i = self.iw(i, j, k);
                    self.w[at_i] = field(at).z;
                }
            }
        }
        for k in 0..nz {
            for j in 0..=ny {
                for i in 0..nx {
                    let at = DVec3::new(i as f64 + 0.5, j as f64, k as f64 + 0.5) * h;
                    let at_i = self.iv(i, j, k);
                    self.v[at_i] = field(at).y;
                }
            }
        }
        self.apply_walls();
        self
    }

    /// The velocity at a cell centre, by averaging the faces around it.
    pub fn velocity_at(&self, i: usize, j: usize, k: usize) -> DVec3 {
        let (nx, ny, nz) = self.counts;
        let (i, j, k) = (i.min(nx - 1), j.min(ny - 1), k.min(nz - 1));
        DVec3::new(
            0.5 * (self.u[self.iu(i, j, k)] + self.u[self.iu((i + 1) % nx, j, k)]),
            0.5 * (self.v[self.iv(i, j, k)] + self.v[self.iv(i, j + 1, k)]),
            0.5 * (self.w[self.iw(i, j, k)] + self.w[self.iw(i, j, (k + 1) % nz)]),
        )
    }

    /// The mean `x` velocity of the whole box.
    pub fn mean_speed(&self) -> Velocity {
        Velocity::from_si(self.u.iter().sum::<f64>() / self.u.len() as f64)
    }

    /// The mean `x` velocity of one layer of cells, which is what a profile is made of.
    pub fn layer_speed(&self, j: usize) -> Velocity {
        let (nx, ny, nz) = self.counts;
        let j = j.min(ny - 1);
        let mut sum = 0.0;
        for k in 0..nz {
            for i in 0..nx {
                sum += self.u[self.iu(i, j, k)];
            }
        }
        Velocity::from_si(sum / (nx * nz) as f64)
    }

    /// Kinetic energy, `½∫ρ|u|²dV`, from the face values.
    pub fn kinetic_energy(&self) -> Energy {
        let cell = self.dx.powi(3);
        let rho = self.fluid.density.to_si();
        // Face values, each carrying the half cell either side of it. `v` at a wall carries only
        // the half inside, which is why its ends are halved.
        let (_, ny, _) = self.counts;
        let sum_u: f64 = self.u.iter().map(|a| a * a).sum();
        let sum_w: f64 = self.w.iter().map(|a| a * a).sum();
        let mut sum_v = 0.0;
        for (idx, val) in self.v.iter().enumerate() {
            let j = (idx / self.counts.0) % (ny + 1);
            let weight = if j == 0 || j == ny { 0.5 } else { 1.0 };
            sum_v += weight * val * val;
        }
        Energy::from_si(0.5 * rho * (sum_u + sum_v + sum_w) * cell)
    }

    /// Total `x` momentum, `ρ∫u dV`.
    ///
    /// Conserved **exactly** in a periodic box with no force and no walls: the advection is in flux
    /// form, so every face's contribution appears twice with opposite signs, and the pressure
    /// gradient of a periodic field sums to zero. That is a machine-precision statement and it is
    /// what a decay rate is too coarse to check.
    pub fn momentum_x(&self) -> f64 {
        self.fluid.density.to_si() * self.dx.powi(3) * self.u.iter().sum::<f64>()
    }

    /// The largest `|∇·u|` anywhere, times the cell — a velocity, so it compares to the flow.
    ///
    /// After the projection this is the pressure solve's residual and nothing else. Weaker than
    /// electromagnetism's divergence identity, which holds exactly; here it holds to whatever the
    /// solve was asked for, and reporting it is the difference between knowing that and assuming
    /// it.
    pub fn divergence(&self) -> f64 {
        let (nx, ny, nz) = self.counts;
        let mut worst: f64 = 0.0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    worst = worst.max(self.divergence_at(i, j, k).abs());
                }
            }
        }
        worst * self.dx
    }

    fn divergence_at(&self, i: usize, j: usize, k: usize) -> f64 {
        let (nx, _, nz) = self.counts;
        ((self.u[self.iu((i + 1) % nx, j, k)] - self.u[self.iu(i, j, k)])
            + (self.v[self.iv(i, j + 1, k)] - self.v[self.iv(i, j, k)])
            + (self.w[self.iw(i, j, (k + 1) % nz)] - self.w[self.iw(i, j, k)]))
            / self.dx
    }

    /// The cell Reynolds number, `|u|dx/ν`.
    ///
    /// A property of the mesh and the flow in it, not of the step. See [`CELL_REYNOLDS_LIMIT`].
    pub fn cell_reynolds(&self) -> f64 {
        self.peak_speed() * self.dx / self.fluid.kinematic_viscosity.to_si()
    }

    /// The fastest face velocity anywhere.
    pub fn peak_speed(&self) -> f64 {
        self.u
            .iter()
            .chain(&self.v)
            .chain(&self.w)
            .fold(0.0f64, |a, b| a.max(b.abs()))
    }

    /// The viscous limit, `dx²/(6ν)` — the same Fourier number conduction has.
    pub fn viscous_limit(&self) -> Time {
        Time::from_si(self.dx * self.dx / (6.0 * self.fluid.kinematic_viscosity.to_si()))
    }

    /// The advective limit, `dx/|u|max`. Infinite for a fluid at rest.
    pub fn courant_limit(&self) -> Time {
        let speed = self.peak_speed();
        Time::from_si(if speed > 0.0 {
            self.dx / speed
        } else {
            f64::INFINITY
        })
    }

    /// Whether the last pressure solve met its tolerance.
    pub fn converged(&self) -> bool {
        self.converged
    }

    /// The relative residual the last pressure solve reached.
    pub fn residual(&self) -> f64 {
        self.residual
    }

    // --- indexing ----------------------------------------------------------

    fn iu(&self, i: usize, j: usize, k: usize) -> usize {
        let (nx, ny, _) = self.counts;
        i + nx * (j + ny * k)
    }
    fn iv(&self, i: usize, j: usize, k: usize) -> usize {
        let (nx, ny, _) = self.counts;
        i + nx * (j + (ny + 1) * k)
    }
    fn iw(&self, i: usize, j: usize, k: usize) -> usize {
        let (nx, ny, _) = self.counts;
        i + nx * (j + ny * k)
    }

    /// The `u` value one cell above or below `j`, through the wall if there is one.
    ///
    /// A no-slip wall is enforced by reflection: the ghost value is `2·u_wall − u_inside`, so the
    /// interpolated value **at** the wall is `u_wall` exactly. Setting the ghost to `u_wall`
    /// instead puts the no-slip condition half a cell into the fluid, which is a first-order error
    /// dressed as a boundary condition and would put the Poiseuille profile visibly off.
    fn u_at(&self, i: usize, j: isize, k: usize) -> f64 {
        let (_, ny, _) = self.counts;
        match self.walls {
            Walls::None => {
                let jj = ((j % ny as isize) + ny as isize) as usize % ny;
                self.u[self.iu(i, jj, k)]
            }
            Walls::Sliding { low, high } => {
                if j < 0 {
                    2.0 * low - self.u[self.iu(i, 0, k)]
                } else if j >= ny as isize {
                    2.0 * high - self.u[self.iu(i, ny - 1, k)]
                } else {
                    self.u[self.iu(i, j as usize, k)]
                }
            }
        }
    }

    fn w_at(&self, i: usize, j: isize, k: usize) -> f64 {
        let (_, ny, _) = self.counts;
        match self.walls {
            Walls::None => {
                let jj = ((j % ny as isize) + ny as isize) as usize % ny;
                self.w[self.iw(i, jj, k)]
            }
            Walls::Sliding { .. } => {
                if j < 0 {
                    -self.w[self.iw(i, 0, k)]
                } else if j >= ny as isize {
                    -self.w[self.iw(i, ny - 1, k)]
                } else {
                    self.w[self.iw(i, j as usize, k)]
                }
            }
        }
    }

    /// Zero the through-wall velocity, which is the only condition `v` has.
    fn apply_walls(&mut self) {
        if let Walls::Sliding { .. } = self.walls {
            let (nx, ny, nz) = self.counts;
            for k in 0..nz {
                for i in 0..nx {
                    let (a, b) = (self.iv(i, 0, k), self.iv(i, ny, k));
                    self.v[a] = 0.0;
                    self.v[b] = 0.0;
                }
            }
        } else {
            // Periodic in y: the two faces **are** the same face, so the high one is a copy of the
            // low one and not an average with it.
            //
            // Averaging was the first version, and it is a half-step lag rather than a boundary
            // condition: the update writes `v[0]` and leaves `v[ny]` stale, so the mean moves
            // `v[0]` only half as far as the physics did. It cost 4.7% of a Taylor-Green decay
            // rate — visible only because that rate has a closed form to be 4.7% away from.
            let (nx, ny, nz) = self.counts;
            for k in 0..nz {
                for i in 0..nx {
                    let (a, b) = (self.iv(i, 0, k), self.iv(i, ny, k));
                    self.v[b] = self.v[a];
                }
            }
        }
    }
    /// `v` at a `y` index that may be outside, resolved by the wall rule.
    fn v_at(&self, i: usize, j: isize, k: usize) -> f64 {
        let (_, ny, _) = self.counts;
        match self.walls {
            Walls::None => {
                let jj = (((j % ny as isize) + ny as isize) % ny as isize) as usize;
                self.v[self.iv(i, jj, k)]
            }
            // A wall has no flow through it, so the face itself is zero and a `v` beyond it is the
            // reflection of the one inside.
            Walls::Sliding { .. } => {
                if j < 0 {
                    -self.v[self.iv(i, 1, k)]
                } else if j > ny as isize {
                    -self.v[self.iv(i, ny - 1, k)]
                } else {
                    self.v[self.iv(i, j as usize, k)]
                }
            }
        }
    }

    /// Which `v` faces are free to move: all of them when periodic, the interior when walled.
    fn v_interior(&self) -> (usize, usize) {
        let (_, ny, _) = self.counts;
        match self.walls {
            Walls::None => (0, ny),
            Walls::Sliding { .. } => (1, ny),
        }
    }

    /// Advection, diffusion and the body force, into a provisional velocity.
    ///
    /// # Flux form, and what it buys
    ///
    /// The advection is `div(uu)` rather than `u.grad u`. They are the same thing for a
    /// divergence-free field and they are **not** the same discretisation: in flux form every
    /// face's contribution appears twice with opposite signs, so total momentum changes only by
    /// what the boundaries and the force do — exactly, rather than to within a truncation error.
    ///
    /// Central differences, second order, no dissipation of their own. The price is
    /// [`CELL_REYNOLDS_LIMIT`]: with nothing damping what advection sharpens, a mesh too coarse for
    /// the viscosity goes unstable and no time step rescues it. Upwinding would trade that for a
    /// numerical viscosity often larger than the real one, which is how a scheme comes to report a
    /// Reynolds number it is not running at.
    fn advance(&mut self, dt: f64) {
        let (nx, ny, nz) = self.counts;
        let h = self.dx;
        let nu = self.fluid.kinematic_viscosity.to_si();
        let (mut du, mut dv, mut dw) = (
            vec![0.0; self.u.len()],
            vec![0.0; self.v.len()],
            vec![0.0; self.w.len()],
        );
        let left = |i: usize| (i + nx - 1) % nx;
        let right = |i: usize| (i + 1) % nx;
        let back = |k: usize| (k + nz - 1) % nz;
        let front = |k: usize| (k + 1) % nz;

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let jj = j as isize;

                    // ---- u, on the x face at (i, j+1/2, k+1/2) -----------------------------
                    let uc = self.u_at(i, jj, k);
                    let ue = 0.5 * (uc + self.u_at(right(i), jj, k));
                    let uw = 0.5 * (self.u_at(left(i), jj, k) + uc);
                    let duudx = (ue * ue - uw * uw) / h;

                    let v_up = 0.5 * (self.v_at(i, jj + 1, k) + self.v_at(left(i), jj + 1, k));
                    let v_dn = 0.5 * (self.v_at(i, jj, k) + self.v_at(left(i), jj, k));
                    let u_up = 0.5 * (uc + self.u_at(i, jj + 1, k));
                    let u_dn = 0.5 * (self.u_at(i, jj - 1, k) + uc);
                    let duvdy = (u_up * v_up - u_dn * v_dn) / h;

                    let w_f = 0.5
                        * (self.w[self.iw(i, j, front(k))] + self.w[self.iw(left(i), j, front(k))]);
                    let w_b = 0.5 * (self.w[self.iw(i, j, k)] + self.w[self.iw(left(i), j, k)]);
                    let u_f = 0.5 * (uc + self.u_at(i, jj, front(k)));
                    let u_b = 0.5 * (self.u_at(i, jj, back(k)) + uc);
                    let duwdz = (u_f * w_f - u_b * w_b) / h;

                    let lap = (self.u_at(right(i), jj, k)
                        + self.u_at(left(i), jj, k)
                        + self.u_at(i, jj + 1, k)
                        + self.u_at(i, jj - 1, k)
                        + self.u_at(i, jj, front(k))
                        + self.u_at(i, jj, back(k))
                        - 6.0 * uc)
                        / (h * h);
                    let at = self.iu(i, j, k);
                    du[at] = dt * (-(duudx + duvdy + duwdz) + nu * lap + self.force.x);

                    // ---- w, on the z face at (i+1/2, j+1/2, k) -----------------------------
                    let wc = self.w[self.iw(i, j, k)];
                    let wf = 0.5 * (wc + self.w[self.iw(i, j, front(k))]);
                    let wb = 0.5 * (self.w[self.iw(i, j, back(k))] + wc);
                    let dwwdz = (wf * wf - wb * wb) / h;

                    let u_e = 0.5 * (self.u_at(right(i), jj, k) + self.u_at(right(i), jj, back(k)));
                    let u_w2 = 0.5 * (self.u_at(i, jj, k) + self.u_at(i, jj, back(k)));
                    let w_e = 0.5 * (wc + self.w[self.iw(right(i), j, k)]);
                    let w_w = 0.5 * (self.w[self.iw(left(i), j, k)] + wc);
                    let dwudx = (w_e * u_e - w_w * u_w2) / h;

                    let v_up2 = 0.5 * (self.v_at(i, jj + 1, k) + self.v_at(i, jj + 1, back(k)));
                    let v_dn2 = 0.5 * (self.v_at(i, jj, k) + self.v_at(i, jj, back(k)));
                    let w_up = 0.5 * (wc + self.w_at(i, jj + 1, k));
                    let w_dn = 0.5 * (self.w_at(i, jj - 1, k) + wc);
                    let dwvdy = (w_up * v_up2 - w_dn * v_dn2) / h;

                    let lap = (self.w[self.iw(right(i), j, k)]
                        + self.w[self.iw(left(i), j, k)]
                        + self.w_at(i, jj + 1, k)
                        + self.w_at(i, jj - 1, k)
                        + self.w[self.iw(i, j, front(k))]
                        + self.w[self.iw(i, j, back(k))]
                        - 6.0 * wc)
                        / (h * h);
                    let at = self.iw(i, j, k);
                    dw[at] = dt * (-(dwwdz + dwudx + dwvdy) + nu * lap + self.force.z);
                }
            }
        }

        let (lo, hi) = self.v_interior();
        for k in 0..nz {
            for j in lo..hi {
                for i in 0..nx {
                    let jj = j as isize;
                    let vc = self.v_at(i, jj, k);
                    let vu = 0.5 * (vc + self.v_at(i, jj + 1, k));
                    let vd = 0.5 * (self.v_at(i, jj - 1, k) + vc);
                    let dvvdy = (vu * vu - vd * vd) / h;

                    let u_e = 0.5 * (self.u_at(right(i), jj, k) + self.u_at(right(i), jj - 1, k));
                    let u_w = 0.5 * (self.u_at(i, jj, k) + self.u_at(i, jj - 1, k));
                    let v_e = 0.5 * (vc + self.v_at(right(i), jj, k));
                    let v_w = 0.5 * (self.v_at(left(i), jj, k) + vc);
                    let dvudx = (v_e * u_e - v_w * u_w) / h;

                    let w_f = 0.5 * (self.w_at(i, jj, front(k)) + self.w_at(i, jj - 1, front(k)));
                    let w_b = 0.5 * (self.w_at(i, jj, k) + self.w_at(i, jj - 1, k));
                    let v_f = 0.5 * (vc + self.v_at(i, jj, front(k)));
                    let v_b = 0.5 * (self.v_at(i, jj, back(k)) + vc);
                    let dvwdz = (v_f * w_f - v_b * w_b) / h;

                    let lap = (self.v_at(right(i), jj, k)
                        + self.v_at(left(i), jj, k)
                        + self.v_at(i, jj + 1, k)
                        + self.v_at(i, jj - 1, k)
                        + self.v_at(i, jj, front(k))
                        + self.v_at(i, jj, back(k))
                        - 6.0 * vc)
                        / (h * h);
                    let at = self.iv(i, j, k);
                    dv[at] = dt * (-(dvvdy + dvudx + dvwdz) + nu * lap + self.force.y);
                }
            }
        }

        for (a, b) in self.u.iter_mut().zip(&du) {
            *a += b;
        }
        for (a, b) in self.v.iter_mut().zip(&dv) {
            *a += b;
        }
        for (a, b) in self.w.iter_mut().zip(&dw) {
            *a += b;
        }
        self.apply_walls();
    }

    /// Make the provisional velocity divergence-free by subtracting a pressure gradient.
    ///
    /// Solve `lap(phi) = div(u*)/dt` and set `u = u* - dt grad(phi)`. The operator is the same
    /// finite-volume Laplacian `Conductor` and `Puck` solve, with no flux through a wall and
    /// periodic elsewhere, and it is singular by a constant — the pressure of an incompressible
    /// flow is defined only up to one. The right-hand side has its mean removed so the system is
    /// consistent, and the answer has its mean removed so it is a particular one.
    fn project(&mut self, dt: f64) -> bool {
        let (nx, ny, nz) = self.counts;
        let cells = nx * ny * nz;
        let mut b = vec![0.0; cells];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    b[i + nx * (j + ny * k)] = self.divergence_at(i, j, k) / dt;
                }
            }
        }
        let mean: f64 = b.iter().sum::<f64>() / cells as f64;
        for v in b.iter_mut() {
            *v -= mean;
        }

        let mut x = std::mem::take(&mut self.p);
        if x.len() != cells {
            x = vec![0.0; cells];
        }
        let ax = self.laplacian(&x);
        let mut r: Vec<f64> = b.iter().zip(&ax).map(|(bi, axi)| bi - axi).collect();
        let mut p = r.clone();
        let mut rr: f64 = r.iter().map(|v| v * v).sum();
        let scale = b
            .iter()
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt()
            .max(f64::MIN_POSITIVE);
        let budget = ITERATION_BUDGET * cells + 64;
        let mut iterations = 0;
        while rr.sqrt() / scale > self.tolerance && iterations < budget {
            let ap = self.laplacian(&p);
            let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
            if pap.abs() <= 0.0 {
                break;
            }
            let alpha = rr / pap;
            for (xi, pi) in x.iter_mut().zip(&p) {
                *xi += alpha * pi;
            }
            for (ri, api) in r.iter_mut().zip(&ap) {
                *ri -= alpha * api;
            }
            let rr_next: f64 = r.iter().map(|v| v * v).sum();
            let beta = rr_next / rr;
            for (pi, ri) in p.iter_mut().zip(&r) {
                *pi = ri + beta * *pi;
            }
            rr = rr_next;
            iterations += 1;
        }
        let phi_mean: f64 = x.iter().sum::<f64>() / cells as f64;
        for v in x.iter_mut() {
            *v -= phi_mean;
        }
        self.p = x;
        self.residual = rr.sqrt() / scale;
        self.converged = self.residual <= self.tolerance;

        let h = self.dx;
        let phi = self.p.clone();
        let idx = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let at = self.iu(i, j, k);
                    self.u[at] -= dt * (phi[idx(i, j, k)] - phi[idx((i + nx - 1) % nx, j, k)]) / h;
                    let at = self.iw(i, j, k);
                    self.w[at] -= dt * (phi[idx(i, j, k)] - phi[idx(i, j, (k + nz - 1) % nz)]) / h;
                }
            }
        }
        let (lo, hi) = self.v_interior();
        for k in 0..nz {
            for j in lo..hi {
                for i in 0..nx {
                    let up = idx(i, j % ny, k);
                    let dn = idx(i, (j + ny - 1) % ny, k);
                    let at = self.iv(i, j, k);
                    self.v[at] -= dt * (phi[up] - phi[dn]) / h;
                }
            }
        }
        self.apply_walls();
        self.converged
    }

    /// The finite-volume Laplacian at cell centres: no flux through a wall, periodic elsewhere.
    fn laplacian(&self, x: &[f64]) -> Vec<f64> {
        let (nx, ny, nz) = self.counts;
        let h2 = self.dx * self.dx;
        let idx = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
        let mut y = vec![0.0; x.len()];
        let walled = matches!(self.walls, Walls::Sliding { .. });
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let c = idx(i, j, k);
                    let mut acc =
                        x[idx((i + 1) % nx, j, k)] + x[idx((i + nx - 1) % nx, j, k)] - 2.0 * x[c];
                    acc +=
                        x[idx(i, j, (k + 1) % nz)] + x[idx(i, j, (k + nz - 1) % nz)] - 2.0 * x[c];
                    // The `y` direction, where a wall means the neighbour is absent rather than
                    // mirrored: a zero-flux face contributes nothing at all.
                    if walled {
                        if j + 1 < ny {
                            acc += x[idx(i, j + 1, k)] - x[c];
                        }
                        if j > 0 {
                            acc += x[idx(i, j - 1, k)] - x[c];
                        }
                    } else {
                        acc += x[idx(i, (j + 1) % ny, k)] + x[idx(i, (j + ny - 1) % ny, k)]
                            - 2.0 * x[c];
                    }
                    y[c] = acc / h2;
                }
            }
        }
        y
    }
}

impl Domain for Channel {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    fn max_stable_dt(&self, _now: Time) -> Time {
        Time::from_si(
            self.viscous_limit()
                .to_si()
                .min(self.courant_limit().to_si()),
        )
    }

    fn step(&mut self, _t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let h = dt.to_si();
        let limit = self.max_stable_dt(Time::from_si(0.0)).to_si();
        if h > limit {
            return Err(Violation {
                quantity: "flow step".into(),
                site: self.name.clone(),
                before: limit,
                after: h,
                scale: limit,
                tolerance: 0.0,
            });
        }
        let re = self.cell_reynolds();
        if re > CELL_REYNOLDS_LIMIT {
            return Err(Violation {
                quantity: "cell Reynolds number".into(),
                site: self.name.clone(),
                before: CELL_REYNOLDS_LIMIT,
                after: re,
                scale: CELL_REYNOLDS_LIMIT,
                tolerance: 0.0,
            });
        }

        // The work the drive did, counted in the same statement that does it. Measured as the
        // change in kinetic energy across a step with no other source: viscosity takes energy out
        // and the drive puts it in, and at steady state the two are equal — so this keeps climbing
        // while the kinetic energy stops, which is the true statement about a pumped channel.
        let before = self.kinetic_energy().to_si();
        self.advance(h);
        if !self.project(h) {
            return Err(Violation::at(
                self.name.clone(),
                "pressure residual",
                self.residual,
            ));
        }
        if self.force != DVec3::ZERO {
            self.driven += self.kinetic_energy().to_si() - before;
        }
        Ok(())
    }

    /// The kinetic energy the box is holding.
    /// The kinetic energy it holds, less the work the drive has put in.
    ///
    /// Two contributions rather than their difference: `Ledger::add` raises an entry's *scale* to
    /// the largest thing added to it and the audit judges a change against that, so pre-summing a
    /// near-zero net would leave no scale at all. A channel starting from rest holds exactly zero,
    /// and the first `2.9e-12` J the drive did was judged a hundred-percent change and stopped a
    /// correct run on its first step.
    fn ledger(&self) -> Ledger {
        Ledger::new()
            .with(quantity::ENERGY, self.kinetic_energy().to_si())
            .with(quantity::ENERGY, -self.driven)
    }

    fn readings(&self) -> Vec<Reading> {
        let mut out = vec![
            Reading::new(&self.name, "mean speed", self.mean_speed().to_si(), "m/s"),
            Reading::new(&self.name, "peak speed", self.peak_speed(), "m/s"),
            Reading::new(
                &self.name,
                "kinetic energy",
                self.kinetic_energy().to_si(),
                "J",
            ),
            Reading::new(&self.name, "divergence", self.divergence(), "m/s"),
            Reading::new(&self.name, "cell Reynolds", self.cell_reynolds(), ""),
        ];
        // **Only for a channel that is driven**, the way a block reports `melted` only if it can
        // melt. At steady state the kinetic energy stops moving and this keeps climbing, which is
        // the pump's power made visible — and without it a reader would see a flow holding still
        // and no sign that anything was paying for it.
        if self.force != DVec3::ZERO {
            out.push(Reading::new(&self.name, "work driven in", self.driven, "J"));
        }
        out
    }

    /// **Speed**, so a flow can be looked at.
    ///
    /// The honest caveat is in the unit and in this sentence rather than in a refusal to draw one:
    /// velocity is a *vector* and this is its magnitude, so a picture of it shows where the fluid
    /// is moving fast and not which way it is going. The components are in the JSON, and
    /// `layer_speed` is what a profile should be read from.
    ///
    /// Drawn at all because a domain nobody can see is a domain nobody trusts, and this crate's own
    /// documentation says why that matters here more than anywhere: "it looks like a fluid" is the
    /// easiest wrong answer in computational physics to accept, and the answer to that is closed
    /// forms **and** a picture, not one instead of the other.
    fn as_field(&self) -> Option<&dyn dualis_core::ScalarField> {
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
            u: self.u.clone(),
            v: self.v.clone(),
            w: self.w.clone(),
            p: self.p.clone(),
        }));
        self.saved_driven = self.driven;
    }

    fn restore(&mut self) {
        if let Some(s) = self.saved.take() {
            self.u = s.u;
            self.v = s.v;
            self.w = s.w;
            self.p = s.p;
            self.driven = self.saved_driven;
        }
    }
}
