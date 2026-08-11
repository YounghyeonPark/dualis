//! A rectangular box with conducting walls, marched by the Yee update.

use dualis_core::conserved::quantity;
use dualis_core::{Domain, Exchange, Kind, Ledger, Reading, Violation};
use dualis_units::{Energy, Frequency, Length, LengthVec, Time};
use glam::DVec3;

use crate::{cavity_frequency, Medium};

/// The Courant number a three-dimensional Yee grid is stable up to.
///
/// `dt ≤ dx / (c√3)`, so the number is `1/√3`. The same shape as the acoustic wave equation's, and
/// for the same reason: three axes each contribute a term to the amplification factor and the
/// worst case is the mode that is sharpest along all three at once.
pub const COURANT_3D: f64 = 0.577_350_269_189_625_8;

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
    /// Whether `H` has been advanced the half step the leapfrog opens with.
    started: bool,
    dissipated: f64,
    saved: Option<Box<Saved>>,
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
        Cavity {
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
            started: false,
            dissipated: 0.0,
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

    /// Sum `f²` over one component array, with the boundary half-cells weighted as such.
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
        let e = self.medium.permittivity();
        // Ex is cell-centred in x and node-centred in y and z; the others rotate that.
        let sum = self.weighted_square(&self.ex, (nx, ny + 1, nz + 1), (false, true, true))
            + self.weighted_square(&self.ey, (nx + 1, ny, nz + 1), (true, false, true))
            + self.weighted_square(&self.ez, (nx + 1, ny + 1, nz), (true, true, false));
        Energy::from_si(0.5 * e * sum * v)
    }

    /// Magnetic energy, `½∫μH² dV`.
    pub fn magnetic_energy(&self) -> Energy {
        let (nx, ny, nz) = self.counts;
        let v = self.dx.powi(3);
        let mu = self.medium.permeability();
        let sum = self.weighted_square(&self.hx, (nx + 1, ny, nz), (true, false, false))
            + self.weighted_square(&self.hy, (nx, ny + 1, nz), (false, true, false))
            + self.weighted_square(&self.hz, (nx, ny, nz + 1), (false, false, true));
        Energy::from_si(0.5 * mu * sum * v)
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
        let f = dt / (self.medium.permeability() * self.dx);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..=nx {
                    let curl = (self.ez[self.iez(i, j + 1, k)] - self.ez[self.iez(i, j, k)])
                        - (self.ey[self.iey(i, j, k + 1)] - self.ey[self.iey(i, j, k)]);
                    let idx = self.ihx(i, j, k);
                    self.hx[idx] -= f * curl;
                }
            }
        }
        for k in 0..nz {
            for j in 0..=ny {
                for i in 0..nx {
                    let curl = (self.ex[self.iex(i, j, k + 1)] - self.ex[self.iex(i, j, k)])
                        - (self.ez[self.iez(i + 1, j, k)] - self.ez[self.iez(i, j, k)]);
                    let idx = self.ihy(i, j, k);
                    self.hy[idx] -= f * curl;
                }
            }
        }
        for k in 0..=nz {
            for j in 0..ny {
                for i in 0..nx {
                    let curl = (self.ey[self.iey(i + 1, j, k)] - self.ey[self.iey(i, j, k)])
                        - (self.ex[self.iex(i, j + 1, k)] - self.ex[self.iex(i, j, k)]);
                    let idx = self.ihz(i, j, k);
                    self.hz[idx] -= f * curl;
                }
            }
        }
    }

    /// Advance `E` by `dt`, from `H`. Interior edges only; the walls are held by
    /// [`Cavity::enforce_walls`].
    fn advance_electric(&mut self, dt: f64) -> f64 {
        let (nx, ny, nz) = self.counts;
        let eps = self.medium.permittivity();
        let sigma = self.medium.conductivity.max(0.0);
        // The semi-implicit form, which is stable for any conductivity rather than only for
        // `dt < 2ε/σ`. An explicit `E -= (σ/ε)E dt` is the obvious version and it explodes in a
        // good conductor at any step this grid would use.
        let decay = (1.0 - sigma * dt / (2.0 * eps)) / (1.0 + sigma * dt / (2.0 * eps));
        let gain = (dt / (eps * self.dx)) / (1.0 + sigma * dt / (2.0 * eps));
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
                    let before = self.ez[idx];
                    self.ez[idx] = decay * before + gain * curl;
                    let mid = 0.5 * (before + self.ez[idx]);
                    lost += sigma * mid * mid * dt * volume;
                }
            }
        }
        if sigma > 0.0 {
            lost.max(0.0)
        } else {
            0.0
        }
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
        self.advance_magnetic(h);
        self.dissipated += self.advance_electric(h);
        self.enforce_walls();
        Ok(())
    }

    /// The field energy, plus whatever a lossy medium has turned into heat.
    ///
    /// Arranged so the total is constant: in a lossless box nothing leaves, and in a lossy one what
    /// the field lost is on the books beside it rather than gone.
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.energy().to_si() + self.dissipated)
    }

    fn readings(&self) -> Vec<Reading> {
        vec![
            Reading::new(&self.name, "field energy", self.energy().to_si(), "J"),
            Reading::new(&self.name, "electric", self.electric_energy().to_si(), "J"),
            Reading::new(&self.name, "magnetic", self.magnetic_energy().to_si(), "J"),
            Reading::new(&self.name, "div B", self.magnetic_divergence(), ""),
            Reading::new(&self.name, "dissipated", self.dissipated, "J"),
        ]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
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
