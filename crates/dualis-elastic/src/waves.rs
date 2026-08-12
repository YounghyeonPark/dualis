//! Elastic waves: the same body, given inertia.
//!
//! [`Block`](crate::Block) solves `∇·σ = 0` — an equilibrium, in zero time, with the density it was
//! given going unused. This is the other half of the same operator:
//!
//! ```text
//!   ρ ∂²u/∂t² = ∇·σ,   σ = λ tr(ε) I + 2μ ε
//! ```
//!
//! A separate type rather than a mode on `Block`, and that is not squeamishness: `Block` is
//! [`Kind::QuasiStatic`] and its own documentation says there is no state to roll forward. That is
//! true about it. A body with velocity in it has a different lifecycle, a different stability limit
//! and a different thing to conserve, and pretending one type is both would make every one of those
//! depend on which method was called last.
//!
//! What they share is the element, which is the part worth sharing: the same 24×24 trilinear
//! stiffness assembled the same way. The static tests check that operator against four exact moduli;
//! this one checks the same operator against two exact speeds.
//!
//! # What it is checked against
//!
//! Two speeds, and their ratio:
//!
//! ```text
//!   c_p = √((λ+2μ)/ρ)      compression, sides held
//!   c_s = √(μ/ρ)           shear
//!   c_p/c_s = √(2(1−ν)/(1−2ν))
//! ```
//!
//! The ratio is the sharp one, because `E` and `ρ` cancel out of it: a scheme with the wrong
//! stiffness *and* the wrong mass can still get one speed right by accident and cannot get the ratio
//! right. [`hold`](Waves::hold) is what makes each speed separately measurable — freezing two
//! displacement components leaves the one-dimensional problem the closed form is about.
//!
//! # The time stepping, and why its dispersion is removed rather than tolerated
//!
//! Central differences, which is the same leapfrog `Room` and `Cavity` use:
//!
//! ```text
//!   u^{n+1} = 2u^n − u^{n−1} − dt² M⁻¹ K u^n
//! ```
//!
//! `M` is **lumped** — each element gives an eighth of its mass to each of its eight nodes — so
//! `M⁻¹` is a division and there is no solve anywhere in a step.
//!
//! On a single eigenmode this recurrence oscillates at exactly `Ω = 2·arcsin(ω dt/2)` per step rather
//! than at `ω dt`, which is the leapfrog's own dispersion and is a property of the time scheme and not
//! of the elasticity. So the tests measure `Ω`, invert that relation to recover `ω`, and compare
//! *that* — leaving only the spatial discretisation, which is second order. A test that compared the
//! raw period would be measuring the two errors added together and calling the sum an accuracy.
//!
//! # What a step costs, measured
//!
//! The 24×24 element loop, and nothing else is close:
//!
//! ```text
//!   16³ elements   4,913 nodes    2.47 ms/step
//!   32³ elements  35,937 nodes   19.69 ms/step
//! ```
//!
//! Linear in nodes, so the two `Vec` allocations a step makes are **under 1%** and are left alone —
//! measured before deciding, because a scratch buffer threaded through `&mut self` would have cost
//! clarity for nothing.
//!
//! What that buys is worth knowing before sizing a run. At a millimetre in aluminium the limit is
//! `1.15e-7 s`, so a wave crossing a 32 mm body takes about **90 steps** — under two seconds of
//! compute, which is the case this is for. A body left *ringing* for a millisecond is 8,700 steps and
//! **three minutes**, which is the case to reach for a modal analysis instead.

use dualis_core::conserved::quantity;
use dualis_core::{
    units::{Energy, Length, LengthVec, Time, Velocity},
    Domain, Exchange, Kind, Ledger, Reading, ScalarField, Violation,
};
use dualis_units::Frequency;

use crate::element::{lame, stiffness, CORNERS, DOF};
use crate::Elastic;

/// One of the three axes.
///
/// # Why a type and not a `usize`
///
/// Because the alternative was measured and it was worse. Five methods on [`Waves`] took an axis as
/// an index, and they disagreed about what to do with a fourth one: four returned silently — leaving
/// a body with three free components where the caller believed it had one, and therefore a wave speed
/// that is not the one being asked about — and `mode_frequency` clamped to 2, which is a different
/// wrong answer to the same mistake.
///
/// Three responses to one bad input in one type is not a policy. An enum makes the input
/// unrepresentable, which is better than any of the three.
///
/// [`Face::axis`](crate::Face::axis) predates this and still returns an index; it is on a published
/// type and changing it is a separate decision from getting this one right.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Axis {
    /// x.
    X,
    /// y.
    Y,
    /// z.
    Z,
}

impl Axis {
    /// 0, 1 or 2 — for indexing the three components of a displacement.
    pub fn index(self) -> usize {
        match self {
            Axis::X => 0,
            Axis::Y => 1,
            Axis::Z => 2,
        }
    }

    /// All three, in order.
    pub const ALL: [Axis; 3] = [Axis::X, Axis::Y, Axis::Z];
}

/// A rectangular body of elastic material, carrying waves.
///
/// Cells are **cubes**, for the reason every grid in this workspace is: an anisotropic element makes
/// the stability limit and the truncation error different along each axis, which would resolve one
/// direction better than another for no physical reason.
///
/// Nothing is held and nothing is moving until something releases a mode or holds a face, so a fresh
/// body sits still forever — which is the correct answer to `ρü = ∇·σ` with `u = 0`.
#[derive(Clone, Debug)]
pub struct Waves {
    name: String,
    /// Elements along each axis.
    counts: (usize, usize, usize),
    /// Nodes along each axis, one more than the elements.
    nodes: (usize, usize, usize),
    dx: f64,
    /// What the block was built from, and element zero's material. See [`Waves::material`].
    material: Elastic,
    /// The materials any element may be, in the order they were introduced. One entry for a block
    /// nobody filled.
    ///
    /// A palette rather than one material per element, because a 24×24 element stiffness is 4.6 kB and
    /// a laminate has two of them however many elements it has. The same shape `Solid3D` uses for the
    /// same reason.
    materials: Vec<Elastic>,
    /// Which palette entry each element is. `u16`, so a block may be made of 65 536 different things,
    /// which is more than a caller will ever write down and a quarter of the memory of a `usize`.
    which: Vec<u16>,
    /// The 24×24 element stiffness of each palette entry.
    kes: Vec<Vec<f64>>,
    /// Displacement now and one step ago, three per node.
    u: Vec<f64>,
    prev: Vec<f64>,
    /// Lumped nodal mass, one per node: `(elements touching it) × ρ dx³ / 8`.
    mass: Vec<f64>,
    /// Whether each degree of freedom is held at zero.
    held: Vec<bool>,
    /// `2/√λ` where `λ` bounds the eigenvalues of `M⁻¹K`, so leapfrog is stable at or under it.
    limit: f64,
    saved: Option<Box<(Vec<f64>, Vec<f64>)>>,
}

impl Waves {
    /// A body of `counts` cubic elements of side `cell`.
    /// Arguments in the same order as [`Block::new`](crate::Block::new) — name, shape, spacing,
    /// material. They disagreed until a stabilisation pass noticed, which is the sort of thing only a
    /// caller writing both in one file ever sees.
    pub fn new(
        name: impl Into<String>,
        counts: (usize, usize, usize),
        cell: Length,
        material: Elastic,
    ) -> Waves {
        let counts = (counts.0.max(1), counts.1.max(1), counts.2.max(1));
        let nodes = (counts.0 + 1, counts.1 + 1, counts.2 + 1);
        let n = nodes.0 * nodes.1 * nodes.2;
        let dx = cell.to_si();
        let (lambda, mu) = lame(material.youngs_modulus.to_si(), material.poisson_ratio);
        let mut w = Waves {
            name: name.into(),
            counts,
            nodes,
            dx,
            material,
            materials: vec![material],
            which: vec![0; counts.0 * counts.1 * counts.2],
            kes: vec![stiffness(dx, lambda, mu)],
            u: vec![0.0; 3 * n],
            prev: vec![0.0; 3 * n],
            mass: vec![0.0; n],
            held: vec![false; 3 * n],
            limit: 0.0,
            saved: None,
        };
        w.resolve();
        w
    }

    /// Elements along each axis.
    pub fn elements(&self) -> (usize, usize, usize) {
        self.counts
    }

    /// Nodes along each axis, one more than the elements.
    pub fn node_counts(&self) -> (usize, usize, usize) {
        self.nodes
    }

    /// The element side.
    pub fn cell(&self) -> Length {
        Length::from_si(self.dx)
    }

    /// The material.
    pub fn material(&self) -> Elastic {
        self.material
    }

    /// Every material this block is made of, in the order they were introduced. Length one unless
    /// [`Waves::fill`] has been called.
    pub fn materials(&self) -> &[Elastic] {
        &self.materials
    }

    /// What one **element** is made of. Elements, not nodes: a material is a property of the volume and
    /// a node sits between up to eight of them.
    pub fn material_at(&self, e_x: usize, e_y: usize, e_z: usize) -> Elastic {
        let (ex, ey, _) = self.counts;
        self.materials[self.which[e_x + ex * (e_y + ey * e_z)] as usize]
    }

    /// Make every element the predicate accepts out of `material`, and report how many changed.
    ///
    /// The predicate takes **element** indices. A `Waves` of `(n, n, n)` has `n³` elements and
    /// `(n+1)³` nodes, and the off-by-one between those is the mistake this signature is shaped to
    /// avoid making silently.
    ///
    /// ```
    /// # use dualis_elastic::{Elastic, Waves};
    /// # use dualis_units::Length;
    /// let mut w = Waves::new("laminate", (1, 1, 8), Length::mm(1.0), Elastic::aluminium_6061());
    /// // Alternating layers, one element thick, perpendicular to z.
    /// let changed = w.fill(Elastic::steel(), |_, _, e_z| e_z % 2 == 1);
    /// assert_eq!(changed, 4);
    /// assert_eq!(w.materials().len(), 2);
    /// ```
    ///
    /// # A fill that changes nothing is not an error
    ///
    /// Unlike a scene's region, which is refused when it selects no elements. The difference is that a
    /// region is a bound somebody typed into a file and an empty one is almost always a typo, whereas
    /// `fill` is called from code with a predicate the caller can read. The **count is returned** so a
    /// caller who wants the region's behaviour can have it, and `a_layered_wave.rs` checks the count
    /// rather than trusting the predicate.
    ///
    /// A block filled with what it already held is unchanged bit for bit: the palette does not grow,
    /// because an equal material resolves to the entry it already has.
    pub fn fill(
        &mut self,
        material: Elastic,
        which: impl Fn(usize, usize, usize) -> bool,
    ) -> usize {
        let index = match self.materials.iter().position(|m| *m == material) {
            Some(i) => i as u16,
            None => {
                let (lambda, mu) = lame(material.youngs_modulus.to_si(), material.poisson_ratio);
                self.materials.push(material);
                self.kes.push(stiffness(self.dx, lambda, mu));
                (self.materials.len() - 1) as u16
            }
        };
        let (ex, ey, ez) = self.counts;
        let mut changed = 0;
        for e_z in 0..ez {
            for e_y in 0..ey {
                for e_x in 0..ex {
                    if !which(e_x, e_y, e_z) {
                        continue;
                    }
                    let e = e_x + ex * (e_y + ey * e_z);
                    if self.which[e] != index {
                        self.which[e] = index;
                        changed += 1;
                    }
                }
            }
        }
        if changed > 0 {
            // The mass and the stability limit both depend on which element is what, and a step sized
            // before a fill would be sized for the wrong block.
            self.resolve();
        }
        changed
    }

    /// Hold one displacement component at zero on **every** node.
    ///
    /// What makes each wave speed separately measurable, and therefore what the closed forms are
    /// about. Freezing two of the three components leaves a one-dimensional problem:
    ///
    /// ```text
    ///   hold x and y, let z move, vary along z    →   ρü = (λ+2μ)u''   →   c_p
    ///   hold y and z, let x move, vary along z    →   ρü = μu''        →   c_s
    /// ```
    ///
    /// The first is the *constrained* compression — the sides cannot bulge because they are held —
    /// which is why it is `λ+2μ` and not `E`. That distinction is the commonest way to get a wave
    /// speed wrong by 20%.
    ///
    pub fn hold(&mut self, axis: Axis) -> &mut Waves {
        let axis = axis.index();
        let n = self.nodes.0 * self.nodes.1 * self.nodes.2;
        for node in 0..n {
            self.held[3 * node + axis] = true;
            self.u[3 * node + axis] = 0.0;
            self.prev[3 * node + axis] = 0.0;
        }
        self.resolve();
        self
    }

    /// Clamp both ends of an axis: every displacement component zero on those two faces.
    ///
    /// What makes a standing mode standing. A free end is a different boundary with a different
    /// closed form — quarter-wave rather than half-wave — and this type offers the one the tests use.
    pub fn clamp_ends(&mut self, axis: Axis) -> &mut Waves {
        let axis = axis.index();
        let (nx, ny, nz) = self.nodes;
        let last = [nx - 1, ny - 1, nz - 1][axis];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let at = [i, j, k][axis];
                    if at == 0 || at == last {
                        let node = i + nx * (j + ny * k);
                        for c in 0..3 {
                            self.held[3 * node + c] = true;
                            self.u[3 * node + c] = 0.0;
                            self.prev[3 * node + c] = 0.0;
                        }
                    }
                }
            }
        }
        self.resolve();
        self
    }

    /// Release a standing half-wave: `amplitude · sin(nπ x_axis / L)` in the `along` component.
    ///
    /// Set as **both** the current and the previous displacement, so the body starts at rest at the
    /// mode's extreme. That is the initial condition the closed form is about: a cosine in time, so
    /// the first quarter period is the whole of the information and no start-up transient is
    /// introduced.
    ///
    /// `mode` counts half-waves, so `1` is the fundamental of a clamped-clamped span.
    pub fn release_mode(&mut self, mode: usize, vary: Axis, along: Axis, amplitude: Length) {
        let (vary, along) = (vary.index(), along.index());
        let (nx, ny, nz) = self.nodes;
        let span = self.counts_of(vary) as f64 * self.dx;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let x = [i, j, k][vary] as f64 * self.dx;
                    let v =
                        amplitude.to_si() * (mode as f64 * std::f64::consts::PI * x / span).sin();
                    let dof = 3 * (i + nx * (j + ny * k)) + along;
                    if !self.held[dof] {
                        self.u[dof] = v;
                        self.prev[dof] = v;
                    }
                }
            }
        }
    }

    /// The amplitude of one mode currently present, by projection onto its shape.
    ///
    /// The counterpart to [`release_mode`](Waves::release_mode), and what makes a frequency
    /// measurable rather than merely visible: a sine is orthogonal to the others on this grid, so for
    /// a body holding one mode this is exact and for one holding several it is the right coefficient.
    pub fn mode_amplitude(&self, mode: usize, vary: Axis, along: Axis) -> f64 {
        let (vary, along) = (vary.index(), along.index());
        let (nx, ny, nz) = self.nodes;
        let span = self.counts_of(vary) as f64 * self.dx;
        let (mut num, mut den) = (0.0, 0.0);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let x = [i, j, k][vary] as f64 * self.dx;
                    let shape = (mode as f64 * std::f64::consts::PI * x / span).sin();
                    num += self.u[3 * (i + nx * (j + ny * k)) + along] * shape;
                    den += shape * shape;
                }
            }
        }
        if den > 0.0 {
            num / den
        } else {
            0.0
        }
    }

    /// The displacement of one node.
    pub fn displacement_at(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        let (nx, ny, nz) = self.nodes;
        let node = i.min(nx - 1) + nx * (j.min(ny - 1) + ny * k.min(nz - 1));
        [self.u[3 * node], self.u[3 * node + 1], self.u[3 * node + 2]]
    }

    /// The closed-form frequency of a clamped-clamped standing mode, `n·c/(2L)`.
    ///
    /// Public because it is what a caller sizing a grid needs before there is anything to measure,
    /// and because stating it here means the tests can compare against something written once.
    /// `speed` is the caller's choice of [`Elastic::p_wave_speed`] or
    /// [`Elastic::s_wave_speed`](Elastic::s_wave_speed) — the type cannot know which wave a given
    /// [`hold`](Waves::hold) has left free.
    pub fn mode_frequency(&self, mode: usize, vary: Axis, speed: Velocity) -> Frequency {
        let span = self.counts_of(vary.index()) as f64 * self.dx;
        Frequency::from_si(mode as f64 * speed.to_si() / (2.0 * span))
    }

    /// The kinetic energy, from the central-difference velocity.
    ///
    /// `(u^n − u^{n−1})/dt` is a **backward** difference and so is the velocity half a step ago, not
    /// now. That half-step offset is why the total energy of a leapfrog swings instead of sitting
    /// still, and the swing is `2 sin(ωΔt/2)` rather than a defect — `Room` records the same thing.
    pub fn kinetic_energy(&self, dt: Time) -> Energy {
        let h = dt.to_si();
        if h <= 0.0 {
            return Energy::from_si(0.0);
        }
        let mut total = 0.0;
        for node in 0..self.mass.len() {
            let mut v2 = 0.0;
            for c in 0..3 {
                let d = (self.u[3 * node + c] - self.prev[3 * node + c]) / h;
                v2 += d * d;
            }
            total += 0.5 * self.mass[node] * v2;
        }
        Energy::from_si(total)
    }

    /// The strain energy, `½ uᵀKu`.
    pub fn strain_energy(&self) -> Energy {
        let ku = self.apply(&self.u);
        Energy::from_si(0.5 * self.u.iter().zip(&ku).map(|(a, b)| a * b).sum::<f64>())
    }

    /// Kinetic plus strain — what a body with no boundary doing work on it holds constant.
    pub fn total_energy(&self, dt: Time) -> Energy {
        Energy::from_si(self.kinetic_energy(dt).to_si() + self.strain_energy().to_si())
    }

    fn counts_of(&self, axis: usize) -> usize {
        [self.counts.0, self.counts.1, self.counts.2][axis]
    }

    /// Rebuild the lumped mass and the stability limit.
    ///
    /// Called by the constructor and by every mutator, on the rule this workspace learned from
    /// `Puck::repack`: a mutator that leaves a cached number stale reports something that was true
    /// about the previous object, and it looks exactly like a number.
    fn resolve(&mut self) {
        let (nx, ny, nz) = self.nodes;
        let (ex, ey, ez) = self.counts;
        let volume = self.dx.powi(3);
        self.mass = vec![0.0; nx * ny * nz];
        // Row sums of |K|, accumulated the same way the mass is, so the two are about the same
        // assembly and the ratio below is a genuine bound rather than two guesses divided.
        let mut rows = vec![0.0; 3 * nx * ny * nz];
        let mut map = [0usize; 8];
        for e_z in 0..ez {
            for e_y in 0..ey {
                for e_x in 0..ex {
                    let e = e_x + ex * (e_y + ey * e_z);
                    let ke = &self.kes[self.which[e] as usize];
                    // An element's mass goes to its eight corners in equal shares, so a node between
                    // two materials gets a share of each. That is what makes a lumped mass right at an
                    // interface rather than an average of the two densities.
                    let per = self.materials[self.which[e] as usize].density.to_si() * volume / 8.0;
                    for (c, corner) in CORNERS.iter().enumerate() {
                        map[c] =
                            (e_x + corner[0]) + nx * ((e_y + corner[1]) + ny * (e_z + corner[2]));
                        self.mass[map[c]] += per;
                    }
                    for a in 0..DOF {
                        let row = 3 * map[a / 3] + a % 3;
                        for b in 0..DOF {
                            rows[row] += ke[a * DOF + b].abs();
                        }
                    }
                }
            }
        }
        // Gershgorin on `M⁻¹K`: every eigenvalue is inside a disc of radius the row sum over the
        // mass, so `λ_max ≤ max_i rows_i / m_i`, and leapfrog needs `dt² λ_max ≤ 4`.
        let mut worst = 0.0f64;
        for node in 0..self.mass.len() {
            for c in 0..3 {
                let dof = 3 * node + c;
                if !self.held[dof] && self.mass[node] > 0.0 {
                    worst = worst.max(rows[dof] / self.mass[node]);
                }
            }
        }
        self.limit = if worst > 0.0 {
            2.0 / worst.sqrt()
        } else {
            f64::INFINITY
        };
    }

    /// `K·x`, assembled element by element — the same assembly [`Block`](crate::Block) uses.
    fn apply(&self, x: &[f64]) -> Vec<f64> {
        let (nx, ny, _) = self.nodes;
        let (ex, ey, ez) = self.counts;
        let mut y = vec![0.0; x.len()];
        let mut map = [0usize; 8];
        for e_z in 0..ez {
            for e_y in 0..ey {
                for e_x in 0..ex {
                    let ke = &self.kes[self.which[e_x + ex * (e_y + ey * e_z)] as usize];
                    for (c, corner) in CORNERS.iter().enumerate() {
                        map[c] =
                            (e_x + corner[0]) + nx * ((e_y + corner[1]) + ny * (e_z + corner[2]));
                    }
                    for a in 0..DOF {
                        let row = 3 * map[a / 3] + a % 3;
                        let mut acc = 0.0;
                        for b in 0..DOF {
                            acc += ke[a * DOF + b] * x[3 * map[b / 3] + b % 3];
                        }
                        y[row] += acc;
                    }
                }
            }
        }
        y
    }
}

impl Domain for Waves {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// `2/√λ_max(M⁻¹K)`, bounded by Gershgorin and computed from the assembled rows.
    ///
    /// Not `dx/(c_p√3)`. That form is the right answer for a three-point stencil per axis and this is
    /// a trilinear element with a lumped mass, whose off-diagonal reach is the whole 24-node
    /// neighbourhood — so the bound is computed from the operator rather than assumed from a
    /// stencil it is not.
    ///
    /// Measured on a 4×4×4 cube of aluminium at a millimetre it comes out **1.229×** the Courant form,
    /// so borrowing `dx/(c_p√3)` would have been safe and 23% wasteful. That is the direction worth
    /// knowing and it is not the direction that was assumed before measuring it.
    fn max_stable_dt(&self, _now: Time) -> Time {
        Time::from_si(self.limit)
    }

    fn step(&mut self, _t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let h = dt.to_si();
        if h > self.limit * (1.0 + 1e-12) {
            return Err(Violation {
                quantity: "leapfrog stability".to_string(),
                site: format!("{} (explicit elastodynamics)", self.name),
                before: self.limit,
                after: h,
                scale: self.limit,
                tolerance: 1e-12,
            });
        }
        let ku = self.apply(&self.u);
        let mut next = vec![0.0; self.u.len()];
        for node in 0..self.mass.len() {
            for c in 0..3 {
                let dof = 3 * node + c;
                if self.held[dof] {
                    continue;
                }
                next[dof] = 2.0 * self.u[dof] - self.prev[dof] - h * h * ku[dof] / self.mass[node];
            }
        }
        self.prev = std::mem::replace(&mut self.u, next);
        Ok(())
    }

    /// The strain energy alone.
    ///
    /// **Not the total**, and the reason is the half-step: the kinetic term needs a `dt` to exist at
    /// all, and `ledger` is not given one. A ledger that guessed a step would report an energy that
    /// depended on a number nobody passed. The audit therefore sees a quantity that legitimately
    /// oscillates, so a scene using this domain has to say so with
    /// `conservation_tolerance_for(ENERGY, ..)` — which is the same bargain `runtime/gpu` asks for and
    /// is better than a total that is quietly wrong by a half step.
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.strain_energy().to_si())
    }

    fn books_balance(&self) -> bool {
        false
    }

    fn readings(&self) -> Vec<Reading> {
        let peak = self.u.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        vec![
            Reading::new(
                &self.name,
                "strain energy",
                self.strain_energy().to_si(),
                "J",
            ),
            Reading::new(&self.name, "peak displacement", peak * 1e9, "nm"),
        ]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    /// The **magnitude** of the displacement, so the analysis layer can draw a wave.
    ///
    /// `Domain::as_field` nominates one scalar and a displacement is a vector, so this is `|u|`:
    /// the one scalar that is nonzero wherever the body is moving, whichever way it is moving. A
    /// component would have been a choice — and the wrong one for a shear wave drawn on the axis it
    /// travels along, which is zero everywhere.
    ///
    /// It costs the sign. A compression and an extension look the same, which for a standing mode
    /// makes a node dark and both antinodes bright — twice the spatial frequency of the mode. That is
    /// the honest trade for one scalar and it is stated rather than discovered from a picture.
    ///
    /// [`Block`](crate::Block) still offers nothing here, which is a real gap and older than this
    /// type: a static solve's displacement is just as drawable. It is on a published type, so
    /// adding it is a separate decision.
    fn as_field(&self) -> Option<&dyn ScalarField> {
        Some(self)
    }

    fn checkpoint(&mut self) {
        self.saved = Some(Box::new((self.u.clone(), self.prev.clone())));
    }

    fn restore(&mut self) {
        if let Some(s) = self.saved.take() {
            self.u = s.0;
            self.prev = s.1;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }
}

impl ScalarField for Waves {
    /// **Metres**, because that is what a displacement is. A view converting to nanometres is making
    /// a presentation choice and this is not the place for it — `readings` reports nanometres because
    /// a number in a table has to be readable, and a field has a legend.
    fn unit(&self) -> &'static str {
        "m"
    }

    /// `|u|` trilinearly interpolated between nodes, clamped at the faces.
    ///
    /// Clamped rather than extrapolated: outside the body there is no displacement defined, and
    /// continuing the gradient would draw material moving where there is none.
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let (nx, ny, nz) = self.nodes;
        let q = p.to_si() / self.dx;
        // NaN spelled out rather than folded into a comparison: a visualiser can hand one over and it
        // must not reach the cast below.
        if q.is_nan() {
            return 0.0;
        }
        let axis = |v: f64, n: usize| -> (usize, f64) {
            let last = n.saturating_sub(1);
            if v <= 0.0 {
                return (0, 0.0);
            }
            if v >= last as f64 {
                return (last, 0.0);
            }
            let i = v.floor();
            (i as usize, v - i)
        };
        let (i, fx) = axis(q.x, nx);
        let (j, fy) = axis(q.y, ny);
        let (k, fz) = axis(q.z, nz);
        let (i1, j1, k1) = (
            (i + 1).min(nx - 1),
            (j + 1).min(ny - 1),
            (k + 1).min(nz - 1),
        );
        let mag = |a: usize, b: usize, c: usize| {
            let n = 3 * (a + nx * (b + ny * c));
            (self.u[n] * self.u[n] + self.u[n + 1] * self.u[n + 1] + self.u[n + 2] * self.u[n + 2])
                .sqrt()
        };
        let lerp = |lo: f64, hi: f64, t: f64| lo * (1.0 - t) + hi * t;
        let z0 = lerp(
            lerp(mag(i, j, k), mag(i1, j, k), fx),
            lerp(mag(i, j1, k), mag(i1, j1, k), fx),
            fy,
        );
        let z1 = lerp(
            lerp(mag(i, j, k1), mag(i1, j, k1), fx),
            lerp(mag(i, j1, k1), mag(i1, j1, k1), fx),
            fy,
        );
        lerp(z0, z1, fz)
    }
}
