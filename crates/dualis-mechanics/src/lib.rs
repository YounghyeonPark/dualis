//! dualis-mechanics: motion under force, as a domain on the `dualis-core` kernel.
//!
//! Five systems, chosen for what they prove about the kernel rather than for
//! completeness:
//!
//! - [`NBody`] is conservative and has closed-form answers. Gravity between point
//!   masses gives Kepler's third law to check against, and pairwise forces that are
//!   equal and opposite give **exactly** conserved momentum — which makes it the
//!   second customer for [`audit`](dualis_core::audit) after energy, and the test
//!   that the ledger convention was not built around one quantity.
//! - [`ContactSystem`] is not conservative. A penalty contact dissipates through its
//!   dashpot, and that energy becomes heat published on the [`Exchange`] — the same
//!   channel `dualis-thermal` already consumes from `dualis-optics`. Neither
//!   mechanics nor thermal names the other, which is what says the seam is general
//!   and not something optics arranged.
//!
//! # Momentum is exact; energy is not, and that distinction is the interesting one
//!
//! Newton's third law is enforced by construction: forces are accumulated once per
//! pair, added to one body and subtracted from the other, so their sum is zero to
//! the last bit. Momentum is therefore conserved to within the one division and one
//! multiplication that turn a force into an acceleration and back — about 1e-16
//! relative, and it stays there for as long as the simulation runs.
//!
//! Energy is different, and no explicit integrator fixes it. A symplectic method
//! conserves a *nearby* Hamiltonian rather than the true one, so the real energy
//! oscillates by `O(dt²)` forever instead of drifting. That is the best available
//! behaviour and it is still not exact. So [`NBody`] reports momentum in its ledger
//! and exposes energy as a method instead: [`audit`](dualis_core::audit) applies one
//! tolerance to every quantity it sees, and a tolerance loose enough for an
//! integrated energy would blunt the momentum check that is genuinely sharp.

// Every public item carries a doc comment. Denied rather than warned: a public physics API
// whose `Length::mm` shows a blank summary in rustdoc is documented in the sense that a
// paragraph exists somewhere, and not in the sense a reader needs.
#![deny(missing_docs)]
pub mod collision;
pub mod rigid;
pub mod tree;

pub use collision::{resolve, Rolling, Sphere, Surface};
pub use rigid::{Inertia, RigidBody};
pub use tree::TreeNBody;

use dualis_core::conserved::quantity;
use dualis_core::{velocity_verlet, Domain, Exchange, Kind, Ledger, Newtonian, State, Violation};
use dualis_units::{
    Damping, Energy, Length, LengthVec, Mass, Momentum, Qty, Stiffness, Time, Velocity, VelocityVec,
};
use glam::DVec3;

/// Newton's gravitational constant, m³·kg⁻¹·s⁻².
pub const GRAVITATION: Qty<3, -1, -2, 0, 0, 0, 0> = Qty::from_si(6.674_30e-11);

/// The bus channel dissipated mechanical energy leaves on, in joules.
///
/// Deliberately the same name `dualis-thermal` consumes: friction produces heat, and
/// heat is heat. That the two crates agree on a string and nothing else is the whole
/// interface.
pub const HEAT: &str = quantity::ENERGY;

/// Momentum has three components and a [`Ledger`] entry holds one number, so it is
/// audited per axis. The magnitude would not do: a ledger sums its contributions, and
/// the sum of magnitudes is not the magnitude of the sum.
pub mod conserved {
    /// The `x` component of linear momentum, kg·m·s⁻¹.
    pub const MOMENTUM_X: &str = "momentum_x";
    /// The `y` component of linear momentum, kg·m·s⁻¹.
    pub const MOMENTUM_Y: &str = "momentum_y";
    /// The `z` component of linear momentum, kg·m·s⁻¹.
    pub const MOMENTUM_Z: &str = "momentum_z";
}

/// A point mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Body {
    /// How much of it there is.
    pub mass: Mass,
    /// Where it is.
    pub position: LengthVec,
    /// How fast and which way it is going.
    pub velocity: VelocityVec,
}

impl Body {
    /// A point mass with a position and a velocity.
    pub fn new(mass: Mass, position: LengthVec, velocity: VelocityVec) -> Body {
        Body {
            mass,
            position,
            velocity,
        }
    }

    /// ½mv².
    pub fn kinetic_energy(&self) -> Energy {
        dualis_units::vector::kinetic_energy(self.mass, self.velocity)
    }

    /// `mv`. Audited per axis rather than as a magnitude — see [`conserved`].
    pub fn momentum(&self) -> dualis_units::MomentumVec {
        self.mass * self.velocity
    }
}

/// Coordinates of a whole system, as the kernel's integrators want them.
///
/// A newtype over `Vec<DVec3>` in SI metres. The dimension is checked at [`NBody`]'s
/// boundary and dropped here, for the same reason the ray intersection kernels work
/// in raw metres: an integrator does arithmetic on a state vector and cannot know
/// what any of its components mean.
#[derive(Clone, Debug, PartialEq)]
pub struct Coords(pub Vec<DVec3>);

impl State for Coords {
    fn axpy(&mut self, a: f64, other: &Self) {
        // Fixed order, one element at a time — see the kernel's note on why a
        // reduction that could be farmed out is not one.
        for (mine, theirs) in self.0.iter_mut().zip(other.0.iter()) {
            *mine += *theirs * a;
        }
    }

    fn scale(&mut self, a: f64) {
        for v in self.0.iter_mut() {
            *v *= a;
        }
    }

    fn zeros_like(&self) -> Self {
        Coords(vec![DVec3::ZERO; self.0.len()])
    }
}

/// Point masses attracting each other under gravity.
///
/// Conservative, so it is stepped with [`velocity_verlet`] rather than Runge-Kutta:
/// over a long run the symplectic method's energy error stays bounded while RK4's
/// drains the system, which the kernel's integrator tests demonstrate on a harmonic
/// oscillator and which matters far more here, where a drifting orbit is the answer.
pub struct NBody {
    name: String,
    masses: Vec<f64>,
    positions: Coords,
    velocities: Coords,
    /// Softening length, to keep a close pass from producing an infinite force.
    ///
    /// Not a fudge factor to be tuned away: a point mass really does have a singular
    /// potential, and any simulation with a finite step will eventually place two
    /// bodies close enough for it to matter. Zero is allowed, and honest, for systems
    /// where that cannot happen.
    softening: f64,
    saved: Option<(Coords, Coords)>,
}

impl NBody {
    /// Every body attracting every other, summed exactly.
    pub fn new(name: impl Into<String>, bodies: &[Body]) -> NBody {
        NBody {
            name: name.into(),
            masses: bodies.iter().map(|b| b.mass.to_si()).collect(),
            positions: Coords(bodies.iter().map(|b| b.position.to_si()).collect()),
            velocities: Coords(bodies.iter().map(|b| b.velocity.to_si()).collect()),
            softening: 0.0,
            saved: None,
        }
    }

    /// Soften the singularity: `1/(r² + ε²)` instead of `1/r²`.
    ///
    /// Two point masses passing arbitrarily close exchange arbitrarily large momentum in one
    /// step, and no fixed-step integrator survives that. Softening trades the near field for
    /// the far field staying right. Zero is exact Newtonian gravity, and correct only while
    /// nothing gets close.
    pub fn with_softening(mut self, softening: Length) -> NBody {
        self.softening = softening.to_si().abs();
        self
    }

    /// How many bodies.
    pub fn count(&self) -> usize {
        self.masses.len()
    }

    /// One body, reassembled from the separate arrays the integrator works on.
    pub fn body(&self, index: usize) -> Body {
        Body {
            mass: Mass::from_si(self.masses[index]),
            position: LengthVec::from_si(self.positions.0[index]),
            velocity: VelocityVec::from_si(self.velocities.0[index]),
        }
    }

    /// Total mass.
    pub fn total_mass(&self) -> Mass {
        Mass::from_si(self.masses.iter().sum())
    }

    /// Centre of mass. It moves at a constant velocity, always, and a bug in the
    /// force symmetry shows up here first.
    pub fn barycentre(&self) -> LengthVec {
        let total: f64 = self.masses.iter().sum();
        if total <= 0.0 {
            return LengthVec::ZERO;
        }
        let weighted: DVec3 = self
            .masses
            .iter()
            .zip(self.positions.0.iter())
            .map(|(m, p)| *p * *m)
            .fold(DVec3::ZERO, |a, b| a + b);
        LengthVec::from_si(weighted / total)
    }

    /// Total momentum, summed in index order.
    pub fn momentum(&self) -> dualis_units::MomentumVec {
        let p = self
            .masses
            .iter()
            .zip(self.velocities.0.iter())
            .map(|(m, v)| *v * *m)
            .fold(DVec3::ZERO, |a, b| a + b);
        dualis_units::MomentumVec::from_si(p)
    }

    /// Total angular momentum about the origin.
    pub fn angular_momentum(&self) -> DVec3 {
        self.masses
            .iter()
            .zip(self.positions.0.iter().zip(self.velocities.0.iter()))
            .map(|(m, (r, v))| r.cross(*v * *m))
            .fold(DVec3::ZERO, |a, b| a + b)
    }

    /// Σ½mv².
    pub fn kinetic_energy(&self) -> Energy {
        let k: f64 = self
            .masses
            .iter()
            .zip(self.velocities.0.iter())
            .map(|(m, v)| 0.5 * m * v.length_squared())
            .sum();
        Energy::from_si(k)
    }

    /// `-Σ G mᵢmⱼ / rᵢⱼ` over pairs. Negative, because gravity is a bound state.
    pub fn potential_energy(&self) -> Energy {
        let g = GRAVITATION.to_si();
        let mut u = 0.0;
        for i in 0..self.masses.len() {
            for j in (i + 1)..self.masses.len() {
                let d = self.positions.0[j] - self.positions.0[i];
                let r = (d.length_squared() + self.softening * self.softening).sqrt();
                if r > 0.0 {
                    u -= g * self.masses[i] * self.masses[j] / r;
                }
            }
        }
        Energy::from_si(u)
    }

    /// Kinetic plus potential. Bounded rather than constant under a symplectic
    /// integrator — see the module docs.
    pub fn energy(&self) -> Energy {
        self.kinetic_energy() + self.potential_energy()
    }

    /// Circular orbit speed at a radius about a central mass: `√(GM/r)`.
    ///
    /// A closed form, and the thing an orbit test needs to be built from rather than
    /// tuned towards.
    pub fn circular_speed(central: Mass, radius: Length) -> Velocity {
        Velocity::from_si((GRAVITATION.to_si() * central.to_si() / radius.to_si()).sqrt())
    }

    /// Orbital period at a radius about a central mass: `2π√(r³/GM)`.
    ///
    /// Kepler's third law. Independent of the orbiting body's mass, which is the part
    /// that is not obvious and the part worth testing.
    pub fn orbital_period(central: Mass, radius: Length) -> Time {
        let r = radius.to_si();
        Time::from_si(
            std::f64::consts::TAU * (r * r * r / (GRAVITATION.to_si() * central.to_si())).sqrt(),
        )
    }
}

impl Newtonian for NBody {
    type Coords = Coords;

    /// Pairwise gravity, accumulated once per pair.
    ///
    /// The loop runs `i < j` and applies `+f` to `i` and `-f` to `j` from the same
    /// computed vector, so the total force is zero to the last bit rather than to
    /// within a summation error. That is Newton's third law held structurally, and it
    /// is why the momentum audit is sharp enough to be worth running.
    fn acceleration(&self, x: &Coords, _t: Time) -> Coords {
        let n = self.masses.len();
        let g = GRAVITATION.to_si();
        let soft2 = self.softening * self.softening;
        let mut force = vec![DVec3::ZERO; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let d = x.0[j] - x.0[i];
                let r2 = d.length_squared() + soft2;
                if r2 <= 0.0 {
                    continue;
                }
                let inv_r = r2.sqrt().recip();
                // G m_i m_j / r^2 along the unit separation.
                let f = d * (g * self.masses[i] * self.masses[j] * inv_r * inv_r * inv_r);
                force[i] += f;
                force[j] -= f;
            }
        }
        for (a, m) in force.iter_mut().zip(self.masses.iter()) {
            *a /= *m;
        }
        Coords(force)
    }
}

impl Domain for NBody {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// A twentieth of the shortest orbital timescale present, estimated from the
    /// closest pair and the largest mass.
    ///
    /// Gravity has no CFL condition, but it does have a shortest period, and stepping
    /// past a fraction of it turns a bound orbit into an escape. This is an estimate
    /// rather than a bound, and it is stated as such: a genuinely close pass needs a
    /// smaller step than any formula computed at the start of a window can know.
    fn max_stable_dt(&self, _now: Time) -> Time {
        let g = GRAVITATION.to_si();
        let mut fastest = 0.0f64;
        for i in 0..self.masses.len() {
            for j in (i + 1)..self.masses.len() {
                let r = (self.positions.0[j] - self.positions.0[i]).length();
                let r = (r * r + self.softening * self.softening).sqrt();
                if r <= 0.0 {
                    continue;
                }
                let m = self.masses[i] + self.masses[j];
                // Angular frequency of the two-body orbit at this separation.
                fastest = fastest.max((g * m / (r * r * r)).sqrt());
            }
        }
        if fastest <= 0.0 {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(1.0 / (20.0 * fastest))
    }

    fn step(&mut self, t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let mut x = std::mem::replace(&mut self.positions, Coords(Vec::new()));
        let mut v = std::mem::replace(&mut self.velocities, Coords(Vec::new()));
        velocity_verlet(self, &mut x, &mut v, t, dt);
        self.positions = x;
        self.velocities = v;
        Ok(())
    }

    /// Momentum, per axis. Energy is deliberately absent — see the module docs.
    /// Momentum per axis, accumulated **one body at a time**.
    ///
    /// The per-body `add` is the whole point and not a style choice. `Ledger` records the
    /// largest single contribution as an entry's `scale`, and [`audit`](dualis_core::audit)
    /// judges a change against that rather than against the total — because a correct
    /// system's total momentum is usually exactly zero, and a relative tolerance on zero
    /// means nothing.
    ///
    /// Handing the pre-summed total to one `with` call sets the scale to `|total|`, which for
    /// a symmetric system is `0.0`. `audit` then skips the entry entirely at its
    /// `scale < 1e-300` guard, and the momentum audit silently does not run. That is what
    /// this did until an audit set the tolerance to `0.0` — a setting that must reject any
    /// change at all — and the test still passed.
    fn ledger(&self) -> Ledger {
        let mut ledger = Ledger::new();
        for (m, v) in self.masses.iter().zip(self.velocities.0.iter()) {
            let p = *v * *m;
            ledger.add(conserved::MOMENTUM_X, p.x);
            ledger.add(conserved::MOMENTUM_Y, p.y);
            ledger.add(conserved::MOMENTUM_Z, p.z);
        }
        ledger
    }

    fn checkpoint(&mut self) {
        self.saved = Some((self.positions.clone(), self.velocities.clone()));
    }

    fn restore(&mut self) {
        if let Some((x, v)) = self.saved.clone() {
            self.positions = x;
            self.velocities = v;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    /// Opted in so a caller can read the bodies back out — a test asserting an orbit, a
    /// renderer drawing one. Every other domain with state to show already does this; these
    /// were simply never asked, which is what a library with no consumer looks like.
    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }
}

/// A flat, immovable floor at a height along a normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ground {
    /// Outward unit normal, pointing away from the solid side.
    pub normal: DVec3,
    /// Signed offset along the normal: a point `p` is penetrating when
    /// `p·normal < offset`.
    pub offset: Length,
}

impl Ground {
    /// The plane `z = 0` with the solid side below.
    pub fn floor() -> Ground {
        Ground {
            normal: DVec3::Z,
            offset: Length::ZERO,
        }
    }

    /// How far a point has sunk in, or zero if it has not.
    pub fn penetration(&self, p: LengthVec) -> Length {
        let depth = self.offset - p.along(self.normal);
        if depth.to_si() > 0.0 {
            depth
        } else {
            Length::ZERO
        }
    }
}

/// Particles under uniform gravity, resting on a penalty contact.
///
/// The contact is a spring and a dashpot: while a particle is penetrating, it is
/// pushed out in proportion to how far it has sunk and damped in proportion to how
/// fast it is moving. That is the simplest contact that is not a discontinuity, and
/// it has two properties worth having:
///
/// - It has a **real stability limit**. The contact spring's period is `2π√(m/k)`
///   whether or not anybody is interested in it, and an explicit integrator must
///   resolve it. Stiffer contact is more realistic and strictly more expensive, with
///   no way around the trade — which is what
///   [`Schedule::Multirate`](dualis_core::Schedule::Multirate) is for.
/// - It **dissipates**, and says where the energy went. The dashpot's work is
///   published as [`HEAT`], so the books close and a thermal domain can pick it up.
pub struct ContactSystem {
    name: String,
    masses: Vec<f64>,
    positions: Coords,
    velocities: Coords,
    gravity: DVec3,
    ground: Ground,
    stiffness: f64,
    damping: f64,
    /// Joules taken out by the dashpot and handed to the bus.
    dissipated: f64,
    saved: Option<(Coords, Coords, f64)>,
}

impl ContactSystem {
    /// Bodies falling under gravity onto a ground plane, meeting it through a penalty
    /// contact.
    pub fn new(
        name: impl Into<String>,
        bodies: &[Body],
        gravity: dualis_units::AccelerationVec,
        ground: Ground,
        stiffness: Stiffness,
        damping: Damping,
    ) -> ContactSystem {
        ContactSystem {
            name: name.into(),
            masses: bodies.iter().map(|b| b.mass.to_si()).collect(),
            positions: Coords(bodies.iter().map(|b| b.position.to_si()).collect()),
            velocities: Coords(bodies.iter().map(|b| b.velocity.to_si()).collect()),
            gravity: gravity.to_si(),
            ground,
            stiffness: stiffness.to_si(),
            damping: damping.to_si(),
            dissipated: 0.0,
            saved: None,
        }
    }

    /// One body.
    pub fn body(&self, index: usize) -> Body {
        Body {
            mass: Mass::from_si(self.masses[index]),
            position: LengthVec::from_si(self.positions.0[index]),
            velocity: VelocityVec::from_si(self.velocities.0[index]),
        }
    }

    /// How many bodies.
    pub fn count(&self) -> usize {
        self.masses.len()
    }

    /// Height of a body above the ground plane.
    pub fn height(&self, index: usize) -> Length {
        LengthVec::from_si(self.positions.0[index]).along(self.ground.normal) - self.ground.offset
    }

    /// Kinetic plus gravitational plus the energy stored in whatever contact springs
    /// are currently compressed.
    pub fn mechanical_energy(&self) -> Energy {
        let mut e = 0.0;
        for (i, m) in self.masses.iter().enumerate() {
            let (p, v) = (self.positions.0[i], self.velocities.0[i]);
            e += 0.5 * m * v.length_squared();
            // Gravitational potential, measured from the ground plane.
            e -= m * self.gravity.dot(p);
            let depth = self.ground.penetration(LengthVec::from_si(p)).to_si();
            e += 0.5 * self.stiffness * depth * depth;
        }
        Energy::from_si(e)
    }

    /// Joules the dashpot has removed over the run.
    pub fn dissipated_energy(&self) -> Energy {
        Energy::from_si(self.dissipated)
    }

    /// Period of the contact spring for the lightest body: `2π√(m/k)`.
    pub fn contact_period(&self) -> Time {
        let lightest = self.lightest_mass();
        if self.stiffness <= 0.0 || !lightest.is_finite() {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(std::f64::consts::TAU * (lightest / self.stiffness).sqrt())
    }

    fn lightest_mass(&self) -> f64 {
        self.masses.iter().cloned().fold(f64::MAX, f64::min)
    }

    /// Damping ratio of the contact oscillator, `c / (2√(km))`.
    ///
    /// Under 1 the contact bounces; at 1 it is critically damped and stops dead; above
    /// 1 it is overdamped and settles slowly without bouncing.
    pub fn damping_ratio(&self) -> f64 {
        let m = self.lightest_mass();
        if self.stiffness <= 0.0 || !m.is_finite() || m <= 0.0 {
            return 0.0;
        }
        self.damping / (2.0 * (self.stiffness * m).sqrt())
    }

    /// Coefficient of restitution this contact implies: `exp(-ζπ/√(1-ζ²))`.
    ///
    /// How bouncy a penalty contact is **is not a free parameter**. It follows from the
    /// stiffness, the dashpot and the mass, because the contact is a damped oscillator
    /// and a bounce is half of one cycle of it. So "make it bounce 60%" is really
    /// "choose a damping ratio of 0.16", and picking `k` and `c` independently means
    /// picking a restitution without noticing.
    ///
    /// This is a *velocity* ratio; a drop returns to `restitution²` of its height.
    /// Zero for a critically damped or overdamped contact, which does not bounce.
    pub fn restitution(&self) -> f64 {
        let zeta = self.damping_ratio();
        if zeta >= 1.0 {
            return 0.0;
        }
        (-zeta * std::f64::consts::PI / (1.0 - zeta * zeta).sqrt()).exp()
    }

    /// Force on body `i` at a given position and velocity, in SI newtons.
    fn force_on(&self, i: usize, p: DVec3, v: DVec3) -> DVec3 {
        let mut f = self.gravity * self.masses[i];
        let depth = self.ground.penetration(LengthVec::from_si(p)).to_si();
        if depth > 0.0 {
            // Spring out along the normal, plus a dashpot on the normal velocity.
            let closing = -v.dot(self.ground.normal);
            f += self.ground.normal * (self.stiffness * depth + self.damping * closing);
        }
        f
    }
}

impl Domain for ContactSystem {
    fn name(&self) -> &str {
        &self.name
    }

    /// A hundredth of the contact spring's period.
    ///
    /// Unlike [`NBody::max_stable_dt`] the period this is derived from is a genuine
    /// bound rather than an estimate: the stiffness and the masses are known up front,
    /// and what they imply does not depend on where anything is.
    ///
    /// The factor of a hundred is an **accuracy** choice, not a stability one, and it
    /// was measured rather than guessed. Bare stability needs only a handful of steps
    /// per period; at a twentieth, a contact lasts about ten steps and the energy books
    /// close only to about 7% over a second of bouncing. At a hundredth the contact is
    /// resolved by fifty steps and that falls to about 2%. This follows the same
    /// reasoning as `LumpedMass`, which reports a tenth of its time constant rather
    /// than the twice-that stability would allow: a scheduler honouring the limit
    /// should get an answer, not merely a finite one.
    ///
    /// **Why it is only 2% and not better.** Semi-implicit Euler is symplectic, and a
    /// symplectic method's energy error is supposed to stay bounded rather than drift.
    /// That guarantee needs a smooth potential, and a penalty contact is not one: it
    /// switches on the moment a body touches. Each transition shifts the shadow
    /// Hamiltonian the method actually conserves, so the true energy takes an `O(h)`
    /// step at every bounce and those steps accumulate. Resolving the contact more
    /// finely shrinks each one, which is why the factor is worth what it costs, but no
    /// step size recovers the bound. Getting past that needs an integrator that knows
    /// where the transition is.
    fn max_stable_dt(&self, _now: Time) -> Time {
        self.contact_period() / 100.0
    }

    /// Semi-implicit Euler, which is what a contact solver wants.
    ///
    /// Not velocity-Verlet: the dashpot force depends on velocity, so the system is
    /// not [`Newtonian`] in the kernel's sense and there is no energy for a symplectic
    /// method to preserve. Updating the velocity from the force and *then* the position
    /// from the new velocity is stable where fully explicit Euler is not, at first
    /// order and no extra force evaluation.
    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let h = dt.to_si();
        if h <= 0.0 {
            return Ok(());
        }
        let mut dissipated_now = 0.0;
        for i in 0..self.masses.len() {
            let (p, v) = (self.positions.0[i], self.velocities.0[i]);
            let f = self.force_on(i, p, v);
            let v_next = v + f * (h / self.masses[i]);
            let p_next = p + v_next * h;

            // Work done by the dashpot alone, over the distance actually travelled.
            // Positive, always: a dashpot cannot add energy, and computing it from the
            // force rather than from an energy difference is what keeps the books
            // closable when the spring is also doing work.
            let depth = self.ground.penetration(LengthVec::from_si(p)).to_si();
            if depth > 0.0 {
                let closing = -v.dot(self.ground.normal);
                let damper = self.damping * closing;
                // -F_damp · dx, with the force along the normal and dx = v_next * h.
                let work = -damper * v_next.dot(self.ground.normal) * h;
                dissipated_now += work;
            }

            self.velocities.0[i] = v_next;
            self.positions.0[i] = p_next;
        }
        self.dissipated += dissipated_now;
        bus.publish(HEAT, dissipated_now);
        Ok(())
    }

    /// Mechanical energy still in the system, and **only** that.
    ///
    /// Not plus the dissipated total. Those joules were published and taken, so
    /// whatever took them is holding them and reporting them; counting them here as
    /// well makes the simulation's total grow by the heat it produced. That is the
    /// mirror image of the mistake `dualis-thermal` made in the other direction — a
    /// consumer that subtracts what it absorbed as well as reporting what it stored —
    /// and the rule both come down to is the same: **a ledger says what you are
    /// holding, not what has passed through you.**
    ///
    /// [`Exchange::audit_transfers`](dualis_core::Exchange::audit_transfers) is what
    /// makes this safe: joules are only ever briefly nobody's, between a publish and a
    /// take inside one sweep, and the bus refuses to end a step with any left on it.
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.mechanical_energy().to_si())
    }

    fn checkpoint(&mut self) {
        self.saved = Some((
            self.positions.clone(),
            self.velocities.clone(),
            self.dissipated,
        ));
    }

    fn restore(&mut self) {
        if let Some((x, v, d)) = self.saved.clone() {
            self.positions = x;
            self.velocities = v;
            self.dissipated = d;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    /// Opted in so a coupling test can ask what mechanical energy is *left*.
    ///
    /// Without it the only reachable number is the one the system started with, and a test
    /// that wants to say "the joules that crossed the bus are the joules that went missing"
    /// has to weaken itself to an inequality against the starting energy. See
    /// `crates/dualis/tests/friction_heats.rs`, which made exactly that compromise.
    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }
}

/// Momentum magnitude, for reporting.
pub fn momentum_magnitude(p: dualis_units::MomentumVec) -> Momentum {
    p.length()
}

#[cfg(test)]
mod tests {
    use super::*;
    use dualis_core::{audit, Schedule, Simulation};
    use dualis_units::{AccelerationVec, G0};

    /// Earth's mass and a low orbit, so the numbers can be checked against the ones
    /// everybody knows.
    fn earth() -> Mass {
        Mass::from_si(5.972_2e24)
    }

    fn orbit_radius() -> Length {
        // 400 km up, which is roughly where the space station is.
        Length::from_si(6_371e3 + 400e3)
    }

    /// Kepler's third law, against the closed form, and the part of it that is not
    /// obvious: the period does not depend on the orbiting mass at all. A 1 kg
    /// satellite and a 1000 kg one at the same radius take the same time.
    #[test]
    fn keplers_third_law_holds_and_ignores_the_orbiting_mass() {
        let m = earth();
        // Low orbit is about 92 minutes.
        let low = NBody::orbital_period(m, orbit_radius());
        assert!(
            (low.to_si() / 60.0 - 92.6).abs() < 0.5,
            "low orbit period {} min",
            low.to_si() / 60.0
        );
        // Geostationary radius gives a sidereal day.
        let geo = NBody::orbital_period(m, Length::from_si(42_164e3));
        assert!(
            (geo.to_si() / 3600.0 - 23.93).abs() < 0.05,
            "geostationary period {} h",
            geo.to_si() / 3600.0
        );
        // T^2 goes as r^3.
        let ratio_t = (geo / low).powi(2);
        let ratio_r = (Length::from_si(42_164e3) / orbit_radius()).powi(3);
        assert!(
            (ratio_t / ratio_r - 1.0).abs() < 1e-9,
            "T^2/r^3 is not constant: {ratio_t} vs {ratio_r}"
        );
        // And the speed at low orbit is 7.7 km/s.
        let v = NBody::circular_speed(m, orbit_radius());
        assert!(
            (v.to_si() / 1e3 - 7.67).abs() < 0.02,
            "{} km/s",
            v.to_si() / 1e3
        );
    }

    /// A circular orbit integrated for one period comes back to where it started, and
    /// stays at its radius throughout. Built from the closed-form speed rather than
    /// tuned until it looked round.
    #[test]
    fn a_circular_orbit_stays_circular() {
        let m = earth();
        let r = orbit_radius();
        let v = NBody::circular_speed(m, r);
        let period = NBody::orbital_period(m, r);

        // The radius error is second order in the step, so with N steps per orbit it
        // sits near `(2pi/N)^2` — 9.9e-6 at N = 2000, and measured at half that.
        // Rather than assert a threshold and hope, run two step counts and check the
        // error actually falls as the square.
        let worst_radius_error = |steps: u32| {
            // The central mass is held still by being enormous, not by being fixed.
            let mut system = NBody::new(
                "orbit",
                &[
                    Body::new(m, LengthVec::ZERO, VelocityVec::ZERO),
                    Body::new(
                        Mass::kg(1000.0),
                        LengthVec::from_si(DVec3::X * r.to_si()),
                        VelocityVec::from_si(DVec3::Y * v.to_si()),
                    ),
                ],
            );
            let dt = period / steps as f64;
            let mut bus = Exchange::new();
            let mut worst = 0.0f64;
            let mut t = Time::ZERO;
            for _ in 0..steps {
                system.step(t, dt, &mut bus).unwrap();
                t += dt;
                let radius = (system.body(1).position - system.body(0).position).length();
                worst = worst.max((radius / r - 1.0).abs());
            }
            (worst, system)
        };

        let (coarse, _) = worst_radius_error(1000);
        let (fine, system) = worst_radius_error(2000);
        assert!(fine < 1e-5, "the orbit should stay round: {fine:e}");
        assert!(
            (coarse / fine - 4.0).abs() < 0.4,
            "halving the step should quarter the radius error, got {}",
            coarse / fine
        );

        // And after one period it is back where it began.
        let closed = (system.body(1).position - LengthVec::from_si(DVec3::X * r.to_si())).length();
        assert!(
            closed / r < 1e-4,
            "one period should close the orbit: off by {:e} of a radius",
            closed / r
        );
    }

    /// Momentum is exact, and it is exact because Newton's third law is structural
    /// rather than emergent. This is the claim the ledger convention was generalised
    /// for: a second conserved quantity, audited by the same machinery as energy.
    #[test]
    fn momentum_is_conserved_to_machine_precision() {
        let mut system = NBody::new(
            "three-body",
            &[
                Body::new(
                    Mass::kg(3e10),
                    LengthVec::m(0.0, 0.0, 0.0),
                    VelocityVec::m_per_s(0.0, 0.1, 0.0),
                ),
                Body::new(
                    Mass::kg(2e10),
                    LengthVec::m(100.0, 0.0, 0.0),
                    VelocityVec::m_per_s(0.0, -0.15, 0.02),
                ),
                Body::new(
                    Mass::kg(1e10),
                    LengthVec::m(30.0, 80.0, -10.0),
                    VelocityVec::m_per_s(-0.05, 0.0, 0.0),
                ),
            ],
        );
        let before = system.momentum().to_si();
        let scale = before.length();
        assert!(scale > 0.0, "the test needs a nonzero momentum to check");
        // The centre of mass does not start at the origin — it starts wherever these
        // three masses put it — so the straight line it travels along has to be
        // measured from there.
        let start_barycentre = system.barycentre();

        let mut bus = Exchange::new();
        let dt = system.max_stable_dt(Time::ZERO);
        let mut t = Time::ZERO;
        for _ in 0..5000 {
            system.step(t, dt, &mut bus).unwrap();
            t += dt;
        }
        let after = system.momentum().to_si();
        let error = (after - before).length() / scale;
        // 6e-13 over five thousand steps, and two things make that figure pessimistic
        // rather than loose. The forces cancel exactly, but turning each into an
        // acceleration and back multiplies and divides by a mass, and those bits do
        // not come back — so the drift is a random walk at about one ulp a step. And
        // the bodies here move partly against each other, so the net momentum is
        // several times smaller than any individual one, which inflates every relative
        // figure measured against it.
        assert!(
            error < 1e-11,
            "momentum drifted by {error:e} of itself over 5000 steps"
        );

        // Which means the barycentre travels in a straight line at a constant speed:
        // exactly `x0 + (p/M) t`, with no gravitational term surviving at all.
        let expected =
            start_barycentre + LengthVec::from_si(before / system.total_mass().to_si() * t.to_si());
        let travelled = (expected - start_barycentre).length();
        assert!(travelled.to_si() > 1.0, "it should have gone somewhere");
        let drift = (system.barycentre() - expected).length();
        assert!(
            drift / travelled < 1e-11,
            "the barycentre wandered {drift:?} off a {travelled:?} straight line"
        );
    }

    /// The audit sees momentum as a conserved quantity in its own right, and would
    /// fire if the pairwise forces stopped cancelling.
    #[test]
    fn the_kernel_audits_momentum_the_way_it_audits_energy() {
        let mut system = NBody::new(
            "pair",
            &[
                Body::new(
                    Mass::kg(1e12),
                    LengthVec::m(0.0, 0.0, 0.0),
                    VelocityVec::m_per_s(0.0, 0.5, 0.0),
                ),
                Body::new(
                    Mass::kg(1e12),
                    LengthVec::m(1000.0, 0.0, 0.0),
                    VelocityVec::m_per_s(0.0, -0.3, 0.0),
                ),
            ],
        );
        let before = system.ledger();
        assert_eq!(before.get(conserved::MOMENTUM_Y), Some(2e11));

        let mut bus = Exchange::new();
        let dt = system.max_stable_dt(Time::ZERO);
        for _ in 0..200 {
            system.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        audit("nbody", &before, &system.ledger(), 1e-12)
            .expect("pairwise forces must cancel exactly");

        // A ledger with a deliberately broken momentum is caught, so the check above
        // is not passing for want of teeth.
        let broken = Ledger::new()
            .with(conserved::MOMENTUM_X, 0.0)
            .with(conserved::MOMENTUM_Y, 2e11 * 1.001)
            .with(conserved::MOMENTUM_Z, 0.0);
        let err = audit("nbody", &before, &broken, 1e-12).expect_err("0.1% is not rounding");
        assert_eq!(err.quantity, "momentum_y");
    }

    /// Angular momentum is conserved too, for a central force — and that is what makes
    /// an orbit a plane curve rather than a space curve.
    #[test]
    fn angular_momentum_is_conserved() {
        let m = earth();
        let r = orbit_radius();
        let v = NBody::circular_speed(m, r);
        let mut system = NBody::new(
            "orbit",
            &[
                Body::new(m, LengthVec::ZERO, VelocityVec::ZERO),
                Body::new(
                    Mass::kg(1000.0),
                    LengthVec::from_si(DVec3::X * r.to_si()),
                    VelocityVec::from_si(DVec3::Y * v.to_si()),
                ),
            ],
        );
        let before = system.angular_momentum();
        let mut bus = Exchange::new();
        let dt = NBody::orbital_period(m, r) / 500.0;
        for _ in 0..500 {
            system.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        let after = system.angular_momentum();
        assert!(
            (after - before).length() / before.length() < 1e-12,
            "angular momentum drifted"
        );
        // It points along +z and the motion stays in the xy plane.
        assert!(before.normalize().dot(DVec3::Z) > 1.0 - 1e-12);
        assert!(system.body(1).position.along(DVec3::Z).to_si().abs() < 1e-6);
    }

    /// A circular orbit's total energy is exactly `-GMm/2r`, by the virial theorem —
    /// kinetic is half the magnitude of potential and opposite in sign. A closed form
    /// worth checking because it ties the two energy terms together.
    #[test]
    fn a_circular_orbit_satisfies_the_virial_theorem() {
        let m = earth();
        let r = orbit_radius();
        let satellite = Mass::kg(1000.0);
        let v = NBody::circular_speed(m, r);
        let system = NBody::new(
            "orbit",
            &[
                Body::new(m, LengthVec::ZERO, VelocityVec::ZERO),
                Body::new(
                    satellite,
                    LengthVec::from_si(DVec3::X * r.to_si()),
                    VelocityVec::from_si(DVec3::Y * v.to_si()),
                ),
            ],
        );
        let expected = -GRAVITATION.to_si() * m.to_si() * satellite.to_si() / (2.0 * r.to_si());
        assert!(
            (system.energy().to_si() / expected - 1.0).abs() < 1e-12,
            "total energy {:e}, virial {:e}",
            system.energy().to_si(),
            expected
        );
        // 2K = -U.
        let k = system.kinetic_energy().to_si();
        let u = system.potential_energy().to_si();
        assert!(
            (2.0 * k / -u - 1.0).abs() < 1e-12,
            "2K/-U = {}",
            2.0 * k / -u
        );
    }

    /// Softening is a real choice, not a tuning knob: with none, a head-on pass drives
    /// the force to infinity and the integration leaves the domain of arithmetic.
    #[test]
    fn softening_keeps_a_close_pass_finite() {
        let bodies = [
            Body::new(
                Mass::kg(1e14),
                LengthVec::m(-10.0, 0.0, 0.0),
                VelocityVec::m_per_s(3.0, 0.0, 0.0),
            ),
            Body::new(
                Mass::kg(1e14),
                LengthVec::m(10.0, 0.0, 0.0),
                VelocityVec::m_per_s(-3.0, 0.0, 0.0),
            ),
        ];
        let mut softened = NBody::new("softened", &bodies).with_softening(Length::m(1.0));
        let mut bus = Exchange::new();
        let dt = Time::s(0.05);
        for _ in 0..400 {
            softened.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        assert!(
            softened.body(0).position.to_si().is_finite()
                && softened.body(1).position.to_si().is_finite(),
            "a softened pass stays finite"
        );
        // Momentum is still exactly zero: they started symmetric and stay so.
        assert!(softened.momentum().length().to_si() < 1e-6);
    }

    // -----------------------------------------------------------------------
    // Contact
    // -----------------------------------------------------------------------

    fn ball(mass_kg: f64, height_m: f64) -> Body {
        Body::new(
            Mass::kg(mass_kg),
            LengthVec::m(0.0, 0.0, height_m),
            VelocityVec::ZERO,
        )
    }

    fn dropped(damping: f64) -> ContactSystem {
        ContactSystem::new(
            "ball",
            &[ball(0.1, 1.0)],
            AccelerationVec::from_si(-DVec3::Z * G0.to_si()),
            Ground::floor(),
            Stiffness::from_si(1e5),
            Damping::from_si(damping),
        )
    }

    /// The contact stiffness sets the step, and the relationship is the one the trade
    /// is made on: a hundred times stiffer contact is ten times more expensive,
    /// because the period goes as one over the square root of `k`.
    #[test]
    fn contact_stiffness_sets_the_step_it_costs() {
        let soft = dropped(0.0);
        let stiff = ContactSystem::new(
            "stiff",
            &[ball(0.1, 1.0)],
            AccelerationVec::from_si(-DVec3::Z * G0.to_si()),
            Ground::floor(),
            Stiffness::from_si(1e7),
            Damping::from_si(0.0),
        );
        // 0.1 kg on 1e5 N/m gives sqrt(m/k) = 1 ms exactly, so the period is 2*pi ms
        // and the limit is a hundredth of it.
        assert!(
            (soft.contact_period().in_us() - 6283.2).abs() < 1.0,
            "period {} us",
            soft.contact_period().in_us()
        );
        assert!(
            (soft.max_stable_dt(Time::ZERO).in_us() - 62.83).abs() < 0.1,
            "limit {} us",
            soft.max_stable_dt(Time::ZERO).in_us()
        );
        let ratio = soft.max_stable_dt(Time::ZERO) / stiff.max_stable_dt(Time::ZERO);
        assert!(
            (ratio - 10.0).abs() < 1e-9,
            "100x stiffer is 10x dearer: {ratio}"
        );
        // A contact with no stiffness has no limit, and says so.
        let free = ContactSystem::new(
            "free",
            &[ball(0.1, 1.0)],
            AccelerationVec::ZERO,
            Ground::floor(),
            Stiffness::from_si(0.0),
            Damping::from_si(0.0),
        );
        assert!(!free.max_stable_dt(Time::ZERO).to_si().is_finite());
    }

    /// Free fall against the closed form: `h - ½gt²`, before anything touches.
    #[test]
    fn free_fall_matches_the_closed_form() {
        let mut system = dropped(0.0);
        let mut bus = Exchange::new();
        let dt = Time::ms(0.1);
        let steps = 2000; // 0.2 s, well before impact
        for _ in 0..steps {
            system.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        let elapsed = dt.to_si() * steps as f64;
        let expected = 1.0 - 0.5 * G0.to_si() * elapsed * elapsed;
        let got = system.height(0).to_si();
        // Semi-implicit Euler overshoots the drop by half a step of velocity, which
        // at 0.1 ms is 20 micrometres over 0.2 s.
        assert!(
            (got - expected).abs() < 1e-4,
            "after {elapsed} s: got {got:.6} m, closed form {expected:.6} m"
        );
        // Nothing was dissipated, because nothing has touched yet.
        assert_eq!(system.dissipated_energy().to_si(), 0.0);
        assert_eq!(bus.peek(HEAT), 0.0);
    }

    /// An undamped contact is conservative, so the ball comes back to the height it
    /// was dropped from. That is the sharpest available check on the contact force:
    /// a spring that did the wrong amount of work would return it somewhere else.
    #[test]
    fn an_undamped_bounce_returns_to_its_drop_height() {
        let mut system = dropped(0.0);
        let mut bus = Exchange::new();
        let dt = system.max_stable_dt(Time::ZERO) / 10.0;
        let mut highest_after_bounce = 0.0f64;
        let mut touched = false;
        let mut t = Time::ZERO;
        for _ in 0..400_000 {
            system.step(t, dt, &mut bus).unwrap();
            t += dt;
            let h = system.height(0).to_si();
            if h <= 0.0 {
                touched = true;
            } else if touched {
                highest_after_bounce = highest_after_bounce.max(h);
            }
            // Drop everything the contact offered; an undamped one offers nothing.
            bus.take(HEAT);
        }
        assert!(touched, "it should have reached the floor");
        assert!(
            (highest_after_bounce - 1.0).abs() < 0.02,
            "an undamped bounce should return to 1 m, reached {highest_after_bounce:.4} m"
        );
        assert!(
            system.dissipated_energy().to_si().abs() < 1e-12,
            "nothing dissipates without a dashpot"
        );
    }

    /// The bounce height follows the closed-form restitution of a damped oscillator,
    /// and the energy that went missing is exactly the heat that was published.
    ///
    /// Two damping ratios, because one would not distinguish a correct model from a
    /// coincidence: `c = 2` against a critical 200 N·s/m is barely damped and returns
    /// 94% of the height, while `c = 22` returns half of it. Both figures come from
    /// `exp(-2ζπ/√(1-ζ²))` rather than from having been run first.
    #[test]
    fn a_damped_bounce_matches_its_closed_form_restitution() {
        for damping in [2.0f64, 22.0] {
            let mut system = dropped(damping);
            let expected_height = system.restitution().powi(2);
            let start = system.mechanical_energy().to_si();
            let mut bus = Exchange::new();
            let dt = system.max_stable_dt(Time::ZERO) / 10.0;
            let mut collected = 0.0;
            let mut highest_after_bounce = 0.0f64;
            let mut touched = false;
            let mut t = Time::ZERO;
            for _ in 0..400_000 {
                system.step(t, dt, &mut bus).unwrap();
                t += dt;
                collected += bus.take(HEAT);
                let h = system.height(0).to_si();
                if h <= 0.0 {
                    touched = true;
                } else if touched {
                    highest_after_bounce = highest_after_bounce.max(h);
                }
            }
            assert!(touched, "c = {damping}: it should have reached the floor");
            assert!(
                (highest_after_bounce / expected_height - 1.0).abs() < 0.02,
                "c = {damping}: restitution {:.4} predicts {expected_height:.4} m, \
                 reached {highest_after_bounce:.4} m",
                system.restitution()
            );

            // Everything published was accounted for, and it is what the system lost.
            assert!(collected > 0.0, "a dashpot dissipates");
            assert!(
                (collected - system.dissipated_energy().to_si()).abs() < 1e-12,
                "the bus and the books must agree"
            );
            let now = system.mechanical_energy().to_si();
            assert!(
                ((start - now) / collected - 1.0).abs() < 0.02,
                "c = {damping}: lost {:.6} J mechanically and published {:.6} J of heat",
                start - now,
                collected
            );
        }
    }

    /// Restitution is a consequence of the stiffness, the dashpot and the mass, not a
    /// setting — and the two ends of the range behave as they must.
    #[test]
    fn restitution_follows_from_the_damping_ratio() {
        assert!(
            (dropped(0.0).restitution() - 1.0).abs() < 1e-15,
            "no dashpot, no loss"
        );
        // Critical damping for 0.1 kg on 1e5 N/m is 2*sqrt(k m) = 200 N s/m.
        let critical = dropped(200.0);
        assert!((critical.damping_ratio() - 1.0).abs() < 1e-12);
        assert_eq!(
            critical.restitution(),
            0.0,
            "a critically damped contact sticks"
        );
        assert_eq!(
            dropped(400.0).restitution(),
            0.0,
            "and so does an overdamped one"
        );
        // Monotone in between.
        let mut previous = 1.0;
        for c in [1.0f64, 5.0, 22.0, 60.0, 120.0, 190.0] {
            let e = dropped(c).restitution();
            assert!(e < previous && e > 0.0, "c = {c} gave {e}");
            previous = e;
        }
        assert!((dropped(22.0).damping_ratio() - 0.11).abs() < 1e-12);
    }

    /// The domain under the kernel's scheduler, with the audit switched on.
    ///
    /// The tolerance is looser than the default 1e-9, and that is a property of the
    /// integrator rather than a concession: semi-implicit Euler is first order, so the
    /// mechanical energy it reports is first-order accurate too, and demanding nine
    /// digits of a first-order method would be demanding something no explicit
    /// integrator provides. The way to buy accuracy is a smaller reported limit, which
    /// is exactly the trade `max_stable_dt` documents.
    #[test]
    fn the_domain_balances_its_books_under_the_scheduler() {
        struct Sink {
            taken: f64,
        }
        impl Domain for Sink {
            fn name(&self) -> &str {
                "sink"
            }
            fn step(&mut self, _t: Time, _dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
                self.taken += bus.take(HEAT);
                Ok(())
            }
            fn ledger(&self) -> Ledger {
                // What it took, because it is holding it. The mechanics domain reports
                // only the mechanical energy it still has, so this is the other half of
                // the total and the two together are constant.
                Ledger::new().with(quantity::ENERGY, self.taken)
            }
            fn checkpoint(&mut self) {}
            fn restore(&mut self) {}
            fn supports_restore(&self) -> bool {
                true
            }
        }

        let mut sim = Simulation::new(Schedule::Multirate)
            .conservation_tolerance(1e-2)
            .with(dropped(2.0))
            .with(Sink { taken: 0.0 });

        for _ in 0..40 {
            sim.advance(Time::ms(20.0))
                .expect("mechanical energy plus heat must be conserved");
        }
        assert!((sim.time().to_si() - 0.8).abs() < 1e-9);
        assert!(sim.bus().total_consumed(quantity::ENERGY) > 0.0);
        assert!(sim.bus().unclaimed().next().is_none());
    }

    /// The scheduler subcycles the contact domain, because its stability limit is
    /// hundreds of microseconds and a video frame is not.
    #[test]
    fn contact_subcycles_within_a_frame() {
        let mut sim = Simulation::new(Schedule::Multirate)
            .conservation_tolerance(1e-2)
            .with(dropped(0.0));
        // Nothing consumes, so an undamped contact is the only thing that can run
        // without leaving heat on the bus.
        let report = sim.advance(Time::ms(16.7)).unwrap();
        assert_eq!(report.substeps[0].0, "ball");
        // 16.7 ms at a 62.8 us limit is 266 substeps — a video frame costs that many
        // contact steps whether or not anything is touching, which is what stiff
        // contact costs and why the number is worth seeing.
        assert_eq!(report.substeps[0].1, 266);
    }

    /// Penetration is one-sided: a body above the ground is not attracted to it.
    #[test]
    fn the_ground_only_pushes() {
        let ground = Ground::floor();
        assert_eq!(
            ground.penetration(LengthVec::m(0.0, 0.0, 1.0)),
            Length::ZERO
        );
        assert_eq!(
            ground.penetration(LengthVec::m(0.0, 0.0, 0.0)),
            Length::ZERO
        );
        assert!((ground.penetration(LengthVec::m(0.0, 0.0, -0.003)).in_mm() - 3.0).abs() < 1e-9);
        // A tilted ground works the same way along its own normal.
        let ramp = Ground {
            normal: DVec3::new(0.0, 1.0, 1.0).normalize(),
            offset: Length::ZERO,
        };
        assert_eq!(ramp.penetration(LengthVec::m(0.0, 1.0, 1.0)), Length::ZERO);
        assert!(ramp.penetration(LengthVec::m(0.0, -1.0, -1.0)).to_si() > 0.0);
    }
}

#[cfg(test)]
mod ledger_scale_probe {
    use super::*;

    /// **The audit has to have something to judge against.** A ledger whose entry scale is
    /// zero is skipped entirely by `audit`, so a conserved quantity that is correctly zero
    /// would be policed by nothing at all.
    ///
    /// This is the regression test for exactly that: the momentum ledger once handed the
    /// pre-summed total to a single `with`, which set the scale to `|total|` — zero for any
    /// symmetric system, which is every system worth checking.
    #[test]
    fn the_momentum_ledger_carries_a_usable_scale() {
        let bodies = [
            Body::new(
                Mass::kg(2.0e30),
                LengthVec::m(-1.0e11, 0.0, 0.0),
                VelocityVec::m_per_s(0.0, -1.0e4, 0.0),
            ),
            Body::new(
                Mass::kg(2.0e30),
                LengthVec::m(1.0e11, 0.0, 0.0),
                VelocityVec::m_per_s(0.0, 1.0e4, 0.0),
            ),
        ];
        let system = NBody::new("pair", &bodies);
        let ledger = system.ledger();

        // The total is exactly zero, which is the correct physics and the reason the scale
        // cannot come from it.
        assert_eq!(ledger.get(conserved::MOMENTUM_Y), Some(0.0));

        // The scale is one body's momentum, which is what a relative tolerance needs.
        // Relative, not absolute: one ulp of 2e34 is about 4e18, so an absolute bound of 1.0
        // would be asking for sub-ulp agreement. Written the wrong way round first, which is
        // the "precision finer than the representation" item on this workspace's own list.
        let scale = ledger.scale_of(conserved::MOMENTUM_Y).unwrap();
        assert!(
            (scale / 2.0e34 - 1.0).abs() < 1e-12,
            "scale should be one body's 2e34 kg m/s, got {scale}"
        );
        assert!(
            scale > 1e-300,
            "below this `audit` skips the entry and polices nothing"
        );
    }
}
