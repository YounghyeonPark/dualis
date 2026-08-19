//! Bodies that hit each other, and bodies that roll.
//!
//! [`ContactSystem`](crate::ContactSystem) is particles against a plane: no
//! orientation, so nothing can spin up and nothing can roll.
//! [`RigidBody`](crate::RigidBody) rotates but never touches anything. This is where
//! the two meet, and the thing that joins them is friction — a force applied *away
//! from* a body's centre, which is the only way linear and angular motion talk to each
//! other.
//!
//! # Why the closed forms here are unusually good
//!
//! Rolling has an exact answer that most people have met and few would guess from the
//! equations. A solid sphere released on an incline accelerates at
//!
//! ```text
//! a = g sinθ / (1 + I/mr²) = 5/7 g sinθ
//! ```
//!
//! and a hoop, with all its mass at the rim, manages only `½ g sinθ`. Neither depends
//! on the mass or the radius. If it slips instead, the answer changes to
//! `g(sinθ − μcosθ)` and the changeover is at exactly `μ = (2/7) tanθ` for a sphere.
//!
//! Three exact statements and a sharp boundary between two of them, none of it obvious
//! from the code that produces it — which is what makes this module testable in the way
//! the rest of the workspace is.
//!
//! # Impulses conserve what forces conserve, and for the same reason
//!
//! A collision here is resolved by an impulse applied equal and opposite to the two
//! bodies at their contact point. That is Newton's third law again, structurally, so
//! linear momentum is conserved to the last bit — and because both impulses act at the
//! *same point*, angular momentum about any origin is too. Neither is enforced
//! afterwards; both fall out of applying one vector twice with opposite signs.

use glam::DVec3;
use pantometry_core::conserved::quantity;
use pantometry_core::{Domain, Exchange, Kind, Ledger, Violation};
use pantometry_units::{Energy, Length, LengthVec, Mass, Time, VelocityVec};

use crate::rigid::Inertia;

/// A sphere with a place, a velocity and a spin.
///
/// A sphere because its inertia is the same about every axis, which removes the
/// orientation from the dynamics entirely: a rolling ball's contact behaviour does not
/// depend on which way up it is. Anything else would need the full orientation carried
/// through the collision, and that is a different module.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sphere {
    /// How much of it there is.
    pub mass: Mass,
    /// Its radius, which is where contact happens.
    pub radius: Length,
    /// Centre position.
    pub position: LengthVec,
    /// Velocity of the centre.
    pub velocity: VelocityVec,
    /// Angular velocity in the world frame, rad/s. A sphere is a spherical top, so
    /// there is no body frame to distinguish.
    pub spin: DVec3,
    /// Moment of inertia about any axis through the centre.
    pub inertia: f64,
}

impl Sphere {
    /// A uniform solid ball: `⅖mr²`.
    pub fn solid(mass: Mass, radius: Length, position: LengthVec) -> Sphere {
        Sphere {
            mass,
            radius,
            position,
            velocity: VelocityVec::ZERO,
            spin: DVec3::ZERO,
            inertia: Inertia::solid_sphere(mass, radius).principal.x,
        }
    }

    /// A thin shell: `⅔mr²`. Rolls more slowly than a solid ball of the same mass,
    /// because more of it has to be spun up.
    pub fn shell(mass: Mass, radius: Length, position: LengthVec) -> Sphere {
        Sphere {
            mass,
            radius,
            position,
            velocity: VelocityVec::ZERO,
            spin: DVec3::ZERO,
            inertia: Inertia::hollow_sphere(mass, radius).principal.x,
        }
    }

    /// A hoop rolling in a plane: `mr²`, all the mass at the rim.
    ///
    /// Not a sphere, but it behaves as one for rolling in two dimensions and it is the
    /// other end of the range the rolling closed form covers.
    pub fn hoop(mass: Mass, radius: Length, position: LengthVec) -> Sphere {
        Sphere {
            mass,
            radius,
            position,
            velocity: VelocityVec::ZERO,
            spin: DVec3::ZERO,
            inertia: mass.to_si() * radius.to_si().powi(2),
        }
    }

    /// Set the centre-of-mass velocity.
    pub fn with_velocity(mut self, velocity: VelocityVec) -> Sphere {
        self.velocity = velocity;
        self
    }

    /// Set the angular velocity, rad/s in the world frame.
    pub fn with_spin(mut self, spin: DVec3) -> Sphere {
        self.spin = spin;
        self
    }

    /// The dimensionless `I/mr²`, which is the only thing about a body's shape that
    /// rolling cares about: ⅖ for a ball, ⅔ for a shell, 1 for a hoop.
    pub fn inertia_ratio(&self) -> f64 {
        let mr2 = self.mass.to_si() * self.radius.to_si().powi(2);
        if mr2 <= 0.0 {
            return 0.0;
        }
        self.inertia / mr2
    }

    /// Velocity of the material point at `offset` from the centre: `v + ω × r`.
    ///
    /// The quantity a contact actually sees. A rolling wheel's contact point is
    /// stationary even though the wheel is moving, and that is this expression coming
    /// out zero.
    pub fn point_velocity(&self, offset: DVec3) -> DVec3 {
        self.velocity.to_si() + self.spin.cross(offset)
    }

    /// Translational plus rotational: `½mv² + ½Iω²`.
    ///
    /// Both halves, because a rolling sphere keeps two sevenths of its energy in the spin and
    /// a model that forgets it gets the wrong answer for a ball down a slope.
    pub fn kinetic_energy(&self) -> Energy {
        let v = self.velocity.to_si().length();
        Energy::from_si(
            0.5 * self.mass.to_si() * v * v + 0.5 * self.inertia * self.spin.length_squared(),
        )
    }

    /// `mv` for the centre of mass. The spin carries angular momentum, not linear.
    pub fn momentum(&self) -> pantometry_units::MomentumVec {
        self.mass * self.velocity
    }

    /// Angular momentum about a point: the spin plus the orbital part.
    pub fn angular_momentum_about(&self, origin: LengthVec) -> DVec3 {
        let r = (self.position - origin).to_si();
        self.inertia * self.spin + r.cross(self.momentum().to_si())
    }

    /// Apply an impulse at a point offset from the centre.
    ///
    /// Changes both the velocity and the spin, and the spin only because the impulse
    /// was applied off-centre. Hitting a ball through its middle cannot make it turn,
    /// which is the whole content of the cross product below.
    pub fn apply_impulse(&mut self, impulse: DVec3, offset: DVec3) {
        let m = self.mass.to_si();
        if m > 0.0 {
            self.velocity += VelocityVec::from_si(impulse / m);
        }
        if self.inertia > 0.0 {
            self.spin += offset.cross(impulse) / self.inertia;
        }
    }
}

/// How bouncy and how grippy a contact is.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Surface {
    /// Coefficient of restitution, 0 (sticks) to 1 (perfectly elastic).
    pub restitution: f64,
    /// Coulomb friction coefficient. The tangential impulse is capped at `μ` times the
    /// normal one, which is what makes friction a limit rather than a force.
    pub friction: f64,
}

impl Surface {
    /// Perfectly elastic and perfectly slippery: momentum and energy both conserved.
    pub fn frictionless_elastic() -> Surface {
        Surface {
            restitution: 1.0,
            friction: 0.0,
        }
    }

    /// A surface with a coefficient of restitution and one of friction.
    ///
    /// Restitution is clamped to `0..=1`: above one a bounce would return more energy than it
    /// arrived with, which is the one thing this workspace refuses everywhere.
    pub fn new(restitution: f64, friction: f64) -> Surface {
        Surface {
            restitution: restitution.clamp(0.0, 1.0),
            friction: friction.max(0.0),
        }
    }
}

/// Resolve a collision between two spheres by an impulse, returning the energy lost.
///
/// The normal impulse comes from the restitution; the tangential one opposes the
/// sliding at the contact and is capped by Coulomb's limit, so a grippy contact stops
/// the sliding outright and a slippery one only slows it.
///
/// Both bodies receive the same impulse with opposite signs at the same point, so
/// linear and angular momentum are conserved by construction rather than by correction.
pub fn resolve(a: &mut Sphere, b: &mut Sphere, surface: Surface) -> Energy {
    let separation = (b.position - a.position).to_si();
    let distance = separation.length();
    if distance <= 0.0 {
        return Energy::from_si(0.0);
    }
    let normal = separation / distance;

    // Contact point offsets from each centre.
    let ra = normal * a.radius.to_si();
    let rb = -normal * b.radius.to_si();

    let relative = b.point_velocity(rb) - a.point_velocity(ra);
    let approach = relative.dot(normal);
    if approach >= 0.0 {
        // Already separating; resolving again would add energy.
        return Energy::from_si(0.0);
    }

    let before = (a.kinetic_energy() + b.kinetic_energy()).to_si();
    let (ma, mb) = (a.mass.to_si(), b.mass.to_si());
    let inv_mass = 1.0 / ma + 1.0 / mb;

    // Normal impulse. The angular terms drop out for a central impulse, because
    // `r x n` is zero when `r` is along `n`.
    let jn = -(1.0 + surface.restitution) * approach / inv_mass;
    a.apply_impulse(-normal * jn, ra);
    b.apply_impulse(normal * jn, rb);

    // Tangential impulse, opposing whatever sliding is left at the contact.
    if surface.friction > 0.0 {
        let relative = b.point_velocity(rb) - a.point_velocity(ra);
        let tangential = relative - normal * relative.dot(normal);
        let speed = tangential.length();
        if speed > 0.0 {
            let direction = tangential / speed;
            // The impulse that would stop the sliding entirely, given that an impulse at
            // the rim also changes the spin.
            let effective = inv_mass
                + a.radius.to_si().powi(2) / a.inertia.max(f64::MIN_POSITIVE)
                + b.radius.to_si().powi(2) / b.inertia.max(f64::MIN_POSITIVE);
            let stopping = speed / effective;
            // Coulomb: friction cannot deliver more than mu times the normal impulse,
            // which is what separates a grip from a slide.
            let jt = stopping.min(surface.friction * jn);
            a.apply_impulse(direction * jt, ra);
            b.apply_impulse(-direction * jt, rb);
        }
    }

    let after = (a.kinetic_energy() + b.kinetic_energy()).to_si();
    Energy::from_si((before - after).max(0.0))
}

/// A ball on a slope, with friction deciding whether it rolls or slides.
///
/// One body and one plane, which is enough for the closed forms that make this module
/// worth having and no more than that.
pub struct Rolling {
    name: String,
    ball: Sphere,
    /// Incline angle from horizontal, radians.
    slope: f64,
    surface: Surface,
    gravity: f64,
    /// Distance travelled along the slope.
    travelled: f64,
    dissipated: f64,
    saved: Option<(Sphere, f64, f64)>,
}

impl Rolling {
    /// A ball released from rest on a slope.
    ///
    /// The coordinate is one-dimensional: `x` runs down the slope and the ball's spin is
    /// about the horizontal axis across it. Gravity's component along the slope is
    /// `g sinθ` and the normal load is `mg cosθ`, which is all the geometry there is.
    pub fn new(
        name: impl Into<String>,
        ball: Sphere,
        slope_radians: f64,
        surface: Surface,
    ) -> Rolling {
        Rolling {
            name: name.into(),
            ball,
            slope: slope_radians,
            surface,
            gravity: pantometry_units::G0.to_si(),
            travelled: 0.0,
            dissipated: 0.0,
            saved: None,
        }
    }

    /// The sphere as it is now, part way down the slope.
    pub fn ball(&self) -> Sphere {
        self.ball
    }

    /// Speed down the slope.
    pub fn speed(&self) -> f64 {
        self.ball.velocity.to_si().x
    }

    /// Spin rate about the axis across the slope.
    pub fn spin(&self) -> f64 {
        self.ball.spin.z
    }

    /// Velocity of the material point touching the ground: zero when rolling.
    ///
    /// The definition of rolling without slipping, and the thing the friction is trying
    /// to bring to zero.
    pub fn slip_velocity(&self) -> f64 {
        self.speed() - self.ball.radius.to_si() * self.spin()
    }

    /// How far along the slope it has travelled.
    pub fn distance(&self) -> Length {
        Length::from_si(self.travelled)
    }

    /// Energy lost to friction while slipping. Zero once it rolls without slipping, which
    /// is why a rolling ball does not heat up.
    pub fn dissipated_energy(&self) -> Energy {
        Energy::from_si(self.dissipated)
    }

    /// Acceleration of a body that rolls without slipping: `g sinθ / (1 + I/mr²)`.
    ///
    /// ⅚ of `g sinθ`... no: ⁵⁄₇ for a solid ball, ³⁄₅ for a shell, ½ for a hoop. It does
    /// not depend on the mass or the radius, only on how the mass is distributed, which
    /// is why a marble and a cannonball reach the bottom together and a hoop does not.
    pub fn rolling_acceleration(&self) -> f64 {
        self.gravity * self.slope.sin() / (1.0 + self.ball.inertia_ratio())
    }

    /// Acceleration of a body that slides with friction: `g(sinθ − μcosθ)`.
    pub fn sliding_acceleration(&self) -> f64 {
        self.gravity * (self.slope.sin() - self.surface.friction * self.slope.cos())
    }

    /// The friction needed to roll rather than slide, `tanθ · (I/mr²)/(1 + I/mr²)`.
    ///
    /// `(2/7) tanθ` for a solid ball. Below this the ball cannot spin up fast enough to
    /// keep its contact point still and it slides; above it, adding more friction
    /// changes nothing, because friction is a limit and not a force.
    pub fn friction_to_roll(&self) -> f64 {
        let k = self.ball.inertia_ratio();
        self.slope.tan() * k / (1.0 + k)
    }

    /// Whether the ball is currently rolling rather than sliding.
    pub fn is_rolling(&self) -> bool {
        self.slip_velocity().abs() < 1e-9 * self.speed().abs().max(1.0)
    }
}

impl Domain for Rolling {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// No stiffness and no wave, so nothing here sets a stability limit — the equations
    /// are smooth and first order in the state. The one scale that exists is how long it
    /// takes friction to bring the slip to zero, and that is resolved by any step short
    /// enough to be interesting.
    fn max_stable_dt(&self, _now: Time) -> Time {
        Time::from_si(f64::INFINITY)
    }

    /// Semi-implicit Euler on the slope coordinate, with friction resolved as a limit
    /// each step.
    ///
    /// While the contact is sliding, friction is `μmg cosθ` and points against the slip.
    /// Once the slip reaches zero it becomes whatever is needed to keep it there, capped
    /// by the same `μmg cosθ` — which is what makes Coulomb friction a constraint that
    /// switches on rather than a force that is always present.
    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let h = dt.to_si();
        if h <= 0.0 {
            return Ok(());
        }
        let m = self.ball.mass.to_si();
        let r = self.ball.radius.to_si();
        let i = self.ball.inertia;
        let normal_load = m * self.gravity * self.slope.cos();
        let along = m * self.gravity * self.slope.sin();
        let max_friction = self.surface.friction * normal_load;

        let slip = self.slip_velocity();

        // One rule rather than a sliding branch and a rolling branch: find the friction
        // that would bring the slip to exactly zero by the end of this step, and use it
        // if Coulomb allows.
        //
        // The slip closes at `d(slip)/dt = F/m + f(1/m + r²/I)`, so
        //
        //   f_stop = (−slip/h − F/m) / (1/m + r²/I)
        //
        // With `slip = 0` that reduces to `−F k/(1+k)` with `k = I/mr²`, which is the
        // classical rolling friction — so the rolling case is not special-cased, it is
        // what the rule gives when there is nothing left to stop.
        //
        // Branching instead was the first attempt and it chatters. Once the slip crosses
        // zero between two steps the sign of the friction flips, the next step overshoots
        // the other way, and the contact oscillates at a few times `h` forever instead of
        // settling — so a ball that should have started rolling never quite does.
        let mobility = 1.0 / m + if i > 0.0 { r * r / i } else { 0.0 };
        let f_stop = if mobility > 0.0 {
            (-slip / h - along / m) / mobility
        } else {
            0.0
        };
        let friction = if f_stop.abs() <= max_friction {
            f_stop
        } else {
            f_stop.signum() * max_friction
        };

        let accel = (along + friction) / m;
        let angular_accel = if i > 0.0 { -friction * r / i } else { 0.0 };

        let v0 = self.speed();
        let new_v = v0 + accel * h;
        let new_w = self.spin() + angular_accel * h;
        self.ball.velocity = VelocityVec::from_si(DVec3::new(new_v, 0.0, 0.0));
        self.ball.spin = DVec3::new(0.0, 0.0, new_w);
        self.travelled += new_v * h;
        self.ball.position += LengthVec::from_si(DVec3::new(new_v * h, 0.0, 0.0));

        // Friction dissipates only where the surfaces slide against each other, and the
        // work is done over the slip that actually occurred — the average across the
        // step, not the value it started at, since the step is often the one that brings
        // the slip to zero. A rolling contact slides nothing and so heats nothing, which
        // is why rolling is efficient and why a rolling wheel does not get warm.
        let slip_after = new_v - r * new_w;
        let mean_slip = 0.5 * (slip + slip_after);
        if mean_slip.abs() > 0.0 {
            let heat = (friction * mean_slip).abs() * h;
            self.dissipated += heat;
            bus.publish(quantity::ENERGY, heat);
        }
        Ok(())
    }

    /// Kinetic energy plus what has been dissipated, less the potential given up. The
    /// sum holds however the ball is moving.
    fn ledger(&self) -> Ledger {
        let dropped = self.travelled * self.slope.sin() * self.ball.mass.to_si() * self.gravity;
        Ledger::new().with(
            quantity::ENERGY,
            self.ball.kinetic_energy().to_si() + self.dissipated - dropped,
        )
    }

    fn checkpoint(&mut self) {
        self.saved = Some((self.ball, self.travelled, self.dissipated));
    }

    fn restore(&mut self) {
        if let Some((ball, travelled, dissipated)) = self.saved {
            self.ball = ball;
            self.travelled = travelled;
            self.dissipated = dissipated;
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

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pantometry_core::{Schedule, Simulation};

    fn ball_at(x: f64) -> Sphere {
        Sphere::solid(Mass::kg(0.2), Length::m(0.05), LengthVec::m(x, 0.0, 0.0))
    }

    /// Two equal masses meeting head-on, perfectly elastic: they swap velocities. The
    /// oldest closed form in mechanics and the first thing a collision routine should be
    /// held to.
    #[test]
    fn equal_masses_swap_velocities() {
        let mut a = ball_at(-1.0).with_velocity(VelocityVec::m_per_s(3.0, 0.0, 0.0));
        let mut b = ball_at(1.0).with_velocity(VelocityVec::m_per_s(-1.0, 0.0, 0.0));
        let lost = resolve(&mut a, &mut b, Surface::frictionless_elastic());

        assert!(
            (a.velocity.to_si().x + 1.0).abs() < 1e-12,
            "{:?}",
            a.velocity
        );
        assert!(
            (b.velocity.to_si().x - 3.0).abs() < 1e-12,
            "{:?}",
            b.velocity
        );
        assert!(
            lost.to_si().abs() < 1e-12,
            "an elastic collision loses nothing"
        );
        // Neither picked up any spin, because the impulse went through both centres.
        assert_eq!(a.spin, DVec3::ZERO);
        assert_eq!(b.spin, DVec3::ZERO);
    }

    /// Momentum and angular momentum survive a collision exactly, and for the same
    /// structural reason forces do: one impulse, applied twice with opposite signs at
    /// one point.
    #[test]
    fn an_impulse_conserves_both_momenta_exactly() {
        for (restitution, friction) in [(1.0, 0.0), (0.5, 0.3), (0.0, 1.0)] {
            let mut a = Sphere::solid(
                Mass::kg(0.3),
                Length::m(0.04),
                LengthVec::m(-0.05, 0.0, 0.0),
            )
            .with_velocity(VelocityVec::m_per_s(2.0, 0.7, 0.0))
            .with_spin(DVec3::new(0.0, 0.0, 12.0));
            let mut b = Sphere::solid(Mass::kg(0.7), Length::m(0.06), LengthVec::m(0.05, 0.0, 0.0))
                .with_velocity(VelocityVec::m_per_s(-0.5, -0.2, 0.0))
                .with_spin(DVec3::new(0.0, 0.0, -3.0));

            let origin = LengthVec::ZERO;
            let p0 = (a.momentum() + b.momentum()).to_si();
            let l0 = a.angular_momentum_about(origin) + b.angular_momentum_about(origin);
            let e0 = (a.kinetic_energy() + b.kinetic_energy()).to_si();

            let lost = resolve(&mut a, &mut b, Surface::new(restitution, friction));

            let p1 = (a.momentum() + b.momentum()).to_si();
            let l1 = a.angular_momentum_about(origin) + b.angular_momentum_about(origin);
            assert!(
                (p1 - p0).length() / p0.length() < 1e-14,
                "e = {restitution}, mu = {friction}: momentum moved"
            );
            assert!(
                (l1 - l0).length() / l0.length() < 1e-14,
                "e = {restitution}, mu = {friction}: angular momentum moved"
            );
            // And the energy the collision reports losing is the energy it lost.
            let e1 = (a.kinetic_energy() + b.kinetic_energy()).to_si();
            assert!(
                (lost.to_si() - (e0 - e1)).abs() / e0 < 1e-12,
                "the reported loss should be the actual loss"
            );
            assert!(lost.to_si() >= -1e-15, "a collision cannot add energy");
        }
    }

    /// A glancing blow with friction imparts spin, and a central one cannot. That
    /// distinction is the entire reason a body needs an orientation.
    #[test]
    fn friction_turns_a_glancing_blow_into_spin() {
        let make = || {
            (
                ball_at(-0.05).with_velocity(VelocityVec::m_per_s(2.0, 1.0, 0.0)),
                ball_at(0.05),
            )
        };
        // Slippery: the tangential motion passes straight through and nothing turns.
        let (mut a, mut b) = make();
        resolve(&mut a, &mut b, Surface::frictionless_elastic());
        assert_eq!(a.spin, DVec3::ZERO);
        assert_eq!(b.spin, DVec3::ZERO);
        assert!(
            (a.velocity.to_si().y - 1.0).abs() < 1e-12,
            "a frictionless contact cannot change the tangential velocity"
        );

        // Grippy: the same blow spins both, in opposite senses.
        let (mut a, mut b) = make();
        resolve(&mut a, &mut b, Surface::new(1.0, 0.6));
        assert!(a.spin.length() > 1.0, "a should be turning, got {}", a.spin);
        assert!(b.spin.length() > 1.0, "b should be turning, got {}", b.spin);
        assert!(
            a.spin.z * b.spin.z > 0.0,
            "rubbing past each other turns them the same way, not opposite ways"
        );
        assert!(a.velocity.to_si().y < 1.0, "friction slowed the sliding");
    }

    /// A dead collision sticks: the two leave with a common velocity, and the energy
    /// that went missing is exactly what a perfectly inelastic collision must lose.
    #[test]
    fn a_dead_collision_sticks_and_loses_the_right_amount() {
        let mut a = Sphere::solid(
            Mass::kg(1.0),
            Length::m(0.05),
            LengthVec::m(-0.05, 0.0, 0.0),
        )
        .with_velocity(VelocityVec::m_per_s(4.0, 0.0, 0.0));
        let mut b = Sphere::solid(Mass::kg(3.0), Length::m(0.05), LengthVec::m(0.05, 0.0, 0.0));
        let lost = resolve(&mut a, &mut b, Surface::new(0.0, 0.0));

        // Common velocity is the momentum over the total mass: 4 kg m/s over 4 kg.
        assert!((a.velocity.to_si().x - 1.0).abs() < 1e-12);
        assert!((b.velocity.to_si().x - 1.0).abs() < 1e-12);
        // The closed form for the loss: (1/2) m1 m2/(m1+m2) * (relative speed)^2.
        let expected = 0.5 * (1.0 * 3.0 / 4.0) * 16.0;
        assert!(
            (lost.to_si() - expected).abs() < 1e-12,
            "lost {} J, closed form says {expected} J",
            lost.to_si()
        );
    }

    /// Bodies already moving apart are left alone. Resolving a separating pair would add
    /// energy from nowhere, and it is the commonest way a collision loop explodes.
    #[test]
    fn a_separating_pair_is_not_resolved() {
        let mut a = ball_at(-0.05).with_velocity(VelocityVec::m_per_s(-1.0, 0.0, 0.0));
        let mut b = ball_at(0.05).with_velocity(VelocityVec::m_per_s(1.0, 0.0, 0.0));
        let (va, vb) = (a.velocity, b.velocity);
        let lost = resolve(&mut a, &mut b, Surface::frictionless_elastic());
        assert_eq!(a.velocity, va);
        assert_eq!(b.velocity, vb);
        assert_eq!(lost.to_si(), 0.0);
    }

    /// **The closed form this module is worth having for.** A rolling body accelerates
    /// at `g sinθ/(1 + I/mr²)`, and the fraction depends only on how its mass is
    /// arranged: ⁵⁄₇ for a solid ball, ³⁄₅ for a shell, ½ for a hoop.
    ///
    /// So a marble and a cannonball reach the bottom together and a hoop arrives last,
    /// however heavy any of them is — which is not obvious and is exactly what the
    /// integration should reproduce.
    #[test]
    fn a_rolling_body_accelerates_by_its_shape_and_nothing_else() {
        let slope = 20f64.to_radians();
        let g = pantometry_units::G0.to_si();
        let grippy = Surface::new(0.0, 1.0);

        let cases = [
            (
                "solid ball",
                Sphere::solid(Mass::kg(0.2), Length::m(0.05), LengthVec::ZERO),
                5.0 / 7.0,
            ),
            (
                "shell",
                Sphere::shell(Mass::kg(0.2), Length::m(0.05), LengthVec::ZERO),
                3.0 / 5.0,
            ),
            (
                "hoop",
                Sphere::hoop(Mass::kg(0.2), Length::m(0.05), LengthVec::ZERO),
                0.5,
            ),
        ];
        for (name, ball, fraction) in cases {
            let mut system = Rolling::new("slope", ball, slope, grippy);
            let predicted = g * slope.sin() * fraction;
            assert!(
                (system.rolling_acceleration() / predicted - 1.0).abs() < 1e-12,
                "{name}: the closed form should be {fraction} of g sin(theta)"
            );

            // Integrate and measure the acceleration that actually came out.
            let mut bus = Exchange::new();
            let dt = Time::ms(0.1);
            for _ in 0..10_000 {
                system.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            let elapsed = 1.0;
            let measured = system.speed() / elapsed;
            assert!(
                (measured / predicted - 1.0).abs() < 1e-6,
                "{name}: integrated {measured:.5} m/s^2 against a predicted \
                 {predicted:.5}"
            );
            // Rolling: the contact point is not moving.
            assert!(system.is_rolling(), "{name} should be rolling");
            // And rolling dissipates nothing, because the contact does no sliding.
            assert!(
                system.dissipated_energy().to_si() < 1e-12,
                "{name}: rolling should not heat anything"
            );
        }

        // Neither the mass nor the radius enters.
        let heavy = Rolling::new(
            "heavy",
            Sphere::solid(Mass::kg(50.0), Length::m(0.4), LengthVec::ZERO),
            slope,
            grippy,
        );
        let light = Rolling::new(
            "light",
            Sphere::solid(Mass::kg(0.01), Length::m(0.005), LengthVec::ZERO),
            slope,
            grippy,
        );
        assert!(
            (heavy.rolling_acceleration() - light.rolling_acceleration()).abs() < 1e-12,
            "a marble and a cannonball reach the bottom together"
        );
    }

    /// Too little friction and it slides instead, at `g(sinθ − μcosθ)` — a different
    /// closed form, and a slower one, with the changeover at exactly `μ = (2/7)tanθ`.
    #[test]
    fn too_little_friction_slides_at_the_other_closed_form() {
        let slope = 30f64.to_radians();
        let ball = || Sphere::solid(Mass::kg(0.2), Length::m(0.05), LengthVec::ZERO);

        // (2/7) tan(30) = 0.1650.
        let threshold = Rolling::new("s", ball(), slope, Surface::new(0.0, 0.0)).friction_to_roll();
        assert!(
            (threshold - 2.0 / 7.0 * slope.tan()).abs() < 1e-12,
            "the rolling threshold is (2/7) tan(theta)"
        );
        assert!((threshold - 0.1650).abs() < 1e-3, "got {threshold}");

        // Below it, the ball slides and takes the sliding acceleration.
        let mut sliding = Rolling::new("slide", ball(), slope, Surface::new(0.0, 0.08));
        let predicted = sliding.sliding_acceleration();
        let mut bus = Exchange::new();
        for _ in 0..5000 {
            sliding.step(Time::ZERO, Time::ms(0.1), &mut bus).unwrap();
        }
        let measured = sliding.speed() / 0.5;
        assert!(
            (measured / predicted - 1.0).abs() < 1e-6,
            "sliding: {measured:.5} against a predicted {predicted:.5}"
        );
        assert!(!sliding.is_rolling(), "it should still be slipping");
        assert!(
            sliding.dissipated_energy().to_si() > 0.0,
            "sliding makes heat"
        );

        // And a weakly-frictional slide gets *down the slope* faster than rolling would,
        // which is the opposite of what "friction slows things down" suggests. Rolling
        // has to spend 2/7 of `g sinθ` on spinning the ball up, and at mu = 0.08 the
        // friction only takes back 0.08 g cosθ — less than that. A ball on ice beats a
        // ball on tarmac to the bottom, and arrives without turning.
        assert!(
            predicted > sliding.rolling_acceleration(),
            "a weak slide should outrun a roll: {predicted:.4} against {:.4}",
            sliding.rolling_acceleration()
        );

        // Above the threshold, more friction changes nothing: friction is a limit.
        let accel_of = |mu: f64| {
            let mut r = Rolling::new("r", ball(), slope, Surface::new(0.0, mu));
            let mut bus = Exchange::new();
            for _ in 0..2000 {
                r.step(Time::ZERO, Time::ms(0.1), &mut bus).unwrap();
            }
            r.speed() / 0.2
        };
        let (just_over, far_over) = (accel_of(0.2), accel_of(2.0));
        assert!(
            (just_over / far_over - 1.0).abs() < 1e-9,
            "past the threshold, adding grip does nothing: {just_over} against {far_over}"
        );
    }

    /// A ball started sliding on a grippy slope spins up until it rolls, and then stops
    /// dissipating. The transition is the interesting part: friction is doing work right
    /// up until the contact point stops moving, and none afterwards.
    #[test]
    fn sliding_becomes_rolling_and_the_heating_stops() {
        let slope = 15f64.to_radians();
        let mut system = Rolling::new(
            "slope",
            Sphere::solid(Mass::kg(0.2), Length::m(0.05), LengthVec::ZERO)
                // Dropped onto the slope already moving, with no spin: it must slide
                // before it can roll.
                .with_velocity(VelocityVec::m_per_s(3.0, 0.0, 0.0)),
            slope,
            Surface::new(0.0, 0.5),
        );
        assert!(!system.is_rolling(), "it starts out sliding");

        let mut bus = Exchange::new();
        let dt = Time::ms(0.05);
        let mut heat_while_sliding = 0.0;
        let mut steps_to_roll = 0;
        for i in 0..40_000 {
            system.step(Time::ZERO, dt, &mut bus).unwrap();
            heat_while_sliding += bus.take(quantity::ENERGY);
            if system.is_rolling() && steps_to_roll == 0 {
                steps_to_roll = i + 1;
            }
        }
        assert!(steps_to_roll > 0, "it should have settled into rolling");
        assert!(heat_while_sliding > 0.0, "sliding makes heat");

        // Once rolling, nothing more is published: a rolling contact does no work.
        let after_rolling = system.dissipated_energy().to_si();
        for _ in 0..1000 {
            system.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        assert!(
            (system.dissipated_energy().to_si() - after_rolling).abs() < 1e-12,
            "a rolling wheel does not get hot"
        );
        assert!(bus.peek(quantity::ENERGY).abs() < 1e-30);
    }

    /// The domain balances under the scheduler: kinetic plus heat is the potential given
    /// up, whether it rolls, slides, or does one and then the other.
    #[test]
    fn the_books_balance_however_it_moves() {
        for mu in [0.0f64, 0.08, 0.5] {
            struct Sink;
            impl Domain for Sink {
                fn name(&self) -> &str {
                    "sink"
                }
                fn step(
                    &mut self,
                    _t: Time,
                    _dt: Time,
                    bus: &mut Exchange,
                ) -> Result<(), Violation> {
                    bus.take(quantity::ENERGY);
                    Ok(())
                }
                fn checkpoint(&mut self) {}
                fn restore(&mut self) {}
                fn supports_restore(&self) -> bool {
                    true
                }
            }

            let mut sim = Simulation::new(Schedule::Staggered)
                // First-order integration of a smooth problem: the energy books close to
                // the step size, and 0.1 ms over a second is a part in ten thousand.
                .conservation_tolerance(1e-3)
                .with(Rolling::new(
                    "slope",
                    Sphere::solid(Mass::kg(0.2), Length::m(0.05), LengthVec::ZERO)
                        .with_velocity(VelocityVec::m_per_s(2.0, 0.0, 0.0)),
                    25f64.to_radians(),
                    Surface::new(0.0, mu),
                ))
                .with(Sink);
            for _ in 0..100 {
                sim.advance(Time::ms(1.0))
                    .unwrap_or_else(|e| panic!("mu = {mu}: {e}"));
            }
        }
    }

    /// Degenerate bodies do not divide by zero.
    #[test]
    fn degenerate_bodies_are_handled() {
        // A point mass has no inertia, so an impulse cannot spin it.
        let mut point = Sphere {
            mass: Mass::kg(1.0),
            radius: Length::ZERO,
            position: LengthVec::ZERO,
            velocity: VelocityVec::ZERO,
            spin: DVec3::ZERO,
            inertia: 0.0,
        };
        point.apply_impulse(DVec3::X, DVec3::Y);
        assert_eq!(point.spin, DVec3::ZERO);
        assert!((point.velocity.to_si().x - 1.0).abs() < 1e-15);
        assert_eq!(point.inertia_ratio(), 0.0);

        // Two spheres at the same place cannot be resolved, and say so by doing nothing.
        let mut a = ball_at(0.0);
        let mut b = ball_at(0.0);
        assert_eq!(
            resolve(&mut a, &mut b, Surface::frictionless_elastic()).to_si(),
            0.0
        );

        // A flat slope goes nowhere.
        let mut flat = Rolling::new("flat", ball_at(0.0), 0.0, Surface::new(0.0, 0.5));
        let mut bus = Exchange::new();
        for _ in 0..100 {
            flat.step(Time::ZERO, Time::ms(1.0), &mut bus).unwrap();
        }
        assert!(flat.speed().abs() < 1e-15);
        assert_eq!(flat.rolling_acceleration(), 0.0);
    }
}
