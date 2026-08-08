//! Free rotation of a rigid body, and the third quantity the kernel audits.
//!
//! Everything else in this crate is a point mass. A point has no orientation, so it
//! cannot tumble, cannot store energy in a spin, and cannot show the one thing about
//! rotation that surprises people: **a body spun about its middle axis will not stay
//! that way.**
//!
//! # Angular momentum, and why it is worth the crate
//!
//! Energy came first, momentum second, and both were audited by the same machinery.
//! Angular momentum is the third, and it is the one that could most easily have been
//! an accident of how the first two were written: it is conserved for a different
//! reason (rotational rather than translational symmetry), it lives in the world frame
//! while the equations of motion live in the body frame, and getting the frames
//! confused produces something that looks conserved and is not.
//!
//! # Euler's equations, in the body frame
//!
//! In a frame aligned with the principal axes the inertia tensor is diagonal and the
//! equations are three lines:
//!
//! ```text
//! I₁ ω̇₁ = (I₂ - I₃) ω₂ ω₃
//! I₂ ω̇₂ = (I₃ - I₁) ω₃ ω₁
//! I₃ ω̇₃ = (I₁ - I₂) ω₁ ω₂
//! ```
//!
//! The right-hand sides vanish when all three moments are equal, which is why a sphere
//! spins forever about whatever axis it was given, and why nothing else does.
//!
//! # What integrating this costs
//!
//! The state is an angular velocity and an orientation, and `ω̇` depends on `ω`, so
//! this is not a [`Newtonian`](dualis_core::Newtonian) system and there is no
//! symplectic method for it. It goes through
//! [`Integrator::Rk4`](dualis_core::Integrator), which is fourth-order accurate and
//! dissipates slowly — so `|L|` drifts, and the audit tolerance is the integrator's
//! rather than the physics'.
//!
//! The quaternion also has to be renormalised every step. Integrating `q̇ = ½ω⊗q`
//! numerically walks off the unit sphere, and an orientation of norm 1.0001 is not an
//! orientation. Renormalising is a projection and is not itself conservative, which is
//! one more reason the tolerance is what it is.

use dualis_core::{Domain, Dynamics, Exchange, Integrator, Kind, Ledger, State, Violation};
use dualis_units::{AngularMomentum, Energy, Length, Mass, MomentOfInertia, Time};
use glam::{DQuat, DVec3};

/// Names for the axes of the angular momentum audit. See the note on
/// [`conserved`](crate::conserved) for why a vector is audited per component.
pub mod conserved {
    /// The `x` component of angular momentum, kg·m²·s⁻¹.
    pub const ANGULAR_MOMENTUM_X: &str = "angular_momentum_x";
    /// The `y` component.
    pub const ANGULAR_MOMENTUM_Y: &str = "angular_momentum_y";
    /// The `z` component.
    pub const ANGULAR_MOMENTUM_Z: &str = "angular_momentum_z";
}

/// Principal moments of inertia, in the body's own frame.
///
/// Stored diagonal, which costs no generality: every real inertia tensor is symmetric
/// and therefore diagonalisable, and the frame that diagonalises it is what "the body
/// frame" means. What it costs is that a body must be *given* in that frame — an
/// arbitrary tensor would need diagonalising first, which is not implemented here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Inertia {
    /// The three principal moments, kg·m², along the body frame axes.
    pub principal: DVec3,
}

impl Inertia {
    /// From three moments in SI, along the body's x, y and z axes.
    pub fn principal(ix: MomentOfInertia, iy: MomentOfInertia, iz: MomentOfInertia) -> Inertia {
        Inertia {
            principal: DVec3::new(ix.to_si(), iy.to_si(), iz.to_si()),
        }
    }

    /// A uniform solid sphere: `⅖mr²` about every axis.
    ///
    /// The degenerate case, and the useful one to test against: with all three moments
    /// equal, Euler's equations have no right-hand side and the body cannot tumble.
    pub fn solid_sphere(mass: Mass, radius: Length) -> Inertia {
        let i = 0.4 * mass.to_si() * radius.to_si().powi(2);
        Inertia {
            principal: DVec3::splat(i),
        }
    }

    /// A thin spherical shell: `⅔mr²`. Larger than a solid sphere of the same mass,
    /// because all of it is at the radius.
    pub fn hollow_sphere(mass: Mass, radius: Length) -> Inertia {
        let i = 2.0 / 3.0 * mass.to_si() * radius.to_si().powi(2);
        Inertia {
            principal: DVec3::splat(i),
        }
    }

    /// A uniform solid cuboid of the given side lengths: `m(b²+c²)/12` about each axis.
    ///
    /// Three unequal moments, which is what it takes to tumble.
    pub fn solid_box(mass: Mass, sides: DVec3) -> Inertia {
        let m = mass.to_si() / 12.0;
        let (a, b, c) = (sides.x, sides.y, sides.z);
        Inertia {
            principal: DVec3::new(
                m * (b * b + c * c),
                m * (a * a + c * c),
                m * (a * a + b * b),
            ),
        }
    }

    /// A uniform solid cylinder about its own axis (z) and across it.
    pub fn solid_cylinder(mass: Mass, radius: Length, length: Length) -> Inertia {
        let (m, r, l) = (mass.to_si(), radius.to_si(), length.to_si());
        let across = m * (3.0 * r * r + l * l) / 12.0;
        Inertia {
            principal: DVec3::new(across, across, 0.5 * m * r * r),
        }
    }

    /// A thin rod along z: `ml²/12` across it and nothing along it.
    pub fn rod(mass: Mass, length: Length) -> Inertia {
        let across = mass.to_si() * length.to_si().powi(2) / 12.0;
        Inertia {
            principal: DVec3::new(across, across, 0.0),
        }
    }

    /// Whether two of the three moments are equal, which makes the body a symmetric
    /// top and its motion a closed form.
    pub fn is_symmetric_top(&self) -> bool {
        let p = self.principal;
        let scale = p.max_element().max(1e-300);
        let close = |a: f64, b: f64| (a - b).abs() / scale < 1e-12;
        close(p.x, p.y) || close(p.y, p.z) || close(p.z, p.x)
    }

    /// Whether all three are equal — a spherical top, which cannot tumble at all.
    pub fn is_spherical_top(&self) -> bool {
        let p = self.principal;
        let scale = p.max_element().max(1e-300);
        (p.x - p.y).abs() / scale < 1e-12 && (p.y - p.z).abs() / scale < 1e-12
    }

    /// The principal axes sorted by moment, as `(smallest, intermediate, largest)`
    /// axis indices.
    ///
    /// The middle one is the unstable axis, and knowing which it is is the whole of
    /// [`Inertia::intermediate_axis`].
    pub fn sorted_axes(&self) -> (usize, usize, usize) {
        let mut order = [0usize, 1, 2];
        order.sort_by(|a, b| self.principal[*a].total_cmp(&self.principal[*b]));
        (order[0], order[1], order[2])
    }

    /// Index of the intermediate principal axis: the one rotation is unstable about.
    pub fn intermediate_axis(&self) -> usize {
        self.sorted_axes().1
    }
}

/// A body with an orientation as well as a place.
pub struct RigidBody {
    name: String,
    inertia: DVec3,
    /// Angular velocity in the **body** frame, rad/s.
    omega_body: DVec3,
    /// Body-to-world rotation.
    orientation: DQuat,
    /// Applied torque in the **world** frame, N·m.
    ///
    /// A `DVec3` rather than a dimensioned type, for the reason `dualis-units`
    /// documents: a torque is newton-metres and so is an energy, and SI cannot tell
    /// them apart because the radian it would need is dimensionless. Angular velocity
    /// is bare for the same reason, so at least the module is consistent about which
    /// quantities the type system cannot help with.
    torque_world: DVec3,
    saved: Option<(DVec3, DQuat)>,
}

/// The integrated state: angular velocity in the body frame, then the quaternion.
#[derive(Clone, Debug, PartialEq)]
pub struct Spin {
    /// Angular velocity in the *body* frame, rad/s. Body frame because that is where the
    /// inertia tensor is diagonal, which is what makes Euler's equations simple.
    pub omega: DVec3,
    /// Body-to-world rotation. Renormalised each step, since integrating a quaternion drifts
    /// off the unit sphere and a non-unit quaternion is a rotation plus a scaling.
    pub orientation: DQuat,
}

impl State for Spin {
    fn axpy(&mut self, a: f64, other: &Self) {
        self.omega += other.omega * a;
        self.orientation = DQuat::from_xyzw(
            self.orientation.x + other.orientation.x * a,
            self.orientation.y + other.orientation.y * a,
            self.orientation.z + other.orientation.z * a,
            self.orientation.w + other.orientation.w * a,
        );
    }

    fn scale(&mut self, a: f64) {
        self.omega *= a;
        self.orientation = DQuat::from_xyzw(
            self.orientation.x * a,
            self.orientation.y * a,
            self.orientation.z * a,
            self.orientation.w * a,
        );
    }

    fn zeros_like(&self) -> Self {
        Spin {
            omega: DVec3::ZERO,
            orientation: DQuat::from_xyzw(0.0, 0.0, 0.0, 0.0),
        }
    }
}

impl RigidBody {
    /// A body with the given inertia, at rest in its authored orientation.
    pub fn new(name: impl Into<String>, inertia: Inertia) -> RigidBody {
        RigidBody {
            name: name.into(),
            inertia: inertia.principal,
            omega_body: DVec3::ZERO,
            orientation: DQuat::IDENTITY,
            torque_world: DVec3::ZERO,
            saved: None,
        }
    }

    /// Apply a constant torque, in the world frame.
    ///
    /// Held until changed, rather than consumed by a step: a motor keeps pushing. A
    /// torque that should act once is an impulse, and belongs in
    /// [`collision`](crate::collision) instead.
    pub fn with_torque(mut self, torque: DVec3) -> RigidBody {
        self.torque_world = torque;
        self
    }

    /// Apply a torque, N·m in the world frame. Constant until set again.
    pub fn set_torque(&mut self, torque: DVec3) {
        self.torque_world = torque;
    }

    /// The torque currently applied, world frame.
    pub fn torque(&self) -> DVec3 {
        self.torque_world
    }

    /// Set the angular velocity, given in the **body** frame.
    pub fn spinning_body_frame(mut self, omega: DVec3) -> RigidBody {
        self.omega_body = omega;
        self
    }

    /// Set the angular velocity, given in the **world** frame.
    pub fn spinning(mut self, omega: DVec3) -> RigidBody {
        self.omega_body = self.orientation.inverse() * omega;
        self
    }

    /// Point the body somewhere, keeping its *world* angular velocity unchanged.
    ///
    /// Which means the body-frame `omega` is rewritten. Turning a body should not change how
    /// it is spinning in space, and storing the velocity in the body frame makes that a
    /// conversion rather than a no-op.
    pub fn with_orientation(mut self, orientation: DQuat) -> RigidBody {
        let world = self.orientation * self.omega_body;
        self.orientation = orientation.normalize();
        self.omega_body = self.orientation.inverse() * world;
        self
    }

    /// The principal moments this body was given.
    pub fn inertia(&self) -> Inertia {
        Inertia {
            principal: self.inertia,
        }
    }

    /// Body-to-world rotation, normalised.
    pub fn orientation(&self) -> DQuat {
        self.orientation
    }

    /// Angular velocity in the body frame.
    pub fn angular_velocity_body(&self) -> DVec3 {
        self.omega_body
    }

    /// Angular velocity in the world frame.
    pub fn angular_velocity(&self) -> DVec3 {
        self.orientation * self.omega_body
    }

    /// Angular momentum in the **world** frame, which is the conserved one.
    ///
    /// In the body frame `L = Iω` componentwise and it is *not* constant — it moves as
    /// the body tumbles. Rotating it out to the world frame is what makes it constant,
    /// and confusing the two is the classic way to produce a quantity that looks
    /// conserved while the body does something else entirely.
    pub fn angular_momentum(&self) -> DVec3 {
        self.orientation * (self.inertia * self.omega_body)
    }

    /// `|L|`, for a quick check. The audit uses the components, because a ledger sums its
    /// entries and the sum of magnitudes is not the magnitude of the sum.
    pub fn angular_momentum_magnitude(&self) -> AngularMomentum {
        AngularMomentum::from_si(self.angular_momentum().length())
    }

    /// `½ ω·Iω`, which is also conserved for a free body.
    ///
    /// Two conserved scalars — this and `|L|` — confine the body-frame angular velocity
    /// to the intersection of an ellipsoid and a sphere. That intersection is a closed
    /// curve for rotation near the largest or smallest axis and a figure of eight
    /// through the intermediate one, which is the whole geometric content of the
    /// instability.
    pub fn rotational_energy(&self) -> Energy {
        let w = self.omega_body;
        Energy::from_si(0.5 * (self.inertia * w).dot(w))
    }

    fn spin_state(&self) -> Spin {
        Spin {
            omega: self.omega_body,
            orientation: self.orientation,
        }
    }
}

impl Dynamics for RigidBody {
    type S = Spin;

    /// Euler's equations with an applied torque, plus the quaternion kinematics.
    ///
    /// `I ω̇ = τ − ω × (Iω)`, with the torque rotated into the body frame using the
    /// orientation *carried in the state* rather than the one on `self`. Runge-Kutta
    /// evaluates the derivative at intermediate orientations, and using the body's
    /// stored orientation for all four stages would apply the torque about the wrong
    /// axes — a mistake that is invisible for a torque along the spin axis and wrong
    /// for every other one.
    fn derivative(&self, s: &Spin, _t: Time) -> Spin {
        let (i1, i2, i3) = (self.inertia.x, self.inertia.y, self.inertia.z);
        let w = s.omega;
        let torque_body = s.orientation.normalize().inverse() * self.torque_world;
        // A zero moment means no freedom about that axis — a thin rod's own length —
        // rather than an infinite acceleration.
        let d = |num: f64, i: f64| if i > 0.0 { num / i } else { 0.0 };
        let domega = DVec3::new(
            d(torque_body.x + (i2 - i3) * w.y * w.z, i1),
            d(torque_body.y + (i3 - i1) * w.z * w.x, i2),
            d(torque_body.z + (i1 - i2) * w.x * w.y, i3),
        );

        // q̇ = ½ q ⊗ ω, with ω in the body frame as a pure quaternion.
        let q = s.orientation;
        let wq = DQuat::from_xyzw(w.x, w.y, w.z, 0.0);
        let qdot = q.mul_quat(wq);
        Spin {
            omega: domega,
            orientation: DQuat::from_xyzw(0.5 * qdot.x, 0.5 * qdot.y, 0.5 * qdot.z, 0.5 * qdot.w),
        }
    }
}

impl Domain for RigidBody {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// A two-hundredth of the rotation period.
    ///
    /// An accuracy choice like the others here, and measured rather than guessed. At
    /// fifty steps a cycle the world-frame angular momentum drifts about `3e-7` of
    /// itself per second of tumbling; fourth order means quadrupling the count divides
    /// that by 256, and two hundred brings it under `1e-9` where the audit can be
    /// meaningfully tight.
    ///
    /// It is worth noticing that this is far weaker than
    /// [`NBody`](crate::NBody)'s momentum, which holds to `1e-13` over thousands of
    /// steps. The difference is structural rather than a matter of effort: linear
    /// momentum there is exact *by construction*, because equal and opposite forces
    /// cancel bit for bit. Nothing here cancels. Angular momentum is only as good as
    /// the integrator, plus a quaternion renormalisation each step that is a
    /// projection and not a conservation law.
    fn max_stable_dt(&self, _now: Time) -> Time {
        let rate = self.omega_body.length();
        if rate <= 0.0 {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(std::f64::consts::TAU / (200.0 * rate))
    }

    fn step(&mut self, t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let next = Integrator::Rk4.step(self, &self.spin_state(), t, dt);
        self.omega_body = next.omega;
        // Back onto the unit sphere. Integrating the kinematics numerically walks off
        // it, and an orientation of norm 1.0001 is not one.
        self.orientation = next.orientation.normalize();
        Ok(())
    }

    /// Angular momentum, per world axis.
    ///
    /// Only conserved when no torque is applied. A driven body's angular momentum grows
    /// at exactly `τ`, which is the definition of a torque rather than a leak — so the
    /// audit is meaningful for a free body and expected to fire for a driven one.
    fn ledger(&self) -> Ledger {
        let l = self.angular_momentum();
        Ledger::new()
            .with(conserved::ANGULAR_MOMENTUM_X, l.x)
            .with(conserved::ANGULAR_MOMENTUM_Y, l.y)
            .with(conserved::ANGULAR_MOMENTUM_Z, l.z)
    }

    fn checkpoint(&mut self) {
        self.saved = Some((self.omega_body, self.orientation));
    }

    fn restore(&mut self) {
        if let Some((w, q)) = self.saved {
            self.omega_body = w;
            self.orientation = q;
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

/// Precession rate of a symmetric top's angular velocity about its symmetry axis, in
/// the body frame: `ω₃(I₃ - I₁)/I₁`.
///
/// A closed form, and the only free-rotation case that has one. Positive for an oblate
/// body (a coin), negative for a prolate one (a rugby ball), and zero for a sphere.
pub fn symmetric_top_precession(inertia: Inertia, omega_body: DVec3) -> f64 {
    let p = inertia.principal;
    // The symmetry axis is the odd one out; assume z if x and y match.
    if (p.x - p.y).abs() / p.max_element().max(1e-300) < 1e-12 {
        if p.x <= 0.0 {
            return 0.0;
        }
        omega_body.z * (p.z - p.x) / p.x
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dualis_core::{audit, Schedule, Simulation};

    fn run(body: &mut RigidBody, duration: Time, steps: u32) {
        let dt = duration / steps as f64;
        let mut bus = Exchange::new();
        let mut t = Time::ZERO;
        for _ in 0..steps {
            body.step(t, dt, &mut bus).unwrap();
            t += dt;
        }
    }

    /// The moments of inertia against the textbook, which is where they come from and
    /// the only honest place to check them.
    #[test]
    fn the_standard_shapes_have_their_textbook_moments() {
        let m = Mass::kg(2.0);
        let r = Length::m(0.5);

        let sphere = Inertia::solid_sphere(m, r);
        assert!((sphere.principal.x - 0.4 * 2.0 * 0.25).abs() < 1e-15);
        assert!(sphere.is_spherical_top());

        // A shell has more inertia than a solid of the same mass, by 5/3.
        let shell = Inertia::hollow_sphere(m, r);
        assert!((shell.principal.x / sphere.principal.x - 5.0 / 3.0).abs() < 1e-12);

        // A cuboid's three moments differ, which is what lets it tumble.
        let brick = Inertia::solid_box(m, DVec3::new(0.1, 0.3, 0.6));
        assert!(!brick.is_symmetric_top());
        assert!((brick.principal.x - 2.0 / 12.0 * (0.09 + 0.36)).abs() < 1e-15);
        // The long axis is the easiest to spin about.
        assert!(brick.principal.z < brick.principal.y);
        assert!(brick.principal.y < brick.principal.x);

        // A cylinder is a symmetric top: two moments equal, one different.
        let can = Inertia::solid_cylinder(m, r, Length::m(1.0));
        assert!(can.is_symmetric_top() && !can.is_spherical_top());
        assert!((can.principal.z - 0.5 * 2.0 * 0.25).abs() < 1e-15);

        // A thin rod has no inertia about its own length.
        let rod = Inertia::rod(m, Length::m(1.2));
        assert_eq!(rod.principal.z, 0.0);
        assert!((rod.principal.x - 2.0 * 1.44 / 12.0).abs() < 1e-15);
    }

    /// A sphere spins about whatever axis it was given and never leaves it, because
    /// Euler's equations have no right-hand side when the three moments agree.
    #[test]
    fn a_sphere_cannot_tumble() {
        let axis = DVec3::new(0.3, -0.7, 0.5).normalize();
        let mut body = RigidBody::new("ball", Inertia::solid_sphere(Mass::kg(1.0), Length::m(0.2)))
            .spinning_body_frame(axis * 4.0);
        let before = body.angular_momentum();
        run(&mut body, Time::s(20.0), 4000);

        // The body-frame angular velocity has not moved at all.
        let after_body = body.angular_velocity_body();
        assert!(
            (after_body - axis * 4.0).length() / 4.0 < 1e-12,
            "a spherical top's body-frame omega is a constant, got {after_body}"
        );
        // And the world-frame angular momentum is exactly where it started.
        assert!(
            (body.angular_momentum() - before).length() / before.length() < 1e-12,
            "angular momentum should not have moved"
        );
    }

    /// A symmetric top precesses at exactly `ω₃(I₃-I₁)/I₁` in its own frame — the one
    /// closed form free rotation offers, and the thing to check the integrator against.
    ///
    /// The transverse part of `ω` goes round a circle at that rate while `ω₃` stays
    /// put, so after one precession period the body-frame `ω` is back where it began.
    #[test]
    fn a_symmetric_top_precesses_at_the_closed_form_rate() {
        let inertia = Inertia::solid_cylinder(Mass::kg(3.0), Length::m(0.4), Length::m(0.2));
        // Oblate: a short fat cylinder, so I3 > I1 and the precession is positive.
        assert!(inertia.principal.z > inertia.principal.x);

        let omega0 = DVec3::new(0.5, 0.0, 6.0);
        let rate = symmetric_top_precession(inertia, omega0);
        assert!(rate > 0.0, "an oblate top precesses forwards, got {rate}");
        // omega3 (I3 - I1) / I1, worked out by hand from the moments.
        let (i1, i3) = (inertia.principal.x, inertia.principal.z);
        assert!((rate - 6.0 * (i3 - i1) / i1).abs() < 1e-12);

        let period = Time::s(std::f64::consts::TAU / rate.abs());
        let mut body = RigidBody::new("top", inertia).spinning_body_frame(omega0);

        // A quarter of the precession period should have turned the transverse part of
        // omega by ninety degrees, from +x to +y.
        run(&mut body, period / 4.0, 4000);
        let w = body.angular_velocity_body();
        assert!(
            (w.z - 6.0).abs() < 1e-9,
            "the axial component is constant, got {}",
            w.z
        );
        assert!(
            w.x.abs() < 5e-3 && (w.y - 0.5).abs() < 5e-3,
            "a quarter period should take omega from +x to +y, got {w}"
        );

        // And after a whole period it is back.
        let mut body = RigidBody::new("top", inertia).spinning_body_frame(omega0);
        run(&mut body, period, 16000);
        let w = body.angular_velocity_body();
        assert!(
            (w - omega0).length() / omega0.length() < 1e-4,
            "one precession period should close the circle, got {w}"
        );
    }

    /// The intermediate axis theorem, which is the reason a rigid body is worth
    /// modelling at all rather than being a point with a spin attached.
    ///
    /// Rotation about the largest or smallest principal axis is stable: a small
    /// disturbance stays small. About the middle one it is not — the disturbance grows
    /// exponentially until the body flips over, then grows again, forever. A tennis
    /// racket does it, and so does a wing nut in orbit.
    ///
    /// Nothing about that is obvious from the equations until you perturb them, which
    /// is exactly why it belongs in a test.
    #[test]
    fn rotation_about_the_intermediate_axis_is_unstable() {
        // Three clearly different moments.
        let inertia = Inertia::solid_box(Mass::kg(1.0), DVec3::new(0.1, 0.4, 0.9));
        let (small, middle, large) = inertia.sorted_axes();
        assert_eq!(inertia.intermediate_axis(), middle);
        assert_ne!(small, large);

        // Spin about an axis with a one-percent wobble on the next one round, and track
        // the *signed* projection of the spin onto the axis it started on. A flip is a
        // change of sign; a magnitude would saturate at "perpendicular" and could not
        // tell a body that toppled from one that merely wandered.
        let spin_about = |axis: usize| {
            let mut w = DVec3::ZERO;
            w[axis] = 10.0;
            w[(axis + 1) % 3] = 0.1;
            let mut body = RigidBody::new("brick", inertia).spinning_body_frame(w);
            let mut lowest = 1.0f64;
            let dt = Time::ms(0.5);
            let mut bus = Exchange::new();
            let mut t = Time::ZERO;
            for _ in 0..40_000 {
                body.step(t, dt, &mut bus).unwrap();
                t += dt;
                let now = body.angular_velocity_body();
                lowest = lowest.min(now[axis] / now.length());
            }
            lowest
        };

        // The stable axes hold: the spin never leaves the axis it started on by more
        // than the wobble it was given.
        for (name, axis) in [("smallest", small), ("largest", large)] {
            let lowest = spin_about(axis);
            assert!(
                lowest > 0.99,
                "the {name} axis should be stable, but the spin fell to {lowest:.4} of \
                 it"
            );
        }
        // The middle one turns the body right over: the spin ends up pointing the
        // opposite way along the axis it started on.
        let lowest = spin_about(middle);
        assert!(
            lowest < -0.9,
            "the intermediate axis should flip the body over, but the spin only reached \
             {lowest:.4} of it"
        );
    }

    /// Both scalars a free body conserves, checked against the integrator rather than
    /// assumed: `|L|` and the rotational energy. Runge-Kutta dissipates, so the figures
    /// are what fourth order at fifty steps a cycle buys.
    #[test]
    fn a_free_body_conserves_its_angular_momentum_and_energy() {
        let inertia = Inertia::solid_box(Mass::kg(2.0), DVec3::new(0.2, 0.5, 0.8));
        let mut body =
            RigidBody::new("brick", inertia).spinning_body_frame(DVec3::new(3.0, 7.0, 2.0));
        let l0 = body.angular_momentum();
        let e0 = body.rotational_energy().to_si();

        run(&mut body, Time::s(60.0), 60_000);

        let l = body.angular_momentum();
        assert!(
            (l - l0).length() / l0.length() < 1e-9,
            "angular momentum drifted by {:e}",
            (l - l0).length() / l0.length()
        );
        let e = body.rotational_energy().to_si();
        assert!(
            (e / e0 - 1.0).abs() < 1e-9,
            "rotational energy drifted by {:e}",
            (e / e0 - 1.0).abs()
        );
        // The orientation stays a rotation, which is what the renormalisation is for.
        assert!((body.orientation().length() - 1.0).abs() < 1e-12);
    }

    /// Angular momentum is constant in the **world** frame and moves in the body frame.
    /// Reporting the body-frame vector as conserved is the classic error, and it would
    /// pass a careless test on a symmetric top while failing on anything that tumbles.
    #[test]
    fn angular_momentum_is_constant_in_the_world_frame_only() {
        let inertia = Inertia::solid_box(Mass::kg(1.0), DVec3::new(0.15, 0.45, 0.75));
        let mut body =
            RigidBody::new("brick", inertia).spinning_body_frame(DVec3::new(1.0, 6.0, 1.5));
        let world0 = body.angular_momentum();
        let body0 = inertia.principal * body.angular_velocity_body();

        run(&mut body, Time::s(3.0), 6000);

        let world = body.angular_momentum();
        assert!(
            (world - world0).length() / world0.length() < 1e-9,
            "the world-frame vector is the conserved one"
        );
        let body_now = inertia.principal * body.angular_velocity_body();
        assert!(
            (body_now - body0).length() / body0.length() > 0.1,
            "the body-frame vector should have moved appreciably, and did not: \
             {body_now} against {body0}"
        );
        // Their magnitudes agree, since a rotation does not change a length. That is
        // exactly what makes the mistake survive a magnitude-only check.
        assert!((body_now.length() / world.length() - 1.0).abs() < 1e-9);
    }

    /// The kernel audits angular momentum with the same machinery as energy and linear
    /// momentum, and the tolerance is the integrator's.
    #[test]
    fn the_kernel_audits_angular_momentum() {
        let inertia = Inertia::solid_box(Mass::kg(1.5), DVec3::new(0.2, 0.4, 0.7));
        let body = RigidBody::new("brick", inertia).spinning_body_frame(DVec3::new(2.0, 5.0, 1.0));
        let before = body.ledger();

        // Two things set this tolerance, and both are worth knowing.
        //
        // The integrator: two hundred steps a cycle drifts the largest component by
        // 1.3e-9 of itself per window. That measurement doubles as a check on the
        // order — at fifty steps a cycle the same run drifted 2.9e-7, and quadrupling
        // the count cut it by 220, which is the 256 fourth order predicts.
        //
        // And the audit's own shape: a vector is audited component by component, so the
        // *smallest* component is the binding constraint. The absolute error is set by
        // the whole of `|L|`, while the scale it is measured against is only that
        // component — here `L_z` is a sixth of `L_x`, and its relative drift is
        // correspondingly larger. The 2.2e-8 below is `L_z`'s, not the vector's.
        let mut sim = Simulation::new(Schedule::Multirate)
            .conservation_tolerance(1e-7)
            .with(body);
        for _ in 0..50 {
            sim.advance(Time::s(0.2))
                .expect("a free body's angular momentum must hold");
        }
        audit("rigid", &before, &sim.ledger(), 1e-5).expect("and hold end to end");

        // A ledger with a deliberately wrong axis is caught, so the check has teeth.
        let get = |q: &str| before.get(q).unwrap();
        let broken = Ledger::new()
            .with(
                conserved::ANGULAR_MOMENTUM_X,
                get(conserved::ANGULAR_MOMENTUM_X) * 1.001,
            )
            .with(
                conserved::ANGULAR_MOMENTUM_Y,
                get(conserved::ANGULAR_MOMENTUM_Y),
            )
            .with(
                conserved::ANGULAR_MOMENTUM_Z,
                get(conserved::ANGULAR_MOMENTUM_Z),
            );
        let err = audit("rigid", &before, &broken, 1e-5).expect_err("0.1% is not rounding");
        assert_eq!(err.quantity, "angular_momentum_x");
    }

    /// Frames are given explicitly and converted correctly: the same physical spin
    /// described either way produces the same angular momentum.
    #[test]
    fn body_and_world_frames_describe_the_same_spin() {
        let inertia = Inertia::solid_box(Mass::kg(1.0), DVec3::new(0.2, 0.4, 0.6));
        let tilt = DQuat::from_axis_angle(DVec3::new(1.0, 1.0, 0.0).normalize(), 0.7);
        let world_omega = DVec3::new(0.0, 0.0, 5.0);

        let by_world = RigidBody::new("a", inertia)
            .with_orientation(tilt)
            .spinning(world_omega);
        assert!(
            (by_world.angular_velocity() - world_omega).length() < 1e-12,
            "a world-frame spin should read back in the world frame"
        );

        let by_body = RigidBody::new("b", inertia)
            .with_orientation(tilt)
            .spinning_body_frame(tilt.inverse() * world_omega);
        assert!(
            (by_body.angular_momentum() - by_world.angular_momentum()).length() < 1e-12,
            "the same spin, described either way, is the same spin"
        );
        // A tilted body's angular momentum is not parallel to its angular velocity
        // unless it is a spherical top — which is the reason a wheel out of balance
        // shakes its bearing.
        let l = by_world.angular_momentum().normalize();
        let w = by_world.angular_velocity().normalize();
        assert!(
            l.dot(w) < 0.999,
            "L and omega should not be parallel for an asymmetric body: {}",
            l.dot(w)
        );
    }

    /// A body at rest reports no limit rather than a zero step, and a rod with no
    /// inertia about its length does not divide by it.
    #[test]
    fn degenerate_bodies_are_handled() {
        let still = RigidBody::new(
            "still",
            Inertia::solid_sphere(Mass::kg(1.0), Length::m(1.0)),
        );
        assert!(!still.max_stable_dt(Time::ZERO).to_si().is_finite());
        assert_eq!(still.angular_momentum(), DVec3::ZERO);
        assert_eq!(still.rotational_energy().to_si(), 0.0);

        let mut rod = RigidBody::new("rod", Inertia::rod(Mass::kg(1.0), Length::m(2.0)))
            .spinning_body_frame(DVec3::new(1.0, 0.0, 3.0));
        let mut bus = Exchange::new();
        for _ in 0..1000 {
            rod.step(Time::ZERO, Time::ms(1.0), &mut bus).unwrap();
        }
        assert!(
            rod.angular_velocity_body().is_finite(),
            "a zero moment is no freedom, not an infinite acceleration"
        );
        // The rod has no angular momentum about its own length, whatever it was given.
        assert!(rod.angular_momentum().is_finite());
    }
}
