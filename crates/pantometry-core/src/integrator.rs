//! Time evolution for systems that have no closed form.
//!
//! [`Motion`](crate::motion::Motion) is a function of `t`: ask for the world at
//! 0.7 s and you get it, without having computed 0.6 s first. That is worth a
//! great deal — an exposure can be sampled at seven instants for motion blur, and
//! frame 7 of a recording does not depend on having rendered frame 6 — and it is
//! why drift, oscillation and spin are written the way they are.
//!
//! It is also not available in general. Three bodies under gravity have no closed
//! form, and neither do contact, heat conduction, or a stiff reaction network.
//! Those systems have to be rolled forward, and frame 7 genuinely does depend on
//! frame 6.
//!
//! # Reproducibility survives the trade, under three rules
//!
//! - **Fixed steps only.** [`Integrator::step`] takes `dt` and uses it. An
//!   adaptive step chosen from the local error makes the floating-point path
//!   depend on the values, so two runs that should agree diverge at the first
//!   place one of them decided to halve the step. Where stability demands a
//!   smaller step, take a fixed *number* of substeps — see
//!   [`substeps_for`].
//! - **No wall clock.** Nothing here reads a timer.
//! - **Ordered reduction.** Summing forces in parallel changes the answer,
//!   because floating-point addition is not associative. That rule belongs to the
//!   domains, but it is the reason [`State::axpy`] is a sequential operation on a
//!   whole state rather than a per-element one to be farmed out.
//!
//! # Symplectic versus accurate
//!
//! [`Integrator::Rk4`] is fourth-order accurate and loses energy steadily.
//! [`velocity_verlet`] is second-order and does not: its energy error oscillates
//! within a bound instead of drifting, because it preserves the geometric
//! structure of a Newtonian system rather than merely fitting its derivative.
//! Over ten steps RK4 wins; over ten million, it has quietly cooled the system
//! down. For anything conservative — orbits, molecules, an undamped spring — use
//! the symplectic one and let [`crate::conserved::audit`] confirm it.
//!
//! The test module proves exactly this on a harmonic oscillator, against the
//! closed-form energy.

use pantometry_units::Time;

/// A state vector that an integrator can do arithmetic on.
///
/// Deliberately tiny: an integrator needs to scale a state and add a multiple of
/// another to it, and nothing else. Anything a domain wants to keep alongside its
/// numbers — a mesh, a material table, a name — stays out of the state and lives
/// on the [`Dynamics`] instead, where it costs nothing to carry.
pub trait State: Clone {
    /// `self += a * other`, the one operation every explicit integrator is made
    /// of. Implementations must visit their elements in a fixed order.
    fn axpy(&mut self, a: f64, other: &Self);

    /// `self *= a`.
    fn scale(&mut self, a: f64);

    /// A zero of the same shape as this state.
    fn zeros_like(&self) -> Self;
}

/// A first-order system: `ds/dt = f(s, t)`.
///
/// The derivative has a different dimension from the state, which is a thing this
/// crate's unit types cannot express through a trait — so [`State`] is in raw SI
/// numbers and the dimensions live in the domain's own types on either side.
pub trait Dynamics {
    /// The state this system evolves.
    type S: State;

    /// `f(s, t)`. Must be a pure function: two calls with the same arguments have
    /// to give the same answer, or none of the reproducibility above holds.
    fn derivative(&self, s: &Self::S, t: Time) -> Self::S;
}

/// A Newtonian system: `d²x/dt² = a(x, t)`, with no dependence on velocity.
///
/// The restriction is what buys the symplectic integrator. A velocity-dependent
/// force — drag, friction, a magnetic field — is not conservative, so there is no
/// energy for a symplectic method to preserve; express those through [`Dynamics`]
/// and [`Integrator::Rk4`] instead.
pub trait Newtonian {
    /// The configuration space — positions, not positions and velocities.
    type Coords: State;

    /// `a(x, t)`. Must not depend on velocity; see the note on this trait for why.
    fn acceleration(&self, x: &Self::Coords, t: Time) -> Self::Coords;
}

/// Explicit fixed-step integrators.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Integrator {
    /// First order. Cheap, and wrong fast enough to be visible — useful mainly as
    /// the thing a better integrator is compared against.
    Euler,
    /// Second-order midpoint.
    Midpoint,
    /// Classical fourth-order Runge-Kutta. Accurate per step, and it dissipates:
    /// see the module docs before using it on anything conservative.
    Rk4,
}

impl Integrator {
    /// One step of `dt`, from `t`.
    pub fn step<D: Dynamics>(&self, system: &D, s: &D::S, t: Time, dt: Time) -> D::S {
        let h = dt.to_si();
        match self {
            Integrator::Euler => {
                let mut next = s.clone();
                next.axpy(h, &system.derivative(s, t));
                next
            }
            Integrator::Midpoint => {
                let k1 = system.derivative(s, t);
                let mut mid = s.clone();
                mid.axpy(h / 2.0, &k1);
                let k2 = system.derivative(&mid, t + dt / 2.0);
                let mut next = s.clone();
                next.axpy(h, &k2);
                next
            }
            Integrator::Rk4 => {
                let k1 = system.derivative(s, t);
                let mut y = s.clone();
                y.axpy(h / 2.0, &k1);
                let k2 = system.derivative(&y, t + dt / 2.0);
                let mut y = s.clone();
                y.axpy(h / 2.0, &k2);
                let k3 = system.derivative(&y, t + dt / 2.0);
                let mut y = s.clone();
                y.axpy(h, &k3);
                let k4 = system.derivative(&y, t + dt);

                // (k1 + 2 k2 + 2 k3 + k4) / 6, accumulated in one fixed order.
                let mut slope = k1;
                slope.axpy(2.0, &k2);
                slope.axpy(2.0, &k3);
                slope.axpy(1.0, &k4);
                slope.scale(1.0 / 6.0);
                let mut next = s.clone();
                next.axpy(h, &slope);
                next
            }
        }
    }

    /// `n` steps of `dt`, which is how a domain subcycles inside a larger step.
    pub fn advance<D: Dynamics>(&self, system: &D, s: &D::S, t: Time, dt: Time, n: u32) -> D::S {
        let mut state = s.clone();
        let mut now = t;
        for _ in 0..n {
            state = self.step(system, &state, now, dt);
            now += dt;
        }
        state
    }
}

/// One velocity-Verlet step, in place. Symplectic, second order, and the right
/// default for anything whose energy is supposed to stay put.
///
/// The classic kick-drift-kick form: half a kick from the acceleration where we
/// are, a full drift, then half a kick from the acceleration where we arrived. It
/// needs one acceleration evaluation per step, since the second half-kick's value
/// is reused as the next step's first.
pub fn velocity_verlet<N: Newtonian>(
    system: &N,
    x: &mut N::Coords,
    v: &mut N::Coords,
    t: Time,
    dt: Time,
) {
    let h = dt.to_si();
    let a0 = system.acceleration(x, t);
    v.axpy(h / 2.0, &a0);
    x.axpy(h, v);
    let a1 = system.acceleration(x, t + dt);
    v.axpy(h / 2.0, &a1);
}

/// How many equal substeps of at most `limit` it takes to cover `dt`.
///
/// The deterministic answer to a stability limit: rather than shrinking the step
/// to whatever the local error asks for, take an integer number of equal ones.
/// Two runs of the same scene take the same count, so the arithmetic follows the
/// same path.
pub fn substeps_for(dt: Time, limit: Time) -> u32 {
    let (dt, limit) = (dt.to_si(), limit.to_si());
    // An infinite limit means "no limit", a NaN one means the domain does not know,
    // and both answer 1 — as does a step that is not going anywhere.
    if !limit.is_finite() || limit <= 0.0 || dt <= 0.0 {
        return 1;
    }
    let n = (dt / limit).ceil();
    if !n.is_finite() || n < 1.0 {
        1
    } else {
        n.min(u32::MAX as f64) as u32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A one-dimensional harmonic oscillator: `x'' = -x`, whose energy
    /// `(x² + v²)/2` is exactly constant and whose solution is a cosine. Every
    /// claim below is checked against that closed form rather than against
    /// another integrator.
    #[derive(Clone, Debug, PartialEq)]
    struct Pair(f64, f64);

    impl State for Pair {
        fn axpy(&mut self, a: f64, other: &Self) {
            self.0 += a * other.0;
            self.1 += a * other.1;
        }
        fn scale(&mut self, a: f64) {
            self.0 *= a;
            self.1 *= a;
        }
        fn zeros_like(&self) -> Self {
            Pair(0.0, 0.0)
        }
    }

    /// As a first-order system: state is (x, v), derivative is (v, -x).
    struct Spring;

    impl Dynamics for Spring {
        type S = Pair;
        fn derivative(&self, s: &Pair, _t: Time) -> Pair {
            Pair(s.1, -s.0)
        }
    }

    /// As a Newtonian system: coordinate is x, acceleration is -x.
    #[derive(Clone, Debug)]
    struct Scalar(f64);

    impl State for Scalar {
        fn axpy(&mut self, a: f64, other: &Self) {
            self.0 += a * other.0;
        }
        fn scale(&mut self, a: f64) {
            self.0 *= a;
        }
        fn zeros_like(&self) -> Self {
            Scalar(0.0)
        }
    }

    impl Newtonian for Spring {
        type Coords = Scalar;
        fn acceleration(&self, x: &Scalar, _t: Time) -> Scalar {
            Scalar(-x.0)
        }
    }

    fn energy(x: f64, v: f64) -> f64 {
        (x * x + v * v) / 2.0
    }

    /// Order of accuracy, measured rather than asserted: halving the step should
    /// cut Euler's error by 2, the midpoint rule's by 4 and RK4's by 16.
    #[test]
    fn each_integrator_shows_its_order() {
        let exact = |t: f64| t.cos();
        let error_at = |method: Integrator, steps: u32| {
            let dt = Time::s(1.0 / steps as f64);
            let end = method.advance(&Spring, &Pair(1.0, 0.0), Time::ZERO, dt, steps);
            (end.0 - exact(1.0)).abs()
        };
        for (method, expected_ratio) in [
            (Integrator::Euler, 2.0),
            (Integrator::Midpoint, 4.0),
            (Integrator::Rk4, 16.0),
        ] {
            let coarse = error_at(method, 200);
            let fine = error_at(method, 400);
            let ratio = coarse / fine;
            assert!(
                (ratio - expected_ratio).abs() / expected_ratio < 0.15,
                "{method:?}: halving the step changed the error by {ratio:.2}, \
                 expected about {expected_ratio}"
            );
        }
    }

    /// The reason the symplectic integrator exists. Over 200 000 steps RK4 —
    /// which is *more* accurate per step — has visibly drained the oscillator,
    /// while velocity-Verlet's energy error is still bounded and oscillating.
    ///
    /// This is the whole argument for `velocity_verlet` in one assertion, and it
    /// is why "use the higher-order method" is the wrong instinct for anything
    /// whose conservation is being audited.
    #[test]
    fn only_the_symplectic_integrator_keeps_its_energy() {
        // The step has to be a real fraction of the period for RK4's dissipation to
        // show at all: its energy loss per step goes as dt⁶, so at dt = 0.05 the
        // drift over these many steps is only 4e-5 and the comparison would be
        // measuring nothing. At dt = 0.2 — still thirty steps per period — it is
        // visible, and that is the honest regime in which the two differ.
        const STEPS: u32 = 200_000;
        let dt = Time::s(0.2);

        let rk4 = Integrator::Rk4.advance(&Spring, &Pair(1.0, 0.0), Time::ZERO, dt, STEPS);
        let rk4_drift = (energy(rk4.0, rk4.1) - 0.5).abs() / 0.5;

        let (mut x, mut v) = (Scalar(1.0), Scalar(0.0));
        let mut t = Time::ZERO;
        let mut worst = 0.0f64;
        for _ in 0..STEPS {
            velocity_verlet(&Spring, &mut x, &mut v, t, dt);
            t += dt;
            // Verlet's velocity is half a step out of phase with its position, so
            // the energy read from the pair wobbles by O(dt²) — bounded, which is
            // the property being tested, and not drifting.
            worst = worst.max((energy(x.0, v.0) - 0.5).abs() / 0.5);
        }

        assert!(
            rk4_drift > 0.05,
            "RK4 should have leaked energy over {STEPS} steps, drift {rk4_drift:.3e}"
        );
        assert!(
            worst < 0.05,
            "velocity-Verlet's energy error should stay bounded, worst {worst:.3e}"
        );
        assert!(
            worst < rk4_drift / 5.0,
            "the symplectic method should hold energy far better: {worst:.3e} vs {rk4_drift:.3e}"
        );
    }

    /// Integration is reproducible to the last bit, including when it is done in
    /// two halves instead of one run.
    #[test]
    fn integration_is_bit_reproducible() {
        let dt = Time::s(0.01);
        let once = Integrator::Rk4.advance(&Spring, &Pair(1.0, 0.0), Time::ZERO, dt, 500);
        let twice = {
            let half = Integrator::Rk4.advance(&Spring, &Pair(1.0, 0.0), Time::ZERO, dt, 250);
            Integrator::Rk4.advance(&Spring, &half, Time::s(2.5), dt, 250)
        };
        assert_eq!(once, twice, "restarting mid-run must change nothing");

        let again = Integrator::Rk4.advance(&Spring, &Pair(1.0, 0.0), Time::ZERO, dt, 500);
        assert_eq!(once, again);
    }

    /// Substep counts are integers derived from a fixed limit, never from the
    /// local error — that is what keeps two runs on the same arithmetic path.
    #[test]
    fn substep_counts_are_deterministic_integers() {
        assert_eq!(substeps_for(Time::s(1.0), Time::s(0.3)), 4);
        assert_eq!(substeps_for(Time::s(1.0), Time::s(0.5)), 2);
        assert_eq!(substeps_for(Time::s(1.0), Time::s(2.0)), 1);
        // A quasi-static domain reports no limit at all.
        assert_eq!(substeps_for(Time::s(1.0), Time::from_si(f64::INFINITY)), 1);
        // Degenerate limits do not produce a zero or an overflowing count.
        assert_eq!(substeps_for(Time::s(1.0), Time::ZERO), 1);
        assert_eq!(substeps_for(Time::s(1.0), Time::s(-1.0)), 1);
        assert_eq!(substeps_for(Time::ZERO, Time::s(0.1)), 1);
        // And covering the step is guaranteed: n * limit >= dt.
        for (dt, limit) in [(1.0, 0.3), (7.3, 0.11), (1e-3, 1e-7)] {
            let n = substeps_for(Time::s(dt), Time::s(limit));
            assert!(n as f64 * limit >= dt - 1e-12, "{n} x {limit} < {dt}");
        }
    }
}
