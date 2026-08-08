//! What `Schedule::Multirate` does with a coupled quantity, and what the audit cannot see.
//!
//! Found from outside, by an application that compared a lumped plate under a lamp against the
//! closed form of its own scheme. The total energy is right to 1e-12 every step; the *time* at
//! which it arrives is not, and `Ledger` has no representation for when.

use dualis::prelude::*;

/// A steady source with a finite tank, quasi-static: it publishes once per outer step.
struct Lamp {
    watts: f64,
    reserve: f64,
}

impl Domain for Lamp {
    fn name(&self) -> &str {
        "lamp"
    }
    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }
    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let j = (self.watts * dt.to_si()).min(self.reserve);
        self.reserve -= j;
        bus.publish(HEAT, j);
        Ok(())
    }
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.reserve)
    }
}

fn plate() -> LumpedMass {
    LumpedMass::new(
        "plate",
        Substance::aluminium_6061(),
        Volume::from_si(60e-3 * 60e-3 * 3e-3),
        Length::mm(1.5),
        Temperature::celsius(20.0),
        Environment::still_air(
            Temperature::celsius(20.0),
            Area::from_si(2.0 * 60e-3 * 60e-3),
        ),
    )
}

/// Run to 600 s in `outer` steps and report the rise.
fn rise_after(schedule: Schedule, outer: Time, steps: usize) -> f64 {
    let mut sim = Simulation::new(schedule)
        .conservation_tolerance(1e-6)
        .with(Lamp {
            watts: 2.0,
            reserve: f64::INFINITY,
        })
        .with(plate());
    for _ in 0..steps {
        sim.advance(outer).expect("the books close");
    }
    sim.domain_as::<LumpedMass>("plate")
        .expect("the plate is there")
        .rise()
        .to_si()
}

/// **Subcycling buys nothing across a coupling, and the audit stays green throughout.**
///
/// `Simulation::sweep` steps one domain to completion before the next. A quasi-static source
/// is not subcycled, so it publishes a whole outer step's joules once; the consumer then
/// subcycles, and `Exchange::take` on its *first* substep empties the channel. Every joule of
/// the interval is therefore deposited at its beginning and then decays for the rest of it.
///
/// So refining the substep does not refine the answer. Taking the limit of the recursion
/// `u ← u·gⁿ + (P·dt/C)·g^(n−1)` with `g = 1 − h/τ` as `n → ∞` gives
/// `u·e^{−dt/τ} + (P·dt/C)·e^{−dt/τ}`, which is not the solution — it is first order in the
/// **outer** step with no dependence on the substep at all.
///
/// This test pins the consequence rather than the mechanism: at a fixed outer step, `Multirate`
/// is not reliably better than `Staggered`, and at a coarse enough step it is worse. Both pass
/// the conservation audit with a residual around 1e-12, because the total that crossed is
/// exactly right and only its distribution in time is wrong. That is the time-domain twin of
/// the reason `Exchange::audit_transfers` had to become a per-face check in space.
#[test]
fn multirate_does_not_improve_a_coupled_quantity_and_the_audit_cannot_tell() {
    // The closed form of the continuous problem, for reference.
    let p = plate();
    let c = p.heat_capacity().to_si();
    let tau = p.time_constant().to_si();
    let ha = c / tau;
    let analytic = 2.0 / ha * (1.0 - (-600.0 / tau).exp());

    // At a coarse outer step, multirate subcycles the plate and staggered does not.
    let coarse = Time::from_si(300.0);
    let (stag, multi) = (
        rise_after(Schedule::Staggered, coarse, 2),
        rise_after(Schedule::Multirate, coarse, 2),
    );
    // Both audit clean. Neither is close, and the subcycling one is *further* away.
    assert!(
        (multi - analytic).abs() > (stag - analytic).abs(),
        "multirate should be the worse of the two at a 300 s outer step: \
         staggered {stag:.6} K, multirate {multi:.6} K, analytic {analytic:.6} K"
    );

    // And refining only the outer step is what helps: halving it roughly halves the error,
    // which is first order in the *outer* step even though the substep shrank too.
    let errs: Vec<f64> = [(300.0, 2usize), (150.0, 4), (75.0, 8)]
        .iter()
        .map(|(dt, n)| (rise_after(Schedule::Multirate, Time::from_si(*dt), *n) - analytic).abs())
        .collect();
    for w in errs.windows(2) {
        let ratio = w[0] / w[1];
        assert!(
            (1.6..2.6).contains(&ratio),
            "halving the outer step should roughly halve the error: {:?}, ratio {ratio:.3}",
            errs
        );
    }
}
