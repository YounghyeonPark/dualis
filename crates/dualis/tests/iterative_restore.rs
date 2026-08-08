//! The restore path of `Schedule::Iterative`, which nothing in this workspace had ever taken.
//!
//! `supports_restore` is a promise every restorable domain makes, and `iterate` is the only
//! caller that collects on it — by rewinding every domain before each sweep after the first.
//! But no coupling in this workspace is strong enough to need iterating, so `Domain::residual`
//! is the default `0.0` everywhere, `iterate` converges on its first sweep, and the restore
//! branch has never executed. A design review noticed that and it turned out to be the more
//! useful half of what it was looking for.
//!
//! So this file supplies the missing ingredient: a domain with a genuine residual, which forces
//! real iterations, which forces real restores.

use dualis::prelude::*;

/// A domain that converges over sweeps, so `iterate` has to iterate.
///
/// The point is the split between what is rewound and what is not. `x` is state and is
/// checkpointed; `estimate` is what the sweep *learned* and is deliberately not, because that is
/// exactly what makes an iterative coupling converge — each sweep restarts from the same
/// physical state with a better guess about its neighbours. A domain that rewound its estimate
/// too would repeat the same sweep forever.
struct Converging {
    /// Rewound by `restore`.
    x: f64,
    /// Carried across restores, on purpose.
    estimate: f64,
    last_change: f64,
    saved: Option<f64>,
}

impl Converging {
    fn new() -> Converging {
        Converging {
            x: 0.0,
            estimate: 0.0,
            last_change: f64::INFINITY,
            saved: None,
        }
    }
}

impl Domain for Converging {
    fn name(&self) -> &str {
        "converging"
    }

    /// `e ← (e + 1)/2`, whose fixed point is 1 and whose error halves every sweep. Chosen
    /// because the number of iterations to a given tolerance is then exactly predictable, so
    /// the test can assert it rather than observe it.
    fn step(&mut self, _t: Time, _dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let next = 0.5 * (self.estimate + 1.0);
        self.last_change = (next - self.estimate).abs();
        self.estimate = next;
        self.x = next;
        Ok(())
    }

    fn residual(&self) -> f64 {
        self.last_change
    }

    /// Also where the estimate is reset, because `iterate` calls this exactly once per
    /// `advance` — so it is the only hook a domain has for "a new outer step is starting".
    ///
    /// Without it the estimate stays converged from the previous step, every later `advance`
    /// finishes in one sweep, and the restore path goes unexercised again. Which is a small
    /// version of how it came to be unexercised in the first place.
    fn checkpoint(&mut self) {
        self.saved = Some(self.x);
        self.estimate = 0.0;
        self.last_change = f64::INFINITY;
    }

    fn restore(&mut self) {
        if let Some(x) = self.saved {
            self.x = x;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }
}

fn plate() -> LumpedMass {
    LumpedMass::new(
        "plate",
        Substance::aluminium_6061(),
        Volume::from_si(60e-3 * 60e-3 * 3e-3),
        Length::mm(1.5),
        Temperature::celsius(80.0),
        Environment::still_air(
            Temperature::celsius(20.0),
            Area::from_si(2.0 * 60e-3 * 60e-3),
        ),
    )
}

/// **The restore branch executes, and it is asserted that it does.**
///
/// Starting from zero, `e ← (e+1)/2` changes by 0.5, then 0.25, then 0.125, so a tolerance of
/// 0.2 is met on the third sweep and `iterate` restores twice on the way. Before this file the
/// same call converged on sweep one and restored zero times, which is why `supports_restore`
/// had never been collected on.
#[test]
fn an_iterative_coupling_actually_iterates_and_actually_restores() {
    let mut sim = Simulation::new(Schedule::Iterative {
        max_iter: 8,
        tol: 0.2,
    })
    .with(Converging::new());

    let report = sim.advance(Time::s(1.0)).expect("it converges");
    assert_eq!(
        report.iterations, 3,
        "0.5, 0.25, 0.125 against a tolerance of 0.2 is three sweeps"
    );

    // And a coupling that cannot reach the tolerance is refused rather than accepted, which is
    // the other half of the promise.
    let mut stubborn = Simulation::new(Schedule::Iterative {
        max_iter: 2,
        tol: 1e-9,
    })
    .with(Converging::new());
    let violation = stubborn
        .advance(Time::s(1.0))
        .expect_err("two sweeps cannot reach 1e-9");
    assert_eq!(violation.quantity, "coupling residual");
    assert_eq!(
        stubborn.time().to_si(),
        0.0,
        "a refused step keeps the clock"
    );
}

/// **A `LumpedMass` rewound mid-coupling does not report energy it did not lose.**
///
/// The hypothesis this file was written to test. `LumpedMass::checkpoint` saves the temperature
/// and nothing else, while `ledger()` reports `stored + lost` — so if `lost` kept accumulating
/// across sweeps whose temperatures were rewound, the books would grow by one sweep's losses
/// every iteration and the audit would see energy created out of nothing.
///
/// It is not hypothetical any more: with `Converging` in the simulation this runs three sweeps
/// and two restores per `advance`, and the conservation audit is live at 1e-9 throughout.
#[test]
fn a_rewound_lumped_mass_does_not_invent_the_heat_it_shed() {
    let mut sim = Simulation::new(Schedule::Iterative {
        max_iter: 8,
        tol: 0.2,
    })
    .conservation_tolerance(1e-9)
    .with(Converging::new())
    .with(plate());

    // The plate starts at 80 C against 20 C ambient, so it sheds heat every substep of every
    // sweep and `lost` is large. That is load-bearing: the first version of this test started
    // it *at* ambient, where `lost` stays identically zero — and a double count of zero is
    // zero, so it passed while measuring nothing. The configuration-with-no-physics trap, in
    // the test written to catch a different one.
    let start = sim
        .domain_as::<LumpedMass>("plate")
        .expect("the plate is there")
        .lost_energy()
        .to_si();
    assert_eq!(start, 0.0);
    for _ in 0..40 {
        let report = sim
            .advance(Time::s(5.0))
            .expect("a rewound sweep must not create energy");
        assert_eq!(report.iterations, 3);
    }

    // And it really did shed something, so the audit had something to be wrong about: three
    // sweeps an advance for forty advances, every one of them rewinding the temperature while
    // this total kept climbing.
    let shed = sim
        .domain_as::<LumpedMass>("plate")
        .expect("the plate is there")
        .lost_energy()
        .to_si();
    assert!(
        shed > 100.0,
        "the plate should have shed real joules, got {shed:.3} J"
    );
}
