//! One tolerance for the whole simulation meant the loosest quantity set what every other one was
//! checked against.
//!
//! This is the failure, demonstrated: two quantities on one bus, one carried by a scheme that
//! cannot do better than `1e-6` and one carried by a scheme exact to `1e-15`. Under a single
//! number the run is either refused for being correct or passed while leaking.

use pantometry_core::conserved::quantity;
use pantometry_core::units::Time;
use pantometry_core::{
    audit, audit_with, Domain, Exchange, Kind, Ledger, Schedule, Simulation, Tolerances, Violation,
};

/// A domain holding two quantities, each of which it can be told to leak.
///
/// Nothing is published, so the bus never sees anything and the only thing under test is the
/// whole-simulation ledger audit — which is what the tolerance applies to.
struct Leaky {
    energy: f64,
    momentum: f64,
    energy_leak: f64,
    momentum_leak: f64,
}

impl Domain for Leaky {
    fn name(&self) -> &str {
        "leaky"
    }
    fn kind(&self) -> Kind {
        Kind::Evolving
    }
    fn step(&mut self, _t: Time, _dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        self.energy *= 1.0 - self.energy_leak;
        self.momentum *= 1.0 - self.momentum_leak;
        Ok(())
    }
    fn ledger(&self) -> Ledger {
        Ledger::new()
            .with(quantity::ENERGY, self.energy)
            .with(quantity::MOMENTUM, self.momentum)
    }
}

fn leaky(energy_leak: f64, momentum_leak: f64) -> Leaky {
    Leaky {
        energy: 1000.0,
        momentum: 1000.0,
        energy_leak,
        momentum_leak,
    }
}

/// **The failure a single tolerance produces, in both directions.**
///
/// One scheme loses momentum at `1e-7` a step — a Barnes-Hut tree, where each body sees its own
/// approximation of the rest and their mutual forces no longer cancel. Another leaks energy at
/// `1e-8`, which is a bug.
///
/// A single number cannot separate them:
///
/// - at `1e-9`, the correct momentum drift is refused;
/// - at `1e-6`, the energy leak passes.
///
/// Both are demonstrated here rather than asserted, because "one number is not enough" is the
/// kind of claim that is easy to state and easy to state wrongly.
#[test]
fn a_single_tolerance_cannot_separate_two_schemes() {
    let strict = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .with(leaky(0.0, 1e-7));
    let mut strict = strict;
    let refused = strict
        .advance(Time::from_si(1.0))
        .expect_err("1e-9 refuses a momentum drift the scheme cannot avoid");
    assert_eq!(refused.quantity, "momentum");

    let mut loose = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-6)
        .with(leaky(1e-8, 1e-7));
    loose
        .advance(Time::from_si(1.0))
        .expect("and 1e-6 lets a real energy leak straight through");

    // The leak was real and the loose run simply could not see it.
    let after = loose.ledger().get(quantity::ENERGY).unwrap_or(0.0);
    assert!(
        (1000.0 - after) / 1000.0 > 5e-9,
        "the energy really did leak: {after}"
    );
}

/// **A tolerance per quantity separates them.**
///
/// The same two schemes, with momentum given the accuracy its scheme can reach and energy kept at
/// the accuracy its scheme can reach. The correct momentum drift passes and the energy leak is
/// caught — which is the whole point, and neither is possible under one number.
#[test]
fn a_tolerance_per_quantity_catches_the_leak_and_allows_the_drift() {
    let mut sim = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .conservation_tolerance_for(quantity::MOMENTUM, 1e-6)
        .with(leaky(0.0, 1e-7));
    sim.advance(Time::from_si(1.0))
        .expect("momentum at 1e-6 allows what the tree cannot avoid");

    let mut sim = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .conservation_tolerance_for(quantity::MOMENTUM, 1e-6)
        .with(leaky(1e-8, 1e-7));
    let caught = sim
        .advance(Time::from_si(1.0))
        .expect_err("and energy at 1e-9 still sees the leak");
    assert_eq!(caught.quantity, "energy");

    // **The violation carries the tolerance that actually applied**, not the default. A reader
    // being told "1e-9" when the number in force was `1e-6` is worse than being told nothing.
    assert!(
        (caught.tolerance - 1e-9).abs() < 1e-30,
        "the violation should name the energy tolerance, said {}",
        caught.tolerance
    );
}

/// **The order the two setters are called in does not matter.**
///
/// `conservation_tolerance` sets the default and `conservation_tolerance_for` sets an override,
/// and a builder where the second silently undid the first would be a trap — the natural reading
/// of `.conservation_tolerance(x)` after an override is "change the default", not "reset
/// everything".
#[test]
fn the_default_and_the_overrides_are_independent() {
    let a = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-12)
        .conservation_tolerance_for(quantity::MOMENTUM, 1e-6);
    let b = Simulation::new(Schedule::Staggered)
        .conservation_tolerance_for(quantity::MOMENTUM, 1e-6)
        .conservation_tolerance(1e-12);

    for sim in [&a, &b] {
        assert_eq!(sim.tolerances().for_quantity(quantity::ENERGY), 1e-12);
        assert_eq!(sim.tolerances().for_quantity(quantity::MOMENTUM), 1e-6);
        assert_eq!(sim.tolerances().default_tolerance(), 1e-12);
    }
    assert_eq!(a.tolerances(), b.tolerances());

    // And what was overridden is readable, so a report can say what a run was checked against.
    let named: Vec<(&str, f64)> = a.tolerances().overrides().collect();
    assert_eq!(named, vec![(quantity::MOMENTUM, 1e-6)]);
}

/// **`audit` and `audit_with` agree when the tolerances are uniform.**
///
/// The compatibility claim, checked rather than assumed: the old signature is the new one with a
/// `Tolerances::uniform`, so every existing caller is unaffected.
#[test]
fn the_uniform_case_is_unchanged() {
    let before = Ledger::new()
        .with(quantity::ENERGY, 100.0)
        .with(quantity::MOMENTUM, 50.0);
    let after = Ledger::new()
        .with(quantity::ENERGY, 100.0 * (1.0 - 1e-7))
        .with(quantity::MOMENTUM, 50.0);

    for tol in [1e-9, 1e-6, 1e-3] {
        let old = audit("site", &before, &after, tol);
        let new = audit_with("site", &before, &after, &Tolerances::uniform(tol));
        assert_eq!(old.is_err(), new.is_err(), "at {tol}");
        if let (Err(a), Err(b)) = (old, new) {
            assert_eq!(a.quantity, b.quantity);
            assert_eq!(a.tolerance, b.tolerance);
        }
    }

    // A quantity nobody named falls back to the default rather than to something permissive.
    let t = Tolerances::uniform(1e-9).with(quantity::MOMENTUM, 1.0);
    assert_eq!(t.for_quantity(quantity::CHARGE), 1e-9);
    assert_eq!(t.for_quantity("a channel invented by a domain"), 1e-9);
}

/// **The overrides iterate in name order, whatever order they were set in.**
///
/// Determinism, and the reason `Tolerances` holds a `BTreeMap`. A violation's message and a
/// report's contents must not depend on the order a builder was called in.
#[test]
fn the_overrides_are_ordered_not_insertion_ordered() {
    let one = Tolerances::uniform(1e-9)
        .with(quantity::PHOTONS, 3.0)
        .with(quantity::CHARGE, 1.0)
        .with(quantity::MOMENTUM, 2.0);
    let other = Tolerances::uniform(1e-9)
        .with(quantity::MOMENTUM, 2.0)
        .with(quantity::PHOTONS, 3.0)
        .with(quantity::CHARGE, 1.0);

    let names: Vec<&str> = one.overrides().map(|(q, _)| q).collect();
    assert_eq!(names, vec!["charge", "momentum", "photons"]);
    assert_eq!(
        one.overrides().collect::<Vec<_>>(),
        other.overrides().collect::<Vec<_>>()
    );
    assert_eq!(one, other);
}
