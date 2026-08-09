//! Electricity publishing heat and thermal consuming it, over a bus neither crate names.
//!
//! Lives in the facade rather than in either domain, because a domain that dev-depended on
//! another would break the property the crate split exists to hold — and `invariant-guard`
//! greps the manifest, so a dev-dependency is caught the same as a real one. The README says
//! this is where the cross-domain integration tests live; this is the sixth pair.

use dualis::prelude::*;
use dualis_electrical::Winding;

/// **A winding warms a plate, and the joules that arrive are `I²R·t` to the bit.**
///
/// The coupling that motivated the crate: electrical publishes, thermal consumes, over the
/// channel that already existed, with neither crate naming the other. Run under a live
/// conservation audit at 1e-9.
///
/// The rise is checked against `Q/C` with the capacity computed here from the substance, so
/// this is not the simulation compared with itself.
#[test]
fn a_winding_warms_a_plate_by_the_joules_it_dissipated() {
    let volume = Volume::from_si(60e-3 * 60e-3 * 3e-3);
    let start = Temperature::celsius(20.0);
    let seconds = 5.0;

    let coil = Winding::of_resistance("coil", Resistance::ohm(2.0), 0.0, start)
        .driven_at(Current::a(1.5))
        .with_reserve(1_000.0);
    let watts = coil.dissipation().to_si();
    assert!(
        (watts - 4.5).abs() < 1e-12,
        "I^2 R = 1.5^2 * 2 = 4.5, got {watts}"
    );

    let plate = LumpedMass::new(
        "plate",
        Substance::aluminium_6061(),
        volume,
        Length::mm(1.5),
        start,
        // Ambient *at* the starting temperature and no radiation would still leave convection
        // shedding heat as it warms, so this instead gives it a vanishing area: the plate keeps
        // every joule, and the closed form below is exact rather than approximate.
        Environment::still_air(start, Area::from_si(0.0)),
    );

    let mut sim = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .with(coil)
        .with(plate);
    let dt = 0.01;
    for _ in 0..(seconds / dt) as usize {
        sim.advance(Time::s(dt)).expect("the books close");
    }

    let capacity = Substance::aluminium_6061()
        .heat_capacity(volume)
        .expect("aluminium has a specific heat")
        .to_si();
    let want = watts * seconds / capacity;
    let got = sim
        .domain_as::<LumpedMass>("plate")
        .expect("the plate is there")
        .temperature()
        .to_si()
        - start.to_si();
    assert!(
        (got / want - 1.0).abs() < 1e-9,
        "rise {got:.9} K against {want:.9} K"
    );

    // The winding's books say the same thing from the other side.
    let coil = sim.domain_as::<Winding>("coil").expect("the coil is there");
    assert!((coil.dissipated_energy().to_si() - watts * seconds).abs() < 1e-9);
    assert!((coil.reserve().to_si() - (1_000.0 - watts * seconds)).abs() < 1e-9);
}

/// **A source with an unlimited supply is refused by the domain, because the audit cannot.**
///
/// This test was written expecting the *audit* to catch it and it did not: the ledger reported
/// `inf` before and `inf` after, `inf` equals itself, and a winding pouring joules into a plate
/// ran green. An infinite reserve is not a large number; it is a number that cannot be
/// subtracted, so it switches the check off rather than failing it.
///
/// The same shape as `dualis-world`'s lamp, where a skipped `with_reserve` left a scene
/// auditing clean at tolerance zero with the lamp doing nothing. So `Winding::step` refuses
/// while the reserve is infinite, and this pins that it is the domain saying so and by name.
#[test]
fn a_winding_with_no_reserve_is_caught_by_the_audit() {
    let coil = Winding::of_resistance(
        "coil",
        Resistance::ohm(2.0),
        0.0,
        Temperature::celsius(20.0),
    )
    .driven_at(Current::a(1.5));
    let plate = LumpedMass::new(
        "plate",
        Substance::aluminium_6061(),
        Volume::from_si(1e-4),
        Length::mm(1.5),
        Temperature::celsius(20.0),
        Environment::still_air(Temperature::celsius(20.0), Area::from_si(0.0)),
    );
    let mut sim = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .with(coil)
        .with(plate);

    let violation = sim
        .advance(Time::s(1.0))
        .expect_err("joules from an infinite tank are joules from nowhere");
    assert_eq!(violation.site, "coil", "the winding names itself");
    assert!(
        violation.quantity.contains("with_reserve"),
        "{}",
        violation.quantity
    );
    assert_eq!(sim.time().to_si(), 0.0, "a refused step keeps the clock");
}
