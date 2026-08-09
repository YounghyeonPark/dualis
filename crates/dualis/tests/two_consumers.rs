//! Two domains taking from one channel, which every scene and test here had avoided.
//!
//! `Exchange::take` empties a channel, so the second consumer gets zero — and the conservation
//! audit reports it clean, because everything published *was* consumed. Two plates under one
//! lamp warm at the rate of one plate and every total agrees to the bit.
//!
//! Found by trying to couple more than two domains. Pairwise coupling cannot reach it, and
//! pairwise is all this workspace had ever built: thirteen scenes, maximum two domains each.

use dualis::prelude::*;

fn plate(name: &str) -> LumpedMass {
    LumpedMass::new(
        name,
        Substance::aluminium_6061(),
        Volume::cm3(10.0),
        Length::mm(1.5),
        Temperature::celsius(20.0),
        // A real surface, so the plate has a finite time constant and the schedule has a limit
        // to respect. Zero area gives an infinite one, which multirate turns into a hang.
        Environment::still_air(Temperature::celsius(20.0), Area::cm2(60.0)),
    )
}

fn coil() -> Winding {
    Winding::of_resistance(
        "coil",
        Resistance::ohm(2.0),
        0.0,
        Temperature::celsius(20.0),
    )
    .driven_at(Current::a(1.0))
    .with_reserve(1_000.0)
}

/// **A second consumer is refused, and named.**
///
/// The alternative is what the library did until now: silently give the first plate everything
/// and the second nothing, with the audit green at any tolerance.
#[test]
fn a_second_consumer_of_one_channel_is_refused() {
    let mut sim = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .with(coil())
        .with(plate("first"))
        .with(plate("second"));

    let violation = sim
        .advance(Time::s(1.0))
        .expect_err("two plates cannot share a channel that empties on the first take");
    assert_eq!(violation.quantity, "energy");
    assert!(
        violation.site.starts_with("second"),
        "the *second* taker is the one that found nothing: {}",
        violation.site
    );
    assert!(
        violation.site.contains("already emptied"),
        "{}",
        violation.site
    );
    assert_eq!(sim.time().to_si(), 0.0, "a refused step keeps the clock");
}

/// **One consumer subcycling many times is not that, and is not refused.**
///
/// The distinction the check exists to draw. Under `Multirate` a plate steps many times per
/// interval and takes on each, which is one consumer collecting its own interval in pieces —
/// counted per *turn*, not per call. Getting this wrong would refuse every multirate scene in
/// the repository, which is how it was caught: four of them failed on the first attempt.
#[test]
fn one_consumer_taking_on_every_substep_is_fine() {
    let mut sim = Simulation::new(Schedule::Multirate)
        .conservation_tolerance(1e-9)
        .with(coil())
        .with(plate("only"));

    // 300 s against a time constant of 576 s, so max_stable_dt is 57.6 s and this needs six.
    let report = sim
        .advance(Time::s(300.0))
        .expect("one consumer, many substeps");
    let substeps: std::collections::BTreeMap<_, _> = report.substeps.into_iter().collect();
    assert!(
        substeps["only"] > 1,
        "the plate should have subcycled; it took {} steps",
        substeps["only"]
    );
    assert!(
        (sim.domain_as::<LumpedMass>("only")
            .unwrap()
            .absorbed_energy()
            .to_si()
            - 600.0)
            .abs()
            < 1e-9,
        "every joule should still arrive"
    );
}
