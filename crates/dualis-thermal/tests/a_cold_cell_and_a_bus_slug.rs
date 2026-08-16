//! The two ways a cooled face could be handed a step it could not take.
//!
//! Both were found by review, both were reproduced before they were fixed, and both had the
//! same signature: the block does not diverge loudly, it **overshoots the air it is losing to**
//! while the conservation audit stays perfectly happy — because the heat really did leave.

use dualis_core::conserved::quantity;
use dualis_core::units::{Area, Length, Temperature, Time};
use dualis_core::{Domain, Exchange, Substance};
use dualis_thermal::{Environment, Face, Solid3D, HEAT};

/// **A cell colder than its surroundings is not handed a step sized for the wrong curve.**
///
/// The radiative derivative `4εσT³` dominates the secant only for `T ≥ T∞`; below ambient it
/// reverses without bound. Measured before the fix: a 20 mm cell at 20 °C inside a 700 °C
/// radiant enclosure had a true step ratio of **13**, and one accepted step took it to 8847 °C.
#[test]
fn a_cell_colder_than_its_surroundings_is_not_overdriven() {
    let glass = Substance::from_name("borosilicate").expect("the catalogue has it");
    let mut block = Solid3D::new(
        "part",
        glass,
        (1, 1, 1),
        Length::mm(20.0),
        Temperature::celsius(20.0),
    )
    .losing_from(
        Face::ZMax,
        Environment {
            // A radiant enclosure far hotter than the part: the loss runs the other way.
            ambient: Temperature::celsius(700.0),
            convection_w_per_m2_k: 0.0,
            area: Area::from_si(20.0e-3 * 20.0e-3),
        },
    );

    let limit = block.max_stable_dt(Time::ZERO);
    let mut bus = Exchange::new();
    block
        .step(Time::ZERO, limit, &mut bus)
        .expect("the reported limit must be a step this domain can take");

    let after = block.temperature_at(0, 0, 0).to_si();
    let ambient = Temperature::celsius(700.0).to_si();
    assert!(
        after > Temperature::celsius(20.0).to_si() - 1e-9 && after <= ambient,
        "a cell warming toward its enclosure must not pass it: {after} K against {ambient} K"
    );
}

/// **Heat off the bus arrives after the limit was computed, and the step is re-checked against
/// the state the film actually meets.**
///
/// `stability_ratio` reads the temperatures the step *started* with. Bus heat is spread over
/// every cell before the film is applied, so a slug can make a cell far stiffer than the limit
/// was sized for. Measured before the fix, on a cube at ambient given its rise in one accepted
/// step: 300 K took a corner to −137 °C and 900 K to −7007 °C, with the ledger balancing
/// perfectly throughout.
///
/// The domain refuses now. Refusing rather than clamping, because a clamp would quietly return
/// a different answer near the limit and this workspace would rather stop than be plausible.
#[test]
fn a_bus_slug_cannot_shed_a_cell_past_ambient() {
    let pla = Substance::from_name("pla").expect("the catalogue has it");
    let side = 4.0 * 5.0e-3;
    let build = || {
        let mut block = Solid3D::new(
            "part",
            pla.clone(),
            (4, 4, 4),
            Length::mm(5.0),
            Temperature::celsius(20.0),
        );
        for face in Face::ALL {
            block = block.losing_from(
                face,
                Environment {
                    ambient: Temperature::celsius(20.0),
                    convection_w_per_m2_k: 0.0,
                    area: Area::from_si(side * side),
                },
            );
        }
        block
    };

    // A slug big enough that the film, met at the post-slug temperature, would overshoot.
    let mut block = build();
    let limit = block.max_stable_dt(Time::ZERO);
    let capacity = block.heat_capacity().to_si();
    let mut bus = Exchange::new();
    bus.publish(HEAT, 900.0 * capacity);

    match block.step(Time::ZERO, limit, &mut bus) {
        // Either the step is refused by name...
        Err(v) => assert_eq!(v.quantity, "surface loss ratio", "{v}"),
        // ...or it was taken, and then no cell may sit below the air it is losing to.
        Ok(()) => {
            let ambient = Temperature::celsius(20.0).to_si();
            for k in 0..4 {
                for j in 0..4 {
                    for i in 0..4 {
                        let t = block.temperature_at(i, j, k).to_si();
                        assert!(
                            t >= ambient - 1e-6,
                            "cell ({i},{j},{k}) was shed to {t} K, below the {ambient} K air"
                        );
                    }
                }
            }
        }
    }
}

/// **The books still balance while all of that is refused or survived.** The audit could never
/// see either defect — the heat genuinely left — so this is not the check that catches them; it
/// is the check that the fixes did not break the thing that was right.
#[test]
fn the_books_balance_through_a_refusal() {
    let glass = Substance::from_name("borosilicate").expect("the catalogue has it");
    let mut block = Solid3D::new(
        "part",
        glass,
        (2, 2, 2),
        Length::mm(10.0),
        Temperature::celsius(400.0),
    )
    .losing_from(
        Face::ZMax,
        Environment::still_air(Temperature::celsius(20.0), Area::from_si(4.0e-4)),
    );

    let opening = block.ledger().get(quantity::ENERGY).unwrap_or(0.0);
    let dt = block.max_stable_dt(Time::ZERO);
    let mut bus = Exchange::new();
    for _ in 0..200 {
        block.step(Time::ZERO, dt, &mut bus).expect("a stable step");
    }
    let closing = block.ledger().get(quantity::ENERGY).unwrap_or(0.0);
    let moved = block.lost_energy().to_si();
    assert!(moved > 0.0);
    assert!(
        (closing - opening).abs() / moved < 1e-12,
        "the ledger moved by {:e} J against {moved:e} J shed",
        closing - opening
    );
}
