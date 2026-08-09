//! `Winding` against closed forms, and against the class its own audit cannot see.
//!
//! The conservation audit checks that the joules published are the joules consumed. It cannot
//! check that `I²R` is the *right* number of joules: both sides of the bus agree perfectly about
//! whatever is published, so a resistance wrong by a factor of two balances the books exactly and
//! reports a winding that runs at half the temperature it will. Every check here is therefore
//! against arithmetic done in this file, or against a limit the physics has to satisfy.

use dualis_core::{Domain, Exchange};
use dualis_electrical::{Winding, COPPER_ALPHA, COPPER_RESISTIVITY_20C, HEAT};
use dualis_units::{Conductance, Current, Length, Resistance, Temperature, Time, Voltage};

/// **`R = ρL/A`, against the resistivity written out here.**
///
/// The one number this crate cannot derive — copper's resistivity — is a constant, so the check
/// is that the geometry is used correctly around it: double the length and the resistance
/// doubles, double the area and it halves. A transposed `L` and `A` passes any single-point
/// check and fails both of these.
#[test]
fn resistance_is_resistivity_times_length_over_area() {
    let (length, area) = (Length::m(24.0), 0.35e-6);
    let coil = Winding::of_copper("coil", length, area, Temperature::celsius(20.0));

    let want = COPPER_RESISTIVITY_20C * 24.0 / 0.35e-6;
    assert!(
        (coil.resistance().to_si() / want - 1.0).abs() < 1e-12,
        "{} against {want}",
        coil.resistance().to_si()
    );

    let twice_as_long = Winding::of_copper("b", Length::m(48.0), area, Temperature::celsius(20.0));
    let twice_as_fat = Winding::of_copper("c", length, 2.0 * area, Temperature::celsius(20.0));
    assert!((twice_as_long.resistance().to_si() / want - 2.0).abs() < 1e-12);
    assert!((twice_as_fat.resistance().to_si() / want - 0.5).abs() < 1e-12);
}

/// **The temperature coefficient, at the two points where it is a known number.**
///
/// `ρ(T) = ρ₂₀(1 + α(T − 20 °C))`, so the resistance at 20 °C is exactly the reference value and
/// the *slope* is `α` per kelvin. Checked at 20 °C, where the answer is 1×, and at 120 °C, where
/// a hundred kelvin of copper is `1 + 0.393`, a factor a motor designer knows by heart.
#[test]
fn copper_gains_four_tenths_of_a_percent_per_kelvin() {
    let make = |c: f64| {
        Winding::of_copper("coil", Length::m(10.0), 1e-6, Temperature::celsius(c))
            .resistance()
            .to_si()
    };
    let r20 = make(20.0);

    assert_eq!(make(20.0), r20, "20 C is the reference point exactly");
    assert!((make(21.0) / r20 - (1.0 + COPPER_ALPHA)).abs() < 1e-12);

    // A hundred kelvin hotter: 39.3% more resistance. This is the number that makes a winding's
    // losses at temperature meaningfully different from its losses on the bench.
    assert!(
        (make(120.0) / r20 - 1.393).abs() < 1e-12,
        "at 120 C: {}",
        make(120.0) / r20
    );

    // Below the reference it falls, which is the sign of the effect and not just its size.
    assert!(make(-30.0) < r20);
}

/// **Constant current and constant voltage are the same winding at the same operating point.**
///
/// Drive it at `I`, read off the voltage it develops, then drive an identical winding *from*
/// that voltage. Every quantity must match to the bit: same current, same voltage, same
/// dissipation. It is the one identity that ties the two drive modes together, and a
/// `V²/R`-vs-`I²R` mixup — using one formula with the other's stored value — breaks it while
/// leaving each mode self-consistent.
#[test]
fn the_two_drive_modes_meet_at_the_same_operating_point() {
    let build = || Winding::of_copper("coil", Length::m(24.0), 0.35e-6, Temperature::celsius(75.0));
    let by_current = build().driven_at(Current::a(3.0));
    let by_voltage = build().driven_from(by_current.voltage());

    // Not bit-for-bit, and it cannot be. `V` here is `I·R`, so `V²/R` re-does a multiply and a
    // divide that `I²R` performs in a different order, and floating-point multiplication is not
    // associative — the two land one unit in the last place apart. The tolerance is a few
    // machine epsilons because that is the size of the effect, not because it was loosened
    // until this passed: at 1e-15 relative it is about five eps, and a genuine formula mixup
    // would be off by the resistance, a factor of 1.15 here.
    let (i, v) = (
        by_current.dissipation().to_si(),
        by_voltage.dissipation().to_si(),
    );
    assert!(
        (i / v - 1.0).abs() < 1e-15,
        "{i} against {v}, a relative difference of {:e}",
        (i / v - 1.0).abs()
    );
    assert!((by_voltage.current().to_si() - 3.0).abs() < 1e-12);

    // And `P = VI`, computed here from the two readings rather than from either formula.
    let p = by_current.voltage().to_si() * by_current.current().to_si();
    assert!((by_current.dissipation().to_si() / p - 1.0).abs() < 1e-12);
}

/// **The feedback has opposite signs in the two modes, which is the whole runaway question.**
///
/// A winding held at constant *current* dissipates `I²R` and gets worse as it warms. The same
/// winding on a constant *voltage* dissipates `V²/R` and gets better. That sign is why a
/// current-driven coil can run away and a voltage-driven one cannot, and it is a claim about the
/// physics rather than about the arithmetic — so it is checked by warming one of each and
/// looking at which way the power went.
#[test]
fn current_drive_rises_with_temperature_and_voltage_drive_falls() {
    let mut hot_i = Winding::of_copper("i", Length::m(10.0), 1e-6, Temperature::celsius(20.0))
        .driven_at(Current::a(2.0));
    let cold_watts = hot_i.dissipation().to_si();
    hot_i.at_temperature(Temperature::celsius(120.0));
    assert!(
        (hot_i.dissipation().to_si() / cold_watts - 1.393).abs() < 1e-12,
        "constant current should rise by exactly the resistance ratio"
    );

    let mut hot_v = Winding::of_copper("v", Length::m(10.0), 1e-6, Temperature::celsius(20.0))
        .driven_from(Voltage::v(5.0));
    let cold_watts = hot_v.dissipation().to_si();
    hot_v.at_temperature(Temperature::celsius(120.0));
    assert!(
        (hot_v.dissipation().to_si() / cold_watts - 1.0 / 1.393).abs() < 1e-12,
        "constant voltage should fall by exactly the resistance ratio"
    );
}

/// **Nothing is published when nothing is driven.**
///
/// A winding built and never given a current is at zero amps, so it dissipates nothing and
/// publishes nothing — and an amount published with no consumer is a `Violation`, so "nothing"
/// has to mean *nothing* rather than a very small number.
#[test]
fn an_undriven_winding_publishes_nothing_at_all() {
    let mut coil = Winding::of_copper("coil", Length::m(10.0), 1e-6, Temperature::celsius(20.0));
    assert_eq!(coil.dissipation().to_si(), 0.0);

    let mut bus = Exchange::new();
    coil.step(Time::s(0.0), Time::s(1.0), &mut bus).unwrap();
    assert_eq!(bus.peek(HEAT), 0.0);

    // A short across an ideal voltage source is not a physical answer, and an infinity here
    // would reach the audit as a NaN a step later with its origin lost.
    let shorted = Winding::of_resistance(
        "short",
        Resistance::ohm(0.0),
        0.0,
        Temperature::celsius(20.0),
    )
    .driven_from(Voltage::v(12.0));
    assert_eq!(shorted.dissipation().to_si(), 0.0);
    assert_eq!(shorted.current().to_si(), 0.0);
}

/// **The runaway threshold, and that it is a threshold rather than a number.**
///
/// `dP/dT > dQ_out/dT` is the condition, and for a constant-current winding both sides are
/// closed forms: `I²R₂₀α` against `g`. So `I_crit = √(g/(R₂₀α))` is exact, and the check is that
/// the *inequality flips* across it — a formula reproduced from itself would pass any single
/// point.
#[test]
fn the_runaway_current_is_where_the_feedback_overtakes_the_heat_path() {
    let coil = Winding::of_copper("coil", Length::m(62.0), 0.35e-6, Temperature::celsius(20.0))
        .driven_at(Current::a(1.0));
    let g = Conductance::w_per_k(0.203);
    let crit = coil.runaway_current(g).expect("a current drive has one");

    // Against the formula written out here.
    let r_20 = 1.724e-8 * 62.0 / 0.35e-6;
    let want = (0.203_f64 / (r_20 * 0.00393)).sqrt();
    assert!(
        (crit.to_si() / want - 1.0).abs() < 1e-12,
        "{} A against {want} A",
        crit.to_si()
    );

    // And that it separates. `dP/dT` is `I²R₂₀α`, computed here from two dissipations a kelvin
    // apart rather than from the same expression the threshold uses — so this compares the
    // slope the model actually has against the conductance, not a formula with itself.
    let slope_at = |amps: f64| {
        let w = Winding::of_copper("c", Length::m(62.0), 0.35e-6, Temperature::celsius(20.0))
            .driven_at(Current::a(amps));
        (w.dissipation_at(Temperature::celsius(21.0)).to_si()
            - w.dissipation_at(Temperature::celsius(20.0)).to_si())
            / 1.0
    };
    let below = crit.to_si() * 0.9;
    let above = crit.to_si() * 1.1;
    assert!(slope_at(below) < 0.203, "below: {} W/K", slope_at(below));
    assert!(slope_at(above) > 0.203, "above: {} W/K", slope_at(above));
    // At the threshold itself the two are equal, which is what makes it one.
    assert!(
        (slope_at(crit.to_si()) / 0.203 - 1.0).abs() < 1e-12,
        "at the threshold: {} W/K",
        slope_at(crit.to_si())
    );

    // **The joints move it.** Reaching air through 0.9 and 2.4 W/K of joints and then 0.294 W/K
    // of convection is a series conductance of 0.203, not 0.294 — and the threshold falls from
    // 4.95 A to 4.11 A. A lumped model reports 17% of margin that is not there.
    let surface_only = coil
        .runaway_current(Conductance::w_per_k(0.294))
        .expect("still a current drive");
    let series = 1.0 / (1.0 / 0.9 + 1.0 / 2.4 + 1.0 / 0.294);
    let with_joints = coil
        .runaway_current(Conductance::w_per_k(series))
        .expect("still a current drive");
    assert!(
        (surface_only.to_si() - 4.949).abs() < 1e-3,
        "{}",
        surface_only.to_si()
    );
    assert!(
        (with_joints.to_si() - 4.111).abs() < 1e-3,
        "{}",
        with_joints.to_si()
    );
    assert!(with_joints.to_si() < surface_only.to_si());

    // A voltage-driven winding cannot run away: `V²/R` falls as it warms, so there is no
    // threshold and reporting one would be worse than reporting none.
    let from_volts = Winding::of_copper("v", Length::m(62.0), 0.35e-6, Temperature::celsius(20.0))
        .driven_from(Voltage::v(12.0));
    assert!(from_volts.runaway_current(g).is_none());
    assert!(
        from_volts
            .dissipation_at(Temperature::celsius(120.0))
            .to_si()
            < from_volts
                .dissipation_at(Temperature::celsius(20.0))
                .to_si()
    );
}

/// `dissipation()` is `dissipation_at(its own temperature)`, to the bit.
///
/// The pure function and the method must not be two implementations that happen to agree —
/// which is the state `linearised_loss_conductance` was extracted to prevent in `dualis-thermal`.
#[test]
fn the_pure_function_and_the_method_are_the_same_arithmetic() {
    for c in [-40.0, 20.0, 75.0, 180.0] {
        let w = Winding::of_copper("coil", Length::m(24.0), 0.35e-6, Temperature::celsius(c))
            .driven_at(Current::a(3.0));
        assert_eq!(
            w.dissipation().to_si().to_bits(),
            w.dissipation_at(Temperature::celsius(c)).to_si().to_bits()
        );
        assert_eq!(
            w.resistance().to_si().to_bits(),
            w.resistance_at(Temperature::celsius(c)).to_si().to_bits()
        );
    }
}
