//! What crosses a gap of nothing, checked against the closed form that says how much.
//!
//! Void broke conduction, which is right, and left two parts exchanging *nothing* — which is
//! wrong for any real clearance, because two surfaces that see each other radiate whether or not
//! there is anything between them. A vacuum gap has no other path at all, so radiation is not a
//! correction there: it is the whole answer.
//!
//! The form is the parallel-plate series, `σA(T₁⁴ − T₂⁴)/(1/ε₁ + 1/ε₂ − 1)`, exact for two
//! facing surfaces that see only each other. Everything below is against that or against a
//! conservation identity, and the one assumption it makes — a view factor of one — is checked
//! for the case it is *not* right in as well, because knowing where a closed form stops being
//! valid is worth as much as the closed form.

use pantometry_core::conserved::quantity;
use pantometry_core::units::{Length, Temperature, Time};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::Solid3D;

const SIGMA: f64 = 5.670_374_419e-8;

/// Two single-cell plates with `gap` empty cells between them, hot on the left.
fn plates(gap: usize, emissivity: f64, hot_c: f64, cold_c: f64) -> Solid3D {
    let mut metal = Substance::copper();
    if let Some(thermal) = metal.thermal.as_mut() {
        thermal.emissivity = emissivity;
    }
    let n = 2 + gap;
    let mut block = Solid3D::new(
        "plates",
        metal,
        (1, 1, n),
        Length::mm(5.0),
        Temperature::celsius(cold_c),
    );
    block.set_temperature(0, 0, 0, Temperature::celsius(hot_c));
    block.empty(|_, _, k| k > 0 && k + 1 < n)
}

/// **One step's exchange is the parallel-plate closed form, to rounding.**
///
/// The claim the whole feature rests on. A single step from a known pair of temperatures moves a
/// known number of joules, and the number is `σA(T₁⁴ − T₂⁴)/(1/ε₁ + 1/ε₂ − 1)·dt` — computed here
/// from the constants rather than from anything the domain did.
///
/// The tolerance is the enthalpy round trip, not a physical budget: the step is explicit and
/// reads the temperatures it started with, so the joules that move are exactly that product.
#[test]
fn one_step_moves_what_the_parallel_plate_form_says() {
    let emissivity = 0.8;
    let (hot_c, cold_c) = (500.0, 20.0);
    let mut block = plates(1, emissivity, hot_c, cold_c);

    let before_hot = block.temperature_at(0, 0, 0).to_si();
    let before_cold = block.temperature_at(0, 0, 2).to_si();
    let capacity = block.heat_capacity().to_si() / 2.0; // two solid cells, equal

    let dt = 1e-4;
    let mut bus = Exchange::new();
    block
        .step(Time::ZERO, Time::from_si(dt), &mut bus)
        .expect("a stable step");

    let moved = (before_cold - block.temperature_at(0, 0, 2).to_si()).abs() * capacity;
    let area = 5.0e-3 * 5.0e-3;
    let expect =
        SIGMA * area * (before_hot.powi(4) - before_cold.powi(4)) / (2.0 / emissivity - 1.0) * dt;

    assert!(
        (moved - expect).abs() / expect < 1e-9,
        "the gap should carry {expect:e} J and carried {moved:e} J"
    );
}

/// **What leaves one plate arrives in the other, to the bit.**
///
/// The exchange is written as one antisymmetric statement, so a gap conserves the way a
/// conduction face does rather than to a tolerance. Checked over a long run, where a per-step
/// error would accumulate into something a single step could hide.
#[test]
fn a_gap_conserves_exactly() {
    let mut block = plates(1, 0.9, 600.0, 20.0);
    let opening = block.ledger().get(quantity::ENERGY).unwrap_or(0.0);

    let dt = block.max_stable_dt(Time::ZERO);
    let mut bus = Exchange::new();
    for _ in 0..5000 {
        block.step(Time::ZERO, dt, &mut bus).expect("a stable step");
    }
    let closing = block.ledger().get(quantity::ENERGY).unwrap_or(0.0);

    // The plates really did exchange, or this proves nothing.
    assert!(
        block.temperature_at(0, 0, 2).to_si() > Temperature::celsius(30.0).to_si(),
        "the cold plate should have warmed: {:.2} K",
        block.temperature_at(0, 0, 2).to_si()
    );
    let scale = opening.abs().max(1.0);
    assert!(
        (closing - opening).abs() / scale < 1e-12,
        "an isolated pair conserves: {opening:e} to {closing:e}"
    );
}

/// **Two plates alone come to the same temperature and stop.**
///
/// Radiation between two bodies that see only each other has one fixed point and it is equality:
/// `T₁⁴ − T₂⁴` is zero only there. The mean is the conserved quantity — equal capacities, so the
/// answer is the average of the two starts — and neither overshoots it.
#[test]
fn an_isolated_pair_settles_at_the_mean() {
    let (hot_c, cold_c) = (400.0, 0.0);
    let mut block = plates(1, 0.9, hot_c, cold_c);
    let dt = block.max_stable_dt(Time::ZERO);
    let mut bus = Exchange::new();
    for _ in 0..400_000 {
        block.step(Time::ZERO, dt, &mut bus).expect("a stable step");
    }

    let hot = block.temperature_at(0, 0, 0).to_si();
    let cold = block.temperature_at(0, 0, 2).to_si();
    // Copper's specific heat is constant here, so equal cells means the mean in kelvin.
    let mean = 0.5 * (Temperature::celsius(hot_c).to_si() + Temperature::celsius(cold_c).to_si());
    assert!(
        (hot - cold).abs() < 1.0,
        "they should have levelled: {hot:.2} K against {cold:.2} K"
    );
    assert!(
        (0.5 * (hot + cold) - mean).abs() < 1e-6,
        "and the mean is conserved: {:.4} K against {mean:.4} K",
        0.5 * (hot + cold)
    );
}

/// **A wider gap carries the same, and that is the boundary condition, not a bug.**
///
/// This assertion used to be labelled as a known defect: the parallel-plate form has a view factor
/// of one, two facing plates only see all of each other when the gap is narrow, so charging every
/// width as one looked like an approximation waiting to be fixed. Measuring it says otherwise.
///
/// The gap here is bounded sideways by the **block's own outer faces**, and those are insulated —
/// implemented, literally, as a mirror in `Solid3D::mirrored`. A mirror puts an image of each
/// plate beyond it, and the images tile the plane, so the pair *is* two infinite parallel plates
/// and `F̄ = 1` is exact at every width. The old comment described a different geometry from the
/// one the model has.
///
/// What is worth knowing is how different the answer would be for somebody who meant the gap to be
/// open to space — two small plates floating in vacuum, where most of what leaves one does miss
/// the other. That is the ordinary view factor, and `GapPatch::view_factor` computes it: for these
/// 5 mm cells it is **0.1998** at one cell of clearance and **0.0124** at five, a factor of
/// sixteen between two cases this model calls identical. The number is reported rather than
/// charged, because which reading is right is the user's statement about their boundary and not
/// something the grid can infer.
#[test]
fn a_wide_gap_carries_the_same_because_the_sides_are_mirrors() {
    let step_once = |gap: usize| {
        let mut block = plates(gap, 0.8, 500.0, 20.0);
        let cold_index = 1 + gap;
        let before = block.temperature_at(0, 0, cold_index).to_si();
        let mut bus = Exchange::new();
        block
            .step(Time::ZERO, Time::from_si(1e-4), &mut bus)
            .expect("a stable step");
        block.temperature_at(0, 0, cold_index).to_si() - before
    };

    let narrow = step_once(1);
    let wide = step_once(5);
    assert!(
        (narrow - wide).abs() / narrow < 1e-12,
        "mirrored sides make every width the infinite-plate answer: {narrow:e} against {wide:e}"
    );

    // And the open-gap factors, which say how much that reading is worth. Both computed from the
    // closed form rather than read back, and both far from one — which is exactly why the reading
    // has to be stated instead of assumed.
    let factor_of = |gap: usize| {
        pantometry_thermal::GapPatch {
            pairs: 1,
            span: (5e-3, 5e-3),
            distance: gap as f64 * 5e-3,
            rectangular: true,
        }
        .view_factor()
    };
    assert!(
        (factor_of(1) - 0.199_82).abs() < 5e-6 && (factor_of(5) - 0.012_4).abs() < 5e-5,
        "open to space these would differ by sixteen: {:.5} and {:.5}",
        factor_of(1),
        factor_of(5)
    );
}

/// **The limit knows about the gap**, so a pair hot enough to radiate hard is stepped short
/// enough to survive it.
///
/// The radiative conductance is `4σAT³/(1/ε₁+1/ε₂−1)` and grows with the cube of temperature, so
/// a pair at 1500 K is far stiffer than the same pair at 300 K. A limit that ignored gaps would
/// hand the hot pair the cold pair's step, and an explicit exchange past its limit overshoots —
/// the hot plate ends below the cold one, every step, while the pair conserves perfectly.
#[test]
fn a_hot_pair_gets_a_shorter_limit_and_does_not_overshoot() {
    let cool = plates(1, 0.9, 40.0, 20.0).max_stable_dt(Time::ZERO).to_si();
    let hot = plates(1, 0.9, 1200.0, 20.0)
        .max_stable_dt(Time::ZERO)
        .to_si();
    assert!(
        hot < cool,
        "a hotter pair radiates harder and must take shorter steps: {hot:e} against {cool:e}"
    );

    // And at the reported limit the exchange does not cross over.
    let mut block = plates(1, 0.9, 1200.0, 20.0);
    let dt = block.max_stable_dt(Time::ZERO);
    let mut bus = Exchange::new();
    for _ in 0..200 {
        block.step(Time::ZERO, dt, &mut bus).expect("a stable step");
        assert!(
            block.temperature_at(0, 0, 0).to_si() >= block.temperature_at(0, 0, 2).to_si(),
            "the hot plate must not end up colder than the one it is heating"
        );
    }
}
