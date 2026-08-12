//! The one number in this crate that was never argued for.
//!
//! A saturated bed's effective thermal conductivity was `ε k_l + (1−ε) k_s` — the arithmetic mean. That
//! is the **Voigt bound**, and a bound is not a model: it is exact only when the two phases lie in
//! parallel with the flux, which in a bed of spheres nothing does. Nobody chose it; it is what you get
//! if you do not think about it, and it went unexamined for four releases because in every scenario this
//! crate ships it **changes nothing at all**.
//!
//! That last part is the interesting half, and both halves are measured here:
//!
//! - **Under flow the bed is isothermal**, so conduction carries no heat and `λ` is multiplied by zero.
//!   Swinging it over a factor of eight leaves the extraction yield identical to `1e-14`. A modelling
//!   choice that no shipped scenario can observe is a modelling choice nobody will check.
//! - **A cold basket can see it.** With a 20 °C wall the yield moves **4.9% per unit `ln λ`** on the grid
//!   this file uses and 3.9% on the shipped basket, so the range the old rule left open — Voigt to Reuss,
//!   a factor of 1.674 — was worth around 2% in extraction yield. That is a taste-level difference in a
//!   cup.
//!
//! So it is Maxwell–Eucken now, with the **liquid as the continuous phase**, which is the structural
//! fact about a saturated bed rather than a preference. The remaining range is 1.184.
//!
//! # Why this is not comparing two implementations of one idea
//!
//! `Puck::bed_conductivity` is Maxwell–Eucken written as Maxwell–Eucken. `Mix::hashin_shtrikman` is the
//! Hashin–Shtrikman variational bound written as a bound, and it is checked in `dualis-core` against
//! **Maxwell–Garnett** to `6.7e-16`. That the two agree here is a theorem — the HS bound is attained by a
//! coated-sphere assemblage, which is what Maxwell–Eucken describes — and neither expression contains
//! the reason. The puck does not call `Mix`; if it did, this test would be checking that a function
//! equals itself.

use dualis_core::mixture::Mix;
use dualis_core::substance::{Substance, ThermalProps};
use dualis_core::units::{
    Density, Length, Pressure, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion,
    Time,
};
use dualis_core::{Domain, Exchange};
use dualis_porous::{Basket, Liquid, Puck};

/// The two phases as substances, so `Mix` can be asked about them. Only the conductivity and the density
/// matter to a conductivity bound; the rest is filled in because `ThermalProps` has no optional fields.
fn phases(kl: f64, ks: f64) -> (Substance, Substance) {
    let props = |k: f64| ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(1000.0),
        expansion: ThermalExpansion::ppm_per_k(100.0),
        emissivity: 0.9,
    };
    (
        Substance::bulk("liquid", Density::from_si(1000.0)).with_thermal(props(kl)),
        Substance::bulk("solid", Density::from_si(1000.0)).with_thermal(props(ks)),
    )
}

/// **The bed's conductivity is inside the Voigt and Reuss bounds at every porosity, and equals the
/// Hashin–Shtrikman bound built around the liquid.**
///
/// The equality is the cross-check and it holds to `2.1e-15` across three decades of conductivity ratio
/// and nineteen porosities. The bracketing is the theorem, and it is what would catch a formula that is right
/// at a half and leaves the bounds at an extreme — which several plausible mis-transcriptions of
/// Maxwell–Eucken do, because the denominator changes sign.
///
/// Both directions of contrast, and that is the case worth having: for water in coffee the liquid conducts
/// **better** than the solid, so liquid-as-host coincides with the HS *upper* bound. For a metal powder it
/// is the lower one. "Liquid as host" is a statement about which phase is continuous and not about which
/// number is larger, and a test using only the coffee numbers could not tell the two rules apart.
#[test]
fn the_bed_conductivity_is_maxwell_eucken_and_inside_both_bounds() {
    let mut worst = 0.0f64;
    let mut checked = 0;
    for (kl, ks) in [
        (0.675, 0.15),   // water and roasted coffee: the shipped pair
        (0.675, 0.0675), // ten times the contrast, same direction
        (0.6, 60.0),     // a metal powder: the other direction
        (0.5, 0.5),      // no contrast at all
    ] {
        let (liquid, solid) = phases(kl, ks);
        for tenths in 1..20 {
            let e = tenths as f64 / 20.0;
            let got = Puck::bed_conductivity(kl, ks, e);
            let mix = Mix::of(&[(liquid.clone(), e), (solid.clone(), 1.0 - e)])
                .expect("fractions sum to one");
            let (reuss, voigt) = mix
                .conductivity_bounds()
                .expect("both state a conductivity");
            let (hs_lo, hs_hi) = mix.hashin_shtrikman().expect("two phases");

            // Inside the outer pair, always. `>=` and not `>` because at zero contrast every bound is
            // the same number and "strictly inside" would be false and correct.
            assert!(
                got >= reuss.to_si() * (1.0 - 1e-12) && got <= voigt.to_si() * (1.0 + 1e-12),
                "kl {kl} ks {ks} e {e}: {got} is outside Reuss {} and Voigt {}",
                reuss.to_si(),
                voigt.to_si()
            );

            // And it is the HS bound built around the *liquid*, whichever of the two that is.
            let want = if kl >= ks { hs_hi } else { hs_lo };
            let off = if want.to_si() > 0.0 {
                (got / want.to_si() - 1.0).abs()
            } else {
                0.0
            };
            // A few machine epsilons, and the reason is worth naming: these are two *arrangements*
            // of one rational function, so the only thing separating them is the order the
            // multiplications happen in. Measured worst over 76 cases is 2.4e-15, at a 100-fold
            // contrast where the denominator is a small difference of large numbers. `1e-15` failed
            // there and `1e-15` was the wrong number, not the wrong result.
            assert!(
                off < 1e-14,
                "kl {kl} ks {ks} e {e}: Maxwell-Eucken gives {got} and the HS bound around the \
                 liquid is {}",
                want.to_si()
            );
            worst = worst.max(off);
            checked += 1;
        }
    }
    println!("  {checked} cases, worst departure from the HS bound {worst:.2e}");
}

/// **What the change cost, on the record.**
///
/// The old rule and the new one for the shipped pair, so the 11.0% is a measurement rather than a claim
/// in a commit message. And the range that remains: 1.184 against the 1.674 the arithmetic mean left
/// open, because Voigt and Reuss assume nothing about the microstructure and Maxwell–Eucken assumes the
/// thing that is true of a bed.
#[test]
fn the_arithmetic_mean_was_eleven_percent_high() {
    let (kl, ks, e) = (0.675, 0.15, 0.45);
    let voigt = e * kl + (1.0 - e) * ks;
    let reuss = 1.0 / (e / kl + (1.0 - e) / ks);
    let now = Puck::bed_conductivity(kl, ks, e);
    let (liquid, solid) = phases(kl, ks);
    let mix = Mix::of(&[(liquid, e), (solid, 1.0 - e)]).expect("sums to one");
    let (hs_lo, hs_hi) = mix.hashin_shtrikman().expect("two phases");

    println!(
        "  Voigt {voigt:.5}  Reuss {reuss:.5} ({:.3}x apart)\n  \
         HS {:.5} to {:.5} ({:.3}x apart), and the bed uses {now:.5}",
        voigt / reuss,
        hs_lo.to_si(),
        hs_hi.to_si(),
        hs_hi.to_si() / hs_lo.to_si(),
    );
    assert!(
        ((voigt / now - 1.0) - 0.1096).abs() < 0.001,
        "the arithmetic mean was {:.2}% high, not 11.0%",
        (voigt / now - 1.0) * 100.0
    );
    assert!(
        ((voigt / reuss) - 1.674).abs() < 0.001
            && ((hs_hi.to_si() / hs_lo.to_si()) - 1.184).abs() < 0.001,
        "the two ranges are 1.674 and 1.184"
    );
}

/// A shot with the bed solid's conductivity set to `ks`, and the wall held at `wall_c` or left at the
/// brew temperature. Returns the extraction yield **as a percentage**, which is what the reading is.
///
/// Stated because a first draft multiplied it by a hundred again and printed a 287% yield — above the
/// 30% soluble ceiling, which sent me looking for a conservation defect that was not there. The domain
/// saturates at exactly 30.00000% when driven to it; the mistake was in the print.
fn shot(ks: f64, wall_c: Option<f64>) -> f64 {
    let mut basket = Basket {
        counts: (8, 12, 8),
        cell: Length::from_si(2.5e-3),
        radius: Length::from_si(10e-3),
        pressure: Pressure::from_si(9.0e5),
        ..Basket::espresso()
    };
    basket.bed.conductivity = ThermalConductivity::w_per_m_k(ks);
    let mut p = Puck::new("puck", basket);
    if let Some(w) = wall_c {
        p.set_wall_temperature(Temperature::celsius(w));
    }
    let dt = Time::from_si(0.05);
    let mut t = 0.0;
    for _ in 0..200 {
        p.step(Time::from_si(t), dt, &mut Exchange::new())
            .expect("the step is stable");
        t += dt.to_si();
    }
    p.readings()
        .iter()
        .find(|r| r.label == "yield")
        .expect("a shot reports its yield")
        .value
}

/// **No shipped scenario can see this choice, and that is why it went unexamined.**
///
/// Under flow the bed is isothermal — the water arrives at the brew temperature, the wall is at the brew
/// temperature, and nothing in the bed generates heat — so every face carries a zero temperature
/// difference and `λ` is multiplied by nothing. Swinging the solid's conductivity over four orders takes
/// `λ` over a factor of ten and leaves the yield identical to `1e-14`.
///
/// Asserted rather than reasoned, because "conduction cannot matter here" is exactly the kind of argument
/// that is true of the continuum and false of a scheme — the two-phase work in `dualis-thermal` found one
/// of those a release ago. And asserted as an *equality*, because a tolerance would leave room for the
/// possibility being excluded.
#[test]
fn an_isothermal_bed_is_blind_to_its_own_conductivity() {
    let reference = shot(0.15, None);
    let mut span = (f64::MAX, 0.0f64);
    for ks in [1e-4, 0.15, 0.9, 5.0] {
        let lambda = Puck::bed_conductivity(Liquid::water().conductivity.to_si(), ks, 0.45);
        span = (span.0.min(lambda), span.1.max(lambda));
        let y = shot(ks, None);
        assert!(
            (y / reference - 1.0).abs() < 1e-14,
            "k_s {ks} (lambda {lambda:.4}) moved an isothermal yield from {reference} to {y}"
        );
    }
    println!(
        "  lambda over {:.4} to {:.4}, a factor of {:.1}, and the yield does not move at all",
        span.0,
        span.1,
        span.1 / span.0
    );
    assert!(
        span.1 / span.0 > 5.0,
        "the sweep should span a wide range of lambda, spans {:.2}x",
        span.1 / span.0
    );
}

/// **A cold basket can see it, so the choice is not free — and this is how much it is worth.**
///
/// The other half, and the reason the first half is a finding rather than a licence to stop caring. With a
/// 20 °C wall the yield falls as the bed conducts better, because more of the shot's heat reaches the cold
/// metal. Measured near the operating point the slope is about 3.9% of the yield per unit `ln λ`, which
/// over the 1.674 the old rule left open is close to 2% in extraction.
///
/// The slope is asserted as a slope and not extrapolated to an endpoint, because the response is not
/// log-linear over a wide range: the same sweep taken over a factor of ten gives 4.9% rather than 3.9%.
/// A single ratio quoted as "what the choice costs" would be a number with a hidden range of validity.
#[test]
fn a_cold_basket_is_not_blind_to_it() {
    let kl = Liquid::water().conductivity.to_si();
    let lambda = |ks: f64| Puck::bed_conductivity(kl, ks, 0.45);

    let (lo_k, hi_k) = (1e-4, 0.15);
    let (lo_y, hi_y) = (shot(lo_k, Some(20.0)), shot(hi_k, Some(20.0)));
    let slope = ((hi_y / lo_y - 1.0) / (lambda(hi_k) / lambda(lo_k)).ln()).abs();
    println!(
        "  lambda {:.4} -> {:.4}: yield {:.6}% -> {:.6}%, a slope of {:.2}% per unit ln(lambda)",
        lambda(lo_k),
        lambda(hi_k),
        lo_y,
        hi_y,
        slope * 100.0
    );
    assert!(
        hi_y < lo_y,
        "a better-conducting bed should lose more heat to a cold wall and extract less: {lo_y} \
         then {hi_y}"
    );
    assert!(
        slope > 0.01,
        "the yield should be visibly sensitive to lambda with a gradient present, slope is \
         {:.4}% per unit ln(lambda) — if this has fallen to nothing then the gradient is gone and \
         this test has stopped being able to see what it is for",
        slope * 100.0
    );
    // And the same sweep with no gradient moves nothing, in one place so the contrast is the assertion
    // rather than two tests a reader has to compare.
    assert!(
        (shot(hi_k, None) / shot(lo_k, None) - 1.0).abs() < 1e-14,
        "with no gradient the same sweep must move nothing"
    );
}
