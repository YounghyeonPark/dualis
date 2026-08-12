//! The algebra of a composite: what is exact, what is bounded, and what is refused.
//!
//! `a_composite.rs` in `dualis-thermal` checks the bounds against a resolved geometry that attains them.
//! This file checks the arithmetic that produces them, and it is where the *exact* claims live — because
//! an exact claim can be asserted as an equality, and an equality is a much stronger test than a
//! tolerance. Density, volumetric heat capacity and latent heat are conservation statements, so every
//! assertion about them here is `1e-15` or better.
//!
//! The one non-trivial closed form is the dilute limit of the Hashin–Shtrikman bound. As the guest phase
//! vanishes, HS with the host as matrix has to reduce to **Maxwell–Garnett** — an independently derived
//! result for a dilute suspension of spheres — and that reduction is what says the HS algebra is the HS
//! algebra rather than a plausible rational function of the right general shape.

use dualis_core::mixture::Mix;
use dualis_core::substance::{FusionProps, Substance, ThermalProps};
use dualis_core::units::{
    Density, LatentHeat, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion, Volume,
};

/// A wax, as a phase-change filler. The same numbers `substances_from_a_file.rs` marches.
fn wax() -> Substance {
    Substance::bulk("n-octadecane", Density::from_si(814.0))
        .with_thermal(ThermalProps {
            conductivity: ThermalConductivity::w_per_m_k(0.358),
            specific_heat: SpecificHeat::j_per_kg_k(1934.0),
            expansion: ThermalExpansion::ppm_per_k(800.0),
            emissivity: 0.9,
        })
        .with_fusion(FusionProps::new(
            Temperature::from_si(301.3),
            LatentHeat::kj_per_kg(244.0),
        ))
}

/// **Fractions that do not sum to one are refused, not normalised.**
///
/// Normalising silently is the tempting behaviour and it is wrong: 45% and 50% is somebody's
/// transcription mistake, and turning it into 47.4% and 52.6% answers a question nobody asked and gives
/// no sign that it happened. A mixture is data a person wrote down, and this format treats it the way
/// `deny_unknown_fields` treats a mistyped key.
#[test]
fn a_mixture_that_does_not_add_up_is_refused() {
    let (a, b) = (Substance::aluminium_6061(), Substance::copper());
    let err = Mix::of(&[(a.clone(), 0.45), (b.clone(), 0.5)]).expect_err("refused");
    assert!(
        err.contains("sum to 1") && err.contains("0.95"),
        "the message should say what they summed to: {err}"
    );

    for bad in [-0.5, 0.0, f64::NAN, f64::INFINITY] {
        assert!(
            Mix::of(&[(a.clone(), bad), (b.clone(), 1.0 - bad)]).is_err(),
            "a fraction of {bad} was accepted"
        );
    }
    assert!(Mix::of(&[]).is_err(), "an empty mixture was accepted");

    // Binary arithmetic: three thirds do not sum to one exactly, and refusing them would be refusing
    // the correct answer. `1e-9` of slack, not zero.
    let third = 1.0 / 3.0;
    assert!(Mix::of(&[(a.clone(), third), (b.clone(), third), (a, third)]).is_ok());
}

/// **A mixture of one is the substance it is made of.**
///
/// Not a special case in the code and it should not be one in the answers: every bound collapses to the
/// substance's own value and the exact properties are its own. This is the identity that says the
/// weighting is a weighting rather than something that only works at a half.
#[test]
fn a_mixture_of_one_reproduces_its_only_part() {
    for s in [
        Substance::aluminium_6061(),
        Substance::borosilicate_crown(),
        wax(),
    ] {
        let m = Mix::of(&[(s.clone(), 1.0)]).expect("one part sums to one");
        let t = s.thermal.expect("all three state thermal properties");
        assert_eq!(m.density(), s.density);
        assert_eq!(m.specific_heat(), Some(t.specific_heat));
        let (lo, hi) = m.conductivity_bounds().expect("bounds");
        assert_eq!(lo, t.conductivity, "{}: Reuss", s.name);
        assert_eq!(hi, t.conductivity, "{}: Voigt", s.name);
        assert_eq!(m.mass_fraction(0), Some(1.0));
    }
}

/// **Density and volumetric heat capacity are exact, and the specific heat is mass-weighted.**
///
/// The equalities are `1e-15`, because these are conservation statements: a cubic metre of composite
/// holds each part's share of the joules, and weighs each part's share of the kilograms. Anything looser
/// would be hiding the possibility that they are not exact, which is the thing being asserted.
///
/// The second half is the trap, quantified. Volume-weighting `c_p` is only 0.08% wrong for aluminium and
/// borosilicate, whose densities are within 8% of each other — so a caller who tested the distinction
/// there would conclude it does not matter. For copper and FR-4 it is **46%**. The test asserts both, so
/// the *range* of the error is on the record and not just its existence.
#[test]
fn the_exact_properties_are_exact_and_the_specific_heat_is_mass_weighted() {
    let cases: [(Substance, Substance, f64, f64); 3] = [
        (
            Substance::aluminium_6061(),
            Substance::borosilicate_crown(),
            0.5,
            7.8945e-4,
        ),
        (
            Substance::aluminium_6061(),
            Substance::copper(),
            0.5,
            0.27253,
        ),
        (Substance::copper(), Substance::fr4(), 0.5, 0.46337),
    ];
    for (a, b, fa, expect_trap) in cases {
        let m = Mix::of(&[(a.clone(), fa), (b.clone(), 1.0 - fa)]).expect("sums to one");
        let (ta, tb) = (a.thermal.expect("a"), b.thermal.expect("b"));

        // Density: exact.
        let want_rho = fa * a.density.to_si() + (1.0 - fa) * b.density.to_si();
        assert!(
            (m.density().to_si() / want_rho - 1.0).abs() < 1e-15,
            "{} + {}: density {} against {want_rho}",
            a.name,
            b.name,
            m.density().to_si()
        );

        // Volumetric heat capacity: exact, and this is the additive quantity.
        let want_volumetric = fa * a.density.to_si() * ta.specific_heat.to_si()
            + (1.0 - fa) * b.density.to_si() * tb.specific_heat.to_si();
        let got_volumetric = m.density().to_si() * m.specific_heat().expect("c_p").to_si();
        assert!(
            (got_volumetric / want_volumetric - 1.0).abs() < 1e-15,
            "{} + {}: rho c {got_volumetric} against {want_volumetric}",
            a.name,
            b.name
        );

        // And the mass fractions are what turn one into the other.
        let (wa, wb) = (
            m.mass_fraction(0).expect("a"),
            m.mass_fraction(1).expect("b"),
        );
        assert!((wa + wb - 1.0).abs() < 1e-15, "mass fractions sum to one");
        let by_mass = wa * ta.specific_heat.to_si() + wb * tb.specific_heat.to_si();
        assert!(
            (m.specific_heat().expect("c_p").to_si() / by_mass - 1.0).abs() < 1e-15,
            "{} + {}: c_p is the mass-weighted mean",
            a.name,
            b.name
        );

        // The trap, sized.
        let by_volume = fa * ta.specific_heat.to_si() + (1.0 - fa) * tb.specific_heat.to_si();
        let trap = (by_volume / m.specific_heat().expect("c_p").to_si() - 1.0).abs();
        println!(
            "  {:>22} + {:<14} c_p {:.1} correct, {:.1} volume-weighted — {:.2}% out",
            a.name,
            b.name,
            m.specific_heat().expect("c_p").to_si(),
            by_volume,
            trap * 100.0
        );
        // Half a percent of the measured figure, which is documentation precision: these three numbers
        // are quoted in `mixture.rs`'s own docs and this test's job is to keep them true.
        assert!(
            (trap / expect_trap - 1.0).abs() < 0.005,
            "{} + {}: the volume-weighting error is {:.4}%, was measured at {:.2}%",
            a.name,
            b.name,
            trap * 100.0,
            expect_trap * 100.0
        );

        // A volume of it weighs and holds what its parts do.
        let v = Volume::from_si(1e-6);
        assert!(
            (m.mass_of(v).to_si() / (want_rho * 1e-6) - 1.0).abs() < 1e-15,
            "the mass of a volume follows from the density"
        );
        assert!(
            (m.heat_capacity(v).expect("capacity").to_si() / (want_volumetric * 1e-6) - 1.0).abs()
                < 1e-15,
            "as does the capacity"
        );
    }
}

/// **The four bounds are always in order, and the outer pair collapses when the phases match.**
///
/// `Reuss ≤ HS− ≤ HS+ ≤ Voigt` is the theorem, checked across nine volume fractions and three pairs
/// spanning contrasts from 1.2× to 335×. Nine fractions rather than one because a bound that is ordered
/// at a half and crosses over at a tenth is a bound with a sign error in it, and that is exactly the kind
/// of mistake that a symmetric test cannot see.
///
/// Two identical phases are the degenerate case and all four have to be that conductivity: a "mixture"
/// of copper and copper is copper, and a bound that opened up there would be reporting uncertainty about
/// a homogeneous material.
#[test]
fn the_bounds_are_ordered_at_every_fraction() {
    let pairs = [
        (Substance::aluminium_6061(), Substance::copper()),
        (Substance::aluminium_6061(), Substance::borosilicate_crown()),
        (Substance::copper(), Substance::fr4()),
    ];
    for (a, b) in pairs {
        for tenths in 1..10 {
            let f = tenths as f64 / 10.0;
            let m = Mix::of(&[(a.clone(), f), (b.clone(), 1.0 - f)]).expect("sums to one");
            let (reuss, voigt) = m.conductivity_bounds().expect("bounds");
            let (hs_lo, hs_hi) = m.hashin_shtrikman().expect("two parts");
            assert!(
                reuss.to_si() < hs_lo.to_si()
                    && hs_lo.to_si() <= hs_hi.to_si()
                    && hs_hi.to_si() < voigt.to_si(),
                "{} {f} + {}: {} <= {} <= {} <= {} is out of order",
                a.name,
                b.name,
                reuss.to_si(),
                hs_lo.to_si(),
                hs_hi.to_si(),
                voigt.to_si()
            );
            // And both bounds sit between the two conductivities, which no weighted mean can leave.
            let (ka, kb) = (
                a.thermal.expect("a").conductivity.to_si(),
                b.thermal.expect("b").conductivity.to_si(),
            );
            let (small, large) = (ka.min(kb), ka.max(kb));
            assert!(reuss.to_si() > small && voigt.to_si() < large);
        }
    }

    // Two of the same thing.
    let cu = Substance::copper();
    let k = cu.thermal.expect("copper").conductivity;
    let same = Mix::of(&[(cu.clone(), 0.3), (cu, 0.7)]).expect("sums to one");
    let (reuss, voigt) = same.conductivity_bounds().expect("bounds");
    let (hs_lo, hs_hi) = same.hashin_shtrikman().expect("two parts");
    for (what, got) in [
        ("Reuss", reuss),
        ("Voigt", voigt),
        ("HS-", hs_lo),
        ("HS+", hs_hi),
    ] {
        assert!(
            (got.to_si() / k.to_si() - 1.0).abs() < 1e-15,
            "{what}: a mixture of copper and copper conducts {} and not {}",
            got.to_si(),
            k.to_si()
        );
    }
}

/// **The Hashin–Shtrikman bound reduces to Maxwell–Garnett as the guest phase vanishes.**
///
/// The independent check on the algebra, and the reason it is worth having: the HS expression is a
/// rational function and several plausible wrong versions of it are also rational functions with the
/// right limits at zero and one. Maxwell–Garnett is a separately derived result for a dilute suspension
/// of spheres in a matrix,
///
/// ```text
///   k = k_h [ 1 + 3φ(k_g − k_h) / (k_g + 2k_h − φ(k_g − k_h)) ]
/// ```
///
/// and the HS bound taken with the matrix as host has to *be* it — not approach it, be it, because the
/// coated-sphere assemblage that attains the bound is what Maxwell–Garnett describes. So this is an
/// equality at every fraction rather than a limit, which is a far stronger statement than a dilute
/// check would be, and it holds to `1e-14` from 0.1% filler to 90%.
#[test]
fn hashin_shtrikman_is_maxwell_garnett_with_the_matrix_as_host() {
    let matrix = Substance::fr4();
    let filler = Substance::copper();
    let kh = matrix.thermal.expect("matrix").conductivity.to_si();
    let kg = filler.thermal.expect("filler").conductivity.to_si();

    let mut worst = 0.0f64;
    for phi in [0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9] {
        let m =
            Mix::of(&[(matrix.clone(), 1.0 - phi), (filler.clone(), phi)]).expect("sums to one");
        // The *lower* HS bound is the one built around the low-conductivity phase, which is the matrix
        // here — so it is the one Maxwell–Garnett describes.
        let (hs_lo, _) = m.hashin_shtrikman().expect("two parts");
        let mg = kh * (1.0 + 3.0 * phi * (kg - kh) / (kg + 2.0 * kh - phi * (kg - kh)));
        let off = (hs_lo.to_si() / mg - 1.0).abs();
        println!(
            "  phi {phi:<6} HS- {:.6} against Maxwell-Garnett {mg:.6} — off {off:.2e}",
            hs_lo.to_si()
        );
        assert!(
            off < 1e-14,
            "phi {phi}: HS- is {} and Maxwell-Garnett is {mg}",
            hs_lo.to_si()
        );
        worst = worst.max(off);
    }
    println!("  worst {worst:.2e} over three decades of filler fraction");
}

/// **A phase-change composite's latent heat is exact, and it is the mass fraction that makes it so.**
///
/// Wax at 814 kg/m³ filling 80% of the volume of an aluminium matrix is 54.67% of the mass, so a wax
/// storing 244 kJ/kg gives a composite storing 133.4. Using the volume fraction instead would claim
/// 195.2 — **46.7% high**, on the single property a thermal buffer is bought for.
///
/// Both are asserted: the correct value as an equality, and the size of the error the other way, so the
/// trap is on the record rather than merely avoided.
#[test]
fn a_phase_change_composite_dilutes_its_latent_heat_by_mass() {
    let m = Mix::of(&[(wax(), 0.8), (Substance::aluminium_6061(), 0.2)]).expect("sums to one");
    let (point, latent) = m
        .fusion()
        .expect("one part melts")
        .expect("so there is a fusion");

    assert_eq!(
        point,
        Temperature::from_si(301.3),
        "the melting point is the melting part's, undiluted — diluting a temperature is meaningless"
    );
    let w = m.mass_fraction(0).expect("the wax");
    assert!(
        (latent.to_si() / (w * 244_000.0) - 1.0).abs() < 1e-15,
        "latent heat is the mass fraction times the wax's: {} against {}",
        latent.to_si(),
        w * 244_000.0
    );
    println!(
        "  wax 80% by volume is {:.2}% by mass: {:.1} kJ/kg, against {:.1} if volume were used \
         — {:.1}% high",
        w * 100.0,
        latent.to_si() / 1000.0,
        0.8 * 244.0,
        (0.8 * 244_000.0 / latent.to_si() - 1.0) * 100.0
    );
    assert!(
        ((0.8 * 244_000.0 / latent.to_si() - 1.0) - 0.4634).abs() < 0.001,
        "the volume-fraction mistake is 46.3% here"
    );

    // Nothing melts: `None`, and that is not an error.
    let inert = Mix::of(&[
        (Substance::aluminium_6061(), 0.5),
        (Substance::copper(), 0.5),
    ])
    .expect("sums to one");
    assert_eq!(inert.fusion(), Ok(None));

    // Two melt: refused, because two plateaux is not one melting point and there is no ordering of the
    // two that would be right.
    let both = Mix::of(&[(wax(), 0.5), (Substance::ice(), 0.5)]).expect("sums to one");
    let err = both.fusion().expect_err("refused");
    assert!(
        err.contains("octadecane") && err.contains("ice") && err.contains("melting point"),
        "the message should name both and say what the problem is: {err}"
    );
}

/// **`as_substance` refuses a conductivity no microstructure realises, and accepts both bounds.**
///
/// The un-skippable part of the design. A caller has to choose a conductivity because the library cannot,
/// and the one thing the library *can* do is refuse an impossible choice — so the check lives in the
/// constructor rather than in a separate validator somebody may not call.
///
/// Both endpoints are allowed, because both are attained: `a_composite.rs` measures a laminate hitting
/// each. A hair outside either is refused. And the emissivity is an argument for the reason the doc
/// gives — it is a surface property and a mixture has no surface — so an impossible one is refused here
/// too rather than clamped.
#[test]
fn a_composite_substance_refuses_an_impossible_conductivity() {
    let m = Mix::of(&[
        (Substance::aluminium_6061(), 0.5),
        (Substance::borosilicate_crown(), 0.5),
    ])
    .expect("sums to one");
    let (reuss, voigt) = m.conductivity_bounds().expect("bounds");

    // Both ends, and the midpoint of the tighter pair, which is what a foam would use.
    let (hs_lo, hs_hi) = m.hashin_shtrikman().expect("two parts");
    for k in [
        reuss,
        voigt,
        ThermalConductivity::from_si(0.5 * (hs_lo.to_si() + hs_hi.to_si())),
    ] {
        let s = m
            .as_substance("Al/BK7", k, 0.9)
            .unwrap_or_else(|e| panic!("{} should be allowed: {e}", k.to_si()));
        let t = s.thermal.expect("it states thermal properties");
        assert_eq!(t.conductivity, k, "the choice is carried through");
        assert_eq!(s.density, m.density());
        assert_eq!(Some(t.specific_heat), m.specific_heat());
        // And the substance's own accessor agrees, so a domain reading it gets the same diffusivity.
        let alpha = k.to_si() / (m.density().to_si() * m.specific_heat().expect("c_p").to_si());
        assert!((s.diffusivity().expect("alpha").to_si() / alpha - 1.0).abs() < 1e-15);
        assert!(s.check().is_ok(), "{}", s.check().unwrap_err());
    }

    for outside in [
        reuss.to_si() * 0.999,
        voigt.to_si() * 1.001,
        0.0,
        -1.0,
        f64::NAN,
    ] {
        let err = m
            .as_substance("Al/BK7", ThermalConductivity::from_si(outside), 0.9)
            .expect_err("refused");
        assert!(err.contains("no microstructure"), "{outside}: {err}");
    }

    for bad in [-0.1, 1.1] {
        assert!(
            m.as_substance("Al/BK7", voigt, bad).is_err(),
            "an emissivity of {bad} was accepted"
        );
    }

    // A part with no thermal properties makes the whole thing unknown rather than partially known.
    let unknown = Mix::of(&[
        (Substance::bulk("mystery", Density::g_per_cm3(2.0)), 0.5),
        (Substance::copper(), 0.5),
    ])
    .expect("sums to one");
    assert_eq!(unknown.conductivity_bounds(), None);
    assert_eq!(unknown.specific_heat(), None);
    assert!(unknown.as_substance("x", voigt, 0.9).is_err());
    // But the density is still exact, because density is not optional on a substance.
    assert!((unknown.density().to_si() - 0.5 * (2000.0 + 8960.0)).abs() < 1e-12);
}

/// **Mixing in two steps gives the same answer as mixing in one — for the exact properties and the outer
/// bounds, and not for Hashin–Shtrikman.**
///
/// A composite of a composite is a real thing: a filled resin, itself a fraction of a laminate. So the
/// weighting has to be associative, and this is the test that catches a fraction applied at the wrong
/// level — which is the mistake that has the same shape as the volume-and-mass confusion and would
/// otherwise show up as a plausible number.
///
/// It holds exactly for density, volumetric heat capacity and both outer bounds, because all four are
/// weighted means of the same form. It does **not** hold for HS, and that is not a defect: HS is derived
/// for a *two*-phase isotropic medium, so applying it to a phase that is itself a composite is asking a
/// different question. Asserted as an inequality so the distinction stays visible.
#[test]
fn mixing_in_two_steps_matches_mixing_in_one() {
    let (a, b, c) = (
        Substance::aluminium_6061(),
        Substance::borosilicate_crown(),
        Substance::copper(),
    );
    // A and B half and half, then that at 40% with C at 60%. So A and B are 20% each of the whole.
    let inner = Mix::of(&[(a.clone(), 0.5), (b.clone(), 0.5)]).expect("sums to one");
    let (reuss_in, voigt_in) = inner.conductivity_bounds().expect("bounds");
    let flat = Mix::of(&[(a, 0.2), (b, 0.2), (c.clone(), 0.6)]).expect("sums to one");

    // The inner mixture as a substance, at each of its own bounds, so the outer answer can be compared
    // with the flat one's corresponding bound.
    for (k_inner, pick) in [(reuss_in, "Reuss"), (voigt_in, "Voigt")] {
        let composite = inner
            .as_substance("Al/BK7", k_inner, 0.9)
            .expect("a bound is allowed");
        let nested = Mix::of(&[(composite, 0.4), (c.clone(), 0.6)]).expect("sums to one");

        assert!(
            (nested.density().to_si() / flat.density().to_si() - 1.0).abs() < 1e-15,
            "{pick}: density is associative"
        );
        let volumetric = |m: &Mix| m.density().to_si() * m.specific_heat().expect("c_p").to_si();
        assert!(
            (volumetric(&nested) / volumetric(&flat) - 1.0).abs() < 1e-15,
            "{pick}: volumetric heat capacity is associative"
        );

        let (reuss_n, voigt_n) = nested.conductivity_bounds().expect("bounds");
        let (reuss_f, voigt_f) = flat.conductivity_bounds().expect("bounds");
        let want = if pick == "Reuss" { reuss_f } else { voigt_f };
        let got = if pick == "Reuss" { reuss_n } else { voigt_n };
        assert!(
            (got.to_si() / want.to_si() - 1.0).abs() < 1e-15,
            "{pick}: nesting gives {} and flattening gives {}",
            got.to_si(),
            want.to_si()
        );
    }

    // And HS is not associative, which is a property of what HS assumes rather than of the arithmetic.
    let composite = inner
        .as_substance("Al/BK7", voigt_in, 0.9)
        .expect("allowed");
    let nested = Mix::of(&[(composite, 0.4), (c, 0.6)]).expect("sums to one");
    let (nested_lo, _) = nested.hashin_shtrikman().expect("two parts");
    // The flat mixture has three parts, so it has no HS pair at all — which is the honest answer.
    assert_eq!(
        flat.hashin_shtrikman(),
        None,
        "HS is a two-phase result and a three-part mixture is not two phases"
    );
    println!(
        "  nesting gives an HS- of {:.4}; flattening gives no HS pair, which is the correct answer",
        nested_lo.to_si()
    );
}
