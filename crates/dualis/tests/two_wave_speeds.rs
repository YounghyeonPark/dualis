//! A solid has two wave speeds and `AcousticProps` has room for one.
//!
//! `dualis-elastic` computes `c_p = √((λ+2μ)/ρ)` and `c_s = √(μ/ρ)` from the two constants a
//! material is stated with. `dualis-core`'s catalogue carries a single `sound_speed` per substance,
//! entered independently of the elastic constants beside it. This is where the two meet, and it is
//! the fourth test of that shape — after `fields_and_rays`, `loss_and_lumps` and `a_slit`.
//!
//! # What it found: one field with two meanings
//!
//! A fluid has one speed because it has no shear modulus to carry a second. A solid has two
//! longitudinal ones — the **bulk** wave `√((λ+2μ)/ρ)`, constrained by the material around it, and the
//! **rod** wave `√(E/ρ)`, free to bulge sideways — and one number cannot say which is meant.
//!
//! The six entries turn out not to mean the same one, and the split is clean:
//!
//! ```text
//!                     rod   stated    bulk    vs bulk   vs rod
//!   ice              3150     3840    3834      +0.1%   +21.9%     bulk
//!   Al 6061          5052     6320    6149      +2.8%   +25.1%     bulk
//!   304 stainless    4912     5790    5623      +3.0%   +17.9%     bulk
//!   Cu ETP           3614     4760    4483      +6.2%   +31.7%     bulk
//!   N-BK7            5716     5680    6048      −6.1%    −0.6%     rod
//!   elec. steel      5113     5100    5853     −12.9%    −0.3%     rod
//! ```
//!
//! Four are the bulk wave and two are the rod wave, and every one of the six is within 6.2% of
//! whichever it is. Nothing in [`AcousticProps`](dualis_core::substance::AcousticProps) records which,
//! because for the fluid the type was designed around there is only one.
//!
//! # Why the four are not exact, and ice is
//!
//! Read as a bulk wave, a stated speed implies a modulus. Against the catalogued `E`:
//!
//! ```text
//!   ice               9.1 against  9.1 GPa      exact
//!   Al 6061          72.8 against 68.9          +5.7%
//!   304 stainless   204.7 against 193.0         +6.1%
//!   Cu ETP          131.9 against 117.0        +12.7%
//! ```
//!
//! **A tensile test and an ultrasonic measurement are not the same measurement.** The static modulus
//! includes whatever the specimen does inelastically at low stress; the dynamic one is taken at
//! megahertz and small strain, and for annealed metals it comes out higher. Copper is the softest of
//! the three and the furthest apart, which is the direction that says so.
//!
//! Ice agrees exactly because both of its numbers came from the same acoustic experiment — its
//! elastic constants *are* back-calculated from velocity. So the one entry that agrees is the one
//! where agreement was never independent, and the three that disagree are the three where it was.
//!
//! # What is therefore asserted
//!
//! The `ν`-only identity, at machine precision, which is about the arithmetic and not the catalogue.
//! Ice's two routes agreeing to 0.14%, which is the one place they can be compared. And that every
//! entry sits within 7% of one of the two speeds — a bound that a shear speed, a fluid's, or a
//! slipped decimal all fail, and which is as tight as independently sourced data allows.
//!
use dualis::prelude::*;
use dualis_elastic::Elastic;

/// The substances that carry both an elastic description and a sound speed.
fn both() -> Vec<(&'static str, Substance)> {
    vec![
        ("ice", Substance::ice()),
        ("Al 6061", Substance::aluminium_6061()),
        ("304 stainless", Substance::stainless_304()),
        ("Cu ETP", Substance::copper()),
        ("N-BK7", Substance::borosilicate_crown()),
        ("electrical steel", Substance::electrical_steel()),
    ]
}

fn elastic_of(s: &Substance) -> Elastic {
    let m = s.mechanical.expect("these all have mechanical properties");
    Elastic::new(m.youngs_modulus, m.poisson_ratio, s.density)
        .expect("a catalogue entry is a representable material")
}

/// **The speed ratio is Poisson's ratio and nothing else, for every material in the catalogue.**
///
/// `c_p/c_s = √(2(1−ν)/(1−2ν))`. Both `E` and `ρ` cancel out of it, which is what makes it the
/// sharpest statement available here: it is checked against `speed_ratio`, which is written from `ν`
/// alone and shares no arithmetic with the two speeds it is compared to.
///
/// A machine-precision equality, because it is an algebraic identity rather than a discretisation.
#[test]
fn the_ratio_of_the_two_speeds_is_poissons_ratio_alone() {
    for (name, s) in both() {
        let e = elastic_of(&s);
        let (cp, cs) = (e.p_wave_speed().to_si(), e.s_wave_speed().to_si());
        let from_nu = e.speed_ratio();
        let off = (cp / cs / from_nu - 1.0).abs();
        println!(
            "  {name:17} c_p {cp:6.0}  c_s {cs:6.0}  ratio {:.4} against sqrt(2(1-v)/(1-2v)) \
             {from_nu:.4} — off {off:.1e}",
            cp / cs
        );
        assert!(
            off < 1e-14,
            "{name}: the ratio is an identity in nu: {:.9} against {from_nu:.9}",
            cp / cs
        );
        // And it is greater than √2 for every real material, which is the ν > 0 half of the range.
        assert!(
            cp / cs > std::f64::consts::SQRT_2,
            "{name}: a positive Poisson ratio puts the ratio above sqrt(2): {:.4}",
            cp / cs
        );
    }
}

/// **A shear wave is slower than a pressure wave, always, and a fluid has no shear wave at all.**
///
/// `c_p > c_s` for every `ν > −1`, which is why the P in P-wave means *primary*: it is the first
/// arrival at a seismometer, and the gap between the two is how the distance to the source is got.
///
/// The fluid half is the one that says the two speeds are about the material and not the solver:
/// water has no `mechanical` entry at all, so there is no `μ` to make a shear wave out of, and the
/// catalogue declines rather than reporting zero.
#[test]
fn a_fluid_has_no_shear_modulus_to_carry_a_second_wave() {
    for (name, s) in both() {
        let e = elastic_of(&s);
        assert!(
            e.p_wave_speed().to_si() > e.s_wave_speed().to_si(),
            "{name}: the pressure wave is the primary arrival"
        );
    }
    let water = Substance::water();
    assert!(
        water.mechanical.is_none(),
        "water has no shear modulus, and the catalogue says so by declining rather than by zero"
    );
    assert!(
        water.acoustic.is_some(),
        "but it does have a sound speed, which is the whole asymmetry: one number is enough \
         for a fluid"
    );
}

/// **Every catalogued `sound_speed` is one of the two longitudinal speeds, to within 7%.**
///
/// The cross-domain check, and it is a *classification* rather than an equality because measuring it
/// showed the field carries two different quantities. Four entries are the bulk wave and two are the
/// rod wave; the test works out which each is nearer and requires the gap to be small.
///
/// That is a real constraint. A shear speed would be 40% below the rod speed, a fluid's further, and
/// a slipped decimal further still — all fail. What it cannot do is be tighter, because the numbers
/// either side of it come from different experiments: see this file's header for the static-versus-
/// dynamic modulus gap that puts copper 6.2% out.
///
/// Two earlier drafts asserted `rod ≤ stated ≤ bulk` with 0.1% and then 2% of slack. Both were false
/// — of six entries, three sit *above* the bulk speed and two *below* the rod speed — and the second
/// draft failed on aluminium, which is not an edge case.
#[test]
fn every_catalogued_sound_speed_is_one_of_the_two_longitudinal_speeds() {
    let mut worst: f64 = 0.0;
    let mut kinds = Vec::new();
    for (name, s) in both() {
        let e = elastic_of(&s);
        let m = s.mechanical.expect("mechanical");
        let rod = (m.youngs_modulus.to_si() / s.density.to_si()).sqrt();
        let bulk = e.p_wave_speed().to_si();
        let stated = s.acoustic.expect("acoustic").sound_speed.to_si();
        let (d_rod, d_bulk) = ((stated / rod - 1.0), (stated / bulk - 1.0));
        let (kind, gap) = if d_bulk.abs() < d_rod.abs() {
            ("bulk", d_bulk)
        } else {
            ("rod", d_rod)
        };
        println!(
            "  {name:17} rod {rod:6.0}  stated {stated:6.0}  bulk {bulk:6.0}   \
             {:+6.1}% of bulk {:+7.1}% of rod   -> {kind} to {:.1}%",
            d_bulk * 100.0,
            d_rod * 100.0,
            gap.abs() * 100.0
        );
        assert!(
            gap.abs() < 0.07,
            "{name}: {stated:.0} m/s is {:.1}% from the nearer of {rod:.0} and {bulk:.0} — which is \
             a speed that is not a longitudinal wave in this material, or an elastic constant that \
             is not this material's",
            gap.abs() * 100.0
        );
        worst = worst.max(gap.abs());
        kinds.push(kind);
    }
    println!(
        "  worst gap {:.1}%, and the field means both things",
        worst * 100.0
    );

    // **The split is the finding, so it is asserted.** If a future edit made all six the same kind,
    // this file's whole point would have quietly gone away — and the honest response to that would be
    // to record which wave the field means, not to delete the test.
    assert!(
        kinds.contains(&"bulk") && kinds.contains(&"rod"),
        "one field, two meanings: {kinds:?}"
    );
}

/// **Ice is the one entry where the two crates agree to a tenth of a percent, and that is the check.**
///
/// The others are bracketed; this one is *equal*, which is what says the elastic route and the
/// acoustic route are computing the same physical quantity rather than two plausible ones. `E`, `ν`
/// and `ρ` go in one side, a measured velocity comes out the other, and 3834 against 3840 m/s is
/// three independent numbers meeting a fourth.
///
/// It matters that this is ice and not a metal: ice is isotropic polycrystalline and its elastic
/// constants come from the same acoustic measurement its sound speed does, so the two are *supposed*
/// to agree. Where they do not — copper, electrical steel — the numbers came from different
/// experiments on differently processed material, and that is the honest reading rather than an
/// error in either crate.
#[test]
fn for_ice_the_two_routes_to_a_wave_speed_agree() {
    let ice = Substance::ice();
    let e = elastic_of(&ice);
    let computed = e.p_wave_speed().to_si();
    let stated = ice.acoustic.expect("acoustic").sound_speed.to_si();
    let off = (computed / stated - 1.0).abs();
    println!(
        "  ice: E, nu and rho give {computed:.0} m/s; the table states {stated:.0} — off {:.2}%",
        off * 100.0
    );
    assert!(
        off < 3e-3,
        "the elastic and acoustic descriptions of ice are the same wave: {computed:.1} against \
         {stated:.1} m/s"
    );
    // And the shear wave is about half of it, which is the ν = 0.33 signature.
    let ratio = computed / e.s_wave_speed().to_si();
    println!("  and c_p/c_s = {ratio:.3}, which nu = 0.33 alone predicts");
    assert!(
        (ratio - 1.985).abs() < 1e-3,
        "nu = 0.33 gives 1.985: {ratio:.4}"
    );
}
