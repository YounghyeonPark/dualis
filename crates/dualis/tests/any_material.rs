//! A material that is not in the catalogue, applied to every domain that takes one.
//!
//! The catalogue holds nine entries and there are hundreds of thousands of materials. Enumeration
//! never closes that gap, so the question is whether a caller with a datasheet can express theirs and
//! get the same treatment — and until `with_thermal`, `with_mechanical`, `with_acoustic` and
//! `with_fusion` existed, they could not: `Substance::bulk` sets all four blocks to `None` and there
//! was no way to fill any of them except a struct literal naming every field.
//!
//! That literal path is the one that breaks. Adding `fusion` for latent heat stopped every literal
//! outside the crate from compiling, and the callers it cost were exactly the ones the catalogue was
//! least able to help.
//!
//! # What the library can still do for a number it did not choose
//!
//! Not check that it is right. `check` checks that it is **possible** — bounds on each field, and one
//! cross-check that is not a bound: a substance stating both a sound speed and elastic constants has
//! three independent numbers describing one thing, and a longitudinal wave has to sit near the rod
//! speed `sqrt(E/rho)` or the bulk speed `sqrt((lambda+2mu)/rho)` those give.
//!
//! An impossible material otherwise produces an answer that is plausible and wrong, which is the
//! failure this workspace is organised around not having.

use dualis::prelude::*;
use dualis::units::LatentHeat;
use dualis_core::substance::{AcousticProps, FusionProps, MechanicalProps, ThermalProps};
use dualis_elastic::{Axis, Elastic, Waves};

/// Ti-6Al-4V, from a datasheet. Not in the catalogue and it does not need to be.
fn titanium() -> Substance {
    Substance::bulk("Ti-6Al-4V", Density::g_per_cm3(4.43))
        .with_thermal(ThermalProps {
            conductivity: ThermalConductivity::w_per_m_k(6.7),
            specific_heat: SpecificHeat::j_per_kg_k(526.0),
            expansion: ThermalExpansion::ppm_per_k(8.6),
            emissivity: 0.30,
        })
        .with_mechanical(MechanicalProps {
            youngs_modulus: Pressure::from_si(113.8e9),
            poisson_ratio: 0.342,
            yield_strength: Pressure::from_si(880.0e6),
        })
        .with_acoustic(AcousticProps {
            sound_speed: Velocity::m_per_s(6100.0),
        })
}

/// Paraffin wax, which is what a phase-change thermal buffer actually is.
fn paraffin() -> Substance {
    Substance::bulk("paraffin RT44", Density::g_per_cm3(0.80))
        .with_thermal(ThermalProps {
            conductivity: ThermalConductivity::w_per_m_k(0.2),
            specific_heat: SpecificHeat::j_per_kg_k(2000.0),
            expansion: ThermalExpansion::ppm_per_k(300.0),
            emissivity: 0.90,
        })
        .with_fusion(FusionProps::new(
            Temperature::celsius(44.0),
            LatentHeat::kj_per_kg(250.0),
        ))
}

/// **Every entry in the catalogue passes its own validator.**
///
/// The first thing `check` has to be is true of the nine materials this crate is answerable for. If it
/// were not, the bounds would be wrong rather than the materials — and a validator that the shipped
/// data fails is a validator nobody will keep.
#[test]
fn the_catalogue_passes_its_own_check() {
    for s in [
        Substance::aluminium_6061(),
        Substance::borosilicate_crown(),
        Substance::stainless_304(),
        Substance::copper(),
        Substance::fr4(),
        Substance::electrical_steel(),
        Substance::pla(),
        Substance::water(),
        Substance::ice(),
    ] {
        let name = s.name.clone();
        assert!(s.check().is_ok(), "{name}: {}", s.check().unwrap_err());
    }
    // And so do the two materials this file invents, which is the point of the file.
    for s in [titanium(), paraffin()] {
        assert!(s.check().is_ok(), "{}", s.check().unwrap_err());
    }
    // A substance known only by its density passes too: unknown is not wrong, and every block it
    // could not fill in is honestly absent rather than zero.
    assert!(Substance::bulk("mystery", Density::g_per_cm3(2.0))
        .check()
        .is_ok());
}

/// **An impossible material is refused, and the message names the field.**
///
/// Each of these is a transcription mistake somebody makes: a percentage where a fraction belongs, a
/// sign, a column read one across. Every one of them would otherwise run and produce an answer.
///
/// The last is the interesting one. Its numbers are individually plausible — an aluminium density with
/// a steel modulus — and only the **cross-check** catches it: those two do not give a sound speed
/// anywhere near the one stated beside them.
#[test]
fn an_impossible_material_is_refused_by_the_field_that_is_wrong() {
    let base = || Substance::bulk("suspect", Density::g_per_cm3(2.7));
    let cases: Vec<(&str, Substance, &str)> = vec![
        (
            "emissivity as a percentage",
            base().with_thermal(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(167.0),
                specific_heat: SpecificHeat::j_per_kg_k(896.0),
                expansion: ThermalExpansion::ppm_per_k(23.6),
                emissivity: 9.0,
            }),
            "emissivity",
        ),
        (
            "a negative conductivity",
            base().with_thermal(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(-167.0),
                specific_heat: SpecificHeat::j_per_kg_k(896.0),
                expansion: ThermalExpansion::ppm_per_k(23.6),
                emissivity: 0.1,
            }),
            "conductivity",
        ),
        (
            "an incompressible Poisson ratio",
            base().with_mechanical(MechanicalProps {
                youngs_modulus: Pressure::from_si(68.9e9),
                poisson_ratio: 0.5,
                yield_strength: Pressure::from_si(276.0e6),
            }),
            "poisson_ratio",
        ),
        (
            "a latent heat of nothing",
            base().with_fusion(FusionProps::new(
                Temperature::celsius(44.0),
                LatentHeat::kj_per_kg(0.0),
            )),
            "latent_heat",
        ),
        (
            "a shear speed where the longitudinal one belongs",
            base()
                .with_mechanical(MechanicalProps {
                    youngs_modulus: Pressure::from_si(68.9e9),
                    poisson_ratio: 0.33,
                    yield_strength: Pressure::from_si(276.0e6),
                })
                .with_acoustic(AcousticProps {
                    // 6061's shear speed, not its longitudinal one.
                    sound_speed: Velocity::m_per_s(3100.0),
                }),
            "sound_speed",
        ),
        (
            "a modulus from the row below",
            base()
                .with_mechanical(MechanicalProps {
                    youngs_modulus: Pressure::from_si(200.0e9),
                    poisson_ratio: 0.29,
                    yield_strength: Pressure::from_si(276.0e6),
                })
                .with_acoustic(AcousticProps {
                    sound_speed: Velocity::m_per_s(6320.0),
                }),
            "sound_speed",
        ),
    ];
    for (why, s, field) in cases {
        let Err(message) = s.check() else {
            panic!("{why} must be refused");
        };
        println!("  {why:44} -> {message}");
        assert!(
            message.contains(field),
            "{why}: the message must name {field}, says {message}"
        );
    }

    // Several wrong at once are all reported, because a material read off the wrong column usually is.
    let doomed = Substance::bulk("doomed", Density::from_si(-1.0)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(0.0),
        specific_heat: SpecificHeat::j_per_kg_k(-1.0),
        expansion: ThermalExpansion::ppm_per_k(1.0),
        emissivity: -0.5,
    });
    let message = doomed.check().expect_err("four problems");
    let named = ["density", "conductivity", "specific_heat", "emissivity"]
        .iter()
        .filter(|f| message.contains(**f))
        .count();
    println!("  four problems at once, {named} named: {message}");
    assert_eq!(named, 4, "all of them, not the first: {message}");
}

/// **A material nobody in this crate chose runs in three domains and gets the same treatment.**
///
/// Titanium through conduction and through elastic waves, and paraffin through a phase change. The
/// point is not that the numbers are special — it is that nothing along the way asks whether the
/// substance came from the catalogue.
///
/// Each is checked against a closed form rather than against a picture, on the same footing as a
/// catalogue material: the stability limit, the wave speed ratio, and the melting plateau.
#[test]
fn a_datasheet_material_works_in_every_domain_that_takes_one() {
    let ti = titanium();

    // Conduction. The limit is `dx^2/(6 alpha)` and titanium's diffusivity is a third of stainless
    // steel's, so the step it allows is a number about the material and not about the catalogue.
    let block = dualis::thermal::Solid3D::new(
        "ti",
        ti.clone(),
        (5, 5, 5),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    );
    let alpha = ti.diffusivity().expect("it has thermal properties").to_si();
    let limit = block.max_stable_dt(Time::from_si(0.0)).to_si();
    let closed = 1e-6 / (6.0 * alpha);
    println!(
        "  Ti: alpha {alpha:.4e} m2/s, conduction limit {limit:.4e} s against dx^2/6alpha {closed:.4e}"
    );
    assert!(
        (limit / closed - 1.0).abs() < 1e-12,
        "a datasheet material gets the same limit formula: {limit:.6e} against {closed:.6e}"
    );

    // Elastic waves. The ratio is `sqrt(2(1-nu)/(1-2nu))` and nu = 0.342 is a value no catalogue entry
    // has, so this cannot be passing on a memorised number.
    let e = Elastic::from_substance(&ti).expect("it has mechanical properties");
    let ratio = e.p_wave_speed().to_si() / e.s_wave_speed().to_si();
    println!(
        "  Ti: c_p {:.0} m/s, c_s {:.0}, ratio {ratio:.6} against sqrt(2(1-v)/(1-2v)) {:.6}",
        e.p_wave_speed().to_si(),
        e.s_wave_speed().to_si(),
        e.speed_ratio()
    );
    assert!(
        (ratio / e.speed_ratio() - 1.0).abs() < 1e-14,
        "the identity holds for a material this crate never saw"
    );
    let mut w = Waves::new("ti bar", (1, 1, 16), Length::mm(1.0), e);
    w.hold(Axis::X);
    w.hold(Axis::Y);
    w.clamp_ends(Axis::Z);
    w.release_mode(1, Axis::Z, Axis::Z, Length::from_si(1e-9));
    let dt = Time::from_si(w.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    for n in 0..200 {
        w.step(
            Time::from_si(n as f64 * dt.to_si()),
            dt,
            &mut Exchange::new(),
        )
        .expect("a datasheet material is as stable as a catalogue one");
    }
    assert!(
        w.strain_energy().to_si() > 0.0,
        "and it is actually ringing"
    );

    // A phase change. Paraffin's plateau is `rho L V / P`, exactly as ice's is.
    let wax = paraffin();
    let mut buffer = dualis::thermal::Solid3D::new(
        "wax",
        wax.clone(),
        (1, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(44.0),
    );
    let volume = Volume::from_si(1e-9);
    let latent = wax.latent_energy(volume).expect("it melts").to_si();
    let power = 1e-3;
    let (dt, mut t, mut left_at) = (0.05, 0.0, None);
    while t < 2.0 * latent / power {
        buffer.deposit(0, 0, 0, Energy::from_si(power * dt));
        t += dt;
        if left_at.is_none()
            && buffer.temperature_at(0, 0, 0).to_si() > Temperature::celsius(44.0).to_si()
        {
            left_at = Some(t);
        }
    }
    let left = left_at.expect("it warms eventually");
    println!(
        "  paraffin: held at 44 C for {left:.2} s against rho L V / P = {:.2}",
        latent / power
    );
    assert!(
        (left - latent / power).abs() <= dt * 1.001,
        "the plateau is the latent heat over the power, for a wax as much as for ice"
    );
}

/// **A material survives being written to JSON and read back, so a file of them is a supported path.**
///
/// The route to "every material" that does not involve this crate learning any of them. `Substance`
/// has derived `Serialize` and `Deserialize` all along; what was missing was any statement that the
/// round trip is exact, and a format nobody checks is a format that drifts.
///
/// Asserted from the parsed value rather than from a second serialisation — comparing
/// `to_string(parse(x))` with `to_string(parse(to_string(parse(x))))` puts serialiser output on both
/// sides and cannot see a field the parser dropped. `dualis-world` learned that one the hard way.
#[test]
fn a_material_round_trips_through_json_so_a_file_can_hold_a_catalogue() {
    for original in [titanium(), paraffin(), Substance::ice()] {
        let text = serde_json::to_string(&original).expect("it serialises");
        let back: Substance = serde_json::from_str(&text).expect("and parses");
        assert_eq!(back, original, "the round trip is exact");
        assert!(back.check().is_ok(), "and still valid");
    }

    // Hand-written, which is what a file actually is — and the absent blocks stay absent rather than
    // arriving as zeros.
    let text = r#"{
        "name": "PEEK",
        "density": 1320.0,
        "thermal": { "conductivity": 0.25, "specific_heat": 1340.0,
                     "expansion": 4.7e-5, "emissivity": 0.9 }
    }"#;
    let peek: Substance = serde_json::from_str(text).expect("a partial material parses");
    assert_eq!(peek.name, "PEEK");
    assert!(peek.check().is_ok(), "{}", peek.check().unwrap_err());
    assert!(
        peek.mechanical.is_none() && peek.acoustic.is_none() && peek.fusion.is_none(),
        "what the file did not say is absent, not zero"
    );
    assert!(
        peek.diffusivity().is_some(),
        "and what it did say is usable: {:?}",
        peek.diffusivity()
    );

    // A key this format does not know is refused rather than dropped, which is the same rule the
    // scene format has and for the same reason: a discarded key makes the file say something the code
    // will not do.
    let typo = r#"{ "name": "x", "density": 1000.0, "thermalz": {} }"#;
    assert!(
        serde_json::from_str::<Substance>(typo).is_err(),
        "an unknown key must not be silently dropped"
    );
}
