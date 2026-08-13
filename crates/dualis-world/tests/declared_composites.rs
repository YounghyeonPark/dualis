//! A scene that mixes its own materials.
//!
//! `declared_materials.rs` covers a scene bringing a substance. This covers the other half, which needed
//! `Mix` to exist: a scene bringing a **composite** — a motor that is copper, steel, magnets and air, a
//! board that is FR-4 and copper, a buffer that is wax in a metal matrix. Each wants to be one material a
//! coarse grid can hold, and its properties are not the properties of its main constituent.
//!
//! # What the format can check and what it cannot
//!
//! It cannot compute the conductivity. No single value exists for a composite without the microstructure,
//! and that is physics rather than a gap — a laminate of the same two materials in the same proportions
//! conducts 38 times more one way than the other.
//!
//! What it *can* do is refuse an impossible one, and the refusal is where the value is: a scene file has
//! nowhere else to learn what range is achievable, so the message carries the Voigt and Reuss bounds and
//! the tighter Hashin–Shtrikman pair. A caller who guessed is told the four numbers that matter.
//!
//! It cannot check the emissivity at all, because that is a property of the surface and a mixture has no
//! surface. A half-copper board is 0.05 bare and 0.9 under solder mask.

use dualis::core::mixture::Mix;
use dualis::prelude::*;
use dualis_world::{Scene, World};

/// A scene whose block is made of `material`, with `materials` and `composites` blocks spliced in.
fn scene_text(materials: &str, composites: &str, material: &str) -> String {
    format!(
        r#"{{
        "title": "a declared composite",
        "duration_s": 1.0,
        "frames": 2,
        "materials": {materials},
        "composites": {composites},
        "domains": [
            {{ "kind": "block", "name": "slab", "cells": [2, 2, 8], "cell_mm": 2.0,
               "initial_c": 20.0, "material": "{material}" }}
        ]
    }}"#
    )
}

/// n-Octadecane, so a composite can be made partly of something the catalogue has never heard of.
const WAX: &str = r#"{
    "octadecane": {
        "name": "n-octadecane",
        "density": 814.0,
        "thermal": { "conductivity": 0.358, "specific_heat": 1934.0,
                     "expansion": 8.0e-4, "emissivity": 0.9 },
        "fusion": { "melting_point": 301.3, "latent_heat": 244000.0 }
    }
}"#;

/// Wax in an aluminium matrix, four fifths wax by volume. Bounds for that pair: Reuss 0.4473, Voigt
/// 33.6864 W/m·K.
fn buffer(conductivity: f64) -> String {
    format!(
        r#"{{
        "buffer": {{
            "parts": [
                {{ "material": "octadecane", "volume_fraction": 0.8 }},
                {{ "material": "aluminium",  "volume_fraction": 0.2 }}
            ],
            "conductivity_w_per_m_k": {conductivity},
            "emissivity": 0.9
        }}
    }}"#
    )
}

fn build(materials: &str, composites: &str, material: &str) -> Result<World, String> {
    let scene: Scene = serde_json::from_str(&scene_text(materials, composites, material))
        .map_err(|e| format!("parse: {e}"))?;
    World::build(scene)
}

/// The refusal, or a panic naming what was accepted.
fn refusal(materials: &str, composites: &str, material: &str) -> String {
    match build(materials, composites, material) {
        Err(e) => e,
        Ok(_) => panic!("built, and should not have"),
    }
}

/// **A composite of a catalogue material and a declared one reaches the domain with the mixture's own
/// numbers.**
///
/// The whole point, and building is not the assertion — a composite that parsed and was then ignored
/// would build perfectly well and give a block of aluminium. What pins it is the stability limit, which
/// is `ρ c dx²/(n k)` and therefore a function of the mixture's density, its specific heat and the chosen
/// conductivity together.
///
/// The density and the specific heat are the *exact* mixture rules — volume-additive and mass-weighted —
/// so this also checks that the format did not confuse them, which for this pair is worth a factor of
/// two: wax at 814 kg/m³ filling 80% of the volume is 54.7% of the mass.
#[test]
fn a_declared_composite_reaches_the_domain_with_the_mixtures_numbers() {
    let world = build(WAX, &buffer(5.0), "buffer").expect("it builds");

    // The mixture, computed here from the two constituents rather than read back from the scene.
    let (rho_w, c_w, f_w) = (814.0, 1934.0, 0.8);
    let (rho_a, c_a, f_a) = (2700.0, 896.0, 0.2);
    let rho = f_w * rho_w + f_a * rho_a;
    let volumetric = f_w * rho_w * c_w + f_a * rho_a * c_a;
    let alpha = 5.0 / volumetric;
    let dx = 2e-3;
    // A 2×2×8 block: one face in x, one in y, at most two in z. Four, not six — this domain sums the
    // actual face conductances, so the limit is a property of the shape.
    let expected = dx * dx / (4.0 * alpha);

    let block = world
        .simulation()
        .domain_as::<dualis::thermal::Solid3D>("slab")
        .expect("the block is there");
    let got = block
        .max_stable_dt(dualis::units::Time::from_si(0.0))
        .to_si();
    println!(
        "  mixture: rho {rho:.1} kg/m3, rho·c {volumetric:.0} J/m3/K, alpha {alpha:.3e} m2/s\n  \
         limit {got:e} s against {expected:e}"
    );
    assert!(
        (got / expected - 1.0).abs() < 1e-12,
        "the block's limit is {got:e} s and the mixture's own numbers give {expected:e}"
    );

    // And it is not either constituent, stated as the comparison the failure would produce.
    for (name, r, c) in [("wax", rho_w, c_w), ("aluminium", rho_a, c_a)] {
        let theirs = dx * dx / (4.0 * (5.0 / (r * c)));
        assert!(
            (got / theirs - 1.0).abs() > 0.1,
            "the limit matches {name}'s {theirs:e} rather than the mixture's"
        );
    }
}

/// **A conductivity no microstructure could have is refused, and the message carries all four bounds.**
///
/// The valuable half. A scene file cannot compute an effective conductivity and has nowhere to look one
/// up, so being told "no microstructure of these parts conducts that, and here is the range that is
/// achievable and the narrower range if it is isotropic" is the difference between a refusal and a useful
/// refusal.
///
/// Both ends, because a value below the Reuss bound and one above the Voigt bound are different mistakes
/// — too conservative and too optimistic — and a check on one is not a check on the other.
#[test]
fn an_impossible_conductivity_is_refused_with_the_bounds_in_the_message() {
    // Reuss 0.4473, Voigt 33.6864 for wax 0.8 / aluminium 0.2.
    for (chosen, what) in [(0.4, "below Reuss"), (40.0, "above Voigt")] {
        let error = refusal(WAX, &buffer(chosen), "buffer");
        println!(
            "  {chosen} W/m/K ({what}): {}",
            &error[..error.len().min(150)]
        );
        assert!(
            error.contains("composites.buffer") && error.contains("no microstructure"),
            "the message should name the composite and say what is wrong: {error}"
        );
        for bound in ["0.44", "33.68"] {
            assert!(
                error.contains(bound),
                "the message should carry the achievable range, is missing {bound}: {error}"
            );
        }
        assert!(
            error.contains("Hashin") && error.contains("narrower"),
            "and the tighter pair, which is where a caller should look: {error}"
        );
    }

    // Both bounds themselves are allowed, because both are attained by a real laminate. The endpoints
    // come from `Mix` rather than being typed in: what is under test here is the boundary condition of
    // the acceptance check, and a hand-typed 33.68644 against a true 33.6864 fails it for being 1.2e-6
    // too high — which is the assertion catching the test rather than the code, as it did on the first
    // run.
    let mix = Mix::of(&[
        (
            serde_json::from_str::<Substance>(
                &serde_json::to_string(
                    &serde_json::from_str::<serde_json::Value>(WAX).expect("wax parses")
                        ["octadecane"],
                )
                .expect("re-serialises"),
            )
            .expect("the wax is a substance"),
            0.8,
        ),
        (Substance::aluminium_6061(), 0.2),
    ])
    .expect("fractions sum to one");
    let (reuss, voigt) = mix.conductivity_bounds().expect("both conduct");
    for edge in [reuss.to_si(), voigt.to_si()] {
        assert!(
            build(WAX, &buffer(edge), "buffer").is_ok(),
            "{edge} W/m/K is a bound and bounds are reachable"
        );
    }
}

/// **A material used only inside a composite counts as used.**
///
/// The interaction a first draft got wrong. `Palette::unused` refuses a declared substance nothing names,
/// because a block with no `material` silently runs as aluminium — but a wax declared purely to be mixed
/// *is* named, by the composite, and refusing it would be the rule firing on the case it exists to
/// protect.
///
/// The composite itself is still subject to the rule, and that is the other half of this test: a composite
/// nothing names is dead weight in exactly the way a substance is.
#[test]
fn a_material_used_only_by_a_composite_is_not_dead_weight() {
    // The wax is named by the composite and by nothing else. The block asks for the composite.
    assert!(
        build(WAX, &buffer(5.0), "buffer").is_ok(),
        "a material mixed into a composite is used"
    );

    // A composite nothing names is refused, as a substance would be.
    let error = refusal(WAX, &buffer(5.0), "aluminium");
    println!("  unused composite: {}", &error[..error.len().min(140)]);
    assert!(
        error.contains("buffer") && error.contains("used by nothing"),
        "an unnamed composite is dead weight: {error}"
    );

    // And a material named by neither is still refused, so the rule did not simply stop applying.
    let extra = WAX.replace(
        "\"octadecane\": {",
        "\"spare\": { \"name\": \"spare\", \"density\": 1000.0 },\n    \"octadecane\": {",
    );
    let error = refusal(&extra, &buffer(5.0), "buffer");
    assert!(
        error.contains("spare"),
        "a material named by nothing at all is still dead weight: {error}"
    );
}

/// **A composite may not be made of another composite, and the refusal says to flatten it.**
///
/// Not an ordering problem, though it is that too — a `BTreeMap` has no declaration order and resolving
/// nested composites would need a dependency walk. The reason worth stating is that nesting changes what
/// the bounds *mean*: Hashin–Shtrikman is a two-phase result and a mixture of three things has no HS pair
/// at all, which `a_mixture.rs` records rather than papers over. Flattening asks the caller for the thing
/// the bounds are actually about.
#[test]
fn a_composite_of_composites_is_refused_with_the_remedy() {
    let nested = r#"{
        "inner": {
            "parts": [
                { "material": "octadecane", "volume_fraction": 0.5 },
                { "material": "aluminium",  "volume_fraction": 0.5 }
            ],
            "conductivity_w_per_m_k": 5.0,
            "emissivity": 0.9
        },
        "outer": {
            "parts": [
                { "material": "inner",  "volume_fraction": 0.5 },
                { "material": "copper", "volume_fraction": 0.5 }
            ],
            "conductivity_w_per_m_k": 50.0,
            "emissivity": 0.9
        }
    }"#;
    let error = refusal(WAX, nested, "outer");
    println!("  nested: {}", &error[..error.len().min(160)]);
    assert!(
        error.contains("flatten"),
        "the refusal should say what to do instead: {error}"
    );
}

/// **Fractions that do not sum to one, a name that collides, and an impossible emissivity are all
/// refused.**
///
/// The mechanical refusals, each one a mistake somebody makes in a file with no completion behind it. The
/// fractions one is the important one and is `Mix`'s own rule: they are **not** normalised, because 45%
/// and 50% is a transcription mistake rather than a request for 47.4% and 52.6%.
#[test]
fn the_mechanical_mistakes_in_a_composite_are_refused() {
    let bad_fractions = buffer(5.0).replace("0.8 }", "0.7 }");
    let error = refusal(WAX, &bad_fractions, "buffer");
    assert!(
        error.contains("composites.buffer") && error.contains("sum to 1"),
        "fractions are not normalised for the caller: {error}"
    );

    // A name that is already a catalogue material.
    let shadow = buffer(5.0).replace("\"buffer\"", "\"copper\"");
    let error = refusal(WAX, &shadow, "copper");
    assert!(
        error.contains("catalogue"),
        "a composite may not redefine a catalogue name: {error}"
    );

    // A name that is also a declared material — one name, two things.
    let collide = buffer(5.0).replace("\"buffer\"", "\"octadecane\"");
    let error = refusal(WAX, &collide, "octadecane");
    assert!(
        error.contains("two things"),
        "one name cannot be both a material and a composite: {error}"
    );

    // An emissivity outside 0..=1.
    let bad_emissivity = buffer(5.0).replace("\"emissivity\": 0.9", "\"emissivity\": 1.4");
    assert!(
        refusal(WAX, &bad_emissivity, "buffer").contains("emissivity"),
        "an emissivity above one is not a fraction of a blackbody's"
    );

    // A part naming nothing at all.
    let missing = buffer(5.0).replace("\"octadecane\"", "\"octadecain\"");
    let error = refusal(WAX, &missing, "buffer");
    assert!(
        error.contains("composites.buffer.parts[0]") && error.contains("unknown material"),
        "the site should name the composite and the index: {error}"
    );
}

/// **A scene with no `composites` key is unchanged, and one survives the round trip.**
///
/// The format promise is that a file which loads today loads tomorrow, and this key is new. Twenty-one
/// scene files do not have it.
///
/// The round trip is asserted from the *parsed* value rather than by comparing two serialisations, because
/// comparing `to_string(parse(x))` with itself puts serialiser output on both sides and cannot see a field
/// the parser dropped. This crate learned that one the hard way.
#[test]
fn the_key_is_optional_and_a_composite_round_trips() {
    let text = scene_text("{}", "{}", "copper").replace("\"composites\": {},", "");
    let scene: Scene = serde_json::from_str(&text).expect("parses without the key");
    assert!(scene.composites.is_empty());
    assert!(World::build(scene).is_ok(), "and builds");

    let text = scene_text(WAX, &buffer(5.0), "buffer");
    let scene: Scene = serde_json::from_str(&text).expect("parses");
    let written = serde_json::to_string(&scene).expect("serialises");
    let back: Scene = serde_json::from_str(&written).expect("and parses again");
    let spec = back.composites.get("buffer").expect("still declared");
    assert_eq!(spec.parts.len(), 2);
    assert_eq!(spec.parts[0].material, "octadecane");
    assert_eq!(spec.parts[0].volume_fraction, 0.8);
    assert_eq!(spec.conductivity_w_per_m_k, 5.0);
    assert!(World::build(back).is_ok(), "and the reparsed scene builds");

    // An empty map is not written out, so twenty-one existing files re-serialise to what they were.
    let plain: Scene = serde_json::from_str(&scene_text("{}", "{}", "copper")).expect("parses");
    assert!(
        !serde_json::to_string(&plain)
            .expect("serialises")
            .contains("composites"),
        "an empty map should not be written"
    );
}
