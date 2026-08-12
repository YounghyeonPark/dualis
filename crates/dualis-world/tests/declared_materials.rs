//! A scene that brings its own materials.
//!
//! The catalogue holds nine substances and the world holds hundreds of thousands, so a format that can
//! only *name* a material can only ever describe nine kinds of thing. `Substance` has been
//! `Deserialize` all along; what was missing was a place in the file to put one.
//!
//! The physics of a declared substance is checked in `dualis`'s `substances_from_a_file.rs`, against
//! Neumann's exact solution and beside ice through the same harness. This file is about the other half,
//! which is the half a format gets wrong: **what happens when the declaration is a mistake.**
//!
//! Three of them are refused, and the third is the one worth reading. `material` on a block is optional
//! and defaults to aluminium, so a scene that declares a substance and then does not name it — a field
//! left off, a name misspelled at the use site where the declaration is spelled right — runs as a block
//! of aluminium. It runs, it audits, it renders, and it answers a question about the wrong material
//! with nothing anywhere saying so. That is the same failure a region selecting no cells has, one level
//! up, and it gets the same treatment.

use dualis::prelude::*;
use dualis_world::{Scene, World};

/// A block of a declared substance, as a scene file. `{extra}` is where each test puts its mistake.
fn scene_text(materials: &str, block_material: &str) -> String {
    format!(
        r#"{{
        "title": "a declared substance",
        "duration_s": 1.0,
        "frames": 2,
        "materials": {materials},
        "domains": [
            {{ "kind": "block", "name": "slab", "cells": [2, 2, 8], "cell_mm": 2.0,
               "initial_c": 20.0{block_material} }}
        ]
    }}"#
    )
}

/// n-Octadecane, the wax a phase-change thermal buffer is made of, as a scene would write it.
const WAX: &str = r#"{
    "octadecane": {
        "name": "n-octadecane",
        "density": 814.0,
        "thermal": { "conductivity": 0.358, "specific_heat": 1934.0,
                     "expansion": 8.0e-4, "emissivity": 0.9 },
        "fusion": { "melting_point": 301.3, "latent_heat": 244000.0 }
    }
}"#;

fn build(materials: &str, block_material: &str) -> Result<World, String> {
    let scene: Scene = serde_json::from_str(&scene_text(materials, block_material))
        .map_err(|e| format!("parse: {e}"))?;
    World::build(scene)
}

/// The refusal, or a panic naming what was accepted. `expect_err` needs `Debug` on the success type and
/// a `World` has no business having one.
fn refusal(materials: &str, block_material: &str) -> String {
    match build(materials, block_material) {
        Err(e) => e,
        Ok(_) => panic!("built, and should not have: {block_material}"),
    }
}

/// **A substance the library has never heard of can be written in a scene and used.**
///
/// The whole point, and the assertion is not that it builds — it is that the *properties reached the
/// domain*. A declaration that parsed and was then ignored would build perfectly well and give a block
/// of aluminium, so building proves nothing on its own.
///
/// What pins it is the stability limit, which is `ρ c dx²/(n k)` and therefore a function of all three
/// declared numbers with no pair of mistakes in them cancelling. A temperature or a cell count could not
/// say that.
///
/// `n = 4`, and getting that wrong is how this test first failed. A 2×2×8 block has more than one cell
/// along all three axes, so the naive answer is the three-dimensional `6α` — but this domain sums the
/// **actual face conductances**, and with only two cells along x and y every cell has one face in each of
/// those and at most two in z. Four faces, not six. The limit is a property of the shape, which is the
/// whole point of that change, and a test asserting `6α` here would have been asserting the bug it
/// replaced.
#[test]
fn a_declared_substance_reaches_the_domain_with_its_own_numbers() {
    let world = build(WAX, r#", "material": "octadecane""#).expect("it builds");

    let alpha = 0.358 / (814.0 * 1934.0);
    let dx = 2e-3;
    let expected = dx * dx / (4.0 * alpha);

    let block = world
        .simulation()
        .domain_as::<dualis::thermal::Solid3D>("slab")
        .expect("the block is there");
    let got = block
        .max_stable_dt(dualis::units::Time::from_si(0.0))
        .to_si();
    assert!(
        (got / expected - 1.0).abs() < 1e-12,
        "the block's limit is {got:e} s, and the declared wax's own numbers give {expected:e}"
    );

    // And the substance is not aluminium, stated as the comparison rather than left implied — the
    // failure being excluded is precisely "the declaration was ignored and this is the default".
    // 167 W/m/K over 2700 kg/m³ and 896 J/kg/K, on the same four faces: a limit 303× shorter.
    let aluminium = dx * dx / (4.0 * 167.0 / (2700.0 * 896.0));
    assert!(
        got > 100.0 * aluminium,
        "the limit {got:e} s is near aluminium's {aluminium:e}, so the declaration did not take"
    );
}

/// **An impossible declared substance is refused, and the message names the material and the field.**
///
/// `Substance::check` is what does the work and is tested on its own in `any_material.rs`; what is
/// tested here is that the scene format *calls* it, before any domain is built. A negative conductivity
/// otherwise reaches the sweep as a `NaN` — which this domain does refuse, but by then the message can
/// only say that a step was unstable, not that `materials.octadecane` has a minus sign in it.
#[test]
fn an_impossible_declared_substance_is_refused_at_the_declaration() {
    let bad = WAX.replace("0.358", "-0.358");
    let error = refusal(&bad, r#", "material": "octadecane""#);
    assert!(
        error.contains("materials.octadecane") && error.contains("conductivity"),
        "the message has to name both the declaration and the field: {error}"
    );

    // A property key this format does not know is refused rather than dropped, which is `Substance`'s
    // own `deny_unknown_fields` doing its job through the scene. Dropped, `thermalz` would leave a
    // substance whose conductivity is *unknown* — and unknown is not the same failure as mistyped.
    let typo = WAX.replace("\"thermal\":", "\"thermalz\":");
    assert!(
        build(&typo, r#", "material": "octadecane""#).is_err(),
        "a mistyped property key must not be silently dropped"
    );
}

/// **A declaration that shadows a catalogue name is refused.**
///
/// Two files that both say `"copper"` have to mean the same copper. If one of them can redefine it,
/// then a run's material is a property of the file it was launched from, and no comparison between two
/// runs means anything — including the comparisons this workspace's own tests make.
///
/// Refused rather than resolved by precedence in either direction. Catalogue-wins makes the
/// declaration silently dead, which is the failure mode this whole file is about; scene-wins is the
/// ambiguity above. There is no third behaviour that is not one of those two, so the answer is to
/// refuse.
#[test]
fn a_declaration_may_not_redefine_a_catalogue_material() {
    let shadow = WAX.replace("octadecane", "copper");
    let error = refusal(&shadow, r#", "material": "copper""#);
    assert!(
        error.contains("copper") && error.contains("catalogue"),
        "the message should say what the collision is: {error}"
    );

    // Every one of the nine, not just the one that happened to be tried. A partial guard here would be
    // a guard that depends on which name somebody reached for first.
    for name in dualis_world::MATERIALS {
        let shadow = WAX.replace("octadecane", name);
        assert!(
            build(&shadow, &format!(", \"material\": {name:?}")).is_err(),
            "{name:?} could be redefined"
        );
    }
}

/// **A declared material that nothing uses is refused, because the alternative is a run on aluminium.**
///
/// The one that is not obvious. Nothing is wrong with the declaration and nothing is wrong with the
/// domain; the mistake is the *absence* of a connection between them, and absence is the shape of
/// failure this workspace keeps finding.
///
/// Both spellings of it are here. Leaving `material` off the block is the one that costs a wrong
/// answer — aluminium conducts 466× better than this wax, so the block reaches steady state in a
/// millisecond where the wax would take a second, and the picture looks like a working simulation
/// throughout. Misspelling it at the use site is caught by the resolver anyway, and is here so that
/// the two failures are known to produce different messages: one says the name is unknown, the other
/// says the declaration is unused.
#[test]
fn a_declaration_nothing_uses_is_refused() {
    // The block says nothing, so it is aluminium and the wax is dead weight.
    let error = refusal(WAX, "");
    assert!(
        error.contains("octadecane") && error.contains("aluminium"),
        "the message has to say what the scene would have run as: {error}"
    );

    // The use site is misspelled. Caught by the resolver, with a different message, and the list of
    // what would have worked includes the declared name — which is the only way a caller finds out
    // that their declaration is fine and their *reference* is not.
    let error = refusal(WAX, r#", "material": "octadecain""#);
    assert!(
        error.contains("octadecain") && error.contains("(declared)"),
        "an unknown name should list the declared ones: {error}"
    );

    // And the ordinary case still builds, so this is not a test that everything is refused.
    assert!(build(WAX, r#", "material": "octadecane""#).is_ok());
}

/// **A scene with no `materials` key is unchanged, and every catalogue name still resolves.**
///
/// The format promise is that a file which loads today loads tomorrow, and this key is new. Nineteen
/// scene files in this repository do not have it.
///
/// `water` is in the loop, and it is the reason this whole thing is one list now: it was in the
/// catalogue and not in the format's private copy of the catalogue, so no scene could name it for
/// eleven releases. Nothing could have noticed — an absent name is not a wrong answer, it is a
/// substance that never appears.
#[test]
fn the_catalogue_still_resolves_and_the_key_is_optional() {
    assert!(
        build("{}", r#", "material": "copper""#).is_ok(),
        "an empty declaration block is not a declaration nothing uses"
    );
    // No `materials` key at all, which is what every scene written before this looks like.
    let text = scene_text("{}", r#", "material": "copper""#).replace("\"materials\": {},", "");
    let scene: Scene = serde_json::from_str(&text).expect("parses without the key");
    assert!(World::build(scene).is_ok(), "and builds");

    for name in dualis_world::MATERIALS {
        assert!(
            build("{}", &format!(", \"material\": {name:?}")).is_ok(),
            "{name:?} is in MATERIALS and a scene cannot use it"
        );
    }
    assert_eq!(
        dualis_world::MATERIALS.len(),
        9,
        "and there are nine of them, `water` included"
    );
}

/// **A declared substance survives the round trip through the file the run writes.**
///
/// A scene is re-serialised into the output, and a declaration that did not survive that would produce
/// a file describing a run nobody can reproduce — with the material silently back to whatever the name
/// resolves to elsewhere, which is nothing.
///
/// Asserted from the *parsed* value rather than by comparing two serialisations, because comparing
/// `to_string(parse(x))` with itself puts serialiser output on both sides and cannot see a field the
/// parser dropped. This crate learned that one the hard way.
#[test]
fn a_declaration_survives_being_written_back_out() {
    let text = scene_text(WAX, r#", "material": "octadecane""#);
    let scene: Scene = serde_json::from_str(&text).expect("parses");
    let written = serde_json::to_string(&scene).expect("serialises");
    let back: Scene = serde_json::from_str(&written).expect("and parses again");

    let original = scene.materials.get("octadecane").expect("declared");
    let survived = back.materials.get("octadecane").expect("still declared");
    assert_eq!(original, survived, "the declaration round trips exactly");
    assert!(
        World::build(back).is_ok(),
        "and the reparsed scene still builds"
    );

    // A scene with no declarations does not grow an empty key, so nineteen existing files re-serialise
    // to what they were.
    let plain: Scene =
        serde_json::from_str(&scene_text("{}", r#", "material": "copper""#)).expect("parses");
    assert!(
        !serde_json::to_string(&plain)
            .expect("serialises")
            .contains("materials"),
        "an empty map should not be written"
    );
}
