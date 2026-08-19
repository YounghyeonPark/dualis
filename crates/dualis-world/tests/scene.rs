//! The app holds itself to the library's standard: every number checked against something
//! the code did not compute.

use dualis::electrical::Winding;
use dualis::prelude::*;
use dualis_world::{DomainSpec, Scene, World};

fn scene_from(text: &str) -> Scene {
    serde_json::from_str(text).expect("it parsed once already")
}

fn room_scene(duration_s: f64, frames: usize, cells: usize) -> Scene {
    serde_json::from_str(&format!(
        r#"{{
          "title": "test", "schedule": "multirate",
          "duration_s": {duration_s}, "frames": {frames},
          "domains": [{{ "kind": "room", "name": "room", "width_m": 4.4,
            "height_m": 3.1, "cells_across": {cells},
            "release": {{ "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 }} }}]
        }}"#
    ))
    .expect("the test scene parses")
}

/// A scene survives being written down and read back.
///
/// The point of a scene format is that it leaves the process. If a round trip loses a field,
/// every saved world is quietly a different world.
/// **Asserted from the hand-written text, not from the serialiser's own output.**
///
/// The previous version compared `to_string(parse(hand_written))` against
/// `to_string(parse(to_string(parse(hand_written))))` — both sides serialiser output, so the
/// hand-written spelling never entered any comparison. A key the parser silently dropped was
/// already gone before the first `to_string`, and the assertion could not see it. Regressing
/// the helper to the pre-`release` spelling left all ten tests passing.
///
/// This goes the other way: parse the text a person would type, serialise it, and require
/// specific substrings to be present. A dropped key then shows as a missing key in the bytes.
/// `deny_unknown_fields` on the scene types is the real guard; this is the test that would
/// have noticed before it existed.
#[test]
fn a_hand_written_scene_survives_being_parsed() {
    let text = r#"{
      "title": "round trip", "schedule": "staggered",
      "duration_s": 0.25, "frames": 4, "conservation_tolerance": 1e-7,
      "domains": [
        { "kind": "room", "name": "room", "width_m": 4.4, "height_m": 3.1,
          "cells_across": 41,
          "release": { "as": "pulse", "x_m": 1.0, "y_m": 0.8,
                       "radius_m": 0.2, "amplitude_pa": 3.0 } },
        { "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 21,
          "area_mm2": 100.0, "initial_c": 20.0,
          "exposes": { "name": "face", "face_area_mm2": 0.5 } }
      ]
    }"#;
    let scene: Scene = serde_json::from_str(text).expect("the hand-written scene parses");
    let out = serde_json::to_string(&scene).unwrap();

    // Every value the text states must come back out. A `#[serde(default)]` field that was
    // dropped and refilled would fail here on the value, not merely on the key.
    for want in [
        r#""title":"round trip""#,
        r#""schedule":"staggered""#,
        r#""conservation_tolerance":1e-7"#,
        r#""as":"pulse""#,
        r#""amplitude_pa":3.0"#,
        r#""radius_m":0.2"#,
        r#""exposes":{"name":"face""#,
    ] {
        assert!(
            out.contains(want),
            "{want} did not survive the round trip:
{out}"
        );
    }
    assert_eq!(scene.domains.len(), 2);
    assert!((scene.duration_s - 0.25).abs() < 1e-15);

    // And a second pass is byte-stable, which is what catches a field that serialises to
    // something it cannot read back.
    let again: Scene = serde_json::from_str(&out).unwrap();
    assert_eq!(serde_json::to_string(&again).unwrap(), out);
}

/// A key this format does not know is refused, not discarded.
///
/// The pre-`release` spelling, which used to parse into the default and produce a
/// byte-identical run — so editing it was a no-op that reported success.
#[test]
fn an_unknown_key_is_refused_rather_than_dropped() {
    let stale = r#"{
      "title": "the old spelling", "duration_s": 0.01, "frames": 2,
      "domains": [
        { "kind": "room", "name": "room", "width_m": 4.4, "height_m": 3.1,
          "cells_across": 41, "mode": [3, 2], "amplitude_pa": 7.0 }
      ]
    }"#;
    let err = serde_json::from_str::<Scene>(stale).expect_err("`mode` is not a field any more");
    let message = err.to_string();
    assert!(
        message.contains("mode") && message.contains("release"),
        "the refusal should name the key and what is expected: {message}"
    );
}

/// The standing mode oscillates at the frequency the closed form gives, and the gap
/// between the two converges at **second** order.
///
/// The room is released at its (1,1) antinode with amplitude 1 Pa. A standing mode is
/// separable, so every point follows `cos(2 pi f t)` and the peak of the field is
/// `|cos(2 pi f t)|`. That is a closed form the integration never sees.
///
/// Worst departure over a run of 0.02 s, against grid resolution, before and after the
/// leapfrog startup was fixed:
///
/// ```text
///   cells     first order      second order
///      31        0.0528            0.00238
///      61        0.0265            0.00059
///     121        0.0151            0.00007
///     241        0.0076            0.00002
///     481        0.0039            0.0000032
/// ```
///
/// The left column halves on refinement and the right one quarters. The cause of the left
/// column was `Room::released_from` leaving the velocity at `t = 0` when a staggered
/// leapfrog carries it at `t = -h/2`; the first velocity update then travelled a whole step
/// where it was owed half. `O(h)`, permanent, and enough to drag a second-order scheme to
/// first. `Tube` had it too.
///
/// **The rate is asserted and not the size**, because only the rate separates a coarse
/// scheme from a wrong one — the same lesson as the wall-weighting defect, which was 1.4%
/// and looked like coarseness. Measured across three doublings rather than one: the
/// per-doubling ratio bounces between 3.9 and 8.1 because "worst over forty sampled frames"
/// is a maximum and therefore noisy, while the span 31 -> 241 is a stable 127x. Second order
/// over three doublings is 64x and first order is 8x, so 40x separates them with room on
/// both sides.
#[test]
fn the_room_rings_at_the_closed_form_frequency_and_converges_at_second_order() {
    let worst_at = |cells: usize| {
        let probe = Room::of_air("probe", Length::m(4.4), Length::m(3.1), cells);
        let f = probe.mode_frequency(1, 1).to_si();
        let mut world = World::build(room_scene(0.02, 40, cells)).expect("the scene builds");
        world
            .run()
            .expect("a rigid room conserves")
            .iter()
            .map(|frame| {
                let peak = frame.panels[0]
                    .values()
                    .iter()
                    .fold(0.0f64, |m, v| m.max(v.abs()));
                let want = (2.0 * std::f64::consts::PI * f * frame.time_s).cos().abs();
                (peak - want).abs()
            })
            .fold(0.0f64, f64::max)
    };

    // The size first: a 61-cell grid tracks the closed form to a tenth of a percent of the
    // amplitude. Measured 0.00059, so this uses 30% of its budget.
    let mid = worst_at(61);
    assert!(mid < 0.002, "61 cells departed by {mid:.5} Pa");

    let (coarse, fine) = (worst_at(31), worst_at(241));
    let fall = coarse / fine;
    assert!(
        fall > 40.0,
        "31 -> 241 cells is three doublings: second order is 64x, first order 8x.          Got {coarse:.5} -> {fine:.5}, a factor of {fall:.1}"
    );
}

/// A scene that cannot describe a simulation is refused before anything runs.
#[test]
fn a_scene_that_makes_no_sense_is_refused() {
    let mut empty = room_scene(0.01, 4, 41);
    empty.domains.clear();
    assert!(World::build(empty).is_err());

    let mut backwards = room_scene(0.01, 4, 41);
    backwards.duration_s = -1.0;
    assert!(World::build(backwards).is_err());

    let mut still = room_scene(0.01, 4, 41);
    still.frames = 0;
    assert!(World::build(still).is_err());
}

/// **A region that selects no cells is refused, and a block with no material is still aluminium.**
///
/// Two halves of the same worry, in opposite directions.
///
/// A mistyped region is the silent failure this key can have. `"to": [9, 9, 9]` where `[9, 9, 18]`
/// was meant selects nothing, and a block of one material runs, audits, renders and answers the
/// wrong question with nothing anywhere saying the coating was not applied. So an empty region — or
/// one naming a cell the block does not have — is a refusal with the convention spelled out, because
/// half-open bounds are exactly what a person gets wrong here.
///
/// And the other way: `material` and `regions` were **added** to a format that carries a
/// compatibility promise, so every block written before they existed has to keep meaning what it
/// meant. Absence is aluminium, and that is asserted from the old spelling rather than assumed.
#[test]
fn a_region_that_selects_nothing_is_refused_and_absence_is_aluminium() {
    let block = |extra: &str| {
        format!(
            r#"{{ "title": "t", "schedule": "multirate", "duration_s": 1e-4, "frames": 2,
              "domains": [{{ "kind": "block", "name": "b", "cells": [4, 4, 8],
                "cell_mm": 1.0, "initial_c": 20.0{extra} }}] }}"#
        )
    };

    // The pre-`material` spelling, which is what every shipped block scene was.
    let old: Scene = serde_json::from_str(&block("")).expect("an old block still parses");
    let world = World::build(old).expect("and still builds");
    let b = world
        .simulation()
        .domain_as::<dualis::thermal::Solid3D>("b")
        .expect("it is a block");
    assert_eq!(b.substances(), 1, "no regions means one material");
    assert_eq!(
        b.substance_at(0, 0, 0).name,
        "Al 6061",
        "absence of `material` is aluminium, as it was before the key existed"
    );

    // A region that is fine, to establish that the refusals below are about the bounds.
    let good = block(
        r#", "material": "copper",
             "regions": [{ "material": "fr4", "from": [0, 0, 4], "to": [4, 4, 8] }]"#,
    );
    let world = World::build(serde_json::from_str(&good).expect("parses")).expect("builds");
    let b = world
        .simulation()
        .domain_as::<dualis::thermal::Solid3D>("b")
        .expect("it is a block");
    assert_eq!(b.substances(), 2, "the region is a second material");
    assert_eq!(b.substance_at(0, 0, 0).name, "Cu ETP");
    assert_eq!(b.substance_at(0, 0, 7).name, "FR-4");

    for (why, spec) in [
        (
            "an empty range",
            r#", "regions": [{ "material": "fr4", "from": [0, 0, 4], "to": [4, 4, 4] }]"#,
        ),
        (
            "a range that runs backwards",
            r#", "regions": [{ "material": "fr4", "from": [0, 0, 6], "to": [4, 4, 2] }]"#,
        ),
        (
            "a cell the block does not have",
            r#", "regions": [{ "material": "fr4", "from": [0, 0, 4], "to": [4, 4, 9] }]"#,
        ),
        (
            "an unknown material",
            r#", "regions": [{ "material": "unobtainium", "from": [0, 0, 4], "to": [4, 4, 8] }]"#,
        ),
    ] {
        let scene: Scene =
            serde_json::from_str(&block(spec)).expect("it parses; the check is later");
        let Err(err) = World::build(scene) else {
            panic!("{why} must be refused");
        };
        assert!(
            err.contains("b/regions[0]"),
            "{why}: the message should say which region: {err}"
        );
    }
}

/// **A region can be made of nothing, and a scene cannot redefine what nothing means.**
///
/// The format could say what a box was made of and could not say it was *empty*, so a clearance
/// between two parts had to be spelled as some substance with a low conductivity — which still
/// conducts, still stores heat and still sets a stability limit. `"void"` is the one reserved
/// spelling here, and it is reserved rather than resolved: a scene that declared a material by
/// that name would be solid in one file and empty in another.
#[test]
fn a_region_can_be_nothing_and_nothing_is_a_reserved_word() {
    let scene = |extra: &str| {
        format!(
            r#"{{ "title": "t", "schedule": "multirate", "duration_s": 1e-4, "frames": 2{extra},
              "domains": [{{ "kind": "block", "name": "b", "cells": [4, 4, 4],
                "cell_mm": 5.0, "initial_c": 20.0,
                "regions": [{{ "material": "void", "from": [0, 0, 1], "to": [4, 4, 3] }}] }}] }}"#
        )
    };

    let world = World::build(serde_json::from_str(&scene("")).expect("parses")).expect("builds");
    let b = world
        .simulation()
        .domain_as::<dualis::thermal::Solid3D>("b")
        .expect("it is a block");
    assert_eq!(b.void_cells(), 32, "two layers of sixteen cells are empty");
    assert!(b.is_void(0, 0, 1) && !b.is_void(0, 0, 0));
    assert!(
        b.temperature_at(0, 0, 1).to_si().is_nan(),
        "there is nothing there to have a temperature"
    );

    // And the name cannot be taken. The message has to say why, or somebody will simply rename
    // their material and wonder where the gap went.
    let declared = r#", "materials": { "void": { "name": "not nothing", "density": 1.0,
        "thermal": { "conductivity": 1.0, "specific_heat": 1.0,
                     "expansion": 0.0, "emissivity": 0.5 } } }"#;
    let Err(err) = World::build(serde_json::from_str(&scene(declared)).expect("parses")) else {
        panic!("a scene may not redefine nothing");
    };
    assert!(
        err.contains("materials.void") && err.contains("nothing"),
        "the refusal should name the key and say what the word is for: {err}"
    );
}

/// **A region can start at its own temperature**, which is the difference between a scene that
/// can pose a transient and one that can only pose a source.
///
/// Found by writing scene `23`: a hot part under a cold lid is the commonest thermal question
/// there is, and the format could not state it in either direction. `regions` said what a box was
/// made of and `initial_c` was a property of the whole block, while the bus — deliberately — has
/// no location at all, so heat arriving there spreads over everything that can hold it. Neither
/// could put 300 C in one corner.
///
/// Refused on a void region, because nothing has no temperature to start at, and that would
/// otherwise be a value quietly discarded.
#[test]
fn a_region_can_start_hot_and_nothing_cannot() {
    let block = |regions: &str| {
        format!(
            r#"{{ "title": "t", "schedule": "multirate", "duration_s": 1e-6, "frames": 2,
              "domains": [{{ "kind": "block", "name": "b", "cells": [4, 4, 4],
                "cell_mm": 5.0, "initial_c": 20.0, "regions": {regions} }}] }}"#
        )
    };

    let hot = block(
        r#"[{ "material": "copper", "from": [0, 0, 0], "to": [4, 4, 1], "initial_c": 300.0 }]"#,
    );
    let world = World::build(serde_json::from_str(&hot).expect("parses")).expect("builds");
    let b = world
        .simulation()
        .domain_as::<dualis::thermal::Solid3D>("b")
        .expect("it is a block");
    assert!(
        (b.temperature_at(0, 0, 0).to_si() - 573.15).abs() < 1e-9,
        "the region starts where it says: {:.4} K",
        b.temperature_at(0, 0, 0).to_si()
    );
    assert!(
        (b.temperature_at(0, 0, 3).to_si() - 293.15).abs() < 1e-9,
        "and the rest starts at the block's own figure: {:.4} K",
        b.temperature_at(0, 0, 3).to_si()
    );

    // Omitting it leaves the block's figure, which is what every scene written before the key
    // existed relies on.
    let cold = block(r#"[{ "material": "copper", "from": [0, 0, 0], "to": [4, 4, 1] }]"#);
    let world = World::build(serde_json::from_str(&cold).expect("parses")).expect("builds");
    let b = world
        .simulation()
        .domain_as::<dualis::thermal::Solid3D>("b")
        .expect("it is a block");
    assert!((b.temperature_at(0, 0, 0).to_si() - 293.15).abs() < 1e-9);

    let void = block(
        r#"[{ "material": "void", "from": [0, 0, 0], "to": [4, 4, 1], "initial_c": 300.0 }]"#,
    );
    let Err(err) = World::build(serde_json::from_str(&void).expect("parses")) else {
        panic!("nothing cannot start at a temperature");
    };
    assert!(
        err.contains("b/regions[0]"),
        "the refusal should say which region: {err}"
    );
}

/// **A scene can say where the heat is made, and is refused when it says nowhere.**
///
/// The gap this key closes is structural rather than cosmetic. Every other source in this format
/// hands watts to the bus, and the bus carries an amount and no location *by design* — heat
/// arriving there spreads to a uniform rise over everything that can hold it, which is the only
/// choice that adds no information and the wrong answer for every real thing that dissipates.
///
/// The refusals matter as much as the key. A box that selects no cell, or only void, is a scene
/// that states watts and generates none — and it would run, conserve, and answer a question about
/// a different object.
#[test]
fn a_scene_can_state_where_its_heat_is_made() {
    let block = |extra: &str| {
        format!(
            r#"{{ "title": "t", "schedule": "multirate", "duration_s": 1e-3, "frames": 2,
              "domains": [{{ "kind": "block", "name": "b", "cells": [4, 4, 4],
                "cell_mm": 5.0, "initial_c": 20.0{extra} }}] }}"#
        )
    };

    let good = block(r#", "dissipation": [{ "watts": 30.0, "from": [0, 0, 0], "to": [2, 2, 2] }]"#);
    let world = World::build(serde_json::from_str(&good).expect("parses")).expect("builds");
    let b = world
        .simulation()
        .domain_as::<dualis::thermal::Solid3D>("b")
        .expect("it is a block");
    assert!(
        (b.generated_power().to_si() - 30.0).abs() < 1e-12,
        "the watts stated are the watts held: {} W",
        b.generated_power().to_si()
    );

    // **The total is the box's, not the cell's.** Eight cells or one, the block generates 30 W —
    // which is what lets `verify` refine the grid without changing the physics.
    let one = block(r#", "dissipation": [{ "watts": 30.0, "from": [0, 0, 0], "to": [1, 1, 1] }]"#);
    let world = World::build(serde_json::from_str(&one).expect("parses")).expect("builds");
    let b = world
        .simulation()
        .domain_as::<dualis::thermal::Solid3D>("b")
        .expect("it is a block");
    assert!((b.generated_power().to_si() - 30.0).abs() < 1e-12);

    // Absence generates nothing, which every scene written before this key existed relies on.
    let none = World::build(serde_json::from_str(&block("")).expect("parses")).expect("builds");
    assert_eq!(
        none.simulation()
            .domain_as::<dualis::thermal::Solid3D>("b")
            .expect("it is a block")
            .generated_power()
            .to_si(),
        0.0
    );

    for (why, spec) in [
        (
            "an empty range",
            r#", "dissipation": [{ "watts": 30.0, "from": [0, 0, 0], "to": [0, 4, 4] }]"#,
        ),
        (
            "a cell the block does not have",
            r#", "dissipation": [{ "watts": 30.0, "from": [0, 0, 0], "to": [4, 4, 5] }]"#,
        ),
        (
            "watts that are not a number",
            r#", "dissipation": [{ "watts": null, "from": [0, 0, 0], "to": [4, 4, 4] }]"#,
        ),
    ] {
        let parsed: Result<Scene, _> = serde_json::from_str(&block(spec));
        let refused = match parsed {
            Err(_) => true, // `null` is not a number and never reaches the builder
            Ok(scene) => World::build(scene).is_err(),
        };
        assert!(refused, "{why} must be refused");
    }

    // A box entirely inside a clearance states watts that would be generated nowhere.
    let in_void = block(
        r#", "regions": [{ "material": "void", "from": [0, 0, 0], "to": [4, 4, 2] }],
             "dissipation": [{ "watts": 30.0, "from": [0, 0, 0], "to": [4, 4, 2] }]"#,
    );
    let Err(err) = World::build(serde_json::from_str(&in_void).expect("parses")) else {
        panic!("watts generated in nothing must be refused");
    };
    assert!(
        err.contains("b/dissipation[0]") && err.contains("nowhere"),
        "the refusal should say which box and why: {err}"
    );
}

/// Two domains that do not interact still run, and both are captured.
#[test]
fn a_scene_can_hold_more_than_one_domain() {
    let mut scene = room_scene(0.01, 3, 31);
    scene.domains.push(
        serde_json::from_str::<DomainSpec>(
            r#"{ "kind": "bar", "name": "bar", "length_mm": 20.0,
                 "cells": 21, "area_mm2": 100.0, "initial_c": 20.0 }"#,
        )
        .unwrap(),
    );
    let mut world = World::build(scene).unwrap();
    let frames = world
        .run()
        .expect("neither domain publishes, so nothing can go missing");
    assert_eq!(frames[0].panels.len(), 2);
    assert_eq!(frames[0].panels[1].name, "bar");
    // An isolated bar at a uniform temperature has nowhere to send heat, so it stays put.
    let bar = &frames[frames.len() - 1].panels[1];
    for v in bar.values() {
        assert!((v - 293.15).abs() < 1e-9, "bar drifted to {v} K");
    }
}

/// **Two domains that actually talk**, driven entirely from data.
///
/// Everything before this ran domains side by side: they stepped, they were drawn, and none
/// of them ever put anything on the bus. So the part of the library the whole architecture
/// exists for — domains meeting on `Exchange` and the kernel auditing the crossing — had only
/// ever been exercised by tests written inside the workspace, which is the condition that
/// produced the API frictions in the first place.
///
/// A heater pays joules onto the channel; a bar takes them and warms. Neither names the
/// other, and the scene names neither: it says "heater" and "bar" and the coupling is the
/// channel they happen to share.
///
/// The heater is defined in *this* crate, which is the other half of the point. A consumer
/// implementing `Domain` from outside needs the trait, `Exchange`, `Ledger`, `Kind`,
/// `Violation` and a channel constant, and all six come out of `dualis::prelude`. Nothing
/// private was required.
///
/// `multirate` and not `staggered`, and the difference is not cosmetic: half a second is
/// thirty-eight times the bar's diffusion limit, so a staggered scene diverges. See
/// `a_schedule_a_scene_cannot_survive_is_named_and_refused` — the kernel catches it, but
/// only when the step is taken, which is a thing a scene format could check earlier.
#[test]
fn a_heater_and_a_bar_meet_on_the_bus() {
    let scene: Scene = serde_json::from_str(
        r#"{
          "title": "coupled", "schedule": "multirate",
          "duration_s": 4.0, "frames": 8,
          "conservation_tolerance": 1e-9,
          "domains": [
            { "kind": "heater", "name": "element", "watts": 2.0, "reserve_j": 6.0 },
            { "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 21,
              "area_mm2": 100.0, "initial_c": 20.0 }
          ]
        }"#,
    )
    .expect("the coupled scene parses");

    let mut world = World::build(scene).expect("it builds");
    // The audit is the assertion. At 1e-9 nothing may appear or vanish between the tank and
    // the bar, and anything left unclaimed on the channel at the end of a step is refused.
    let frames = world.run().expect("the books close across the bus");

    // The heater has no field, so it contributes no panel: one domain, one drawing.
    assert_eq!(frames[0].panels.len(), 1);
    assert_eq!(frames[0].panels[0].name, "bar");

    // 2 W for 4 s is 8 J of capacity against a 6 J tank, so it runs dry and the bar receives
    // exactly the tank — which is the number to check, because it is set by the scene rather
    // than by the physics.
    let heater = world
        .simulation()
        .domain_as::<dualis_world::heater::Heater>("element")
        .expect("the heater is still there");
    assert!(
        heater.reserve().to_si() < 1e-12,
        "the tank should be empty, {} J left",
        heater.reserve().to_si()
    );
    let crossed = world.simulation().bus().total_consumed(quantity::ENERGY);
    assert!(
        (crossed - 6.0).abs() < 1e-9,
        "every joule in the tank should have crossed, got {crossed}"
    );

    // And the bar is warmer by the amount those joules buy. An independent number: 6 J into
    // 20 mm x 1 cm^2 of aluminium is 6 / (rho V c_p), and the bar is insulated so none of it
    // leaves. Computed here from the substance rather than read off the bar.
    //
    // **Read from the domain, not from the panel.** This assertion used to average the render
    // panel, which is `ScalarField` *sampled* at evenly spaced points including both ends — and
    // the bar's grid is cell-centred, so the end samples sit half a cell outside the outermost
    // centres and the average is not the state. It passed at 1e-6 only because the field was
    // nearly uniform by the end of the run; changing *when* the heat arrives changed the
    // profile enough to expose it at 4.1e-6. That is `FRICTION.md` finding 10, in a test
    // written after finding 10 was documented.
    let capacity = Substance::aluminium_6061()
        .heat_capacity(Volume::from_si(20e-3 * 100e-6))
        .expect("aluminium has a specific heat");
    let expected_rise = 6.0 / capacity.to_si();
    let bar = world
        .simulation()
        .domain_as::<Bar1D>("bar")
        .expect("the bar is still there");
    let mean_rise = bar.mean_temperature().to_si() - Temperature::celsius(20.0).to_si();
    assert!(
        (mean_rise / expected_rise - 1.0).abs() < 1e-9,
        "the bar should have risen {expected_rise:.6} K, got {mean_rise:.6}"
    );
}

/// A scene can ask for a schedule its domains cannot survive, and the refusal says which.
///
/// The same scene as above with `staggered` in place of `multirate`. Half a second is
/// thirty-eight times the bar's explicit-diffusion limit, and without subcycling the bar
/// would fill with oscillating nonsense — the classic silently-wrong result.
///
/// It does not. The step is refused, and the refusal names the quantity (`Fourier number`),
/// the site (`bar (explicit conduction)`), the limit and the value. That is the difference
/// this library is for: an application that chose the wrong schedule from a config file is
/// told which domain broke and by how much, rather than drawing something plausible.
///
/// **Also a finding.** Nothing checks this until the first step is taken. A scene format
/// could ask every domain for its `max_stable_dt` at build time and refuse there, where the
/// message could name the file. `Domain::max_stable_dt` is public, so an application can do
/// it; this one does not yet.
#[test]
fn a_schedule_a_scene_cannot_survive_is_named_and_refused() {
    let scene: Scene = serde_json::from_str(
        r#"{
          "title": "too big a step", "schedule": "staggered",
          "duration_s": 4.0, "frames": 8,
          "domains": [
            { "kind": "heater", "name": "element", "watts": 2.0, "reserve_j": 6.0 },
            { "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 21,
              "area_mm2": 100.0, "initial_c": 20.0 }
          ]
        }"#,
    )
    .unwrap();

    let violation = World::build(scene)
        .expect("it builds; the trouble only shows when it runs")
        .run()
        .expect_err("half a second is far past the bar's diffusion limit");

    assert_eq!(violation.quantity, "Fourier number");
    assert!(
        violation.site.starts_with("bar"),
        "the refusal should name the domain, got {:?}",
        violation.site
    );
    // 0.5 is the explicit limit and the scene asked for about 38.
    assert!((violation.before - 0.5).abs() < 1e-12);
    assert!(violation.after > 30.0, "got {}", violation.after);
}

/// **A beam that heats where it lands**, over a boundary two domains share, from data.
///
/// The last part of the kernel with no outside consumer. A plain channel carries an amount;
/// an `Interface` carries an amount *and a place*, and is audited face by face — because a
/// flux redistributed to the wrong part of a boundary keeps the total exactly right, so a
/// total-only check would pass the one bug the spatial coupling exists to prevent.
///
/// The scene builds a bar that exposes a boundary and a beam that publishes onto it. Neither
/// names the other; they name the same boundary.
///
/// The check is a shape, not a total. A total would be satisfied by spreading the joules
/// evenly, which is exactly the failure mode: the middle cell must end up hotter than the
/// ends by the ratio the Gaussian says, computed here rather than read off the bar.
#[test]
fn a_beam_heats_the_bar_where_it_lands() {
    let cells = 41;
    let scene: Scene = serde_json::from_str(&format!(
        r#"{{
          "title": "a beam on a bar", "schedule": "multirate",
          "duration_s": 0.2, "frames": 4,
          "conservation_tolerance": 1e-9,
          "domains": [
            {{ "kind": "beam", "name": "beam", "onto": "face", "faces": {cells},
               "face_area_mm2": 0.5, "watts": 4.0, "reserve_j": 0.4,
               "waist_fraction": 0.12 }},
            {{ "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": {cells},
               "area_mm2": 100.0, "initial_c": 20.0,
               "exposes": {{ "name": "face", "face_area_mm2": 0.5 }} }}
          ]
        }}"#
    ))
    .expect("the spatial scene parses");

    let mut world = World::build(scene).expect("it builds");
    // The audit is face by face, and anything left unclaimed on a spatial channel at the end
    // of a step is refused just as it is on a plain one.
    let frames = world.run().expect("the books close over the boundary");

    let beam = world
        .simulation()
        .domain_as::<dualis_world::beam::Beam>("beam")
        .expect("the beam is still there");
    assert!(
        beam.reserve().to_si() < 1e-12,
        "0.4 J at 4 W runs out in 0.1 s, {} J left",
        beam.reserve().to_si()
    );

    let profile = frames[frames.len() - 1].panels[0].values();
    assert_eq!(profile.len(), cells);

    // The shape. Conduction has had 0.2 s to smear it, so the peak is lower than the
    // Gaussian's and the check is a bound rather than an equality: the middle must still be
    // clearly hotter than the ends, and a uniform spread would put it at exactly 1.0.
    // Ambient in kelvin, because the panel is kelvin: `ScalarField::at` returns what the cells
    // hold and the celsius a picture wants is a view's conversion. Subtracting 20 here — which
    // is what this said while the app did the conversion inside its own sampler — leaves 273.15
    // in both terms and turns a ratio of 4 into a ratio of 1.002. The assertion would still have
    // passed at `> 1.0` and measured nothing.
    let middle = profile[cells / 2] - 293.15;
    let edge = profile[0] - 293.15;
    assert!(
        middle / edge > 3.0,
        "the beam should land in the middle: centre rose {middle:.4} K, end {edge:.4} K, \
         ratio {:.2} — a flat spread would give 1.0",
        middle / edge
    );

    // And the total is exactly what was paid, independent of where it went. 0.4 J into
    // 20 mm x 1 cm^2 of aluminium, insulated, computed from the substance rather than the bar.
    //
    // **Read from the bar and not from the panel**, which is a trap worth naming. A panel is
    // `ScalarField` *sampled* at evenly spaced points including both ends; the bar's grid is
    // cell-centred, so the two end samples sit half a cell outside the outermost centres.
    // Averaging the samples is not averaging the cells, and it comes out about 1/2n low —
    // 1.2% at 41 cells, which is exactly what this assertion caught the first time. An
    // application that reported a mean temperature from its own render buffer would be wrong
    // by that much and would have no way to know. See `FRICTION.md`, finding 10.
    let capacity = Substance::aluminium_6061()
        .heat_capacity(Volume::from_si(20e-3 * 100e-6))
        .expect("aluminium has a specific heat");
    let bar = world
        .simulation()
        .domain_as::<Bar1D>("bar")
        .expect("the bar is still there");
    let mean_rise = bar.mean_temperature().to_si() - Temperature::celsius(20.0).to_si();
    assert!(
        (mean_rise / (0.4 / capacity.to_si()) - 1.0).abs() < 1e-6,
        "the bar holds every joule the beam paid: {mean_rise:.5} K"
    );

    // The sampled panel is close but not equal, and the gap is the sampling and not a leak.
    let sampled_rise = profile.iter().sum::<f64>() / cells as f64 - 293.15;
    let gap = (sampled_rise / mean_rise - 1.0).abs();
    assert!(
        gap > 1e-4 && gap < 0.05,
        "the panel average should differ from the cell average by about 1/2n, got {gap:.4}"
    );
}

/// Two sides that disagree about the boundary are refused, and told both numbers.
///
/// The face count is stated twice in a scene — once as the bar's `cells` and once as the
/// beam's `faces` — because nothing derives one from the other. So it can be stated wrongly,
/// and this is what happens when it is.
///
/// Not a silent renormalisation onto whichever discretisation happened to be first. A flux
/// padded or truncated to fit would put energy on the wrong part of the boundary while
/// keeping the total right, which is the failure the spatial channel exists to prevent and
/// the one a conservation audit cannot see.
#[test]
fn a_boundary_the_two_sides_cut_differently_is_refused() {
    let scene: Scene = serde_json::from_str(
        r#"{
          "title": "mismatched", "schedule": "multirate",
          "duration_s": 0.1, "frames": 2,
          "domains": [
            { "kind": "beam", "name": "beam", "onto": "face", "faces": 20,
              "face_area_mm2": 0.5, "watts": 4.0, "reserve_j": 0.4,
              "waist_fraction": 0.12 },
            { "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 41,
              "area_mm2": 100.0, "initial_c": 20.0,
              "exposes": { "name": "face", "face_area_mm2": 0.5 } }
          ]
        }"#,
    )
    .unwrap();

    let violation = World::build(scene)
        .expect("it builds; the two sides only meet when they run")
        .run()
        .expect_err("20 faces cannot be published onto a boundary cut into 41");

    // The refusal names the boundary and the channel, and carries both counts.
    assert!(
        violation.site.contains("face") && violation.site.contains("energy"),
        "the refusal should name the boundary and channel, got {:?}",
        violation.site
    );
    assert!(
        (violation.before - 41.0).abs() < 0.5 || (violation.after - 41.0).abs() < 0.5,
        "one side of the report should be the 41 faces the bar offered: \
         before {} after {}",
        violation.before,
        violation.after
    );
}

/// **Every scene that ships is run**, because a scene in this repository is a claim.
///
/// The same rule the library's examples follow: one that compiles and then produces nonsense
/// is worse than none at all, and the only way to know is to run it. Running one is not a
/// weak check — the conservation audit is live for the whole run, so a scene that leaked
/// energy or left it unclaimed on a channel would fail here rather than draw something
/// plausible.
///
/// Each also gets one number asserted, chosen to be a property of the physics rather than of
/// the file: what would change if the scene were edited, and what would change if the library
/// broke.
#[test]
fn every_scene_that_ships_runs_and_says_something_true() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes");
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .expect("the scenes directory is there")
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".json"))
        .collect();
    names.sort();
    assert!(
        names.len() >= 5,
        "expected the shipped scenes, found {names:?}"
    );

    for name in &names {
        let text = std::fs::read_to_string(dir.join(name)).unwrap();
        let scene: Scene =
            serde_json::from_str(&text).unwrap_or_else(|e| panic!("{name} does not parse: {e}"));
        let title = scene.title.clone();
        let mut world = World::build(scene).unwrap_or_else(|e| panic!("{name}: {e}"));
        let frames = world
            .run()
            .unwrap_or_else(|v| panic!("{name} ({title}) stopped conserving: {v}"));

        let last = frames.last().expect("a run produces frames");
        // The crate's only guard against a scene that checks nothing. It used to be an
        // out-of-bounds panic inside a loop over ten files, naming neither the scene nor the
        // domain; a shipped scene whose every domain lacked a field would have reported that.
        //
        // A scene may legitimately draw nothing — a thermal network's nodes have capacities and
        // not positions, so `as_field` declines — but then it owes an arm in the match below,
        // and this list is how it says so. Adding a name here without an arm reintroduces
        // exactly the hole: a scene that runs, draws nothing, checks nothing, and passes.
        const NOTHING_TO_DRAW: [&str; 3] = [
            "11-motor-thermal-network.json",
            "12-winding-heats-a-motor.json",
            "13-winding-that-heats-itself.json",
        ];
        // Zero for a scene with no panel, which never reads it — the arms below that do are all
        // in the drawn branch.
        let peak = last.panels.first().map_or(0.0, |p| {
            p.values().iter().fold(0.0f64, |m, v| m.max(v.abs()))
        });
        if NOTHING_TO_DRAW.contains(&name.as_str()) {
            assert!(
                last.panels.is_empty(),
                "{name}: listed as undrawable but it produced a panel — drop it from the list"
            );
        } else {
            assert!(
                !last.panels.is_empty(),
                "{name} ({title}): no domain produced a panel, so there is nothing to draw \
                 and nothing to check"
            );
            assert!(peak.is_finite(), "{name}: the field went to {peak}");
        }

        match name.as_str() {
            // The electro-thermal feedback, closed between frames because it cannot be closed
            // inside the step loop. The claim is the *amplification*: a winding whose resistance
            // follows its own temperature settles hotter than one held at ambient, by
            // 1/(1 − g) where g = I²R₂₀·α·R_th is the loop gain.
            //
            // Checked as a ratio against the same scene without `tracks`, because a ratio
            // cancels most of what the two runs share and isolates the feedback.
            "13-winding-that-heats-itself.json" => {
                let net = world
                    .simulation()
                    .domain_as::<ThermalNetwork>("motor")
                    .expect("the motor is a network");
                let node = |n: &str| net.node_named(n).expect("the node is there");
                let hot = net.temperature(node("winding")).to_si() - 273.15;
                let housing = net.temperature(node("housing")).to_si() - 273.15;

                // The same scene with the feedback removed.
                let mut open = scene_from(&text);
                if let Some(DomainSpec::Winding { tracks, .. }) = open.domains.get_mut(0) {
                    *tracks = None;
                }
                let mut open = World::build(open).expect("the open-loop scene builds");
                open.run().expect("it conserves too");
                let cold = open
                    .simulation()
                    .domain_as::<ThermalNetwork>("motor")
                    .expect("the motor is a network")
                    .temperature(node("winding"))
                    .to_si()
                    - 273.15;

                let measured = (hot - 25.0) / (cold - 25.0);

                // The loop gain, computed here. `R_th` is the series path winding → stator →
                // housing → ambient, and the last term is **not** just convection: at 74.6 °C
                // the housing's linearised radiative conductance `4εσAT³` is 0.036 W/K against
                // 0.294 for convection, which is 12% of the path's weakest link and moves the
                // prediction from 1.310 to 1.280. Leaving it out is a 2.2% error, which is
                // larger than the agreement being asserted — so it is in.
                let (rho_20, alpha, sigma, emissivity, area) =
                    (1.724e-8, 0.00393, 5.670_374_419e-8, 0.09, 0.042);
                let r_20 = rho_20 * 62.0 / 0.35e-6;
                let radiative = 4.0 * emissivity * sigma * area * (housing + 273.15).powi(3);
                let r_th = 1.0 / 0.9 + 1.0 / 2.4 + 1.0 / (7.0 * area + radiative);
                let gain = 2.0 * 2.0 * r_20 * alpha * r_th;
                let want = 1.0 / (1.0 - gain);

                assert!(
                    (measured / want - 1.0).abs() < 5e-3,
                    "{name}: amplification {measured:.4} against {want:.4} from a loop gain of \
                     {gain:.4}"
                );
                // And it is a real effect rather than a rounding one: 16 K of winding.
                assert!(
                    hot - cold > 10.0,
                    "{name}: the feedback is worth only {:.2} K",
                    hot - cold
                );
            }
            // The first scene with more than two domains, and the reason it exists: the crate
            // split's claim is that domains compose, and it had been verified in pairs and
            // never beyond. Building this found two things pairwise coupling cannot reach —
            // a second consumer of one channel silently getting nothing, and the fact that a
            // world's tolerance is set by its loosest domain.
            //
            // Four crates on one bus: optics publishes twice (spatially onto the bar's face,
            // and as a plain amount from the lamp), thermal consumes both, acoustics and
            // mechanics conserve alongside without exchanging anything. All at 1e-9.
            "14-a-world.json" => {
                assert_eq!(world.scene().domains.len(), 5);

                // The bar took from *both* couplings, which is the composition being tested.
                let bar = world
                    .simulation()
                    .domain_as::<Bar1D>("bar")
                    .expect("the bar is there");
                let took = bar.absorbed_energy().to_si();
                assert!(
                    took > 0.5,
                    "{name}: the bar should have absorbed from beam and lamp, got {took:.4} J"
                );

                // Both producers reached it, and neither reached it whole.
                //
                // The beam's 0.4 J is a flux onto a face and arrives entire. The lamp's 2.0 J
                // is *spent*, not delivered: it lands on an aluminium coating that reflects
                // most of it, and only the absorbed fraction becomes heat. Asserting 2.4 J
                // here — which I did first — reads the reserve as though a mirror were a
                // blackbody, and the run said 0.708 instead.
                assert!(
                    took > 0.4,
                    "{name}: {took:.4} J is no more than the beam alone, so the lamp gave nothing"
                );
                assert!(
                    took < 2.4,
                    "{name}: {took:.4} J means the mirror absorbed everything it reflected"
                );

                // The room is still ringing and the orbits still moving, which is what makes
                // this a world rather than a coupling with spectators.
                assert!(peak > 0.1, "{name}: the room went quiet, peak {peak:.4} Pa");
            }
            // The scene that computes its own watts. `11` states 12 W; this one derives them
            // from 62 m of 0.35 mm² copper at 1.75 A, and the point is that the number is now
            // wrong if the geometry is wrong. Checked against `I²R` written out here, with
            // copper's resistivity and coefficient as literals rather than read from the
            // library — otherwise this compares the library with itself.
            "12-winding-heats-a-motor.json" => {
                let rho_20 = 1.724e-8;
                let r_20 = rho_20 * 62.0 / 0.35e-6;
                let r_90 = r_20 * (1.0 + 0.00393 * 70.0);
                let watts = 1.75 * 1.75 * r_90;

                let coil = world
                    .simulation()
                    .domain_as::<Winding>("coil")
                    .expect("the coil is a winding");
                assert!(
                    (coil.dissipation().to_si() / watts - 1.0).abs() < 1e-9,
                    "{name}: dissipating {:.4} W against {watts:.4} W",
                    coil.dissipation().to_si()
                );

                // Evaluated hot, and that is worth 27.5% — the whole reason the temperature is
                // a parameter rather than an omission.
                let cold = 1.75 * 1.75 * r_20;
                assert!(
                    (watts / cold - 1.2751).abs() < 1e-3,
                    "{name}: hot/cold is {:.4}",
                    watts / cold
                );

                // Every joule it spent reached the network, which is the coupling itself.
                let net = world
                    .simulation()
                    .domain_as::<ThermalNetwork>("motor")
                    .expect("the motor is a network");
                assert!(
                    (net.absorbed_energy().to_si() / coil.dissipated_energy().to_si() - 1.0).abs()
                        < 1e-12,
                    "{name}: {} J absorbed against {} J dissipated",
                    net.absorbed_energy().to_si(),
                    coil.dissipated_energy().to_si()
                );

                // And it lands where `11` put it with a stated 12 W, which is the scene's
                // argument: the guess was reasonable, and this one would have caught it if not.
                let winding = net.node_named("winding").expect("there is a winding");
                let hot = net.temperature(winding).to_si() - 273.15;
                assert!((hot - 54.85).abs() < 0.5, "{name}: winding at {hot:.2} C");
            }
            // The ordering a network exists to produce, and the reason a `lump` will not do:
            // heat enters the winding and leaves through the housing, so the temperatures must
            // fall along the chain and every drop must be positive. A transposed index or a
            // one-sided link keeps the ledger exact and breaks this.
            "11-motor-thermal-network.json" => {
                let net = world
                    .simulation()
                    .domain_as::<ThermalNetwork>("motor")
                    .expect("the motor is a thermal network");
                let temps: Vec<(&str, f64)> = net
                    .handles()
                    .map(|(n, label)| (label, net.temperature(n).to_si()))
                    .collect();
                assert_eq!(temps.len(), 3);
                for pair in temps.windows(2) {
                    assert!(
                        pair[0].1 > pair[1].1,
                        "{}: {} at {:.2} K is not above {} at {:.2} K",
                        name,
                        pair[0].0,
                        pair[0].1,
                        pair[1].0,
                        pair[1].1
                    );
                }

                // And the winding is meaningfully hotter than the housing, not hotter by a
                // rounding error — 12 W across 0.9 and 2.4 W/K is 13.3 + 5.0 K at steady state,
                // and this run stops around half a time constant in, so a good part of it.
                let (hot, cold) = (temps[0].1, temps[2].1);
                assert!(
                    hot - cold > 10.0,
                    "{name}: the winding is only {:.2} K above the housing",
                    hot - cold
                );
            }
            // A standing mode keeps its shape and rides |cos|, so it can never exceed the
            // amplitude it was released at. A scheme going unstable shows up here first.
            "01-room-mode.json" | "02-room-higher-mode.json" => {
                assert!(
                    peak <= 1.0 + 1e-9,
                    "{name}: a standing mode cannot exceed its release amplitude, got {peak}"
                );
            }
            // A pulse released from rest splits into waves going every way, so no part of it
            // keeps the full height. It must have spread and it must not have blown up.
            "03-room-pulse.json" => {
                assert!(
                    peak < 0.5 && peak > 0.02,
                    "{name}: a spread pulse should be well under its release height and \
                     still visible, got {peak}"
                );
            }
            // Six joules into 20 mm x 1 cm^2 of aluminium, insulated, is 1.24 K — computed
            // from the substance and not from the bar.
            //
            // Against the *mean*, and read from the domain rather than the panel. Two traps
            // in one assertion, both met while writing it. The peak is 1.30 K, not 1.24,
            // because heat arriving on a plain channel has no place and `Bar1D` puts it in
            // cell 0 — four seconds of conduction have not finished levelling it, and that
            // gradient is the physics rather than an error. And a mean taken over the panel
            // would be about 1/2n low, for the reason in `FRICTION.md` finding 10.
            "04-heater-and-bar.json" => {
                let capacity = Substance::aluminium_6061()
                    .heat_capacity(Volume::from_si(20e-3 * 100e-6))
                    .unwrap();
                let want = 6.0 / capacity.to_si();
                let bar = world
                    .simulation()
                    .domain_as::<Bar1D>("bar")
                    .expect("the bar is still there");
                let mean = bar.mean_temperature().to_si() - Temperature::celsius(20.0).to_si();
                assert!(
                    (mean / want - 1.0).abs() < 1e-6,
                    "{name}: the bar holds every joule: wanted {want:.4} K, got {mean:.4}"
                );
                assert!(
                    peak - 20.0 > mean,
                    "{name}: the fed end should still be above the mean, peak {:.4} against mean {mean:.4}",
                    peak - 20.0
                );
            }
            // The beam lands in the middle, so the middle must be hotter than the ends. A
            // flat spread would make these equal — which is what the scene showed at its
            // first duration of 1.5 s, because 20 mm of aluminium has a diffusion time
            // constant of about 0.59 s and had levelled itself twice over. The scene is
            // 0.2 s now, which spans the beam being on and the spot starting to spread.
            "05-beam-on-bar.json" => {
                let v = last.panels[0].values();
                let (middle, end) = (v[v.len() / 2] - 293.15, v[0] - 293.15);
                assert!(
                    middle > 2.0 * end,
                    "{name}: the beam landed in the middle: {middle:.4} K against {end:.4} K"
                );
            }
            // Kepler's third law, from the picture. The satellites are on circular orbits,
            // so `v = sqrt(GM/r)` and the fastest is the innermost: 7546 m/s at 7000 km
            // against Earth's mass, computed here and not read off the domain.
            "06-orbits.json" => {
                let want = (6.674_30e-11f64 * 5.972e24 / 7.0e6).sqrt();
                assert!(
                    (peak / want - 1.0).abs() < 0.02,
                    "{name}: the innermost should be at {want:.1} m/s, got {peak:.1}"
                );
                let (positions, bounds) = match &last.panels[0].data {
                    dualis_world::PanelData::Points {
                        positions, bounds, ..
                    } => (positions, bounds),
                    _ => panic!("{name}: an orbit is bodies, not a field"),
                };
                // The frame holds the widest orbit, in all three axes.
                assert!(
                    bounds[3] > 2.0e7,
                    "{name}: the frame should hold the widest orbit"
                );
                // And the satellites are genuinely out of one plane, which is what the third
                // axis is for. A flat system would leave every z at zero and the projection
                // would be an expensive way to draw a circle.
                let out_of_plane = positions.iter().map(|p| p[2].abs()).fold(0.0f64, f64::max);
                assert!(
                    out_of_plane > 5.0e6,
                    "{name}: the orbits should be inclined, largest |z| is {out_of_plane:.3e}"
                );
            }
            // The dashpot takes the height away. Restitution is about 0.51, so after a
            // second the ball is on the floor and not moving — and the audit having passed
            // at all means a lump took every joule the contact published.
            "07-bouncing-ball.json" => {
                assert!(
                    peak < 0.2,
                    "{name}: a second is enough bounces to settle, still at {peak:.3} m/s"
                );
            }
            // Equipartition, twice. The mean square speed is `3(N-1)k_B T / N m`, so the
            // liquid at T* = 1.4 must be quicker than the crystal at 0.15 by about
            // sqrt(1.4/0.15) = 3.06. Peaks are noisier than means, so this is a band.
            "08-atoms-crystal.json" => {
                assert!(
                    (0.5..4.0).contains(&peak),
                    "{name}: a cold crystal's fastest atom, got {peak:.3}"
                );
            }
            "09-atoms-liquid.json" => {
                assert!(
                    peak > 3.0,
                    "{name}: at nine times the temperature the atoms are quicker, got {peak:.3}"
                );
            }
            // Every joule the lamp paid arrived in the mirror.
            //
            // Against what the lamp actually spent, not against its 12 J budget: at 3.6 W of
            // absorbed light the three-second run only gets through about 10.9 J, so asserting
            // the budget would be asserting the scene's arithmetic rather than the coupling's.
            // The absorbed *fraction* — the optics' own answer — is checked separately.
            "10-lamp-on-a-mirror.json" => {
                let capacity = Substance::aluminium_6061()
                    .heat_capacity(Volume::from_si(20e-3 * 100e-6))
                    .unwrap();
                let lamp = world
                    .simulation()
                    .domain_as::<dualis_world::light::Light>("lamp")
                    .expect("the lamp is still there");
                let paid = 12.0 - lamp.reserve().to_si();
                assert!(
                    paid > 8.0,
                    "{name}: the lamp should have spent most of its budget, got {paid:.3} J"
                );
                let want = paid / capacity.to_si();
                let bar = world
                    .simulation()
                    .domain_as::<Bar1D>("mirror")
                    .expect("the mirror bar is still there");
                let mean = bar.mean_temperature().to_si() - Temperature::celsius(20.0).to_si();
                assert!(
                    (mean / want - 1.0).abs() < 1e-6,
                    "{name}: {paid:.4} J is {want:.4} K, got {mean:.4}"
                );
            }
            // **The claim only a three-dimensional model can make.** One cell of a 9x9x9
            // aluminium block starts 60 K hot; six milliseconds later the spot has spread.
            //
            // Three things are checked and each rules out a different wrong model.
            //
            // The mean is exactly `20 + 60/729`, because the faces are insulated and nothing is
            // on the bus — so the block holds every joule it started with and the audit at 1e-9
            // is not the only thing saying so.
            //
            // The spread is **isotropic**: the neighbour one cell away along z is exactly as
            // warm as the one along x. A model that resolved a plane and stacked it, or one that
            // used the wrong spacing on one axis, fails here and passes everything else.
            //
            // And the spot is still a spot. A block that had levelled completely would satisfy
            // both of the above trivially, which is the vacuous version of this scene.
            "15-a-hot-spot-in-a-block.json" => {
                let block = world
                    .simulation()
                    .domain_as::<dualis::thermal::Solid3D>("block")
                    .expect("the block is still there");

                let ambient = Temperature::celsius(20.0).to_si();
                let levelled = ambient + 60.0 / 729.0;
                let mean = block.mean_temperature().to_si();
                assert!(
                    (mean - levelled).abs() < 1e-9,
                    "{name}: insulated, so the mean is fixed at {levelled:.9} K, got {mean:.9}"
                );

                let hot = block.temperature_at(4, 4, 4).to_si();
                let along = |d: (usize, usize, usize)| block.temperature_at(d.0, d.1, d.2).to_si();
                let (x_arm, y_arm, z_arm) = (along((5, 4, 4)), along((4, 5, 4)), along((4, 4, 5)));
                assert!(
                    (x_arm - z_arm).abs() < 1e-9 && (x_arm - y_arm).abs() < 1e-9,
                    "{name}: the spread must be isotropic: x {x_arm:.9}, y {y_arm:.9}, \
                     z {z_arm:.9}"
                );
                assert!(
                    z_arm - ambient > 1.0,
                    "{name}: the third axis should have carried real heat, only {:.4} K",
                    z_arm - ambient
                );

                assert!(
                    hot > z_arm && z_arm > along((4, 4, 6)),
                    "{name}: it should still fall away from the spot: {hot:.4} > {z_arm:.4} > {:.4}",
                    along((4, 4, 6))
                );
                assert!(
                    hot - ambient > 2.0 && hot - ambient < 30.0,
                    "{name}: the spot should be well spread and still visible, {:.3} K above",
                    hot - ambient
                );

                // And the panel is a volume rather than a plane, which is what the capture layer
                // could not express before this domain existed.
                let panel = &last.panels[0];
                assert_eq!(panel.grid(), Some((9, 9, 9)));
                assert_eq!(panel.values().len(), 729, "{name}: a slice is not a solid");
                assert!(
                    panel.slice(0) != panel.slice(4),
                    "{name}: every slice is identical, so z was never sampled"
                );
            }
            // **A heater melting a block of ice**, and the plateau it holds at while it does. The
            // scene a domain with no latent heat cannot express at all: without it the block would
            // sail through zero and be at 190 °C by the end.
            //
            // Three numbers, and none of them is the temperature — which is the point. The melt rate
            // is `P/ρL` and nothing else: a hundred watts against 306 mJ per cubic millimetre is
            // 327 mm³/s, and it is a *straight line* because a phase change has no rate constant in
            // it. Then the plateau ends when the latent heat is paid, and the run's leftover joules
            // warm what is now all liquid by exactly what its capacity says.
            "20-melting-a-block-of-ice.json" => {
                let block = world
                    .simulation()
                    .domain_as::<dualis::thermal::Solid3D>("ice")
                    .expect("the ice is still there");
                let cells = 11.0 * 11.0 * 11.0;
                let cell = 1e-9;
                // Ice, written out here rather than read back from the substance the scene named.
                let (rho, cp, latent) = (917.0, 2050.0, 333_550.0);
                let capacity = cells * cell * rho * cp;
                let fusion = cells * cell * rho * latent;

                // Every cubic millimetre is liquid by the end, and the reading exists at all — a
                // block that cannot melt does not report this column, so its presence is the check
                // that the scene's `"material": "ice"` reached the domain.
                let melted = last
                    .readings
                    .iter()
                    .find(|r| r.label == "melted")
                    .expect("a block that can melt reports how much has");
                assert!(
                    (melted.value - cells).abs() < 1e-9,
                    "{name}: all {cells} mm3 should be liquid, got {}",
                    melted.value
                );

                // The plateau's length is what the two heats cost, and its end is where the
                // temperature leaves zero. 25.0 J of warming then 407.1 J of melting, at 100 W.
                let warming = capacity * 10.0;
                let done_at = (warming + fusion) / 100.0;
                let leftover = (5.0 - done_at) * 100.0;
                let want = leftover / capacity;
                let got = block.mean_temperature().to_si() - 273.15;
                println!(
                    "  {name}: melting took {done_at:.3} s of the 5, and the remaining \
                     {leftover:.1} J warmed it {want:.4} K — measured {got:.4}, off {:.2e}",
                    (got / want - 1.0).abs()
                );
                // `1e-4`, which is the resolution the delivery has: the heater pays in whole steps
                // and the last frame lands where it lands. A first draft allowed 5%, on numbers that
                // agree to five digits.
                assert!(
                    (got / want - 1.0).abs() < 1e-4,
                    "{name}: the leftover joules warm the water by {want:.4} K, got {got:.4}"
                );

                // And the rate in between is `P/ρL`, straight. Read off two frames well inside the
                // plateau, where nothing else is happening.
                let at = |n: usize| {
                    let f = &frames[n];
                    (
                        f.time_s,
                        f.readings
                            .iter()
                            .find(|r| r.label == "melted")
                            .expect("melted")
                            .value,
                    )
                };
                let (t1, v1) = at(2);
                let (t2, v2) = at(8);
                let slope = (v2 - v1) / (t2 - t1);
                let closed = 100.0 / (rho * latent * cell);
                println!(
                    "  and it melted at {slope:.2} mm3/s against P/rho L = {closed:.2} mm3/s — off \
                     {:.2e}",
                    (slope / closed - 1.0).abs()
                );
                // **Machine precision, and that is not luck.** Inside the plateau every joule the
                // heater pays goes to melting and none to warming, so the discrete slope *is* `P/ρL`
                // — the scheme has no discretisation error here at all, and 3.33e-16 is the last bit
                // of the division. A tolerance of a percent would have passed for a scheme that was
                // merely close.
                assert!(
                    (slope / closed - 1.0).abs() < 1e-12,
                    "{name}: the melt rate is the power over the latent heat: {slope:.2} against \
                     {closed:.2} mm3/s"
                );
            }
            // **A wall of glass halfway down a block of aluminium**, and the heat piling up
            // against it. The scene a single-material block cannot express.
            //
            // What is checked is where the *gradient* is, which is the only thing a two-material
            // block says that a one-material block does not. Aluminium's 150x conductivity leaves
            // it nearly isothermal over its own 9 mm, and borosilicate's diffusion length in half
            // a second is half a cell — so the whole temperature drop sits on one face, and the
            // largest cell-to-cell step along z must be exactly the one at the material interface.
            //
            // A block with the coating quietly not applied would put its largest step next to the
            // hot spot instead, which is the failure this catches and a peak temperature would not.
            "19-a-coating-stops-the-heat.json" => {
                let block = world
                    .simulation()
                    .domain_as::<dualis::thermal::Solid3D>("joint")
                    .expect("the joint is still there");
                let ambient = Temperature::celsius(20.0).to_si();
                let at = |k: usize| block.temperature_at(4, 4, k).to_si();

                // The books, not the mean: with two capacities it is `Σ Cᵢ Tᵢ` that is fixed, and
                // the mean temperature is not. 60 K in one aluminium cell of 2.4192 mJ/K.
                let opening = 2700.0 * 896.0 * 1e-9 * 60.0;
                let held = block
                    .ledger()
                    .get(quantity::ENERGY)
                    .expect("energy is on the books");
                assert!(
                    (held - opening).abs() < 1e-9 * opening,
                    "{name}: insulated, so it still holds {opening:.9} J, got {held:.9}"
                );

                let steps: Vec<f64> = (0..17).map(|k| at(k) - at(k + 1)).collect();
                let (worst, _) = steps.iter().enumerate().fold((0, 0.0f64), |acc, (k, s)| {
                    if *s > acc.1 {
                        (k, *s)
                    } else {
                        acc
                    }
                });
                assert_eq!(
                    worst, 8,
                    "{name}: the gradient belongs on the interface face, not at cell {worst}: \
                     {steps:?}"
                );

                // The metal is nearly isothermal and the glass is nearly untouched, which is the
                // same statement read from either side.
                let (near, far) = (at(0) - ambient, at(8) - ambient);
                assert!(
                    (near - far).abs() < 0.2 * far,
                    "{name}: 150x the conductivity should level the metal: {near:.3} K at the far \
                     face against {far:.3} K at the interface"
                );
                assert!(
                    at(12) - ambient < 0.02 * far,
                    "{name}: the glass should have stopped it: {:.4} K four cells in against \
                     {far:.3} K in the metal",
                    at(12) - ambient
                );

                let panel = &last.panels[0];
                assert_eq!(panel.grid(), Some((9, 9, 18)));
            }
            // **The mode a floor plan does not have.** A 4.4 x 3.1 x 2.4 m room released in its
            // oblique (1,1,1) mode, which needs all three axes at once.
            //
            // The peak of a standing mode rides `|cos(2 pi f t)|` exactly, and `f` is the
            // rigid-wall closed form — computed here from the **quantised** dimensions, because
            // the grid makes the ceiling 3.2 m rather than the 3.1 m the file asks for and a
            // closed form about the wrong room is not a check.
            //
            // 97.46 Hz, so 0.02 s is 1.949 periods and the peak should be back near 0.949.
            "16-a-room-with-a-ceiling.json" => {
                let hall = world
                    .simulation()
                    .domain_as::<dualis::acoustic::Hall>("hall")
                    .expect("the hall is still there");
                let (lx, ly, lz) = (
                    hall.width().to_si(),
                    hall.height().to_si(),
                    hall.depth().to_si(),
                );
                let f = 343.0 / 2.0
                    * ((1.0f64 / lx).powi(2) + (1.0 / ly).powi(2) + (1.0 / lz).powi(2)).sqrt();
                let want = (2.0 * std::f64::consts::PI * f * last.time_s).cos().abs();
                assert!(
                    (peak - want).abs() < 0.02,
                    "{name}: a standing mode rides |cos(2 pi f t)|: {peak:.4} against                      {want:.4} at {f:.2} Hz"
                );
                assert!(
                    peak <= 1.0 + 1e-9,
                    "{name}: and can never exceed its release amplitude, got {peak}"
                );

                // The vertical mode is the thing `DomainSpec::Room` cannot express at all. It is
                // not a smaller number there; it is absent.
                let vertical = hall.mode_frequency((0, 0, 1)).to_si();
                assert!(
                    (vertical - 343.0 / (2.0 * lz)).abs() < 1e-9
                        && (60.0..85.0).contains(&vertical),
                    "{name}: the floor-to-ceiling mode is c/2Lz, got {vertical:.2} Hz"
                );

                // And the panel is a volume, sampled at the grid's own node count.
                assert_eq!(last.panels[0].grid(), Some(hall.nodes()));
            }
            // **A resistance that no formula gives.** A 12 x 5 x 5 mm copper busbar with a
            // notch three cells deep across the middle, driven at 1 mV.
            //
            // Two bounds, both provable, and neither is a value: `rho L/A` for the full section
            // is a floor, because removing conductor cannot help; and a naive series estimate
            // that treats the notched slice as a shorter bar of reduced section is *also* a
            // floor, because the current has to spread back out and spreading costs. The excess
            // over the second is the spreading resistance, and it has no closed form for this
            // shape -- which is the entire reason to solve rather than to state.
            //
            // Measured 12.392 uohm, against 8.275 for the full section and 9.310 for the naive
            // series. A bound rather than the measurement, because the measurement is what the
            // code produced and a test that asserts it checks nothing.
            "17-a-busbar-with-a-notch.json" => {
                let bar = world
                    .simulation()
                    .domain_as::<dualis::electrical::Conductor>("busbar")
                    .expect("the busbar is still there");
                assert!(bar.converged(), "residual {:.3e}", bar.residual());

                let (rho, dx) = (1.724e-8, 1e-3);
                let full = rho * (12.0 * dx) / (25.0 * dx * dx);
                let naive = rho * (11.0 * dx) / (25.0 * dx * dx) + rho * dx / (10.0 * dx * dx);
                let got = bar.resistance().to_si();
                assert!(
                    got > full * 1.2,
                    "{name}: a notch must cost: {got:.4e} against rho L/A = {full:.4e}"
                );
                assert!(
                    got > naive,
                    "{name}: spreading costs more than a series estimate:                      {got:.4e} against {naive:.4e}"
                );

                // Tellegen, through the scene format: the power from the field equals the power
                // at the terminals.
                let terminal = bar.drive().to_si() * bar.current().to_si();
                assert!(
                    (bar.dissipation().to_si() / terminal - 1.0).abs() < 1e-9,
                    "{name}: field power against terminal power"
                );

                // And every joule it paid, the heatsink took -- which the audit at 1e-9 also
                // says, from the other side.
                let sink = last
                    .readings
                    .iter()
                    .find(|r| r.domain == "heatsink" && r.label == "absorbed")
                    .expect("the heatsink reports what it absorbed");
                assert!(
                    (sink.value / bar.dissipated_energy().to_si() - 1.0).abs() < 1e-9,
                    "{name}: {} J absorbed against {} J paid",
                    sink.value,
                    bar.dissipated_energy().to_si()
                );

                // The panel is the potential, a volume, at the conductor's own cell count.
                assert_eq!(last.panels[0].grid(), Some(bar.counts()));
                assert_eq!(last.panels[0].unit, "V");
            }
            // **Two baskets, identical but for the ring against the wall**, and the whole
            // difference between them is a permeability that nobody stated.
            //
            // The flow ratio is exact rather than approximate, and the reason is worth stating:
            // every column of cells from inlet to outlet is an independent series chain, so the
            // basket is columns in **parallel** and its conductance is their sum. Widening a
            // column changes only its own term. So
            //
            //     Q_gap / Q_even  =  (1 - f) + f * m(0.60)/m(0.45),    m(e) = e^3/(1-e)^2
            //
            // with `f` the fraction of the cross-section the ring covers. Kozeny-Carman gives
            // the mobility ratio as 4.482; `f` is **counted** rather than estimated, because
            // estimating it from `2 pi r / pi r^2` was wrong by a third -- a staircase ring on
            // a 15-cell radius is not a circle -- and the test was asserting a bound that the
            // right answer failed.
            //
            // The pair is what makes it a channel rather than merely a faster bed: more liquid,
            // and *less* coffee in it. Both directions are asserted, because either alone is
            // ambiguous -- a bed that is simply coarser all through would give the first and not
            // the second.
            "18-an-espresso-shot.json" => {
                let read = |domain: &str, label: &str| {
                    last.readings
                        .iter()
                        .find(|r| r.domain == domain && r.label == label)
                        .unwrap_or_else(|| panic!("{name}: {domain} reports {label}"))
                        .value
                };
                let (even_g, bad_g) = (read("even", "delivered"), read("wall gap", "delivered"));
                let (even_tds, bad_tds) = (read("even", "TDS"), read("wall gap", "TDS"));
                let (even_ring, bad_ring) = (
                    read("even", "ring over core"),
                    read("wall gap", "ring over core"),
                );

                let puck = world
                    .simulation()
                    .domain_as::<dualis::porous::Puck>("even")
                    .expect("the basket is still there");
                let (nx, _, nz) = puck.counts();
                let (mut packed, mut ring) = (0usize, 0usize);
                for kk in 0..nz {
                    for i in 0..nx {
                        if !puck.is_packed(i, 0, kk) {
                            continue;
                        }
                        packed += 1;
                        if i == 0
                            || kk == 0
                            || i + 1 == nx
                            || kk + 1 == nz
                            || !puck.is_packed(i - 1, 0, kk)
                            || !puck.is_packed(i + 1, 0, kk)
                            || !puck.is_packed(i, 0, kk - 1)
                            || !puck.is_packed(i, 0, kk + 1)
                        {
                            ring += 1;
                        }
                    }
                }
                let f = ring as f64 / packed as f64;
                let mobility = |e: f64| e.powi(3) / (1.0 - e).powi(2);
                let predicted = (1.0 - f) + f * mobility(0.60) / mobility(0.45);
                let measured = bad_g / even_g;
                assert!(
                    (measured / predicted - 1.0).abs() < 1e-6,
                    "{name}: columns in parallel give the flow ratio exactly: {measured:.6}x \
                     against {predicted:.6}x, with the ring {:.1}% of the section",
                    f * 100.0
                );
                assert!(
                    bad_tds < even_tds,
                    "{name}: and carry less coffee in it: {bad_tds:.3}% against {even_tds:.3}%"
                );

                // The diagnosis, which is the reading that separates the two hypotheses.
                assert!(
                    (even_ring - 1.0).abs() < 0.02,
                    "{name}: an evenly packed basket extracts its ring and its core alike: \
                     {even_ring:.4}"
                );
                assert!(
                    bad_ring > 1.05,
                    "{name}: and the gap's ring must outrun the core it starved: {bad_ring:.4} \
                     against {even_ring:.4}"
                );

                // Darcy in closed form, through the scene format. Nothing in the file is a flow
                // rate; this is what the permeability the grind gives actually produces.
                let mu = dualis::porous::Liquid::water()
                    .viscosity(dualis::units::Temperature::celsius(93.0))
                    .to_si();
                let k = dualis::porous::Grind::sieved(dualis::units::Length::from_si(250e-6))
                    .permeability(0.45)
                    .to_si();
                let ny = puck.counts().1;
                let dx = puck.spacing().to_si();
                let closed = k * (packed as f64 * dx * dx) * 9.0e5 / (mu * ny as f64 * dx);
                let measured =
                    puck.flow_rate().to_si() / dualis::porous::Liquid::water().density.to_si();
                assert!(
                    (measured / closed - 1.0).abs() < 1e-9,
                    "{name}: Q = kA dp / (mu L): {measured:.6e} against {closed:.6e}"
                );

                // The panel is a volume at the basket's own cells, and it is the extraction
                // rather than the temperature -- a bed under flow is isothermal, so a
                // temperature panel would be a flat rectangle that renders and says nothing.
                assert_eq!(last.panels[0].grid(), Some(puck.counts()));
                assert_eq!(last.panels[0].unit, "");
            }
            "21-a-wax-thermal-buffer.json" => {
                // n-octadecane, and every number here is in the scene file rather than in the
                // library. That is the point of the scene: nothing in `dualis` knows this substance.
                let (rho, cp, latent, melting_c) = (814.0, 1934.0, 244_000.0, 301.3 - 273.15);
                // `melted` is reported in mm3 and the cells are 2 mm, so the two conversions are
                // different numbers and mixing them up is easy: one cell is 8 mm3, and a cubic
                // millimetre is 1e-9 m3 whatever the cell size. A first draft used 8e-9 as the
                // mm3-to-m3 factor and reported the block as 402% melted.
                const MM3: f64 = 1e-9;
                let total_mm3 = 11.0 * 11.0 * 11.0 * 8.0;
                let mass = total_mm3 * MM3 * rho;

                let reading = |f: &dualis_world::Frame, label: &str| {
                    f.readings
                        .iter()
                        .find(|r| r.label == label)
                        .unwrap_or_else(|| panic!("{name}: no {label} reading"))
                        .value
                };

                // **The plateau is at the declared melting point**, and that number appears nowhere
                // except the file. Ice's plateau is 0.000 and every catalogue material has none at
                // all, so this single assertion is what says the declaration reached the domain
                // rather than being parsed and dropped in favour of the aluminium default.
                let held = world
                    .simulation()
                    .domain_as::<dualis::thermal::Solid3D>("wax")
                    .expect("the wax is still there")
                    .mean_temperature()
                    .to_si()
                    - 273.15;
                assert!(
                    (held - melting_c).abs() < 1e-9,
                    "{name}: the plateau should sit at the declared {melting_c} C, is {held}"
                );

                // Half melted and stopped there, because the reserve ran out. Both halves matter:
                // a block that melts entirely would not show the plateau ending, and one that melts
                // nothing would look identical to a declaration that never took.
                let melted = reading(last, "melted");
                let fraction = melted / total_mm3;
                assert!(
                    (0.45..0.55).contains(&fraction),
                    "{name}: about half should be liquid, {fraction:.4} is"
                );

                // The rate is `P/rho L` with the density and the latent heat coming out of the file.
                // Read between two frames well inside the plateau and after the sensible warming is
                // done, where every joule arriving goes to melting and none to warming.
                let slope = (reading(&frames[5], "melted") - reading(&frames[2], "melted"))
                    / (frames[5].time_s - frames[2].time_s);
                let closed = 20.0 / (rho * latent * MM3);
                println!(
                    "  {name}: melted at {slope:.4} mm3/s against P/rho L = {closed:.4} — off {:.2e}",
                    (slope / closed - 1.0).abs()
                );
                // **Machine precision, for the reason `20` states.** Inside the plateau every joule
                // the heater pays goes to melting and none to warming, so the discrete slope *is*
                // `P/rho L` and the scheme has no discretisation error here at all. Measured 8.3e-15.
                // A tolerance of 1e-4 would have been eleven orders too loose and would have passed a
                // density read from the wrong column.
                assert!(
                    (slope / closed - 1.0).abs() < 1e-12,
                    "{name}: {slope:.4} mm3/s against a closed-form {closed:.4}"
                );

                // The energy accounts, and this is the check a wrong specific heat could not pass:
                // the sensible term uses `c_p` and the latent term does not, so the two are only
                // consistent with the reserve if both declared numbers are right.
                let sensible = mass * cp * (melting_c - 20.0);
                let fusion = melted * MM3 * rho * latent;
                let absorbed = reading(last, "absorbed");
                println!(
                    "  and {absorbed:.4} J absorbed = {sensible:.4} warming + {fusion:.4} melting                      = {:.4}, off {:.2e}",
                    sensible + fusion,
                    ((sensible + fusion) / absorbed - 1.0).abs()
                );
                // 1.13e-14, and exact for the same reason: in the enthalpy method energy *is* the
                // state, so the split between sensible and latent is bookkeeping rather than
                // approximation. `1e-12` and not `1e-3`, because the looser number would pass a
                // specific heat that was 8% wrong — the sensible term is a ninth of the total.
                assert!(
                    ((sensible + fusion) / absorbed - 1.0).abs() < 1e-12,
                    "{name}: {absorbed:.4} J absorbed against {:.4} accounted for",
                    sensible + fusion
                );

                // And once the reserve is empty nothing moves. A block still melting after its
                // heater is exhausted is inventing energy, which the audit would catch on the sum
                // and not on this domain — so it is worth saying here where it is unambiguous.
                let empty: Vec<f64> = frames
                    .iter()
                    .filter(|f| reading(f, "melted") > 0.0 && f.time_s > 70.0)
                    .map(|f| reading(f, "melted"))
                    .collect();
                assert!(
                    empty.len() >= 4,
                    "{name}: not enough frames after the reserve ran dry"
                );
                for v in &empty {
                    assert_eq!(
                        *v, melted,
                        "{name}: melting continued after the heater was empty"
                    );
                }
            }
            "22-wax-in-an-aluminium-matrix.json" => {
                // The same wax as `21`, now four fifths of a composite whose other fifth is aluminium.
                // Every number here comes from the two constituents; nothing is read back from the scene.
                const MM3: f64 = 1e-9;
                let (rho_w, c_w, l_w, phi_w) = (814.0, 1934.0, 244_000.0, 0.8);
                let (rho_a, c_a, phi_a) = (2700.0, 896.0, 0.2);
                let melting_c = 301.3 - 273.15;
                let total_mm3 = 11.0 * 11.0 * 11.0 * 8.0;

                let reading = |f: &dualis_world::Frame, label: &str| {
                    f.readings
                        .iter()
                        .find(|r| r.label == label)
                        .unwrap_or_else(|| panic!("{name}: no {label} reading"))
                        .value
                };

                // The plateau sits at the melting part's own melting point, **undiluted**. A latent heat
                // dilutes with the mass fraction and a temperature does not, and mixing the two rules up
                // is the mistake `Mix::fusion` exists to prevent.
                let held = world
                    .simulation()
                    .domain_as::<dualis::thermal::Solid3D>("buffer")
                    .expect("the buffer is still there")
                    .mean_temperature()
                    .to_si()
                    - 273.15;
                assert!(
                    (held - melting_c).abs() < 1e-9,
                    "{name}: the plateau should sit at the wax's own {melting_c} C, is {held}"
                );

                // **The melting rate is `P/(phi_w rho_w L_w)`**, which is the composite's volumetric latent
                // heat and nothing else. Read between two frames well inside the plateau.
                let slope = (reading(&frames[5], "melted") - reading(&frames[2], "melted"))
                    / (frames[5].time_s - frames[2].time_s);
                let closed = 20.0 / (phi_w * rho_w * l_w * MM3);
                println!(
                    "  {name}: melted at {slope:.4} mm3/s against P/(phi rho L) = {closed:.4} — off                      {:.2e}",
                    (slope / closed - 1.0).abs()
                );
                assert!(
                    (slope / closed - 1.0).abs() < 1e-12,
                    "{name}: {slope:.4} mm3/s against a closed-form {closed:.4}"
                );

                // **And it is exactly `1/phi_w` times scene 21's rate**, which is the statement only
                // having both scenes can make: the ratio of the two has no material property left in it
                // at all. Diluting the wax by a fifth makes a cubic millimetre of composite hold a fifth
                // less latent heat, so the same watts clear it a quarter faster.
                let pure = 20.0 / (rho_w * l_w * MM3);
                println!(
                    "  and {:.6}x scene 21's {pure:.4} mm3/s, against 1/phi = {:.6}",
                    slope / pure,
                    1.0 / phi_w
                );
                assert!(
                    (slope / pure * phi_w - 1.0).abs() < 1e-12,
                    "{name}: the ratio to pure wax should be 1/phi_w = {:.6}, is {:.6}",
                    1.0 / phi_w,
                    slope / pure
                );

                // The energy accounts, and this is where the volume-and-mass rules are both exercised at
                // once: the sensible term uses the volume-additive `rho c` and the latent term uses the
                // mass-diluted `L`. Getting either weighting wrong breaks the sum.
                let volumetric = phi_w * rho_w * c_w + phi_a * rho_a * c_a;
                let sensible = volumetric * total_mm3 * MM3 * (melting_c - 20.0);
                let melted = reading(last, "melted");
                let fusion = melted * MM3 * phi_w * rho_w * l_w;
                let absorbed = reading(last, "absorbed");
                println!(
                    "  and {absorbed:.4} J = {sensible:.4} warming + {fusion:.4} melting = {:.4}, off                      {:.2e}",
                    sensible + fusion,
                    ((sensible + fusion) / absorbed - 1.0).abs()
                );
                assert!(
                    ((sensible + fusion) / absorbed - 1.0).abs() < 1e-12,
                    "{name}: {absorbed:.4} J absorbed against {:.4} accounted for",
                    sensible + fusion
                );

                // About three fifths melted, and stopped: the reserve is the limit rather than the block.
                let fraction = melted / total_mm3;
                assert!(
                    (0.55..0.70).contains(&fraction),
                    "{name}: about three fifths should be liquid, {fraction:.4} is"
                );
                for f in frames.iter().filter(|f| f.time_s > 70.0) {
                    assert_eq!(
                        reading(f, "melted"),
                        melted,
                        "{name}: melting continued after the heater was empty"
                    );
                }
            }
            "23-a-part-radiating-to-its-lid.json" => {
                // A hot part sitting in a sealed housing, surrounded on five sides by **nothing**
                // and facing a cooled lid across a clearance. The geometry is what makes this
                // checkable: the part has no conducting face to anywhere, so the only path its
                // heat has is the parallel-plate exchange across the gap, and the lid's only exit
                // is the film on `z-max`. That is a two-body lumped system whose closed form is a
                // pair of ODEs in nothing but constants.
                const SIGMA: f64 = 5.670_374_419e-8;
                let cell = 8e-3_f64; // metres, from the scene
                let (rho, c, emissivity) = (2700.0_f64, 896.0_f64, 0.85_f64);
                let (h, ambient) = (12.0_f64, 293.15_f64);

                let face = cell * cell;
                let gap_area = 16.0 * face; // the 4x4 part looking up at the lid
                let top_area = 36.0 * face; // the whole 6x6 lid looking at the room
                let part_capacity = 48.0 * cell.powi(3) * rho * c; // 4 x 4 x 3 cells
                let lid_capacity = 72.0 * cell.powi(3) * rho * c; // 6 x 6 x 2 cells
                let series = 2.0 / emissivity - 1.0;

                let reading = |f: &dualis_world::Frame, label: &str| {
                    f.readings
                        .iter()
                        .find(|r| r.label == label)
                        .unwrap_or_else(|| panic!("{name}: no {label} reading"))
                        .value
                        + 273.15
                };

                // **The part really does cool, and only radiation can have done it.** Without the
                // gap exchange this assertion is false by construction: void breaks conduction, so
                // a part with no radiating clearance sits at its starting temperature forever.
                let start = reading(&frames[0], "peak");
                let end = reading(frames.last().expect("frames"), "peak");
                assert!(
                    (start - 573.15).abs() < 0.1,
                    "{name}: the part starts where the region says, {start:.2} K"
                );
                assert!(
                    start - end > 20.0,
                    "{name}: the clearance should have carried tens of kelvin, carried              {:.2} K",
                    start - end
                );

                // **And the amount is the lumped two-body answer**, integrated here from the
                // constants above with RK4 — no property is read back from the library, and the
                // library's own arithmetic is cellwise rather than lumped, so this is a check
                // against the physics and not against a second copy of the code.
                let rates = |t_part: f64, t_lid: f64| {
                    let radiated = SIGMA * gap_area * (t_part.powi(4) - t_lid.powi(4)) / series;
                    let filmed = h * top_area * (t_lid - ambient);
                    (
                        -radiated / part_capacity,
                        (radiated - filmed) / lid_capacity,
                    )
                };
                let (mut t_part, mut t_lid) = (573.15_f64, 293.15_f64);
                let dt = 0.01;
                for _ in 0..60_000 {
                    let (a1, b1) = rates(t_part, t_lid);
                    let (a2, b2) = rates(t_part + 0.5 * dt * a1, t_lid + 0.5 * dt * b1);
                    let (a3, b3) = rates(t_part + 0.5 * dt * a2, t_lid + 0.5 * dt * b2);
                    let (a4, b4) = rates(t_part + dt * a3, t_lid + dt * b3);
                    t_part += dt / 6.0 * (a1 + 2.0 * a2 + 2.0 * a3 + a4);
                    t_lid += dt / 6.0 * (b1 + 2.0 * b2 + 2.0 * b3 + b4);
                }

                // The lumped model's error is the gradient it assumes away, and both gradients here
                // are small on purpose: the part's diffusion time over its 24 mm is `L^2/alpha` =
                // 8 s and the lid's over its 48 mm span is 33 s, against a 600 s run. What is left
                // is that the library radiates cell by cell on `T^4` where this averages first, and
                // Jensen's inequality makes that a second-order term in the spread. Judged on the
                // **cooling**, which is the quantity the scene is about; as a fraction of the
                // absolute temperature it would be twenty times looser and would notice nothing.
                //
                // Measured at **0.20%**, and the bound is 1% rather than the 5% the lag argument
                // alone would buy, because a tolerance twenty-five times the residual is one that
                // would sit still through a real change to the exchange.
                let cooled = start - end;
                let closed = 573.15 - t_part;
                println!(
                    "  {name}: the part shed {cooled:.2} K against a lumped {closed:.2} K — off {:.2}%",
                    100.0 * (cooled / closed - 1.0).abs()
                );
                assert!(
                    (cooled / closed - 1.0).abs() < 0.01,
                    "{name}: {cooled:.2} K shed against a closed-form {closed:.2} K"
                );

                // **The lid warms, and stays far below the part**, which is what says the exchange
                // ran one way. A gap that conducted rather than radiated would have levelled them.
                let lid = reading(frames.last().expect("frames"), "coldest");
                assert!(
                    lid > 293.15 + 1.0 && lid < end - 100.0,
                    "{name}: the lid should warm a little and stay cold, {lid:.2} K against a part at {end:.2} K"
                );
            }
            "24-a-power-module-junction-to-ambient.json" => {
                // A die dissipating 45 W through the stack under it — solder, DBC ceramic, copper,
                // a cold plate. The scene the format could not state at all until `dissipation`
                // existed: every other source here hands watts to the bus, and the bus carries an
                // amount and no location, so a die's heat would have spread over the baseplate as
                // fast as over the die and there would have been no junction temperature to read.
                //
                // **At steady state the answer is a resistance stack**, and every term in it is
                // written out below from the geometry and the conductivities. Nothing is read back
                // from the scene: the layers are named here, the harmonic face mean is the two
                // half-cells in series, and the film carries the last half cell of copper.
                let dx = 1.5e-3;
                let area = (8.0 * dx) * (8.0 * dx);
                let watts = 45.0;
                let ambient = 40.0;
                // k = 0,1,2 copper baseplate | 3 alumina | 4 copper | 5 solder | 6,7 silicon
                let k_of = |layer: usize| match layer {
                    3 => 24.0,      // Al2O3 96%
                    5 => 58.0,      // SAC305
                    6 | 7 => 148.0, // Si
                    _ => 401.0,     // Cu ETP
                };

                // Centre to centre, half a cell of each material per face. This *is* the harmonic
                // mean the sweep uses, written the other way round — as two resistances in series
                // rather than as one averaged conductivity — so agreeing with it is a statement
                // about the physics and not a repeat of the implementation.
                let mut resistance = 0.0;
                for upper in (1..8).rev() {
                    resistance +=
                        (dx / 2.0) / (k_of(upper) * area) + (dx / 2.0) / (k_of(upper - 1) * area);
                }
                // The bottom cell's own half, then the film.
                resistance += (dx / 2.0) / (k_of(0) * area);
                resistance += 1.0 / (3000.0 * area);

                let junction = ambient + watts * resistance;
                let peak = frames
                    .last()
                    .expect("frames")
                    .readings
                    .iter()
                    .find(|r| r.label == "peak")
                    .expect("a block reports its peak")
                    .value;
                println!(
                    "  {name}: junction {peak:.2} C against a {resistance:.4} K/W stack giving                      {junction:.2} C — off {:.2e}",
                    (peak - junction).abs() / (junction - ambient)
                );
                assert!(
                    (peak - junction).abs() / (junction - ambient) < 2e-3,
                    "{name}: {peak:.3} C against the stack's {junction:.3} C"
                );

                // **And it is steady**, or the agreement above is a coincidence of run length.
                // The last two frames are 15 s apart and the junction must have stopped moving.
                let earlier = frames[frames.len() - 2]
                    .readings
                    .iter()
                    .find(|r| r.label == "peak")
                    .expect("a peak")
                    .value;
                assert!(
                    (peak - earlier).abs() < 0.05,
                    "{name}: still climbing at the end, {earlier:.3} then {peak:.3} C"
                );

                // **The heat leaves where the scene said it does.** A stack whose film was on the
                // wrong face would still reach a steady temperature, and a wrong one; the die must
                // be the hottest thing and the baseplate the coldest.
                let block = world
                    .simulation()
                    .domain_as::<dualis::thermal::Solid3D>("module")
                    .expect("the module is a block");
                assert!(
                    block.temperature_at(4, 4, 7).to_si() > block.temperature_at(4, 4, 0).to_si(),
                    "{name}: the die should be hotter than the baseplate"
                );
                assert!(
                    (block.generated_power().to_si() - watts).abs() < 1e-12,
                    "{name}: the scene states {watts} W and the block has {} W",
                    block.generated_power().to_si()
                );
            }
            "25-what-140-kelvin-does-to-the-solder.json" => {
                // The same power module as `24`, with a body on the same grid whose stress-free
                // strain is the block's temperature. The first scene in this format that couples
                // two physics on one mesh, and the first that reaches `dualis-elastic` at all.
                //
                // **The failure mode of a real module is not its temperature.** It is the solder,
                // fatigued by silicon at 2.6e-6 per kelvin sitting on solder at 2.15e-5 through
                // every power cycle. The platform could compute the temperature to four figures
                // and could say nothing about that until a structure could follow a block.
                let reading = |f: &dualis_world::Frame, label: &str| {
                    f.readings
                        .iter()
                        .find(|r| r.label == label)
                        .unwrap_or_else(|| panic!("{name}: no {label} reading"))
                        .value
                };

                // **The first frame is exact and it is not zero.** The block starts uniform at
                // 40 C and the body is the size it was drawn at 217 C — its solder's reflow — so
                // every element is already strained before anything is switched on. The peak is
                // the solder's own, `alpha (T - T_ref)`, computed here from the scene's declared
                // coefficient and nothing read back.
                let (reference_c, start_c, solder_alpha) = (217.0, 40.0, 2.15e-5);
                let closed = solder_alpha * (start_c - reference_c);
                let at_start = reading(&frames[0], "free strain");
                assert!(
                    (at_start - closed).abs() < 1e-12,
                    "{name}: cold, the solder wants {closed:e} and has {at_start:e}"
                );
                // An equilibrium rather than a body handed a strain and not yet asked. Judged
                // against the solve's own tolerance and not against zero: an unsolved body reports
                // an **infinite** residual, so this separates the two cases without asserting that
                // conjugate gradients terminate exactly.
                let residual = reading(&frames[0], "residual");
                assert!(
                    residual.is_finite() && residual < 1e-9,
                    "{name}: the first frame must be an equilibrium: residual {residual:e}"
                );

                // **And the last frame is exact too**, against the temperature the block reached.
                // This is the claim that says the coupling is live rather than frozen at build.
                let block = world
                    .simulation()
                    .domain_as::<dualis::thermal::Solid3D>("module")
                    .expect("the module is a block");
                // Scanned over every element, exactly as the domain does, with the coefficient
                // the scene's layers imply. The peak is **not** the hottest cell of any layer: all
                // the strains are negative, so the largest magnitude belongs to whatever sits
                // furthest *below* the reference — and which layer that is depends on the product
                // of the coefficient and the distance, not on either alone. Getting it by hand was
                // wrong twice before it was scanned.
                let copper_alpha = dualis::prelude::Substance::copper()
                    .thermal
                    .expect("copper expands")
                    .expansion
                    .to_si();
                let alpha_at = |k: usize| match k {
                    3 => 7.4e-6,       // Al2O3 96%
                    5 => solder_alpha, // SAC305
                    6 | 7 => 2.6e-6,   // Si
                    _ => copper_alpha, // Cu ETP, the block's bulk
                };
                let mut want = 0.0f64;
                for k in 0..8 {
                    for j in 0..8 {
                        for i in 0..8 {
                            let t = block.temperature_at(i, j, k).to_si() - 273.15;
                            let e = alpha_at(k) * (t - reference_c);
                            if e.abs() > want.abs() {
                                want = e;
                            }
                        }
                    }
                }
                let at_end = reading(frames.last().expect("frames"), "free strain");
                assert!(
                    (at_end - want).abs() < 1e-9,
                    "{name}: hot, the worst element wants {want:e} and has {at_end:e}"
                );

                // **The module is most stressed cold and relaxes as it heats**, which is the
                // answer a module assembled at reflow has to give: warming it back towards the
                // temperature it was built at relieves the assembly strain. Monotone across the
                // run, not merely lower at the end.
                let energy: Vec<f64> = frames.iter().map(|f| reading(f, "strain energy")).collect();
                for pair in energy.windows(2) {
                    assert!(
                        pair[1] <= pair[0] + 1e-12,
                        "{name}: heating towards the reference must relieve, not add: {energy:?}"
                    );
                }
                println!(
                    "  {name}: strain energy {:.4} J cold falling to {:.4} J hot, {:.1}x",
                    energy[0],
                    energy[energy.len() - 1],
                    energy[0] / energy[energy.len() - 1]
                );
                assert!(
                    energy[0] / energy[energy.len() - 1] > 5.0,
                    "{name}: it should relax by a lot: {:?} to {:?}",
                    energy[0],
                    energy[energy.len() - 1]
                );

                // **And the thermal answer is untouched by being watched.** The coupling is
                // one-way — nothing writes back into the block — so this run is still the same
                // physics scene `24` settles: climbing monotonically towards the 181.190771 C its
                // resistance stack fixes, and not past it. A structure that fed back would leave
                // that path, and the shape of the climb is a stronger statement than one number
                // because it holds at every frame rather than at the last.
                let junction: Vec<f64> = frames.iter().map(|f| reading(f, "peak")).collect();
                for pair in junction.windows(2) {
                    assert!(
                        pair[1] >= pair[0] - 1e-9,
                        "{name}: the junction only climbs while the module warms: {junction:?}"
                    );
                }
                let peak = *junction.last().expect("frames");
                assert!(
                    peak < 181.190_771 && peak > 179.0,
                    "{name}: 60 s is short of steady, so it sits just under scene 24's                      181.190771 C: {peak:.6} C"
                );
            }
            "26-poiseuille-in-a-cooling-channel.json" => {
                // The first scene to reach `dualis-fluid`, and the crate's own docs say why one
                // has to be written carefully: "it looks like a fluid" is the easiest wrong answer
                // in computational physics to accept, because a scheme with the wrong viscosity
                // still makes plausible vortices. So this is written around an exact solution.
                //
                // **And against the *discrete* parabola, not the continuum one.** A no-slip wall
                // imposed by reflecting the first cell makes the linear interpolation vanish at
                // the wall, and a parabola is not its own linear interpolation — so the settled
                // mean is `(gh²/12ν)(1 + 2/n²)`, exactly, and the `2/n²` is a closed form rather
                // than an error term. `dualis-fluid`'s own test derives it; this scene checks that
                // the *format* delivers the same problem the library was given.
                let (g, nu, cells) = (0.02, 1.004e-6, 16.0);
                let gap = cells * 0.125e-3;
                let continuum = dualis::fluid::poiseuille_mean_speed(g, gap, nu);
                let discrete = continuum * (1.0 + 2.0 / (cells * cells));

                let mean = frames
                    .last()
                    .expect("frames")
                    .readings
                    .iter()
                    .find(|r| r.label == "mean speed")
                    .expect("a channel reports its mean speed")
                    .value;

                // The tolerance is the **startup transient**, and it is predicted rather than
                // chosen. A channel from rest approaches its profile as `exp(−π²νt/h²)`, whose
                // time constant here is `h²/(π²ν)` = 0.405 s; the scene runs 4 s, which is 9.9 of
                // them, so what is left is `e^-9.9` = 5.0e-5 of the answer. Measured 4.8e-5.
                let constant = gap * gap / (std::f64::consts::PI.powi(2) * nu);
                let left = (-4.0f64 / constant).exp();
                println!(
                    "  {name}: mean {mean:.9} against a discrete {discrete:.9} m/s, off \
                     {:.2e}, with {left:.1e} of the startup left",
                    (mean / discrete - 1.0).abs()
                );
                assert!(
                    (mean / discrete - 1.0).abs() < 3.0 * left,
                    "{name}: {mean:.9} against the discrete parabola's {discrete:.9} m/s"
                );

                // **And it is nearer the discrete answer than the continuum one**, which is the
                // half that says the `2/n²` is real. They differ by 0.78% here and the run agrees
                // with one of them to 5e-5.
                assert!(
                    (mean - discrete).abs() < (mean - continuum).abs() / 100.0,
                    "{name}: the discrete parabola is the one this scheme solves: {mean:.9}                      against {discrete:.9} and {continuum:.9}"
                );

                // **The flow is incompressible and the advection is resolved**, or the profile
                // above would be a coincidence of a scheme that had stopped meaning anything.
                let last = frames.last().expect("frames");
                let of = |label: &str| {
                    last.readings
                        .iter()
                        .find(|r| r.label == label)
                        .unwrap_or_else(|| panic!("{name}: no {label}"))
                        .value
                };
                assert!(
                    of("divergence").abs() < 1e-9,
                    "{name}: divergence {:e} is not zero",
                    of("divergence")
                );
                assert!(
                    of("cell Reynolds") < dualis::fluid::CELL_REYNOLDS_LIMIT,
                    "{name}: cell Reynolds {:.3} is over the limit",
                    of("cell Reynolds")
                );
                // The pump is doing work and the run says how much — a driven channel whose books
                // showed nothing arriving would be one the audit had been talked out of.
                assert!(
                    of("work driven in") > 0.0,
                    "{name}: the drive must have put energy in: {:e} J",
                    of("work driven in")
                );
            }
            other => panic!("{other} ships but nothing checks it; add a claim for it"),
        }
    }
}

/// **The lamp's colour changes how much of it becomes heat**, and that is the whole point of
/// carrying spectra around.
///
/// The mirror is aluminium-like: about 90.5% reflective at 380 nm and 96.8% at 700 nm, so it
/// absorbs roughly three times as much in the blue as in the red. A blackbody at 6500 K puts
/// far more of its visible output at the blue end than one at 2800 K, so the same hundred
/// watts leaves more heat behind.
///
/// A flat reflectance would make this ratio exactly 1 and the entire spectral apparatus —
/// `Spectrum`, `SpectralPower`, `SurfaceOptics::absorptance` — would be an expensive way to
/// multiply by a constant. Asserting the *difference between two colour temperatures* is what
/// makes this a test of the optics rather than of a number.
#[test]
fn a_hotter_lamp_leaves_more_heat_on_a_mirror_that_is_worse_in_the_blue() {
    use dualis_world::light::{aluminium_mirror, Light};

    let at = |k: f64| Light::new("lamp", 100.0, k, aluminium_mirror()).absorbed_fraction();
    let (warm, cool) = (at(2800.0), at(6500.0));

    // Both are small: a good mirror absorbs a few percent, which is exactly why a hundred
    // watts on one is a thermal problem and not a catastrophe.
    assert!(
        (0.02..0.12).contains(&warm) && (0.02..0.12).contains(&cool),
        "a mirror absorbs a few percent, got {warm:.4} and {cool:.4}"
    );
    assert!(
        cool > warm,
        "6500 K is bluer than 2800 K and the mirror is worse in the blue: \
         {cool:.4} against {warm:.4}"
    );
}

/// **A scene can set a tolerance per quantity, and a typo in one is refused.**
///
/// The kernel gained `conservation_tolerance_for` because one number meant the loosest quantity in
/// a simulation set what every other one was checked against. Reaching it from data is this
/// crate's job, and the interesting half is the refusal: a channel name is matched against the
/// kernel's constants rather than passed through, so a misspelling cannot quietly leave a quantity
/// at the default.
///
/// That is the same failure `aluminum` for `aluminium` produced once already in this format — a
/// one-character difference that turned off the check the library exists for.
#[test]
fn a_scene_can_set_a_tolerance_per_quantity() {
    let text = r#"{
      "title": "per quantity", "duration_s": 0.01, "frames": 2,
      "conservation_tolerance": 1e-9,
      "tolerance_for": { "momentum": 1e-6, "photons": 0.5 },
      "domains": [
        { "kind": "room", "name": "room", "width_m": 4.4, "height_m": 3.1,
          "cells_across": 21 }
      ]
    }"#;
    let scene: Scene = serde_json::from_str(text).expect("it parses");
    let world = World::build(scene).expect("it builds");
    let tol = world.simulation().tolerances();

    assert_eq!(tol.default_tolerance(), 1e-9);
    assert_eq!(
        tol.for_quantity("energy"),
        1e-9,
        "unnamed keeps the default"
    );
    assert_eq!(tol.for_quantity("momentum"), 1e-6);
    assert_eq!(tol.for_quantity("photons"), 0.5);

    // A name the kernel does not have is refused, and the message says what is known.
    let typo = text.replace("\"momentum\"", "\"momentom\"");
    let scene: Scene = serde_json::from_str(&typo).expect("the JSON is still valid JSON");
    let Err(err) = World::build(scene) else {
        panic!("a misspelt channel must not be ignored");
    };
    assert!(
        err.contains("momentom") && err.contains("momentum"),
        "the refusal should name both what was written and what is known: {err}"
    );

    // And a scene that says nothing is unchanged — the feature costs nothing to ignore.
    let plain = text.replace(
        "\"tolerance_for\": { \"momentum\": 1e-6, \"photons\": 0.5 },",
        "",
    );
    let scene: Scene = serde_json::from_str(&plain).expect("it parses");
    let world = World::build(scene).expect("it builds");
    assert_eq!(world.simulation().tolerances().overrides().count(), 0);
    assert_eq!(
        world.simulation().tolerances().for_quantity("momentum"),
        1e-9
    );
}

/// **A scene from a newer build is refused, not half-run.**
///
/// `deny_unknown_fields` already catches a key that was *added*. It cannot catch a key whose
/// **meaning changed** — same name, same type, different semantics — and that is the whole reason
/// a format carries a version.
///
/// The failure being prevented is one this format has already had in a smaller form: a key that
/// is not read leaves its field at a default, and the run proceeds quietly doing something other
/// than what the file says. Refusing forward rather than attempting a best-effort downgrade is
/// deliberate: a plausible run that is not the one written down is worse than no run.
#[test]
fn a_scene_from_the_future_is_refused() {
    let text = r#"{
      "format": 99, "title": "from a newer dualis", "duration_s": 0.01, "frames": 2,
      "domains": [
        { "kind": "room", "name": "room", "width_m": 4.4, "height_m": 3.1, "cells_across": 21 }
      ]
    }"#;
    let scene: Scene = serde_json::from_str(text).expect("it is still valid JSON");
    assert_eq!(scene.format, 99);
    let Err(err) = World::build(scene) else {
        panic!("a format this build cannot read must not run");
    };
    assert!(err.contains("99") && err.contains("upgrade"), "{err}");

    // Zero is not a version this format ever had, and is what an uninitialised field looks like.
    let zero = text.replace("\"format\": 99", "\"format\": 0");
    let scene: Scene = serde_json::from_str(&zero).unwrap();
    assert!(World::build(scene).is_err(), "0 is not a version");
}

/// **Every scene that ships loads, and absence of the key means version 1.**
///
/// The seventeen shipped files were written before the field existed and none of them has it.
/// Defaulting to 1 is what makes that true rather than a special case — and it is why the default
/// is *absence means the original format* rather than *absence means whatever is current*.
#[test]
fn the_shipped_scenes_are_format_one_by_omission() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes");
    let mut seen = 0;
    for entry in std::fs::read_dir(&dir).expect("the scenes are there") {
        let path = entry.expect("readable").path();
        if path.extension().and_then(|e| e.to_str()) != Some("json") {
            continue;
        }
        let text = std::fs::read_to_string(&path).expect("readable");
        assert!(
            !text.contains("\"format\""),
            "{}: written before the field existed",
            path.display()
        );
        let scene: Scene = serde_json::from_str(&text).expect("it parses");
        assert_eq!(scene.format, 1, "{}", path.display());
        scene.check_version().expect("format 1 is readable");
        seen += 1;
    }
    assert!(seen >= 17, "only {seen} scenes were checked");
}

/// **What this build writes, it can read back — including the version.**
///
/// The narrow promise the number carries: within one version, a file that loads today loads
/// tomorrow. A round trip is the cheapest test of it and the one that would catch a serialiser
/// that emitted a version its own reader refuses.
#[test]
fn what_it_writes_it_reads() {
    let text = r#"{
      "title": "round trip", "duration_s": 0.25, "frames": 4,
      "domains": [
        { "kind": "room", "name": "room", "width_m": 4.4, "height_m": 3.1, "cells_across": 21 }
      ]
    }"#;
    let scene: Scene = serde_json::from_str(text).expect("it parses");
    let written = serde_json::to_string(&scene).expect("it serialises");

    assert!(
        written.contains(&format!("\"format\":{}", dualis_world::FORMAT)),
        "a file this build writes states its version: {written}"
    );
    let back: Scene = serde_json::from_str(&written).expect("it reads its own output");
    assert_eq!(back.format, dualis_world::FORMAT);
    back.check_version().expect("its own output is readable");
    World::build(back).expect("and runnable");
}
