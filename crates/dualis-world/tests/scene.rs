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
                    "{name}: the fed end should still be above the mean,                      peak {:.4} against mean {mean:.4}",
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
                let (even_ring, bad_ring) =
                    (read("even", "ring over core"), read("wall gap", "ring over core"));

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
                let measured = puck.flow_rate().to_si() / dualis::porous::Liquid::water().density.to_si();
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
