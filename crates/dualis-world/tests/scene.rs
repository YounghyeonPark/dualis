//! The app holds itself to the library's standard: every number checked against something
//! the code did not compute.

use dualis::prelude::*;
use dualis_world::{DomainSpec, Scene, World};

fn room_scene(duration_s: f64, frames: usize, cells: usize) -> Scene {
    serde_json::from_str(&format!(
        r#"{{
          "title": "test", "schedule": "multirate",
          "duration_s": {duration_s}, "frames": {frames},
          "domains": [{{ "kind": "room", "name": "room", "width_m": 4.4,
            "height_m": 3.1, "cells_across": {cells}, "mode": [1, 1],
            "amplitude_pa": 1.0 }}]
        }}"#
    ))
    .expect("the test scene parses")
}

/// A scene survives being written down and read back.
///
/// The point of a scene format is that it leaves the process. If a round trip loses a field,
/// every saved world is quietly a different world.
#[test]
fn a_scene_round_trips_through_json() {
    let scene = room_scene(0.01, 4, 41);
    let text = serde_json::to_string(&scene).unwrap();
    let back: Scene = serde_json::from_str(&text).unwrap();

    assert_eq!(back.title, scene.title);
    assert_eq!(back.frames, scene.frames);
    assert_eq!(back.domains.len(), 1);
    assert_eq!(back.domains[0].name(), "room");
    assert!((back.duration_s - scene.duration_s).abs() < 1e-15);
    // Serialising the deserialised copy must give the same bytes, which catches a field that
    // round trips into a *default* rather than into itself.
    assert_eq!(serde_json::to_string(&back).unwrap(), text);
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
        assert!((v - 20.0).abs() < 1e-9, "bar drifted to {v} C");
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
    let al = Substance::aluminium_6061();
    let volume = Volume::from_si(20e-3 * 100e-6);
    let capacity = al
        .heat_capacity(volume)
        .expect("aluminium has a specific heat");
    let expected_rise = 6.0 / capacity.to_si();
    let mean_rise = frames[frames.len() - 1].panels[0]
        .values()
        .iter()
        .sum::<f64>()
        / frames[0].panels[0].values().len() as f64
        - 20.0;
    assert!(
        (mean_rise / expected_rise - 1.0).abs() < 1e-6,
        "the bar should have risen {expected_rise:.4} K, got {mean_rise:.4}"
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
    let middle = profile[cells / 2] - 20.0;
    let edge = profile[0] - 20.0;
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
    let sampled_rise = profile.iter().sum::<f64>() / cells as f64 - 20.0;
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
        let peak = last.panels[0]
            .values()
            .iter()
            .fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(peak.is_finite(), "{name}: the field went to {peak}");

        match name.as_str() {
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
                let (middle, end) = (v[v.len() / 2] - 20.0, v[0] - 20.0);
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
                // And every satellite is still in orbit rather than having been flung out.
                let (_, _, bounds) = match &last.panels[0].data {
                    dualis_world::PanelData::Points {
                        positions,
                        values,
                        bounds,
                    } => (positions, values, bounds),
                    _ => panic!("{name}: an orbit is bodies, not a field"),
                };
                assert!(
                    bounds[2] > 2.0e7,
                    "{name}: the frame should hold the widest orbit"
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
            other => panic!("{other} ships but nothing checks it; add a claim for it"),
        }
    }
}
