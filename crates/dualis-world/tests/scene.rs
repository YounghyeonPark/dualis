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
                    .values
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
    for v in &bar.values {
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
        .values
        .iter()
        .sum::<f64>()
        / frames[0].panels[0].values.len() as f64
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
