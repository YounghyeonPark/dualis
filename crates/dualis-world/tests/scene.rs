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
