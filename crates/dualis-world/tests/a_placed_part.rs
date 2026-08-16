//! A scene can say where a part sits — and cannot say that two parts touch.
//!
//! `Pose` has been in the kernel since placement arrived, and no file could reach it: every
//! domain sat at the origin because the format had nowhere to say otherwise. `poses` is that
//! place. What it moves is where a domain is captured and drawn, through the kernel's rigid
//! motion, so distances and angles survive exactly.
//!
//! What it deliberately does not do is make placed parts interact, and the last two tests pin
//! that rather than leaving it to be discovered. Domains meet on the bus and over an
//! `Interface` matched by face count; no geometry is consulted anywhere, so parts drawn in
//! contact exchange nothing. `ARCHITECTURE.md` names that as the gap that arrives first in
//! practice, and placing is the half that can be done without answering the other half.
//!
//! Writing those tests turned up a second obstacle, larger and not in that list: **two domains
//! that consume the same channel cannot share a scene at all.** Two blocks with no heater
//! anywhere are refused on the second one's first step, because the kernel's guard counts takes
//! rather than amounts. Today a scene holds one heat-consuming domain, which is a harder limit
//! on assembly than contact is.

use dualis_world::{Scene, World};

fn scene(json: &str) -> Scene {
    serde_json::from_str(json).expect("the test scene parses")
}

/// A block, with whatever `poses` block the test wants beside it.
fn block(poses: &str) -> String {
    format!(
        r#"{{
  "title": "a placed block",
  "duration_s": 1.0,
  "frames": 1,
  {poses}
  "domains": [
    {{ "kind": "block", "name": "block", "cells": [2, 2, 2], "cell_mm": 10.0,
      "initial_c": 20.0 }}
  ]
}}"#
    )
}

/// The world-space corners of the one domain's extent, low corner first.
fn extent_of(json: &str) -> ([f64; 3], [f64; 3]) {
    let s = scene(json);
    let placed = s.placements();
    let p = placed.get("block").expect("the block is placed");
    let e = p.extent.expect("a block is a field with an extent");
    let lo = p.pose.point_to_world(e.min).to_si();
    let hi = p.pose.point_to_world(e.max).to_si();
    ([lo.x, lo.y, lo.z], [hi.x, hi.y, hi.z])
}

/// **An absent `poses` block leaves every domain exactly where it was.** The guarantee that let
/// this key be added without a format bump: no existing file changes meaning.
#[test]
fn a_scene_without_poses_is_unmoved() {
    let (lo, hi) = extent_of(&block(""));
    assert_eq!(lo, [0.0, 0.0, 0.0]);
    // Two cells of 10 mm on each axis.
    assert_eq!(hi, [0.02, 0.02, 0.02]);
}

/// **A stated position moves the part there, exactly.** A translation is exact arithmetic —
/// there is no scheme in the way — so this is an equality and not a tolerance.
#[test]
fn a_stated_position_moves_the_part() {
    let (lo, hi) = extent_of(&block(
        r#""poses": { "block": { "at_m": [0.5, -0.25, 1.0] } },"#,
    ));
    assert_eq!(lo, [0.5, -0.25, 1.0]);
    assert_eq!(hi, [0.52, -0.23, 1.02]);
}

/// **A turn is a rigid motion: it moves the part without changing its size.**
///
/// Ninety degrees about z carries `+x` onto `+y`, so a 20 mm span along x arrives as a 20 mm
/// span along y. The tolerance is the quaternion's rounding on a right angle — `sin(π/4)` is
/// irrational, so the corner does not land on an exact decimal — and 1e-12 is orders above the
/// few ulps that costs.
#[test]
fn a_turn_moves_the_part_without_resizing_it() {
    let (lo, hi) = extent_of(&block(
        r#""poses": { "block": { "at_m": [0.0, 0.0, 0.0],
                                "turn": { "axis": [0.0, 0.0, 1.0], "degrees": 90.0 } } },"#,
    ));
    // The low corner is the origin and stays put; the high corner's x span becomes a y span.
    assert!(lo.iter().all(|v| v.abs() < 1e-12), "{lo:?}");
    assert!(
        (hi[0] + 0.02).abs() < 1e-12
            && (hi[1] - 0.02).abs() < 1e-12
            && (hi[2] - 0.02).abs() < 1e-12,
        "a right angle about z should carry +x onto +y, got {hi:?}"
    );

    // And the diagonal — the thing a rigid motion must preserve — is the same length as unposed.
    let (lo0, hi0) = extent_of(&block(""));
    let span = |a: [f64; 3], b: [f64; 3]| {
        ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2) + (b[2] - a[2]).powi(2)).sqrt()
    };
    assert!(
        (span(lo, hi) - span(lo0, hi0)).abs() < 1e-12,
        "a turn changed the part's size"
    );
}

/// **A pose naming a domain the scene does not define is refused**, with the file's own
/// vocabulary quoted back. It would otherwise place nothing, silently — the same shape as a
/// `tracks` pointing at a node that is not there.
#[test]
fn a_pose_for_a_domain_that_is_not_there_is_refused() {
    let why = World::build(scene(&block(
        r#""poses": { "blcok": { "at_m": [0.0, 0.0, 0.0] } },"#,
    )))
    .err()
    .expect("a misspelled domain must refuse");
    assert!(why.contains("blcok") && why.contains("block"), "{why}");
}

/// **A turn with no direction is refused rather than normalised into a NaN**, and an unknown key
/// does not parse — `deny_unknown_fields` holds here as it does everywhere else in this format,
/// which is why `poses` is a map rather than a flattened field.
#[test]
fn a_degenerate_pose_is_refused() {
    let why = World::build(scene(&block(
        r#""poses": { "block": { "at_m": [0.0, 0.0, 0.0],
                                "turn": { "axis": [0.0, 0.0, 0.0], "degrees": 45.0 } } },"#,
    )))
    .err()
    .expect("a zero axis must refuse");
    assert!(why.contains("no direction"), "{why}");

    let unknown: Result<Scene, _> = serde_json::from_str(&block(
        r#""poses": { "block": { "at_m": [0.0, 0.0, 0.0], "scale": 2.0 } },"#,
    ));
    assert!(
        unknown.is_err(),
        "an unknown pose key must not parse — a scale is not a rigid motion"
    );
}

/// **Two parts that consume the same channel cannot share a scene at all** — which is a larger
/// obstacle to assembly than the one this file was written to record.
///
/// Two blocks, no heater, nothing published anywhere. The kernel refuses on the second block's
/// first step: `Exchange::take` empties a channel, so a second consumer would silently get
/// nothing and warm at the wrong rate, and the kernel would rather say so than let two plates
/// under one lamp warm as one. The guard counts **takes**, not amounts, so it fires even when
/// the channel was empty and nothing could have been mis-split.
///
/// Recorded rather than worked around, because it is a fact somebody assembling parts meets
/// immediately: today a scene holds **one** heat-consuming domain. Whether the guard should look
/// at the amount rather than the count is a question about the kernel, and the kernel is not
/// this file's to change.
#[test]
fn two_parts_that_take_from_one_channel_are_refused() {
    let s = scene(
        r#"{
  "title": "two blocks, face to face",
  "duration_s": 10.0,
  "frames": 2,
  "poses": { "cold": { "at_m": [0.02, 0.0, 0.0] } },
  "domains": [
    { "kind": "block", "name": "hot", "cells": [2, 2, 2], "cell_mm": 10.0,
      "initial_c": 500.0 },
    { "kind": "block", "name": "cold", "cells": [2, 2, 2], "cell_mm": 10.0,
      "initial_c": 20.0 }
  ]
}"#,
    );

    // They really are face to face: the hot block ends where the cold one begins.
    let placed = s.placements();
    let hot = placed.get("hot").unwrap();
    let cold = placed.get("cold").unwrap();
    let hot_hi = hot.pose.point_to_world(hot.extent.unwrap().max).to_si();
    let cold_lo = cold.pose.point_to_world(cold.extent.unwrap().min).to_si();
    assert!(
        (hot_hi.x - cold_lo.x).abs() < 1e-12,
        "the blocks should share a face at x = {}, and the cold one starts at {}",
        hot_hi.x,
        cold_lo.x
    );

    // And the scene builds — the refusal is the kernel's, at the first step, not the builder's.
    let mut world = World::build(s).expect("two blocks describe a simulation");
    let violation = world.run().expect_err("two heat consumers must be refused");
    assert_eq!(violation.quantity, "energy");
    assert!(
        violation.site.contains("already emptied"),
        "the refusal should name the emptied channel: {}",
        violation.site
    );
}

/// **Two parts that can share a scene, placed in contact, exchange nothing.**
///
/// A block and a room, drawn face to face and coupled by nothing: domains meet on the bus and
/// over an `Interface` matched by face count, and no geometry is consulted anywhere. The block
/// ends exactly as hot as it started.
///
/// This is `ARCHITECTURE.md`'s "two parts have no way to touch", measured. When that gap closes
/// this test is the one that must change, and its failure will be the notification.
#[test]
fn two_placed_parts_are_drawn_touching_and_exchange_nothing() {
    let s = scene(
        r#"{
  "title": "a block against a room",
  "duration_s": 0.01,
  "frames": 2,
  "poses": { "room": { "at_m": [0.02, 0.0, 0.0] } },
  "domains": [
    { "kind": "block", "name": "block", "cells": [2, 2, 2], "cell_mm": 10.0,
      "initial_c": 500.0 },
    { "kind": "room", "name": "room", "width_m": 0.4, "height_m": 0.2,
      "cells_across": 9 }
  ]
}"#,
    );

    let placed = s.placements();
    let b = placed.get("block").unwrap();
    let r = placed.get("room").unwrap();
    let block_hi = b.pose.point_to_world(b.extent.unwrap().max).to_si();
    let room_lo = r.pose.point_to_world(r.extent.unwrap().min).to_si();
    assert!(
        (block_hi.x - room_lo.x).abs() < 1e-12,
        "the two should share a face at x = {}, and the room starts at {}",
        block_hi.x,
        room_lo.x
    );

    let mut world = World::build(s).expect("the scene builds");
    let frames = world.run().expect("and runs");
    let last = frames.last().unwrap();
    let mean = last
        .readings
        .iter()
        .find(|r| r.domain == "block" && r.label == "mean")
        .map(|r| r.value)
        .expect("the block reports a mean");
    assert!(
        (mean - 500.0).abs() < 1e-9,
        "the block changed temperature, which would mean geometry started coupling domains: {mean}"
    );
}
