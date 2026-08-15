//! The `verify` battery, held to its own standard.
//!
//! The battery's whole claim is that it measures rather than asserts, so its tests are the
//! places a measurement has a known answer: a refinement that must preserve a total exactly, a
//! window error whose order the kernel documentation derives and measures, a mode whose
//! convergence rate is the scheme's own second order, and a hazard the schedule documentation
//! states in words that a test can arrange in numbers.

use dualis_world::verify::{verify, Order, SweepOutcome};
use dualis_world::{DomainSpec, Scene};

fn scene(json: &str) -> Scene {
    serde_json::from_str(json).expect("the test scene parses")
}

/// The sweep behind an outcome, unwrapped where the test knows it must have run.
fn ran(outcome: &SweepOutcome) -> &dualis_world::verify::Sweep {
    match outcome {
        SweepOutcome::Ran(s) => s,
        SweepOutcome::Skipped(why) => panic!("the sweep was skipped: {why}"),
        SweepOutcome::Failed(why) => panic!("the sweep failed: {why}"),
    }
}

/// The order a sweep measured for one reading, unwrapped where the test knows it must exist.
fn order_of(
    battery: &dualis_world::verify::Battery,
    which: &str,
    domain: &str,
    label: &str,
) -> f64 {
    let sweep = ran(match which {
        "window" => &battery.window,
        _ => &battery.resolution,
    });
    let (_, _, order) = sweep
        .orders
        .iter()
        .find(|(d, l, _)| d == domain && l == label)
        .expect("the reading was matched across the runs");
    match order {
        Order::Measured(p) => *p,
        Order::BelowFloor => {
            panic!("{domain}/{label}: at the floor, where this test expected a measurable order")
        }
        Order::NotAsymptotic => {
            panic!("{domain}/{label}: not asymptotic, where this test expected a measurable order")
        }
    }
}

/// **Two runs of one scene produce identical bytes, and the battery checks it rather than
/// assuming it.** This is the determinism promise at the scale of a whole run — every panel,
/// body and reading of every frame — not just the pinned generator stream.
#[test]
fn two_runs_are_bit_identical_and_the_battery_says_so() {
    let s = scene(
        r#"{
  "title": "a small room, checked",
  "duration_s": 0.005,
  "frames": 3,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 2.0, "height_m": 1.0,
      "cells_across": 21,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } }
  ]
}"#,
    );
    let b = verify(&s, false).expect("the battery runs");
    assert!(b.deterministic, "one scene, two different runs");
    assert!(
        b.findings.is_empty(),
        "findings on a clean scene: {:?}",
        b.findings
    );
    // Both sweeps ran and actually compared something — a sweep with zero shifts would be
    // indistinguishable in the report's shape from "swept and nothing moved".
    assert!(!ran(&b.window).shifts.is_empty());
    assert!(!ran(&b.resolution).shifts.is_empty());
    // The margin rows exist: energy was on the ledger and the room reported a limit.
    assert!(b.base.drift.iter().any(|(q, _, _)| q == "energy"));
    assert!(b.base.stability.iter().any(|(n, _, _)| n == "room"));
}

/// **A scene that reports no readings is a finding, not a verification.**
///
/// Neither mechanics domain implements `Domain::readings`, so an orbit scene gives the sweeps
/// nothing to compare. The battery used to run the sweeps anyway, compare zero numbers, print
/// a heading with nothing under it and exit 0 — a report whose shape is identical to "swept
/// and nothing moved". A user would conclude the orbit is insensitive to the coupling window;
/// in fact nothing was measured. Now: both sweeps are skipped with the reason stated, the
/// report says which domain was not swept, and the empty measurement is a finding that
/// carries the exit code.
#[test]
fn a_scene_with_no_readings_is_a_finding_not_a_verification() {
    let s = scene(
        r#"{
  "title": "an orbit nobody reads",
  "duration_s": 4.0,
  "frames": 4,
  "domains": [
    { "kind": "orbit", "name": "orbit", "central_kg": 5e24,
      "radii_m": [7e6], "satellite_kg": 1000.0 }
  ]
}"#,
    );
    let b = verify(&s, false).expect("the run itself is fine");
    assert!(
        b.findings
            .iter()
            .any(|f| f.contains("nothing here was verified")),
        "an empty measurement must be a finding: {:?}",
        b.findings
    );
    assert!(matches!(b.window, SweepOutcome::Skipped(_)));
    assert!(b.unread.contains(&"orbit".to_string()));
    let report = b.render();
    assert!(
        report.contains("SKIPPED: no domain reports a reading"),
        "the report does not say what was not measured:\n{report}"
    );
}

/// **A sweep run the audit refused is a FAILED finding, not a benign skip.**
///
/// A staggered room just under its CFL limit at base resolution is over it at 2×, because
/// refining halves the limit. The refined run is refused on its first advance — and a scene
/// that passes at its own settings and fails when a knob moves is one knob away from wrong,
/// which is about the most alarming fact this battery can learn. It used to render as
/// `SKIPPED:` and exit 0, the same reading as "not applicable".
#[test]
fn a_sweep_run_the_audit_refused_is_a_finding_not_a_skip() {
    let s = scene(
        r#"{
  "title": "stable at its own resolution and not at the next",
  "schedule": "staggered",
  "duration_s": 0.00045,
  "frames": 3,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.0, "height_m": 2.0,
      "cells_across": 41,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } }
  ]
}"#,
    );
    // Base window 1.5e-4 s sits under the 41-sample limit of 2.06e-4 s, so the base run and
    // the hazard check both pass; at 81 samples the limit is 1.03e-4 s and the refined run is
    // refused.
    let b = verify(&s, false).expect("the base scene runs");
    assert!(matches!(b.resolution, SweepOutcome::Failed(_)));
    assert!(
        b.findings
            .iter()
            .any(|f| f.contains("one knob away from wrong")),
        "a refused sweep run must be a finding: {:?}",
        b.findings
    );
    let report = b.render();
    assert!(
        report.contains("FAILED:"),
        "not rendered as a failure:\n{report}"
    );
}

/// **The report states what it did not measure, row by row.**
///
/// A mixed scene: the room is swept and the orbit — which reports no readings — must appear
/// as "no readings — not swept" rather than being simply absent, because a sweep section
/// listing every other domain looks complete. And a flat orbit's `momentum_z` never rises
/// above the audit's scale floor, so it is never audited; the margins table must say that
/// rather than omit the row, because a table of quantities-with-margin that quietly drops one
/// is the shape of a check that turned itself off.
#[test]
fn the_report_states_what_it_did_not_measure() {
    let s = scene(
        r#"{
  "title": "a room and an unread orbit",
  "duration_s": 0.005,
  "frames": 3,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 2.0, "height_m": 1.0,
      "cells_across": 21,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } },
    { "kind": "orbit", "name": "orbit", "central_kg": 5e24,
      "radii_m": [7e6], "satellite_kg": 1000.0 }
  ]
}"#,
    );
    let b = verify(&s, false).expect("the battery runs");
    let report = b.render();
    assert!(
        report.contains("orbit          no readings — not swept"),
        "the unread domain is invisible:\n{report}"
    );
    assert!(
        report.contains("never above the scale floor — not audited"),
        "the unaudited quantity is invisible:\n{report}"
    );
    // And the room was genuinely swept beside it.
    assert!(!ran(&b.window).shifts.is_empty());
}

/// **A staggered window past a domain's stability limit is named as the likely cause.**
///
/// `ScheduleSpec::Staggered`'s own documentation: a frame interval larger than the limit is
/// silently unstable. The room turns out not to be silent — its step reports the Courant
/// number as a created quantity and the audit refuses the very first advance — but a refusal
/// that only says "Courant number created" leaves the *why* to be worked out. The battery
/// checks the hazard on the built world before running anything, so its error carries both
/// halves: what the kernel refused, and the schedule choice that made it inevitable.
#[test]
fn a_staggered_window_past_the_limit_is_named_as_the_likely_cause() {
    let unstable = r#"{
  "title": "one step across the whole window",
  "schedule": "staggered",
  "duration_s": 0.002,
  "frames": 3,
  "conservation_tolerance": 1e6,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.0, "height_m": 2.0,
      "cells_across": 41,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } }
  ]
}"#;
    let why = verify(&scene(unstable), false).expect_err("the room refuses an over-CFL step");
    assert!(
        why.contains("likely why") && why.contains("silently unstable"),
        "the error does not name the schedule hazard: {why}"
    );

    // The same room under multirate subcycles, and the same battery finds nothing.
    let subcycled = unstable.replace("\"staggered\"", "\"multirate\"");
    let b = verify(&scene(&subcycled), false).expect("multirate subcycles the same window");
    assert!(
        !b.findings.iter().any(|f| f.contains("silently unstable")),
        "multirate flagged as unstable: {:?}",
        b.findings
    );
}

/// **A refinement preserves the problem, exactly where exactness is checkable.**
///
/// A bar doubles its cells and halves each face; the beam that lands on it doubles its faces
/// and halves each area. The exposed boundary's total area — `cells × face_area` — is the
/// number both sides must preserve to the bit, because the kernel refuses a face-count
/// mismatch and conserves the flux across whatever count both sides agree on.
#[test]
fn a_bar_and_its_beam_refine_together() {
    let s = scene(
        r#"{
  "title": "a beam on a bar",
  "duration_s": 1.0,
  "frames": 2,
  "domains": [
    { "kind": "beam", "name": "beam", "onto": "face", "faces": 10,
      "face_area_mm2": 2.0, "watts": 5.0, "reserve_j": 100.0, "waist_fraction": 0.2 },
    { "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 10, "area_mm2": 4.0,
      "initial_c": 20.0, "exposes": { "name": "face", "face_area_mm2": 2.0 } }
  ]
}"#,
    );
    let fine = s.refined().expect("a bar and a beam both refine");
    let mut checked = 0;
    for spec in &fine.domains {
        match spec {
            DomainSpec::Beam {
                faces,
                face_area_mm2,
                watts,
                ..
            } => {
                assert_eq!(*faces, 20);
                assert_eq!(*face_area_mm2, 1.0);
                assert_eq!(*watts, 5.0, "refinement moved a physical parameter");
                checked += 1;
            }
            DomainSpec::Bar { cells, exposes, .. } => {
                assert_eq!(*cells, 20);
                let b = exposes.as_ref().expect("the boundary survives refinement");
                assert_eq!(b.face_area_mm2, 1.0);
                checked += 1;
            }
            other => panic!("an unexpected domain appeared: {other:?}"),
        }
    }
    assert_eq!(checked, 2);
}

/// **A room refines by intervals, not by count.** Its samples sit on the walls, so `n` samples
/// are `n − 1` intervals and halving the spacing is `2n − 1`. `2n` would be a spacing 60/121
/// of the old one — near half, not half — and the residue would be read as physics.
#[test]
fn a_room_refines_by_intervals() {
    let s = scene(
        r#"{
  "title": "a room",
  "duration_s": 0.01,
  "frames": 2,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.0, "height_m": 2.0,
      "cells_across": 41 }
  ]
}"#,
    );
    let fine = s.refined().expect("a room refines");
    match &fine.domains[0] {
        DomainSpec::Room { cells_across, .. } => assert_eq!(*cells_across, 81),
        other => panic!("{other:?}"),
    }
}

/// **What cannot be refined honestly refuses, naming the feature.** A hot spot is one cell of
/// excess temperature: at half the cell it holds an eighth of the joules, so a "refined" scene
/// would be a different experiment whose difference the sweep would report as discretisation
/// error. The same for a conductor's one-cell notch and a puck's one-cell channel ring.
#[test]
fn one_cell_features_refuse_to_refine() {
    let hot = scene(
        r#"{
  "title": "a warmed cell",
  "duration_s": 1.0,
  "frames": 2,
  "domains": [
    { "kind": "block", "name": "block", "cells": [4, 4, 4], "cell_mm": 2.0,
      "initial_c": 20.0, "hot_spot": { "at": [2, 2, 2], "above_k": 10.0 } }
  ]
}"#,
    );
    let why = hot.refined().expect_err("a hot spot must refuse");
    assert!(
        why.contains("hot spot"),
        "the refusal does not say why: {why}"
    );

    let notched = scene(
        r#"{
  "title": "a notched conductor",
  "duration_s": 1.0,
  "frames": 2,
  "domains": [
    { "kind": "conductor", "name": "busbar", "cells": [8, 2, 2], "cell_mm": 1.0,
      "resistivity_ohm_m": 1.724e-8, "volts": 0.01, "blocked": [[4, 0, 0]] }
  ]
}"#,
    );
    let why = notched.refined().expect_err("a blocked cell must refuse");
    assert!(why.contains("notch"), "{why}");

    let channelled = scene(
        r#"{
  "title": "a channelled puck",
  "duration_s": 5.0,
  "frames": 2,
  "domains": [
    { "kind": "puck", "name": "puck", "cells": [9, 6, 9], "cell_mm": 2.0,
      "radius_mm": 7.0, "grind_um": 250.0, "porosity": 0.45, "bar": 9.0,
      "brew_c": 93.0, "channel_porosity": 0.6 }
  ]
}"#,
    );
    let why = channelled
        .refined()
        .expect_err("a channel ring must refuse");
    assert!(why.contains("channel"), "{why}");
}

/// **A scene of pure systems has no finer statement, and says so instead of "verifying" a
/// no-op.** Refining nothing and comparing a run against itself would print a zero shift and
/// look like the cleanest convergence in the file — the empty-panel failure shape, one level
/// up.
#[test]
fn a_scene_of_pure_systems_has_nothing_to_refine() {
    let s = scene(
        r#"{
  "title": "a heater warming a lump",
  "duration_s": 10.0,
  "frames": 5,
  "domains": [
    { "kind": "heater", "name": "element", "watts": 5.0, "reserve_j": 1000.0 },
    { "kind": "lump", "name": "plate", "volume_cm3": 0.5, "thickness_mm": 2.0,
      "initial_c": 20.0, "ambient_c": 20.0, "area_cm2": 400.0 }
  ]
}"#,
    );
    let why = s.refined().expect_err("nothing here is a discretisation");
    assert!(
        why.contains("nothing in this scene is a discretisation"),
        "{why}"
    );
}

/// **The window error of the once-per-frame feedback is first order, and the battery sees it.**
///
/// The `tracks` loop is closed between frames because it cannot be closed inside the step
/// loop, and its own documentation states the price: the winding's temperature lags by up to
/// one frame interval. A resistance evaluated at a stale temperature dissipates the wrong
/// watts for the whole window, so the accumulated error is first order in the window — the
/// error class the conservation audit structurally cannot see, since every joule the winding
/// spends does arrive. The battery's window sweep must read an order near 1 here; reading 2
/// would mean it cannot see the one error class it most exists to see. (The other candidate,
/// a subcycling consumer taking a deposited interval, stopped being first order when
/// `Exchange::take_share` fixed it — measured by `multirate_timing.rs`, and measured again
/// here in the negative: this scene is the one with a first-order knob left.)
///
/// The band is earned by measurement against the theory's limit, **and it holds only at this
/// pin**. First order is the *limit* as the window shrinks: at 20 s windows against a ~90 s
/// time constant the order read 0.88, at 10 s it read 1.37 — outside this band, one octave
/// from the asymptotic regime — and at the 5 s windows this test uses it reads 1.01, with the
/// instantaneous readings at 1.00 confirming the regime. So `[0.8, 1.2]` around a measured
/// 1.01: the width is against benign code changes moving the arithmetic path, not against
/// platforms — results are bit-identical across platforms by convention 3, so that costs no
/// width at all — and it is narrow enough that second order cannot pass.
#[test]
fn the_window_error_of_the_tracks_feedback_is_first_order() {
    let s = scene(
        r#"{
  "title": "a winding tracking its own heat",
  "duration_s": 60.0,
  "frames": 12,
  "domains": [
    { "kind": "winding", "name": "coil", "length_m": 10.0, "cross_section_mm2": 0.1,
      "amps": 3.0, "at_c": 20.0, "reserve_j": 40000.0, "tracks": "motor/coil" },
    { "kind": "network", "name": "motor",
      "nodes": [ { "name": "coil", "material": "copper", "volume_cm3": 2.0,
                   "thickness_mm": 2.0, "initial_c": 20.0,
                   "loses_to": { "ambient_c": 20.0, "area_cm2": 100.0 } } ],
      "links": [],
      "absorbing": "coil" }
  ]
}"#,
    );
    let b = verify(&s, true).expect("the battery runs");
    let p = order_of(&b, "window", "coil", "spent");
    assert!(
        (0.8..=1.2).contains(&p),
        "the tracks-feedback error should read first order, measured {p:.3}"
    );
}

/// **A room mode with commensurate dimensions converges at the scheme's second order.**
///
/// Width 4 m at 0.1 m spacing and height 2 m — exactly 20 cells — stay commensurate through
/// both refinements, so the quantised geometry does not move and what remains is the scheme:
/// second order in space, second order in time, first-order startup fixed. The measured order
/// of the peak reading should sit near 2. Under the incommensurate default scene it reads
/// near 1, because the *height itself* converges at first order — measured and then confirmed
/// against `f(1,1)` by hand — which is why this test pins the commensurate case: it is the one
/// with a clean known answer.
///
/// Measured: 2.24. Above 2 rather than below because the substep count is a ceiling, so the
/// effective time step shrinks slightly faster than the spacing; 2.4 is that ceiling argument
/// with room, and 1.6 is where a first-order term stops being subdominant. What the band can
/// and cannot refuse is a property of a three-point measurement worth stating exactly: with
/// error `A·h + B·h²`, `p ≥ 1.6` admits a first-order term contributing up to ~40% of the
/// coarse-grid difference — but that term *doubles* its share per refinement, so it is caught
/// one octave later, and a **dominant** first-order defect (both shipped acoustic defects
/// were) reads near 1 and fails here outright.
#[test]
fn a_commensurate_room_mode_converges_at_second_order() {
    let s = scene(
        r#"{
  "title": "a commensurate room",
  "duration_s": 0.02,
  "frames": 5,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.0, "height_m": 2.0,
      "cells_across": 41,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } }
  ]
}"#,
    );
    let b = verify(&s, true).expect("the battery runs");
    let p = order_of(&b, "resolution", "room", "peak");
    assert!(
        (1.6..=2.4).contains(&p),
        "a second-order scheme on unmoving geometry, measured {p:.3}"
    );
}
