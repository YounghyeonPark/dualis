//! The stage is a statement, and every domain answers it one of three ways.
//!
//! Gravity used to be an assumption hardcoded into the bounce's constructor: no scene could
//! be written on the Moon, none in free fall, and `9.80665` appeared in no file anywhere a
//! reader could see it. `environment` makes the stage a statement — and a stated condition
//! admits exactly three answers: consumed (the bounce falls at the stated rate, checked
//! against the closed form), refused (an orbit cannot stand in a uniform field, named), or
//! dismissed with the measurement that earns it (a molecular fluid, reported as a build
//! note). The fourth answer, silence, is what these tests exist to make impossible.

use pantometry_world::{PanelData, Scene, World};

fn scene(json: &str) -> Scene {
    serde_json::from_str(json).expect("the test scene parses")
}

// The tolerance is the integrator's, not a concession: `ContactSystem` marches semi-implicit
// Euler, whose energy on free flight oscillates at `O(dt)` — measured 1.3e-4 relative at this
// stiffness's substep — so 1e-3 is that with ~8x headroom, and the default 1e-6 refuses a
// correct run.
fn bounce(env: &str) -> String {
    format!(
        r#"{{
  "title": "a ball over a floor",
  "duration_s": 0.35,
  "frames": 2,
  "conservation_tolerance": 1e-3,
  {env}
  "domains": [
    {{ "kind": "bounce", "name": "ball", "drop_m": 1.0, "mass_kg": 1.0,
      "stiffness": 1e4, "damping": 0.0 }}
  ]
}}"#
    )
}

/// The ball's height in the last frame.
fn final_z(scene_json: &str) -> f64 {
    let mut world = World::build(scene(scene_json)).expect("the scene builds");
    let frames = world.run().expect("the scene runs");
    let last = frames.last().expect("frames were captured");
    let panel = last
        .panels
        .iter()
        .find(|p| p.name == "ball")
        .expect("the ball has a panel");
    match &panel.data {
        PanelData::Points { positions, .. } => positions[0][2],
        other => panic!("a ball is bodies, got {other:?}"),
    }
}

/// **An absent environment means what every existing scene means** — the old hardcoded
/// standard gravity, to the bit. This equality is what let the key be added without a format
/// bump: stating the assumption changed nothing, only made it statable.
#[test]
fn an_absent_environment_is_the_old_assumption_to_the_bit() {
    let assumed = final_z(&bounce(""));
    let stated = final_z(&bounce(
        r#""environment": { "gravity_m_per_s2": 9.80665 },"#,
    ));
    assert_eq!(
        assumed, stated,
        "stating the default must not move a single bit"
    );
}

/// **The stated gravity is the gravity it falls at, checked two ways against closed forms.**
///
/// The scheme is in the way of the naive check, instructively: `ContactSystem` marches
/// semi-implicit Euler, whose free-flight position is `h − ½g(t² + t·dt)` — a first-order
/// term the substep sets, so `z` does not land on `h − ½gt²` to rounding and a test
/// asserting it would be asserting the wrong closed form. Two forms survive the scheme:
///
/// - **The ratio.** Both the `½gt²` term and the `½gt·dt` error scale linearly in `g`, and
///   the substep is set by the contact stiffness, not by gravity — so two stages differing
///   only in `g` drop in exactly the ratio of their gravities, scheme error included. That
///   is checked to 1e-12, which rounding earns.
/// - **The band.** The error term is bounded by `½·g·t·dt` with `dt` at the contact limit
///   (~2e-4 s at this stiffness), giving ≤ 1e-4 m on the Moon drop; 1e-3 m is that with an
///   order of headroom, and tight enough that a wrong `g` — the next candidate up or down —
///   misses it by two orders.
#[test]
fn a_stated_gravity_is_the_gravity_it_falls_at() {
    let earth = final_z(&bounce(
        r#""environment": { "gravity_m_per_s2": 9.80665 },"#,
    ));
    let moon = final_z(&bounce(r#""environment": { "gravity_m_per_s2": 1.62 },"#));

    let ratio = (1.0 - moon) / (1.0 - earth);
    let exact_ratio = 1.62 / 9.80665;
    assert!(
        (ratio - exact_ratio).abs() < 1e-12,
        "two stages must drop in the ratio of their gravities: {ratio} against {exact_ratio}"
    );

    let exact = 1.0 - 0.5 * 1.62 * 0.35 * 0.35;
    assert!(
        (moon - exact).abs() < 1e-3,
        "the Moon drop should land within the scheme's own band of the closed form: {moon} \
         against {exact}"
    );

    // And zero is a legitimate stage: a drop tower's ball floats on its dashpot, exactly —
    // no field, no motion, no scheme error to have.
    let floating = final_z(&bounce(r#""environment": { "gravity_m_per_s2": 0.0 },"#));
    assert_eq!(floating, 1.0, "free fall means the ball does not move");
}

/// **An orbit refuses a stated uniform field, by name.** Its bodies fall only toward each
/// other; building it under stated gravity would run a free-fall experiment while the file
/// says Earth, which is the wrong-experiment failure this format exists to make impossible.
#[test]
fn an_orbit_refuses_a_stated_uniform_gravity() {
    let orbit = |env: &str| {
        format!(
            r#"{{
  "title": "an orbit on a stage",
  "duration_s": 4.0,
  "frames": 2,
  {env}
  "domains": [
    {{ "kind": "orbit", "name": "orbit", "central_kg": 5e24,
      "radii_m": [7e6], "satellite_kg": 1000.0 }}
  ]
}}"#
        )
    };
    let why = World::build(scene(&orbit(
        r#""environment": { "gravity_m_per_s2": 9.80665 },"#,
    )))
    .err()
    .expect("an orbit under uniform gravity must refuse");
    assert!(
        why.contains("orbit") && why.contains("cannot stand"),
        "{why}"
    );

    // A stated free-fall stage is exactly what an orbit is in, and builds.
    World::build(scene(&orbit(
        r#""environment": { "gravity_m_per_s2": 0.0 },"#,
    )))
    .expect("free fall is an orbit's own stage");
}

/// **A molecular fluid dismisses the stated gravity, and the dismissal is readable.** The
/// dismissal is correct physics — `m·g·σ` across one molecular diameter is ~1.3e-13 of the
/// Lennard-Jones well for argon — and the point of the note is that "correctly ignored" is
/// something a person reads rather than something that silently happened.
#[test]
fn a_molecular_fluid_dismisses_gravity_in_writing() {
    let atoms = scene(
        r#"{
  "title": "atoms on a stage",
  "duration_s": 0.5,
  "frames": 2,
  "environment": { "gravity_m_per_s2": 9.80665 },
  "domains": [
    { "kind": "atoms", "name": "atoms", "cells": 2, "density": 0.8,
      "temperature": 1.5, "seed": 7 }
  ]
}"#,
    );
    let world = World::build(atoms).expect("the dismissal is not a refusal");
    assert_eq!(world.notes().len(), 1);
    assert!(
        world.notes()[0].contains("dismissed") && world.notes()[0].contains("1.3e-13"),
        "{:?}",
        world.notes()
    );
}

/// **A stage that does not describe one is refused** — negative gravity is a direction flip
/// nobody means, and an unknown key is a typo that would otherwise fall back to a default and
/// run something other than what the file says.
#[test]
fn a_nonsense_stage_is_refused() {
    let why = World::build(scene(&bounce(
        r#""environment": { "gravity_m_per_s2": -9.8 },"#,
    )))
    .err()
    .expect("negative gravity must refuse");
    assert!(why.contains("magnitude"), "{why}");

    let unknown: Result<Scene, _> = serde_json::from_str(&bounce(
        r#""environment": { "gravity_m_per_s2": 9.8, "wind": 3.0 },"#,
    ));
    assert!(unknown.is_err(), "an unknown stage key must not parse");
}
