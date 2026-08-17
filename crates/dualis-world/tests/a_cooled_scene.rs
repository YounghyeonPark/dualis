//! A scene can say which faces of a block are cooled, and the block then has a steady state.
//!
//! The capability arrived in `dualis-thermal` first and no file could reach it, which is the
//! oldest shape in `FRICTION.md`: the library can do a thing and the consumer cannot ask for
//! it. Until this key existed, every three-dimensional thermal scene in this repository warmed
//! for as long as it ran and settled nowhere — honest for a pulse and no answer at all to the
//! question a designer asks, which is *what temperature does this run at*.

use dualis_world::{Scene, World};

fn scene(json: &str) -> Scene {
    serde_json::from_str(json).expect("the test scene parses")
}

/// A heater and a block, with whatever `cooling` array the test wants.
fn heated(cooling: &str, watts: f64, seconds: f64, frames: usize) -> String {
    format!(
        r#"{{
  "title": "a part with something to lose heat to",
  "duration_s": {seconds},
  "frames": {frames},
  "domains": [
    {{ "kind": "heater", "name": "element", "watts": {watts}, "reserve_j": 1e9 }},
    {{ "kind": "block", "name": "part", "cells": [4, 4, 4], "cell_mm": 10.0,
      "initial_c": 20.0, "material": "aluminium"{cooling} }}
  ]
}}"#
    )
}

/// The block's mean temperature at the end of a run.
fn final_mean(json: &str) -> f64 {
    let mut world = World::build(scene(json)).expect("the scene builds");
    let frames = world.run().expect("the scene runs");
    frames
        .last()
        .expect("frames were captured")
        .readings
        .iter()
        .find(|r| r.domain == "part" && r.label == "mean")
        .map(|r| r.value)
        .expect("the block reports a mean")
}

/// **A scene with no `cooling` is the insulated block every earlier file is** — it warms without
/// limit, and the guarantee that let this key be added without changing any existing scene.
#[test]
fn a_scene_without_cooling_warms_without_limit() {
    let short = final_mean(&heated("", 40.0, 200.0, 4));
    let long = final_mean(&heated("", 40.0, 400.0, 8));
    // Twice the time, twice the rise: nothing is leaving.
    let (a, b) = (short - 20.0, long - 20.0);
    assert!(
        (b / a - 2.0).abs() < 1e-6,
        "an insulated block's rise is linear in time: {a:.3} K then {b:.3} K"
    );
}

/// **A cooled scene reaches a steady state, and it is the one the energy balance names.**
///
/// At steady state the block sheds exactly what the heater delivers, so `Q = h_total·A·ΔT` — and
/// the closed form has to use the *total* film, convective **plus** radiative, because the
/// domain models both. Aluminium's emissivity is 0.09, so radiation is a small correction here
/// rather than the factor of two it is for a black surface; the check computes it rather than
/// assuming either.
#[test]
fn a_cooled_scene_settles_where_the_energy_balance_says() {
    // Six faces of a 40 mm cube: 0.0096 m², 96 cm², 16 cm² each.
    let cooling = r#",
      "cooling": [
        { "face": "x-min", "ambient_c": 20.0, "convection_w_per_m2_k": 25.0, "area_cm2": 16.0 },
        { "face": "x-max", "ambient_c": 20.0, "convection_w_per_m2_k": 25.0, "area_cm2": 16.0 },
        { "face": "y-min", "ambient_c": 20.0, "convection_w_per_m2_k": 25.0, "area_cm2": 16.0 },
        { "face": "y-max", "ambient_c": 20.0, "convection_w_per_m2_k": 25.0, "area_cm2": 16.0 },
        { "face": "z-min", "ambient_c": 20.0, "convection_w_per_m2_k": 25.0, "area_cm2": 16.0 },
        { "face": "z-max", "ambient_c": 20.0, "convection_w_per_m2_k": 25.0, "area_cm2": 16.0 }
      ]"#;

    // Long enough to be at steady state, and the constant is worth computing rather than
    // guessing: C = 6.4e-5 m^3 × 2700 × 896 = 155 J/K and hA = 25 × 0.0096 = 0.24 W/K, so
    // τ = 645 s. A first pass at this test used a thousand seconds — 1.5 τ — and read a block
    // still visibly climbing.
    let settled = final_mean(&heated(cooling, 8.0, 6000.0, 12));
    let hotter = final_mean(&heated(cooling, 8.0, 8000.0, 16));
    assert!(
        (hotter - settled).abs() < 0.5,
        "nine time constants in, it should have stopped moving: {settled:.2} then {hotter:.2}"
    );

    // The balance, with the radiative film evaluated at the temperature it settled at.
    let area = 6.0 * 16.0e-4;
    let emissivity = 0.09_f64;
    let sigma = 5.670_374_419e-8_f64;
    let ambient = 293.15_f64;
    let t = settled + 273.15;
    let shed =
        25.0 * area * (t - ambient) + emissivity * sigma * area * (t.powi(4) - ambient.powi(4));
    assert!(
        (shed - 8.0).abs() / 8.0 < 0.05,
        "at {settled:.2} C the block sheds {shed:.3} W and the heater delivers 8 W"
    );
}

/// **The same part exposed the same way is the same problem on a finer grid.**
///
/// `area_cm2` is the whole face's, so refining the block does not change what the file says it
/// exposes — which is what lets `verify`'s resolution sweep compare two grids of one problem
/// rather than two problems. A per-cell area would have had to be halved and this would drift.
#[test]
fn refining_a_cooled_block_keeps_the_same_problem() {
    let cooling = r#",
      "cooling": [
        { "face": "z-max", "ambient_c": 20.0, "convection_w_per_m2_k": 40.0, "area_cm2": 16.0 }
      ]"#;
    let coarse = scene(&heated(cooling, 4.0, 600.0, 6));
    let fine = coarse.refined().expect("a block refines");

    let mean = |s: Scene| {
        let mut world = World::build(s).expect("builds");
        let frames = world.run().expect("runs");
        frames
            .last()
            .unwrap()
            .readings
            .iter()
            .find(|r| r.domain == "part" && r.label == "mean")
            .map(|r| r.value)
            .unwrap()
    };
    let a = mean(coarse);
    let b = mean(fine);
    // Same physical problem, so the means agree to discretisation — a per-cell area would put
    // the finer grid at half the loss and this would be off by tens of kelvin.
    assert!(
        (a - b).abs() < 1.0,
        "the same part on two grids should agree: {a:.2} C against {b:.2} C"
    );
}

/// **A face cooled twice is refused**, because a face cannot lose heat to two different airs and
/// silently summing them would answer a question nobody asked.
#[test]
fn a_face_cooled_twice_is_refused() {
    let cooling = r#",
      "cooling": [
        { "face": "z-max", "ambient_c": 20.0, "convection_w_per_m2_k": 25.0, "area_cm2": 16.0 },
        { "face": "z-max", "ambient_c": 5.0, "convection_w_per_m2_k": 60.0, "area_cm2": 16.0 }
      ]"#;
    let why = World::build(scene(&heated(cooling, 4.0, 10.0, 2)))
        .err()
        .expect("one face, two airs");
    assert!(why.contains("already cooled"), "{why}");
}

/// **A face with no area, a negative film, and an unknown face name are all refused.** The first
/// is the dangerous one: it would run as an insulated block and answer a different question with
/// nothing saying so.
#[test]
fn a_cooling_entry_that_cools_nothing_is_refused() {
    let no_area = r#",
      "cooling": [
        { "face": "z-max", "ambient_c": 20.0, "convection_w_per_m2_k": 25.0, "area_cm2": 0.0 }
      ]"#;
    let why = World::build(scene(&heated(no_area, 4.0, 10.0, 2)))
        .err()
        .expect("no area");
    assert!(why.contains("loses nothing"), "{why}");

    let uphill = r#",
      "cooling": [
        { "face": "z-max", "ambient_c": 20.0, "convection_w_per_m2_k": -5.0, "area_cm2": 16.0 }
      ]"#;
    let why = World::build(scene(&heated(uphill, 4.0, 10.0, 2)))
        .err()
        .expect("negative film");
    assert!(why.contains("uphill"), "{why}");

    let nonsense: Result<Scene, _> = serde_json::from_str(&heated(
        r#",
      "cooling": [
        { "face": "top", "ambient_c": 20.0, "convection_w_per_m2_k": 25.0, "area_cm2": 16.0 }
      ]"#,
        4.0,
        10.0,
        2,
    ));
    assert!(
        nonsense.is_err(),
        "a face this format does not have must not parse"
    );
}
