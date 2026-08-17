//! A scene can name STL files and get an assembly: parts where their own coordinates put them,
//! nothing between them, and a rasterisation report for what the grid cost.
//!
//! The library could do this a commit ago and no file could ask for it — the oldest shape in
//! `FRICTION.md`. What a scene could state before was a box of one material with boxes of other
//! materials cut into it, which is a layered wall and not an assembly.

use dualis_world::{Scene, World};

fn scene(json: &str) -> Scene {
    serde_json::from_str(json).expect("the test scene parses")
}

/// Write a binary STL of an axis-aligned brick, in **millimetres**, and give back its path.
///
/// Millimetres because that is what an STL means: the format is unitless and every CAD tool
/// writes mm, so `Mesh::from_stl` scales by `1e-3`. Writing metres here made every brick a
/// thousand times too small, which voxelised to no cells at all — the first version of this
/// fixture did exactly that, and the error it produced ("every part voxelised to no cells") was
/// the right one for a reason it did not name.
///
/// Written by the test rather than committed as a fixture: the geometry under test should be
/// visible beside the assertion, not one binary file away from it. Binary rather than ASCII
/// because that is what a CAD tool exports.
fn brick_stl(dir: &std::path::Path, name: &str, low: [f32; 3], high: [f32; 3]) -> String {
    let v = |i: usize| {
        [
            if i & 1 == 0 { low[0] } else { high[0] },
            if i & 2 == 0 { low[1] } else { high[1] },
            if i & 4 == 0 { low[2] } else { high[2] },
        ]
    };
    // Twelve triangles, wound so each face's normal points out of the brick.
    let tris: [[usize; 3]; 12] = [
        [0, 4, 6],
        [0, 6, 2],
        [1, 3, 7],
        [1, 7, 5],
        [0, 1, 5],
        [0, 5, 4],
        [2, 6, 7],
        [2, 7, 3],
        [0, 2, 3],
        [0, 3, 1],
        [4, 5, 7],
        [4, 7, 6],
    ];

    let mut bytes = vec![0u8; 80];
    bytes.extend_from_slice(&(tris.len() as u32).to_le_bytes());
    for t in tris {
        // The normal is left at zero: STL stores one and every reader worth the name recomputes
        // it from the winding, which is the only copy that cannot disagree with the vertices.
        bytes.extend_from_slice(&[0u8; 12]);
        for idx in t {
            for component in v(idx) {
                bytes.extend_from_slice(&component.to_le_bytes());
            }
        }
        bytes.extend_from_slice(&[0u8; 2]);
    }

    let path = dir.join(name);
    std::fs::write(&path, bytes).expect("the fixture is writable");
    path.to_string_lossy().replace('\\', "/")
}

/// A scratch directory of this test's own, removed by the OS eventually and unique per run.
fn scratch(tag: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("dualis-assembly-{tag}"));
    std::fs::create_dir_all(&dir).expect("a scratch directory");
    dir
}

/// **Two parts from two files land where their files put them, with nothing between them.**
///
/// The left brick spans 0–20 mm and the right 25–40 mm on a 5 mm grid, so there is one empty
/// column between them — and it is *empty*, not the block's bulk material, which is what makes
/// this an assembly in air rather than two inclusions in a billet.
#[test]
fn two_parts_from_files_assemble_with_a_gap_of_nothing() {
    let dir = scratch("gap");
    let left = brick_stl(&dir, "left.stl", [0.0, 0.0, 0.0], [20.0, 10.0, 10.0]);
    let right = brick_stl(&dir, "right.stl", [25.0, 0.0, 0.0], [40.0, 10.0, 10.0]);

    let json = format!(
        r#"{{
  "title": "two parts in air",
  "duration_s": 1.0,
  "frames": 1,
  "domains": [
    {{ "kind": "block", "name": "assembly", "cells": [8, 2, 2], "cell_mm": 5.0,
      "initial_c": 20.0,
      "parts": [
        {{ "stl": "{left}", "material": "copper" }},
        {{ "stl": "{right}", "material": "aluminium" }}
      ] }}
  ]
}}"#
    );

    let world = World::build(scene(&json)).expect("the assembly builds");
    let block = world
        .simulation()
        .domain_as::<dualis::thermal::Solid3D>("assembly")
        .expect("the block is there");

    // Columns 0..4 are the left part, 5..8 the right, and column 4 is the gap.
    for i in 0..8 {
        let void = block.is_void(i, 0, 0);
        assert_eq!(
            void,
            i == 4,
            "column {i} should be {}",
            if i == 4 { "empty" } else { "part of something" }
        );
    }
    assert_eq!(
        block.substance_at(0, 0, 0).name,
        dualis::prelude::Substance::copper().name
    );

    // And the build reported what voxelising cost, per part.
    assert_eq!(world.notes().len(), 2, "{:?}", world.notes());
    assert!(
        world.notes()[0].contains("volume") && world.notes()[0].contains("boundary cells"),
        "{}",
        world.notes()[0]
    );
}

/// **Heat crosses a join and does not cross a gap.**
///
/// The physics that makes void worth having, through the file path. The left part starts hot and
/// nothing else feeds the block, so the only way the far end warms is conduction across whatever
/// is between the two parts.
///
/// **No heater**, and that is a correction rather than a simplification: the first version of
/// this test drove both cases with one and read the gapped assembly as *hotter*. It has one cell
/// fewer, so the same joules spread over less capacity — the test was measuring heat capacity
/// and calling it conduction. With nothing arriving, the far part's temperature is about the
/// path and nothing else.
#[test]
fn heat_crosses_a_join_and_not_a_gap() {
    let dir = scratch("join");
    let left = brick_stl(&dir, "left.stl", [0.0, 0.0, 0.0], [20.0, 10.0, 10.0]);
    let touching = brick_stl(&dir, "touching.stl", [20.0, 0.0, 0.0], [40.0, 10.0, 10.0]);
    let apart = brick_stl(&dir, "apart.stl", [25.0, 0.0, 0.0], [40.0, 10.0, 10.0]);

    let run_with = |right: &str| {
        let json = format!(
            r#"{{
  "title": "a join or a gap",
  "duration_s": 60.0,
  "frames": 2,
  "domains": [
    {{ "kind": "block", "name": "assembly", "cells": [8, 2, 2], "cell_mm": 5.0,
      "initial_c": 20.0,
      "parts": [
        {{ "stl": "{left}", "material": "copper" }},
        {{ "stl": "{right}", "material": "copper" }}
      ] }}
  ]
}}"#
        );
        let mut world = World::build(scene(&json)).expect("builds");
        {
            let block = world
                .simulation_mut()
                .domain_as_mut::<dualis::thermal::Solid3D>("assembly")
                .expect("the block is there");
            for k in 0..2 {
                for j in 0..2 {
                    for i in 0..4 {
                        block.set_temperature(i, j, k, dualis::units::Temperature::celsius(300.0));
                    }
                }
            }
        }
        world.run().expect("runs");
        let block = world
            .simulation()
            .domain_as::<dualis::thermal::Solid3D>("assembly")
            .expect("the block is there");
        block.temperature_at(7, 0, 0).to_si()
    };

    let start = dualis::units::Temperature::celsius(20.0).to_si();
    let joined = run_with(&touching);
    let gapped = run_with(&apart);
    assert!(
        joined > start + 50.0,
        "a shared face should carry the heat: {joined:.2} K from {start:.2} K"
    );
    assert!(
        (gapped - start).abs() < 1e-9,
        "and nothing carries nothing: {gapped:.6} K should still be {start:.2} K"
    );
}

/// **A part that reaches outside the block is refused with both boxes named**, rather than
/// cropped into a different shape.
#[test]
fn a_part_that_does_not_fit_is_refused() {
    let dir = scratch("toobig");
    let big = brick_stl(&dir, "big.stl", [0.0, 0.0, 0.0], [60.0, 10.0, 10.0]);
    let json = format!(
        r#"{{
  "title": "a part too long for its block",
  "duration_s": 1.0,
  "frames": 1,
  "domains": [
    {{ "kind": "block", "name": "assembly", "cells": [8, 2, 2], "cell_mm": 5.0,
      "initial_c": 20.0,
      "parts": [ {{ "stl": "{big}", "material": "copper" }} ] }}
  ]
}}"#
    );
    let why = World::build(scene(&json))
        .err()
        .expect("a part longer than its block must refuse");
    assert!(why.contains("cut off"), "{why}");
}

/// **A missing file names itself**, and a set of parts that voxelised to nothing is refused
/// rather than run as an empty box.
#[test]
fn a_missing_file_and_an_empty_assembly_are_both_refused() {
    let json = |parts: &str| {
        format!(
            r#"{{
  "title": "nothing to assemble",
  "duration_s": 1.0,
  "frames": 1,
  "domains": [
    {{ "kind": "block", "name": "assembly", "cells": [4, 2, 2], "cell_mm": 5.0,
      "initial_c": 20.0, "parts": [{parts}] }}
  ]
}}"#
        )
    };

    let why = World::build(scene(&json(
        r#"{ "stl": "definitely-not-here.stl", "material": "copper" }"#,
    )))
    .err()
    .expect("a missing mesh must refuse");
    assert!(why.contains("definitely-not-here.stl"), "{why}");

    // A brick far smaller than one cell rasterises to nothing at all.
    let dir = scratch("tiny");
    let tiny = brick_stl(&dir, "tiny.stl", [1.0, 1.0, 1.0], [1.2, 1.2, 1.2]);
    let why = World::build(scene(&json(&format!(
        r#"{{ "stl": "{tiny}", "material": "copper" }}"#
    ))))
    .err()
    .expect("an assembly of nothing must refuse");
    assert!(why.contains("no cells"), "{why}");
}
