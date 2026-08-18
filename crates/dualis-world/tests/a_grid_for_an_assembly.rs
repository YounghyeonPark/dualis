//! Choosing a cell size by measuring, which is the third of `ARCHITECTURE.md`'s three assembly
//! gaps and the one it left open on purpose.
//!
//! The reason it was left open is worth keeping in front of the tests: it "would be the first place
//! in the workspace that *guesses* — so it has to guess visibly". Measuring every candidate removes
//! the guess from the numbers entirely; what is left is the **rule** for which row to recommend,
//! and these check that the rule does what its one sentence says on geometry whose answer can be
//! counted by hand.

use dualis::shape::Mesh;
use dualis_world::fit::propose;

/// A binary STL of an axis-aligned brick, in **millimetres**, which is what `Mesh::from_stl` reads.
fn brick(low: [f32; 3], high: [f32; 3]) -> Mesh {
    let v = |i: usize| {
        [
            if i & 1 == 0 { low[0] } else { high[0] },
            if i & 2 == 0 { low[1] } else { high[1] },
            if i & 4 == 0 { low[2] } else { high[2] },
        ]
    };
    let quads = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    let mut bytes = vec![0u8; 80];
    bytes.extend_from_slice(&(12u32).to_le_bytes());
    for q in quads {
        for tri in [[q[0], q[1], q[2]], [q[0], q[2], q[3]]] {
            bytes.extend_from_slice(&[0u8; 12]); // the normal, which the reader recomputes
            for c in tri {
                for x in v(c) {
                    bytes.extend_from_slice(&x.to_le_bytes());
                }
            }
            bytes.extend_from_slice(&[0u8; 2]);
        }
    }
    Mesh::from_stl(&bytes).expect("the fixture is a mesh")
}

/// **The ladder is a statement about the assembly, not about millimetres.**
///
/// A 40 x 20 x 10 mm brick has a thinnest dimension of 10 mm, so the rows give that 1, 2, 4 … cells
/// and the cell sizes are 10, 5, 2.5 mm — every one of them a sentence a person can check. The grid
/// each implies is the extent divided by the cell, which for a brick on cell boundaries is exact.
#[test]
fn the_rows_resolve_the_thinnest_feature_by_powers_of_two() {
    let fit = propose(
        &[("brick.stl".to_string(), brick([0.0; 3], [40.0, 20.0, 10.0]))],
        100_000,
    )
    .expect("a brick has an extent");

    assert!(
        (fit.extent_m[0] - 0.040).abs() < 1e-12
            && (fit.extent_m[1] - 0.020).abs() < 1e-12
            && (fit.extent_m[2] - 0.010).abs() < 1e-12,
        "{:?}",
        fit.extent_m
    );

    for (n, want_mm) in [(0usize, 10.0), (1, 5.0), (2, 2.5)] {
        let c = &fit.candidates[n];
        assert!(
            (c.cell_m * 1e3 - want_mm).abs() < 1e-12,
            "row {n} should be {want_mm} mm, is {:.4}",
            c.cell_m * 1e3
        );
        let across = 1usize << n;
        assert_eq!(
            c.counts,
            (4 * across, 2 * across, across),
            "row {n}: the grid is the extent over the cell"
        );
    }
}

/// **A part finer than the grid disappears, and the table says so.**
///
/// The failure this whole tool exists to prevent: rasterising a rib below the cell size does not
/// error, it produces nothing, and the run is perfectly well behaved about a different object. Two
/// parts, a 1 mm plate with a 0.2 mm shim on top of it, and `sound()` is what says a row is
/// unusable rather than merely rough. The two are not the same judgement and the tool must not
/// blur them: a rough grid is a choice, a grid missing a part is an answer about something else.
#[test]
fn a_part_finer_than_the_grid_vanishes_and_the_row_is_not_sound() {
    // The shim's thickness is what sets the ladder, so row 0 gives it exactly one cell — the case
    // where a rasteriser's tie-breaking decides whether it exists at all.
    let parts = vec![
        (
            "plate.stl".to_string(),
            brick([0.0, 0.0, 0.0], [6.0, 6.0, 1.0]),
        ),
        (
            "shim.stl".to_string(),
            brick([0.0, 0.0, 1.0], [6.0, 6.0, 1.2]),
        ),
    ];
    let fit = propose(&parts, 2_000_000).expect("two bricks have an extent");

    // Row 0 gives the thinnest feature one cell, which is exactly the case where a 0.4 mm shim is
    // one cell tall and may or may not survive the rasteriser's tie-breaking. Whatever it does,
    // `sound()` has to agree with `filled`.
    for c in &fit.candidates {
        let vanished = c.parts.iter().any(|p| p.filled == 0);
        assert_eq!(
            c.sound(),
            !vanished && c.parts.iter().all(|p| p.ambiguous_rows == 0),
            "a row is sound exactly when every part is present and nothing was undecidable: {c:?}"
        );
    }

    // And by row two, which is four cells across the thinnest feature by the ladder's own
    // definition, it is certainly there. Indexed rather than searched by cell size: an STL stores
    // f32, so a 0.2 mm shim is 0.20000005 mm and a nominal comparison misses by 1e-11.
    assert!(fit.candidates.len() > 2, "{:?}", fit.candidates.len());
    let fine = &fit.candidates[2];
    assert!(
        fine.parts.iter().all(|p| p.filled > 0),
        "both parts are there once the thin one is resolved: {fine:?}"
    );
}

/// **The recommendation is the coarsest sound row under the stated uncertainty**, which is what its
/// one sentence says and is checked here against the table it came from rather than against a
/// remembered number.
///
/// The check is deliberately circular-proof: it finds the recommended row, asserts it satisfies the
/// rule, and then asserts every *earlier* row fails it. A rule that returned the finest row, or the
/// first row, or a row that happened to be good, would pass one of those and not both.
#[test]
fn the_recommendation_is_the_coarsest_row_that_qualifies() {
    let fit = propose(
        &[("brick.stl".to_string(), brick([0.0; 3], [40.0, 20.0, 10.0]))],
        5_000_000,
    )
    .expect("a brick has an extent");

    let uncertainty = 0.5;
    let picked = fit
        .recommended(uncertainty)
        .expect("some row resolves a brick to under half boundary");

    assert!(picked.sound() && picked.worst_boundary() <= uncertainty);
    for c in &fit.candidates {
        if std::ptr::eq(c, picked) {
            break;
        }
        assert!(
            !(c.sound() && c.worst_boundary() <= uncertainty),
            "a coarser row qualified and was not picked: {c:?}"
        );
    }

    // An uncertainty nothing can meet is a `None` rather than a shrug — a real answer that the
    // parts are finer than the budget allows.
    assert!(
        fit.recommended(0.0).is_none() || fit.recommended(0.0).unwrap().worst_boundary() == 0.0,
        "nothing may be recommended above the bar it was given"
    );
}

/// **Boundary fraction falls as the grid refines and volume error does not**, which is the whole
/// reason the rule steers by the first.
///
/// `Loss::volume_error` documents this about itself — a sphere at three resolutions gives
/// `+4.9%, +5.8%, −2.3%` — and a tool that extrapolated from one row to the next would be
/// confidently wrong. Here it is measured on the assembly's own numbers rather than quoted.
#[test]
fn the_rule_steers_by_the_quantity_that_actually_falls() {
    let fit = propose(
        &[("brick.stl".to_string(), brick([0.0; 3], [37.0, 23.0, 11.0]))],
        20_000_000,
    )
    .expect("a brick has an extent");

    let boundary: Vec<f64> = fit.candidates.iter().map(|c| c.worst_boundary()).collect();
    assert!(
        boundary.len() >= 5,
        "enough rows to see a trend: {boundary:?}"
    );
    for pair in boundary.windows(2) {
        assert!(
            pair[1] <= pair[0] + 1e-12,
            "refining must not make the rasterisation *less* certain: {boundary:?}"
        );
    }
    assert!(
        boundary[0] > 0.9 && *boundary.last().expect("rows") < 0.2,
        "and it should fall a long way: {boundary:?}"
    );
}

/// **The fragment it prints is a scene that builds.**
///
/// A table nobody can act on is a table. This asserts the suggested `cells`, `cell_mm` and `parts`
/// go into a domain that `World::build_with` accepts and that the part is really there — which is
/// the round trip from "I have a CAD file" to "it is running" with nothing in between to get wrong.
#[test]
fn the_fragment_it_suggests_is_a_scene_that_builds() {
    let mesh = brick([0.0; 3], [40.0, 20.0, 10.0]);
    // Re-emit the same brick as STL so the scene has something to read.
    let bytes = {
        let stl = |low: [f32; 3], high: [f32; 3]| {
            let v = |i: usize| {
                [
                    if i & 1 == 0 { low[0] } else { high[0] },
                    if i & 2 == 0 { low[1] } else { high[1] },
                    if i & 4 == 0 { low[2] } else { high[2] },
                ]
            };
            let quads = [
                [0, 2, 3, 1],
                [4, 5, 7, 6],
                [0, 1, 5, 4],
                [2, 6, 7, 3],
                [0, 4, 6, 2],
                [1, 3, 7, 5],
            ];
            let mut b = vec![0u8; 80];
            b.extend_from_slice(&(12u32).to_le_bytes());
            for q in quads {
                for tri in [[q[0], q[1], q[2]], [q[0], q[2], q[3]]] {
                    b.extend_from_slice(&[0u8; 12]);
                    for c in tri {
                        for x in v(c) {
                            b.extend_from_slice(&x.to_le_bytes());
                        }
                    }
                    b.extend_from_slice(&[0u8; 2]);
                }
            }
            b
        };
        stl([0.0; 3], [40.0, 20.0, 10.0])
    };

    let fit = propose(&[("brick.stl".to_string(), mesh)], 200_000).expect("it fits");
    let picked = fit.recommended(0.6).expect("a row qualifies");
    let fragment = fit.scene_fragment(picked, "aluminium");

    let json = format!(
        "{{ \"title\": \"from the fitter\", \"duration_s\": 1.0, \"frames\": 2,\n  \
         \"domains\": [{{ \"kind\": \"block\", \"name\": \"assembly\", \"initial_c\": 20.0,\n{fragment} }}] }}"
    );
    let scene: dualis_world::Scene =
        serde_json::from_str(&json).unwrap_or_else(|e| panic!("{e}\n{json}"));
    let files = dualis_world::Uploaded::new().with("brick.stl", bytes);
    let world = dualis_world::World::build_with(scene, &files)
        .unwrap_or_else(|e| panic!("the suggested grid should build: {e}\n{json}"));

    let block = world
        .simulation()
        .domain_as::<dualis::thermal::Solid3D>("assembly")
        .expect("it is a block");
    assert!(
        block.void_cells() < picked.total_cells,
        "the part has to be in there: {} void of {}",
        block.void_cells(),
        picked.total_cells
    );
}

/// **Nothing to measure is refused with what was missing**, rather than answered with an empty
/// table somebody would read as "no problems".
#[test]
fn an_empty_assembly_is_refused() {
    let err = propose(&[], 1000).expect_err("no parts is no assembly");
    assert!(err.contains("at least one part"), "{err}");

    let flat = propose(
        &[("flat.stl".to_string(), brick([0.0; 3], [10.0, 10.0, 0.0]))],
        1000,
    );
    assert!(
        flat.is_err(),
        "a part with no thickness has nothing to resolve"
    );
}
