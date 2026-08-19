//! Two CAD parts, rasterised onto one grid, touching — and heat crossing between them.
//!
//! `ARCHITECTURE.md` named this as the gap that arrives first in practice: "each rasterises onto
//! its own grid from its own bounding box. Placing them in one grid needs a shared origin and a
//! rule for a cell both of them claim, and neither exists."
//!
//! The shared origin is [`Voxels::onto`]. The rule for a contested cell turns out to need no new
//! mechanism at all: a domain fills from each part in turn, so a cell both claim is the later
//! one's, which is the same **last writer wins** that `regions` has documented since it arrived —
//! "applied in order, so a later region overwrites an earlier one where they overlap ... it is
//! how a coating *on* a layer is written and there is no other way to mean it".
//!
//! And once they share a grid they are already coupled: `Solid3D`'s stencil crosses the face
//! between two cells of different materials with the harmonic mean of their conductivities, which
//! this workspace measured as *exact* for a layered wall at every resolution. No interface, no
//! bus channel, no new physics — the parts touch because they are cells of one array.

use pantometry::prelude::*;
use pantometry::shape::{Mesh, Triangle, Voxels};
use pantometry::thermal::Solid3D;

/// An axis-aligned box as twelve triangles, wound outward, in metres.
///
/// Written out rather than taken from a file so the test states its own geometry: an STL fixture
/// would put the thing under test one parser away from the assertion.
fn brick(low: [f64; 3], high: [f64; 3]) -> Vec<Triangle> {
    let v = |i: usize| {
        // Through `LengthVec` because the vertex type is glam's and the facade does not
        // re-export it — a test reaches the same value the units crate already carries.
        LengthVec::m(
            if i & 1 == 0 { low[0] } else { high[0] },
            if i & 2 == 0 { low[1] } else { high[1] },
            if i & 4 == 0 { low[2] } else { high[2] },
        )
        .to_si()
    };
    // Each face as two triangles, wound so the normal points out of the box.
    let quads = [
        // -x, +x
        ([0, 4, 6], [0, 6, 2]),
        ([1, 3, 7], [1, 7, 5]),
        // -y, +y
        ([0, 1, 5], [0, 5, 4]),
        ([2, 6, 7], [2, 7, 3]),
        // -z, +z
        ([0, 2, 3], [0, 3, 1]),
        ([4, 5, 7], [4, 7, 6]),
    ];
    let mut out = Vec::new();
    for (a, b) in quads {
        out.push(Triangle {
            a: v(a[0]),
            b: v(a[1]),
            c: v(a[2]),
        });
        out.push(Triangle {
            a: v(b[0]),
            b: v(b[1]),
            c: v(b[2]),
        });
    }
    out
}

/// **Two parts rasterised onto one grid occupy the cells their own coordinates put them in, and
/// meet at the plane they share.**
///
/// The left brick spans `x = 0..20 mm` and the right `20..40 mm`, so on a 5 mm grid the left owns
/// columns 0 to 3 and the right 4 to 7, with no cell claimed twice and no gap between them. Both
/// were rasterised onto the *same* origin, which is the whole difference from `Voxels::of`.
#[test]
fn two_parts_land_where_their_own_coordinates_put_them() {
    let cell = Length::mm(5.0);
    let origin = LengthVec::m(0.0, 0.0, 0.0);
    let counts = (8, 2, 2);

    let left = Mesh::new(brick([0.0, 0.0, 0.0], [0.020, 0.010, 0.010]));
    let right = Mesh::new(brick([0.020, 0.0, 0.0], [0.040, 0.010, 0.010]));

    let a = Voxels::onto(&left, origin, counts, cell).expect("the left part fits");
    let b = Voxels::onto(&right, origin, counts, cell).expect("the right part fits");

    for k in 0..2 {
        for j in 0..2 {
            for i in 0..8 {
                let in_a = a.contains(i, j, k);
                let in_b = b.contains(i, j, k);
                assert!(
                    !(in_a && in_b),
                    "cell ({i},{j},{k}) is claimed by both parts"
                );
                assert!(
                    in_a || in_b,
                    "cell ({i},{j},{k}) belongs to neither, so the two do not touch"
                );
                assert_eq!(
                    in_a,
                    i < 4,
                    "the left part should own columns 0..4, not {i}"
                );
            }
        }
    }
}

/// **A part that would not fit the grid is refused rather than cropped.**
///
/// Silently cutting a corner off gives a run that audits clean and answers about a different
/// shape — the failure this workspace keeps finding. The message names both boxes so a caller
/// can see which way to move or grow.
#[test]
fn a_part_that_does_not_fit_is_refused_with_both_boxes() {
    let cell = Length::mm(5.0);
    let too_big = Mesh::new(brick([0.0, 0.0, 0.0], [0.060, 0.010, 0.010]));
    let why = Voxels::onto(&too_big, LengthVec::m(0.0, 0.0, 0.0), (8, 2, 2), cell)
        .expect_err("a part longer than the grid must refuse");
    assert!(
        why.contains("cut off") && why.contains("0.0600"),
        "the refusal should name both boxes: {why}"
    );
}

/// **Heat crosses between two assembled parts, and it crosses at the rate the two materials give
/// — not at either one's.**
///
/// This is the payoff. Two bricks of *different* metals filled into one `Solid3D` from their own
/// meshes: copper on the left starting hot, stainless on the right starting cold. Nothing couples
/// them but the grid, and the stencil's harmonic-mean face conductivity is what carries the heat.
///
/// The check is a **conservation** statement plus a **direction** one, both of which the domain
/// does not compute for itself: the total enthalpy is unchanged to the bit because no face is
/// exposed and nothing is on the bus, and the cold part warms while the hot part cools. A pair
/// that did not touch would leave both exactly where they started, which is what this looked like
/// before `onto` existed.
#[test]
fn heat_crosses_between_two_assembled_parts() {
    let cell = Length::mm(5.0);
    let origin = LengthVec::m(0.0, 0.0, 0.0);
    let counts = (8, 2, 2);

    let left = Mesh::new(brick([0.0, 0.0, 0.0], [0.020, 0.010, 0.010]));
    let right = Mesh::new(brick([0.020, 0.0, 0.0], [0.040, 0.010, 0.010]));
    let a = Voxels::onto(&left, origin, counts, cell).expect("fits");
    let b = Voxels::onto(&right, origin, counts, cell).expect("fits");

    // Both parts are filled from their own meshes. The right one could have been left as the
    // block's bulk and the test would still pass -- which is why it is not: an assembly is two
    // parts each put where its own file says, and a test that only placed one would be checking
    // half of that.
    let steel = Substance::from_name("stainless_304").expect("the catalogue has it");
    let mut block = Solid3D::new(
        "assembly",
        steel.clone(),
        counts,
        cell,
        Temperature::celsius(20.0),
    );
    // The left part is copper; the right keeps the block's own steel. Filled from the mesh's own
    // predicate, which is the only thing that says where the part is.
    block.fill(Substance::copper(), |i, j, k| a.contains(i, j, k));
    block.fill(steel.clone(), |i, j, k| b.contains(i, j, k));
    for k in 0..counts.2 {
        for j in 0..counts.1 {
            for i in 0..counts.0 {
                if a.contains(i, j, k) {
                    block.set_temperature(i, j, k, Temperature::celsius(300.0));
                }
            }
        }
    }

    let hot_before = (0..counts.2)
        .flat_map(|k| (0..counts.1).map(move |j| (j, k)))
        .flat_map(|(j, k)| (0..4).map(move |i| (i, j, k)))
        .map(|(i, j, k)| block.temperature_at(i, j, k).to_si())
        .fold(f64::MIN, f64::max);
    let cold_before = block.temperature_at(7, 0, 0).to_si();
    let opening = block.ledger().get(quantity::ENERGY).unwrap_or(0.0);

    let dt = block.max_stable_dt(Time::ZERO);
    let mut bus = Exchange::new();
    for _ in 0..400 {
        block.step(Time::ZERO, dt, &mut bus).expect("a stable step");
    }

    let cold_after = block.temperature_at(7, 0, 0).to_si();
    let closing = block.ledger().get(quantity::ENERGY).unwrap_or(0.0);

    assert!(
        cold_after > cold_before + 1.0,
        "the far end of the cold part should have warmed: {cold_before:.2} K to {cold_after:.2} K"
    );
    assert!(
        cold_after < hot_before,
        "and not past the hot part it is drawing from"
    );
    // Nothing leaves an assembly with no exposed face: the enthalpy is unchanged to the bit,
    // which is the `Σ Cᵢ Tᵢ` the stencil conserves by construction.
    let scale = opening.abs().max(closing.abs()).max(1.0);
    assert!(
        (closing - opening).abs() / scale < 1e-12,
        "an unexposed assembly conserves its enthalpy: {opening:e} to {closing:e}"
    );
}

/// **A cell two parts both claim goes to the later one, which is the rule `regions` already
/// documents.** Assembly needed no new rule, and a rule invented here would have been a second
/// answer to a question this format had already answered.
#[test]
fn a_contested_cell_belongs_to_the_part_filled_last() {
    let cell = Length::mm(5.0);
    let origin = LengthVec::m(0.0, 0.0, 0.0);
    let counts = (4, 2, 2);

    // Deliberately overlapping: both cover x = 5..10 mm.
    let first = Mesh::new(brick([0.0, 0.0, 0.0], [0.010, 0.010, 0.010]));
    let second = Mesh::new(brick([0.005, 0.0, 0.0], [0.020, 0.010, 0.010]));
    let a = Voxels::onto(&first, origin, counts, cell).expect("fits");
    let b = Voxels::onto(&second, origin, counts, cell).expect("fits");
    assert!(
        a.contains(1, 0, 0) && b.contains(1, 0, 0),
        "the overlap is real"
    );

    let mut block = Solid3D::new(
        "assembly",
        Substance::from_name("stainless_304").expect("the catalogue has it"),
        counts,
        cell,
        Temperature::celsius(20.0),
    );
    block.fill(Substance::copper(), |i, j, k| a.contains(i, j, k));
    block.fill(Substance::aluminium_6061(), |i, j, k| b.contains(i, j, k));

    assert_eq!(
        block.substance_at(1, 0, 0).name,
        Substance::aluminium_6061().name,
        "the contested cell should be the later part's"
    );
    assert_eq!(
        block.substance_at(0, 0, 0).name,
        Substance::copper().name,
        "and the uncontested one the earlier part's"
    );
}
