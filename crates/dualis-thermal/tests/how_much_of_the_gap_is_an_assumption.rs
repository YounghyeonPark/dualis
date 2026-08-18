//! What a clearance's geometry says about how much of its answer is a boundary condition.
//!
//! A gap is charged `σA(T₁⁴−T₂⁴)/(1/ε₁+1/ε₂−1)` with a view factor of **one**, and that is exact
//! when the sides of the gap are mirrors. The block's own outer faces are mirrors: an insulated
//! boundary is implemented as one, and a mirror puts an image of each surface beyond it, which
//! extends both to infinity. So for every clearance bounded laterally by the block's edge — which
//! is every clearance in this repository — the exchange is the infinite-parallel-plate result and
//! there is nothing approximate about it.
//!
//! What is *not* stated anywhere is how different the answer would be if somebody meant the gap to
//! be open to space, which is what two small plates floating in vacuum are. That is `F₁₂`, the
//! ordinary view factor, and the two agree only when the gap is narrow compared with the surfaces.
//! [`GapPatch::view_factor`] computes it so the difference is a number rather than a caveat.

use dualis_core::units::{Length, Temperature};
use dualis_core::Substance;
use dualis_thermal::{GapPatch, Solid3D};

/// **The view factor is the published closed form**, checked against values computed from it
/// independently and against the two limits it has to have.
///
/// `F(1,1) = 0.19982` for two squares one side-length apart is the number in every radiation text,
/// and it is not a number this code could produce by accident. The limits are the honest half: a
/// gap that closes sees everything, and a gap that opens sees an area over a sphere.
#[test]
fn the_view_factor_is_the_parallel_rectangle_form() {
    let patch = |span: f64, distance: f64| GapPatch {
        pairs: 1,
        span: (span, span),
        distance,
        rectangular: true,
    };

    // The textbook value, to five figures.
    let square = patch(1.0, 1.0).view_factor();
    assert!(
        (square - 0.199_82).abs() < 5e-6,
        "two squares one side apart see 0.19982 of each other, computed {square:.6}"
    );

    // **Closing the gap goes to one.** Not asserted at a point but as a limit, because the whole
    // claim of the model's `F = 1` is that it is the narrow-gap answer.
    let mut previous = 0.0;
    for n in 1..=6 {
        let f = patch(10f64.powi(n), 1.0).view_factor();
        assert!(f > previous, "closing the gap must see more: {f}");
        previous = f;
    }
    assert!(
        previous > 0.999,
        "a gap narrow against its plates sees essentially all of them: {previous:.6}"
    );

    // **Opening it goes as `A/(πc²)`**, the solid angle a small patch subtends, and it goes there
    // *monotonically*. The second half is the one that found something: the published form is a
    // small difference of larger terms, so at `X` of 1e-3 it has cancelled four digits away and
    // the ratio turns around — measured 0.999993 at 1/300, 1.00009 at 1/1000, 1.007 at 1/3000,
    // growing where it should converge. A closed form that is quietly worse the further apart two
    // surfaces are is worth catching, and the fix is the series below that crossover.
    let mut worst = 0.0f64;
    for c in [10.0, 30.0, 100.0, 300.0, 1_000.0, 10_000.0, 1e6] {
        let far = patch(1.0, c);
        let solid_angle = 1.0 / (std::f64::consts::PI * c * c);
        let off = (far.view_factor() / solid_angle - 1.0).abs();
        // The physical term is `(X²+Y²)/3`, so at a separation of `c` the deviation is `2/(3c²)`.
        // Asserted against that rather than against a tolerance, which is what makes it a check of
        // the *form* and not of one number.
        let expected = 2.0 / (3.0 * c * c);
        assert!(
            off <= expected * 1.05 + 1e-9,
            "at c = {c} the deviation should be 2/(3c^2) = {expected:e}, is {off:e}"
        );
        if c >= 300.0 {
            worst = worst.max(off);
        }
    }
    assert!(
        worst < 1e-5,
        "a distant patch must keep converging rather than turn around: {worst:e}"
    );

    // The degenerate ends are the form's own limits and not sentinels, which matters because this
    // is public and a caller can build a patch the block would never produce. A surface with no
    // extent sees nothing; two surfaces with nothing between them see each other entirely. Neither
    // may be a `NaN`, which is the value that would reach a plot and be believed.
    assert_eq!(patch(0.0, 1.0).view_factor(), 0.0);
    assert_eq!(patch(1.0, 0.0).view_factor(), 1.0);
    assert!(patch(1.0, f64::NAN).view_factor().is_finite());
    assert!(patch(f64::NAN, 1.0).view_factor().is_finite());
}

/// **A patch is the sheet of facing area, not the pairs it is made of.**
///
/// The grouping is what makes a view factor meaningful: `F` for a 4x4 sheet at a two-cell
/// separation is nothing like `F` for one cell pair at the same separation, and charging the
/// second where the first is true would be a different wrong answer from the one being measured.
#[test]
fn a_patch_is_a_sheet_and_not_its_pairs() {
    // A 6x6x7 block: a 4x4 part on the bottom, two empty layers, a lid across the whole top.
    let mut block = Solid3D::new(
        "housing",
        Substance::copper(),
        (6, 6, 7),
        Length::mm(8.0),
        Temperature::celsius(20.0),
    );
    block = block
        .empty(|i, j, k| k < 3 && (i == 0 || i == 5 || j == 0 || j == 5))
        .empty(|_, _, k| k == 3 || k == 4);

    let patches = block.gap_patches();
    assert_eq!(patches.len(), 1, "one clearance, one patch: {patches:?}");
    let p = patches[0];
    assert_eq!(
        p.pairs, 16,
        "the 4x4 part faces the lid across sixteen pairs"
    );
    assert!(p.rectangular, "and the sheet fills its bounding box");
    assert!(
        (p.span.0 - 0.032).abs() < 1e-12 && (p.span.1 - 0.032).abs() < 1e-12,
        "four 8 mm cells across: {:?}",
        p.span
    );
    assert!(
        (p.distance - 0.016).abs() < 1e-12,
        "two 8 mm cells apart: {}",
        p.distance
    );

    // **And this is the number the caveat was hiding.** X = Y = 32/16 = 2, so a gap open to space
    // would carry 0.415 of what the block charges — the block is right for its own mirrored
    // boundary and this says how much that reading is worth.
    let f = p.view_factor();
    assert!(
        (f - 0.415_25).abs() < 5e-5,
        "F(2,2) is 0.41525, computed {f:.6}"
    );
}

/// **Two clearances at the same width are two patches**, not one sheet that describes neither.
///
/// The failure this prevents is an average: two 1x1 gaps at opposite corners have the same
/// separation as one 2x2 gap and a completely different view factor, and a grouping that only
/// sorted by width would report the wrong one for both.
#[test]
fn separate_clearances_stay_separate() {
    let mut block = Solid3D::new(
        "two",
        Substance::copper(),
        (5, 1, 3),
        Length::mm(5.0),
        Temperature::celsius(20.0),
    );
    // Cells 0 and 4 are pillars, 1..4 are solid at top and bottom with nothing between.
    block = block.empty(|i, _, k| k == 1 && (i == 0 || i == 4));

    let mut patches = block.gap_patches();
    patches.sort_by_key(|p| p.pairs);
    assert_eq!(patches.len(), 2, "two pillars, two patches: {patches:?}");
    for p in &patches {
        assert_eq!(p.pairs, 1);
        assert!((p.distance - 0.005).abs() < 1e-12);
        // One 5 mm cell facing across 5 mm: X = Y = 1, the textbook square.
        assert!((p.view_factor() - 0.199_82).abs() < 5e-6, "{p:?}");
    }
}

/// **A block with no void has no patches**, which is every block that existed before there was a
/// void at all — and the guard that says this cost nothing to a scene without one.
#[test]
fn a_solid_block_has_no_clearances() {
    let block = Solid3D::new(
        "solid",
        Substance::copper(),
        (4, 4, 4),
        Length::mm(5.0),
        Temperature::celsius(20.0),
    );
    assert!(block.gap_patches().is_empty());
}
