//! What a *reader* sees where there is nothing, which is not what a void cell contains.
//!
//! Void arrived and broke conduction correctly, and every output kept drawing it. The cause is one
//! line: `ScalarField::at` read `self.cells` directly, and an emptied cell still holds whatever it
//! held when it was emptied. So a clearance exported as a piece of the block sitting at its start
//! temperature forever — a plausible number that never moves, which is worse than an obviously
//! wrong one because nothing about it looks wrong.
//!
//! Measured on scene `23`, a part under a lid with a two-layer gap: the glTF carried **252** points
//! for a grid with **120** solid cells, so two thirds of what left this workspace for Blender was
//! material that is not there. The tests below are on the sampler, which is where the defect was;
//! `dualis-view` has its own for what the writers do with a sample that is not a number.

use dualis_core::units::{Length, LengthVec, Temperature, Time};
use dualis_core::{ScalarField, Substance};
use dualis_thermal::Solid3D;

/// A 4x1x1 bar of 10 mm cells with the two middle cells emptied.
fn bar() -> Solid3D {
    let mut block = Solid3D::new(
        "bar",
        Substance::copper(),
        (4, 1, 1),
        Length::mm(10.0),
        Temperature::celsius(20.0),
    );
    block.set_temperature(0, 0, 0, Temperature::celsius(300.0));
    block.empty(|i, _, _| i == 1 || i == 2)
}

/// **A sample inside a clearance is not a number**, and one on a solid cell centre is that cell.
///
/// The two halves are one claim: the mask has to be exact where there is material, or a fix for
/// the gap would be a blur everywhere else. Cell centres sit at `(i + 1/2)·dx`.
#[test]
fn a_sample_in_the_gap_has_no_value_and_one_on_the_metal_is_exact() {
    let block = bar();
    let at = |mm: f64| {
        block.at(
            LengthVec::from_si(glam::DVec3::new(mm * 1e-3, 5e-3, 5e-3)),
            Time::ZERO,
        )
    };

    assert!(
        (at(5.0) - Temperature::celsius(300.0).to_si()).abs() < 1e-12,
        "the hot cell's centre reads the hot cell: {:.6} K",
        at(5.0)
    );
    assert!(
        (at(35.0) - Temperature::celsius(20.0).to_si()).abs() < 1e-12,
        "and the far cell reads the far cell: {:.6} K",
        at(35.0)
    );
    for mm in [15.0, 20.0, 25.0] {
        assert!(
            at(mm).is_nan(),
            "there is nothing at {mm} mm to have a temperature, and it read {:.4} K",
            at(mm)
        );
    }
}

/// **The sample belongs to the cell it is in**, so the gap's edge is where the material ends and
/// not half a cell further in.
///
/// This is the half that was still wrong after the mask alone. `Scene::capture` samples a field on
/// a grid spread across the block's *extent*, which does not land on cell centres — for a six-cell
/// axis it lands at `0, 1/5, 2/5 …` of the span. Those samples straddle the boundary, the masked
/// weights answer with whichever solid corner they lean towards, and a whole layer of the
/// clearance came back solid. Measured before the guard: 156 points for 120 solid cells.
#[test]
fn a_sample_that_straddles_the_boundary_belongs_to_the_cell_it_is_in() {
    let block = bar();
    let at = |mm: f64| {
        block.at(
            LengthVec::from_si(glam::DVec3::new(mm * 1e-3, 5e-3, 5e-3)),
            Time::ZERO,
        )
    };

    // 9 mm is inside the hot cell, 11 mm is inside the empty one, and they are 2 mm apart.
    assert!(at(9.0).is_finite(), "9 mm is still metal");
    assert!(at(11.0).is_nan(), "11 mm is already nothing");

    // Exactly on the face, the round is to the lower index, which is the solid one here. What
    // matters is that it is *one or the other* rather than a blend of both.
    let face = at(10.0);
    assert!(
        face.is_nan() || (face - Temperature::celsius(300.0).to_si()).abs() < 1e-12,
        "a face sample is one cell or the other, not an average: {face:.4} K"
    );
}

/// **A block with no void samples exactly as it did**, which is what says this is a mask and not a
/// new interpolation.
///
/// The renormalised weights have to collapse to plain trilinear when every corner counts, or every
/// scene in the repository would move by a little and the determinism digest would say so without
/// saying why. Checked against the trilinear form written out here from the cell values.
#[test]
fn a_solid_block_interpolates_exactly_as_before() {
    let mut block = Solid3D::new(
        "solid",
        Substance::copper(),
        (2, 2, 2),
        Length::mm(10.0),
        Temperature::celsius(20.0),
    );
    // Eight distinct values, so any weight that went missing shows up.
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..2 {
                let c = 100.0 * (i + 2 * j + 4 * k) as f64;
                block.set_temperature(i, j, k, Temperature::celsius(c));
            }
        }
    }

    // A point at 30% along each axis between the two cell centres, which are at 5 and 15 mm.
    let (fx, fy, fz) = (0.3, 0.6, 0.8);
    let mm = |f: f64| (5.0 + 10.0 * f) * 1e-3;
    let got = block.at(
        LengthVec::from_si(glam::DVec3::new(mm(fx), mm(fy), mm(fz))),
        Time::ZERO,
    );

    let cell = |i: usize, j: usize, k: usize| {
        Temperature::celsius(100.0 * (i + 2 * j + 4 * k) as f64).to_si()
    };
    let lerp = |lo: f64, hi: f64, t: f64| lo * (1.0 - t) + hi * t;
    let z0 = lerp(
        lerp(cell(0, 0, 0), cell(1, 0, 0), fx),
        lerp(cell(0, 1, 0), cell(1, 1, 0), fx),
        fy,
    );
    let z1 = lerp(
        lerp(cell(0, 0, 1), cell(1, 0, 1), fx),
        lerp(cell(0, 1, 1), cell(1, 1, 1), fx),
        fy,
    );
    let want = lerp(z0, z1, fz);

    assert!(
        (got - want).abs() < 1e-9,
        "with nothing masked the mask must be invisible: {got:.9} against {want:.9}"
    );
}

/// **The derivatives know about void too**, which `at` did and its three siblings did not.
///
/// `gradient`, `laplacian` and `rate` are the rest of the `ScalarField` surface, they are what the
/// SDK's own examples call, and all three read the raw cell array through `mirrored`. So a
/// consumer asking for the slope at the edge of a clearance was handed the difference between
/// their metal and a frozen number from before the cell was emptied — a plausible gradient across
/// a gap that conducts nothing.
///
/// The closed form is the boundary condition: a face with nothing behind it carries no heat, which
/// is the same statement as the block's own insulated outer face. So the difference is **one
/// sided**, and the check is that a cell against a clearance reports the same slope as the same
/// cell against the edge of the grid — geometry that differs in every way except the one that
/// matters.
#[test]
fn a_gradient_at_a_clearance_is_the_one_a_block_face_gives() {
    let build = |n: usize, gap: bool| {
        let mut block = Solid3D::new(
            "bar",
            Substance::copper(),
            (n, 1, 1),
            Length::mm(10.0),
            Temperature::celsius(20.0),
        );
        for i in 0..n {
            block.set_temperature(i, 0, 0, Temperature::celsius(20.0 + 30.0 * i as f64));
        }
        if gap {
            block.empty(|i, _, _| i >= 3)
        } else {
            block
        }
    };
    let slope_at = |block: &Solid3D, i: usize| {
        block
            .gradient(
                LengthVec::from_si(glam::DVec3::new((i as f64 + 0.5) * 10e-3, 5e-3, 5e-3)),
                Time::ZERO,
                Length::mm(10.0),
            )
            .x
    };

    // Cell 2 is the last cell either way: against the edge of a three-cell block, and against a
    // clearance in a six-cell one.
    let edge = slope_at(&build(3, false), 2);
    let clearance = slope_at(&build(6, true), 2);
    assert!(
        (edge - clearance).abs() < 1e-12,
        "a clearance is a boundary: {edge:e} K/m at a face against {clearance:e} K/m at a gap"
    );

    // And it is the one-sided value, computed here: (T2 - T1)/(2 dx), because the mirror puts T2
    // on both sides of the stencil and the 2 dx denominator stays.
    let closed = 30.0 / (2.0 * 10e-3);
    assert!(
        (clearance - closed).abs() < 1e-9,
        "one-sided over the full stencil width: {clearance:e} against {closed:e} K/m"
    );
}

/// **Inside a clearance there is no slope, no curvature and no rate**, rather than three numbers.
///
/// `rate` is the one that was quietly confident: `mobility` is `1/C`, an empty cell has no `C`, so
/// the product was `inf * 0.0` and the function's own `is_finite` guard turned that `NaN` into
/// **zero** — a cell that is not warming, which is exactly what a solid cell in equilibrium reports.
#[test]
fn nothing_has_no_derivatives_either() {
    let block = bar();
    let inside = LengthVec::from_si(glam::DVec3::new(15e-3, 5e-3, 5e-3));
    assert!(
        block
            .gradient(inside, Time::ZERO, Length::mm(10.0))
            .is_nan(),
        "no slope where there is nothing"
    );
    assert!(
        block
            .laplacian(inside, Time::ZERO, Length::mm(10.0))
            .is_nan(),
        "no curvature either"
    );
    assert!(
        block.rate(inside, Time::ZERO, Time::from_si(1e-6)).is_nan(),
        "and it is not warming at zero, it is not warming"
    );

    // The metal beside it still answers, or the guard has swallowed the block as well.
    let metal = LengthVec::from_si(glam::DVec3::new(5e-3, 5e-3, 5e-3));
    assert!(block
        .gradient(metal, Time::ZERO, Length::mm(10.0))
        .is_finite());
    assert!(block
        .rate(metal, Time::ZERO, Time::from_si(1e-6))
        .is_finite());
}

/// **A solid block's derivatives are untouched**, the same claim as for `at` and for the same
/// reason: this is a mask, and every scene in the repository would otherwise move by a little.
#[test]
fn a_solid_block_differentiates_exactly_as_before() {
    let mut block = Solid3D::new(
        "solid",
        Substance::copper(),
        (5, 1, 1),
        Length::mm(10.0),
        Temperature::celsius(20.0),
    );
    for i in 0..5 {
        block.set_temperature(i, 0, 0, Temperature::celsius(20.0 + 30.0 * i as f64));
    }
    let at = |mm: f64| LengthVec::from_si(glam::DVec3::new(mm * 1e-3, 5e-3, 5e-3));

    // A linear ramp: the interior slope is the ramp's own, and the curvature is zero.
    let slope = block.gradient(at(25.0), Time::ZERO, Length::mm(10.0)).x;
    assert!(
        (slope - 30.0 / 10e-3).abs() < 1e-9,
        "30 K per 10 mm cell: {slope:e} K/m"
    );
    assert!(
        block
            .laplacian(at(25.0), Time::ZERO, Length::mm(10.0))
            .abs()
            < 1e-6,
        "a straight line has no curvature"
    );
}
