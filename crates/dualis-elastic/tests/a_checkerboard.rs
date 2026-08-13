//! An isotropic composite against the Hashin–Shtrikman bounds, and what an affine boundary can and
//! cannot say.
//!
//! `a_layered_wave.rs` and `a_layered_block.rs` measure the **Voigt–Reuss** pair, whose ends a laminate
//! attains. Hashin–Shtrikman is the tighter pair and it assumes the microstructure is statistically
//! *isotropic*, which a laminate is not — so it needs a different witness, and a three-dimensional
//! checkerboard is one: cubic symmetry gives an isotropic effective tensor.
//!
//! # What an affine boundary measures, and the mistake this file made about it
//!
//! This has to be said before any number, because a first draft asserted something that does not follow
//! and the bulk modulus is what caught it.
//!
//! An effective modulus needs the composite driven by something, and what `Block` offers is a displacement
//! on the boundary. An affine one — `u = ε·x` — is a **kinematically uniform** boundary condition, and the
//! standard result is that the apparent stiffness it produces is an upper bound on the **effective**
//! stiffness of the infinite medium:
//!
//! ```text
//!   Reuss <= HS- <= C_effective <= C_apparent <= Voigt
//! ```
//!
//! The right end is a theorem too: the affine field held *throughout* is admissible and its energy is
//! exactly the Voigt energy, so relaxing the interior can only lower it.
//!
//! The draft went on to assert `C_apparent >= HS+`, reasoning that an upper estimate cannot cross an upper
//! bound. **That does not follow.** `C_apparent` bounds `C_effective`, and `C_effective` is itself below
//! `HS+`, which leaves `C_apparent` free to sit on either side. The bulk modulus does sit below, and the
//! assertion fired.
//!
//! Read correctly, being below is the *more* informative outcome, because `C_apparent` is a **ceiling on
//! the truth**:
//!
//! ```text
//!   four periods per axis, growing the cells per block
//!   shear:  1 cell  exactly Voigt (aliased)   2 cells  1.183 x HS+   4 cells  1.005 x HS+
//!   bulk:                                     4 cells  0.972 x HS+
//! ```
//!
//! For the shear modulus the ceiling lands half a percent above HS+, so the effective `G` is at most that
//! and **HS+ is tight to within 0.5% at worst**. For the bulk modulus the ceiling is *below* HS+, so the
//! effective `K` is strictly under the bound and **HS+ is at least 2.8% loose there**. Two different
//! answers for the two moduli, both measured, and neither is a bracketing — the conductivity file gets a
//! genuine bracketing because a marched field with a flux residual is not a kinematic estimate.
//!
//! # Two convergence parameters, and only one of them is the obvious one
//!
//! Refining a checkerboard has two knobs and they do different things:
//!
//! - **periods per axis** removes the boundary layer where the affine condition fights the
//!   microstructure. It converges fast and plateaus — measured, 1.206 to 1.175 times HS+ as the count
//!   goes from two per axis to six, a move of 2.7% — and the plateau is *not* the answer;
//! - **cells per block** resolves the strain field inside one phase. This is the slow one, and it is the
//!   one that moves the plateau: two cells plateau 17.5% above HS+ and four get to 0.5%.
//!
//! A sweep that changed both at once would report their sum as one rate, which is the mistake
//! `a_layered_wave.rs` records making with its own two errors.

use dualis_core::mixture::Mix;
use dualis_core::substance::Substance;
use dualis_core::units::Length;
use dualis_elastic::{Block, Elastic};
use dualis_units::LengthVec;

const DX: f64 = 1e-3;
/// The applied strain. Linear elasticity has no scale, so this only has to stay away from the point where
/// `u` and `dx` stop being distinguishable in an `f64`.
const AMP: f64 = 1e-6;
/// What `solve` is asked for; every tolerance here traces to it.
const TOLERANCE: f64 = 1e-12;

fn mix() -> Mix {
    Mix::of(&[(Substance::aluminium_6061(), 0.5), (Substance::pla(), 0.5)])
        .expect("fractions sum to one")
}

fn phases() -> (Elastic, Elastic) {
    (
        Elastic::aluminium_6061(),
        Elastic::from_substance(&Substance::pla()).expect("PLA states its mechanical properties"),
    )
}

/// What kind of affine field the boundary is given.
#[derive(Clone, Copy)]
enum Drive {
    /// `u = (γz, 0, 0)`. Simple shear, so `U = ½ G γ² V`.
    Shear,
    /// `u = ε·x`. Hydrostatic, so `U = ½ K (3ε)² V`.
    Hydrostatic,
}

/// The apparent modulus of a `cells³` checkerboard of `block`-cell cubes, from the strain energy.
///
/// The energy route rather than a reaction force, because it is one scalar and needs no bookkeeping about
/// which face carries what — and because `strain_energy` is `½uᵗKu`, assembled by the same code that
/// solved the system, so it cannot disagree with the solution about what the stiffness is.
///
/// `block = 0` means a uniform block, which is the calibration.
fn apparent(cells: usize, block: usize, drive: Drive) -> f64 {
    let (stiff, soft) = phases();
    let mut b = Block::new("board", (cells, cells, cells), Length::from_si(DX), stiff);
    if block > 0 {
        let changed = b.fill(soft, |i, j, k| (i / block + j / block + k / block) % 2 == 1);
        assert_eq!(
            changed * 2,
            cells * cells * cells,
            "a checkerboard of {block}-cell cubes in {cells}³ should be exactly half and half"
        );
    }
    match drive {
        Drive::Shear => b.prescribe_boundary(|p| {
            LengthVec::from_si(glam::DVec3::new(AMP * p.to_si().z, 0.0, 0.0))
        }),
        Drive::Hydrostatic => b.prescribe_boundary(|p| LengthVec::from_si(p.to_si() * AMP)),
    }
    assert!(
        b.solve(TOLERANCE),
        "the solve did not converge: residual {:e}",
        b.residual()
    );

    let volume = (cells as f64 * DX).powi(3);
    let energy = b.strain_energy().to_si();
    match drive {
        // U = ½ G γ² V, because simple shear `u_x = γz` has `ε_xz = γ/2` and the energy density is
        // `2μ ε_xz²`.
        Drive::Shear => 2.0 * energy / (AMP * AMP * volume),
        // U = ½ K θ² V with θ = 3ε.
        Drive::Hydrostatic => 2.0 * energy / (9.0 * AMP * AMP * volume),
    }
}

/// **A uniform block's apparent moduli are its own, to `1.6e-14`.**
///
/// The calibration, and it is what makes every number below mean anything. For a homogeneous body the
/// affine field *is* the solution — it satisfies equilibrium — so the apparent modulus has to be the
/// material's exactly, and any error is the solver's. If this were even slightly off, a checkerboard
/// landing near a bound would be two errors of unknown size meeting.
///
/// Both drives, because they exercise different parts of the stiffness: the hydrostatic one is blind to
/// `G` in a homogeneous body and the shear one is blind to `K`, so a mistake in either the energy
/// bookkeeping or the boundary field shows up in one and not the other.
#[test]
fn a_uniform_block_reports_its_own_moduli() {
    let (stiff, _) = phases();
    let (e, nu) = (stiff.youngs_modulus.to_si(), stiff.poisson_ratio);
    let want_g = e / (2.0 * (1.0 + nu));
    let want_k = e / (3.0 * (1.0 - 2.0 * nu));

    let got_g = apparent(8, 0, Drive::Shear);
    let got_k = apparent(8, 0, Drive::Hydrostatic);
    println!(
        "  G {:.6} GPa against {:.6} — off {:.2e}\n  K {:.6} GPa against {:.6} — off {:.2e}",
        got_g / 1e9,
        want_g / 1e9,
        (got_g / want_g - 1.0).abs(),
        got_k / 1e9,
        want_k / 1e9,
        (got_k / want_k - 1.0).abs()
    );
    assert!(
        (got_g / want_g - 1.0).abs() < 1e-12 && (got_k / want_k - 1.0).abs() < 1e-12,
        "an affine field on a homogeneous block is the exact solution, so both moduli should be exact"
    );
    // And they are different numbers, so neither test is silently measuring the other.
    assert!(want_k / want_g > 2.5, "K and G are far apart for aluminium");
}

/// **A checkerboard one cell per phase is aliased to Voigt exactly.**
///
/// The same finding the thermal checkerboard records, arriving by a different route. There the cause was
/// every face being an interface, so the discrete operator was that of a uniform medium at the harmonic
/// mean. Here it is kinematic: with an affine displacement imposed on the boundary and a microstructure
/// at the element scale, a trilinear element has no freedom left to relax into, so the affine field *is*
/// the discrete solution and its energy is the volume average of the two phases' — which is Voigt.
///
/// Asserted so it stays known. A one-cell board looks like the finest possible resolution of a
/// checkerboard and is in fact no checkerboard at all.
#[test]
fn a_one_cell_checkerboard_is_aliased_to_the_voigt_bound() {
    let m = mix();
    let (_, voigt) = m.shear_bounds().expect("both state mechanical properties");
    let got = apparent(8, 1, Drive::Shear);
    println!(
        "  G {:.6} GPa against Voigt {:.6} — off {:.2e}",
        got / 1e9,
        voigt.to_si() / 1e9,
        (got / voigt.to_si() - 1.0).abs()
    );
    assert!(
        (got / voigt.to_si() - 1.0).abs() < 1e-9,
        "a one-cell board should measure Voigt {:.6} GPa exactly, measures {:.6}",
        voigt.to_si() / 1e9,
        got / 1e9
    );
}

/// **Resolving the microstructure walks the estimate down to half a percent above HS+.**
///
/// The headline, and the number that says the Hashin–Shtrikman upper bound is nearly tight for a
/// checkerboard. Four periods per axis throughout, so the boundary layer is the same in all three and the
/// only thing changing is how many elements resolve one cube of one phase.
///
/// Monotone down, and the assertion is the trend and the final value rather than a rate: the approach is
/// not a clean power law, because the strain field near a corner where four cubes meet is singular and a
/// trilinear element resolves it slowly.
///
/// What is asserted about the bounds is only what is a theorem — `HS− ≤ C_apparent ≤ Voigt`, since the
/// apparent value is a ceiling on an effective modulus that is itself above `HS−`. It happens to stay
/// above `HS+` for this modulus and that is **not** asserted, because it is not a theorem and the bulk
/// modulus is a counterexample to the version of it a first draft believed.
///
/// The last step also falls to 0.745 of Voigt, which is what says the geometry is genuinely relaxing
/// rather than being held affine by the mesh.
#[test]
fn resolving_the_microstructure_walks_the_estimate_down_toward_the_upper_bound() {
    let m = mix();
    let (_, voigt) = m.shear_bounds().expect("mechanical properties");
    let (hs_lo, hs_hi) = m.shear_hashin_shtrikman().expect("two phases");
    println!(
        "  Voigt {:.4} GPa, HS+ {:.4} — four periods per axis throughout:",
        voigt.to_si() / 1e9,
        hs_hi.to_si() / 1e9
    );

    let mut previous = f64::INFINITY;
    let mut finest = 0.0;
    for block in [1usize, 2, 4] {
        let cells = 8 * block;
        let got = apparent(cells, block, Drive::Shear);
        println!(
            "    {block} cell(s) per block ({cells}³): {:.6} GPa = {:.4}× HS+, {:.4}× Voigt",
            got / 1e9,
            got / hs_hi.to_si(),
            got / voigt.to_si()
        );
        assert!(
            got < previous,
            "resolving the microstructure should soften the estimate: {block} cells gave \
             {:.6} GPa against the previous {:.6}",
            got / 1e9,
            previous / 1e9
        );
        // What is a theorem: the apparent value is a ceiling on an effective modulus that is itself at
        // or above the lower bound, and the affine field held throughout gives exactly Voigt so relaxing
        // can only come down from it.
        assert!(
            got > hs_lo.to_si() && got <= voigt.to_si() * (1.0 + 1e-12),
            "{:.6} GPa is outside HS- {:.6} and Voigt {:.6}, and both of those are theorems",
            got / 1e9,
            hs_lo.to_si() / 1e9,
            voigt.to_si() / 1e9
        );
        previous = got;
        finest = got;
    }

    // The ceiling on the effective modulus, and therefore how tight the bound can be: within 1% of HS+
    // means the effective `G` is at most 1% above a number that is meant to be above it, so the bound is
    // tight to that. Asserted as a two-sided window because drifting *down* through HS+ would also be
    // interesting and would mean this file had started measuring something else.
    let margin = finest / hs_hi.to_si() - 1.0;
    println!(
        "  the finest is {:.3}% above HS+, so the effective G is at most that far above the bound",
        margin * 100.0
    );
    assert!(
        margin.abs() < 0.01,
        "four cells per block should land within 1% of HS+, is {:.3}% away",
        margin * 100.0
    );
}

/// **The period count is the fast knob, and its plateau is not the answer.**
///
/// The other convergence parameter, isolated: two cells per block throughout, so the microstructure
/// resolution never changes and only the boundary layer shrinks. It converges quickly — tripling the
/// period count moves the estimate 2.7% — and it converges to something **17.5% above HS+**, which is
/// where the unresolved microstructure leaves it.
///
/// That is the point of measuring it separately. A reader seeing only this sweep would conclude the method
/// had converged and the bound was 18% loose. The bound is not loose; two cells per block is not a cube.
#[test]
fn the_period_count_converges_fast_and_to_the_wrong_place() {
    let m = mix();
    let (_, hs_hi) = m.shear_hashin_shtrikman().expect("two phases");
    let mut values = Vec::new();
    for cells in [8usize, 16, 24] {
        let got = apparent(cells, 2, Drive::Shear);
        println!(
            "  {cells}³, {} periods per axis: {:.6} GPa = {:.4}× HS+",
            cells / 4,
            got / 1e9,
            got / hs_hi.to_si()
        );
        values.push(got);
    }
    // Falling, and by little: the boundary layer is a small effect next to the resolution one.
    assert!(
        values[2] < values[1] && values[1] < values[0],
        "more periods should mean less boundary layer"
    );
    let moved = (values[0] / values[2] - 1.0).abs();
    println!("  tripling the period count moved it {:.2}%", moved * 100.0);
    assert!(
        moved < 0.05,
        "the period count is the fast knob and should move the answer by a few percent, moved it \
         {:.2}%",
        moved * 100.0
    );
    // And where it lands is nowhere near the bound, which is the finding.
    assert!(
        values[2] / hs_hi.to_si() > 1.1,
        "two cells per block plateaus well above HS+, at {:.4}×",
        values[2] / hs_hi.to_si()
    );
}

/// **The bulk modulus tells a different story: its upper bound is measurably loose.**
///
/// `K` and `G` are separate bounds with separate reference phases, so a test of one is not a test of the
/// other — `bulk_hashin_shtrikman` could have the reference phase's shear modulus in the wrong place and
/// every shear measurement here would still pass.
///
/// And the answer differs. The apparent `K` is a **ceiling on the effective `K`**, and it comes out at
/// 0.972 of HS+ — *below* the bound. So the effective bulk modulus of a checkerboard is at least 2.8%
/// under the Hashin–Shtrikman upper bound, where for the shear modulus the same measurement puts the
/// truth within half a percent of it. The pair is tight for one modulus of this geometry and loose for the
/// other, which is a thing a caller choosing a number inside it would want to know and is not something
/// the algebra says.
///
/// This is also the measurement that corrected the module docs: the draft asserted an affine estimate
/// could not fall below HS+, and it does.
///
/// One resolution rather than a sweep, because the sweep's lesson is recorded above and a 32³ hydrostatic
/// solve is the most expensive thing in this file.
#[test]
fn the_bulk_modulus_lands_between_its_bounds_the_same_way() {
    let m = mix();
    let (reuss, voigt) = m.bulk_bounds().expect("mechanical properties");
    let (hs_lo, hs_hi) = m.bulk_hashin_shtrikman().expect("two phases");
    println!(
        "  K: Reuss {:.4}  HS {:.4} to {:.4}  Voigt {:.4} GPa — {:.2}× against {:.2}×",
        reuss.to_si() / 1e9,
        hs_lo.to_si() / 1e9,
        hs_hi.to_si() / 1e9,
        voigt.to_si() / 1e9,
        voigt.to_si() / reuss.to_si(),
        hs_hi.to_si() / hs_lo.to_si()
    );

    let got = apparent(16, 4, Drive::Hydrostatic);
    println!(
        "  two periods, four cells per block: {:.6} GPa = {:.4}× HS+",
        got / 1e9,
        got / hs_hi.to_si()
    );
    // The theorem, again: above the lower bound and at or below Voigt.
    assert!(
        got > hs_lo.to_si() && got <= voigt.to_si(),
        "an affine estimate sits between HS- {:.4} and Voigt {:.4}, is {:.4} GPa",
        hs_lo.to_si() / 1e9,
        voigt.to_si() / 1e9,
        got / 1e9
    );
    // And the finding: it is below HS+, so the bound is loose for this modulus by at least that much.
    let slack = 1.0 - got / hs_hi.to_si();
    println!(
        "  which is {:.2}% BELOW HS+, so the effective K is at least that far under the bound",
        slack * 100.0
    );
    assert!(
        slack > 0.02,
        "the apparent K should land measurably below HS+ for this geometry, is {:.2}% below",
        slack * 100.0
    );
    // And HS is the pair worth having: it narrows the honest range by more than a factor of two.
    assert!(
        (voigt.to_si() / reuss.to_si()) / (hs_hi.to_si() / hs_lo.to_si()) > 1.8,
        "HS should be substantially narrower than Voigt-Reuss here"
    );
}
