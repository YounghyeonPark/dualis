//! A laminate's stiffness, measured statically instead of with a wave.
//!
//! `a_layered_wave.rs` checks `Mix`'s stiffness bounds against Backus averaging by marching a wave and
//! reading a frequency. It gets 0.035% at its finest mesh on the overlapping modulus and costs 68
//! seconds in debug, and both of
//! those are the price of *time being in the problem*: a leapfrog has a dispersion relation, a period has
//! to be fitted over several cycles, and the mesh error is second order.
//!
//! The same modulus falls out of a static solve with none of that. `Block` is elliptic — there is no
//! step, no dispersion and no marching — so a laminate's harmonic modulus comes out to **solver
//! tolerance**. Measured, `4.8e-13` against the wave's `3.5e-4`: **nine orders** sharper. And 0.13 s in
//! debug against 67.5 for the wave file, which is **520× cheaper** for the overlapping claim.
//!
//! ```text
//!   pull a laminate column across its layers   ->  <1/(lambda+2mu)>^-1  = 11.1237 GPa
//! ```
//!
//! # Why a traction and not a prescribed displacement
//!
//! The harmonic mean is the **uniform-stress** average, so the way to realise it is to apply a stress and
//! let the strains come out however they will. A column with a traction on one end and equilibrium
//! everywhere has `σ_zz` uniform through the stack by construction — one dimension, no body force — and
//! each layer then contributes `σ/M_i` of strain. The compliances add, which is what makes the mean
//! harmonic. Prescribing the displacement instead measures the same thing the long way round, since the
//! interior still redistributes to make the stress uniform.
//!
//! Rollers on all four sides and a single element across `x` and `y`, so **every** node is on both an
//! x face and a y face and the lateral displacement is zero *pointwise* rather than only on the boundary.
//! That is what makes it uniaxial strain and therefore `λ+2μ` rather than `E`, and for aluminium those
//! differ by 48%.
//!
//! # What is not here
//!
//! The **arithmetic** mean, and for the reason `a_layered_wave.rs` gives about `C11`: realising it needs
//! the lateral strain zero on average while free locally, and `roller` constrains a face rather than a
//! plane through the interior. Statics gets the sharp half of the pair; the wave gets both halves less
//! sharply. Together they cover it, and neither on its own does.

use dualis_core::mixture::Mix;
use dualis_core::substance::Substance;
use dualis_core::units::{Length, Pressure};
use dualis_elastic::{Block, Elastic, Face};

/// The same pair `a_layered_wave.rs` uses: aluminium and PLA, a 17-fold contrast in the constrained
/// modulus, so the harmonic mean sits nowhere near either constituent and nowhere near their average.
fn phases() -> (Elastic, Elastic) {
    (
        Elastic::aluminium_6061(),
        Elastic::from_substance(&Substance::pla()).expect("PLA states its mechanical properties"),
    )
}

fn mix() -> Mix {
    Mix::of(&[(Substance::aluminium_6061(), 0.5), (Substance::pla(), 0.5)])
        .expect("fractions sum to one")
}

const DX: f64 = 1e-3;
/// The applied stress. Small enough that linear elasticity is the right model and large enough that the
/// solve is nowhere near the noise floor — a strain of about `10⁻⁴` at these moduli.
const PULL: f64 = 1.0e6;
/// What `solve` is asked for. Every tolerance below traces to it, which is the same discipline
/// `a_composite.rs` applies to its marched residual.
const TOLERANCE: f64 = 1e-12;

/// A column of alternating layers, pulled along `z` with its sides on rollers, and the effective
/// constrained modulus that comes out.
fn constrained_modulus(elements: usize, thickness: usize) -> f64 {
    let (stiff, soft) = phases();
    let mut b = Block::new("laminate", (1, 1, elements), Length::from_si(DX), stiff);
    let changed = b.fill(soft, |_, _, e_z| (e_z / thickness) % 2 == 1);
    assert_eq!(
        changed * 2,
        elements,
        "alternating layers should be exactly half of {elements} elements, are {changed}"
    );
    assert_eq!(b.materials().len(), 2, "two materials in the palette");

    // Rollers on all four sides: one element across x and y, so every node is on both pairs of faces and
    // the lateral displacement is zero everywhere rather than only on the boundary.
    for face in [Face::XLow, Face::XHigh, Face::YLow, Face::YHigh] {
        b.roller(face);
    }
    b.clamp(Face::ZLow);
    b.pull(Face::ZHigh, Pressure::from_si(PULL));
    assert!(
        b.solve(TOLERANCE),
        "the solve did not converge: residual {:e}",
        b.residual()
    );

    // `σ = M ε`, with `ε` the mean strain along the pull. The mean is the right one: each layer strains
    // differently and the composite's strain is their volume average, which over equal-thickness layers
    // is what the end-to-end displacement over the length gives.
    PULL / b.mean_strain(2)
}

/// **A laminate pulled across its layers has exactly the harmonic mean of `λ+2μ`.**
///
/// The sharp version of `a_layered_wave.rs`'s `C33`, and the difference is what an elliptic problem buys:
/// worst `4.8e-13` here against `3.5e-4` there, at four resolutions, in 0.13 s of debug.
///
/// It is not quite an algebraic identity the way the thermal laminate's Reuss result is — that one is
/// exact because a harmonic face mean makes the discrete resistance chain equal the continuum one, and a
/// trilinear element has no equivalent. What is left is the conjugate-gradient tolerance, which is why
/// every number below traces to `TOLERANCE` rather than to a mesh.
///
/// Four resolutions because an exact result and a small convergent error look identical at one, and this
/// workspace has read a first-order error as a tolerance before. The error does grow with the element
/// count — `5.9e-15` at eight against `4.8e-13` at sixty-four — which is the conjugate gradient's
/// tolerance accumulating over more degrees of freedom and not a mesh error: a mesh error would not care
/// how many iterations it took. `1e-11` is twenty times the worst of those.
#[test]
fn a_laminate_pulled_across_its_layers_is_the_harmonic_constrained_modulus() {
    let m = mix();
    let (reuss, voigt) = m
        .p_wave_modulus_bounds()
        .expect("both state mechanical properties");
    println!(
        "  bounds: {:.4} to {:.4} GPa, a factor of {:.3}",
        reuss.to_si() / 1e9,
        voigt.to_si() / 1e9,
        voigt.to_si() / reuss.to_si()
    );

    let mut worst = 0.0f64;
    for elements in [8, 16, 32, 64] {
        let got = constrained_modulus(elements, 1);
        let off = (got / reuss.to_si() - 1.0).abs();
        println!(
            "  {elements:2} elements: {:.6} GPa against {:.6} — off {off:.2e}",
            got / 1e9,
            reuss.to_si() / 1e9
        );
        assert!(
            off < 1e-11,
            "{elements} elements: pulling across the layers should give the harmonic mean {:.6} GPa, \
             gives {:.6}",
            reuss.to_si() / 1e9,
            got / 1e9
        );
        worst = worst.max(off);
    }
    println!("  worst {worst:.2e}, and the solver was asked for {TOLERANCE:e}");

    // And it is nowhere near the arithmetic mean, which is the whole point of the bound being a range:
    // 4.85× apart, so a caller who took the average of the two constituents would be 385% out.
    assert!(
        voigt.to_si() / reuss.to_si() > 4.0,
        "the pair should be far apart for this contrast, is {:.3}x",
        voigt.to_si() / reuss.to_si()
    );
}

/// **The layer thickness does not matter, which is what says this is not a long-wavelength result.**
///
/// The wave measurement of the *arithmetic* mean pays an error that grows with layer thickness — 23% at
/// eight layers per wavelength — because it needs the layers to move together and that only happens when
/// a wavelength spans many of them. A static uniform stress needs nothing of the kind: the stress is
/// uniform because equilibrium says so, whatever the layers look like.
///
/// So the same 64-element column with layers 1, 2, 4 and 8 elements thick gives the same modulus to
/// `4.8e-13`, and the four numbers do not trend — 4.83, 4.35, 4.27 and 3.18 times `1e-13`, which is solver
/// noise and not a thickness dependence. That contrast between the two files is the honest statement of
/// when Backus averaging is a limit and when it is an identity, and it is a thing neither file could say
/// alone.
#[test]
fn the_layer_thickness_does_not_change_the_static_answer() {
    let m = mix();
    let (reuss, _) = m.p_wave_modulus_bounds().expect("mechanical properties");
    for thickness in [1usize, 2, 4, 8] {
        let got = constrained_modulus(64, thickness);
        let off = (got / reuss.to_si() - 1.0).abs();
        println!(
            "  layers {thickness} element(s) thick: {:.6} GPa — off {off:.2e}",
            got / 1e9
        );
        assert!(
            off < 1e-11,
            "layers {thickness} thick: a static uniform stress does not care how thick the layers \
             are, and this is off {off:.2e}"
        );
    }
}

/// **A uniform block gives the constrained modulus of its one material, and a filled-with-itself block is
/// unchanged.**
///
/// The two compatibility statements. The first is the calibration: if a uniform column did not give
/// `λ+2μ` then the laminate agreeing with a harmonic mean of `λ+2μ` would be two mistakes cancelling. It
/// is `λ+2μ` and not `E` because the rollers hold the lateral strain at zero, and for aluminium those
/// differ by **48%** — so a version of this test that had the boundary conditions wrong would miss by far
/// more than any tolerance.
///
/// The second can be an exact equality, because an equal material resolves to the palette entry it
/// already has rather than adding a second one.
#[test]
fn a_uniform_column_gives_its_own_constrained_modulus() {
    let (stiff, _) = phases();
    let mut b = Block::new("column", (1, 1, 32), Length::from_si(DX), stiff);
    for face in [Face::XLow, Face::XHigh, Face::YLow, Face::YHigh] {
        b.roller(face);
    }
    b.clamp(Face::ZLow);
    b.pull(Face::ZHigh, Pressure::from_si(PULL));
    assert!(b.solve(TOLERANCE), "the solve converged");
    let got = PULL / b.mean_strain(2);
    let want = stiff.constrained_modulus().to_si();
    println!(
        "  uniform: {:.6} GPa against λ+2μ = {:.6}, and E is {:.6}",
        got / 1e9,
        want / 1e9,
        stiff.youngs_modulus.to_si() / 1e9
    );
    assert!(
        (got / want - 1.0).abs() < 1e-11,
        "a uniform column gives its own constrained modulus: {got:e} against {want:e}"
    );
    // The 48% that says the boundary conditions are the ones intended.
    assert!(
        (want / stiff.youngs_modulus.to_si() - 1.482).abs() < 0.001,
        "λ+2μ over E is 1.482 for this material, is {:.4}",
        want / stiff.youngs_modulus.to_si()
    );

    // Filling with what it already held changes nothing, and un-converges nothing.
    let before = b.mean_strain(2);
    assert_eq!(b.fill(stiff, |_, _, _| true), 0, "nothing changed");
    assert_eq!(b.materials().len(), 1, "and the palette did not grow");
    assert!(b.converged(), "so the solution is still the solution");
    assert_eq!(b.mean_strain(2), before, "bit for bit");
}

/// **A fill invalidates the solution rather than leaving a stale one.**
///
/// The silent failure this guards against: `fill` changes the assembly, so the displacement field on
/// record solves the *previous* problem. Left alone it would still read as converged, `mean_strain` would
/// still return a number, and the number would be the answer to a question the caller has changed — which
/// is the shape of failure this workspace keeps finding, an answer that is present and wrong rather than
/// absent.
#[test]
fn filling_a_solved_block_invalidates_the_solution() {
    let (stiff, soft) = phases();
    let mut b = Block::new("column", (1, 1, 16), Length::from_si(DX), stiff);
    for face in [Face::XLow, Face::XHigh, Face::YLow, Face::YHigh] {
        b.roller(face);
    }
    b.clamp(Face::ZLow);
    b.pull(Face::ZHigh, Pressure::from_si(PULL));
    assert!(b.solve(TOLERANCE) && b.converged());

    let changed = b.fill(soft, |_, _, e_z| e_z % 2 == 1);
    assert_eq!(changed, 8);
    assert!(
        !b.converged(),
        "a block whose materials changed is not solved any more"
    );
    assert!(
        !b.residual().is_finite(),
        "and its residual says so rather than reporting the old one: {:e}",
        b.residual()
    );

    // Solving again gives the laminate's answer, which is the softer one.
    assert!(b.solve(TOLERANCE));
    let laminate = PULL / b.mean_strain(2);
    assert!(
        laminate < 0.5 * stiff.constrained_modulus().to_si(),
        "half PLA by volume should be far softer than aluminium: {:.4} GPa against {:.4}",
        laminate / 1e9,
        stiff.constrained_modulus().to_si() / 1e9
    );
}
