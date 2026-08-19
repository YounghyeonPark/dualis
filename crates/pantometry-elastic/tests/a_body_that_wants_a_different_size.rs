//! A stress-free strain, and the four closed forms that catch four different ways of getting it
//! wrong.
//!
//! An **eigenstrain** is the size a body would take if nothing held it: `αΔT` for thermal
//! expansion, and the same statement for swelling, curing shrinkage or a phase change. The domain
//! is deliberately not told which, so it never has to depend on whatever computes the temperature.
//!
//! Trilinear elements reproduce any linear displacement field exactly, and a uniform eigenstrain
//! produces one — so **every case here is exact at any mesh size**, not converging. That is what
//! makes the tolerances machine epsilon rather than a judgement, and it is why a wrong factor
//! cannot hide behind a coarse grid.

use pantometry_core::units::{Density, Length, Pressure};
use pantometry_elastic::{Block, Elastic, Face};

/// Aluminium-like, with round numbers so every closed form below can be checked by hand.
fn metal() -> Elastic {
    Elastic {
        youngs_modulus: Pressure::from_si(70e9),
        poisson_ratio: 0.33,
        density: Density::from_si(2700.0),
    }
}

fn block(n: usize) -> Block {
    Block::new("bar", (n, n, n), Length::mm(2.0), metal())
}

/// **A body free to expand takes exactly the size the eigenstrain asks for, and carries no
/// stress.**
///
/// The first form, and the one that catches a sign or a factor immediately: with nothing holding
/// it, `ε = ε₀` and `σ = D(ε − ε₀) = 0`. The strain energy is the sharper half of the claim — a
/// body that expanded by the right amount while storing energy would have the displacement right
/// and the stress wrong, and only the energy would say so.
#[test]
fn a_free_body_takes_the_size_it_asks_for_and_stores_nothing() {
    let strain = 1.7e-3; // e.g. 23e-6 per kelvin over 74 K
    let mut b = block(4);
    // Held against rigid-body motion only: one corner pinned, and two rollers to stop rotation.
    b.roller(Face::XLow);
    b.roller(Face::YLow);
    b.roller(Face::ZLow);
    b.stress_free_strain(|_, _, _| strain);
    assert!(b.solve(1e-12), "a free expansion should solve");

    for axis in 0..3 {
        let got = b.mean_strain(axis);
        assert!(
            (got - strain).abs() < 1e-12,
            "axis {axis}: free expansion is the eigenstrain itself, {got:e} against {strain:e}"
        );
    }
    // Three linear strains of `e` make a volumetric strain of `3e` to first order, and the
    // domain's own measure should agree to the second-order term it keeps.
    assert!(
        (b.volumetric_strain() - 3.0 * strain).abs() < 1e-8,
        "volumetric strain is three times the linear one: {:e}",
        b.volumetric_strain()
    );
    assert!(
        b.strain_energy().to_si().abs() < 1e-9,
        "an unresisted expansion stores no energy: {:e} J",
        b.strain_energy().to_si()
    );
}

/// **A body held in every direction develops `σ = −E·ε₀/(1−2ν)` in every direction.**
///
/// The second form, and it is a different combination of `λ` and `μ` from the first: with `ε = 0`,
/// `σ = −D ε₀ = −(3λ+2μ)ε₀`, and `3λ+2μ = 3K = E/(1−2ν)`. A body that had the two Lamé constants
/// transposed passes the free-expansion test and fails this one.
///
/// Read as the reaction on a held face, which is the force the clamp has to supply — the number a
/// designer actually wants, and one that comes from the assembled system rather than from a strain
/// the test could have arranged.
#[test]
fn a_fully_held_body_pushes_back_with_three_k_epsilon() {
    let strain = 1.7e-3;
    let n = 4;
    let mut b = block(n);
    for face in [
        Face::XLow,
        Face::XHigh,
        Face::YLow,
        Face::YHigh,
        Face::ZLow,
        Face::ZHigh,
    ] {
        b.clamp(face);
    }
    b.stress_free_strain(|_, _, _| strain);
    assert!(b.solve(1e-12), "a held expansion should solve");

    let (e, nu) = (70e9, 0.33);
    let closed = -e * strain / (1.0 - 2.0 * nu);
    let area = (n as f64 * 2e-3).powi(2);
    let got = b.normal_reaction(Face::XHigh).to_si() / area;
    assert!(
        (got / closed - 1.0).abs() < 1e-9,
        "a fully held expansion pushes at -E e/(1-2v) = {closed:e} Pa, pushes {got:e} Pa"
    );
}

/// **A bar held along one axis and free on the others develops `σ = −E·ε₀` along it, and nothing
/// across it.**
///
/// The third form, and the classic one — a rail that cannot move develops `EαΔT` and buckles. It
/// is a *different* combination again: the free directions relieve the Poisson coupling entirely,
/// so the answer has no `ν` in it at all. A body whose held directions leaked into its free ones
/// would come out with a `ν` here and would pass both tests above.
#[test]
fn a_bar_held_along_one_axis_develops_e_epsilon_and_no_more() {
    let strain = 1.7e-3;
    let n = 4;
    let mut b = block(n);
    // **Rollers, not clamps.** A clamp holds all three components, so clamping the two x faces
    // would hold the transverse expansion at the ends as well and the answer would not be the bar
    // case at all — measured, it comes out 1.82 times `Ee`, which is a three-dimensional answer to
    // a different question. A roller holds only the normal component, which is what "held along
    // one axis" means.
    b.roller(Face::XLow);
    b.roller(Face::XHigh);
    // Free across, but pinned against sliding.
    b.roller(Face::YLow);
    b.roller(Face::ZLow);
    b.stress_free_strain(|_, _, _| strain);
    assert!(b.solve(1e-12), "a uniaxially held expansion should solve");

    let closed = -70e9 * strain;
    let area = (n as f64 * 2e-3).powi(2);
    let along = b.normal_reaction(Face::XHigh).to_si() / area;
    assert!(
        (along / closed - 1.0).abs() < 1e-9,
        "held along one axis the stress is -E e = {closed:e} Pa, is {along:e} Pa"
    );
    // And across it the body simply grew: the free axes carry the eigenstrain plus what the
    // held axis's Poisson coupling adds, `e·(1 + v·... )` — asserted as *positive and free of
    // reaction* rather than as a formula, because the formula is the previous claim again.
    assert!(
        b.normal_reaction(Face::YHigh).to_si().abs() < 1e-3,
        "a free face carries no reaction: {:e} N",
        b.normal_reaction(Face::YHigh).to_si()
    );
    assert!(
        b.mean_strain(1) > strain,
        "the free axis grows by more than the eigenstrain, because the held one squeezes it: {:e}",
        b.mean_strain(1)
    );
}

/// **Two halves that want different sizes fight, and the answer converges on force balance from
/// the stiff side.**
///
/// The fourth form, and the one a real assembly is made of: a bimaterial strip. Nothing is clamped
/// — the stress is entirely internal, generated by the mismatch — so it is the case a test that
/// only ever clamped things would never see, and it is the mechanism behind solder fatigue in the
/// power module scene `pantometry-world` runs.
///
/// **This one is not exact, and the crate's own docs say why before it is measured.** A uniform
/// eigenstrain gives a linear displacement field, which a trilinear element reproduces exactly; a
/// *mismatched* one makes the strip bend, and "bending only converges, and from the stiff side —
/// a fully integrated trilinear element develops shear strain where it should be flexing". So the
/// claim here is the convergence rather than a tolerance somebody chose.
///
/// The limit is force balance: with free ends the net axial force is zero, and for two equal
/// halves of the same modulus that puts the mean strain at the average of what they wanted.
/// Measured errors `−1.91e-5, −1.12e-5, −5.99e-6` at 4, 8 and 16 elements — ratios 1.70 and 1.87,
/// approaching two — and **negative at every resolution**, which is the stiff side the doc names.
#[test]
fn two_halves_that_want_different_sizes_converge_on_force_balance() {
    let (hot, cold) = (2.0e-3, 0.0);
    let average = 0.5 * (hot + cold);

    let error = |n: usize| {
        let mut b = Block::new(
            "strip",
            (n, n, n),
            Length::from_si(16e-3 / n as f64),
            metal(),
        );
        b.roller(Face::XLow);
        b.roller(Face::YLow);
        b.roller(Face::ZLow);
        // The lower half wants to be bigger; the upper half wants to stay put.
        b.stress_free_strain(|_, _, k| if k < n / 2 { hot } else { cold });
        assert!(
            b.solve(1e-13),
            "a mismatched strip should solve at {n} elements"
        );
        assert!(
            b.strain_energy().to_si() > 0.0,
            "a mismatch stores strain energy, unlike a free expansion: {:e} J",
            b.strain_energy().to_si()
        );
        assert!(
            b.mean_strain(0) < hot && b.mean_strain(0) > cold,
            "neither half gets what it asked for: {:e}",
            b.mean_strain(0)
        );
        b.mean_strain(0) - average
    };

    let (coarse, middle, fine) = (error(4), error(8), error(16));
    println!(
        "  errors {coarse:+.4e} {middle:+.4e} {fine:+.4e}, ratios {:.2} {:.2}",
        coarse / middle,
        middle / fine
    );
    for e in [coarse, middle, fine] {
        assert!(
            e < 0.0,
            "trilinear elements are too stiff in bending, so the mean is short of the average:              {e:e}"
        );
    }
    for (a, b) in [(coarse, middle), (middle, fine)] {
        let ratio = a / b;
        assert!(
            (1.6..2.3).contains(&ratio),
            "first order is a ratio of two per doubling: {ratio:.3}"
        );
    }
    assert!(
        middle / fine > coarse / middle,
        "and it should be approaching two rather than sitting still: {:.3} then {:.3}",
        coarse / middle,
        middle / fine
    );
}

/// **A uniform eigenstrain applies no net force**, which is the identity the load assembly has to
/// satisfy before any of the forms above can be right.
///
/// A body expanding into nothing does not push itself across the room. Checked on the free case by
/// the reactions summing to zero over all six faces, which is the assembled statement rather than
/// the per-element one.
#[test]
fn expanding_into_nothing_pushes_nothing() {
    let mut b = block(3);
    b.roller(Face::XLow);
    b.roller(Face::YLow);
    b.roller(Face::ZLow);
    b.stress_free_strain(|_, _, _| 5e-3);
    assert!(b.solve(1e-12));

    let mut total = 0.0;
    for face in [
        Face::XLow,
        Face::XHigh,
        Face::YLow,
        Face::YHigh,
        Face::ZLow,
        Face::ZHigh,
    ] {
        total += b.normal_reaction(face).to_si().abs();
    }
    assert!(
        total < 1e-6,
        "an unresisted expansion is self-equilibrated: {total:e} N of reaction"
    );
}

/// **A block with no eigenstrain behaves exactly as it did**, which is what says this is an
/// addition and not a change.
///
/// Every scene and every test in this workspace predates the field, so the untouched path has to be
/// bit-for-bit what it was. Checked against the uniaxial modulus, which is one of the four the
/// crate's own docs say a trilinear element reproduces exactly.
#[test]
fn a_block_with_no_eigenstrain_is_untouched() {
    let n = 4;
    let mut b = block(n);
    b.roller(Face::XLow);
    b.roller(Face::YLow);
    b.roller(Face::ZLow);
    b.pull(Face::XHigh, Pressure::from_si(1e6));
    assert!(b.solve(1e-12));

    assert_eq!(b.stress_free_at(0, 0, 0), 0.0);
    let got = 1e6 / b.mean_strain(0);
    assert!(
        (got / 70e9 - 1.0).abs() < 1e-9,
        "uniaxial stress over strain is E: {got:e} against 70e9"
    );
}
