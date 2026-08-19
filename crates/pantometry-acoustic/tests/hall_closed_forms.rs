//! `Hall` against the rigid-wall mode frequencies, which are exact.
//!
//! The wave equation in three dimensions has the same trap the two-dimensional one shipped with:
//! a scheme that is first order where it claims to be second, and looks merely coarse. So the
//! checks here are on **rates** wherever a rate exists, and the one place a rate cannot help — the
//! energy — is checked as an exact equality instead, because rigid walls make it one.

use pantometry_acoustic::Hall;
use pantometry_core::units::{Length, Pressure, Time};
use pantometry_core::{Domain, Exchange};

fn hall(across: usize) -> Hall {
    Hall::of_air(
        "hall",
        Length::m(4.4),
        Length::m(3.1),
        Length::m(2.4),
        across,
    )
}

/// Run a hall at a fraction of its Courant limit, returning the peak at each step.
fn run(h: &mut Hall, dt: Time, steps: usize) -> Vec<(f64, f64)> {
    let mut bus = Exchange::new();
    let mut t = 0.0;
    let mut out = Vec::with_capacity(steps);
    for _ in 0..steps {
        h.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
        out.push((t, h.peak_pressure().to_si()));
    }
    out
}

/// **The mode frequency is the closed form, and the third axis is in it.**
///
/// `(c/2)·√((a/Lx)² + (b/Ly)² + (c/Lz)²)`, computed here rather than taken from the domain.
/// The `(0,0,1)` mode is the one a two-dimensional model does not have at all — not less
/// accurately, at all — so it is checked explicitly.
#[test]
fn the_mode_frequencies_are_the_closed_form() {
    let h = hall(45);
    let (lx, ly, lz) = (h.width().to_si(), h.height().to_si(), h.depth().to_si());
    let want = |a: f64, b: f64, c: f64| {
        343.0 / 2.0 * ((a / lx).powi(2) + (b / ly).powi(2) + (c / lz).powi(2)).sqrt()
    };

    for (mode, a, b, c) in [
        ((1, 0, 0), 1.0, 0.0, 0.0),
        ((0, 1, 0), 0.0, 1.0, 0.0),
        ((0, 0, 1), 0.0, 0.0, 1.0),
        ((1, 1, 1), 1.0, 1.0, 1.0),
        ((3, 2, 1), 3.0, 2.0, 1.0),
    ] {
        let got = h.mode_frequency(mode).to_si();
        assert!(
            (got - want(a, b, c)).abs() < 1e-9,
            "{mode:?}: {got:.6} Hz against {:.6}",
            want(a, b, c)
        );
    }

    // The vertical mode is the one a floor plan cannot have. 2.4 m of ceiling puts it at 71 Hz,
    // well inside the range a room is judged on.
    let vertical = h.mode_frequency((0, 0, 1)).to_si();
    assert!(
        (vertical - 343.0 / (2.0 * lz)).abs() < 1e-9,
        "the vertical mode is c/2Lz: {vertical:.4} Hz"
    );
    assert!((60.0..85.0).contains(&vertical));
}

/// **A released mode rides `|cos(2πft)|`, and the departure converges at second order.**
///
/// A standing mode is separable, so every point follows the same cosine and the peak of the field
/// is `|cos(2πft)|` exactly. That is the closed form; the integration never sees it.
///
/// The **rate** is what is asserted. A first-order scheme also tracks a cosine and its error also
/// gets small — what separates them is that halving the spacing quarters one error and halves the
/// other. `Room` and `Tube` both shipped first order while looking merely coarse, and this domain
/// exists partly to demonstrate that the fix was understood rather than pattern-matched.
///
/// **Measured across three doublings, not one**, and that is not a convenience. "Worst departure
/// over a run" is a maximum and therefore noisy; the per-doubling ratios measured here are
///
/// ```text
///     7 nodes   2.320e-2
///    13 nodes   4.123e-3    5.63x
///    25 nodes   1.402e-3    2.94x
///    49 nodes   3.020e-4    4.64x
///    97 nodes   7.549e-5    4.00x
/// ```
///
/// — bouncing around four rather than sitting on it, so a single doubling proves nothing in
/// either direction. Over 13 to 97 the fall is **54.6x**, against 64x for second order and 8x for
/// first, which separates them with room on both sides. `Room`'s own convergence test reaches the
/// same conclusion by the same route, and for the same reason.
#[test]
fn a_released_mode_tracks_the_cosine_and_converges_at_second_order() {
    let worst_at = |across: usize| {
        let mut h = hall(across).released_in_mode((1, 1, 1), Pressure::from_si(1.0));
        let f = h.mode_frequency((1, 1, 1)).to_si();
        // A fixed fraction of the limit, so the step refines with the grid — which is what makes
        // this a test of the *scheme* rather than of the step size.
        let dt = Time::from_si(h.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
        let steps = (0.006 / dt.to_si()).round() as usize;
        run(&mut h, dt, steps)
            .into_iter()
            .map(|(t, peak)| (peak - (2.0 * std::f64::consts::PI * f * t).cos().abs()).abs())
            .fold(0.0f64, f64::max)
    };

    let (coarse, fine) = (worst_at(13), worst_at(97));
    let fall = coarse / fine;
    assert!(
        fall > 30.0,
        "13 -> 97 nodes is three doublings: second order is 64x, first order 8x. \
         {coarse:.6} -> {fine:.6}, a factor of {fall:.1}"
    );
    // And the size, so the rate is not a ratio of two numbers that are both wrong: 97 nodes
    // across tracks a 1 Pa mode to under a thousandth of its amplitude.
    assert!(fine < 1e-3, "97 nodes departed by {fine:.6} Pa");
}

/// **Rigid surfaces conserve energy exactly, and the startup error is `O(h²)`.**
///
/// Exact rather than to a tolerance: no velocity passes a rigid wall, so the discrete invariant
/// is conserved to rounding from the first step onward. The one departure is the leapfrog
/// startup, which converts a state with velocity at `t = 0` into one with velocity at `t = −h/2`
/// — that is discretisation, and `startup_adjustment` reports it rather than hiding it.
#[test]
fn a_rigid_hall_conserves_its_invariant_and_reports_its_startup() {
    let adjustment_at = |across: usize| {
        // (1,1,1), so every axis carries a velocity. This said (1,1,0) first, and a mutation
        // showed why that was worthless: with no z variation `vz` stays zero, so the z wall
        // weight contributes nothing to either the update or the invariant and a wrong one
        // cancels itself out of both. The mode has to excite what the test claims to check.
        let mut h = hall(across).released_in_mode((1, 1, 1), Pressure::from_si(2.0));
        let dt = Time::from_si(h.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
        run(&mut h, dt, 200);
        let energies: Vec<f64> = {
            let mut e = Vec::new();
            let mut bus = Exchange::new();
            for _ in 0..200 {
                h.step(Time::from_si(0.0), dt, &mut bus).expect("stable");
                e.push(h.energy().to_si());
            }
            e
        };
        let first = energies[0];
        let drift = energies
            .iter()
            .map(|e| (e - first).abs())
            .fold(0.0f64, f64::max);
        (
            drift / first.abs(),
            h.startup_adjustment().to_si().abs() / first.abs(),
        )
    };

    let (drift, coarse_start) = adjustment_at(13);
    assert!(
        drift < 1e-12,
        "after the startup the invariant is exact, drifted {drift:.3e}"
    );

    // The startup adjustment is real and quarters on refinement — O(h²), because h follows dx.
    let (_, fine_start) = adjustment_at(25);
    assert!(
        coarse_start > 1e-6,
        "the startup adjustment should be a real number, got {coarse_start:.3e}"
    );
    let fall = coarse_start / fine_start;
    assert!(
        (2.5..6.0).contains(&fall),
        "the startup error is O(h²): {coarse_start:.3e} -> {fine_start:.3e}, {fall:.2}x"
    );
}

/// **The Courant limit is `dx/(c√3)`, and past it the scheme is refused.**
///
/// `√3` rather than `√2`, because a wave crossing a cube diagonally covers `√3` cells while the
/// stencil sees one. Checked against arithmetic done here, and against the two-dimensional limit
/// it must be tighter than.
#[test]
fn the_courant_limit_is_the_diagonal_of_a_cube() {
    let h = hall(45);
    let dx = h.spacing().to_si();
    let limit = h.max_stable_dt(Time::from_si(0.0)).to_si();
    let exact = dx / (343.0 * 3f64.sqrt());
    assert!(
        (limit / exact - 1.0).abs() < 1e-12,
        "dx/(c√3): {limit:.9e} against {exact:.9e}"
    );
    assert!((h.courant(Time::from_si(limit)) - 1.0).abs() < 1e-12);

    // Tighter than the two-dimensional limit by exactly √(3/2) — 22%, which is what the crate
    // docs have claimed since before there was a three-dimensional domain to check it against.
    let two_d = dx / (343.0 * 2f64.sqrt());
    assert!(
        (two_d / limit - (1.5f64).sqrt()).abs() < 1e-12,
        "3D against 2D should be √1.5 = 1.2247, got {:.6}",
        two_d / limit
    );

    let mut h = hall(45).released_in_mode((1, 0, 0), Pressure::from_si(1.0));
    let err = h
        .step(
            Time::from_si(0.0),
            Time::from_si(limit * 1.05),
            &mut Exchange::new(),
        )
        .expect_err("5% past the limit must be refused");
    assert_eq!(err.quantity, "Courant number");
    assert!(
        (err.after - 1.05).abs() < 1e-9,
        "by how much: {}",
        err.after
    );

    // And exactly at it is accepted, so this is a limit rather than a blanket refusal.
    hall(45)
        .released_in_mode((1, 0, 0), Pressure::from_si(1.0))
        .step(
            Time::from_si(0.0),
            Time::from_si(limit),
            &mut Exchange::new(),
        )
        .expect("exactly at the limit is stable");
}

/// **A hall one node thick is a room, and one node thick twice over is a tube.**
///
/// The reduction that ties the three dimensionalities together. An axis of one node has no faces
/// along it, contributes nothing to the divergence, and must leave the remaining axes behaving as
/// the lower-dimensional scheme does — which is checkable through the mode frequency, since a
/// collapsed axis drops out of the closed form as well.
#[test]
fn collapsing_an_axis_reduces_to_the_lower_dimension() {
    // One node deep: 0 m of depth, so the z term vanishes and the frequencies are the 2D ones.
    let flat = Hall::of_air("flat", Length::m(4.4), Length::m(3.1), Length::m(0.0), 45);
    assert_eq!(flat.nodes().2, 1, "zero depth is one node");
    let (lx, ly) = (flat.width().to_si(), flat.height().to_si());
    let want = 343.0 / 2.0 * ((1.0 / lx).powi(2) + (1.0 / ly).powi(2)).sqrt();
    assert!(
        (flat.mode_frequency((1, 1, 0)).to_si() - want).abs() < 1e-9,
        "a collapsed axis should drop out of the closed form"
    );

    // And it still runs, conserves and oscillates, rather than dividing by a zero extent.
    let mut flat = flat.released_in_mode((1, 1, 0), Pressure::from_si(1.0));
    let dt = Time::from_si(flat.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
    let trace = run(&mut flat, dt, 400);
    let lowest = trace.iter().map(|(_, p)| *p).fold(f64::MAX, f64::min);
    assert!(
        lowest < 0.2,
        "the mode should pass through zero, lowest peak was {lowest:.4}"
    );
    let highest = trace.iter().map(|(_, p)| *p).fold(0.0f64, f64::max);
    assert!(
        highest <= 1.0 + 1e-9,
        "and never exceed its release amplitude: {highest:.6}"
    );
}

/// **The mode count grows as `f³`, which is the qualitative thing 2D gets wrong.**
///
/// A two-dimensional room's modes fill an *area* of the index lattice and grow as `f²`; a real
/// room's fill a volume. That is why resonances merge into a hiss above the Schroeder frequency,
/// and it is not a quantitative correction — it is a different curve.
///
/// Checked against Weyl's estimate, **with its surface and edge terms**:
///
/// ```text
///   N(f) = (4π/3)·V·(f/c)³ + (π/4)·S·(f/c)² + (1/8)·L·(f/c)
/// ```
///
/// The leading term alone is not close enough to assert against — it gives 92 modes below 300 Hz
/// where this room has 134, a 46% gap — because for **rigid** walls the surface term is positive
/// and large at these frequencies. That is a real feature of the Neumann problem rather than a
/// fudge: the Dirichlet version of the same formula subtracts it. Room-acoustics texts quote the
/// three-term version for exactly this reason.
#[test]
fn the_mode_count_grows_as_the_cube_of_frequency() {
    use pantometry_core::units::Frequency;
    let h = hall(45);
    let below = |f: f64| h.modes_below(Frequency::from_si(f)).len() as f64;

    let (low, high) = (below(150.0), below(300.0));
    assert!(low > 10.0, "150 Hz should already hold modes, got {low}");
    let ratio = high / low;
    assert!(
        (5.5..10.0).contains(&ratio),
        "doubling the frequency should roughly octuple the count, got {ratio:.2}"
    );

    // Against the three-term Weyl estimate.
    let (lx, ly, lz) = (h.width().to_si(), h.height().to_si(), h.depth().to_si());
    let v = lx * ly * lz;
    let surface = 2.0 * (lx * ly + ly * lz + lz * lx);
    let edges = 4.0 * (lx + ly + lz);
    let pi = std::f64::consts::PI;
    let weyl = |f: f64| {
        let k = f / 343.0;
        4.0 * pi / 3.0 * v * k.powi(3) + pi / 4.0 * surface * k.powi(2) + edges * k / 8.0
    };
    let predicted = weyl(300.0);
    assert!(
        (0.9..1.1).contains(&(high / predicted)),
        "{high} modes below 300 Hz against Weyl's {predicted:.1}"
    );
    // The leading term alone is not enough, and saying so is the point: it is 46% low here.
    let leading = 4.0 * pi / 3.0 * v * (300.0f64 / 343.0).powi(3);
    assert!(
        high / leading > 1.2,
        "the surface term should matter at this frequency: {high} against {leading:.1}"
    );

    // Ascending, and the first one is the longest axis.
    let modes = h.modes_below(Frequency::from_si(120.0));
    for pair in modes.windows(2) {
        assert!(pair[0].1.to_si() <= pair[1].1.to_si(), "not sorted");
    }
    assert_eq!(modes[0].0, (1, 0, 0), "the lowest mode is along the width");
}
