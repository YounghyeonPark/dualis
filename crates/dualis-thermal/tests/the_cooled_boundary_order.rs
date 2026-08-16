//! How fast a cooled face converges, measured — because a rate is what tells a wrong boundary
//! from a coarse one.
//!
//! This workspace has twice found a wrong boundary condition hiding inside a plausible error,
//! and both times the thing that exposed it was the **order**: refining the grid halved the
//! error instead of quartering it, which a coarse condition cannot do and a wrong one must.
//! When `Solid3D` gained a cooled face, a review measured exactly that signature — first order
//! at the boundary against second in the interior — because the film was charged against the
//! cell centre with the half cell of solid between centre and surface left out of the path.
//! `1/(1/G_film + dx/2kA)` is the fix, and this file is the measurement that says it worked.
//!
//! # Why the measurement is a power iteration and not one step
//!
//! The first attempt at this test sampled the **continuum** mode `cos(λz/L)` at the cell
//! centres and read the decay off a single step. That is exact only when the state is the
//! *discrete* operator's eigenvector, which it is not: the sampled continuum mode carries a
//! little of every faster mode, and the one-step reading is contaminated by all of them. It
//! measured a bias that did not shrink under refinement and was not the boundary's.
//!
//! Marching until the faster modes have died and then reading the ratio between successive
//! steps is the power method, and it converges on the discrete operator's own slowest
//! eigenvalue — which is the number that must approach the continuum one at second order.

use dualis_core::units::{Area, Length, Temperature, Time};
use dualis_core::{Domain, Exchange, Substance};
use dualis_thermal::{Environment, Face, Solid3D};

/// The slab: `cells` deep, insulated everywhere but `ZMax`, which sees a film of `h`.
///
/// Non-radiating on purpose. Radiation would make the problem nonlinear and the closed form
/// below is the **linear** Robin condition, so `ε = 0` is what keeps this a comparison against
/// a closed form rather than against a nearby one.
fn slab(cells: usize, h: f64, length: f64) -> Solid3D {
    let mut glass = Substance::from_name("borosilicate").expect("the catalogue has it");
    if let Some(thermal) = glass.thermal.as_mut() {
        thermal.emissivity = 0.0;
    }
    let dx = length / cells as f64;
    Solid3D::new(
        "slab",
        glass,
        (1, 1, cells),
        Length::from_si(dx),
        Temperature::celsius(120.0),
    )
    .losing_from(
        Face::ZMax,
        Environment {
            ambient: Temperature::celsius(20.0),
            convection_w_per_m2_k: h,
            area: Area::from_si(dx * dx),
        },
    )
}

/// The slowest decay rate the marched operator actually has, in s⁻¹.
///
/// A power iteration: march until the faster modes are gone, then read the ratio between two
/// successive states. The step is a fixed fraction of the limit so the *time* discretisation is
/// the same at every grid and cannot be mistaken for the boundary's contribution — what is
/// being compared across grids is the space.
fn measured_rate(cells: usize, h: f64, length: f64) -> f64 {
    let mut block = slab(cells, h, length);
    let ambient = Temperature::celsius(20.0).to_si();
    let dt = block.max_stable_dt(Time::ZERO).to_si() * 0.2;
    let mut bus = Exchange::new();

    let amplitude = |b: &Solid3D| {
        (0..cells)
            .map(|k| b.temperature_at(0, 0, k).to_si() - ambient)
            .sum::<f64>()
    };

    // Long enough for the second mode to be gone: it decays about nine times faster than the
    // first for a slab of this Biot number, so a few first-mode time constants leave it at
    // machine level. Scaled with the grid because the step is scaled with the grid.
    let settle = 40 * cells * cells;
    for _ in 0..settle {
        block
            .step(Time::ZERO, Time::from_si(dt), &mut bus)
            .expect("a stable step");
    }

    let before = amplitude(&block);
    block
        .step(Time::ZERO, Time::from_si(dt), &mut bus)
        .expect("a stable step");
    let after = amplitude(&block);
    (1.0 - after / before) / dt
}

/// The continuum answer: the slowest mode of a slab insulated on one face and cooled on the
/// other decays at `α λ²/L²`, where `λ tan λ = Bi` and `Bi = hL/k`.
fn exact_rate(h: f64, length: f64) -> f64 {
    let glass = Substance::from_name("borosilicate").expect("the catalogue has it");
    let thermal = glass.thermal.expect("borosilicate conducts");
    let k = thermal.conductivity.to_si();
    let alpha = k / glass.density.to_si() / thermal.specific_heat.to_si();
    let bi = h * length / k;

    let (mut lo, mut hi) = (1e-12, std::f64::consts::FRAC_PI_2 - 1e-12);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid * mid.tan() < bi {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let lambda = 0.5 * (lo + hi);
    alpha * lambda * lambda / (length * length)
}

/// **The cooled boundary converges at second order.**
///
/// The error against `α λ²/L²` has to fall by four per grid doubling. First order halves it,
/// and that is the signature the review measured before the half cell went into the path —
/// 2.08, 2.04, 2.02, against the interior's 4.0.
///
/// The band is `> 3.4` rather than `≈ 4`, and it is earned rather than chosen: the next term in
/// the expansion is `O(dx³)`, so at the coarse end the measured ratio approaches four **from
/// below** by a margin that shrinks with the grid. What it must never do is sit near two, and
/// 3.4 is comfortably clear of it while leaving the coarse pair its contamination.
#[test]
fn the_cooled_boundary_converges_at_second_order() {
    let length = 0.048;
    let h = 40.0;
    let exact = exact_rate(h, length);

    let errors: Vec<f64> = [4usize, 8, 16, 32]
        .iter()
        .map(|n| (measured_rate(*n, h, length) - exact).abs() / exact)
        .collect();
    let ratios: Vec<f64> = errors.windows(2).map(|w| w[0] / w[1]).collect();

    for (n, r) in ratios.iter().enumerate() {
        assert!(
            *r > 3.4,
            "doubling {n} improved the boundary by {r:.3}, which is first order dressed as \
             second — the errors were {errors:?}, ratios {ratios:?}"
        );
    }
    assert!(
        errors[3] < 5e-3,
        "the finest grid should be close, got {:.4e} — errors {errors:?}",
        errors[3]
    );
}

/// **And it converges at second order for a stiff film too**, where the boundary carries most of
/// the resistance and a wrong condition has the most room to hide. `h = 400` puts the Biot
/// number at 17, an order above the lumped criterion.
#[test]
fn a_stiff_film_converges_at_second_order_as_well() {
    let length = 0.048;
    let h = 400.0;
    let exact = exact_rate(h, length);

    let errors: Vec<f64> = [4usize, 8, 16, 32]
        .iter()
        .map(|n| (measured_rate(*n, h, length) - exact).abs() / exact)
        .collect();
    let ratios: Vec<f64> = errors.windows(2).map(|w| w[0] / w[1]).collect();

    for (n, r) in ratios.iter().enumerate() {
        assert!(
            *r > 3.4,
            "doubling {n} improved the boundary by {r:.3} at Bi = 17 — errors {errors:?}"
        );
    }
}
