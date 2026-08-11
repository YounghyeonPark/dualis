//! What a field does inside something lossy, and why a lumped resistance cannot say it.
//!
//! `dualis-electrical`'s `Winding` states a resistance and dissipates `I²R`. That is exactly right
//! for a wire at low frequency and it is the whole of what a lumped model can say: `R` is an input.
//! `Conductor` improved on it by solving `∇·(σ∇φ) = 0`, so a resistance became a property of a
//! *shape* — but it is quasi-static, and a quasi-static solve has no frequency in it at all.
//!
//! A field does. Current crowds toward the surface as the frequency rises, and the depth it crowds
//! into has an exact closed form for **any** conductivity:
//!
//! ```text
//!   p = σ/(ωε)                            the loss tangent
//!   α = ω√(εμ/2) · √(√(1+p²) − 1)         the amplitude decay, exactly
//!   α → √(ωμσ/2) = 1/δ                    the skin depth, as p → ∞
//! ```
//!
//! The second line is not an approximation. The third is the one every handbook quotes, and it is
//! the `p ≫ 1` limit of the second — so a test that checks both checks the formula and its own
//! famous special case.
//!
//! # One run, three frequencies
//!
//! A pulse carries a band, so a single march gives the decay at every frequency in it. Taking the
//! discrete Fourier transform of the recorded field at each depth separates them, and each comes
//! with its own `α` — three independent comparisons from one run, at three different loss tangents.
//!
//! That is also the honest way to measure a frequency-dependent quantity: driving a single
//! sinusoid and waiting for steady state measures one number and hides whether the *dependence* is
//! right.

use dualis::em::{Boundary, Cavity, Medium, Wall};
use dualis::prelude::*;

/// Cells per vacuum wavelength.
const PER_WAVELENGTH: usize = 24;

/// The exact amplitude decay constant of a conducting medium, in nepers per metre.
fn decay_constant(omega: f64, permittivity: f64, permeability: f64, conductivity: f64) -> f64 {
    let p = conductivity / (omega * permittivity);
    omega * (permittivity * permeability / 2.0).sqrt() * ((1.0 + p * p).sqrt() - 1.0).sqrt()
}

/// The handbook skin depth, `√(2/(ωμσ))`.
fn skin_depth(omega: f64, permeability: f64, conductivity: f64) -> f64 {
    (2.0 / (omega * permeability * conductivity)).sqrt()
}

/// The magnitude of one frequency in a time series, by direct transform.
fn amplitude_at(series: &[f64], frequency: f64, dt: f64) -> f64 {
    let (mut re, mut im) = (0.0, 0.0);
    for (n, v) in series.iter().enumerate() {
        let phase = -2.0 * std::f64::consts::PI * frequency * n as f64 * dt;
        re += v * phase.cos();
        im += v * phase.sin();
    }
    (re * re + im * im).sqrt()
}

/// Fit `ln A = c − αz` by least squares, returning `α`.
fn fitted_decay(depths: &[f64], amplitudes: &[f64]) -> f64 {
    let n = depths.len() as f64;
    let (sx, sy): (f64, f64) = (
        depths.iter().sum(),
        amplitudes.iter().map(|a| a.ln()).sum::<f64>(),
    );
    let sxy: f64 = depths.iter().zip(amplitudes).map(|(z, a)| z * a.ln()).sum();
    let sxx: f64 = depths.iter().map(|z| z * z).sum();
    -((n * sxy - sx * sy) / (n * sxx - sx * sx))
}

/// **The decay inside a lossy medium is the exact `α(ω)`, at three frequencies from one run.**
///
/// A pulse enters a conducting half-space and its amplitude falls exponentially. The rate depends on
/// the frequency, so measuring it at three frequencies in one pulse's band checks the *dependence*
/// and not just a number — which is the thing a lumped resistance has none of.
#[test]
fn a_field_decays_in_a_conductor_at_the_exact_rate() {
    let wavelength = 1.0;
    let dx = wavelength / PER_WAVELENGTH as f64;
    let f0 = 299_792_458.0 / wavelength;
    // A loss tangent near a third at the carrier: strong enough to decay inside the box, gentle
    // enough that a decay length is twenty cells rather than two.
    let epsilon = dualis::em::EPSILON_0;
    let mu = dualis::em::MU_0;
    let conductivity = 0.3 * 2.0 * std::f64::consts::PI * f0 * epsilon;

    let interface = 8 * PER_WAVELENGTH;
    let nz = 20 * PER_WAVELENGTH;
    let mut c = Cavity::new("slab", (2, 1, nz), Length::from_si(dx), Medium::vacuum());
    c.set_boundary(Wall::XLow, Boundary::Magnetic);
    c.set_boundary(Wall::XHigh, Boundary::Magnetic);
    c.set_boundary(Wall::ZLow, Boundary::Open);
    c.set_boundary(Wall::ZHigh, Boundary::Open);
    c.fill(
        Medium {
            conductivity,
            ..Medium::vacuum()
        },
        |_, _, k| k >= interface,
    );

    let dt = Time::from_si(c.courant_limit().to_si() * 0.5);
    let k0 = 2.0 * std::f64::consts::PI / wavelength;
    let (centre, width) = (3.0 * wavelength, 2.0 * wavelength);
    c.launch_along_z(dt, move |z| {
        let u = (z - centre) / width;
        (-u * u).exp() * (k0 * (z - centre)).cos()
    });

    // Depths well inside the medium, and clear of the interface so the measurement is of the
    // travelling wave rather than of the reflection standing on top of it.
    let probes: Vec<usize> = (1..=6)
        .map(|n| interface + n * PER_WAVELENGTH / 2)
        .collect();
    let mut series: Vec<Vec<f64>> = vec![Vec::new(); probes.len()];
    let mut bus = Exchange::new();
    // Long enough for the pulse to reach the deepest probe and pass it completely.
    let steps = (1.6 * nz as f64 / 0.2887) as usize;
    for n in 0..steps {
        c.step(Time::from_si(n as f64 * dt.to_si()), dt, &mut bus)
            .expect("stable");
        for (s, p) in series.iter_mut().zip(&probes) {
            s.push(c.electric_at(1, 0, *p).y);
        }
    }

    let depths: Vec<f64> = probes
        .iter()
        .map(|p| (*p - interface) as f64 * dx)
        .collect();
    for factor in [0.85, 1.0, 1.15] {
        let f = f0 * factor;
        let omega = 2.0 * std::f64::consts::PI * f;
        let amplitudes: Vec<f64> = series
            .iter()
            .map(|s| amplitude_at(s, f, dt.to_si()))
            .collect();
        let measured = fitted_decay(&depths, &amplitudes);
        let closed = decay_constant(omega, epsilon, mu, conductivity);
        println!(
            "  f/f0 = {factor:.2}  p = {:.4}:  alpha {:.5} against {:.5} /m — off {:.2}%",
            conductivity / (omega * epsilon),
            measured,
            closed,
            (measured / closed - 1.0).abs() * 100.0
        );
        assert!(
            (measured / closed - 1.0).abs() < 0.08,
            "the decay is the exact alpha: {measured:.5} against {closed:.5}"
        );
    }
}

/// **The famous skin depth is the `p ≫ 1` limit of that formula, and the approach is measured.**
///
/// `δ = √(2/(ωμσ))` is what every handbook quotes and it is an approximation — good when the loss
/// tangent is large and badly wrong when it is not. Measured against the exact `α`:
///
/// ```text
///   p =    1     off by 55.4%
///   p =   10            5.1%
///   p =  100            0.50%
///   p = 1000            0.050%
/// ```
///
/// A tenfold in the loss tangent is a tenfold in the accuracy, which is what a first-order limit
/// looks like. The 55% at `p = 1` is the part worth knowing: a lossy dielectric is not a conductor,
/// and using the skin depth on one is not a small error.
///
/// So this is a limit rather than a value, which is the same shape as `Conductor`'s Maxwell
/// constriction: the closed form is approached and the *approach* is the claim.
#[test]
fn the_skin_depth_is_the_good_conductor_limit() {
    let mu = dualis::em::MU_0;
    let epsilon = dualis::em::EPSILON_0;
    let omega = 2.0 * std::f64::consts::PI * 1e9;
    let mut errors = Vec::new();
    for p in [1.0, 10.0, 100.0, 1000.0] {
        let sigma = p * omega * epsilon;
        let exact = decay_constant(omega, epsilon, mu, sigma);
        let handbook = 1.0 / skin_depth(omega, mu, sigma);
        let error = (handbook / exact - 1.0).abs();
        println!(
            "  p = {p:>6}: 1/delta {handbook:.5} against the exact {exact:.5} — off {:.3}%",
            error * 100.0
        );
        errors.push(error);
    }
    for pair in errors.windows(2) {
        assert!(
            pair[1] < pair[0],
            "the approximation must improve with the loss tangent: {:.4} then {:.4}",
            pair[0],
            pair[1]
        );
    }
    assert!(
        errors[0] > 0.05,
        "and at p = 1 it must be visibly wrong, or this is not a limit: off by {:.2}%",
        errors[0] * 100.0
    );
    assert!(
        errors[3] < 1e-3,
        "while at p = 1000 it is the answer: off by {:.4}%",
        errors[3] * 100.0
    );
}

/// **A lumped resistance has no frequency in it, and that is the gap this measures.**
///
/// `Winding::resistance` is a number. `Conductor` solves for a shape but is quasi-static, so its
/// answer is the same at every frequency too. A field is not: the decay constant above rises as
/// `√ω` once the loss tangent is large, so the depth carrying the current shrinks and the effective
/// resistance of the same piece of metal grows.
///
/// This states the ratio rather than computing a resistance, because a resistance needs a geometry
/// and the point is the *frequency* dependence. A hundredfold in frequency is a tenfold in the
/// crowding, and neither lumped model can produce either.
#[test]
fn the_depth_carrying_the_current_shrinks_as_the_root_of_the_frequency() {
    let mu = dualis::em::MU_0;
    let sigma = 5.96e7; // copper
    let mut depths = Vec::new();
    for f in [50.0, 5_000.0, 500_000.0] {
        let d = skin_depth(2.0 * std::f64::consts::PI * f, mu, sigma);
        println!("  copper at {f:>9} Hz: delta {:.4} mm", d * 1e3);
        depths.push(d);
    }
    // A hundredfold in frequency is a tenfold in the depth, exactly.
    for pair in depths.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (ratio - 10.0).abs() < 1e-9,
            "delta goes as 1/sqrt(f): {ratio:.6} for a hundredfold"
        );
    }
    // And at mains frequency it is millimetres, which is why a busbar is a busbar and not a rod.
    assert!(
        (8.0..10.0).contains(&(depths[0] * 1e3)),
        "copper at 50 Hz is about 9 mm: {:.3} mm",
        depths[0] * 1e3
    );
    // The lumped model's own answer, for contrast: `Winding` states a resistance and it is the same
    // number at every one of those frequencies.
    let w = dualis::electrical::Winding::of_copper(
        "coil",
        Length::m(10.0),
        1e-6,
        Temperature::celsius(20.0),
    );
    println!(
        "  and Winding says {:.4} ohm, at every frequency above",
        w.resistance().to_si()
    );
    assert!(
        w.resistance().to_si() > 0.0,
        "the lumped model has an answer; what it does not have is a frequency"
    );
}
