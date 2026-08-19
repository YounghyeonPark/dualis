//! A rectangular waveguide, and the whole dispersion relation from one march.
//!
//! A box with conducting side walls and open ends **is** a waveguide, so this needs nothing the
//! cavity tests did not already have. What it measures is the thing a cavity cannot show: a mode
//! that travels, and the frequency below which it stops.
//!
//! ```text
//!   f_c = (c/2)√((m/a)² + (n/b)²)      the cutoff, from the cross-section alone
//!   β  = √(k² − k_c²)                  above it: propagating, and slower to turn its phase
//!   α  = √(k_c² − k²)                  below it: evanescent, decaying without absorbing
//!   λ_g = λ / √(1 − (f_c/f)²)          so the guide wavelength is always longer than free space
//!   v_p v_g = c²                       the phase outruns light and the pulse does not
//! ```
//!
//! # Why one run gives all of it
//!
//! A pulse carries a band. The mode's profile is `sin(πx/a)`, which is orthogonal to every other
//! mode's, so what propagates downstream is TE₁₀ and nothing else — at every frequency in the band
//! at once. A discrete transform of the recorded field at a row of stations gives, for each
//! frequency, an amplitude and a phase against distance: the phase slope is `β` above cutoff and the
//! amplitude slope is `α` below it.
//!
//! So the check is a **curve** rather than a point, and it spans the cutoff. That matters because
//! `β = √(k² − k_c²)` and `β ≈ k` agree to a percent at high frequency and disagree by everything
//! near the cutoff: a test at one frequency well above it would pass for a solver with no cutoff at
//! all.
//!
//! # This is why a fibre has modes
//!
//! Not the same geometry — a fibre guides by index contrast rather than by conducting walls — but
//! the same statement: a transverse dimension quantises the wavenumber, what is left over goes into
//! propagation, and below a threshold there is nothing left over.

use pantometry_core::{Domain, Exchange};
use pantometry_em::{Boundary, Cavity, Medium, Wall};
use pantometry_units::{Length, Time};

/// Cells across the guide's wide dimension. The cutoff wavelength is twice this.
const ACROSS: usize = 20;
/// Cells along it.
const ALONG: usize = 420;

/// One station's recorded field, transformed at a frequency: `(amplitude, phase)`.
fn transform(series: &[f64], frequency: f64, dt: f64) -> (f64, f64) {
    let (mut re, mut im) = (0.0, 0.0);
    for (n, v) in series.iter().enumerate() {
        let phase = -2.0 * std::f64::consts::PI * frequency * n as f64 * dt;
        re += v * phase.cos();
        im += v * phase.sin();
    }
    ((re * re + im * im).sqrt(), im.atan2(re))
}

/// Fit a slope by least squares.
fn slope(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (sx, sy): (f64, f64) = (x.iter().sum(), y.iter().sum());
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| a * b).sum();
    let sxx: f64 = x.iter().map(|a| a * a).sum();
    (n * sxy - sx * sy) / (n * sxx - sx * sx)
}

/// Make a monotone phase out of one that wraps at ±π.
fn unwrapped(phases: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(phases.len());
    let mut offset = 0.0;
    for (n, p) in phases.iter().enumerate() {
        if n > 0 {
            let jump = p + offset - out[n - 1];
            if jump > std::f64::consts::PI {
                offset -= 2.0 * std::f64::consts::PI;
            } else if jump < -std::f64::consts::PI {
                offset += 2.0 * std::f64::consts::PI;
            }
        }
        out.push(p + offset);
    }
    out
}

/// What one march produced: where the stations are, and what they recorded.
struct Run {
    dt: f64,
    cutoff: f64,
    stations: Vec<f64>,
    series: Vec<Vec<f64>>,
}

impl Run {
    /// The propagation constant at a frequency, from the phase slope along the guide.
    fn beta(&self, frequency: f64) -> f64 {
        let phases: Vec<f64> = self
            .series
            .iter()
            .map(|s| transform(s, frequency, self.dt).1)
            .collect();
        -slope(&self.stations, &unwrapped(&phases))
    }
}

/// Send a TE₁₀ pulse down the guide and record it at a row of stations.
///
/// # The pulse's band and the stations' spacing are both measurements, not decorations
///
/// Two mistakes, both made first and both silent:
///
/// A **band that does not cover the frequency being asked about** returns a transform of rounding.
/// A pulse of half-width 2.5 wavelengths has a fractional bandwidth near 0.13, so at 0.6 of its
/// carrier there is nothing there — and the first version of these tests read the propagation
/// constant correctly at the carrier, to 0.84%, and 87% wrong two octaves away. Each test now
/// marches with a carrier near the frequencies it asks about.
///
/// A **spacing wider than what is being measured** aliases or vanishes. Above cutoff the stations
/// must be closer than half a guide wavelength or the unwrapped phase is wrong; below it they must
/// span a couple of decay lengths and no more, because a field that has fallen by `e⁻²⁵` is at the
/// numerical floor and its slope is noise. Those are different spacings, so each test chooses.
fn march(
    carrier_over_cutoff: f64,
    width_wavelengths: f64,
    first_cell: usize,
    spacing: usize,
    count: usize,
) -> Run {
    // The cutoff wavelength is `2a`, so pick the cell from the guide's width.
    let dx = 1.0 / ACROSS as f64;
    let a = ACROSS as f64 * dx;
    let mut c = Cavity::new(
        "guide",
        (ACROSS, 4, ALONG),
        Length::from_si(dx),
        Medium::vacuum(),
    );
    // The four side walls are conductors, which is what a guide is; the ends are open, which is
    // what makes it a guide rather than a cavity.
    c.set_boundary(Wall::ZLow, Boundary::Open);
    c.set_boundary(Wall::ZHigh, Boundary::Open);

    let cutoff = c.cutoff_frequency((1, 0)).to_si();
    let carrier = cutoff * carrier_over_cutoff;
    let wavelength = 299_792_458.0 / carrier;
    let dt = Time::from_si(c.courant_limit().to_si() * 0.5);

    let k = 2.0 * std::f64::consts::PI / wavelength;
    let (centre, width) = (6.0 * a, width_wavelengths * wavelength);
    let pi = std::f64::consts::PI;
    c.launch_along_z(dt, move |x, _y, z| {
        let u = (z - centre) / width;
        (pi * x / a).sin() * (-u * u).exp() * (k * (z - centre)).cos()
    });

    let probes: Vec<usize> = (0..count).map(|n| first_cell + n * spacing).collect();
    let mut series: Vec<Vec<f64>> = vec![Vec::new(); probes.len()];
    let mut bus = Exchange::new();
    // **Stopped before either end's reflection returns.** Mur's coefficient is tuned to `c`, and a
    // guided mode's phase velocity is `c/√(1−(f_c/f)²)` — 1.67c just above cutoff. So the absorbing
    // faces are badly mismatched exactly where this test is most interesting, and they reflect.
    //
    // Letting the march run long enough for the far end's reflection to come back put a backward wave
    // under the stations and read a guide wavelength *shorter* than free space, which is impossible.
    // Excluding that left the near end's, which travels source → `ZLow` → stations and arrives at
    // about three and a half times the direct path: a window that ends between the two is the whole
    // trick, and it is why this is 0.8 of a length rather than two.
    let steps = (0.8 * ALONG as f64 / 0.2887) as usize;
    for n in 0..steps {
        c.step(Time::from_si(n as f64 * dt.to_si()), dt, &mut bus)
            .expect("stable");
        for (s, p) in series.iter_mut().zip(&probes) {
            s.push(c.electric_at(ACROSS / 2, 2, *p).y);
        }
    }
    Run {
        dt: dt.to_si(),
        cutoff,
        stations: probes.iter().map(|p| *p as f64 * dx).collect(),
        series,
    }
}

/// **Above cutoff the guide wavelength is `λ/√(1−(f_c/f)²)`, over two and a half octaves.**
///
/// Measured from the phase slope, so it is the propagation constant itself rather than a wavelength
/// counted off a picture. Two marches with different carriers, each asked about frequencies inside
/// its own band, so the range spans from `1.4 f_c` — where the guide wavelength is 43% longer than
/// the free-space one — to `3.1 f_c`, where the two agree to 5%.
///
/// That range is the test. `β = √(k² − k_c²)` and `β ≈ k` agree to a percent at the top of it and
/// disagree by everything at the bottom, so a check at one high frequency would pass for a solver
/// with no cutoff at all.
///
/// # It stops at `1.4 f_c`, and the reason is the boundary rather than the scheme
///
/// [`Boundary::Open`] is Mur's condition, whose one coefficient is `(cΔt−Δ)/(cΔt+Δ)` — tuned to `c`.
/// A guided mode's phase velocity is `c/√(1−(f_c/f)²)`, which is **2c** at `1.25 f_c`, so the
/// absorbing faces are worst matched exactly where this measurement is most interesting. Their
/// reflection then travels back under the stations and corrupts the phase fit: measured, `β` comes
/// out 15% low at `1.25 f_c` while it is 2.3% at `1.4` and 0.1% at `2.1`.
///
/// Truncating the march before each reflection returns is what recovers the rest of the band, and
/// it is why the run is 0.8 of a guide length rather than two. Below `1.4 f_c` no window works: the
/// near end's reflection arrives at three and a half times the direct path and the direct pulse is
/// itself slow, so the two overlap.
///
/// What would fix it is a boundary matched to the guide, or a perfectly matched layer. Both are
/// larger than this crate has, and the honest thing is to say where the limit is and which side of
/// it the number came from.
#[test]
fn above_cutoff_the_phase_turns_more_slowly_than_in_free_space() {
    let mut worst: f64 = 0.0;
    let mut widest: f64 = 1.0;
    for (carrier, ratios) in [(1.5, [1.4, 1.5, 1.7]), (2.6, [2.1, 2.6, 3.1])] {
        // Four cells apart is under half a guide wavelength even at the top of the band; twenty of
        // them span four metres, which is several guide wavelengths at the bottom.
        let run = march(carrier, 1.4, ALONG / 3, 4, 20);
        for ratio in ratios {
            let f = run.cutoff * ratio;
            let beta = run.beta(f).abs();
            let k = 2.0 * std::f64::consts::PI * f / 299_792_458.0;
            let kc = 2.0 * std::f64::consts::PI * run.cutoff / 299_792_458.0;
            let closed = (k * k - kc * kc).sqrt();
            let off = (beta / closed - 1.0).abs();
            println!(
                "  f/f_c = {ratio:.2}: beta {beta:.3} against sqrt(k^2 - kc^2) {closed:.3} /m — \
                 off {:.2}%; lambda_g/lambda = {:.3}",
                off * 100.0,
                k / beta
            );
            worst = worst.max(off);
            widest = widest.max(k / beta);
        }
    }
    assert!(
        worst < 0.04,
        "the dispersion relation holds across the band: worst {:.3}%",
        worst * 100.0
    );
    println!("  and near cutoff the phase outran light by {widest:.3}x");
    assert!(
        widest > 1.4,
        "near cutoff the guide wavelength is much longer than free space: {widest:.3}x"
    );
}

/// **Below cutoff the mode is evanescent, decaying at `√(k_c² − k²)` without anything absorbing.**
///
/// The half of the dispersion relation with no counterpart in free space. Nothing here is lossy —
/// vacuum, conducting walls, no conductivity anywhere — and the field still falls off exponentially,
/// because below the cutoff there is no real wavenumber left to travel with.
///
/// # This one needs a driven source, and that is the physics rather than the API
///
/// Every other test here sets a field and lets it go. That cannot work below cutoff: an evanescent
/// field is a **near field**, it does not travel, and with nothing driving it it decays to zero.
/// There is no steady state for it to settle into and so no spatial profile to fit — the first
/// attempt returned decay constants of 0.285, −0.891 and −0.102 where the closed form says 2.7, and
/// a *negative* decay is the measurement saying it was fitting noise.
///
/// So a sinusoid is impressed on one plane, the march runs until the near field has settled, and the
/// amplitude is read along the guide.
#[test]
fn below_cutoff_the_mode_decays_without_being_absorbed() {
    let c0 = 299_792_458.0;
    let dx = 1.0 / ACROSS as f64;
    let a = ACROSS as f64 * dx;
    let mut c = Cavity::new(
        "guide",
        (ACROSS, 4, 120),
        Length::from_si(dx),
        Medium::vacuum(),
    );
    // `ZLow` stays a conductor. Mur's condition is derived for a wave that travels and an
    // evanescent field does not, so an absorbing face behind the source acts as a source of its
    // own: the same measurement read -7%, +11%, -2.5% and -4.1% depending only on which stations
    // it used. A conductor there gives an image field with the same exponential shape, which
    // changes the amplitude and not the slope.
    c.set_boundary(Wall::ZHigh, Boundary::Open);
    let cutoff = c.cutoff_frequency((1, 0)).to_si();
    let dt = Time::from_si(c.courant_limit().to_si() * 0.5);

    let mut biases = Vec::new();
    for ratio in [0.45, 0.6, 0.75] {
        let mut g = c.clone();
        let f = cutoff * ratio;
        let omega = 2.0 * std::f64::consts::PI * f;
        let period = 1.0 / f;
        let pi = std::f64::consts::PI;
        let source = 40usize;

        // Long enough for the near field to settle: many periods, and many transits of the region
        // being measured.
        let settle = (16.0 * period / dt.to_si()) as usize;
        let mut bus = Exchange::new();
        // **Impressed after the step, not before.** The electric update rewrites every interior
        // face from the curl, so a source added first is erased by the very step it was meant to
        // drive — what survives is its effect on `H` for one half step, which is a different source
        // with a different phase. It measured 3.14 against a closed form of 2.81 that way.
        for n in 0..settle {
            let t = n as f64 * dt.to_si();
            g.step(Time::from_si(t), dt, &mut bus).expect("stable");
            g.impress(source, |x, _y| {
                0.02 * (pi * x / a).sin() * (omega * t).sin()
            });
        }
        // Then four whole periods, recorded. **Whole**: a transform over a window that is not an
        // integer number of cycles leaks between bins, and with the amplitude falling by `e⁻²` across
        // the stations that leakage is the difference between reading 13% off and 4%.
        let window = (4.0 * period / dt.to_si()).round() as usize;
        // Adjacent cells, close in. Widening the span to hold a fixed number of decay lengths was
        // tried and made it worse — 12.4% where the tight span gives 5.5% — because the far stations
        // are down among the boundary's leakage and a longer lever on a noisier point is not a
        // better fit.
        let probes: Vec<usize> = (0..14).map(|n| source + 4 + n).collect();
        let mut series: Vec<Vec<f64>> = vec![Vec::new(); probes.len()];
        for n in 0..window {
            let t = (settle + n) as f64 * dt.to_si();
            g.step(Time::from_si(t), dt, &mut bus).expect("stable");
            g.impress(source, |x, _y| {
                0.02 * (pi * x / a).sin() * (omega * t).sin()
            });
            for (s, p) in series.iter_mut().zip(&probes) {
                s.push(g.electric_at(ACROSS / 2, 2, *p).y);
            }
        }

        let stations: Vec<f64> = probes.iter().map(|p| *p as f64 * dx).collect();
        let logs: Vec<f64> = series
            .iter()
            .map(|s| transform(s, f, dt.to_si()).0.ln())
            .collect();
        let alpha = -slope(&stations, &logs);
        let k = omega / c0;
        let kc = 2.0 * std::f64::consts::PI * cutoff / c0;
        let closed = (kc * kc - k * k).sqrt();
        let off = (alpha / closed - 1.0).abs();
        println!(
            "  f/f_c = {ratio:.2}: alpha {alpha:.4} against sqrt(kc^2 - k^2) {closed:.4} /m — off \
             {:.2}%",
            off * 100.0
        );
        assert!(
            off < 0.08,
            "an evanescent mode decays at the closed form: {alpha:.4} against {closed:.4}"
        );
        biases.push(alpha / closed);
    }

    // Reported rather than asserted, because it is not one common bias: 0.4%, 5.5% and 4.2% at the
    // three frequencies. The residual is the source's own near field over the first cells plus the
    // absorbing faces' mismatch — the same limitation `above_cutoff` runs into, from the other side —
    // and a claim that the three errors were *equal* would be a claim about the boundary that is not
    // true.
    let spread = biases.iter().cloned().fold(0.0f64, f64::max)
        - biases.iter().cloned().fold(f64::MAX, f64::min);
    println!(
        "  the three biases span {:.2} percentage points",
        spread * 100.0
    );
}

/// **The phase outruns light and the pulse does not, and their product is exactly `c²`.**
///
/// `v_p = ω/β > c` and `v_g = dω/dβ < c`, with `v_p v_g = c²`. The first alarms people and should
/// not: a phase carries no information.
///
/// `v_g` is measured as a **derivative of the dispersion curve** — `Δω/Δβ` between two nearby
/// frequencies of the same march — rather than from the arrival of an envelope. A pulse's energy
/// centroid was the first attempt and it read 0.96c against a formula's 0.78c, because a pulse in a
/// dispersive guide is a superposition of group velocities and the fastest of them dominates a
/// centroid. `Δω/Δβ` is the group velocity of one frequency, which is what the identity is about.
#[test]
fn the_phase_velocity_and_the_group_velocity_multiply_to_c_squared() {
    let c0 = 299_792_458.0;
    let run = march(2.0, 1.4, ALONG / 3, 4, 20);
    let ratio = 2.0;
    let f = run.cutoff * ratio;
    let step = f * 0.05;

    let (b_lo, b_hi) = (run.beta(f - step).abs(), run.beta(f + step).abs());
    let v_group = 2.0 * std::f64::consts::PI * 2.0 * step / (b_hi - b_lo);
    let v_phase = 2.0 * std::f64::consts::PI * f / run.beta(f).abs();

    let kc = 2.0 * std::f64::consts::PI * run.cutoff / c0;
    let k = 2.0 * std::f64::consts::PI * f / c0;
    let closed_phase = 2.0 * std::f64::consts::PI * f / (k * k - kc * kc).sqrt();
    let closed_group = c0 * c0 / closed_phase;

    println!(
        "  v_p {:.4}c against the closed form {:.4}c",
        v_phase / c0,
        closed_phase / c0
    );
    println!(
        "  v_g {:.4}c against c^2/v_p = {:.4}c",
        v_group / c0,
        closed_group / c0
    );
    println!("  and v_p v_g / c^2 = {:.6}", v_phase * v_group / (c0 * c0));
    assert!(
        v_phase > c0,
        "the phase velocity exceeds c above cutoff: {:.4}c",
        v_phase / c0
    );
    assert!(
        v_group < c0,
        "and the group velocity does not: {:.4}c",
        v_group / c0
    );
    assert!(
        (v_phase * v_group / (c0 * c0) - 1.0).abs() < 0.05,
        "their product is c squared: {:.6}",
        v_phase * v_group / (c0 * c0)
    );
}

/// **The cutoff is the cross-section's, and nothing else's.**
///
/// `f_c = (c/2)√((m/a)² + (n/b)²)` — a property of the geometry. Worth stating because the guide's
/// *length* does not appear in it, and a reader who has just seen a cavity's three-index resonance
/// will expect it to.
#[test]
fn the_cutoff_is_a_property_of_the_cross_section() {
    let dx = 1.0 / ACROSS as f64;
    let short = Cavity::new("a", (ACROSS, 4, 40), Length::from_si(dx), Medium::vacuum());
    let long = Cavity::new("b", (ACROSS, 4, 400), Length::from_si(dx), Medium::vacuum());
    assert!(
        (short.cutoff_frequency((1, 0)).to_si() - long.cutoff_frequency((1, 0)).to_si()).abs()
            < 1e-6,
        "ten times the length is the same cutoff"
    );
    // And it is `c/2a` for the dominant mode.
    let a = ACROSS as f64 * dx;
    let closed = 299_792_458.0 / (2.0 * a);
    println!(
        "  a {:.3} m gives f_c {:.4} GHz against c/2a = {:.4} GHz",
        a,
        short.cutoff_frequency((1, 0)).to_si() / 1e9,
        closed / 1e9
    );
    assert!(
        (short.cutoff_frequency((1, 0)).to_si() / closed - 1.0).abs() < 1e-12,
        "the dominant mode's cutoff is c/2a"
    );
    // A narrower guide cuts off higher, which is the whole reason a waveguide has a band.
    let narrow = Cavity::new(
        "c",
        (ACROSS / 2, 4, 40),
        Length::from_si(dx),
        Medium::vacuum(),
    );
    assert!(
        (narrow.cutoff_frequency((1, 0)).to_si() / short.cutoff_frequency((1, 0)).to_si() - 2.0)
            .abs()
            < 1e-12,
        "half the width is twice the cutoff"
    );
}
