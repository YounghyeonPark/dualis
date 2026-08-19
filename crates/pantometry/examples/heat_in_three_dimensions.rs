//! A point of heat in a block, and the exponent that says how many dimensions it is spreading in.
//!
//! ```text
//! cargo run --example heat_in_three_dimensions            # numbers, checked
//! cargo run --example heat_in_three_dimensions out.svg    # and a picture
//! ```
//!
//! An instantaneous point source of `Q` joules in an infinite medium has an exact solution, and
//! it is the cleanest closed form three-dimensional conduction has:
//!
//! ```text
//!   T(r,t) − T₀ = Q / (ρc·(4παt)^{d/2}) · exp(−r²/4αt)
//! ```
//!
//! The peak sits at `r = 0` and falls as `t^{−d/2}`, where **d is the number of dimensions the
//! heat can spread in**. A bar gives `−1/2`, a plate `−1`, a block `−3/2`. That exponent is the
//! sharpest possible statement that this is a three-dimensional simulation and not a
//! one-dimensional one with more arrays: nothing 1D can produce it, and no tolerance choice makes
//! it appear.
//!
//! # The window, which is most of the difficulty
//!
//! The closed form is for a **point** source in an **infinite** medium, and the simulation has
//! neither. So it only applies in between:
//!
//! ```text
//!   too early   √(4αt) ≲ Δx     the source is still one cell, not a point
//!   too late    √(4αt) ≳ L/4    the walls reflect, and the block is not infinite
//! ```
//!
//! Both ends are visible in the numbers below rather than asserted away, because a fit reported
//! without its window is a fit that can be made to say anything.
//!
//! # Why the step is half the stability limit, which is the thing worth taking away
//!
//! `Solid3D::max_stable_dt` reports `Δx²/6α`, and at **exactly** that step the sharpest mode the
//! grid can hold has an amplification factor of `−1`: it flips sign every step and never decays.
//! That is what marginal stability means, and it is not a bug — the scheme is stable there, which
//! is all the limit claims.
//!
//! It is not *accurate* there. A point source is the sharpest initial condition possible, so it
//! excites that mode as hard as anything can, and the undamped oscillation rides on the answer.
//! Measured against the closed form at `t = 0.1 s`:
//!
//! ```text
//!   step            peak      against the closed form
//!   the limit    1.969 K            1.96×
//!   half         1.025 K            1.005×
//!   a quarter    1.037 K            1.017×
//! ```
//!
//! A factor of two, from a scheme that never diverges and whose conservation audit is exact to
//! the last bit. **Stable is not the same as accurate**, and a stability limit is the only one of
//! the two a domain can compute for you.

use pantometry::prelude::*;

mod common;
use common::svg::{heat, rgb, ticks, Plot};
use common::{check, check_between, check_zero, heading};

/// Cells on a side. Odd, so there is a centre cell to put the heat in.
const N: usize = 21;
/// One millimetre cells, so the block is 21 mm across.
const DX: f64 = 1e-3;
/// Joules delivered to one cell at `t = 0`.
const Q: f64 = 2.0;
/// The fraction of `max_stable_dt` to run at. See the module docs.
const STEP_FRACTION: f64 = 0.5;

fn main() {
    let aluminium = Substance::aluminium_6061();
    let alpha = aluminium.diffusivity().expect("aluminium conducts").to_si();
    let rho_c = aluminium
        .heat_capacity(Volume::from_si(1.0))
        .expect("and has a specific heat")
        .to_si();

    let mut block = Solid3D::new(
        "block",
        aluminium,
        (N, N, N),
        Length::from_si(DX),
        Temperature::celsius(20.0),
    );
    let ambient = Temperature::celsius(20.0).to_si();
    block.deposit(N / 2, N / 2, N / 2, Energy::from_si(Q));

    let limit = block.max_stable_dt(Time::from_si(0.0));
    let dt = Time::from_si(limit.to_si() * STEP_FRACTION);

    heading("A point of heat in 21 mm of aluminium");
    println!(
        "  {:<28} {:.4} mm cells, {} of them",
        "grid",
        DX * 1e3,
        N * N * N
    );
    println!(
        "  {:<28} {:.4} ms, run at {:.0}% of it",
        "stability limit",
        limit.to_si() * 1e3,
        STEP_FRACTION * 100.0
    );
    println!("  {:<28} {:.1} J into one cell", "delivered", Q);

    // The window the closed form is about, in seconds.
    let span = N as f64 * DX;
    let t_min = (2.5 * DX).powi(2) / (4.0 * alpha);
    let t_max = (span / 4.0).powi(2) / (4.0 * alpha);
    println!(
        "  {:<28} {:.1} ms to {:.1} ms  (√4αt from {:.1} to {:.1} mm)",
        "closed form applies",
        t_min * 1e3,
        t_max * 1e3,
        2.5 * DX * 1e3,
        span / 4.0 * 1e3
    );

    // ---------------------------------------------------------------- run
    let mut bus = Exchange::new();
    let mut t = 0.0;
    let mut trace: Vec<(f64, f64, f64)> = Vec::new();
    // The peak-to-mean gap while it still exists. By the end of the run the block has
    // levelled and the ratio is 1, which is true and says nothing about what 3D buys.
    let mut gap_in_window = (0.0, 0.0, 0.0);
    let mut slices: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut want_slice = [t_min, 0.5 * (t_min + t_max), t_max].to_vec();

    // Run until the diffusion length reaches the whole block, so the walls are unambiguously
    // in the answer by the end. Stopping at 1.6x the window looked like enough and was not: the
    // late fit came out at -1.503, better than the fit *inside* the window, because at that point
    // the walls had barely been reached.
    let t_end = span.powi(2) / (4.0 * alpha);
    while t < t_end {
        block
            .step(Time::from_si(t), dt, &mut bus)
            .expect("half the limit is comfortably stable");
        t += dt.to_si();

        let peak = block.peak_temperature().to_si() - ambient;
        // The closed form at this instant, computed here and not by the domain.
        let gaussian = Q / (rho_c * (4.0 * std::f64::consts::PI * alpha * t).powf(1.5));
        trace.push((t, peak, gaussian));
        if t <= t_max {
            gap_in_window = (t, peak, block.mean_temperature().to_si() - ambient);
        }

        if want_slice.first().is_some_and(|w| t >= *w) {
            want_slice.remove(0);
            let mid = N / 2;
            let plane: Vec<f64> = (0..N)
                .flat_map(|j| (0..N).map(move |i| (i, j)).collect::<Vec<_>>())
                .map(|(i, j)| block.temperature_at(i, j, mid).to_si() - ambient)
                .collect();
            slices.push((t, plane));
        }
    }

    // ---------------------------------------------------------------- the exponent
    let inside: Vec<&(f64, f64, f64)> = trace
        .iter()
        .filter(|(t, _, _)| *t >= t_min && *t <= t_max)
        .collect();
    let slope = fit_log_slope(inside.iter().map(|(t, p, _)| (*t, *p)));

    heading("The exponent is the dimensionality");
    println!(
        "  {:<28} {} samples between {:.1} and {:.1} ms",
        "fitted over",
        inside.len(),
        t_min * 1e3,
        t_max * 1e3
    );
    check(
        "peak decays as t^d/2, d = 3",
        slope,
        -1.5,
        0.02,
        "(exponent)",
    );
    println!(
        "  {:<28} a bar would give −0.5, a plate −1.0",
        "for contrast"
    );

    // Outside the window it drifts, in both directions, and that is the window being real.
    let early = fit_log_slope(
        trace
            .iter()
            .filter(|(t, _, _)| *t < t_min)
            .map(|(t, p, _)| (*t, *p)),
    );
    // "Late" means the diffusion length is past half the block, so the walls are genuinely in it.
    let t_late = (span / 2.0).powi(2) / (4.0 * alpha);
    let late = fit_log_slope(
        trace
            .iter()
            .filter(|(t, _, _)| *t > t_late)
            .map(|(t, p, _)| (*t, *p)),
    );
    println!(
        "  {:<28} {early:+.3}  — the source is still one cell",
        "before the window"
    );
    println!(
        "  {:<28} {late:+.3}  — the walls are reflecting",
        "after the window"
    );
    // Each end must be *clearly* further from −3/2 than the fit inside, or the window is
    // decoration. Absolute deviations rather than a ratio, so the numbers read directly.
    check_between(
        "inside the window",
        (slope + 1.5).abs(),
        0.0,
        0.04,
        "from −3/2",
    );
    check_between("before it", (early + 1.5).abs(), 0.08, 10.0, "from −3/2");
    check_between("after it", (late + 1.5).abs(), 0.08, 10.0, "from −3/2");

    // ---------------------------------------------------------------- the amplitude
    heading("And the amplitude, against the same closed form");
    for (t, peak, gaussian) in inside.iter().step_by(inside.len().max(4) / 4) {
        println!(
            "  t = {:>6.1} ms   peak {:>8.4} K   closed form {:>8.4} K   ratio {:.4}",
            t * 1e3,
            peak,
            gaussian,
            peak / gaussian
        );
    }
    let worst = inside
        .iter()
        .map(|(_, p, g)| (p / g - 1.0).abs())
        .fold(0.0f64, f64::max);
    check_between("worst amplitude error", worst * 100.0, 0.0, 6.0, "%");

    // ---------------------------------------------------------------- conservation
    heading("What a lumped model would say instead");
    let capacity = Substance::aluminium_6061()
        .heat_capacity(block.volume())
        .expect("aluminium has a specific heat")
        .to_si();
    let uniform = Q / capacity;
    let mean = block.mean_temperature().to_si() - ambient;
    let peak = block.peak_temperature().to_si() - ambient;
    println!(
        "  {:<28} {uniform:.5} K — every cell, no gradient",
        "lumped rise"
    );
    // **The tolerance is the rounding, and it is worked out rather than tried.** The state is
    // updated once per cell per step — 9261 cells over about 1300 steps — and each update sums
    // seven terms. Accumulated float error over `k` roundings grows like `√k·ε`, which here is
    // `√(1300 × 9261) × 2.2e-16 ≈ 8e-13`. Measured 1.1e-12, the same order.
    //
    // 1e-12 was the first number written here and it failed by 10%, which is the right way round
    // for a tolerance to be wrong: a guess that is slightly too tight announces itself, and one
    // that is slightly too loose does not.
    check("the block holds every joule", mean, uniform, 1e-11, "K");
    let (gap_t, gap_peak, gap_mean) = gap_in_window;
    println!(
        "  {:<28} {gap_peak:.5} K against a mean of {gap_mean:.5} K, at t = {:.0} ms",
        "the peak really there",
        gap_t * 1e3
    );
    println!(
        "  {:<28} {:.0}× — and a lumped model reports that ratio as exactly 1",
        "so the gap is",
        gap_peak / gap_mean
    );
    println!(
        "  {:<28} {peak:.5} K, {:.3}× the mean — the block has levelled",
        "by the end of the run",
        peak / mean
    );
    check_zero(
        "insulated, so nothing left",
        (mean - uniform) * capacity,
        Q,
        1e-11,
        "J",
    );

    // ---------------------------------------------------------------- picture
    if let Some(path) = common::output_path() {
        common::write(&path, &draw(&trace, &slices, t_min, t_max, slope));
    }
}

/// Least-squares slope of `ln(peak)` against `ln(t)` — the exponent, fitted rather than eyeballed.
fn fit_log_slope(points: impl Iterator<Item = (f64, f64)>) -> f64 {
    let pts: Vec<(f64, f64)> = points
        .filter(|(t, p)| *t > 0.0 && *p > 0.0)
        .map(|(t, p)| (t.ln(), p.ln()))
        .collect();
    if pts.len() < 2 {
        return f64::NAN;
    }
    let n = pts.len() as f64;
    let (sx, sy) = pts.iter().fold((0.0, 0.0), |(a, b), (x, y)| (a + x, b + y));
    let (mx, my) = (sx / n, sy / n);
    let num: f64 = pts.iter().map(|(x, y)| (x - mx) * (y - my)).sum();
    let den: f64 = pts.iter().map(|(x, _)| (x - mx).powi(2)).sum();
    num / den
}

/// The decay on log-log with the exact slope beside it, and three mid-plane slices.
fn draw(
    trace: &[(f64, f64, f64)],
    slices: &[(f64, Vec<f64>)],
    t_min: f64,
    t_max: f64,
    slope: f64,
) -> String {
    let (w, h) = (960.0, 660.0);
    let lo_t = trace.first().map_or(1e-3, |(t, _, _)| *t).ln();
    let hi_t = trace.last().map_or(1.0, |(t, _, _)| *t).ln();
    let lo_p = trace
        .iter()
        .map(|(_, p, _)| *p)
        .fold(f64::MAX, f64::min)
        .max(1e-6)
        .ln();
    let hi_p = trace.iter().map(|(_, p, _)| *p).fold(0.0f64, f64::max).ln();

    let mut plot =
        Plot::new(w, h, (lo_t, hi_t), (lo_p, hi_p)).viewport(74.0, 58.0, w - 130.0, 330.0);
    plot.title("A point of heat in a block: the peak falls as t^(-3/2)");
    plot.caption(&format!(
        "fitted {slope:.3} against the exact -1.5, inside the window where the source is a point and the block is infinite"
    ));

    // The window, shaded by its edges.
    for (edge, label) in [(t_min, "source is a point"), (t_max, "walls reflect")] {
        plot.polyline(
            [(edge.ln(), lo_p), (edge.ln(), hi_p)],
            &rgb(150, 150, 150),
            1.0,
        );
        plot.label(edge.ln(), hi_p, label, 10.0, &rgb(110, 110, 110), "middle");
    }

    plot.polyline(
        trace.iter().map(|(t, _, g)| (t.ln(), g.ln())),
        &rgb(150, 160, 175),
        3.0,
    );
    plot.polyline(
        trace.iter().map(|(t, p, _)| (t.ln(), p.ln())),
        &rgb(255, 82, 33),
        2.0,
    );
    plot.axes(
        &ticks(lo_t, hi_t, 5),
        &ticks(lo_p, hi_p, 5),
        |v| format!("{:.0} ms", v.exp() * 1e3),
        |v| format!("{:.2} K", v.exp()),
    );
    plot.label(
        lo_t + (hi_t - lo_t) * 0.62,
        hi_p - (hi_p - lo_p) * 0.12,
        "grey: the closed form   orange: the simulation",
        11.0,
        &rgb(90, 90, 90),
        "start",
    );
    let mut parts = vec![plot.into_body()];

    // Three mid-plane slices, on one colour scale so they can be compared.
    let scale = slices
        .first()
        .map_or(1.0, |(_, v)| v.iter().fold(0.0f64, |m, x| m.max(*x)));
    for (k, (t, plane)) in slices.iter().enumerate() {
        let x = 74.0 + k as f64 * 280.0;
        let mut p =
            Plot::new(w, h, (0.0, N as f64), (0.0, N as f64)).viewport(x, 452.0, 250.0, 170.0);
        p.raster(N, N, (0.0, N as f64), (0.0, N as f64), |i, j| {
            heat(plane[j * N + i] / scale.max(1e-30))
        });
        p.caption(&format!("z-midplane at t = {:.1} ms", t * 1e3));
        parts.push(p.into_body());
    }
    let mut foot = Plot::new(w, h, (0.0, 1.0), (0.0, 1.0)).viewport(74.0, 636.0, w - 130.0, 1.0);
    foot.footnote(
        "one colour scale across the three slices, so the spot is seen to spread and fade rather than be renormalised",
    );
    parts.push(foot.into_body());
    common::svg::document(w, h, parts)
}
