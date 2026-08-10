//! The modes a floor plan does not have, and the count that grows as the cube of frequency.
//!
//! ```text
//! cargo run --example room_in_three_dimensions            # numbers, checked
//! cargo run --example room_in_three_dimensions out.svg    # and a picture
//! ```
//!
//! `cargo run --example room_modes` is the same room in two dimensions, and it is a good model of
//! a floor plan. What it cannot be is a model of a room, and the difference is not accuracy:
//!
//! ```text
//!   f(a,b,c) = (c/2)·√( (a/Lx)² + (b/Ly)² + (c/Lz)² )
//! ```
//!
//! Set `Lz` aside and the whole `c` family disappears. A 2.4 m ceiling puts the first
//! floor-to-ceiling mode at 71 Hz — squarely inside the range a room is judged on — and a
//! two-dimensional model does not have it *at all*. Nor any of the oblique modes, which need all
//! three indices at once and are most of them.
//!
//! # Why the count matters more than any single mode
//!
//! Modes below a frequency fill a **volume** of the `(a,b,c)` lattice rather than an area, so the
//! count grows as `f³` and not `f²`. That is why a real room's resonances crowd together into a
//! statistical hiss above the Schroeder frequency while a floor plan keeps them separable much
//! further up. A two-dimensional model is not a slightly sparser room; it is a room whose modes
//! never stop being countable.
//!
//! Checked against Weyl's estimate with its **surface and edge terms**, not the leading term
//! alone:
//!
//! ```text
//!   N(f) = (4π/3)·V·(f/c)³ + (π/4)·S·(f/c)² + (1/8)·L·(f/c)
//! ```
//!
//! For rigid walls that surface term is positive and large at these frequencies — the leading
//! term alone is 46% low at 300 Hz — and dropping it is the easy way to conclude the count is
//! wrong when it is the estimate that is incomplete.

use dualis::prelude::*;

mod common;
use common::svg::{diverging, rgb, ticks, Plot};
use common::{check, check_between, heading};

/// A small room: 4.4 m by 3.1 m, with a 2.4 m ceiling.
const WIDTH: f64 = 4.4;
const HEIGHT: f64 = 3.1;
const DEPTH: f64 = 2.4;
/// Nodes across the width. Height and depth are quantised to the same spacing.
const ACROSS: usize = 23;
/// The mode to release: oblique, so all three axes carry a velocity.
const MODE: (u32, u32, u32) = (1, 1, 1);

fn main() {
    let hall = Hall::of_air(
        "hall",
        Length::m(WIDTH),
        Length::m(HEIGHT),
        Length::m(DEPTH),
        ACROSS,
    );
    let (lx, ly, lz) = (
        hall.width().to_si(),
        hall.height().to_si(),
        hall.depth().to_si(),
    );
    let (nx, ny, nz) = hall.nodes();

    heading("A room with a ceiling");
    println!(
        "  {:<30} {lx:.3} x {ly:.3} x {lz:.3} m, quantised from {WIDTH} x {HEIGHT} x {DEPTH}",
        "size"
    );
    println!(
        "  {:<30} {nx} x {ny} x {nz} = {} nodes at {:.1} mm",
        "grid",
        nx * ny * nz,
        hall.spacing().to_si() * 1e3
    );

    // ---------------------------------------------------------------- the modes themselves
    heading("The closed form, and the family a floor plan loses");
    let exact = |a: f64, b: f64, c: f64| {
        343.0 / 2.0 * ((a / lx).powi(2) + (b / ly).powi(2) + (c / lz).powi(2)).sqrt()
    };
    for (mode, a, b, c, what) in [
        ((1, 0, 0), 1.0, 0.0, 0.0, "axial, along the width"),
        ((0, 1, 0), 0.0, 1.0, 0.0, "axial, across"),
        (
            (0, 0, 1),
            0.0,
            0.0,
            1.0,
            "axial, floor to ceiling — absent in 2D",
        ),
        ((1, 1, 0), 1.0, 1.0, 0.0, "tangential, in the floor plan"),
        (
            (1, 1, 1),
            1.0,
            1.0,
            1.0,
            "oblique — needs all three, absent in 2D",
        ),
    ] {
        let got = hall.mode_frequency(mode).to_si();
        println!("  {:<9} {:>8.2} Hz   {what}", format!("{mode:?}"), got);
        check(&format!("mode {mode:?}"), got, exact(a, b, c), 1e-12, "Hz");
    }
    let vertical = hall.mode_frequency((0, 0, 1)).to_si();
    check(
        "the vertical mode is c/2Lz",
        vertical,
        343.0 / (2.0 * lz),
        1e-12,
        "Hz",
    );
    check_between(
        "and it is inside the audible bass",
        vertical,
        60.0,
        85.0,
        "Hz",
    );

    // ---------------------------------------------------------------- the count
    heading("How many modes there are, which is the real difference");
    let pi = std::f64::consts::PI;
    let (v, surface, edges) = (
        lx * ly * lz,
        2.0 * (lx * ly + ly * lz + lz * lx),
        4.0 * (lx + ly + lz),
    );
    let weyl = |f: f64| {
        let k = f / 343.0;
        4.0 * pi / 3.0 * v * k.powi(3) + pi / 4.0 * surface * k.powi(2) + edges * k / 8.0
    };

    let mut counts = Vec::new();
    for f in [100.0, 150.0, 200.0, 300.0] {
        let n = hall.modes_below(Frequency::from_si(f)).len() as f64;
        counts.push((f, n));
        println!(
            "  below {f:>5.0} Hz   {n:>5.0} modes   Weyl {:>7.1}   leading term alone {:>7.1}",
            weyl(f),
            4.0 * pi / 3.0 * v * (f / 343.0).powi(3)
        );
    }
    let (f_hi, n_hi) = *counts.last().unwrap();
    check("against Weyl at 300 Hz", n_hi, weyl(f_hi), 0.1, "modes");
    check_between(
        "the leading term alone is not enough",
        n_hi / (4.0 * pi / 3.0 * v * (f_hi / 343.0).powi(3)),
        1.2,
        2.0,
        "x",
    );

    // Doubling the frequency, as a ratio so no constant has to be right for it to mean something.
    //
    // **Against Weyl's own ratio, not against 8.** The pure `f^3` asymptote gives 8 and a
    // two-dimensional room would give 4, but at these frequencies the positive surface term is
    // still a large part of the count, so the growth is genuinely sub-cubic — 5.82 by the
    // three-term estimate. Asserting 8 here would be asserting the asymptote against data that is
    // nowhere near it, and passing only because the band was wide enough to swallow the gap.
    let ratio = counts[3].1 / counts[1].1;
    let weyl_ratio = weyl(counts[3].0) / weyl(counts[1].0);
    println!(
        "  {:<30} {ratio:.2}x from 150 to 300 Hz",
        "doubling the frequency"
    );
    println!(
        "  {:<30} {weyl_ratio:.2}x — the f^3 asymptote is 8, and a floor plan would give 4",
        "Weyl predicts"
    );
    check(
        "the count grows the way Weyl says",
        ratio,
        weyl_ratio,
        0.08,
        "x",
    );
    check_between("which is well clear of a plane's f^2", ratio, 5.0, 8.0, "x");

    // ---------------------------------------------------------------- it rings at that frequency
    heading("And the integration rides the closed form");
    let f = hall.mode_frequency(MODE).to_si();
    let mut hall = hall.released_in_mode(MODE, Pressure::from_si(1.0));
    let dt = Time::from_si(hall.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
    let mut bus = Exchange::new();
    let (mut t, mut worst) = (0.0, 0.0f64);
    let mut trace = Vec::new();
    while t < 2.5 / f {
        hall.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
        let peak = hall.peak_pressure().to_si();
        let want = (2.0 * pi * f * t).cos().abs();
        worst = worst.max((peak - want).abs());
        trace.push((t, peak, want));
    }
    println!(
        "  {:<30} {f:.2} Hz, followed for {:.1} periods at {} nodes across",
        "the (1,1,1) mode",
        t * f,
        nx
    );
    check_between(
        "worst departure from |cos 2 pi f t|",
        worst,
        0.0,
        0.03,
        "Pa",
    );
    println!(
        "  {:<30} {:.3e} J, reported rather than hidden",
        "startup adjustment",
        hall.startup_adjustment().to_si().abs()
    );

    if let Some(path) = common::output_path() {
        let mid = nz / 2;
        let plane: Vec<f64> = (0..ny)
            .flat_map(|j| (0..nx).map(move |i| (i, j)).collect::<Vec<_>>())
            .map(|(i, j)| hall.pressure_at(i, j, mid).to_si())
            .collect();
        common::write(&path, &draw(&trace, &plane, nx, ny, f, &counts, &weyl));
    }
}

/// The peak against the closed form, the mode's shape at mid-height, and the mode count.
fn draw(
    trace: &[(f64, f64, f64)],
    plane: &[f64],
    nx: usize,
    ny: usize,
    f: f64,
    counts: &[(f64, f64)],
    weyl: &dyn Fn(f64) -> f64,
) -> String {
    let (w, h) = (960.0, 640.0);
    let t_end = trace.last().map_or(1.0, |(t, _, _)| *t);

    let mut top =
        Plot::new(w, h, (0.0, t_end * 1e3), (-0.05, 1.15)).viewport(74.0, 56.0, w - 130.0, 250.0);
    top.title("A room with a ceiling: the oblique (1,1,1) mode");
    top.polyline(
        trace.iter().map(|(t, _, c)| (t * 1e3, *c)),
        &rgb(150, 160, 175),
        3.0,
    );
    top.polyline(
        trace.iter().map(|(t, p, _)| (t * 1e3, *p)),
        &rgb(255, 82, 33),
        2.0,
    );
    top.axes(
        &ticks(0.0, t_end * 1e3, 6),
        &ticks(0.0, 1.0, 5),
        |v| format!("{v:.0} ms"),
        |v| format!("{v:.2}"),
    );
    top.caption(&format!(
        "grey |cos(2 pi f t)| at {f:.2} Hz from the closed form, orange the integration"
    ));

    let scale = plane.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-30);
    let mut left =
        Plot::new(w, h, (0.0, nx as f64), (0.0, ny as f64)).viewport(74.0, 372.0, 330.0, 210.0);
    left.raster(nx, ny, (0.0, nx as f64), (0.0, ny as f64), |i, j| {
        diverging(plane[j * nx + i] / scale)
    });
    left.caption("pressure at mid-height: one nodal line each way, and a third out of the page");

    let f_hi = counts.last().map_or(300.0, |(f, _)| *f);
    let n_hi = counts.last().map_or(1.0, |(_, n)| *n);
    let mut right =
        Plot::new(w, h, (0.0, f_hi), (0.0, n_hi * 1.1)).viewport(500.0, 372.0, w - 560.0, 210.0);
    right.polyline(
        (1..=60).map(|k| {
            let f = f_hi * k as f64 / 60.0;
            (f, weyl(f))
        }),
        &rgb(150, 160, 175),
        3.0,
    );
    right.polyline(counts.iter().map(|(f, n)| (*f, *n)), &rgb(255, 82, 33), 2.0);
    right.axes(
        &ticks(0.0, f_hi, 4),
        &ticks(0.0, n_hi * 1.1, 4),
        |v| format!("{v:.0} Hz"),
        |v| format!("{v:.0}"),
    );
    right.caption("modes below a frequency: counted, against Weyl with its surface term");

    common::svg::document(w, h, [top.into_body(), left.into_body(), right.into_body()])
}
