//! Why a small room booms at one note, and why that note is not a note.
//!
//! ```text
//! cargo run --example room_modes            # numbers, checked
//! cargo run --example room_modes out.svg    # and a picture
//! ```
//!
//! A rectangular room with rigid walls resonates at
//!
//! ```text
//! f(nx, ny) = (c/2) √((nx/Lx)² + (ny/Ly)²)
//! ```
//!
//! for every pair of non-negative integers, and two things follow that a one-dimensional
//! duct never shows. **The series is not harmonic**: `f(1,1)` is `√2` times `f(1,0)` in a
//! square room, which is a tritone, so a room does not ring on a note — it rings on a chord
//! that is not one. And **the modes crowd together as they rise**, going as `f²` in two
//! dimensions, which is why a room's colouration is audible at the bottom of the spectrum
//! and inaudible at the top.
//!
//! The closed form is exact. What is checked here is that the *integration* agrees with it:
//! a mode is released, the pressure at a corner is watched, and the zero crossings are
//! counted. Nothing in the solver knows the formula.

use dualis::prelude::*;
use dualis_acoustic::Room;

mod common;
use common::svg::{diverging, document, rgb, ticks, Plot};
use common::{check, check_between, check_zero, heading};

/// A small room, of the size that has an audible problem.
const WIDTH: f64 = 4.4;
const HEIGHT: f64 = 3.1;
const CELLS: usize = 89;

fn room() -> Room {
    Room::of_air(
        "room",
        Length::from_si(WIDTH),
        Length::from_si(HEIGHT),
        CELLS,
    )
}

fn main() {
    let quiet = room();
    let (nx, ny) = quiet.cells();
    let speed = 343.0;

    heading("A 4.4 by 3.1 metre room in air");
    println!(
        "  {nx} by {ny} nodes, cells {:.1} mm across",
        quiet.width().to_si() / (nx - 1) as f64 * 1e3
    );
    // The three lowest modes, from the closed form. The axial pair first, then the one that
    // uses both dimensions at once.
    for (a, b) in [(1u32, 0u32), (0, 1), (1, 1)] {
        let f = quiet.mode_frequency(a, b).to_si();
        let expected = speed / 2.0
            * ((a as f64 / WIDTH).powi(2) + (b as f64 / quiet.height().to_si()).powi(2)).sqrt();
        check(&format!("mode ({a},{b})"), f, expected, 1e-12, "Hz");
    }
    // Not a harmonic series, and that is the whole difference from a pipe. The ratio is
    // sqrt(1 + (Lx/Ly)^2) — irrational for any room, and √2 only if it happens to be square.
    let ratio = quiet.mode_frequency(1, 1).to_si() / quiet.mode_frequency(1, 0).to_si();
    check(
        "f(1,1) / f(1,0)  -- a pipe would give exactly 2",
        ratio,
        (1.0 + (WIDTH / quiet.height().to_si()).powi(2)).sqrt(),
        1e-12,
        "x",
    );
    println!("  irrational, so the two lowest resonances are not in tune with each other");

    // And they crowd: counting how many lie below a frequency should go as f squared.
    let below = |hz: f64| quiet.modes_below(Frequency::from_si(hz)).len() as f64;
    let (low, high) = (below(100.0), below(200.0));
    println!("  {low:.0} modes below 100 Hz, {high:.0} below 200");
    check_between(
        "  quadrupling with the frequency",
        high / low,
        2.5,
        6.0,
        "x",
    );

    heading("The integration converges on it -- first order, which is a defect");
    // Release mode (1,1), time the zero crossings at a corner, and compare. Nothing in the
    // solver knows the formula, so this is an independent measurement of the same number.
    //
    // It comes out low, and refining the grid halves the error each time rather than
    // quartering it. A scheme that is second order in the interior converging at first order
    // means the *boundary* is first order, and here is why: a wall node's control volume is
    // half a cell wide, but `Room` divides its divergence by the whole `dx`. The walls come
    // out twice as heavy as they should be and every mode reads flat. It is a real defect in
    // the domain rather than a property of the discretisation, and it is asserted here as
    // what it is instead of being hidden behind a loose tolerance.
    let measure = |cells: usize| -> (f64, f64, Vec<(f64, f64)>) {
        let reference = Room::of_air(
            "room",
            Length::from_si(WIDTH),
            Length::from_si(HEIGHT),
            cells,
        );
        let exact = reference.mode_frequency(1, 1).to_si();
        let mut sound = Room::of_air(
            "room",
            Length::from_si(WIDTH),
            Length::from_si(HEIGHT),
            cells,
        )
        .released_in_mode(1, 1, Pressure::from_si(1.0));
        let dt = Time::from_si(sound.max_stable_dt(Time::ZERO).to_si() * 0.9);
        let corner = LengthVec::ZERO;
        let mut previous = sound.at(corner, Time::ZERO);
        let mut elapsed = 0.0;
        let mut bus = Exchange::new();
        let mut crossings: Vec<f64> = Vec::new();
        let mut trace: Vec<(f64, f64)> = vec![(0.0, previous)];
        while elapsed < 0.25 {
            sound
                .step(Time::ZERO, dt, &mut bus)
                .expect("inside the CFL limit");
            let now = sound.at(corner, Time::ZERO);
            elapsed += dt.to_si();
            if previous.signum() != now.signum() {
                // Interpolated rather than counted. Counting whole crossings quantises the
                // frequency by one over twice the window -- 2 Hz here, which would swamp the
                // effect being measured.
                let fraction = previous / (previous - now);
                crossings.push(elapsed - dt.to_si() * (1.0 - fraction));
            }
            previous = now;
            if trace.len() < 4000 {
                trace.push((elapsed, now));
            }
        }
        let (first, last) = (crossings[0], crossings[crossings.len() - 1]);
        let f = (crossings.len() - 1) as f64 / (2.0 * (last - first));
        (f, exact, trace)
    };

    let (coarse, exact_coarse, _) = measure(45);
    let (medium, exact_medium, trace) = measure(89);
    let (fine, exact_fine, _) = measure(177);
    let error = |m: f64, e: f64| (1.0 - m / e).abs();
    for (cells, m, e) in [
        (45, coarse, exact_coarse),
        (89, medium, exact_medium),
        (177, fine, exact_fine),
    ] {
        println!(
            "  {cells:>4} cells: {m:8.3} Hz against {e:8.3} Hz    {:+.2}%",
            (m / e - 1.0) * 100.0
        );
    }
    // Always low, never high: the walls only ever add weight.
    assert!(coarse < exact_coarse && medium < exact_medium && fine < exact_fine);
    // And the error halves rather than quartering, which is the diagnosis.
    check(
        "error ratio, 45 to 89 cells  (2 = first order)",
        error(coarse, exact_coarse) / error(medium, exact_medium),
        2.0,
        0.05,
        "x",
    );
    check(
        "error ratio, 89 to 177 cells",
        error(medium, exact_medium) / error(fine, exact_fine),
        2.0,
        0.05,
        "x",
    );
    check_between(
        "  so 89 cells is off by",
        error(medium, exact_medium) * 100.0,
        1.0,
        2.0,
        "%",
    );

    heading("And the energy stays put");
    // Rigid walls and no loss, so the scheme's own conserved quantity should not drift over
    // several thousand steps. Judged against the energy present, since the total is
    // constant rather than zero.
    let mut settled = room().released_in_mode(2, 1, Pressure::from_si(1.0));
    let dt = Time::from_si(settled.max_stable_dt(Time::ZERO).to_si() * 0.9);
    let mut bus = Exchange::new();
    let started = settled.energy().to_si();
    for _ in 0..4000 {
        settled.step(Time::ZERO, dt, &mut bus).unwrap();
    }
    check_zero(
        "energy drift over 4000 steps",
        settled.energy().to_si() - started,
        started,
        1e-9,
        "J",
    );

    heading("What the field reports, read through ScalarField");
    let shape = room().released_in_mode(1, 1, Pressure::from_si(1.0));
    let h = Length::from_si(WIDTH / (nx - 1) as f64);
    // A rigid wall's boundary condition is exactly zero normal gradient, and on this
    // node-centred grid a node sits on the wall — so it is zero to the last bit. Bar1D
    // reports a nonzero slope at *its* insulated end, because its grid is cell-centred and
    // the nearest sample is half a cell inside. Same physics, different sampling.
    check(
        "dp/dx on the left wall",
        shape.gradient(LengthVec::m(0.0, 1.0, 0.0), Time::ZERO, h).x,
        0.0,
        1e-12,
        "Pa/m",
    );
    check(
        "dp/dy on the floor",
        shape.gradient(LengthVec::m(1.0, 0.0, 0.0), Time::ZERO, h).y,
        0.0,
        1e-12,
        "Pa/m",
    );
    // The curvature carries the mode frequency: c2 lap(p) = -omega2 p. Two routes to one
    // number, sharing no code.
    //
    // Sampled exactly on a node, and that matters: `at` interpolates between nodes while the
    // derivatives snap to the nearest one, so asking for both at an arbitrary point divides a
    // curvature computed here by a value computed slightly over there. It cost 1.2% before I
    // noticed, which is the same size as the boundary defect above and would have been easy
    // to blame on it.
    let dx = h.to_si();
    let probe = LengthVec::from_si(glam::DVec3::new(13.0 * dx, 19.0 * dx, 0.0));
    let omega =
        -speed * speed * shape.laplacian(probe, Time::ZERO, h) / shape.at(probe, Time::ZERO);
    check(
        "mode (1,1) from the field's curvature",
        omega.sqrt() / std::f64::consts::TAU,
        quiet.mode_frequency(1, 1).to_si(),
        1e-3,
        "Hz",
    );
    println!("  the curvature is a spatial statement, so the wall defect does not touch it");

    let Some(path) = common::output_path() else {
        println!("\npass a path to write an SVG, e.g. `cargo run --example room_modes out.svg`");
        return;
    };
    common::write(&path, &draw(&trace, quiet.mode_frequency(1, 1).to_si()));
}

/// Four mode shapes and the corner trace underneath them.
fn draw(trace: &[(f64, f64)], f11: f64) -> String {
    let (w, h) = (900.0, 470.0);
    let mut parts = Vec::new();

    let modes = [(1u32, 0u32), (0, 1), (1, 1), (2, 1)];
    for (k, (a, b)) in modes.into_iter().enumerate() {
        let shape = room().released_in_mode(a, b, Pressure::from_si(1.0));
        let f = shape.mode_frequency(a, b).to_si();
        let x = 56.0 + k as f64 * 214.0;
        // Aspect-correct: the room is wider than it is tall, and drawing it square would
        // misrepresent which direction a mode runs in.
        let (pw, ph) = (186.0, 186.0 * HEIGHT / WIDTH);
        let mut panel = Plot::new(w, h, (0.0, WIDTH), (0.0, HEIGHT)).viewport(x, 56.0, pw, ph);
        // The panels are 186 px wide, so this is roughly three pixels a cell -- blocky up
        // close, but a mode is a large-scale shape and the file is a quarter the size.
        let n = 52;
        panel.raster(n, n, (0.0, WIDTH), (0.0, HEIGHT), |i, j| {
            let p = LengthVec::m(
                WIDTH * (i as f64 + 0.5) / n as f64,
                HEIGHT * (j as f64 + 0.5) / n as f64,
                0.0,
            );
            // Diverging, anchored at zero. A pressure is signed, and a ramp built for
            // temperature would put the midpoint wherever this frame's extremes fell.
            diverging(shape.at(p, Time::ZERO))
        });
        panel.label(
            x + pw / 2.0,
            56.0 + ph + 22.0,
            &format!("({a},{b})   {f:.1} Hz"),
            12.5,
            "#3a3a3a",
            "middle",
        );
        parts.push(panel.into_body());
    }

    // The corner pressure over time, with the closed-form period marked. Two independent
    // things agreeing is the picture.
    let span = trace.last().map(|(t, _)| *t).unwrap_or(1.0).min(0.06);
    let mut wave =
        Plot::new(w, h, (0.0, span * 1e3), (-1.15, 1.15)).viewport(56.0, 320.0, 788.0, 104.0);
    wave.axes(
        &ticks(0.0, span * 1e3, 8),
        &[-1.0, 0.0, 1.0],
        |v| format!("{v:.0}"),
        |v| format!("{v:.0}"),
    );
    let period_ms = 1000.0 / f11;
    let mut marker = period_ms;
    while marker < span * 1e3 {
        wave.polyline([(marker, -1.15), (marker, 1.15)], "#00000028", 1.0);
        marker += period_ms;
    }
    wave.polyline(
        trace
            .iter()
            .filter(|(t, _)| *t <= span)
            .map(|(t, p)| (t * 1e3, *p)),
        &rgb(28, 76, 168),
        1.8,
    );
    wave.label(
        56.0,
        24.0,
        "Rigid-wall modes of a 4.4 x 3.1 m room, and the corner pressure of mode (1,1)",
        15.0,
        "#1b1b1b",
        "start",
    );
    wave.label(
        844.0,
        24.0,
        &format!("grid lines every 1/f(1,1) = {period_ms:.1} ms"),
        12.0,
        "#6a6a6a",
        "end",
    );
    wave.label(
        450.0,
        452.0,
        "milliseconds -- the integration crosses the grid lines it was never told about",
        11.5,
        "#6a6a6a",
        "middle",
    );
    parts.push(wave.into_body());

    document(w, h, parts)
}
