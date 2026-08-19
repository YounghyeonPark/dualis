//! What a perfect lens does to a point.
//!
//! ```text
//! cargo run --example airy_pattern            # numbers, checked
//! cargo run --example airy_pattern out.svg    # and a picture
//! ```
//!
//! The diffraction limit is not a blur that better manufacturing removes. A circular
//! aperture turns a point at infinity into a disc with rings around it, and that is what a
//! *flawless* instrument produces — everything an optical designer does is spent getting
//! close to this pattern, never past it.
//!
//! Three panels, because they are three views of one piece of physics: the image of a point,
//! the fraction of its light inside a given radius, and the contrast surviving at a given
//! spatial frequency. The third comes from the first by an autocorrelation, and it is the one
//! printed on a datasheet.
//!
//! Every number here comes from the Bessel series in `pantometry-optics`, checked against the
//! values an optics text quotes — 83.8% inside the first dark ring, 0.391 at half cutoff —
//! rather than against a table hard-coded somewhere.

use pantometry::prelude::*;
use pantometry_optics::diffraction::{
    airy_intensity, mtf_ideal, FIRST_AIRY_ZERO, SECOND_AIRY_ZERO,
};

mod common;
use common::svg::{document, heat, rgb, ticks, Plot};
use common::{check, check_between, heading};

/// f/4 in the green. Image-side numerical aperture is `1/2F`.
const LAMBDA: f64 = 550e-9;
const F_NUMBER: f64 = 4.0;
const NA: f64 = 0.5 / F_NUMBER;

fn main() {
    let lambda = Length::from_si(LAMBDA);

    heading("A perfect f/4 lens in the green");
    // The exact factor is the first zero of J1 over 2 pi, which is 0.60983 and not 0.61.
    // Both roundings are everywhere in the literature and they agree to three parts in ten
    // thousand, so it never matters -- but a tolerance of 1e-12 against "0.61" would fail,
    // and it is better to know which of the two numbers the code actually holds.
    let airy = airy_radius(lambda, NA);
    let exact = FIRST_AIRY_ZERO / std::f64::consts::TAU;
    check(
        "Airy radius  (J1 zero / 2pi) lambda / NA",
        airy.in_um(),
        exact * LAMBDA / NA * 1e6,
        1e-12,
        "um",
    );
    check(
        "  against the textbook 1.22 lambda F",
        airy.in_um() / (1.22 * LAMBDA * F_NUMBER * 1e6),
        1.0,
        4e-4,
        "x",
    );
    check(
        "Rayleigh limit  (the same distance)",
        rayleigh_limit(lambda, NA).in_um(),
        airy.in_um(),
        1e-12,
        "um",
    );
    // Past this a perfect lens transmits nothing at all — not faintly, nothing. There is no
    // deconvolution that gets it back.
    let cutoff = cutoff_frequency(lambda, NA);
    check(
        "cutoff  2 NA / lambda",
        cutoff / 1000.0,
        1.0 / (LAMBDA * F_NUMBER) / 1000.0,
        1e-12,
        "cy/mm",
    );

    heading("Where the light goes");
    // The two fractions every optics text quotes, from `1 - J0² - J1²` rather than a table.
    check(
        "energy inside the 1st dark ring",
        encircled_energy(FIRST_AIRY_ZERO) * 100.0,
        83.8,
        1e-3,
        "%",
    );
    check(
        "energy inside the 2nd",
        encircled_energy(SECOND_AIRY_ZERO) * 100.0,
        91.0,
        1e-3,
        "%",
    );
    // A sixth of a *perfect* image's light is outside the Airy disc. Not aberration — what
    // a finite aperture does.
    check(
        "energy in the rings, therefore",
        100.0 - encircled_energy(FIRST_AIRY_ZERO) * 100.0,
        16.2,
        6e-3,
        "%",
    );

    // Found by searching between the first two zeros rather than asserted at a remembered
    // position, so the number and its location are both the code's own.
    let (ring_v, ring_peak) = (0..2000)
        .map(|i| {
            let v =
                FIRST_AIRY_ZERO + (SECOND_AIRY_ZERO - FIRST_AIRY_ZERO) * (i as f64 + 0.5) / 2000.0;
            (v, airy_intensity(v))
        })
        .fold((0.0, 0.0), |best, c| if c.1 > best.1 { c } else { best });
    check_between(
        "brightest point of the 1st ring",
        ring_peak * 100.0,
        1.5,
        2.0,
        "%",
    );
    check_between(
        "  which sits at, in Airy radii",
        ring_v / FIRST_AIRY_ZERO,
        1.3,
        1.5,
        "",
    );

    heading("Contrast against spatial frequency");
    // The ideal MTF is a disc autocorrelated with itself: (2/pi)(acos s - s sqrt(1-s^2)).
    for (fraction, expected) in [(0.0, 1.0), (0.25, 0.685), (0.5, 0.391), (0.75, 0.1443)] {
        check(
            &format!("MTF at {:>3.0}% of cutoff", fraction * 100.0),
            mtf_at(fraction * cutoff, lambda, NA),
            expected,
            2e-3,
            "",
        );
    }
    check(
        "MTF at cutoff, and beyond",
        mtf_at(1.5 * cutoff, lambda, NA),
        0.0,
        1e-12,
        "",
    );

    heading("And what a quarter wave of error costs");
    // Rayleigh's rule, and where "diffraction limited" comes from: a quarter wave
    // peak-to-valley is about lambda/14 RMS, and it leaves a Strehl ratio near 0.8.
    let strehl = strehl_from_wavefront_error(1.0 / 14.0);
    check_between("Strehl at lambda/14 RMS", strehl, 0.78, 0.84, "");
    // Depth of focus goes as 1/NA², so doubling the aperture to halve the spot quarters the
    // depth. That trade is the reason a fast lens is hard to focus.
    let dof = depth_of_focus(lambda, NA, 1.0);
    let half_na = depth_of_focus(lambda, NA * 2.0, 1.0);
    check(
        "depth of focus at f/4",
        dof.in_um(),
        LAMBDA / (NA * NA) * 1e6,
        1e-12,
        "um",
    );
    check(
        "  and at f/2, four times shallower",
        dof.in_um() / half_na.in_um(),
        4.0,
        1e-12,
        "x",
    );

    let Some(path) = common::output_path() else {
        println!("\npass a path to write an SVG, e.g. `cargo run --example airy_pattern out.svg`");
        return;
    };
    common::write(&path, &draw(airy.in_um(), cutoff));
}

/// Three panels: the pattern as a raster, the encircled energy, and the MTF.
fn draw(airy_um: f64, cutoff: f64) -> String {
    let (w, h) = (900.0, 330.0);
    let extent = 3.2; // Airy radii from the centre to the edge of the raster

    // --- The pattern. Logarithmic, because the first ring is under 2% of the centre and a
    // linear map shows a white dot on black — the rings are the whole point.
    let mut pattern =
        Plot::new(w, h, (-extent, extent), (-extent, extent)).viewport(56.0, 54.0, 216.0, 216.0);
    let n = 96;
    pattern.raster(n, n, (-extent, extent), (-extent, extent), |ix, iy| {
        let centre = |k: usize| -extent + 2.0 * extent * (k as f64 + 0.5) / n as f64;
        let (x, y) = (centre(ix), centre(iy));
        let intensity = airy_intensity((x * x + y * y).sqrt() * FIRST_AIRY_ZERO);
        // Four decades of range mapped onto the ramp. Linear would show a white dot.
        heat(((intensity.max(1e-6).log10() + 4.0) / 4.0).clamp(0.0, 1.0))
    });
    pattern.axes(
        &[-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0],
        &[],
        |v| format!("{v:.0}"),
        |_| String::new(),
    );
    pattern.title("A point source through a perfect f/4 lens");
    pattern.caption(&format!(
        "Airy radius {airy_um:.2} um   cutoff {:.0} cy/mm",
        cutoff / 1000.0
    ));
    pattern.footnote("the pattern, log scale, in Airy radii");

    // --- Encircled energy, with the two quoted fractions marked where they fall.
    let mut energy = Plot::new(w, h, (0.0, 4.0), (0.0, 1.02)).viewport(340.0, 54.0, 216.0, 216.0);
    energy.axes(
        &ticks(0.0, 4.0, 4),
        &ticks(0.0, 1.0, 5),
        |v| format!("{v:.0}"),
        |v| format!("{:.0}%", v * 100.0),
    );
    energy.polyline(
        (0..=320).map(|i| {
            let r = i as f64 / 80.0;
            (r, encircled_energy(r * FIRST_AIRY_ZERO))
        }),
        &rgb(158, 40, 32),
        2.0,
    );
    for (radii, v) in [(1.0, FIRST_AIRY_ZERO), (2.0, SECOND_AIRY_ZERO)] {
        let e = encircled_energy(v);
        energy.polyline([(radii, 0.0), (radii, e)], "#00000038", 1.0);
        energy.text(
            radii + 0.12,
            e - 0.09,
            &format!("{:.1}%", e * 100.0),
            11.5,
            "#7a3a34",
            "start",
        );
    }
    energy.footnote("energy inside a radius, in Airy radii");

    // --- The MTF, with half cutoff marked because it is the number people compare.
    let mut mtf = Plot::new(w, h, (0.0, 1.0), (0.0, 1.02)).viewport(624.0, 54.0, 216.0, 216.0);
    mtf.axes(
        &ticks(0.0, 1.0, 4),
        &ticks(0.0, 1.0, 5),
        |v| format!("{v:.2}"),
        |v| format!("{v:.1}"),
    );
    mtf.polyline(
        (0..=240).map(|i| {
            let s = i as f64 / 240.0;
            (s, mtf_ideal(s))
        }),
        &rgb(48, 108, 186),
        2.0,
    );
    mtf.polyline([(0.5, 0.0), (0.5, mtf_ideal(0.5))], "#00000038", 1.0);
    mtf.text(
        0.54,
        mtf_ideal(0.5) - 0.08,
        "0.391",
        11.5,
        "#2a4a78",
        "start",
    );
    mtf.footnote("contrast, as a fraction of cutoff");

    document(
        w,
        h,
        [pattern.into_body(), energy.into_body(), mtf.into_body()],
    )
}
