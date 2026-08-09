//! Why read noise is the first number on a camera's datasheet.
//!
//! ```text
//! cargo run --example detector_snr            # numbers, checked
//! cargo run --example detector_snr out.svg    # and a picture
//! ```
//!
//! Everything upstream of a detector is deterministic: a lamp has a spectrum, a filter passes
//! a fraction of it, a lens forms a spot. A detector is where that stops. Counting `N`
//! photons carries an uncertainty of `√N` and no care taken in the optics changes it.
//!
//! ```text
//! SNR = S / sqrt(S + D + R²)
//! ```
//!
//! The whole of low-light imaging is in the shape of that denominator, and it has two
//! regimes that behave completely differently. Below `S = R²` the noise is a constant, so
//! the SNR grows *linearly* with exposure and twice the time is twice the quality. Above it
//! the noise is `√S`, and twice the time buys 41%. The crossover for a 6-electron read noise
//! is 36 electrons — a startlingly small number, and the reason the quiet sensor wins in the
//! dark and makes no difference in daylight.
//!
//! The numbers below are not asserted from the formula alone. Frames are actually drawn from
//! the noise model and their sample variance compared with their sample mean, which is the
//! only way to find out whether the distribution is the one that was intended.

use dualis::prelude::*;
use dualis_optics::detector::{Detector, Exposure};
use dualis_optics::spectrum::Spectrum as Spec;

mod common;
use common::svg::{document, rgb, ticks, Plot};
use common::{check, check_between, heading};

/// Enough frames that a 1/sqrt(N) tolerance is a few tenths of a percent.
const FRAMES: usize = 200_000;

fn main() {
    let scientific = Detector::scientific_cmos();
    let consumer = Detector::consumer_cmos();

    heading("Two sensors, and the electron count where each stops being read-limited");
    // S = R² is where shot noise overtakes read noise. Below it, doubling the exposure
    // doubles the SNR; above it, doubling buys 41%.
    for (name, d) in [("scientific", &scientific), ("consumer", &consumer)] {
        check(
            &format!("{name}: crossover, R squared"),
            d.shot_noise_crossover(),
            d.read_noise * d.read_noise,
            1e-12,
            "e-",
        );
    }
    println!("  a 1.4 e- sensor is shot-limited from 2 electrons; a 6 e- one needs 36");

    heading("The two regimes, measured rather than restated");
    // Through the library's own Exposure rather than by writing S/sqrt(S+D+R2) again here.
    // An example that restates the formula it is demonstrating proves only that I can copy.
    let faint = |d: &Detector, s: f64| {
        Exposure {
            signal: s,
            dark: 0.0,
        }
        .snr(d.read_noise)
    };
    let r = consumer.read_noise * consumer.read_noise;
    check(
        "read-limited: 2x signal at S = R2/16",
        faint(&consumer, r / 8.0) / faint(&consumer, r / 16.0),
        2.0,
        0.04,
        "x",
    );
    // Shot-limited: the noise grows as the root, so the SNR does too.
    check(
        "shot-limited: 2x signal at S = 400 R2",
        faint(&consumer, 800.0 * r) / faint(&consumer, 400.0 * r),
        std::f64::consts::SQRT_2,
        0.01,
        "x",
    );
    // The ideal limit both approach: a perfect counter's SNR is exactly root N.
    let ideal = Detector::ideal();
    check(
        "an ideal counter at 10 000 electrons",
        faint(&ideal, 10_000.0),
        100.0,
        1e-12,
        "SNR",
    );

    heading("And the distribution itself, from 200 000 frames");
    // The claim is about a *distribution*, so it is checked by sampling. Variance equals
    // mean is the signature of a Poisson process, and a model with the right mean but the
    // wrong spread would pass every check above and fail this one.
    // Two hundred thousand frames at three rates, spread over eight threads.
    //
    // This used to draw from one `Rng` in a loop and take the variance as
    // `sum(k²)/N − mean²`. Both were worth changing, and the second more than the first: at a
    // mean of 900 that expression subtracts 1.6e11 from itself to arrive at 900, and throws
    // away most of the digits it needed. `Ensemble` folds with Welford inside a block and
    // Chan's merge between blocks, so nothing large is ever subtracted from nothing large.
    //
    // The parallelism is free of consequences because `Rng::for_index` gives frame `i` the same
    // draw wherever it runs — asserted below rather than assumed.
    let tol = 4.0 * (2.0 / FRAMES as f64).sqrt();
    for mean in [4.0f64, 45.0, 900.0] {
        let frames = Ensemble::new(0xD0_5E_11_A5, FRAMES as u64).with_threads(8);
        let e = frames
            .estimate(|_, mut rng| rng.poisson(mean) as f64)
            .expect("two hundred thousand frames");

        // Tolerance from the sample size, not from taste: the standard error of a variance
        // estimate is about sqrt(2/N), which is 0.32% here. Four of those is a fair bound.
        check(
            &format!("mean {mean:>5.0}: sampled mean"),
            e.mean,
            mean,
            tol,
            "e-",
        );
        check(
            &format!("mean {mean:>5.0}: variance equals it"),
            e.standard_deviation() * e.standard_deviation(),
            mean,
            tol,
            "e-",
        );

        // And the count is what a Monte Carlo result is meaningless without.
        assert_eq!(e.samples, FRAMES as u64);
    }

    heading("The same frames, on a different number of threads");
    // Not a physics claim — a claim about this library. A Monte Carlo drawing from a shared
    // generator gives a different answer on eight cores than on one, and the difference looks
    // exactly like statistical noise, so it is never investigated. Here it is compared on the
    // bits, which noise cannot survive.
    let one_frame = |_: u64, mut rng: Rng| rng.poisson(45.0) as f64;
    let sequential = Ensemble::new(0xD0_5E_11_A5, FRAMES as u64)
        .estimate(one_frame)
        .expect("frames");
    for threads in [2usize, 8, 32] {
        let parallel = Ensemble::new(0xD0_5E_11_A5, FRAMES as u64)
            .with_threads(threads)
            .estimate(one_frame)
            .expect("frames");
        assert_eq!(
            sequential.mean.to_bits(),
            parallel.mean.to_bits(),
            "{threads} threads moved the mean"
        );
        println!(
            "  {threads:>2} threads: mean {:.9} e-, identical to the bit",
            parallel.mean
        );
    }
    // The generator switches from inverse transform to a rounded normal at a mean of 30, and
    // the join is checked from both sides above: 4 is below it, 900 well above, 45 just past.

    heading("What that means for an exposure");
    // A real source rather than a bare rate: twenty attowatts of sunlight, integrated
    // against each sensor's own quantum efficiency. The two sensors see different rates
    // from the same star, which is part of the comparison rather than noise in it.
    let star = SpectralPower::new(
        Spec::blackbody(5800.0),
        Power::from_si(2.0e-17),
        VISIBLE_RANGE,
    );
    let target = 10.0;
    let rate_of = |d: &Detector| star.photon_rate_through(&d.quantum_efficiency).to_si();
    let t_sci = scientific.exposure_for_snr(&star, target);
    let t_con = consumer.exposure_for_snr(&star, target);
    println!(
        "  a 20 aW star: {:.1} e-/s on the scientific sensor, {:.1} on the consumer one",
        rate_of(&scientific),
        rate_of(&consumer)
    );
    check_between(
        "  detected rate, scientific",
        rate_of(&scientific),
        5.0,
        500.0,
        "e-/s",
    );
    println!("  target SNR {target:.0}");
    check_between("  scientific cmos needs", t_sci.to_si(), 0.2, 60.0, "s");
    check_between("  consumer cmos needs", t_con.to_si(), 1.0, 6000.0, "s");
    check_between(
        "  so the quiet sensor is faster by",
        t_con.to_si() / t_sci.to_si(),
        1.5,
        200.0,
        "x",
    );
    // And the reason is dark current as much as read noise: 20 e-/s against 0.5 means the
    // consumer sensor is fighting its own thermal signal before the star arrives.
    check_between(
        "  consumer dark current, for comparison",
        consumer.dark_current.to_si(),
        5.0,
        100.0,
        "e-/s",
    );

    let Some(path) = common::output_path() else {
        println!("\npass a path to write an SVG, e.g. `cargo run --example detector_snr out.svg`");
        return;
    };
    common::write(&path, &draw(&scientific, &consumer, &ideal));
}

/// Two panels: SNR against signal on log axes, and the noise budget against exposure.
fn draw(scientific: &Detector, consumer: &Detector, ideal: &Detector) -> String {
    let (w, h) = (880.0, 380.0);
    // Log10 of electrons, from 1 to 10^6.
    let (lo, hi) = (0.0f64, 6.0);

    let snr = |d: &Detector, s: f64| s / (s + d.read_noise * d.read_noise).sqrt();
    let curve = |d: &Detector| -> Vec<(f64, f64)> {
        (0..=300)
            .map(|i| {
                let l = lo + (hi - lo) * i as f64 / 300.0;
                (l, snr(d, 10f64.powf(l)).max(1e-3).log10())
            })
            .collect()
    };

    let mut left = Plot::new(w, h, (lo, hi), (-1.0, 3.0)).viewport(64.0, 54.0, 340.0, 262.0);
    left.axes(
        &[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        &[-1.0, 0.0, 1.0, 2.0, 3.0],
        |v| format!("10^{v:.0}"),
        |v| format!("10^{v:.0}"),
    );
    left.polyline(curve(ideal), &rgb(150, 150, 150), 1.6);
    left.polyline(curve(scientific), &rgb(48, 108, 186), 2.2);
    left.polyline(curve(consumer), &rgb(158, 40, 32), 2.2);

    // The crossover for each: below it the curve is twice as steep as above.
    for (d, colour) in [
        (scientific, rgb(48, 108, 186)),
        (consumer, rgb(158, 40, 32)),
    ] {
        let x = d.shot_noise_crossover().log10();
        let y = snr(d, d.shot_noise_crossover()).log10();
        left.polyline([(x, -1.0), (x, y)], "#00000030", 1.0);
        left.text(
            x + 0.1,
            -0.82,
            &format!("R2 = {:.0} e-", d.shot_noise_crossover()),
            11.0,
            &colour,
            "start",
        );
    }
    left.text(4.2, 1.55, "ideal counter", 11.5, "#8a8a8a", "start");
    left.text(3.4, 2.05, "1.4 e- read", 11.5, &rgb(48, 108, 186), "start");
    left.text(4.3, 0.85, "6 e- read", 11.5, &rgb(158, 40, 32), "start");
    left.title("Signal to noise, and where each sensor changes regime");
    left.caption("both axes log10");
    left.footnote("SNR against signal electrons");

    // --- Right: how the three noise terms grow with exposure, for the consumer sensor.
    // Read noise is flat, dark grows as the root of time, shot as the root of time too but
    // from a different constant. Which one dominates is the whole design question.
    let seconds = 60.0;
    let rate = 40.0;
    let dark_rate = consumer.dark_current.to_si();
    let noise_at = |t: f64| {
        let (s, d) = (rate * t, dark_rate * t);
        (s + d + consumer.read_noise * consumer.read_noise).sqrt()
    };
    let top = noise_at(seconds) * 1.08;
    let mut right = Plot::new(w, h, (0.0, seconds), (0.0, top)).viewport(470.0, 54.0, 350.0, 262.0);
    right.axes(
        &ticks(0.0, seconds, 6),
        &ticks(0.0, top, 5),
        |v| format!("{v:.0}"),
        |v| format!("{v:.0}"),
    );
    let sampled = |f: &dyn Fn(f64) -> f64| -> Vec<(f64, f64)> {
        (0..=240)
            .map(|i| {
                let t = seconds * i as f64 / 240.0;
                (t, f(t))
            })
            .collect()
    };
    right.polyline(sampled(&|_t| consumer.read_noise), &rgb(120, 120, 120), 1.8);
    right.polyline(sampled(&|t| (dark_rate * t).sqrt()), &rgb(158, 40, 32), 1.8);
    right.polyline(sampled(&|t| (rate * t).sqrt()), &rgb(48, 108, 186), 1.8);
    right.polyline(sampled(&noise_at), &rgb(24, 24, 24), 2.4);
    right.text(
        seconds * 0.06,
        consumer.read_noise + 1.4,
        "read",
        11.5,
        "#787878",
        "start",
    );
    right.text(
        seconds * 0.62,
        (dark_rate * seconds * 0.62).sqrt() + 1.4,
        "dark",
        11.5,
        &rgb(158, 40, 32),
        "start",
    );
    right.text(
        seconds * 0.34,
        (rate * seconds * 0.34).sqrt() + 1.4,
        "shot",
        11.5,
        &rgb(48, 108, 186),
        "start",
    );
    right.text(
        seconds * 0.72,
        noise_at(seconds * 0.72) + 1.6,
        "total",
        11.5,
        "#1a1a1a",
        "start",
    );
    right.footnote("noise electrons against exposure seconds, 6 e- sensor on a 40 e-/s source");

    document(w, h, [left.into_body(), right.into_body()])
}
