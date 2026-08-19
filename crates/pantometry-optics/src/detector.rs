//! Turning light into a number, and the noise that comes with it.
//!
//! Everything upstream of here is deterministic: a lamp has a spectrum, a filter passes a
//! fraction of it, a lens forms a spot. A detector is where that stops. Photons arrive at
//! random, and no amount of care in the optics changes the fact that counting `N` of them
//! carries an uncertainty of `√N`.
//!
//! This is the first module in the workspace whose correctness is a statement about a
//! *distribution* rather than a value. Everything so far has compared a computed number to
//! an analytic one — Planck against Wien, an FFT against a Bessel series. Here the claim
//! is that a sampled variance equals a sampled mean, and the test tolerance is set by
//! `1/√N` of the sample count rather than by an integrator's order. That is a different
//! testing regime, and it arrives with a shape that has closed forms all the way through.
//!
//! # The four noises, and why their order matters
//!
//! ```text
//! SNR = S / √(S + D + R²)
//! ```
//!
//! `S` is the signal in electrons, `D` the dark current's contribution, `R` the read noise.
//! The shape of that expression is the whole of low-light imaging:
//!
//! - **Read-limited**, when `S ≪ R²`: the noise is a constant, so `SNR ∝ S` and
//!   doubling the exposure doubles the quality.
//! - **Shot-limited**, when `S ≫ R²`: the noise is `√S`, so `SNR ∝ √S` and doubling the
//!   exposure buys only 41%.
//!
//! The crossover sits at `S = R²`, which for a 5-electron read noise is 25 electrons —
//! a startlingly small number, and the reason a scientific camera's read noise is quoted
//! before anything else. Past that point the camera has stopped mattering and the photons
//! are in charge.
//!
//! # One element, deliberately
//!
//! A single detector: a photodiode, a photomultiplier, a power meter, an integrating
//! sphere. Not a pixel array, and the reason is structural rather than effort. A sensor
//! needs to say *where* each photon landed, and the kernel's coupling carries one number
//! per channel per step — so a pixel array waits on a spatially resolved flux, and that is
//! a larger decision than a detector model.

use pantometry_core::Rng;
use pantometry_units::{Energy, Frequency, Length, Power, Time};

use crate::radiometry::SpectralPower;
use crate::spectrum::Spectrum;

/// The expected electron counts for one source and one exposure time, worked out once.
///
/// Separating this from [`Detector::sample`] is not premature tidiness. Getting the
/// expected signal means integrating the source against the efficiency curve, which is a
/// few thousand spectrum evaluations — and the reason to have a noise model at all is to
/// take many frames of the *same* scene. Recomputing the integral per frame made a
/// hundred thousand exposures take a minute; hoisting it makes them instant, and the
/// separation says plainly which part depends on the light and which on the dice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Exposure {
    /// Expected electrons from the light.
    pub signal: f64,
    /// Expected electrons from the dark current.
    pub dark: f64,
}

impl Exposure {
    /// Signal-to-noise ratio for this exposure and a given read noise: `S/√(S + D + R²)`.
    pub fn snr(&self, read_noise: f64) -> f64 {
        let noise = (self.signal + self.dark + read_noise * read_noise).sqrt();
        if noise <= 0.0 {
            return 0.0;
        }
        self.signal / noise
    }
}

/// What came out of one exposure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Reading {
    /// Electrons in the well, including every noise source.
    pub electrons: f64,
    /// The digitised output.
    pub counts: u32,
    /// Whether the well filled, which makes the reading a lower bound rather than a
    /// measurement.
    pub saturated: bool,
}

/// A single-element photon detector.
#[derive(Clone, Debug, PartialEq)]
pub struct Detector {
    /// Quantum efficiency against wavelength: the fraction of arriving photons that become
    /// countable electrons.
    ///
    /// Spectral, and that matters more than it looks. A silicon detector is near 90% in the
    /// red and under 20% in the blue, so two sources of equal *power* produce very
    /// different signals — and a detector modelled with one number cannot express the
    /// difference.
    ///
    /// Give the curve an explicit zero wherever the physics has one. A
    /// [`Spectrum::Curve`] holds its endpoint rather than extrapolating, which is right for
    /// a measurement that simply stops — but a silicon sensor does not merely stop being
    /// measured past 1100 nm, it stops working, because a photon below the bandgap cannot
    /// lift an electron at all. A curve left hanging at 5% will happily detect infrared the
    /// device is blind to.
    pub quantum_efficiency: Spectrum,
    /// Read noise, in electrons RMS. Added once per exposure however long it is, which is
    /// what makes it dominate at low light and vanish at high.
    pub read_noise: f64,
    /// Dark current: electrons generated with no light at all, per second.
    ///
    /// Thermally driven, so it roughly halves for every 6 K of cooling — which is the whole
    /// reason a long-exposure camera is cooled, and why cooling buys nothing for a short
    /// one.
    pub dark_current: Frequency,
    /// Full well, in electrons. Past this the response stops being a measurement.
    pub well_depth: f64,
    /// Electrons per digital unit.
    pub gain: f64,
    /// Digitiser depth in bits.
    pub bits: u32,
}

impl Detector {
    /// A cooled scientific CMOS: quiet, deep, and near-perfect in the visible.
    pub fn scientific_cmos() -> Detector {
        Detector {
            quantum_efficiency: Spectrum::curve(vec![
                (350.0, 0.30),
                (400.0, 0.60),
                (500.0, 0.82),
                (600.0, 0.80),
                (700.0, 0.65),
                (850.0, 0.35),
                (1000.0, 0.08),
                // Silicon's bandgap: nothing past here, and the zero has to be written
                // down or the curve would hold 8% forever.
                (1100.0, 0.0),
            ]),
            read_noise: 1.4,
            dark_current: Frequency::hz(0.5),
            well_depth: 30_000.0,
            gain: 0.5,
            bits: 16,
        }
    }

    /// An uncooled consumer sensor: noisier, shallower, and warm.
    pub fn consumer_cmos() -> Detector {
        Detector {
            quantum_efficiency: Spectrum::curve(vec![
                (400.0, 0.35),
                (550.0, 0.55),
                (650.0, 0.50),
                (800.0, 0.25),
                (950.0, 0.05),
                (1100.0, 0.0),
            ]),
            read_noise: 6.0,
            dark_current: Frequency::hz(20.0),
            well_depth: 12_000.0,
            gain: 2.0,
            bits: 12,
        }
    }

    /// An idealised counter: every photon counted, nothing added.
    ///
    /// Useful as the thing a real detector is compared against — its signal-to-noise is
    /// exactly `√N`, which is the ceiling counting allows.
    pub fn ideal() -> Detector {
        Detector {
            quantum_efficiency: Spectrum::constant(1.0),
            read_noise: 0.0,
            dark_current: Frequency::from_si(0.0),
            well_depth: f64::INFINITY,
            gain: 1.0,
            bits: 32,
        }
    }

    /// Largest digital value this detector can report.
    pub fn full_scale(&self) -> u32 {
        if self.bits >= 32 {
            u32::MAX
        } else {
            (1u32 << self.bits) - 1
        }
    }

    /// Expected signal electrons from a source over an exposure.
    ///
    /// The photon rate through the quantum efficiency curve, times the time. Not the power
    /// times an efficiency: a detector counts photons, and the conversion depends on where
    /// in the spectrum they were — see
    /// [`SpectralPower::photon_rate_through`](crate::radiometry::SpectralPower::photon_rate_through).
    pub fn signal(&self, source: &SpectralPower, exposure: Time) -> f64 {
        let rate = source.photon_rate_through(&self.quantum_efficiency);
        (rate * exposure).to_si().max(0.0)
    }

    /// Expected dark electrons over an exposure.
    pub fn dark(&self, exposure: Time) -> f64 {
        (self.dark_current * exposure).to_si().max(0.0)
    }

    /// The expected counts for a source and an exposure, computed once so that many frames
    /// of the same scene do not each re-integrate the spectrum.
    pub fn prepare(&self, source: &SpectralPower, exposure: Time) -> Exposure {
        Exposure {
            signal: self.signal(source, exposure),
            dark: self.dark(exposure),
        }
    }

    /// Signal-to-noise ratio, `S/√(S + D + R²)`.
    ///
    /// The closed form the whole module is checked against, and the expression a camera is
    /// really chosen by.
    pub fn snr(&self, source: &SpectralPower, exposure: Time) -> f64 {
        self.prepare(source, exposure).snr(self.read_noise)
    }

    /// Signal at which shot noise overtakes read noise: `R²`.
    ///
    /// Below it, exposure time is worth its weight linearly; above it, only as a square
    /// root. For a 1.4-electron read noise this is two electrons — which is why a good
    /// scientific camera is shot-limited almost everywhere and its specification sheet
    /// stops being interesting.
    pub fn shot_noise_crossover(&self) -> f64 {
        self.read_noise * self.read_noise
    }

    /// Exposure needed to reach a target signal-to-noise ratio, by solving the CCD
    /// equation.
    ///
    /// Worth having in closed form rather than by search, because the answer is the thing
    /// an experiment is planned around, and because it makes the two regimes visible: the
    /// solution is quadratic in the target, so a factor of two in the SNR costs between
    /// two and four times the exposure depending on which regime you are in.
    pub fn exposure_for_snr(&self, source: &SpectralPower, target: f64) -> Time {
        let rate = source.photon_rate_through(&self.quantum_efficiency).to_si();
        let dark = self.dark_current.to_si();
        if rate <= 0.0 || target <= 0.0 {
            return Time::from_si(f64::INFINITY);
        }
        // S = rate*t, D = dark*t. Solving target^2 (S + D + R^2) = S^2 for t:
        //   rate^2 t^2 - target^2 (rate + dark) t - target^2 R^2 = 0
        let a = rate * rate;
        let b = -target * target * (rate + dark);
        let c = -target * target * self.read_noise * self.read_noise;
        let discriminant = b * b - 4.0 * a * c;
        if discriminant < 0.0 {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si((-b + discriminant.sqrt()) / (2.0 * a))
    }

    /// One exposure, with every noise source drawn.
    ///
    /// Shot noise on the signal and on the dark current separately, because they are
    /// independent arrivals and their variances add; read noise as a Gaussian, because it
    /// is the sum of many small electronic contributions rather than a count of anything.
    pub fn measure(&self, source: &SpectralPower, exposure: Time, rng: &mut Rng) -> Reading {
        self.sample(&self.prepare(source, exposure), rng)
    }

    /// One exposure from counts already worked out, which is what a run of frames should
    /// use — see [`Exposure`].
    pub fn sample(&self, exposure: &Exposure, rng: &mut Rng) -> Reading {
        let signal = rng.poisson(exposure.signal) as f64;
        let dark = rng.poisson(exposure.dark) as f64;
        let read = rng.normal(0.0, self.read_noise);

        let collected = signal + dark + read;
        let saturated = collected >= self.well_depth;
        let electrons = collected.min(self.well_depth);

        // Digitising loses everything below one gain step, which contributes a variance of
        // `g²/12` — the standard result for rounding to a grid, and a real noise source
        // that a coarse gain can make dominant.
        let counts = if self.gain > 0.0 {
            (electrons / self.gain)
                .round()
                .max(0.0)
                .min(self.full_scale() as f64) as u32
        } else {
            0
        };
        Reading {
            electrons,
            counts,
            saturated,
        }
    }

    /// Variance the digitiser adds by rounding, `g²/12`.
    pub fn digitisation_variance(&self) -> f64 {
        self.gain * self.gain / 12.0
    }

    /// Energy one electron's worth of signal represents, at a wavelength.
    ///
    /// The bridge back to radiometry: a reading is a count, and turning it into watts needs
    /// to know both the quantum efficiency and the photon energy there.
    pub fn energy_per_electron(&self, wavelength: Length) -> Energy {
        let qe = self
            .quantum_efficiency
            .at(wavelength)
            .max(f64::MIN_POSITIVE);
        Energy::from_si(pantometry_units::photon_energy(wavelength).to_si() / qe)
    }

    /// Smallest power this detector can see at a given signal-to-noise ratio and exposure,
    /// assuming the source is monochromatic.
    ///
    /// The number a detector is actually specified by, and it falls as the square root of
    /// the exposure only once the shot noise is in charge.
    pub fn noise_equivalent_power(
        &self,
        wavelength: Length,
        exposure: Time,
        target_snr: f64,
    ) -> Power {
        let qe = self.quantum_efficiency.at(wavelength);
        if qe <= 0.0 || exposure.to_si() <= 0.0 {
            return Power::from_si(f64::INFINITY);
        }
        // Solve target^2 (S + D + R^2) = S^2 for S, then convert electrons to watts.
        let d = self.dark(exposure);
        let r2 = self.read_noise * self.read_noise;
        let t2 = target_snr * target_snr;
        let electrons = 0.5 * (t2 + (t2 * t2 + 4.0 * t2 * (d + r2)).sqrt());
        let photons = electrons / qe;
        let energy = pantometry_units::photon_energy(wavelength).to_si();
        Power::from_si(photons * energy / exposure.to_si())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spectrum::VISIBLE_RANGE;

    fn nm(v: f64) -> Length {
        Length::nm(v)
    }

    /// A monochromatic source of a stated power, for putting a known number of photons on
    /// a detector.
    fn line(centre_nm: f64, power: Power) -> SpectralPower {
        SpectralPower::new(
            Spectrum::gaussian(nm(centre_nm), nm(1.0), 1.0),
            power,
            (nm(centre_nm - 15.0), nm(centre_nm + 15.0)),
        )
    }

    /// The signal is a photon count through the efficiency curve, not a power times a
    /// number — and the two differ, because a detector counts photons.
    #[test]
    fn the_signal_is_a_photon_count_through_the_efficiency_curve() {
        let detector = Detector::ideal();
        let source = line(550.0, Power::from_si(1e-15));
        let exposure = Time::ms(100.0);

        // A perfect detector's signal is exactly the photon count.
        let expected = source.photon_rate().to_si() * exposure.to_si();
        assert!(
            (detector.signal(&source, exposure) / expected - 1.0).abs() < 1e-9,
            "got {}, expected {expected}",
            detector.signal(&source, exposure)
        );
        // A femtowatt of green for a tenth of a second is about 280 photons.
        assert!(
            (expected - 277.0).abs() < 5.0,
            "a femtowatt of green is {expected:.0} photons in 100 ms"
        );

        // And a real detector's efficiency is spectral: the same power in the blue and in
        // the green gives different counts.
        let cmos = Detector::scientific_cmos();
        let blue = cmos.signal(&line(400.0, Power::from_si(1e-15)), exposure);
        let green = cmos.signal(&line(550.0, Power::from_si(1e-15)), exposure);
        assert!(
            green > blue * 1.15,
            "the same power should count differently by colour: {green:.1} against \
             {blue:.1}"
        );
    }

    /// **The property this module exists to get right.** Photon arrivals are Poisson, so
    /// the variance of a measurement equals its mean — and the noise on a count is its
    /// square root whatever the optics did beforehand.
    ///
    /// The first test in this workspace whose tolerance comes from `1/√N` of a sample count
    /// rather than from a closed form's accuracy.
    #[test]
    fn shot_noise_has_variance_equal_to_its_mean() {
        let detector = Detector::ideal();
        let exposure = Time::ms(10.0);
        const N: usize = 100_000;

        for power in [1e-16f64, 1e-15, 1e-14] {
            let source = line(550.0, Power::from_si(power));
            let prepared = detector.prepare(&source, exposure);
            let mean = prepared.signal;
            assert!(mean > 1.0, "the test needs a countable signal");

            let readings: Vec<f64> = (0..N)
                .map(|i| {
                    let mut rng = Rng::for_index(0xDE7EC7, i as u64);
                    detector.sample(&prepared, &mut rng).electrons
                })
                .collect();
            let measured_mean = readings.iter().sum::<f64>() / N as f64;
            let variance = readings
                .iter()
                .map(|r| (r - measured_mean).powi(2))
                .sum::<f64>()
                / N as f64;

            assert!(
                (measured_mean / mean - 1.0).abs() < 0.02,
                "mean should be the expected signal: {measured_mean:.2} against {mean:.2}"
            );
            assert!(
                (variance / measured_mean - 1.0).abs() < 0.05,
                "variance should equal the mean: {variance:.2} against {measured_mean:.2}"
            );
        }
    }

    /// The CCD equation, checked against sampling rather than restated.
    #[test]
    fn the_snr_equation_matches_what_the_noise_actually_does() {
        let exposure = Time::ms(50.0);
        const N: usize = 60_000;

        for detector in [Detector::scientific_cmos(), Detector::consumer_cmos()] {
            for power in [1e-16f64, 1e-14] {
                let source = line(550.0, Power::from_si(power));
                let prepared = detector.prepare(&source, exposure);
                let predicted = detector.snr(&source, exposure);

                let readings: Vec<f64> = (0..N)
                    .map(|i| {
                        let mut rng = Rng::for_index(99, i as u64);
                        detector.sample(&prepared, &mut rng).electrons
                    })
                    .collect();
                let mean = readings.iter().sum::<f64>() / N as f64;
                let sigma =
                    (readings.iter().map(|r| (r - mean).powi(2)).sum::<f64>() / N as f64).sqrt();
                // The measured ratio uses the signal alone over the total noise, which is
                // what the closed form means.
                let measured = prepared.signal / sigma;
                assert!(
                    (measured / predicted - 1.0).abs() < 0.05,
                    "read noise {}: predicted SNR {predicted:.3}, measured {measured:.3}",
                    detector.read_noise
                );
            }
        }
    }

    /// The two regimes, and the crossover between them at `S = R²`.
    ///
    /// This is the shape of low-light imaging: below the crossover an exposure twice as
    /// long is twice as good, and above it only 41% better. Both halves are asserted,
    /// because only the second is folklore.
    #[test]
    fn exposure_pays_linearly_below_the_crossover_and_as_a_root_above_it() {
        let detector = Detector::consumer_cmos();
        let crossover = detector.shot_noise_crossover();
        assert_eq!(crossover, 36.0, "6 electrons of read noise crosses at 36");

        // A source faint enough that even a long exposure stays read-limited.
        let faint = line(550.0, Power::from_si(3e-18));
        let short = Time::ms(20.0);
        assert!(
            detector.signal(&faint, short) < crossover / 5.0,
            "the faint case should be well read-limited, got {} electrons",
            detector.signal(&faint, short)
        );
        let ratio_read = detector.snr(&faint, short * 4.0) / detector.snr(&faint, short);
        assert!(
            (ratio_read - 4.0).abs() < 0.35,
            "read-limited: four times the exposure should be about four times the SNR, \
             got {ratio_read:.3}"
        );

        // And a bright one that is firmly shot-limited. 3e-13 W of green through a 55%
        // efficiency for 20 ms is 9100 electrons, 250 times the crossover — a factor of ten
        // less light would only have been 25 times it, which is not "firmly".
        let bright = line(550.0, Power::from_si(3e-13));
        assert!(
            detector.signal(&bright, short) > crossover * 100.0,
            "the bright case should be well shot-limited, got {} electrons against a \
             crossover of {crossover}",
            detector.signal(&bright, short)
        );
        let ratio_shot = detector.snr(&bright, short * 4.0) / detector.snr(&bright, short);
        assert!(
            (ratio_shot - 2.0).abs() < 0.05,
            "shot-limited: four times the exposure buys only twice the SNR, got \
             {ratio_shot:.3}"
        );
    }

    /// A perfect counter's signal-to-noise is exactly the square root of the count, which
    /// is the ceiling counting allows and the thing every real detector is compared to.
    #[test]
    fn an_ideal_detector_is_exactly_root_n() {
        let detector = Detector::ideal();
        let exposure = Time::ms(100.0);
        for power in [1e-16f64, 1e-15, 1e-13] {
            let source = line(550.0, Power::from_si(power));
            let n = detector.signal(&source, exposure);
            assert!(
                (detector.snr(&source, exposure) / n.sqrt() - 1.0).abs() < 1e-12,
                "SNR should be root N exactly"
            );
        }
        assert_eq!(detector.shot_noise_crossover(), 0.0, "nothing to cross");
    }

    /// Dark current grows with the exposure and its noise grows with the root of it, so a
    /// long exposure is eventually dark-limited however bright the source is not.
    ///
    /// And it is why cooling matters for a long exposure and not a short one.
    #[test]
    fn dark_current_grows_with_time_and_its_noise_with_the_root() {
        let warm = Detector::consumer_cmos();
        let cold = Detector::scientific_cmos();
        for seconds in [1.0f64, 10.0, 60.0] {
            let t = Time::s(seconds);
            assert!((warm.dark(t) - 20.0 * seconds).abs() < 1e-9);
            assert!(
                cold.dark(t) < warm.dark(t) / 30.0,
                "a cooled sensor should be far darker"
            );
        }

        // Over a minute the warm sensor's dark shot noise alone is 35 electrons, six times
        // its read noise — so cooling has become the only thing that would help.
        let minute = Time::s(60.0);
        let dark_noise = warm.dark(minute).sqrt();
        assert!(
            (dark_noise - 34.6).abs() < 0.5,
            "dark shot noise over a minute is {dark_noise:.1} electrons"
        );
        assert!(dark_noise > warm.read_noise * 5.0);

        // Over a millisecond it is nothing, and the read noise is everything.
        let instant = Time::ms(1.0);
        assert!(warm.dark(instant).sqrt() < warm.read_noise / 20.0);
    }

    /// Saturation is reported rather than hidden, because a full well makes a reading a
    /// lower bound and a caller that treats it as a measurement gets a wrong answer with
    /// no symptom.
    #[test]
    fn a_full_well_is_reported_as_saturated() {
        let detector = Detector::consumer_cmos();
        let exposure = Time::ms(100.0);
        let blinding = line(550.0, Power::from_si(1e-11));
        assert!(detector.signal(&blinding, exposure) > detector.well_depth * 2.0);

        let mut rng = Rng::new(3);
        let reading = detector.measure(&blinding, exposure, &mut rng);
        assert!(reading.saturated);
        assert!((reading.electrons - detector.well_depth).abs() < 1e-9);
        // The digital value is clipped to full scale too, not wrapped.
        assert!(reading.counts <= detector.full_scale());

        // A modest source does not saturate and is not flagged.
        let modest = line(550.0, Power::from_si(1e-15));
        let reading = detector.measure(&modest, exposure, &mut rng);
        assert!(!reading.saturated);
        assert!(reading.electrons < detector.well_depth);
    }

    /// Digitising is its own noise source, with variance `g²/12` — the standard result for
    /// rounding to a grid. A coarse gain can make it dominant, which is why a cheap
    /// sensor's read noise figure can be misleading.
    #[test]
    fn digitising_adds_the_variance_of_rounding() {
        let detector = Detector::consumer_cmos();
        // Gain 2 electrons per count: variance 1/3, or 0.58 electrons RMS.
        assert!((detector.digitisation_variance() - 1.0 / 3.0).abs() < 1e-12);

        // Measured: read a constant electron count through the digitiser many times with
        // the analogue noise removed, and the spread is the rounding alone.
        let quiet = Detector {
            read_noise: 0.0,
            dark_current: Frequency::from_si(0.0),
            gain: 4.0,
            ..Detector::consumer_cmos()
        };
        assert!((quiet.digitisation_variance() - 16.0 / 12.0).abs() < 1e-12);

        const N: usize = 40_000;
        let residuals: Vec<f64> = (0..N)
            .map(|i| {
                let electrons = 1000.0 + i as f64 * 1e-3;
                let counts = (electrons / quiet.gain).round();
                counts * quiet.gain - electrons
            })
            .collect();
        let mean = residuals.iter().sum::<f64>() / N as f64;
        let variance = residuals.iter().map(|r| (r - mean).powi(2)).sum::<f64>() / N as f64;
        assert!(
            (variance / quiet.digitisation_variance() - 1.0).abs() < 0.05,
            "rounding variance {variance:.4} against g^2/12 = {:.4}",
            quiet.digitisation_variance()
        );
    }

    /// Solving the CCD equation for exposure, checked by putting the answer back in.
    #[test]
    fn the_exposure_for_a_target_snr_round_trips() {
        for detector in [Detector::scientific_cmos(), Detector::consumer_cmos()] {
            let source = line(550.0, Power::from_si(1e-16));
            for target in [3.0f64, 10.0, 50.0] {
                let t = detector.exposure_for_snr(&source, target);
                assert!(t.to_si() > 0.0 && t.to_si().is_finite());
                let achieved = detector.snr(&source, t);
                assert!(
                    (achieved / target - 1.0).abs() < 1e-9,
                    "asked for SNR {target}, the exposure it returned gives {achieved}"
                );
            }
        }
        // A dark detector never gets there, and says so rather than looping.
        let dark = Detector::ideal();
        let nothing =
            SpectralPower::new(Spectrum::constant(0.0), Power::from_si(0.0), VISIBLE_RANGE);
        assert!(!dark.exposure_for_snr(&nothing, 10.0).to_si().is_finite());
    }

    /// Noise-equivalent power, and the thing about it that is easy to get wrong: it falls
    /// as the square root of the exposure only once the shot noise is in charge. While the
    /// read noise dominates it falls *linearly*, which is a much better return.
    #[test]
    fn noise_equivalent_power_improves_faster_while_read_limited() {
        let detector = Detector::consumer_cmos();
        let green = nm(550.0);
        let nep = |ms: f64| {
            detector
                .noise_equivalent_power(green, Time::ms(ms), 5.0)
                .to_si()
        };

        // Short exposures are read-limited, so quadrupling the time nearly quarters the
        // detectable power.
        let read_limited = nep(1.0) / nep(4.0);
        assert!(
            read_limited > 3.0,
            "read-limited: four times the exposure should nearly quarter the NEP, got \
             {read_limited:.2}"
        );

        // Long ones are dark-limited, where the return has fallen to a square root. It has
        // to be *well* into that regime: at ten seconds the read noise and the target's own
        // floor still contribute, and the ratio is 2.29 rather than 2 — approaching the
        // limit rather than sitting at it.
        let dark_limited = nep(100_000.0) / nep(400_000.0);
        assert!(
            (dark_limited - 2.0).abs() < 0.15,
            "dark-limited: four times the exposure only halves it, got {dark_limited:.2}"
        );
        assert!(
            read_limited > dark_limited * 1.4,
            "the return really does fall off"
        );

        // A wavelength the detector cannot see needs infinite power.
        assert!(!detector
            .noise_equivalent_power(nm(1400.0), Time::ms(10.0), 5.0)
            .to_si()
            .is_finite());
    }

    /// The bridge back to radiometry: a count becomes watts only with both the quantum
    /// efficiency and the photon energy, and neither alone is enough.
    #[test]
    fn a_count_becomes_watts_through_the_efficiency_and_the_photon_energy() {
        let detector = Detector::scientific_cmos();
        let green = detector.energy_per_electron(nm(550.0));
        // A green photon is 3.6e-19 J and the efficiency there is about 0.81, so each
        // electron represents rather more than one photon's worth.
        let photon = pantometry_units::photon_energy(nm(550.0));
        assert!(green > photon, "an imperfect detector wastes photons");
        assert!(
            (green / photon - 1.0 / 0.81).abs() < 0.05,
            "the ratio should be one over the efficiency, got {}",
            green / photon
        );
        // In the deep red where the efficiency collapses, an electron is worth far more.
        assert!(detector.energy_per_electron(nm(1000.0)) > green * 5.0);
    }

    /// Measurements are reproducible per index, so a run is repeatable whatever order the
    /// exposures were computed in — the same guarantee the rest of the workspace makes.
    #[test]
    fn measurements_are_reproducible_per_index() {
        let detector = Detector::scientific_cmos();
        let source = line(550.0, Power::from_si(1e-15));
        let exposure = Time::ms(20.0);
        let draw = |i: u64| {
            let mut rng = Rng::for_index(1234, i);
            detector.measure(&source, exposure, &mut rng)
        };
        let forward: Vec<Reading> = (0..64).map(draw).collect();
        let backward: Vec<Reading> = (0..64).rev().map(draw).collect();
        assert_eq!(
            forward,
            backward.into_iter().rev().collect::<Vec<_>>(),
            "the order exposures were computed in must not matter"
        );
        // And two different indices really do differ, so this is not passing trivially.
        assert_ne!(forward[0], forward[1]);
    }

    /// A detector in the dark reads its own noise and nothing else, which is the bias frame
    /// every real camera has to have subtracted.
    #[test]
    fn a_dark_frame_is_the_detectors_own_noise() {
        let detector = Detector::scientific_cmos();
        let nothing =
            SpectralPower::new(Spectrum::constant(0.0), Power::from_si(0.0), VISIBLE_RANGE);
        let exposure = Time::s(1.0);
        assert_eq!(detector.signal(&nothing, exposure), 0.0);

        const N: usize = 40_000;
        let prepared = detector.prepare(&nothing, exposure);
        let frames: Vec<f64> = (0..N)
            .map(|i| {
                let mut rng = Rng::for_index(555, i as u64);
                detector.sample(&prepared, &mut rng).electrons
            })
            .collect();
        let mean = frames.iter().sum::<f64>() / N as f64;
        let variance = frames.iter().map(|f| (f - mean).powi(2)).sum::<f64>() / N as f64;

        // The mean is the dark current alone.
        assert!(
            (mean - detector.dark(exposure)).abs() < 0.1,
            "a dark frame's mean is its dark current: {mean:.3} against {}",
            detector.dark(exposure)
        );
        // And the variance is the dark shot noise plus the read noise squared, added
        // because they are independent.
        let expected = detector.dark(exposure) + detector.read_noise * detector.read_noise;
        assert!(
            (variance / expected - 1.0).abs() < 0.05,
            "variance {variance:.3} against {expected:.3}"
        );
    }
}
