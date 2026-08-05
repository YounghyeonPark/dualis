//! Wavelength-dependent quantities: reflectance, transmittance, quantum
//! efficiency, lamp output, a fluorophore's absorption band.
//!
//! Optics is spectral. A single number for "how reflective is this" is a
//! convenience, not a description: a dichroic reflects 99% at 488 nm and 2% at
//! 509 nm, and the whole point of it lives in that difference. Everything in
//! this crate that used to be one number is a [`Spectrum`].
//!
//! Values are read with [`Spectrum::at`] and clamped, never extrapolated — a
//! curve measured from 400 to 700 nm says nothing about 1200 nm, so it holds its
//! endpoint rather than inventing a trend.

use serde::{Deserialize, Serialize};

/// Wavelength range worth evaluating, nm — near-UV through near-IR, which is
/// where silicon detectors and ordinary glass both work.
pub const VISIBLE_RANGE_NM: (f64, f64) = (350.0, 1100.0);

/// A quantity that varies with wavelength.
///
/// Four shapes, because four is what real data comes as: a flat number, a
/// measured curve, a filter's pass bands, and a Gaussian band such as an LED's
/// output or a fluorophore's absorption.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Spectrum {
    /// The same value at every wavelength. Honest for a metal mirror over a
    /// narrow band, and a useful way to isolate one effect from the others.
    Constant { value: f64 },
    /// A measured curve: `(nm, value)` pairs, linearly interpolated between
    /// them and held flat outside them.
    Curve { samples: Vec<(f64, f64)> },
    /// Pass bands with real edges. `in_band` inside, `out_of_band` outside, and
    /// a transition `edge_nm` wide between them.
    ///
    /// A perfect brick wall does not exist: an interference filter rolls off
    /// over a few nanometres, and it leaks — good blocking is OD 6, which is
    /// `out_of_band: 1e-6`, not zero. Both of those matter, because leaked
    /// excitation is the noise floor of a fluorescence image.
    Bands {
        /// `[low, high]` pairs, nm. A long-pass is one open-ended band.
        bands: Vec<[f64; 2]>,
        in_band: f64,
        #[serde(default)]
        out_of_band: f64,
        /// Width of the transition, nm. 0 is a brick wall.
        #[serde(default)]
        edge_nm: f64,
    },
    /// A Gaussian band: an LED's emission, a laser line, a dye's absorption
    /// peak. `fwhm_nm` is the full width at half maximum, so a 20 nm LED is
    /// `fwhm_nm: 20`.
    Gaussian {
        center_nm: f64,
        fwhm_nm: f64,
        peak: f64,
        /// Value far from the peak. Dyes have a real tail; lasers do not.
        #[serde(default)]
        floor: f64,
    },
    /// Planck's law: what a hot body actually emits.
    ///
    /// This is the shape of a tungsten lamp (about 3200 K), of daylight
    /// (5500 K) and of a star. It is not a convenience curve — it is why a
    /// halogen lamp is orange and starves the blue end of a colour camera,
    /// and why the same camera behaves differently under an LED. Scaled so its
    /// peak is `peak`, since a ray tracer only needs relative weights.
    Blackbody {
        temperature_k: f64,
        #[serde(default = "one")]
        peak: f64,
    },
}

fn one() -> f64 {
    1.0
}

/// Wien's displacement constant, in nm·K — where a blackbody's spectral radiance
/// peaks. Divide by the temperature to get nanometres: 3200 K peaks at 906 nm,
/// deep in the infrared, which is exactly why a tungsten lamp spends most of its
/// power on heat rather than light.
pub const WIEN_NM_K: f64 = 2_897_771.9;

impl Spectrum {
    /// A flat spectrum.
    pub fn constant(value: f64) -> Spectrum {
        Spectrum::Constant { value }
    }

    /// A measured curve. Samples are sorted, so they may be given in any order.
    pub fn curve(samples: Vec<(f64, f64)>) -> Spectrum {
        let mut samples = samples;
        samples.sort_by(|a, b| a.0.total_cmp(&b.0));
        Spectrum::Curve { samples }
    }

    /// A brick-wall multi-band filter: the simple case, kept simple.
    pub fn bands(bands: Vec<[f64; 2]>, in_band: f64, out_of_band: f64) -> Spectrum {
        Spectrum::Bands {
            bands,
            in_band,
            out_of_band,
            edge_nm: 0.0,
        }
    }

    /// A hot body's emission: tungsten at 3200 K, daylight at 5500 K, a star at
    /// whatever its surface temperature is.
    pub fn blackbody(temperature_k: f64) -> Spectrum {
        Spectrum::Blackbody {
            temperature_k,
            peak: 1.0,
        }
    }

    /// A realistic interference filter: `edge_nm` of roll-off and OD-6 blocking.
    pub fn interference_bands(bands: Vec<[f64; 2]>, transmission: f64, edge_nm: f64) -> Spectrum {
        Spectrum::Bands {
            bands,
            in_band: transmission,
            out_of_band: 1e-6,
            edge_nm,
        }
    }

    /// The value at a wavelength, in nanometres.
    pub fn at(&self, wavelength_nm: f64) -> f64 {
        match self {
            Spectrum::Constant { value } => *value,
            Spectrum::Curve { samples } => interpolate(samples, wavelength_nm),
            Spectrum::Bands {
                bands,
                in_band,
                out_of_band,
                edge_nm,
            } => {
                // Combine bands by how far *inside* one the wavelength is, not
                // by which gives the larger value. Bands may overlap and the
                // light does not care which one caught it — and `in_band` is
                // allowed to be the smaller of the two, which is exactly the
                // case for a dichroic's reflectance.
                let inside = bands
                    .iter()
                    .map(|[lo, hi]| band_weight(*lo, *hi, *edge_nm, wavelength_nm))
                    .fold(0.0f64, f64::max);
                out_of_band + (in_band - out_of_band) * inside
            }
            Spectrum::Gaussian {
                center_nm,
                fwhm_nm,
                peak,
                floor,
            } => {
                if *fwhm_nm <= 0.0 {
                    return if (wavelength_nm - center_nm).abs() < 1e-9 {
                        *peak
                    } else {
                        *floor
                    };
                }
                // exp(-4 ln2 (d/fwhm)^2) is exactly a half at d = fwhm/2.
                const FOUR_LN2: f64 = 2.772_588_722_239_781;
                let d = (wavelength_nm - center_nm) / fwhm_nm;
                floor + (peak - floor) * (-FOUR_LN2 * d * d).exp()
            }
            Spectrum::Blackbody {
                temperature_k,
                peak,
            } => {
                if *temperature_k <= 0.0 || wavelength_nm <= 0.0 {
                    return 0.0;
                }
                // Normalised to the peak, so the constants in front of Planck's
                // law cancel and only the shape survives.
                let at = |w: f64| planck_shape(w, *temperature_k);
                let peak_w = WIEN_NM_K / temperature_k;
                let denom = at(peak_w);
                if denom <= 0.0 {
                    return 0.0;
                }
                peak * (at(wavelength_nm) / denom)
            }
        }
    }

    /// Largest value anywhere in `range`, found by sampling. Used to check that
    /// reflectance and transmittance never sum past one.
    pub fn max_over(&self, range: (f64, f64), steps: usize) -> f64 {
        let (lo, hi) = range;
        let steps = steps.max(2);
        (0..=steps)
            .map(|i| self.at(lo + (hi - lo) * i as f64 / steps as f64))
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Wavelength of the peak and its value, sampled over `range`.
    pub fn peak_over(&self, range: (f64, f64), steps: usize) -> (f64, f64) {
        let (lo, hi) = range;
        let steps = steps.max(2);
        (0..=steps)
            .map(|i| {
                let w = lo + (hi - lo) * i as f64 / steps as f64;
                (w, self.at(w))
            })
            .fold(
                (lo, f64::NEG_INFINITY),
                |best, s| {
                    if s.1 > best.1 {
                        s
                    } else {
                        best
                    }
                },
            )
    }

    /// Evenly spaced `(nm, value)` samples across `range`, for plotting.
    pub fn sample(&self, range: (f64, f64), steps: usize) -> Vec<(f64, f64)> {
        let (lo, hi) = range;
        let steps = steps.max(2);
        (0..=steps)
            .map(|i| {
                let w = lo + (hi - lo) * i as f64 / steps as f64;
                (w, self.at(w))
            })
            .collect()
    }
}

/// Planck's law up to a constant: `1 / (l^5 (exp(c2 / (l T)) - 1))`, with the
/// wavelength in metres. Only the shape matters here, so the leading constants
/// are dropped and the result is normalised at the call site.
fn planck_shape(wavelength_nm: f64, temperature_k: f64) -> f64 {
    /// Second radiation constant hc/k, in metre-kelvin.
    const C2: f64 = 1.438_776_877e-2;
    let l = wavelength_nm * 1e-9;
    let x = C2 / (l * temperature_k);
    // exp overflows for cold bodies in the ultraviolet; the answer there is zero.
    if x > 700.0 {
        return 0.0;
    }
    let denom = x.exp_m1();
    if denom <= 0.0 {
        return 0.0;
    }
    1.0 / (l.powi(5) * denom)
}

/// Piecewise-linear read of sorted `(nm, value)` samples, held flat outside the
/// sampled range so a curve never extrapolates into nonsense.
fn interpolate(samples: &[(f64, f64)], wavelength_nm: f64) -> f64 {
    match samples.len() {
        0 => 0.0,
        1 => samples[0].1,
        _ => {
            let first = samples[0];
            let last = samples[samples.len() - 1];
            if wavelength_nm <= first.0 {
                return first.1;
            }
            if wavelength_nm >= last.0 {
                return last.1;
            }
            let i = samples.partition_point(|(w, _)| *w <= wavelength_nm).max(1);
            let (w0, v0) = samples[i - 1];
            let (w1, v1) = samples[i];
            let t = if (w1 - w0).abs() < 1e-12 {
                0.0
            } else {
                (wavelength_nm - w0) / (w1 - w0)
            };
            v0 + t * (v1 - v0)
        }
    }
}

/// How far inside a band a wavelength sits, 0..1, ramping linearly across
/// `edge_nm` at each shoulder.
fn band_weight(lo: f64, hi: f64, edge_nm: f64, wavelength_nm: f64) -> f64 {
    if edge_nm <= 0.0 {
        return f64::from(wavelength_nm >= lo && wavelength_nm <= hi);
    }
    let rising = ((wavelength_nm - (lo - edge_nm / 2.0)) / edge_nm).clamp(0.0, 1.0);
    let falling = (((hi + edge_nm / 2.0) - wavelength_nm) / edge_nm).clamp(0.0, 1.0);
    rising.min(falling)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_constant_spectrum_is_flat() {
        let s = Spectrum::constant(0.42);
        for w in [200.0, 550.0, 2000.0] {
            assert_eq!(s.at(w), 0.42);
        }
    }

    /// A curve interpolates between its samples and holds its endpoints, rather
    /// than extrapolating a trend it has no evidence for.
    #[test]
    fn a_curve_interpolates_and_never_extrapolates() {
        let s = Spectrum::curve(vec![(500.0, 0.2), (600.0, 0.8), (400.0, 0.0)]);
        assert_eq!(s.at(400.0), 0.0);
        assert_eq!(s.at(450.0), 0.1);
        assert_eq!(s.at(550.0), 0.5);
        assert_eq!(s.at(600.0), 0.8);
        // Outside the measured range it holds, in both directions.
        assert_eq!(s.at(100.0), 0.0);
        assert_eq!(s.at(2000.0), 0.8);
    }

    /// A brick-wall band is exact at its edges; a real one rolls off across
    /// `edge_nm` and leaks its blocking floor outside.
    #[test]
    fn bands_have_edges_and_leak() {
        let brick = Spectrum::bands(vec![[500.0, 560.0]], 0.95, 0.0);
        assert_eq!(brick.at(499.9), 0.0);
        assert_eq!(brick.at(500.0), 0.95);
        assert_eq!(brick.at(560.0), 0.95);
        assert_eq!(brick.at(560.1), 0.0);

        let real = Spectrum::interference_bands(vec![[500.0, 560.0]], 0.95, 10.0);
        // Half way up the shoulder at the nominal edge.
        assert!((real.at(500.0) - (1e-6 + 0.95 / 2.0)).abs() < 1e-6);
        assert!((real.at(530.0) - 0.95).abs() < 1e-9, "flat in the middle");
        // Blocking is OD 6, not zero: a real filter leaks, and that leak is the
        // noise floor of a fluorescence image.
        assert_eq!(real.at(400.0), 1e-6);
        assert!(real.at(400.0) > 0.0);
    }

    /// Overlapping bands combine by how far inside one the wavelength is, since
    /// the light does not care which band caught it.
    #[test]
    fn overlapping_bands_take_the_best() {
        let s = Spectrum::bands(vec![[400.0, 500.0], [450.0, 600.0]], 0.9, 0.01);
        assert_eq!(s.at(475.0), 0.9);
        assert_eq!(s.at(550.0), 0.9);
        assert_eq!(s.at(650.0), 0.01);
    }

    /// `in_band` may be *lower* than `out_of_band`, which is how a dichroic's
    /// reflectance is written: high everywhere except where it transmits.
    /// Combining by value rather than by insideness would break this.
    #[test]
    fn a_band_may_invert() {
        let reflect = Spectrum::Bands {
            bands: vec![[495.0, 545.0]],
            in_band: 0.05,
            out_of_band: 0.98,
            edge_nm: 10.0,
        };
        assert!(
            (reflect.at(520.0) - 0.05).abs() < 1e-9,
            "got {}",
            reflect.at(520.0)
        );
        assert!((reflect.at(405.0) - 0.98).abs() < 1e-9);
        assert!((reflect.at(610.0) - 0.98).abs() < 1e-9);
        // Halfway across the shoulder, halfway between the two levels.
        assert!((reflect.at(495.0) - (0.98 + (0.05 - 0.98) / 2.0)).abs() < 1e-9);
    }

    /// Planck's law, checked against Wien's displacement law: a blackbody peaks
    /// at 2898/T micrometres, and nowhere else. That closed form is independent
    /// of the implementation, which is what makes it worth testing against.
    #[test]
    fn a_blackbody_peaks_where_wien_says() {
        for temperature_k in [2000.0, 3200.0, 5500.0, 9000.0] {
            let s = Spectrum::blackbody(temperature_k);
            let expected = WIEN_NM_K / temperature_k;
            // Search a wide range so the peak is genuinely found, not assumed.
            let (peak_nm, peak) = s.peak_over((100.0, 4000.0), 4000);
            assert!(
                (peak_nm - expected).abs() < 3.0,
                "{temperature_k} K should peak at {expected:.0} nm, found {peak_nm:.0} nm"
            );
            // Normalised, so the value *at* the peak is exactly one; the
            // sampled maximum only lands within a grid step of it.
            assert!((s.at(expected) - 1.0).abs() < 1e-12);
            assert!(
                (peak - 1.0).abs() < 1e-4,
                "peak {peak} at {temperature_k} K"
            );
        }
    }

    /// The consequences a lamp choice actually has: tungsten is starved in the
    /// blue and dumps most of its output into the infrared, while daylight is
    /// nearly flat across the visible. This is why a halogen lamp looks orange.
    #[test]
    fn tungsten_is_red_and_daylight_is_not() {
        let tungsten = Spectrum::blackbody(3200.0);
        let daylight = Spectrum::blackbody(5500.0);
        let ratio = |s: &Spectrum| s.at(450.0) / s.at(650.0);
        assert!(
            ratio(&tungsten) < 0.45,
            "tungsten should be far weaker in the blue: ratio {}",
            ratio(&tungsten)
        );
        assert!(
            ratio(&daylight) > 0.9,
            "daylight should be nearly even across the visible: ratio {}",
            ratio(&daylight)
        );
        // And tungsten's peak is out in the infrared, which is where its power
        // goes: 906 nm, not anywhere you can see.
        assert!(tungsten.peak_over((300.0, 3000.0), 3000).0 > 850.0);
        // Nothing is emitted at zero or negative wavelength.
        assert_eq!(tungsten.at(0.0), 0.0);
        assert_eq!(Spectrum::blackbody(-5.0).at(550.0), 0.0);
        // A room-temperature body in the ultraviolet is not *zero* — Planck's law
        // has no cutoff — but it is unimaginably small, and the arithmetic must
        // reach that rather than overflowing on the way.
        let cold = Spectrum::blackbody(300.0).at(200.0);
        assert!(cold.is_finite() && cold < 1e-80, "got {cold}");
        // Far enough into the ultraviolet the exponent does overflow, and the
        // answer there is zero rather than an infinity.
        assert_eq!(Spectrum::blackbody(300.0).at(50.0), 0.0);
    }

    /// A Gaussian is exactly half its peak at half its full width, which is what
    /// "full width at half maximum" means.
    #[test]
    fn a_gaussian_is_half_at_its_half_width() {
        let s = Spectrum::Gaussian {
            center_nm: 488.0,
            fwhm_nm: 20.0,
            peak: 1.0,
            floor: 0.0,
        };
        assert!((s.at(488.0) - 1.0).abs() < 1e-12);
        assert!((s.at(498.0) - 0.5).abs() < 1e-9, "got {}", s.at(498.0));
        assert!((s.at(478.0) - 0.5).abs() < 1e-9);
        assert!(s.at(588.0) < 1e-9);
        let (peak_nm, peak) = s.peak_over(VISIBLE_RANGE_NM, 750);
        assert!((peak_nm - 488.0).abs() < 2.0 && (peak - 1.0).abs() < 0.01);
    }
}
