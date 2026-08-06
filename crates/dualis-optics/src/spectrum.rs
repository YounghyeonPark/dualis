//! Wavelength-dependent quantities: reflectance, transmittance, quantum
//! efficiency, lamp output, a fluorophore's absorption band.
//!
//! Optics is spectral. A single number for "how reflective is this" is a
//! convenience, not a description: a dichroic reflects 99% at 488 nm and 2% at
//! 509 nm, and the whole point of it lives in that difference. Everything in this
//! crate that used to be one number is a [`Spectrum`].
//!
//! Values are read with [`Spectrum::at`] and clamped, never extrapolated — a curve
//! measured from 400 to 700 nm says nothing about 1200 nm, so it holds its
//! endpoint rather than inventing a trend.
//!
//! # Nanometres in the data, [`Length`] in the API
//!
//! The enum's fields are nanometres as plain numbers, because they are the
//! serialised form and `"center_nm": 488` is worth reading where `4.88e-7` is not.
//! The methods take [`Length`], because that is where a wavelength meets the rest
//! of the physics and where confusing it with a path length would matter. The
//! conversion happens once, at the constructor and at [`Spectrum::at`].

use std::ops::{Add, Mul};

use dualis_units::Length;
use serde::{Deserialize, Serialize};

/// Wavelength range worth evaluating — near-UV through near-IR, which is where
/// silicon detectors and ordinary glass both work.
pub const VISIBLE_RANGE: (Length, Length) = (Length::from_si(350e-9), Length::from_si(1100e-9));

/// The same range in nanometres, for the internal loops that work in nm.
pub(crate) const VISIBLE_RANGE_NM: (f64, f64) = (350.0, 1100.0);

/// A quantity that varies with wavelength.
///
/// Five shapes because five is what real data comes as — a flat number, a measured
/// curve, a filter's pass bands, a Gaussian band such as an LED's output or a
/// fluorophore's absorption, and a hot body's Planck spectrum — and two more because
/// data is not the only thing a spectrum comes from. A light path is a *product* of the
/// things along it, and a sample with two dyes in it is a *sum*, so
/// [`Product`](Spectrum::Product) and [`Sum`](Spectrum::Sum) let a whole chain be one
/// value rather than a loop the caller has to write.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Spectrum {
    /// The same value at every wavelength. Honest for a metal mirror over a narrow
    /// band, and a useful way to isolate one effect from the others.
    Constant {
        /// The value, at every wavelength.
        value: f64,
    },
    /// A measured curve: `(nm, value)` pairs, linearly interpolated between them
    /// and held flat outside them.
    Curve {
        /// `(nanometres, value)`, in increasing wavelength. Held flat past either end, which
        /// is right for a measurement that simply stops and wrong for a physical cutoff —
        /// see the note on `Detector::quantum_efficiency` for what that cost once.
        samples: Vec<(f64, f64)>,
    },
    /// Pass bands with real edges. `in_band` inside, `out_of_band` outside, and a
    /// transition `edge_nm` wide between them.
    ///
    /// A perfect brick wall does not exist: an interference filter rolls off over a
    /// few nanometres, and it leaks — good blocking is OD 6, which is
    /// `out_of_band: 1e-6`, not zero. Both of those matter, because leaked
    /// excitation is the noise floor of a fluorescence image.
    Bands {
        /// `[low, high]` pairs, nm. A long-pass is one open-ended band.
        bands: Vec<[f64; 2]>,
        /// Transmission inside a band.
        in_band: f64,
        /// Transmission outside every band. Not zero for a real filter: blocking is finite,
        /// and six orders of magnitude of it is what a fluorescence experiment lives on.
        #[serde(default)]
        out_of_band: f64,
        /// Width of the transition, nm. 0 is a brick wall.
        ///
        /// A zero-width edge is exact only up to the wavelength's round trip
        /// through metres, which is a femtometre — so whether a wavelength written
        /// as exactly the band edge falls inside or outside is not something to
        /// rely on. That is not much of a loss: no real filter has a brick wall,
        /// and a physically meaningful edge has a width.
        #[serde(default)]
        edge_nm: f64,
    },
    /// A Gaussian band: an LED's emission, a laser line, a dye's absorption peak.
    /// `fwhm_nm` is the full width at half maximum, so a 20 nm LED is
    /// `fwhm_nm: 20`.
    Gaussian {
        /// Centre wavelength, nm.
        center_nm: f64,
        /// Full width at half maximum, nm — not the standard deviation.
        fwhm_nm: f64,
        /// Value at the centre.
        peak: f64,
        /// Value far from the peak. Dyes have a real tail; lasers do not.
        #[serde(default)]
        floor: f64,
    },
    /// Planck's law: what a hot body actually emits.
    ///
    /// This is the shape of a tungsten lamp (about 3200 K), of daylight (5500 K)
    /// and of a star. It is not a convenience curve — it is why a halogen lamp is
    /// orange and starves the blue end of a colour camera, and why the same camera
    /// behaves differently under an LED. Scaled so its peak is `peak`, since a ray
    /// tracer only needs relative weights.
    Blackbody {
        /// Colour temperature in kelvin. A tungsten lamp is about 3200, the sun about 5800.
        temperature_k: f64,
        /// What the peak is scaled to. Planck's law in absolute units spans many decades,
        /// and a ray tracer only needs relative weights.
        #[serde(default = "one")]
        peak: f64,
    },
    /// Several spectra multiplied together, wavelength by wavelength.
    ///
    /// A light path is a product. Following excitation through a fluorescence microscope
    /// is the lamp, times the excitation filter, times the dichroic's transmission, times
    /// the dye's absorption — and the emission side is four more. Without this the only
    /// way to ask about the whole chain is to evaluate each part at each wavelength by
    /// hand, and the composite cannot be handed to anything.
    ///
    /// An empty product is 1, which is the identity and lets a chain be built up from
    /// nothing.
    Product {
        /// The factors. Flattened on construction, so `Mul` is exactly associative.
        factors: Vec<Spectrum>,
    },
    /// Several spectra added together.
    ///
    /// Two dyes in the same sample, two lamps in the same illuminator, a signal plus the
    /// leak through a filter's blocking. An empty sum is 0.
    Sum {
        /// The terms. Flattened on construction, as for [`Spectrum::Product`].
        terms: Vec<Spectrum>,
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

    /// A measured curve of `(nm, value)` pairs. Samples are sorted, so they may be
    /// given in any order.
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

    /// An emission or absorption band of the given centre and full width at half
    /// maximum.
    pub fn gaussian(center: Length, fwhm: Length, peak: f64) -> Spectrum {
        Spectrum::Gaussian {
            center_nm: center.in_nm(),
            fwhm_nm: fwhm.in_nm(),
            peak,
            floor: 0.0,
        }
    }

    /// Several spectra multiplied together.
    ///
    /// Nested products are flattened, which is not only tidiness: floating-point
    /// multiplication is commutative but *not* associative, so `(a·b)·c` and `a·(b·c)`
    /// could otherwise differ in their last bits. Flattening makes both spell the same
    /// left-to-right fold, so the operators are exactly associative and a chain built in
    /// any order gives identical numbers.
    pub fn product(factors: Vec<Spectrum>) -> Spectrum {
        let mut flat = Vec::with_capacity(factors.len());
        for factor in factors {
            match factor {
                Spectrum::Product { factors } => flat.extend(factors),
                other => flat.push(other),
            }
        }
        Spectrum::Product { factors: flat }
    }

    /// Several spectra added together. Nested sums are flattened, for the same reason.
    pub fn sum(terms: Vec<Spectrum>) -> Spectrum {
        let mut flat = Vec::with_capacity(terms.len());
        for term in terms {
            match term {
                Spectrum::Sum { terms } => flat.extend(terms),
                other => flat.push(other),
            }
        }
        Spectrum::Sum { terms: flat }
    }

    /// The same spectrum scaled so its largest value over `range` is exactly 1.
    pub fn normalized_peak(&self, range: (Length, Length), steps: usize) -> Spectrum {
        let peak = self.max_over(range, steps);
        if peak.abs() < f64::MIN_POSITIVE {
            return self.clone();
        }
        self.clone() * (1.0 / peak)
    }

    /// The same spectrum scaled so its integral over `range` is exactly 1.
    ///
    /// How an emission spectrum is normally carried: normalise the shape, then multiply by
    /// a total to get a distribution. Note that the result's values are per metre — the
    /// integral being 1 means the values had to absorb the wavelength axis — which is a
    /// real change of meaning and not just of scale.
    pub fn normalized_area(&self, range: (Length, Length), steps: usize) -> Spectrum {
        let area = self.integrate(range, steps);
        if area.abs() < f64::MIN_POSITIVE {
            return self.clone();
        }
        self.clone() * (1.0 / area)
    }

    /// The value at a wavelength.
    pub fn at(&self, wavelength: Length) -> f64 {
        self.at_nm(wavelength.in_nm())
    }

    /// The value at a wavelength given in nanometres.
    ///
    /// The dimensioned [`Spectrum::at`] is the one to reach for; this is what it
    /// calls, and what the sampling loops in this crate use so that a tight loop is
    /// not converting metres to nanometres and back on every step.
    pub fn at_nm(&self, wavelength_nm: f64) -> f64 {
        match self {
            Spectrum::Constant { value } => *value,
            Spectrum::Curve { samples } => interpolate(samples, wavelength_nm),
            Spectrum::Bands {
                bands,
                in_band,
                out_of_band,
                edge_nm,
            } => {
                // Combine bands by how far *inside* one the wavelength is, not by
                // which gives the larger value. Bands may overlap and the light
                // does not care which one caught it — and `in_band` is allowed to
                // be the smaller of the two, which is exactly the case for a
                // dichroic's reflectance.
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
                // Normalised to the peak, so the constants in front of Planck's law
                // cancel and only the shape survives.
                let at = |w: f64| planck_shape(w, *temperature_k);
                let peak_w = WIEN_NM_K / temperature_k;
                let denom = at(peak_w);
                if denom <= 0.0 {
                    return 0.0;
                }
                peak * (at(wavelength_nm) / denom)
            }
            // Folded left to right in the stored order, so the answer is a function of
            // the representation and not of how the expression was written.
            Spectrum::Product { factors } => factors
                .iter()
                .fold(1.0, |acc, f| acc * f.at_nm(wavelength_nm)),
            Spectrum::Sum { terms } => terms
                .iter()
                .fold(0.0, |acc, t| acc + t.at_nm(wavelength_nm)),
        }
    }

    /// Largest value anywhere in `range`, found by sampling. Used to check that
    /// reflectance and transmittance never sum past one.
    pub fn max_over(&self, range: (Length, Length), steps: usize) -> f64 {
        let (lo, hi) = (range.0.in_nm(), range.1.in_nm());
        let steps = steps.max(2);
        (0..=steps)
            .map(|i| self.at_nm(lo + (hi - lo) * i as f64 / steps as f64))
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Wavelength of the peak and its value, sampled over `range`.
    pub fn peak_over(&self, range: (Length, Length), steps: usize) -> (Length, f64) {
        let (lo, hi) = (range.0.in_nm(), range.1.in_nm());
        let steps = steps.max(2);
        let (nm, value) = (0..=steps)
            .map(|i| {
                let w = lo + (hi - lo) * i as f64 / steps as f64;
                (w, self.at_nm(w))
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
            );
        (Length::nm(nm), value)
    }

    /// Evenly spaced samples across `range`, for plotting.
    pub fn sample(&self, range: (Length, Length), steps: usize) -> Vec<(Length, f64)> {
        let (lo, hi) = (range.0.in_nm(), range.1.in_nm());
        let steps = steps.max(2);
        (0..=steps)
            .map(|i| {
                let w = lo + (hi - lo) * i as f64 / steps as f64;
                (Length::nm(w), self.at_nm(w))
            })
            .collect()
    }
}

impl Mul for Spectrum {
    type Output = Spectrum;

    /// A light path, written as one.
    fn mul(self, rhs: Spectrum) -> Spectrum {
        Spectrum::product(vec![self, rhs])
    }
}

impl Mul<f64> for Spectrum {
    type Output = Spectrum;

    /// Scaling, as a product with a flat spectrum — so it flattens into a chain like any
    /// other factor rather than needing a variant of its own.
    fn mul(self, k: f64) -> Spectrum {
        Spectrum::product(vec![self, Spectrum::constant(k)])
    }
}

impl Mul<Spectrum> for f64 {
    type Output = Spectrum;

    fn mul(self, s: Spectrum) -> Spectrum {
        s * self
    }
}

impl Add for Spectrum {
    type Output = Spectrum;

    fn add(self, rhs: Spectrum) -> Spectrum {
        Spectrum::sum(vec![self, rhs])
    }
}

/// Planck's law up to a constant: `1 / (l^5 (exp(c2 / (l T)) - 1))`, with the
/// wavelength in metres. Only the shape matters here, so the leading constants are
/// dropped and the result is normalised at the call site.
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

    fn nm(v: f64) -> Length {
        Length::nm(v)
    }

    #[test]
    fn a_constant_spectrum_is_flat() {
        let s = Spectrum::constant(0.42);
        for w in [200.0, 550.0, 2000.0] {
            assert_eq!(s.at(nm(w)), 0.42);
        }
    }

    /// A curve interpolates between its samples and holds its endpoints, rather
    /// than extrapolating a trend it has no evidence for.
    ///
    /// The comparisons are approximate rather than exact, and for a reason worth
    /// recording: a wavelength written as 450 nm is stored in metres and read back
    /// in nanometres, and `450e-9 * 1e9` is 449.999999999999_94. That is a
    /// femtometre, so it changes no physics, but it does mean an interpolated value
    /// is no longer bit-exact. See the note on brick-wall edges in
    /// [`Spectrum::Bands`].
    #[test]
    fn a_curve_interpolates_and_never_extrapolates() {
        let s = Spectrum::curve(vec![(500.0, 0.2), (600.0, 0.8), (400.0, 0.0)]);
        let close = |got: f64, want: f64| assert!((got - want).abs() < 1e-12, "{got} vs {want}");
        close(s.at(nm(400.0)), 0.0);
        close(s.at(nm(450.0)), 0.1);
        close(s.at(nm(550.0)), 0.5);
        close(s.at(nm(600.0)), 0.8);
        // Outside the measured range it holds, in both directions — and out there
        // the value is exact, because it is a stored endpoint rather than a
        // computed one.
        assert_eq!(s.at(nm(100.0)), 0.0);
        assert_eq!(s.at(nm(2000.0)), 0.8);
        // Reading in nanometres directly skips the round trip and is exact.
        assert_eq!(s.at_nm(450.0), 0.1);
    }

    /// A brick-wall band is exact at its edges; a real one rolls off across
    /// `edge_nm` and leaks its blocking floor outside.
    #[test]
    fn bands_have_edges_and_leak() {
        let brick = Spectrum::bands(vec![[500.0, 560.0]], 0.95, 0.0);
        assert_eq!(brick.at(nm(499.9)), 0.0);
        assert_eq!(brick.at(nm(500.0)), 0.95);
        assert_eq!(brick.at(nm(560.0)), 0.95);
        assert_eq!(brick.at(nm(560.1)), 0.0);

        let real = Spectrum::interference_bands(vec![[500.0, 560.0]], 0.95, 10.0);
        // Half way up the shoulder at the nominal edge.
        assert!((real.at(nm(500.0)) - (1e-6 + 0.95 / 2.0)).abs() < 1e-6);
        assert!(
            (real.at(nm(530.0)) - 0.95).abs() < 1e-9,
            "flat in the middle"
        );
        // Blocking is OD 6, not zero: a real filter leaks, and that leak is the
        // noise floor of a fluorescence image.
        assert_eq!(real.at(nm(400.0)), 1e-6);
        assert!(real.at(nm(400.0)) > 0.0);
    }

    /// Overlapping bands combine by how far inside one the wavelength is, since the
    /// light does not care which band caught it.
    #[test]
    fn overlapping_bands_take_the_best() {
        let s = Spectrum::bands(vec![[400.0, 500.0], [450.0, 600.0]], 0.9, 0.01);
        assert_eq!(s.at(nm(475.0)), 0.9);
        assert_eq!(s.at(nm(550.0)), 0.9);
        assert_eq!(s.at(nm(650.0)), 0.01);
    }

    /// `in_band` may be *lower* than `out_of_band`, which is how a dichroic's
    /// reflectance is written: high everywhere except where it transmits. Combining
    /// by value rather than by insideness would break this.
    #[test]
    fn a_band_may_invert() {
        let reflect = Spectrum::Bands {
            bands: vec![[495.0, 545.0]],
            in_band: 0.05,
            out_of_band: 0.98,
            edge_nm: 10.0,
        };
        assert!(
            (reflect.at(nm(520.0)) - 0.05).abs() < 1e-9,
            "got {}",
            reflect.at(nm(520.0))
        );
        assert!((reflect.at(nm(405.0)) - 0.98).abs() < 1e-9);
        assert!((reflect.at(nm(610.0)) - 0.98).abs() < 1e-9);
        // Halfway across the shoulder, halfway between the two levels.
        assert!((reflect.at(nm(495.0)) - (0.98 + (0.05 - 0.98) / 2.0)).abs() < 1e-9);
    }

    /// Planck's law, checked against Wien's displacement law: a blackbody peaks at
    /// 2898/T micrometres, and nowhere else. That closed form is independent of the
    /// implementation, which is what makes it worth testing against.
    #[test]
    fn a_blackbody_peaks_where_wien_says() {
        for temperature_k in [2000.0, 3200.0, 5500.0, 9000.0] {
            let s = Spectrum::blackbody(temperature_k);
            let expected = WIEN_NM_K / temperature_k;
            // Search a wide range so the peak is genuinely found, not assumed.
            let (peak_nm, peak) = s.peak_over((nm(100.0), nm(4000.0)), 4000);
            assert!(
                (peak_nm.in_nm() - expected).abs() < 3.0,
                "{temperature_k} K should peak at {expected:.0} nm, found {:.0} nm",
                peak_nm.in_nm()
            );
            // Normalised, so the value *at* the peak is exactly one; the sampled
            // maximum only lands within a grid step of it.
            assert!((s.at(nm(expected)) - 1.0).abs() < 1e-12);
            assert!(
                (peak - 1.0).abs() < 1e-4,
                "peak {peak} at {temperature_k} K"
            );
        }
    }

    /// The consequences a lamp choice actually has: tungsten is starved in the blue
    /// and dumps most of its output into the infrared, while daylight is nearly flat
    /// across the visible. This is why a halogen lamp looks orange.
    #[test]
    fn tungsten_is_red_and_daylight_is_not() {
        let tungsten = Spectrum::blackbody(3200.0);
        let daylight = Spectrum::blackbody(5500.0);
        let ratio = |s: &Spectrum| s.at(nm(450.0)) / s.at(nm(650.0));
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
        assert!(tungsten.peak_over((nm(300.0), nm(3000.0)), 3000).0.in_nm() > 850.0);
        // Nothing is emitted at zero or negative wavelength.
        assert_eq!(tungsten.at(Length::ZERO), 0.0);
        assert_eq!(Spectrum::blackbody(-5.0).at(nm(550.0)), 0.0);
        // A room-temperature body in the ultraviolet is not *zero* — Planck's law
        // has no cutoff — but it is unimaginably small, and the arithmetic must
        // reach that rather than overflowing on the way.
        let cold = Spectrum::blackbody(300.0).at(nm(200.0));
        assert!(cold.is_finite() && cold < 1e-80, "got {cold}");
        // Far enough into the ultraviolet the exponent does overflow, and the
        // answer there is zero rather than an infinity.
        assert_eq!(Spectrum::blackbody(300.0).at(nm(50.0)), 0.0);
    }

    /// A Gaussian is exactly half its peak at half its full width, which is what
    /// "full width at half maximum" means.
    #[test]
    fn a_gaussian_is_half_at_its_half_width() {
        let s = Spectrum::gaussian(nm(488.0), nm(20.0), 1.0);
        assert!((s.at(nm(488.0)) - 1.0).abs() < 1e-12);
        assert!(
            (s.at(nm(498.0)) - 0.5).abs() < 1e-9,
            "got {}",
            s.at(nm(498.0))
        );
        assert!((s.at(nm(478.0)) - 0.5).abs() < 1e-9);
        assert!(s.at(nm(588.0)) < 1e-9);
        let (peak_nm, peak) = s.peak_over(VISIBLE_RANGE, 750);
        assert!((peak_nm.in_nm() - 488.0).abs() < 2.0 && (peak - 1.0).abs() < 0.01);
    }

    /// A product is the pointwise product, which is the definition and the thing every
    /// other property here rests on.
    #[test]
    fn a_product_multiplies_wavelength_by_wavelength() {
        let a = Spectrum::gaussian(nm(500.0), nm(40.0), 0.9);
        let b = Spectrum::bands(vec![[480.0, 520.0]], 0.8, 0.01);
        let both = a.clone() * b.clone();
        for w in [400.0, 480.0, 500.0, 519.0, 560.0, 700.0] {
            let expected = a.at_nm(w) * b.at_nm(w);
            assert!(
                (both.at_nm(w) - expected).abs() < 1e-15,
                "at {w} nm: {} against {expected}",
                both.at_nm(w)
            );
        }
        // A sum likewise.
        let either = a.clone() + b.clone();
        for w in [450.0, 500.0, 600.0] {
            assert!((either.at_nm(w) - (a.at_nm(w) + b.at_nm(w))).abs() < 1e-15);
        }
    }

    /// Flattening makes the operators exactly associative, which they would not otherwise
    /// be: floating-point multiplication is commutative but not associative, so a nested
    /// representation would let `(a·b)·c` and `a·(b·c)` differ in their last bits.
    ///
    /// With both spelling the same left-to-right fold, a chain assembled in any order
    /// gives numbers identical to the bit — which is what the rest of this workspace
    /// promises about everything else.
    #[test]
    fn products_flatten_and_so_associate_exactly() {
        let a = Spectrum::gaussian(nm(480.0), nm(30.0), 0.7);
        let b = Spectrum::curve(vec![(400.0, 0.1), (600.0, 0.95)]);
        let c = Spectrum::constant(0.6);

        let left = (a.clone() * b.clone()) * c.clone();
        let right = a.clone() * (b.clone() * c.clone());
        assert_eq!(left, right, "the two spellings should be the same value");
        for w in [420.0, 480.0, 550.0] {
            assert_eq!(left.at_nm(w).to_bits(), right.at_nm(w).to_bits());
        }
        // And the representation is flat rather than nested three deep.
        match &left {
            Spectrum::Product { factors } => assert_eq!(factors.len(), 3),
            other => panic!("expected a flat product, got {other:?}"),
        }
    }

    /// The identities, so a chain can be built up from nothing.
    #[test]
    fn empty_products_and_sums_are_the_identities() {
        assert_eq!(Spectrum::product(vec![]).at_nm(550.0), 1.0);
        assert_eq!(Spectrum::sum(vec![]).at_nm(550.0), 0.0);
        // Which means folding a list of filters onto an empty product works.
        let filters = [
            Spectrum::constant(0.9),
            Spectrum::constant(0.8),
            Spectrum::constant(0.5),
        ];
        let chain = filters
            .iter()
            .cloned()
            .fold(Spectrum::product(vec![]), |acc, f| acc * f);
        assert!((chain.at_nm(550.0) - 0.36).abs() < 1e-15);
    }

    /// Scaling is a product with a flat spectrum, so it joins a chain as a factor rather
    /// than wrapping it — and it works from either side.
    #[test]
    fn scaling_joins_the_chain_rather_than_wrapping_it() {
        let s = Spectrum::gaussian(nm(520.0), nm(25.0), 1.0);
        let half = s.clone() * 0.5;
        assert!((half.at_nm(520.0) - 0.5).abs() < 1e-15);
        assert_eq!(2.0 * s.clone(), s.clone() * 2.0);
        // Two scalings and a spectrum flatten to three factors, not nested wrappers.
        match s * 0.5 * 4.0 {
            Spectrum::Product { factors } => assert_eq!(factors.len(), 3),
            other => panic!("expected a flat product, got {other:?}"),
        }
    }

    /// The connection to the integration already in `radiometry`: multiplying and then
    /// integrating is what `integrate_weighted` does, so the two must agree.
    ///
    /// That is also the trap this replaces. `integrate_weighted` takes exactly one weight,
    /// so a chain of seven had to be evaluated by hand at every wavelength; and the
    /// integral of a product is emphatically *not* the product of the integrals, which is
    /// the mistake a shortcut would be.
    #[test]
    fn multiplying_then_integrating_matches_the_weighted_integral() {
        use crate::radiometry::STEPS;
        let lamp = Spectrum::blackbody(3200.0);
        let filter = Spectrum::interference_bands(vec![[500.0, 560.0]], 0.95, 8.0);

        let by_product = (lamp.clone() * filter.clone()).integrate(VISIBLE_RANGE, STEPS);
        let by_weight = lamp.integrate_weighted(&filter, VISIBLE_RANGE, STEPS);
        assert!(
            (by_product / by_weight - 1.0).abs() < 1e-12,
            "{by_product:e} against {by_weight:e}"
        );

        // And the integral of a product is not the product of the integrals — not wrong by
        // a factor, but *not the same kind of quantity*. Each integral carries a metre
        // from the wavelength axis, so their product carries two while the integral of the
        // product carries one. The ratio between them is therefore a length, and it comes
        // out on the order of the range integrated over.
        //
        // Which is a better warning than "it is off by ten": a shortcut that changes the
        // dimensions cannot be rescued by a correction factor.
        let separately =
            lamp.integrate(VISIBLE_RANGE, STEPS) * filter.integrate(VISIBLE_RANGE, STEPS);
        let ratio = separately / by_product;
        let width = (VISIBLE_RANGE.1 - VISIBLE_RANGE.0).to_si();
        assert!(
            ratio / width > 0.3 && ratio / width < 3.0,
            "the ratio should be a length of order the range: {ratio:e} m against a \
             range of {width:e} m"
        );
    }

    /// **The chain this variant exists for.** A fluorescence path is seven spectra
    /// multiplied: lamp, excitation filter, dichroic transmission, dye absorption on the
    /// way in; dye emission, dichroic reflection and emission filter on the way out.
    ///
    /// Before this, the composite could be evaluated but not *held* — there was no way to
    /// hand the whole path to anything, and `integrate_weighted` takes one weight at a
    /// time.
    #[test]
    fn a_fluorescence_path_is_one_spectrum() {
        // GFP-like: excite near 488, emit near 509.
        let lamp = Spectrum::blackbody(5500.0);
        let excitation_filter = Spectrum::interference_bands(vec![[470.0, 495.0]], 0.95, 6.0);
        let dichroic_pass = Spectrum::interference_bands(vec![[460.0, 500.0]], 0.93, 8.0);
        let dye_absorption = Spectrum::gaussian(nm(488.0), nm(40.0), 1.0);

        let dye_emission = Spectrum::gaussian(nm(509.0), nm(35.0), 1.0);
        let dichroic_reflect = Spectrum::Bands {
            bands: vec![[460.0, 500.0]],
            in_band: 0.05,
            out_of_band: 0.97,
            edge_nm: 8.0,
        };
        let emission_filter = Spectrum::interference_bands(vec![[505.0, 560.0]], 0.94, 6.0);

        let excitation = lamp.clone()
            * excitation_filter.clone()
            * dichroic_pass.clone()
            * dye_absorption.clone();
        let emission = dye_emission.clone() * dichroic_reflect.clone() * emission_filter.clone();

        // Seven factors across two chains, and both are flat.
        match (&excitation, &emission) {
            (Spectrum::Product { factors: a }, Spectrum::Product { factors: b }) => {
                assert_eq!(a.len() + b.len(), 7);
            }
            _ => panic!("both chains should be flat products"),
        }

        // The excitation path peaks where the dye absorbs and is dark where it emits.
        let (peak_nm, _) = excitation.peak_over(VISIBLE_RANGE, 1500);
        assert!(
            (peak_nm.in_nm() - 488.0).abs() < 12.0,
            "excitation should peak near the dye's absorption, got {} nm",
            peak_nm.in_nm()
        );
        assert!(
            excitation.at(nm(509.0)) / excitation.at(nm(488.0)) < 0.1,
            "the excitation path must be dark where the dye emits, or the leak swamps \
             the signal"
        );

        // The emission path peaks where the dye emits and rejects the excitation band.
        let (emit_peak, _) = emission.peak_over(VISIBLE_RANGE, 1500);
        assert!(
            (emit_peak.in_nm() - 515.0).abs() < 15.0,
            "emission should peak near the dye, got {} nm",
            emit_peak.in_nm()
        );
        let rejection = emission.at(nm(488.0)) / emission.at(emit_peak);
        assert!(
            rejection < 1e-3,
            "the emission path should reject the excitation by at least a thousand, got \
             {rejection:e}"
        );

        // Every factor evaluated by hand at one wavelength, as the check that the chain is
        // the chain and not something adjacent to it.
        let w = 505.0;
        let by_hand = dye_emission.at_nm(w) * dichroic_reflect.at_nm(w) * emission_filter.at_nm(w);
        assert!((emission.at_nm(w) - by_hand).abs() < 1e-15);
    }

    /// Normalising, both ways, and the change of meaning one of them carries.
    #[test]
    fn normalising_scales_to_a_peak_or_to_an_area() {
        use crate::radiometry::STEPS;
        let band = Spectrum::gaussian(nm(550.0), nm(30.0), 7.3);

        let by_peak = band.normalized_peak(VISIBLE_RANGE, STEPS);
        assert!((by_peak.max_over(VISIBLE_RANGE, STEPS) - 1.0).abs() < 1e-9);
        // The shape is untouched: every ratio survives.
        assert!(
            (by_peak.at_nm(535.0) / by_peak.at_nm(550.0) - band.at_nm(535.0) / band.at_nm(550.0))
                .abs()
                < 1e-12
        );

        let by_area = band.normalized_area(VISIBLE_RANGE, STEPS);
        assert!((by_area.integrate(VISIBLE_RANGE, STEPS) - 1.0).abs() < 1e-9);
        // Which is a distribution rather than a shape: its values are per metre, so they
        // are enormous compared with the peak-normalised version.
        assert!(by_area.at_nm(550.0) > 1e6);

        // A dark spectrum cannot be normalised and is returned unchanged rather than as a
        // field of infinities.
        let dark = Spectrum::constant(0.0);
        assert_eq!(dark.normalized_peak(VISIBLE_RANGE, STEPS), dark);
        assert_eq!(dark.normalized_area(VISIBLE_RANGE, STEPS), dark);
    }

    /// A composed spectrum serialises, so a whole light path can live in a scene file.
    #[test]
    fn a_composed_path_round_trips_through_json() {
        let path = Spectrum::blackbody(3200.0)
            * Spectrum::interference_bands(vec![[500.0, 560.0]], 0.95, 8.0)
            * 0.5
            + Spectrum::constant(1e-6);
        let json = serde_json::to_string(&path).unwrap();
        let back: Spectrum = serde_json::from_str(&json).unwrap();
        assert_eq!(back, path);
        for w in [450.0, 530.0, 650.0] {
            assert_eq!(back.at_nm(w).to_bits(), path.at_nm(w).to_bits());
        }
        // The nesting is visible in the file rather than flattened away, since a sum of a
        // product is genuinely two levels.
        assert!(json.contains("\"sum\""), "{json}");
        assert!(json.contains("\"product\""), "{json}");
    }

    /// A long chain is a flat fold, not a deep recursion — so a path with many elements
    /// costs one pass and cannot overflow a stack.
    #[test]
    fn a_long_chain_stays_flat() {
        let mut chain = Spectrum::constant(1.0);
        for i in 0..2000 {
            chain = chain * Spectrum::constant(1.0 - 1e-6 * i as f64);
        }
        match &chain {
            Spectrum::Product { factors } => assert_eq!(factors.len(), 2001),
            other => panic!("expected a flat product, got {other:?}"),
        }
        assert!(chain.at_nm(550.0) > 0.0 && chain.at_nm(550.0) < 1.0);
    }

    /// The wavelength argument is a length, so a path length cannot be handed over
    /// as a wavelength — and 550 nm written any of three ways is one wavelength.
    #[test]
    fn wavelengths_are_lengths() {
        let s = Spectrum::gaussian(nm(550.0), nm(20.0), 1.0);
        assert_eq!(s.at(Length::nm(550.0)), s.at(Length::um(0.55)));
        assert_eq!(s.at(Length::nm(550.0)), s.at(Length::from_si(550e-9)));
        assert_eq!(s.at_nm(550.0), s.at(Length::nm(550.0)));
    }
}
