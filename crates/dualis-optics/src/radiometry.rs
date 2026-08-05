//! Integrating a spectrum, and the difference between watts and photons.
//!
//! A [`Spectrum`] is a shape. It says a filter passes 0.95 at 530 nm and 1e-6 at
//! 405 nm, and that is enough to answer questions about ratios, but not to answer
//! "how much light gets through" — for that the shape has to be integrated against
//! something, and the answer has to have units.
//!
//! # Watts are not photons, and a detector counts photons
//!
//! A photon at 450 nm carries 1.44 times the energy of one at 650 nm, so a
//! milliwatt of blue light is 1.44 times *fewer* photons than a milliwatt of red.
//! Every silicon detector responds to the count, not the power. Getting this
//! backwards is a systematic error that grows with how broad the source is — it is
//! invisible for a laser and worth tens of percent for a tungsten lamp, which is
//! exactly the case where nobody notices.
//!
//! So [`SpectralPower::total`] returns a [`Power`] and
//! [`SpectralPower::photon_rate`] returns a rate, and the two are different
//! functions of the same distribution. The bridge is `E = hc/λ`, which lives in
//! `dualis-units` where its dimensions can be checked.
//!
//! # Photon counts have no dimension
//!
//! SI gives a count no dimension, so a photon rate is dimensionally a
//! [`Frequency`] — the same limitation that makes an angular velocity
//! indistinguishable from a frequency. That is not only a wart: `rate * exposure`
//! really is the dimensionless number of photons collected, and the type system
//! gets that right.

use dualis_units::{Dimensionless, Frequency, Length, Power, Time, C, PLANCK};

use crate::spectrum::Spectrum;

/// Integration steps that are enough for anything in this crate.
///
/// Trapezoid over a smooth spectrum converges as the square of the step, and a
/// filter edge a few nanometres wide is the sharpest feature here — 2000 intervals
/// across the working range is under half a nanometre each, which resolves it.
/// A brick-wall edge is a discontinuity and no step size fixes that; give a real
/// filter a real edge width.
pub const STEPS: usize = 2000;

impl Spectrum {
    /// ∫ S(λ) dλ over the range, by the trapezoid rule.
    ///
    /// The result is in value-metres, which is rarely what anyone wants on its own
    /// — the useful quantities built on it are [`Spectrum::weighted_mean`]
    /// (dimensionless), [`Spectrum::centroid`] (a length) and
    /// [`SpectralPower::through`] (a power).
    pub fn integrate(&self, range: (Length, Length), steps: usize) -> f64 {
        self.integrate_with(range, steps, |_, s| s)
    }

    /// ∫ S(λ)·W(λ) dλ — the shape of every throughput calculation: a lamp against a
    /// filter, a filter against a detector's quantum efficiency.
    pub fn integrate_weighted(
        &self,
        weight: &Spectrum,
        range: (Length, Length),
        steps: usize,
    ) -> f64 {
        self.integrate_with(range, steps, |nm, s| s * weight.at_nm(nm))
    }

    /// ∫S·W dλ / ∫S dλ — what fraction of this spectrum `W` lets through.
    ///
    /// Dimensionless, and the number to reach for: "this dichroic passes 91% of that
    /// LED" is this function, and it does not depend on how the LED's shape happened
    /// to be normalised.
    pub fn weighted_mean(&self, weight: &Spectrum, range: (Length, Length), steps: usize) -> f64 {
        let total = self.integrate(range, steps);
        if total.abs() < f64::MIN_POSITIVE {
            return 0.0;
        }
        self.integrate_weighted(weight, range, steps) / total
    }

    /// ∫λS dλ / ∫S dλ — where the spectrum's weight actually sits.
    ///
    /// Not the same as the peak: a blackbody's centroid is well to the red of its
    /// Wien peak, because the long-wavelength tail is far heavier than the short one.
    pub fn centroid(&self, range: (Length, Length), steps: usize) -> Length {
        let total = self.integrate(range, steps);
        if total.abs() < f64::MIN_POSITIVE {
            return Length::ZERO;
        }
        let moment = self.integrate_with(range, steps, |nm, s| s * nm * 1e-9);
        Length::from_si(moment / total)
    }

    /// ∫S dλ divided by the peak value — the width of the brick-wall filter that
    /// would pass the same total.
    ///
    /// This is how a band is quoted when its shape is awkward: a Gaussian of 20 nm
    /// FWHM has an equivalent width of 21.3 nm, and a real interference filter's is
    /// what a photon budget should be built on rather than its nominal width.
    pub fn equivalent_width(&self, range: (Length, Length), steps: usize) -> Length {
        let peak = self.max_over(range, steps);
        if peak.abs() < f64::MIN_POSITIVE {
            return Length::ZERO;
        }
        Length::from_si(self.integrate(range, steps) / peak)
    }

    /// Trapezoid over the range, in nanometres, with the integrand built per sample.
    /// The `1e-9` at the end turns the nanometre step into metres, so the result is
    /// in value-metres like every other integral here.
    fn integrate_with(
        &self,
        range: (Length, Length),
        steps: usize,
        integrand: impl Fn(f64, f64) -> f64,
    ) -> f64 {
        let (lo, hi) = (range.0.in_nm(), range.1.in_nm());
        if hi <= lo {
            return 0.0;
        }
        let steps = steps.max(2);
        let h = (hi - lo) / steps as f64;
        let mut sum = 0.0;
        for i in 0..=steps {
            let nm = lo + h * i as f64;
            let value = integrand(nm, self.at_nm(nm));
            // Trapezoid: the endpoints count half.
            sum += if i == 0 || i == steps {
                value / 2.0
            } else {
                value
            };
        }
        sum * h * 1e-9
    }
}

/// A source of known total power and known spectral shape.
///
/// The shape's normalisation does not matter — every method divides by the shape's
/// own integral — so a [`Spectrum::blackbody`] scaled to a peak of 1 and the same
/// blackbody scaled to 1e6 describe the same lamp.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralPower {
    shape: Spectrum,
    total: Power,
    range: (Length, Length),
    steps: usize,
}

impl SpectralPower {
    /// A source emitting `total` watts, distributed as `shape` over `range`.
    pub fn new(shape: Spectrum, total: Power, range: (Length, Length)) -> SpectralPower {
        SpectralPower {
            shape,
            total,
            range,
            steps: STEPS,
        }
    }

    /// Override the integration resolution.
    pub fn with_steps(mut self, steps: usize) -> SpectralPower {
        self.steps = steps.max(2);
        self
    }

    pub fn shape(&self) -> &Spectrum {
        &self.shape
    }

    pub fn range(&self) -> (Length, Length) {
        self.range
    }

    /// The power emitted over the whole range.
    pub fn total(&self) -> Power {
        self.total
    }

    /// Power surviving a spectral transmission — a filter, a coating, a path through
    /// glass, or all three multiplied together.
    pub fn through(&self, transmission: &Spectrum) -> Power {
        self.total
            * self
                .shape
                .weighted_mean(transmission, self.range, self.steps)
    }

    /// Power-weighted mean wavelength. This is what sets the photon rate, which is
    /// why it is worth a name.
    pub fn centroid(&self) -> Length {
        self.shape.centroid(self.range, self.steps)
    }

    /// Photons per second.
    ///
    /// `∫ P(λ)·λ/(hc) dλ`, which reduces exactly to `total × centroid / (hc)` — so a
    /// source's photon rate is set by its *power-weighted mean wavelength*, and two
    /// lamps of equal wattage differ in photon output by the ratio of their
    /// centroids.
    pub fn photon_rate(&self) -> Frequency {
        Frequency::from_si(self.total.to_si() * self.centroid().to_si() / hc())
    }

    /// Photons per second surviving a transmission.
    ///
    /// Not `photon_rate() * through()/total()`: filtering changes the centroid as
    /// well as the total, so the photon rate has to be integrated against the
    /// filtered distribution. A blue-blocking filter removes the photons that were
    /// worth least, and this is the difference.
    pub fn photon_rate_through(&self, transmission: &Spectrum) -> Frequency {
        let norm = self.shape.integrate(self.range, self.steps);
        if norm.abs() < f64::MIN_POSITIVE {
            return Frequency::ZERO;
        }
        // ∫ s·f·λ dλ / ∫ s dλ, in metres.
        let weighted = self.shape.integrate_with(self.range, self.steps, |nm, s| {
            s * transmission.at_nm(nm) * nm * 1e-9
        });
        Frequency::from_si(self.total.to_si() * (weighted / norm) / hc())
    }

    /// Photons collected over an exposure — dimensionless, because a count is.
    pub fn photons_in(&self, exposure: Time) -> Dimensionless {
        self.photon_rate() * exposure
    }

    /// Power absorbed given a spectral absorptance — the number a thermal domain
    /// receives.
    ///
    /// The same integral as [`SpectralPower::through`]; it has its own name because
    /// it is the seam between optics and heat, and because absorptance is a
    /// different curve from transmittance even when it is one minus it.
    pub fn absorbed_by(&self, absorptance: &Spectrum) -> Power {
        self.through(absorptance)
    }
}

/// `hc`, in joule-metres. The product that turns a wavelength into a photon energy.
fn hc() -> f64 {
    PLANCK.to_si() * C.to_si()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spectrum::VISIBLE_RANGE;
    use crate::SurfaceOptics;
    use dualis_units::photon_energy;

    fn nm(v: f64) -> Length {
        Length::nm(v)
    }

    /// A flat spectrum over a known interval integrates to the rectangle, exactly.
    #[test]
    fn a_constant_integrates_to_its_rectangle() {
        let flat = Spectrum::constant(0.5);
        let got = flat.integrate((nm(400.0), nm(500.0)), 100);
        // 0.5 over 100 nm is 0.5 * 100e-9 value-metres.
        assert!((got - 0.5 * 100e-9).abs() < 1e-20, "got {got:e}");
        // A backwards or empty range integrates to nothing rather than to a
        // negative area.
        assert_eq!(flat.integrate((nm(500.0), nm(400.0)), 100), 0.0);
        assert_eq!(flat.integrate((nm(500.0), nm(500.0)), 100), 0.0);
    }

    /// A Gaussian's integral has a closed form — `peak · FWHM · sqrt(pi/(4 ln2))`,
    /// which is 1.0645 · peak · FWHM — and that is what the equivalent width
    /// measures. A 20 nm Gaussian passes as much as a 21.3 nm brick wall.
    #[test]
    fn a_gaussian_matches_its_closed_form_area() {
        let g = Spectrum::gaussian(nm(550.0), nm(20.0), 1.0);
        let expected = (std::f64::consts::PI / (4.0 * 2f64.ln())).sqrt() * 20e-9;
        let got = g.integrate(VISIBLE_RANGE, STEPS);
        assert!(
            (got / expected - 1.0).abs() < 1e-6,
            "got {got:e}, expected {expected:e}"
        );
        let width = g.equivalent_width(VISIBLE_RANGE, STEPS);
        assert!(
            (width.in_nm() - 21.29).abs() < 0.01,
            "equivalent width {} nm",
            width.in_nm()
        );
        // A symmetric band's centroid is its centre.
        assert!((g.centroid(VISIBLE_RANGE, STEPS).in_nm() - 550.0).abs() < 0.01);
    }

    /// The question that could not be asked before: how much of this lamp gets
    /// through this filter. A brick-wall band over a flat source passes exactly the
    /// ratio of the widths.
    #[test]
    fn throughput_is_the_ratio_of_the_widths() {
        let flat = Spectrum::constant(1.0);
        let band = Spectrum::bands(vec![[500.0, 600.0]], 1.0, 0.0);
        let fraction = flat.weighted_mean(&band, (nm(400.0), nm(800.0)), 4000);
        // 100 nm passed out of 400 nm offered.
        assert!((fraction - 0.25).abs() < 1e-3, "got {fraction}");
        // And the normalisation of the source cannot matter.
        let scaled = Spectrum::constant(1e6);
        let same = scaled.weighted_mean(&band, (nm(400.0), nm(800.0)), 4000);
        assert!((same - fraction).abs() < 1e-12);
    }

    /// A monochromatic source's photon rate is exactly its power over `hc/λ`, which
    /// is the definition, so this is the anchor the broadband cases are trusted
    /// against.
    #[test]
    fn a_narrow_source_gives_the_textbook_photon_rate() {
        // A 1 nm line is narrow enough to be monochromatic for this purpose.
        let laser = SpectralPower::new(
            Spectrum::gaussian(nm(532.0), nm(1.0), 1.0),
            Power::mw(1.0),
            (nm(520.0), nm(545.0)),
        );
        let expected = Power::mw(1.0).to_si() / photon_energy(nm(532.0)).to_si();
        let got = laser.photon_rate().to_si();
        assert!(
            (got / expected - 1.0).abs() < 1e-6,
            "got {got:e}, expected {expected:e}"
        );
        // 2.7e15 photons a second in a milliwatt of green.
        assert!((got / 2.678e15 - 1.0).abs() < 0.01, "got {got:e}");
        // Over a 10 ms exposure that is 2.7e13 photons, dimensionless.
        let count = laser.photons_in(Time::ms(10.0));
        assert!((count.to_si() / 2.678e13 - 1.0).abs() < 0.01, "{count:?}");
    }

    /// The systematic error this module exists to prevent: equal power at different
    /// wavelengths is not equal photons. Red beats blue by the ratio of the
    /// wavelengths, 650/450 = 1.44.
    #[test]
    fn equal_power_is_not_equal_photons() {
        let line = |centre_nm: f64| {
            SpectralPower::new(
                Spectrum::gaussian(nm(centre_nm), nm(1.0), 1.0),
                Power::mw(1.0),
                (nm(centre_nm - 12.0), nm(centre_nm + 12.0)),
            )
        };
        let blue = line(450.0).photon_rate().to_si();
        let red = line(650.0).photon_rate().to_si();
        assert!(
            (red / blue - 650.0 / 450.0).abs() < 1e-4,
            "ratio {}",
            red / blue
        );
        assert!(red > blue, "a red photon is cheaper, so a watt buys more");
    }

    /// A blackbody's centroid sits well to the red of its Wien peak, because the
    /// long-wavelength tail carries far more weight than the short one. Anyone who
    /// used the peak as the mean wavelength would underestimate a tungsten lamp's
    /// photon output.
    #[test]
    fn a_blackbody_centroid_is_redder_than_its_peak() {
        let range = (nm(350.0), nm(2500.0));
        let tungsten = Spectrum::blackbody(3200.0);
        let peak = tungsten.peak_over(range, 4000).0;
        let centroid = tungsten.centroid(range, 4000);
        assert!(
            (peak.in_nm() - 905.0).abs() < 5.0,
            "Wien peak {} nm",
            peak.in_nm()
        );
        assert!(
            centroid > peak,
            "centroid {} nm should be redder than the peak {} nm",
            centroid.in_nm(),
            peak.in_nm()
        );
        assert!(
            centroid.in_nm() > 1000.0,
            "centroid {} nm",
            centroid.in_nm()
        );
    }

    /// Filtering changes the mean wavelength as well as the total, so the photon
    /// rate through a filter is not the power fraction times the unfiltered rate.
    /// A blue-passing filter keeps the photons that were worth least per photon.
    #[test]
    fn filtering_moves_the_centroid_not_just_the_total() {
        let lamp = SpectralPower::new(
            Spectrum::blackbody(3200.0),
            Power::w(1.0),
            (nm(350.0), nm(2500.0)),
        );
        let blue = Spectrum::bands(vec![[400.0, 500.0]], 1.0, 0.0);

        let exact = lamp.photon_rate_through(&blue).to_si();
        let naive =
            lamp.photon_rate().to_si() * (lamp.through(&blue).to_si() / lamp.total().to_si());
        // The naive scaling uses the unfiltered centroid (about 1.2 um) where the
        // truth uses the filtered one (about 450 nm), so it overestimates by nearly
        // three times.
        assert!(
            naive / exact > 2.0,
            "the naive estimate should be badly high: naive {naive:e} exact {exact:e}"
        );
        // Sanity: the filtered rate is a small part of the whole, and positive.
        assert!(exact > 0.0 && exact < lamp.photon_rate().to_si());
    }

    /// The seam to the thermal domain, computed end to end: a real dichroic absorbs
    /// a few percent of a lamp, and that is a definite number of watts.
    #[test]
    fn absorbed_power_is_a_definite_number_of_watts() {
        let lamp = SpectralPower::new(
            Spectrum::blackbody(3200.0),
            Power::w(2.0),
            crate::spectrum::VISIBLE_RANGE,
        );
        let dichroic = SurfaceOptics::dichroic(vec![[495.0, 545.0]], 0.95, 10.0);

        // Absorptance sampled as its own spectrum, which is what the coupling hands
        // over.
        let absorptance = Spectrum::curve(
            (0..=75)
                .map(|i| {
                    let w = 350.0 + i as f64 * 10.0;
                    (w, dichroic.absorptance(nm(w)))
                })
                .collect(),
        );
        let absorbed = lamp.absorbed_by(&absorptance);
        assert!(
            absorbed.in_mw() > 1.0 && absorbed.in_mw() < 60.0,
            "a dichroic should absorb a few tens of milliwatts of a 2 W lamp, got {:?}",
            absorbed.in_mw()
        );
        // And it cannot absorb more than arrived.
        assert!(absorbed < lamp.total());

        // A black surface absorbs the lot, to the last part in a million.
        let all = lamp.absorbed_by(&Spectrum::constant(1.0));
        assert!((all / lamp.total() - 1.0).abs() < 1e-9);
        // A perfect mirror absorbs none.
        assert_eq!(lamp.absorbed_by(&Spectrum::constant(0.0)).to_si(), 0.0);
    }

    /// How good the trapezoid actually is here, which is better than "second order"
    /// in both of the cases that come up.
    ///
    /// On a piecewise-linear [`Spectrum::Curve`] it is exact by construction, at any
    /// step count. On a band that decays smoothly to nothing inside the range it is
    /// exact to machine precision at a hundred intervals — because every derivative
    /// vanishes at the endpoints and the Euler-Maclaurin correction terms vanish
    /// with them, which makes the convergence exponential in the step rather than
    /// quadratic. Neither case is the O(h²) a textbook quotes, and knowing that is
    /// what says [`STEPS`] is generous rather than marginal.
    #[test]
    fn the_trapezoid_is_exact_on_a_ramp_and_on_a_decaying_band() {
        // A ramp from 0 at 400 nm to 1 at 800 nm: a triangle of area 2e-7.
        let ramp = Spectrum::curve(vec![(400.0, 0.0), (800.0, 1.0)]);
        let triangle = 0.5 * 400e-9;
        for steps in [4, 17, 200] {
            let got = ramp.integrate((nm(400.0), nm(800.0)), steps);
            // Exact up to the rounding of summing the samples — a handful of ulp,
            // and independent of how many there were.
            assert!(
                (got / triangle - 1.0).abs() < 1e-14,
                "{steps} steps gave {got:e}, expected {triangle:e}"
            );
        }

        // A 40 nm FWHM band is sigma = 17 nm wide. Across the 750 nm working range:
        let g = Spectrum::gaussian(nm(550.0), nm(40.0), 1.0);
        let exact = (std::f64::consts::PI / (4.0 * 2f64.ln())).sqrt() * 40e-9;
        let error_at = |steps: usize| (g.integrate(VISIBLE_RANGE, steps) / exact - 1.0).abs();

        // 100 intervals is h = 7.5 nm, comfortably inside sigma: machine precision.
        assert!(error_at(100) < 1e-12, "off by {:e}", error_at(100));
        // 44 intervals is h ~ sigma, which is where the exponential term is still
        // tiny — a part in a hundred million.
        assert!(error_at(44) < 1e-6, "off by {:e}", error_at(44));
        // 12 intervals is h = 62 nm, nearly four sigma: the grid steps straight over
        // the band and the error collapses to 14%. The exponential convergence above
        // holds only once the step resolves the feature, and below that there is a
        // cliff rather than a gradient — which is why [`STEPS`] is set from the
        // narrowest filter edge worth representing and not from a convergence study
        // on a smooth curve.
        assert!(
            error_at(12) > 0.1,
            "under-resolving should fail badly, not gently: off by {:e}",
            error_at(12)
        );
    }

    /// A spectrum with no area does not divide by zero on its way to an answer.
    #[test]
    fn an_empty_spectrum_has_no_centroid_rather_than_a_nan() {
        let nothing = Spectrum::constant(0.0);
        assert_eq!(nothing.centroid(VISIBLE_RANGE, STEPS), Length::ZERO);
        assert_eq!(nothing.equivalent_width(VISIBLE_RANGE, STEPS), Length::ZERO);
        assert_eq!(
            nothing.weighted_mean(&Spectrum::constant(1.0), VISIBLE_RANGE, STEPS),
            0.0
        );
        let dark = SpectralPower::new(nothing, Power::w(0.0), VISIBLE_RANGE);
        assert_eq!(dark.photon_rate(), Frequency::ZERO);
        assert_eq!(
            dark.photon_rate_through(&Spectrum::constant(1.0)),
            Frequency::ZERO
        );
    }
}
