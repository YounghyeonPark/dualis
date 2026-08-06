//! How much light agrees with itself, and what that does to an image.
//!
//! Everything else in this crate is incoherent: intensities add. That is right for a
//! lamp, for fluorescence and for daylight, and wrong for a laser, wrong for a
//! condenser stopped down, and wrong for the interference that makes a microscope's
//! phase contrast work.
//!
//! # Two images of the same object, and they are not the same picture
//!
//! With coherent light, amplitudes add and *then* get squared:
//!
//! ```text
//! I = |h ⊛ U|²
//! ```
//!
//! With incoherent light the squaring happens first:
//!
//! ```text
//! I = |h|² ⊛ |U|²
//! ```
//!
//! The consequences are not a matter of degree. A coherent system's cutoff frequency is
//! **half** the incoherent one — `NA/λ` against `2NA/λ` — so it resolves worse; but
//! below that cutoff its transfer is flat at 1 rather than sloping away, so it renders
//! what it does pass at full contrast. And it rings: a sharp edge under coherent light
//! shows fringes that no amount of focusing removes, because they are interference and
//! not aberration.
//!
//! # Van Cittert–Zernike, which is the closed form this module is checked against
//!
//! An incoherent source does not stay incoherent at a distance. Light from a small
//! enough patch of sky arrives at two nearby points having travelled almost the same
//! path, and it interferes. The theorem says the coherence between two points is the
//! normalised Fourier transform of the source's intensity distribution — so for a
//! uniform circular source it is
//!
//! ```text
//! |μ(r)| = |2J₁(v)/v|,   v = π D r / (λ z)
//! ```
//!
//! the same function as the Airy amplitude, which
//! [`diffraction`](crate::diffraction) already computes from Bessel series. Two
//! completely different physical questions with the same answer, and the second one is
//! available here to check the first against.
//!
//! That is also how a star's diameter is measured: find the baseline at which the
//! fringes vanish, and the first zero of that function gives the angle.

use std::f64::consts::PI;

use dualis_units::Length;

use crate::diffraction::{bessel_j1, mtf_ideal, FIRST_AIRY_ZERO};

/// How coherent an illumination is, as the ratio of the condenser's aperture to the
/// objective's.
///
/// The number a microscopist actually sets, by opening or closing the condenser
/// diaphragm. Zero is a point source and fully coherent; 1 is a condenser matched to
/// the objective; above about 1 the illumination is effectively incoherent and opening
/// further only wastes light.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Coherence {
    /// `NA_condenser / NA_objective`, conventionally written `S` or `m`.
    pub ratio: f64,
}

impl Coherence {
    /// A point source: fully coherent.
    pub const COHERENT: Coherence = Coherence { ratio: 0.0 };
    /// Matched condenser, the usual compromise.
    pub const MATCHED: Coherence = Coherence { ratio: 1.0 };

    pub fn new(ratio: f64) -> Coherence {
        Coherence {
            ratio: ratio.max(0.0),
        }
    }

    /// Effectively incoherent, which for imaging purposes is anything past about 1.
    pub fn is_incoherent(&self) -> bool {
        self.ratio >= 1.0
    }

    /// Cutoff frequency as a multiple of the coherent cutoff `NA/λ`.
    ///
    /// `1 + S`: a coherent system stops at 1, a matched one at 2, and past that the
    /// objective's own aperture is the limit and nothing more is gained. This is why
    /// "closing the condenser sharpens the image" is false — it raises contrast at low
    /// frequencies and lowers the cutoff.
    pub fn cutoff_multiple(&self) -> f64 {
        1.0 + self.ratio.min(1.0)
    }
}

/// Cutoff frequency for coherent imaging, `NA/λ` in cycles per metre.
///
/// Exactly half the incoherent cutoff that
/// [`cutoff_frequency`](crate::diffraction::cutoff_frequency) returns. The factor of two
/// is the difference between a pupil and its autocorrelation: incoherent imaging
/// transfers the autocorrelation of the pupil, which is twice as wide as the pupil
/// itself.
pub fn coherent_cutoff(wavelength: Length, na: f64) -> f64 {
    if wavelength.to_si() <= 0.0 {
        return 0.0;
    }
    na / wavelength.to_si()
}

/// Coherent transfer at a fraction of the *coherent* cutoff.
///
/// A rectangle: 1 inside the pupil and 0 outside, because coherent imaging transfers the
/// pupil itself rather than its autocorrelation. Flat contrast up to the cutoff and
/// nothing past it, where an incoherent system slopes away from the start and reaches
/// twice as far.
pub fn coherent_transfer(fraction_of_coherent_cutoff: f64) -> f64 {
    f64::from(fraction_of_coherent_cutoff.abs() <= 1.0)
}

/// Complex degree of coherence between two points a distance apart, for light from a
/// uniform circular incoherent source.
///
/// Van Cittert–Zernike: `2J₁(v)/v` with `v = π D r/(λz)`. Real and signed here, since a
/// centred circular source has no phase — the negative lobes are real and mean the
/// fringes come back inverted.
pub fn van_cittert_zernike(
    separation: Length,
    source_diameter: Length,
    distance: Length,
    wavelength: Length,
) -> f64 {
    let (r, d, z, l) = (
        separation.to_si(),
        source_diameter.to_si(),
        distance.to_si(),
        wavelength.to_si(),
    );
    if z <= 0.0 || l <= 0.0 {
        return 1.0;
    }
    let v = PI * d * r / (l * z);
    if v.abs() < 1e-9 {
        // 2J1(v)/v -> 1 at the origin: a point is perfectly coherent with itself.
        return 1.0 - v * v / 8.0;
    }
    2.0 * bessel_j1(v) / v
}

/// Separation at which the fringes first vanish: `1.22 λz/D`.
///
/// The coherence area's radius, and the measurement a stellar interferometer makes. Find
/// the baseline at which the fringes disappear and this relation gives the source's
/// angular diameter — which is how the first star was measured, and how a lamp's
/// effective size can be found without looking at it.
pub fn coherence_radius(source_diameter: Length, distance: Length, wavelength: Length) -> Length {
    let d = source_diameter.to_si();
    if d <= 0.0 {
        return Length::from_si(f64::INFINITY);
    }
    Length::from_si(FIRST_AIRY_ZERO / PI * wavelength.to_si() * distance.to_si() / d)
}

/// Angular diameter of a source whose fringes vanish at this baseline: `1.22 λ/b`.
///
/// The inverse of [`coherence_radius`], and the form the measurement is actually used
/// in: Betelgeuse was measured this way in 1920 with a twenty-foot beam on a telescope,
/// long before anything could resolve it directly.
pub fn angular_diameter_from_baseline(baseline: Length, wavelength: Length) -> f64 {
    let b = baseline.to_si();
    if b <= 0.0 {
        return f64::INFINITY;
    }
    FIRST_AIRY_ZERO / PI * wavelength.to_si() / b
}

/// Transfer at a spatial frequency under partially coherent illumination, as a fraction
/// of the coherent cutoff.
///
/// Interpolated between the two limits rather than computed from the transmission
/// cross-coefficients, which for a general object are a four-dimensional integral and a
/// different module. The two ends are exact — a rectangle at `S = 0` and the disc
/// autocorrelation at `S ≥ 1` — and the shape between them carries the trade that
/// matters: contrast at low frequencies against reach at high ones.
pub fn partially_coherent_transfer(coherence: Coherence, fraction_of_coherent_cutoff: f64) -> f64 {
    let s = coherence.ratio.min(1.0);
    let f = fraction_of_coherent_cutoff.abs();
    let cutoff = 1.0 + s;
    if f > cutoff {
        return 0.0;
    }
    if s <= 0.0 {
        return coherent_transfer(f);
    }
    // The incoherent limit is the disc autocorrelation, expressed against the coherent
    // cutoff by halving the argument.
    let incoherent = mtf_ideal(f / 2.0);
    if s >= 1.0 {
        return incoherent;
    }
    let coherent = f64::from(f <= 1.0);
    coherent * (1.0 - s) + incoherent * s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diffraction::cutoff_frequency;

    fn nm(v: f64) -> Length {
        Length::nm(v)
    }

    /// The factor of two, which is the whole difference between the two imaging modes
    /// and the commonest thing to get backwards.
    #[test]
    fn coherent_imaging_cuts_off_at_half_the_incoherent_frequency() {
        let green = nm(550.0);
        for na in [0.25f64, 0.7, 1.4] {
            let coherent = coherent_cutoff(green, na);
            let incoherent = cutoff_frequency(green, na);
            assert!(
                (incoherent / coherent - 2.0).abs() < 1e-12,
                "NA {na}: incoherent should be exactly twice coherent"
            );
        }
        // 1.4 NA at 550 nm: 2545 cycles per mm coherently, 5091 incoherently.
        assert!((coherent_cutoff(green, 1.4) / 1e3 - 2545.5).abs() < 1.0);
    }

    /// A coherent system passes less but passes it at full contrast; an incoherent one
    /// reaches further but fades all the way. Neither is simply better, which is why the
    /// condenser has a diaphragm.
    #[test]
    fn coherent_transfer_is_flat_and_incoherent_transfer_slopes() {
        // Flat to the cutoff, then nothing.
        for f in [0.0f64, 0.3, 0.7, 0.99] {
            assert_eq!(coherent_transfer(f), 1.0, "flat at {f}");
        }
        assert_eq!(coherent_transfer(1.01), 0.0);

        // The incoherent one has already lost more than half its contrast by the
        // coherent cutoff, and still has as far again to go.
        let at_coherent_cutoff = mtf_ideal(0.5);
        assert!(
            (at_coherent_cutoff - 0.391).abs() < 0.001,
            "0.391 at half the incoherent cutoff"
        );
        assert!(
            at_coherent_cutoff < coherent_transfer(1.0),
            "coherent light renders what it passes at better contrast"
        );
    }

    /// Partially coherent illumination trades reach against contrast, and the endpoints
    /// are the two exact cases.
    #[test]
    fn partial_coherence_trades_reach_for_contrast() {
        // A point source is the rectangle.
        assert_eq!(partially_coherent_transfer(Coherence::COHERENT, 0.8), 1.0);
        assert_eq!(partially_coherent_transfer(Coherence::COHERENT, 1.2), 0.0);

        // A matched condenser is the incoherent disc autocorrelation.
        let matched = partially_coherent_transfer(Coherence::MATCHED, 1.0);
        assert!((matched - mtf_ideal(0.5)).abs() < 1e-12);
        assert_eq!(partially_coherent_transfer(Coherence::MATCHED, 2.01), 0.0);

        // The cutoff opens as the condenser does, from 1 to 2 and no further.
        assert_eq!(Coherence::COHERENT.cutoff_multiple(), 1.0);
        assert_eq!(Coherence::new(0.5).cutoff_multiple(), 1.5);
        assert_eq!(Coherence::MATCHED.cutoff_multiple(), 2.0);
        assert_eq!(
            Coherence::new(3.0).cutoff_multiple(),
            2.0,
            "opening past the objective wastes light and buys nothing"
        );
        assert!(Coherence::new(1.5).is_incoherent());

        // Closing the condenser raises contrast at low frequencies and *lowers* the
        // cutoff. Both halves of that are worth asserting, because only the first half
        // is folklore.
        let low = 0.6;
        assert!(
            partially_coherent_transfer(Coherence::new(0.2), low)
                > partially_coherent_transfer(Coherence::MATCHED, low),
            "a closed condenser has more contrast at low frequencies"
        );
        let high = 1.6;
        assert!(
            partially_coherent_transfer(Coherence::new(0.2), high)
                < partially_coherent_transfer(Coherence::MATCHED, high),
            "and less reach at high ones"
        );
    }

    /// **Van Cittert–Zernike against the Bessel form.** The coherence of light from a
    /// uniform circular source is `2J₁(v)/v` — the same function as the Airy amplitude,
    /// arrived at from a completely different question.
    #[test]
    fn the_coherence_of_a_circular_source_is_the_airy_function() {
        let wavelength = nm(550.0);
        let source = Length::mm(1.0);
        let distance = Length::m(2.0);

        // Perfectly coherent with itself at zero separation.
        assert!(
            (van_cittert_zernike(Length::ZERO, source, distance, wavelength) - 1.0).abs() < 1e-12
        );

        // And the profile is 2J1(v)/v with v = pi D r / (lambda z), checked point by
        // point against the Bessel function `diffraction` computes.
        for r_mm in [0.1f64, 0.3, 0.67, 1.2, 2.0] {
            let r = Length::mm(r_mm);
            let v = PI * source.to_si() * r.to_si() / (wavelength.to_si() * distance.to_si());
            let expected = 2.0 * bessel_j1(v) / v;
            let measured = van_cittert_zernike(r, source, distance, wavelength);
            assert!(
                (measured - expected).abs() < 1e-12,
                "at {r_mm} mm: {measured} against {expected}"
            );
        }

        // Falls monotonically to the first zero, and the negative lobe past it is real:
        // the fringes come back with their bright and dark swapped.
        let radius = coherence_radius(source, distance, wavelength);
        assert!(
            van_cittert_zernike(radius, source, distance, wavelength).abs() < 1e-6,
            "the fringes should vanish at the coherence radius"
        );
        assert!(
            van_cittert_zernike(radius * 1.3, source, distance, wavelength) < -0.01,
            "and come back inverted past it"
        );
    }

    /// The measurement the theorem is used for: a source's size from the baseline at
    /// which its fringes disappear.
    ///
    /// 1.22 λ/b, which is how Betelgeuse was measured in 1920 with a twenty-foot beam,
    /// years before anything could resolve it directly.
    #[test]
    fn a_source_can_be_measured_by_where_its_fringes_vanish() {
        let wavelength = nm(550.0);
        // A 1 mm lamp at 2 m: the coherence radius is 1.34 mm.
        let radius = coherence_radius(Length::mm(1.0), Length::m(2.0), wavelength);
        assert!(
            (radius.in_mm() - 1.342).abs() < 0.005,
            "coherence radius {} mm",
            radius.in_mm()
        );

        // Round trip: the angular diameter recovered from that baseline is the source's.
        let angular = angular_diameter_from_baseline(radius, wavelength);
        let actual = 1e-3 / 2.0; // 1 mm at 2 m
        assert!(
            (angular / actual - 1.0).abs() < 1e-12,
            "recovered {angular:e} rad against an actual {actual:e}"
        );

        // Betelgeuse: 0.047 arcseconds, which is 2.3e-7 radians, needs a 2.9 m baseline
        // at 550 nm.
        let betelgeuse = 0.047 / 3600.0 * PI / 180.0;
        let baseline = FIRST_AIRY_ZERO / PI * wavelength.to_si() / betelgeuse;
        assert!(
            (baseline - 2.95).abs() < 0.1,
            "a twenty-foot beam, near enough: {baseline:.2} m"
        );

        // A point source is coherent everywhere, and says so rather than dividing by
        // zero.
        assert!(!coherence_radius(Length::ZERO, Length::m(1.0), wavelength).is_finite());
        assert!(!angular_diameter_from_baseline(Length::ZERO, wavelength).is_finite());
    }

    /// A bigger source, or a closer one, is less coherent. Both halves of that are the
    /// same statement and both are worth having the sign of.
    #[test]
    fn coherence_shrinks_with_source_size_and_grows_with_distance() {
        let wavelength = nm(550.0);
        let at = |d_mm: f64, z_m: f64| {
            coherence_radius(Length::mm(d_mm), Length::m(z_m), wavelength).in_mm()
        };
        // Twice the source, half the coherence radius.
        assert!((at(1.0, 2.0) / at(2.0, 2.0) - 2.0).abs() < 1e-12);
        // Twice the distance, twice the coherence radius.
        assert!((at(1.0, 4.0) / at(1.0, 2.0) - 2.0).abs() < 1e-12);
        // The sun is half a degree across, so daylight is coherent over about 60
        // micrometres — which is why you cannot see interference in it without a
        // pinhole.
        let sun = 0.53 / 180.0 * PI;
        let solar_coherence = FIRST_AIRY_ZERO / PI * wavelength.to_si() / sun;
        assert!(
            (solar_coherence * 1e6 - 72.6).abs() < 2.0,
            "daylight is coherent over {} um",
            solar_coherence * 1e6
        );
    }

    /// Degenerate arguments give the limit rather than a NaN.
    #[test]
    fn degenerate_geometry_is_handled() {
        let w = nm(550.0);
        // Zero distance: nothing has propagated, so nothing has become coherent.
        assert_eq!(
            van_cittert_zernike(Length::mm(1.0), Length::mm(1.0), Length::ZERO, w),
            1.0
        );
        // A point source is coherent at every separation.
        let point = van_cittert_zernike(Length::m(1.0), Length::ZERO, Length::m(1.0), w);
        assert!((point - 1.0).abs() < 1e-9, "got {point}");
        assert_eq!(coherent_cutoff(Length::ZERO, 1.0), 0.0);
    }
}
