//! Wave optics: what light does that a ray cannot express.
//!
//! Everything else in this crate is geometric. A ray hits a surface, bends, and
//! arrives somewhere — and in that picture a perfect lens focuses a point source to
//! a point, so resolution is limited only by how well the lens is made. That
//! picture is wrong, and for a microscope it is wrong about the thing that matters
//! most: a perfect lens focuses a point to a disc about `0.61 λ/NA` across, because
//! light is a wave and the aperture that collected it was finite. No amount of
//! polishing improves on that.
//!
//! # Why this module is testable
//!
//! Diffraction from a circular aperture has closed-form answers, all of them. The
//! Airy pattern is `[2J₁(v)/v]²`, its first zero is at `v = 3.8317`, the energy
//! inside that zero is `1 - J₀² - J₁²` which is 83.8%, and the incoherent MTF of a
//! perfect circular pupil is `(2/π)(arccos s - s√(1-s²))`. Every claim below is
//! checked against one of those rather than against another implementation — the
//! same reason Planck's law is checked against Wien's displacement law and the
//! symplectic integrator against a harmonic oscillator's energy.
//!
//! # Conventions, stated because they differ between textbooks
//!
//! `NA` is the numerical aperture on the side being asked about, `n sin θ`. The
//! radial coordinate `v` is `2π·NA·r/λ` with `r` a distance in the image. Spatial
//! frequencies are cycles per unit length in the image, and the cutoff is the
//! *incoherent* one, `2NA/λ`. A coherent system cuts off at half that, and quoting
//! the wrong one is a factor-of-two error in every resolution claim.

use dualis_units::Length;

/// Bessel function of the first kind, order zero.
///
/// Abramowitz and Stegun 9.4.1 and 9.4.3 — polynomial below 3 and an asymptotic
/// form above it, accurate to about 1e-7. Written out rather than pulled in as a
/// dependency, since these two functions are all this crate needs of the family and
/// the accuracy is far past what any optical measurement reaches.
pub fn bessel_j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 3.0 {
        let y = (x / 3.0) * (x / 3.0);
        1.0 + y
            * (-2.249_999_7
                + y * (1.265_620_8
                    + y * (-0.316_386_6
                        + y * (0.044_447_9 + y * (-0.003_944_4 + y * 0.000_210_0)))))
    } else {
        let t = 3.0 / ax;
        let f = 0.797_884_56
            + t * (-0.000_000_77
                + t * (-0.005_527_40
                    + t * (-0.000_095_12
                        + t * (0.001_372_37 + t * (-0.000_728_05 + t * 0.000_144_76)))));
        // A&S quotes this leading term as 0.78539816, which is pi/4 to the precision
        // it gives; the constant is exact, so it is used instead.
        let theta = ax - std::f64::consts::FRAC_PI_4
            + t * (-0.041_663_97
                + t * (-0.000_039_54
                    + t * (0.002_625_73
                        + t * (-0.000_541_25 + t * (-0.000_293_33 + t * 0.000_135_58)))));
        f * theta.cos() / ax.sqrt()
    }
}

/// Bessel function of the first kind, order one.
///
/// Abramowitz and Stegun 9.4.4 and 9.4.6. Odd, so `J₁(-x) = -J₁(x)`.
pub fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    let value = if ax < 3.0 {
        let y = (x / 3.0) * (x / 3.0);
        ax * (0.5
            + y * (-0.562_499_85
                + y * (0.210_935_73
                    + y * (-0.039_542_89
                        + y * (0.004_433_19 + y * (-0.000_317_61 + y * 0.000_011_09))))))
    } else {
        let t = 3.0 / ax;
        let f = 0.797_884_56
            + t * (0.000_001_56
                + t * (0.016_596_67
                    + t * (0.000_171_05
                        + t * (-0.002_495_11 + t * (0.001_136_53 + t * (-0.000_200_33))))));
        // A&S's 2.35619449, which is 3*pi/4.
        let theta = ax - 3.0 * std::f64::consts::FRAC_PI_4
            + t * (0.124_996_12
                + t * (0.000_056_50
                    + t * (-0.006_378_79
                        + t * (0.000_743_48 + t * (0.000_798_24 + t * (-0.000_291_66))))));
        f * theta.cos() / ax.sqrt()
    };
    if x < 0.0 {
        -value
    } else {
        value
    }
}

/// First zero of `J₁`, which is where the Airy pattern's first dark ring sits.
pub const FIRST_AIRY_ZERO: f64 = 3.831_705_970_207_512;
/// Second zero of `J₁` — the second dark ring.
pub const SECOND_AIRY_ZERO: f64 = 7.015_586_669_815_619;

/// The Airy pattern's intensity at radial coordinate `v = 2π·NA·r/λ`, normalised to
/// 1 at the centre.
///
/// `[2J₁(v)/v]²`. The limit at the centre is exactly 1, taken analytically rather
/// than left to divide by zero.
pub fn airy_intensity(v: f64) -> f64 {
    if v.abs() < 1e-8 {
        // 2J1(v)/v -> 1 as v -> 0, and the series correction is O(v^2).
        return 1.0 - v * v / 4.0;
    }
    let a = 2.0 * bessel_j1(v) / v;
    a * a
}

/// A single slit's far-field intensity, normalised to 1 on the axis.
///
/// ```text
///   I(θ) = sinc²(π a sinθ / λ),      sinc x = sin x / x
/// ```
///
/// The one-dimensional counterpart of [`airy_intensity`], and worth having beside it because the two
/// differ in the way that matters: a slit's first zero is at `sinθ = λ/a` **exactly**, while a
/// circular aperture's is at `1.22 λ/D`. The 1.22 is the first zero of `J₁` and has no closed form;
/// the 1 here is `sin π = 0`.
///
/// # What this is an approximation to
///
/// Scalar diffraction: it assumes the field in the aperture is the incident field, which is
/// Kirchhoff's boundary condition and is not what Maxwell's equations give. The screen has thickness,
/// the metal carries currents, and the field in the opening is perturbed near its edges over a
/// distance of order `λ`. So the fraction of the aperture that is wrong goes as `λ/a`, and this
/// formula is exact only as `a/λ → ∞`.
///
/// `crates/dualis/tests/a_slit.rs` measures that convergence against a field solution, because a
/// closed form's regime of validity is worth as much as the closed form. Measured, as the largest
/// absolute difference in intensity: 0.0057 at `a = 12λ`, 0.0311 at `3λ`, and **0.2772 at `1λ`** —
/// where it is wrong by more than a quarter of the pattern it is predicting.
pub fn single_slit_intensity(width: Length, wavelength: Length, sin_theta: f64) -> f64 {
    let u = std::f64::consts::PI * width.to_si() * sin_theta / wavelength.to_si();
    if u.abs() < 1e-12 {
        return 1.0;
    }
    (u.sin() / u).powi(2)
}

/// Where a single slit's `m`-th dark fringe falls, as `sinθ = mλ/a`.
///
/// Returns `None` past grazing, which happens for `m λ > a` — a slit narrower than the wavelength
/// has **no** zeros at all, and that is not an edge case to be clamped away: it is the statement
/// that a sub-wavelength opening does not diffract into a pattern, it radiates.
pub fn slit_zero(width: Length, wavelength: Length, order: u32) -> Option<f64> {
    let s = order as f64 * wavelength.to_si() / width.to_si();
    (order > 0 && s <= 1.0).then_some(s)
}

/// Fraction of the total energy inside radius `v`.
///
/// `1 - J₀²(v) - J₁²(v)`, exactly. At the first dark ring this is 0.8378: a sixth of
/// a perfect image's light is *outside* the Airy disc, spread over rings, and that
/// is not aberration — it is what a finite aperture does.
pub fn encircled_energy(v: f64) -> f64 {
    let (j0, j1) = (bessel_j0(v), bessel_j1(v));
    (1.0 - j0 * j0 - j1 * j1).clamp(0.0, 1.0)
}

/// Radius of the Airy disc — the first dark ring — for a wavelength and numerical
/// aperture.
///
/// `0.61 λ/NA`, which is also the Rayleigh resolution criterion: two point sources
/// this far apart put each one's peak on the other's first zero, and that is the
/// conventional line between resolved and not.
pub fn airy_radius(wavelength: Length, na: f64) -> Length {
    if na <= 0.0 {
        return Length::from_si(f64::INFINITY);
    }
    Length::from_si(FIRST_AIRY_ZERO * wavelength.to_si() / (std::f64::consts::TAU * na))
}

/// The Rayleigh criterion, which is [`airy_radius`] under its other name.
pub fn rayleigh_limit(wavelength: Length, na: f64) -> Length {
    airy_radius(wavelength, na)
}

/// The Abbe limit, `λ/(2 NA)`.
///
/// A different question from Rayleigh's, and a smaller number: Abbe asks what
/// grating period the system can still transmit at all, which is the reciprocal of
/// the incoherent cutoff frequency. Rayleigh asks when two points look like two.
/// Quoting one and meaning the other is a 22% error.
pub fn abbe_limit(wavelength: Length, na: f64) -> Length {
    if na <= 0.0 {
        return Length::from_si(f64::INFINITY);
    }
    Length::from_si(wavelength.to_si() / (2.0 * na))
}

/// Incoherent cutoff frequency, `2NA/λ`, in cycles per metre.
///
/// Past this the modulation transfer is exactly zero: the system does not attenuate
/// finer detail, it does not transmit it at all. There is nothing to deconvolve
/// back.
pub fn cutoff_frequency(wavelength: Length, na: f64) -> f64 {
    if wavelength.to_si() <= 0.0 {
        return 0.0;
    }
    2.0 * na / wavelength.to_si()
}

/// Modulation transfer of a perfect circular pupil at a fraction `s` of the cutoff.
///
/// `(2/π)(arccos s - s√(1-s²))`, the autocorrelation of a disc. This is the ceiling:
/// no lens does better, every real one does worse, and the difference between the
/// two is what a lens is judged on.
pub fn mtf_ideal(s: f64) -> f64 {
    if s <= 0.0 {
        return 1.0;
    }
    if s >= 1.0 {
        return 0.0;
    }
    (2.0 / std::f64::consts::PI) * (s.acos() - s * (1.0 - s * s).sqrt())
}

/// Modulation transfer at a spatial frequency in cycles per metre.
pub fn mtf_at(frequency_per_m: f64, wavelength: Length, na: f64) -> f64 {
    let cutoff = cutoff_frequency(wavelength, na);
    if cutoff <= 0.0 {
        return 0.0;
    }
    mtf_ideal(frequency_per_m / cutoff)
}

/// Depth of focus by the Rayleigh quarter-wave convention: `n λ / NA²`, the full
/// depth over which the wavefront error stays under a quarter wave.
///
/// Conventions for this differ by factors of two between textbooks and between
/// vendors, so the one in use is named here rather than assumed. The scaling is the
/// part that is not a convention: **depth goes as one over NA squared**, so
/// doubling the aperture to halve the spot quarters the depth, and that trade is
/// why a high-NA objective has no working range.
pub fn depth_of_focus(wavelength: Length, na: f64, refractive_index: f64) -> Length {
    if na <= 0.0 {
        return Length::from_si(f64::INFINITY);
    }
    Length::from_si(refractive_index * wavelength.to_si() / (na * na))
}

/// Strehl ratio from an RMS wavefront error in waves, by the Maréchal
/// approximation `exp(-(2π σ)²)`.
///
/// The number that says whether a system is diffraction limited: 0.8 is the
/// conventional threshold, and it corresponds to λ/14 RMS — considerably tighter
/// than the λ/4 peak-to-valley figure the same tradition also quotes. Valid while
/// the error is small; past about half a wave it understates the damage.
pub fn strehl_from_wavefront_error(rms_waves: f64) -> f64 {
    let phase = std::f64::consts::TAU * rms_waves;
    (-phase * phase).exp()
}

/// Fresnel number `a²/(λ z)` for an aperture of radius `a` seen from distance `z`.
///
/// Which regime a propagation is in: much greater than one is the near field, where
/// the shadow of an aperture still looks like the aperture; around or below one is
/// the far field, where it looks like a diffraction pattern instead. Reaching for
/// [`airy_intensity`] in the near field is a category error, and this is how to tell.
pub fn fresnel_number(aperture_radius: Length, distance: Length, wavelength: Length) -> f64 {
    let (a, z, l) = (
        aperture_radius.to_si(),
        distance.to_si(),
        wavelength.to_si(),
    );
    if z <= 0.0 || l <= 0.0 {
        return f64::INFINITY;
    }
    a * a / (l * z)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nm(v: f64) -> Length {
        Length::nm(v)
    }

    /// The Bessel approximations, against the tabulated zeros and values. These are
    /// the foundation everything else in the module stands on, so they are checked
    /// first and against numbers that come from outside the code.
    #[test]
    fn the_bessel_functions_hit_their_tabulated_zeros() {
        // J1's zeros, which set the Airy rings.
        for zero in [
            FIRST_AIRY_ZERO,
            SECOND_AIRY_ZERO,
            10.173_468_135_062_722,
            13.323_691_936_314_223,
        ] {
            assert!(
                bessel_j1(zero).abs() < 2e-7,
                "J1({zero}) should vanish, got {}",
                bessel_j1(zero)
            );
        }
        // J0's zeros.
        for zero in [
            2.404_825_557_695_773,
            5.520_078_110_286_311,
            8.653_727_912_911_013,
        ] {
            assert!(
                bessel_j0(zero).abs() < 2e-7,
                "J0({zero}) should vanish, got {}",
                bessel_j0(zero)
            );
        }
        // Values at the origin and a couple of tabulated points.
        assert!((bessel_j0(0.0) - 1.0).abs() < 1e-7);
        assert!(bessel_j1(0.0).abs() < 1e-7);
        assert!((bessel_j0(1.0) - 0.765_197_686_557_966).abs() < 1e-7);
        assert!((bessel_j1(1.0) - 0.440_050_585_744_933).abs() < 1e-7);
        assert!((bessel_j0(5.0) - (-0.177_596_771_314_338)).abs() < 1e-7);
        assert!((bessel_j1(5.0) - (-0.327_579_137_591_465)).abs() < 1e-7);
        // J1 is odd, J0 even.
        assert!((bessel_j1(-2.5) + bessel_j1(2.5)).abs() < 1e-12);
        assert!((bessel_j0(-2.5) - bessel_j0(2.5)).abs() < 1e-12);
    }

    /// The Airy pattern is 1 at the centre, falls monotonically to zero at the first
    /// ring, and the ring is genuinely dark.
    #[test]
    fn the_airy_pattern_has_a_bright_centre_and_dark_rings() {
        assert!((airy_intensity(0.0) - 1.0).abs() < 1e-12);
        assert!(airy_intensity(FIRST_AIRY_ZERO) < 1e-13);
        assert!(airy_intensity(SECOND_AIRY_ZERO) < 1e-13);

        // Monotone decreasing out to the first zero.
        let mut previous = 1.0;
        for i in 1..=200 {
            let v = FIRST_AIRY_ZERO * i as f64 / 200.0;
            let value = airy_intensity(v);
            assert!(value <= previous + 1e-12, "rose at v = {v}");
            previous = value;
        }
        // Half maximum at v = 1.6163, which is the FWHM of the Airy core.
        let half = airy_intensity(1.616_339);
        assert!((half - 0.5).abs() < 1e-4, "got {half}");
        // The first ring's peak is 1.75% of the centre — faint, and the reason a
        // bright star shows rings while a dim one does not.
        let ring = (0..1000)
            .map(|i| {
                let v = FIRST_AIRY_ZERO + (SECOND_AIRY_ZERO - FIRST_AIRY_ZERO) * i as f64 / 1000.0;
                airy_intensity(v)
            })
            .fold(0.0f64, f64::max);
        assert!((ring - 0.0175).abs() < 0.0005, "first ring peak {ring}");
    }

    /// 83.8% of a perfect image's light is inside the Airy disc and the rest is not.
    /// That number is `1 - J₀²(3.8317)`, a closed form, and it is why an encircled
    /// energy spec of 90% needs a radius well past the first dark ring.
    #[test]
    fn the_airy_disc_holds_five_sixths_of_the_light() {
        assert_eq!(encircled_energy(0.0), 0.0);
        let inside = encircled_energy(FIRST_AIRY_ZERO);
        assert!((inside - 0.837_8).abs() < 1e-4, "got {inside}");
        // By the second ring, 91%.
        let two_rings = encircled_energy(SECOND_AIRY_ZERO);
        assert!((two_rings - 0.910_3).abs() < 1e-3, "got {two_rings}");
        // Monotone, and never past 1.
        let mut previous = 0.0;
        for i in 0..=500 {
            let v = 30.0 * i as f64 / 500.0;
            let e = encircled_energy(v);
            assert!(e >= previous - 1e-6, "fell at v = {v}");
            assert!(e <= 1.0);
            previous = e;
        }
        assert!(
            previous > 0.97,
            "most of it is in eventually, got {previous}"
        );
    }

    /// The resolution limits, against the figures a microscope is sold on: a 1.4 NA
    /// oil objective at 550 nm resolves about 240 nm by Rayleigh and 196 nm by Abbe.
    /// The two differ by 22%, and confusing them is the most common way a resolution
    /// claim is wrong.
    #[test]
    fn resolution_limits_match_the_catalogue_figures() {
        let green = nm(550.0);
        let rayleigh = rayleigh_limit(green, 1.4);
        let abbe = abbe_limit(green, 1.4);
        assert!(
            (rayleigh.in_nm() - 239.6).abs() < 1.0,
            "{} nm",
            rayleigh.in_nm()
        );
        assert!((abbe.in_nm() - 196.4).abs() < 1.0, "{} nm", abbe.in_nm());
        assert!(rayleigh > abbe);
        assert!((rayleigh / abbe - 1.22).abs() < 0.01);

        // Both scale as one over NA: halving the aperture doubles the spot.
        assert!((rayleigh_limit(green, 0.7) / rayleigh_limit(green, 1.4) - 2.0).abs() < 1e-12);
        // And linearly in wavelength, which is why blue light resolves better.
        assert!(rayleigh_limit(nm(405.0), 1.4) < rayleigh_limit(nm(640.0), 1.4));
        // A zero aperture resolves nothing, and says so rather than dividing by zero.
        assert!(!rayleigh_limit(green, 0.0).is_finite());
    }

    /// This is what the geometric half of the crate cannot say. A ray tracer's spot
    /// for a perfect lens is a point; the real one is 240 nm across, and no surface
    /// figure improves it.
    #[test]
    fn a_perfect_lens_still_has_a_finite_spot() {
        let spot = airy_radius(nm(550.0), 1.4);
        assert!(spot.in_nm() > 200.0, "a perfect lens is not a point");
        // Even a perfectly figured surface, with zero wavefront error and therefore
        // a Strehl ratio of exactly 1, has this spot.
        assert_eq!(strehl_from_wavefront_error(0.0), 1.0);
    }

    /// The MTF of a perfect circular pupil, against its closed form at the points
    /// where the closed form is simple: 1 at zero frequency, 0 at cutoff, and
    /// 0.3910 at half cutoff.
    #[test]
    fn the_ideal_mtf_matches_the_disc_autocorrelation() {
        assert!((mtf_ideal(0.0) - 1.0).abs() < 1e-12);
        assert_eq!(mtf_ideal(1.0), 0.0);
        assert_eq!(mtf_ideal(1.5), 0.0, "past cutoff is zero, not negative");
        assert_eq!(mtf_ideal(-0.2), 1.0);
        // (2/pi)(acos(0.5) - sqrt(3)/4) = 0.391002 — the figure quoted as "0.391 at
        // half cutoff" in every diffraction-limited MTF chart.
        assert!(
            (mtf_ideal(0.5) - 0.391_002).abs() < 1e-5,
            "{}",
            mtf_ideal(0.5)
        );
        // Monotone decreasing.
        let mut previous = 1.0;
        for i in 0..=100 {
            let value = mtf_ideal(i as f64 / 100.0);
            assert!(value <= previous + 1e-12);
            previous = value;
        }
    }

    /// The cutoff in the units a lens is specified in: a 1.4 NA objective at 550 nm
    /// transmits up to 5091 cycles per millimetre, and nothing finer. Detail past
    /// that is not attenuated, it is absent.
    #[test]
    fn the_cutoff_is_where_transmission_stops_entirely() {
        let green = nm(550.0);
        let cutoff_per_mm = cutoff_frequency(green, 1.4) / 1000.0;
        assert!(
            (cutoff_per_mm - 5090.9).abs() < 1.0,
            "{cutoff_per_mm} cyc/mm"
        );
        // The cutoff period is exactly the Abbe limit, which is what makes them the
        // same statement.
        let period = Length::from_si(1.0 / cutoff_frequency(green, 1.4));
        assert!((period - abbe_limit(green, 1.4)).abs().in_nm() < 1e-6);

        assert_eq!(mtf_at(cutoff_frequency(green, 1.4) * 1.01, green, 1.4), 0.0);
        assert!(mtf_at(cutoff_frequency(green, 1.4) * 0.99, green, 1.4) > 0.0);
        // A camera lens, in the units its charts use: f/2.8 is about NA 0.18, so it
        // cuts off near 650 cyc/mm — and its charts stop at 40, because that is
        // where a sensor's pixels give up long before the optics do.
        let camera = cutoff_frequency(green, 0.18) / 1000.0;
        assert!(camera > 600.0 && camera < 700.0, "{camera} cyc/mm");
        assert!(
            mtf_at(40e3, green, 0.18) > 0.9,
            "a good lens is flat at 40 cyc/mm"
        );
    }

    /// The trade that decides every objective's design: depth of focus goes as one
    /// over NA squared, so doubling the aperture to halve the spot quarters the
    /// depth. A 1.4 NA lens has about 280 nm of it.
    #[test]
    fn depth_of_focus_falls_as_the_square_of_the_aperture() {
        let green = nm(550.0);
        // A dry 0.7 NA lens against an oil 1.4 NA one, as they are really used. The
        // aperture alone would give a factor of four; immersion in n = 1.515 hands
        // back 1.515 of it, so the real ratio is 2.64. Both effects are in the
        // formula and neither is optional.
        let deep = depth_of_focus(green, 0.7, 1.0);
        let shallow = depth_of_focus(green, 1.4, 1.515);
        assert!(
            (deep / shallow - 2.640).abs() < 0.01,
            "ratio {}",
            deep / shallow
        );
        assert!(
            (deep.in_nm() - 1122.4).abs() < 5.0,
            "0.7 NA dry: {} nm",
            deep.in_nm()
        );
        assert!(
            (shallow.in_nm() - 425.2).abs() < 5.0,
            "1.4 NA in oil: {} nm",
            shallow.in_nm()
        );
        // Halving the spot costs three quarters of the depth.
        let spot_ratio = airy_radius(green, 1.4) / airy_radius(green, 0.7);
        let depth_ratio = depth_of_focus(green, 1.4, 1.0) / depth_of_focus(green, 0.7, 1.0);
        assert!((spot_ratio - 0.5).abs() < 1e-12);
        assert!((depth_ratio - 0.25).abs() < 1e-12);
    }

    /// The diffraction-limited threshold, and two discrepancies inside the tradition
    /// that states it.
    ///
    /// The customary "lambda/14 RMS" gives a Strehl of 0.818, comfortably *past* the
    /// 0.8 threshold rather than at it — the exact threshold is lambda/13.3. And the
    /// same tradition's "quarter wave" figure is peak-to-valley, not RMS: a quarter
    /// wave RMS is a Strehl of 0.085, a thoroughly broken system.
    #[test]
    fn the_diffraction_limit_is_a_fourteenth_of_a_wave_rms() {
        let at_fourteenth = strehl_from_wavefront_error(1.0 / 14.0);
        assert!((at_fourteenth - 0.8176).abs() < 1e-3, "got {at_fourteenth}");
        assert!(at_fourteenth > 0.8, "lambda/14 does meet the criterion");
        // Inverting the Marechal approximation for exactly 0.8.
        let exact_threshold = (-0.8f64.ln()).sqrt() / std::f64::consts::TAU;
        assert!(
            (1.0 / exact_threshold - 13.30).abs() < 0.02,
            "the 0.8 threshold is lambda/{}",
            1.0 / exact_threshold
        );
        assert!((strehl_from_wavefront_error(exact_threshold) - 0.8).abs() < 1e-12);
        // A quarter wave RMS is a badly aberrated system, not a limit case.
        assert!(strehl_from_wavefront_error(0.25) < 0.1);
        // Monotone, bounded, and 1 only for a perfect wavefront.
        assert_eq!(strehl_from_wavefront_error(0.0), 1.0);
        let mut previous = 1.0;
        for i in 1..=50 {
            let s = strehl_from_wavefront_error(i as f64 / 100.0);
            assert!(s < previous && s > 0.0);
            previous = s;
        }
    }

    /// Which regime a propagation is in, so that a far-field formula is not applied
    /// in the near field. A 5 mm aperture at 100 mm in green light has a Fresnel
    /// number of 455: deeply near field, where the beam is still a beam and no Airy
    /// pattern has formed.
    ///
    /// The distance it takes to reach the far field is `a²/λ`, and for this aperture
    /// that is 45 metres — which is why an Airy pattern is something you see at a
    /// lens's focus, where the lens has done the propagating, and not down a bench.
    #[test]
    fn the_fresnel_number_separates_the_regimes() {
        let green = nm(550.0);
        let near = fresnel_number(Length::mm(5.0), Length::mm(100.0), green);
        assert!((near - 454.5).abs() < 0.5, "got {near}");
        assert!(near > 1.0, "near field");
        // Ten metres is still not far enough: 4.5.
        let ten_m = fresnel_number(Length::mm(5.0), Length::m(10.0), green);
        assert!(ten_m > 1.0, "still near field at 10 m, got {ten_m}");
        // A hundred metres is.
        let far = fresnel_number(Length::mm(5.0), Length::m(100.0), green);
        assert!(far < 1.0, "far field, got {far}");
        assert!((far - 0.4545).abs() < 0.001, "got {far}");
        // Zero distance is not a regime, and reports infinity rather than a NaN.
        assert!(!fresnel_number(Length::mm(5.0), Length::ZERO, green).is_finite());
    }

    /// **A slit's zeros are at `mλ/a` exactly**, which is the one place a diffraction pattern has a
    /// closed form with no special constant in it.
    #[test]
    fn a_slits_zeros_are_at_m_lambda_over_a() {
        let (a, l) = (Length::um(10.0), Length::nm(500.0));
        for m in 1..=5u32 {
            let s = slit_zero(a, l, m).expect("a 20-wavelength slit has five zeros");
            assert!(
                (s - m as f64 * 0.05).abs() < 1e-15,
                "the {m}th zero is at {m} lambda / a: {s}"
            );
            assert!(
                single_slit_intensity(a, l, s) < 1e-24,
                "and the intensity there is zero: {:.3e}",
                single_slit_intensity(a, l, s)
            );
        }
        // Halfway between the first two zeros, `sinc²(3π/2) = 1/(2.25 π²)` — 4.503%. Not quite the
        // sidelobe's *peak*, which sits slightly inside that angle at 4.72%: `sinc` is falling as
        // `1/u` while it oscillates, so every maximum is a little before the midpoint.
        let between = 1.5 * 0.05;
        let side = single_slit_intensity(a, l, between);
        let closed = 1.0 / (2.25 * std::f64::consts::PI.powi(2));
        assert!(
            (side - closed).abs() < 1e-12,
            "sinc^2(3 pi / 2) = 1/(2.25 pi^2) = {:.4}%: measured {:.4}%",
            closed * 100.0,
            side * 100.0
        );
        // On the axis it is exactly one, taken analytically rather than left to divide by zero.
        assert_eq!(single_slit_intensity(a, l, 0.0), 1.0);
    }

    /// **A sub-wavelength slit has no zeros at all**, and that is a statement rather than an edge
    /// case.
    ///
    /// `sinθ = λ/a > 1` has no solution: the first dark fringe would be past grazing. Such an opening
    /// does not diffract into a pattern, it radiates — and the scalar formula, which happily returns
    /// a number for any angle, is furthest from the truth exactly there.
    #[test]
    fn a_slit_narrower_than_the_wavelength_has_no_zeros() {
        let l = Length::nm(500.0);
        assert!(slit_zero(Length::nm(400.0), l, 1).is_none());
        assert!(
            slit_zero(Length::nm(500.0), l, 1).is_some(),
            "at a = lambda the zero is at grazing"
        );
        assert!(slit_zero(Length::nm(900.0), l, 2).is_none());
        assert!(
            slit_zero(Length::um(10.0), l, 0).is_none(),
            "there is no zeroth zero"
        );
        // And a narrow slit's pattern is nearly flat, which is what "radiates" means. Measured
        // against a wide one at the same angle rather than against a round number: the contrast is
        // the statement, and either figure alone is just a number.
        let narrow = single_slit_intensity(Length::nm(300.0), l, 1.0);
        let wide = single_slit_intensity(Length::um(10.0), l, 0.53);
        println!("  0.6 lambda at grazing {narrow:.4}, 20 lambda off axis {wide:.6}");
        assert!(
            (0.20..0.30).contains(&narrow),
            "a 0.6-wavelength slit still sends a quarter of its light to grazing: {narrow:.4}"
        );
        assert!(
            narrow / wide > 100.0,
            "which is two orders more than a wide slit sends anywhere off axis: {:.0}x",
            narrow / wide
        );
    }
}
