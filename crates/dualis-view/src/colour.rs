//! What colour a thing at a temperature actually is.
//!
//! Every other colour in this crate is **chosen** — `filmstrip`'s ramp is a designer's blue-to-red
//! that says "more" and "less", and for a pressure swinging about zero that is the right answer,
//! because a sound wave has no colour. A body at 1200 K does. Its spectral exitance is Planck's
//! law, the eye's response to that spectrum is the CIE 1931 colour matching functions, and the
//! result is a number rather than a taste: **1200 K is that particular orange and no other.**
//!
//! That distinction is why this module exists in a library that otherwise refuses to build a
//! renderer. `ARCHITECTURE.md` records the decision not to compete on rendering — RTX, materials,
//! shadows and USD belong to tools that are enormous and good at them, and the move is to be
//! reachable from those rather than to rebuild them. Colour computed from a temperature is the
//! part they *cannot* do: in Blender the colour of a hot block is picked by a person, and here it
//! is a consequence of the physics the simulation already computed.
//!
//! # What this does and does not model
//!
//! It models **thermal emission**, and nothing else. It does not model reflectance, because
//! nothing in this workspace holds a visible reflectance: `Substance` carries an emissivity — a
//! broadband infrared number a thermal balance needs — and no colour. So a body cool enough not
//! to glow has no computed appearance here at all, and [`glow_fraction`] is what says which
//! regime a temperature is in rather than a threshold invented in a renderer.
//!
//! # The approximations, and their stated error
//!
//! Planck's law and the sRGB transfer are exact. The CIE 1931 colour matching functions are the
//! **multi-lobe analytic fit** of Wyman, Sloan and Shirley (JCGT 2013) rather than the tabulated
//! standard: seven Gaussians against 471 rows of table, agreeing with the tables to about one
//! part in a hundred of the chromaticity. Every tolerance below traces to that figure, and the
//! test that pins it is [`the_planckian_locus_passes_through_illuminant_a`](self), which checks a
//! **published** coordinate that this code had no hand in producing.

/// Planck's constant, J·s. Exact by the 2019 SI redefinition.
const H: f64 = 6.626_070_15e-34;
/// The speed of light in vacuum, m/s. Exact.
const C: f64 = 299_792_458.0;
/// The Boltzmann constant, J/K. Exact by the 2019 SI redefinition.
const K_B: f64 = 1.380_649e-23;
/// The Stefan–Boltzmann constant, W·m⁻²·K⁻⁴, which follows from the three above.
const SIGMA: f64 = 5.670_374_419e-8;

/// The visible band this module integrates over, in metres. The CIE tables' own range.
const VISIBLE: (f64, f64) = (360e-9, 830e-9);

/// Spectral radiant exitance of a black body, in W·m⁻³ — power per unit area per unit wavelength.
///
/// `M(λ, T) = (2πhc²/λ⁵) / (exp(hc/λk_BT) − 1)`, the per-wavelength form, integrating to `σT⁴`
/// over all wavelengths — which is [`the_planck_curve_integrates_to_stefan_boltzmann`](self).
///
/// Returns zero rather than a NaN or an infinity for a non-positive wavelength or temperature: a
/// body at absolute zero radiates nothing, and that is a value rather than an error.
pub fn planck_exitance(wavelength_m: f64, temperature_k: f64) -> f64 {
    if wavelength_m <= 0.0 || temperature_k <= 0.0 {
        return 0.0;
    }
    let x = H * C / (wavelength_m * K_B * temperature_k);
    // `exp_m1` rather than `exp() - 1`: at long wavelengths x is small and the subtraction loses
    // every digit it has. At 10 µm and 300 K, x is 4.8 and it hardly matters; at 1 mm it does.
    let denom = x.exp_m1();
    if !denom.is_finite() || denom <= 0.0 {
        // x large: the exponential overflowed, and the true value underflows to zero.
        return 0.0;
    }
    let l5 = wavelength_m.powi(5);
    2.0 * std::f64::consts::PI * H * C * C / (l5 * denom)
}

/// One lobe of the analytic colour-matching fit: a Gaussian with a different width on each side.
fn lobe(nm: f64, mu: f64, sigma_low: f64, sigma_high: f64) -> f64 {
    let sigma = if nm < mu { sigma_low } else { sigma_high };
    let t = (nm - mu) / sigma;
    (-0.5 * t * t).exp()
}

/// The CIE 1931 2° colour matching functions `x̄, ȳ, z̄` at a wavelength, in nanometres.
///
/// The multi-lobe analytic fit named in the module documentation. `ȳ` is the photopic luminous
/// efficiency curve, which is why [`glow_fraction`]'s luminous form and this share an
/// implementation rather than carrying two spellings of the same curve.
fn cie_xyz_bar(nm: f64) -> (f64, f64, f64) {
    let x = 1.056 * lobe(nm, 599.8, 37.9, 31.0) + 0.362 * lobe(nm, 442.0, 16.0, 26.7)
        - 0.065 * lobe(nm, 501.1, 20.4, 26.2);
    let y = 0.821 * lobe(nm, 568.8, 46.9, 40.5) + 0.286 * lobe(nm, 530.9, 16.3, 31.1);
    let z = 1.217 * lobe(nm, 437.0, 11.8, 36.0) + 0.681 * lobe(nm, 459.0, 26.0, 13.8);
    (x, y, z)
}

/// Integrate a spectrum against the colour matching functions, giving unnormalised CIE XYZ.
///
/// One nanometre steps across the CIE range, which is the tables' own resolution — finer buys
/// nothing against a fit whose error is a hundred times larger.
fn xyz_of(spectrum: impl Fn(f64) -> f64) -> (f64, f64, f64) {
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    let (lo, hi) = VISIBLE;
    let steps = ((hi - lo) * 1e9).round() as usize;
    for i in 0..=steps {
        let nm = lo * 1e9 + i as f64;
        let power = spectrum(nm * 1e-9);
        let (xb, yb, zb) = cie_xyz_bar(nm);
        x += power * xb;
        y += power * yb;
        z += power * zb;
    }
    (x, y, z)
}

/// The CIE 1931 chromaticity `(x, y)` of a black body at this temperature — its point on the
/// Planckian locus.
///
/// Chromaticity and not brightness, because those are different questions with different answers:
/// what colour a body at `T` is, is physics; how bright it appears is exposure, which belongs to
/// whoever is drawing. [`glow_fraction`] is the physical half of the second question.
pub fn planckian_chromaticity(temperature_k: f64) -> (f64, f64) {
    let (x, y, z) = xyz_of(|l| planck_exitance(l, temperature_k));
    let sum = x + y + z;
    if sum <= 0.0 {
        // A body cold enough that its visible exitance underflows to nothing has no colour, and
        // the equal-energy point is returned rather than a NaN so arithmetic downstream stays
        // arithmetic. **It is not a colour to draw** — `blackbody_srgb` checks the emission
        // itself rather than trusting this pair, because a body that emits nothing came out
        // *white* when it did trust it, which is a plausible-looking answer to a question with
        // no answer. `glow_fraction` is the caller's version of the same check.
        return (1.0 / 3.0, 1.0 / 3.0);
    }
    (x / sum, y / sum)
}

/// The sRGB triple a black body at this temperature glows, at full brightness.
///
/// **Chromaticity, scaled so the brightest channel is 255.** A caller drawing a picture chooses
/// exposure — the sun and a candle flame are the same two colours whatever the shutter — and
/// [`glow_fraction`] is the physical number to scale by when the picture should say which is
/// which.
///
/// Out-of-gamut chromaticities are clamped per channel, which sRGB requires for any real
/// radiator: a 1000 K body is redder than sRGB's red primary, so its blue channel comes out
/// negative and is clamped to zero. That clamp changes the colour — it has to, because the
/// monitor cannot make that colour — and it is stated here rather than hidden in a saturate.
pub fn blackbody_srgb(temperature_k: f64) -> [u8; 3] {
    // The emission is checked here rather than the chromaticity, and the difference is a defect
    // this function had: a body at absolute zero has no chromaticity, the fallback pair is the
    // equal-energy point, and trusting it rendered a body that emits nothing as **white**.
    let (raw_x, raw_y, raw_z) = xyz_of(|l| planck_exitance(l, temperature_k));
    if raw_x + raw_y + raw_z <= 0.0 {
        return [0, 0, 0];
    }
    let (x, y) = planckian_chromaticity(temperature_k);
    if y <= 0.0 {
        return [0, 0, 0];
    }
    // Chromaticity back to XYZ at unit luminance.
    let (big_x, big_y, big_z) = (x / y, 1.0, (1.0 - x - y) / y);

    // The sRGB primaries with the D65 white point, the standard matrix.
    let r = 3.240_454_2 * big_x - 1.537_138_5 * big_y - 0.498_531_4 * big_z;
    let g = -0.969_266_0 * big_x + 1.876_010_8 * big_y + 0.041_556_0 * big_z;
    let b = 0.055_643_4 * big_x - 0.204_025_9 * big_y + 1.057_225_2 * big_z;

    let linear = [r.max(0.0), g.max(0.0), b.max(0.0)];
    let peak = linear.iter().fold(0.0f64, |m, v| m.max(*v));
    if peak <= 0.0 {
        return [0, 0, 0];
    }
    let mut out = [0u8; 3];
    for (i, v) in linear.iter().enumerate() {
        // The sRGB transfer function, exactly as the standard states it.
        let c = v / peak;
        let encoded = if c <= 0.003_130_8 {
            12.92 * c
        } else {
            1.055 * c.powf(1.0 / 2.4) - 0.055
        };
        out[i] = (encoded * 255.0).round().clamp(0.0, 255.0) as u8;
    }
    out
}

/// The fraction of a black body's total radiated power that falls in the visible band.
///
/// `∫₃₆₀^₈₃₀ M_λ dλ / σT⁴`, and it is the honest answer to "does this glow?" — a question a
/// renderer otherwise answers with a threshold somebody picked. Room temperature returns about
/// `1e-30`: an object at 300 K emits no visible light at all, which is why a cool thing in a dark
/// room is invisible rather than dim. It climbs steeply — the Wien tail, `exp(−hc/λk_BT)` — and
/// at 3000 K, a tungsten filament, it is a few percent.
///
/// A caller drawing a scene multiplies its glow by this and gets the right answer for free: a
/// warm block does not glow, a melting one does, and nothing had to decide where the line is.
pub fn glow_fraction(temperature_k: f64) -> f64 {
    if temperature_k <= 0.0 {
        return 0.0;
    }
    let (lo, hi) = VISIBLE;
    let steps = ((hi - lo) * 1e9).round() as usize;
    let dl = 1e-9;
    let mut visible = 0.0;
    for i in 0..=steps {
        visible += planck_exitance(lo + i as f64 * dl, temperature_k) * dl;
    }
    visible / (SIGMA * temperature_k.powi(4))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **Wien's displacement law**, which fixes where the Planck curve peaks: `λ_max·T = b` with
    /// `b = 2.897771955e-3 m·K`. Located by scanning at 1 pm resolution, so the tolerance is that
    /// grid against the peak's width, not a fudge.
    #[test]
    fn the_planck_peak_obeys_wiens_displacement_law() {
        const B: f64 = 2.897_771_955e-3;
        for t in [300.0, 1000.0, 3000.0, 5772.0] {
            let expect = B / t;
            let (mut best, mut best_at) = (0.0, 0.0);
            // Scan a decade either side of the expected peak.
            let (lo, hi) = (expect * 0.1, expect * 10.0);
            let steps = 200_000;
            for i in 0..=steps {
                let l = lo + (hi - lo) * i as f64 / steps as f64;
                let m = planck_exitance(l, t);
                if m > best {
                    best = m;
                    best_at = l;
                }
            }
            let rel = (best_at - expect).abs() / expect;
            assert!(
                rel < 1e-4,
                "at {t} K the peak is at {best_at:e} m and Wien says {expect:e}"
            );
        }
    }

    /// **Stefan–Boltzmann**: the Planck curve integrates to `σT⁴` over all wavelengths. This is
    /// the check on the prefactor `2πhc²` and on the constants together — get any of them wrong
    /// and the total is wrong by that factor.
    ///
    /// Integrated from `0.01·λ_max` to `100·λ_max` on a logarithmic grid, which contains all but
    /// a few parts in 10⁵ of the area; the tolerance is that truncation plus the quadrature, and
    /// 0.1% has room for both.
    #[test]
    fn the_planck_curve_integrates_to_stefan_boltzmann() {
        for t in [300.0f64, 1500.0, 5772.0] {
            let peak = 2.897_771_955e-3 / t;
            let (lo, hi) = (peak * 0.01, peak * 100.0);
            let steps: usize = 400_000;
            let ratio = (hi / lo).powf(1.0 / steps as f64);
            let mut total = 0.0;
            let mut l = lo;
            for _ in 0..steps {
                let next = l * ratio;
                // Trapezoid on each log-spaced interval.
                total += 0.5 * (planck_exitance(l, t) + planck_exitance(next, t)) * (next - l);
                l = next;
            }
            let expect = SIGMA * t.powi(4);
            let rel = (total - expect).abs() / expect;
            assert!(
                rel < 1e-3,
                "at {t} K the curve integrates to {total:e} and σT⁴ is {expect:e} ({rel:e})"
            );
        }
    }

    /// **Equal-energy white lands on `x = y = 1/3`.** The CIE E illuminant is a flat spectrum, and
    /// its chromaticity is one third by *definition* of how the matching functions are normalised
    /// — so this checks the fit's three curves against each other and the integrator at once. The
    /// tolerance is the analytic fit's own error against the tables, about one part in a hundred.
    #[test]
    fn equal_energy_white_lands_on_the_equal_energy_point() {
        let (x, y, z) = xyz_of(|_| 1.0);
        let sum = x + y + z;
        let (cx, cy) = (x / sum, y / sum);
        assert!(
            (cx - 1.0 / 3.0).abs() < 0.01 && (cy - 1.0 / 3.0).abs() < 0.01,
            "equal-energy white is at ({cx:.4}, {cy:.4}) and should be (0.3333, 0.3333)"
        );
    }

    /// **The Planckian locus passes through CIE Illuminant A**, whose chromaticity is
    /// `(0.44757, 0.40745)` — a **published** coordinate, defined as a Planckian radiator at
    /// 2856 K. Nothing in this crate had a hand in producing that pair, which is what makes it a
    /// check rather than a restatement.
    ///
    /// The tolerance is the analytic colour-matching fit's error against the tabulated standard,
    /// which its authors state as about a percent of chromaticity; 0.005 in each coordinate is
    /// that, and the measured miss is well inside it.
    #[test]
    fn the_planckian_locus_passes_through_illuminant_a() {
        let (x, y) = planckian_chromaticity(2856.0);
        assert!(
            (x - 0.447_57).abs() < 0.005 && (y - 0.407_45).abs() < 0.005,
            "2856 K comes out at ({x:.5}, {y:.5}); Illuminant A is (0.44757, 0.40745)"
        );
    }

    /// **Hotter is bluer, monotonically** — the direction the locus runs, and the reason a
    /// photographer's "warm" light is the cooler one.
    ///
    /// Checked on the **chromaticity**, not on the rendered channels, and the difference is a
    /// lesson this test learned by failing: below about 2000 K the blue channel is clamped to
    /// zero — those colours are outside sRGB's gamut entirely — so a blue-to-red ratio is
    /// `0/255` at 1000 K *and* at 1500 K, and monotonicity is untestable there through a monitor
    /// that cannot show either. The CIE `x` coordinate carries the same physics with no gamut in
    /// the way, and falls without exception.
    #[test]
    fn hotter_is_bluer_all_the_way_up() {
        let mut last = f64::INFINITY;
        for t in [1000.0, 1500.0, 2000.0, 3000.0, 4500.0, 6500.0, 10000.0] {
            let (x, _) = planckian_chromaticity(t);
            assert!(
                x < last,
                "{t} K is not bluer than the step below it: x = {x:.4}"
            );
            last = x;
        }
        // And the ends are the colours everybody knows: a 1000 K body is red, out of gamut on the
        // blue side, and a 10000 K body is blue-white with the blue channel saturated.
        let cool = blackbody_srgb(1000.0);
        let hot = blackbody_srgb(10000.0);
        assert_eq!(
            cool[2], 0,
            "a 1000 K body is redder than sRGB's red primary"
        );
        assert!(cool[0] > cool[1], "and should be red-dominant: {cool:?}");
        assert_eq!(hot[2], 255, "a 10000 K body should saturate blue");
    }

    /// **A body at room temperature does not glow, and the number says so rather than a
    /// threshold.** The visible fraction at 300 K is around `1e-30`; at 1000 K, a dull red heat,
    /// it is still under a part in ten thousand; at 3000 K it is percent-scale.
    ///
    /// The steep climb is Wien's tail, `exp(−hc/λk_BT)`, and it is checked as such rather than by
    /// its level: the logarithm of the visible **exitance** is very nearly linear in `1/T` with
    /// slope `−hc/λ_eff·k_B`, so inverting the measured slope must return a wavelength the band
    /// actually contains. That is a closed-form *behaviour*, and it is the stronger check — a
    /// wrong `h`, `c` or `k_B` moves the slope, where a level check could be absorbed by a
    /// prefactor.
    ///
    /// The `σT⁴` of [`glow_fraction`] is divided back out first, and that division is the point
    /// of this comment: with it left in, the `−4 ln T` term shifts the apparent slope and the
    /// inversion returns 944 nm — outside the band, and a check written against *that* number
    /// would have encoded the mistake as the specification. Without it the slope returns 787 to
    /// 798 nm across the pairs below: inside the band and at its red end, which is where the
    /// integral's weight sits when the blue end is exponentially suppressed.
    #[test]
    fn a_cool_body_emits_no_visible_light_and_the_climb_is_wiens_tail() {
        assert!(
            glow_fraction(300.0) < 1e-20,
            "room temperature should emit essentially no visible light, got {:e}",
            glow_fraction(300.0)
        );
        assert!(glow_fraction(1000.0) < 1e-4);
        assert!(glow_fraction(3000.0) > 0.01 && glow_fraction(3000.0) < 0.2);

        for (t1, t2) in [(600.0f64, 900.0f64), (500.0, 800.0), (700.0, 1000.0)] {
            // The visible exitance itself: the fraction with its normalisation put back.
            let e = |t: f64| glow_fraction(t) * SIGMA * t.powi(4);
            let slope = (e(t2).ln() - e(t1).ln()) / (1.0 / t2 - 1.0 / t1);
            let lambda = -H * C / (K_B * slope);
            assert!(
                (700e-9..=830e-9).contains(&lambda),
                "between {t1} K and {t2} K the tail's slope implies {lambda:e} m, which should \
                 sit at the red end of the 360–830 nm band"
            );
        }
    }

    /// **Zero and nonsense are values, not panics.** A field can carry an absolute zero or a
    /// negative kelvin from a caller's mistake, and a renderer must not die on one.
    #[test]
    fn the_degenerate_temperatures_are_black_rather_than_a_panic() {
        assert_eq!(blackbody_srgb(0.0), [0, 0, 0]);
        assert_eq!(blackbody_srgb(-5.0), [0, 0, 0]);
        assert_eq!(glow_fraction(0.0), 0.0);
        assert_eq!(planck_exitance(500e-9, 0.0), 0.0);
        assert_eq!(planck_exitance(-1.0, 300.0), 0.0);
    }
}
