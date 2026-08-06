//! Walking a field from one plane to another.
//!
//! [`wavefront`](crate::wavefront) transforms a pupil straight into an image: one
//! plane, one transform, and no notion of distance. That covers a lens's focus and
//! nothing else. A beam crossing a bench expands, a diffracting edge casts a fringed
//! shadow that changes with how far away the screen is, and neither is expressible as a
//! single transform.
//!
//! This module propagates a complex field a stated distance. The method is the angular
//! spectrum: transform to spatial frequencies, multiply by the phase each plane wave
//! accumulates over the distance, transform back.
//!
//! ```text
//! U(z) = F⁻¹{ F{U(0)} · exp(i 2π z √(1/λ² − fx² − fy²)) }
//! ```
//!
//! No paraxial approximation is made in that transfer function — it is exact for
//! propagation in free space — so it holds at large angles where the Fresnel
//! approximation does not. Evanescent components, where `fx² + fy² > 1/λ²`, decay
//! rather than propagate and are discarded; keeping them would amplify numerical noise
//! into nonsense.
//!
//! # Checked against a Gaussian beam
//!
//! A Gaussian beam is the one field whose free-space propagation has a closed form all
//! the way through: the waist grows as `w(z) = w₀√(1 + (z/z_R)²)` with
//! `z_R = πw₀²/λ`, and far from the waist it opens into a cone of half-angle
//! `λ/(πw₀)`. So the numerical propagator can be checked against an analytic answer at
//! every distance rather than only at a focus, which is what the rest of this crate
//! does with Bessel functions and Planck's law.
//!
//! # The sampling trap, which is worse here than in `wavefront`
//!
//! The transfer function is a phase that winds faster the further the field goes. Once
//! it advances by more than π between neighbouring frequency samples the grid cannot
//! represent it, and the result is not merely inaccurate — it aliases, wrapping energy
//! around the array as though the field were periodic.
//!
//! [`Grid::max_exact_distance`] reports where that happens and
//! [`Grid::propagate`] refuses to go past it. Silently returning a wrapped field is the
//! one behaviour that would be worse than an error, because it looks like a picture.

use std::f64::consts::{PI, TAU};

use dualis_units::Length;

use dualis_core::transform::{fft2, fftshift, ifft2};

/// A complex field sampled on a square grid with a stated physical spacing.
///
/// The spacing is what makes this different from a [`Pupil`](crate::wavefront::Pupil):
/// a pupil is dimensionless in extent and only its diameter in samples matters, while a
/// propagating field has to know how many metres a pixel is.
#[derive(Clone, Debug)]
pub struct Grid {
    samples: usize,
    /// Sample spacing.
    pitch: f64,
    wavelength: f64,
    re: Vec<f64>,
    im: Vec<f64>,
}

impl Grid {
    /// An empty field.
    pub fn new(samples: usize, pitch: Length, wavelength: Length) -> Grid {
        assert!(
            samples.is_power_of_two() && samples >= 8,
            "the transform is radix two, so the grid must be a power of two"
        );
        assert!(pitch.to_si() > 0.0, "a sample must cover some distance");
        assert!(wavelength.to_si() > 0.0, "light must have a wavelength");
        Grid {
            samples,
            pitch: pitch.to_si(),
            wavelength: wavelength.to_si(),
            re: vec![0.0; samples * samples],
            im: vec![0.0; samples * samples],
        }
    }

    /// A field from a function of position, returning `(amplitude, phase in waves)`.
    pub fn from_fn(
        samples: usize,
        pitch: Length,
        wavelength: Length,
        field: impl Fn(Length, Length) -> (f64, f64),
    ) -> Grid {
        let mut grid = Grid::new(samples, pitch, wavelength);
        let centre = samples as f64 / 2.0;
        for y in 0..samples {
            for x in 0..samples {
                let px = Length::from_si((x as f64 + 0.5 - centre) * grid.pitch);
                let py = Length::from_si((y as f64 + 0.5 - centre) * grid.pitch);
                let (a, waves) = field(px, py);
                let phase = TAU * waves;
                let i = y * samples + x;
                grid.re[i] = a * phase.cos();
                grid.im[i] = a * phase.sin();
            }
        }
        grid
    }

    /// A Gaussian beam at its waist: amplitude `exp(-r²/w₀²)`, flat phase.
    ///
    /// The one field whose propagation is a closed form, and therefore the one worth
    /// building in.
    pub fn gaussian_waist(
        samples: usize,
        pitch: Length,
        wavelength: Length,
        waist: Length,
    ) -> Grid {
        let w0 = waist.to_si();
        Grid::from_fn(samples, pitch, wavelength, |x, y| {
            let r2 = x.to_si() * x.to_si() + y.to_si() * y.to_si();
            ((-r2 / (w0 * w0)).exp(), 0.0)
        })
    }

    /// A clear circular aperture, uniformly illuminated.
    pub fn circular_aperture(
        samples: usize,
        pitch: Length,
        wavelength: Length,
        radius: Length,
    ) -> Grid {
        let r = radius.to_si();
        Grid::from_fn(samples, pitch, wavelength, |x, y| {
            let d = (x.to_si() * x.to_si() + y.to_si() * y.to_si()).sqrt();
            (f64::from(d <= r), 0.0)
        })
    }

    pub fn samples(&self) -> usize {
        self.samples
    }

    pub fn pitch(&self) -> Length {
        Length::from_si(self.pitch)
    }

    pub fn wavelength(&self) -> Length {
        Length::from_si(self.wavelength)
    }

    /// The grid's physical width.
    pub fn width(&self) -> Length {
        Length::from_si(self.samples as f64 * self.pitch)
    }

    /// Intensity at a sample.
    pub fn intensity_at(&self, x: usize, y: usize) -> f64 {
        let i = y * self.samples + x;
        self.re[i] * self.re[i] + self.im[i] * self.im[i]
    }

    /// Intensity everywhere.
    pub fn intensity(&self) -> Vec<f64> {
        (0..self.samples * self.samples)
            .map(|i| self.re[i] * self.re[i] + self.im[i] * self.im[i])
            .collect()
    }

    /// Total intensity, which free-space propagation must not change.
    pub fn power(&self) -> f64 {
        self.intensity().iter().sum()
    }

    /// Second-moment radius of the intensity, which is what "the width of a beam"
    /// means for anything that is not a hard edge.
    ///
    /// For a Gaussian `exp(-2r²/w²)` in intensity this returns exactly `w/2`, so a
    /// beam's `w` is twice it. Defined for any profile, which a
    /// full-width-at-half-maximum is not.
    pub fn rms_radius(&self) -> Length {
        let intensity = self.intensity();
        let total: f64 = intensity.iter().sum();
        if total <= 0.0 {
            return Length::ZERO;
        }
        let centre = self.samples as f64 / 2.0;
        let mut second = 0.0;
        for y in 0..self.samples {
            for x in 0..self.samples {
                let px = (x as f64 + 0.5 - centre) * self.pitch;
                let py = (y as f64 + 0.5 - centre) * self.pitch;
                second += intensity[y * self.samples + x] * (px * px + py * py);
            }
        }
        Length::from_si((second / total / 2.0).sqrt())
    }

    /// Beam radius in the usual `1/e²` intensity convention, `w = 2 × rms_radius`.
    pub fn beam_radius(&self) -> Length {
        Length::from_si(self.rms_radius().to_si() * 2.0)
    }

    /// Largest distance this grid can propagate before the transfer function's phase
    /// aliases.
    ///
    /// The phase advances by `2πz(f_max² λ/2)` at the corner of the frequency grid in
    /// the paraxial limit, and it must not advance by more than π between neighbouring
    /// samples. Working that through for a grid of `N` samples at pitch `Δ`:
    /// `z_max = N Δ² / λ`.
    ///
    /// Which is the same as saying the grid must be wide enough to hold the light that
    /// diffracts across it: past `z_max` a ray leaving one edge at the largest angle
    /// the grid can represent has crossed the whole array.
    pub fn max_exact_distance(&self) -> Length {
        Length::from_si(self.samples as f64 * self.pitch * self.pitch / self.wavelength)
    }

    /// Propagate a distance in free space, by the angular spectrum method.
    ///
    /// Fails rather than aliasing. See [`Grid::max_exact_distance`]: a wrapped field
    /// looks like a picture, and there is no way for a caller to notice.
    pub fn propagate(&self, distance: Length) -> Result<Grid, PropagationError> {
        let z = distance.to_si();
        if z == 0.0 {
            return Ok(self.clone());
        }
        let limit = self.max_exact_distance().to_si();
        if z.abs() > limit {
            return Err(PropagationError::TooFar {
                requested: distance,
                limit: Length::from_si(limit),
            });
        }

        let n = self.samples;
        let mut re = self.re.clone();
        let mut im = self.im.clone();
        fft2(&mut re, &mut im, n);

        // Frequency of each sample, in cycles per metre, with the usual FFT layout:
        // 0..n/2 are positive, n/2..n are negative.
        let df = 1.0 / (n as f64 * self.pitch);
        let inv_lambda2 = 1.0 / (self.wavelength * self.wavelength);
        let frequency = |k: usize| {
            let signed = if k < n / 2 {
                k as f64
            } else {
                k as f64 - n as f64
            };
            signed * df
        };

        for ky in 0..n {
            let fy = frequency(ky);
            for kx in 0..n {
                let fx = frequency(kx);
                let under = inv_lambda2 - fx * fx - fy * fy;
                let i = ky * n + kx;
                if under <= 0.0 {
                    // Evanescent: this plane wave does not propagate, it decays. Over
                    // any distance worth simulating it is gone, and keeping it would
                    // amplify rounding noise into a field.
                    re[i] = 0.0;
                    im[i] = 0.0;
                    continue;
                }
                let phase = TAU * z * under.sqrt();
                let (c, s) = (phase.cos(), phase.sin());
                let (a, b) = (re[i], im[i]);
                re[i] = a * c - b * s;
                im[i] = a * s + b * c;
            }
        }

        ifft2(&mut re, &mut im, n);
        Ok(Grid {
            samples: n,
            pitch: self.pitch,
            wavelength: self.wavelength,
            re,
            im,
        })
    }

    /// Intensity along the x axis through the centre, as `(position, intensity)`.
    pub fn profile_x(&self) -> Vec<(Length, f64)> {
        let c = self.samples / 2;
        let centre = self.samples as f64 / 2.0;
        (0..self.samples)
            .map(|x| {
                (
                    Length::from_si((x as f64 + 0.5 - centre) * self.pitch),
                    self.intensity_at(x, c),
                )
            })
            .collect()
    }

    /// The intensity array with the centre in the middle, for display.
    pub fn centred_intensity(&self) -> Vec<f64> {
        fftshift(&self.intensity(), self.samples)
    }
}

/// Why a propagation was refused.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PropagationError {
    /// Past the distance this grid can represent without aliasing.
    TooFar { requested: Length, limit: Length },
}

impl std::fmt::Display for PropagationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PropagationError::TooFar { requested, limit } => write!(
                f,
                "cannot propagate {:.4} m on this grid: past {:.4} m the transfer \
                 function's phase aliases and the field would wrap around the array",
                requested.to_si(),
                limit.to_si()
            ),
        }
    }
}

impl std::error::Error for PropagationError {}

/// Rayleigh range of a Gaussian beam, `πw₀²/λ`.
///
/// The distance over which a beam stays roughly collimated: at `z_R` its radius has
/// grown by `√2` and its area has doubled. Everything about how a beam behaves is
/// really a statement about `z/z_R`.
pub fn rayleigh_range(waist: Length, wavelength: Length) -> Length {
    Length::from_si(PI * waist.to_si().powi(2) / wavelength.to_si())
}

/// Radius of a Gaussian beam at a distance from its waist: `w₀√(1 + (z/z_R)²)`.
pub fn gaussian_radius_at(waist: Length, wavelength: Length, distance: Length) -> Length {
    let zr = rayleigh_range(waist, wavelength).to_si();
    let ratio = distance.to_si() / zr;
    Length::from_si(waist.to_si() * (1.0 + ratio * ratio).sqrt())
}

/// Far-field half-angle of a Gaussian beam, `λ/(πw₀)`, in radians.
///
/// The tradeoff every beam obeys: a tighter waist diverges faster, and the product of
/// the two is fixed by the wavelength. There is no way to have a small spot that stays
/// small.
pub fn gaussian_divergence(waist: Length, wavelength: Length) -> f64 {
    wavelength.to_si() / (PI * waist.to_si())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn um(v: f64) -> Length {
        Length::um(v)
    }

    /// A beam propagated no distance is the beam.
    #[test]
    fn propagating_nowhere_changes_nothing() {
        let beam = Grid::gaussian_waist(128, um(2.0), Length::nm(633.0), um(20.0));
        let same = beam.propagate(Length::ZERO).unwrap();
        for i in 0..beam.samples * beam.samples {
            assert_eq!(beam.re[i].to_bits(), same.re[i].to_bits());
            assert_eq!(beam.im[i].to_bits(), same.im[i].to_bits());
        }
    }

    /// Free space does not absorb, so propagation conserves power. The cheapest check
    /// that the transfer function is a pure phase and not accidentally a filter.
    #[test]
    fn propagation_conserves_power() {
        let beam = Grid::gaussian_waist(256, um(2.0), Length::nm(633.0), um(30.0));
        let start = beam.power();
        assert!(start > 0.0);
        let mut running = beam;
        for _ in 0..5 {
            running = running.propagate(um(200.0)).unwrap();
            assert!(
                (running.power() / start - 1.0).abs() < 1e-9,
                "power changed to {}",
                running.power() / start
            );
        }
    }

    /// Going forward and back returns the field, which catches a sign error in the
    /// transfer function that a magnitude-only test would not.
    #[test]
    fn propagation_is_reversible() {
        // 12 um waist at 633 nm has a Rayleigh range of 0.715 mm, so 500 um is 0.7 of
        // it and the beam grows by 22% on the way — enough that returning to the
        // original is not a trivial pass. The grid reaches 0.809 mm, which covers it.
        let beam = Grid::gaussian_waist(128, um(2.0), Length::nm(633.0), um(12.0));
        let there = beam.propagate(um(500.0)).unwrap();
        let back = there.propagate(-um(500.0)).unwrap();
        for i in 0..beam.samples * beam.samples {
            assert!(
                (back.re[i] - beam.re[i]).abs() < 1e-9,
                "real part at {i}: {} against {}",
                back.re[i],
                beam.re[i]
            );
            assert!((back.im[i] - beam.im[i]).abs() < 1e-9);
        }
        assert!(
            there.beam_radius() > beam.beam_radius() * 1.15,
            "the beam should have expanded on the way: {} um to {} um",
            beam.beam_radius().in_um(),
            there.beam_radius().in_um()
        );
    }

    /// **The test this module is checked by.** A Gaussian beam's radius follows
    /// `w₀√(1 + (z/z_R)²)` at every distance, and the numerical propagator reproduces
    /// it — an analytic answer available all the way along rather than only at a focus.
    #[test]
    fn a_gaussian_beam_expands_by_its_closed_form() {
        let wavelength = Length::nm(633.0);
        let waist = um(12.0);
        let zr = rayleigh_range(waist, wavelength);
        // pi w0^2 / lambda = pi * 1.44e-10 / 633e-9 = 0.715 mm.
        assert!(
            (zr.in_mm() - 0.7147).abs() < 0.002,
            "Rayleigh range {} mm",
            zr.in_mm()
        );

        // The grid has to reach past the furthest distance asked for: 512 samples at
        // 1.5 um gives 1.82 mm, which is 2.5 Rayleigh ranges. Choosing a waist first and
        // a grid second is the wrong order — the reach is what constrains the run.
        let beam = Grid::gaussian_waist(512, um(1.5), wavelength, waist);
        assert!(beam.max_exact_distance() > zr * 2.0);
        // At the waist the measured radius is the one it was built with.
        assert!(
            (beam.beam_radius().in_um() / waist.in_um() - 1.0).abs() < 5e-3,
            "at the waist: {} um",
            beam.beam_radius().in_um()
        );

        for fraction in [0.25f64, 0.5, 1.0, 1.5] {
            let z = Length::from_si(zr.to_si() * fraction);
            let propagated = beam.propagate(z).expect("within the grid's range");
            let measured = propagated.beam_radius();
            let expected = gaussian_radius_at(waist, wavelength, z);
            assert!(
                (measured / expected - 1.0).abs() < 0.02,
                "at z = {fraction} z_R the propagator gives {} um and the closed form \
                 gives {} um",
                measured.in_um(),
                expected.in_um()
            );
        }
        // At one Rayleigh range the radius has grown by exactly root two.
        let at_zr = beam.propagate(zr).unwrap().beam_radius();
        assert!(
            (at_zr / waist - 2f64.sqrt()).abs() < 0.03,
            "root two at z_R, got {}",
            at_zr / waist
        );
    }

    /// The divergence a beam settles into, and the trade behind it: a tighter waist
    /// opens faster, and the product is fixed by the wavelength alone.
    #[test]
    fn a_tighter_waist_diverges_faster() {
        let wavelength = Length::nm(633.0);
        let tight = gaussian_divergence(um(10.0), wavelength);
        let loose = gaussian_divergence(um(40.0), wavelength);
        assert!(
            (tight / loose - 4.0).abs() < 1e-12,
            "divergence goes as 1/w0"
        );
        // The product of waist and divergence is lambda/pi, whatever the waist.
        for w in [um(5.0), um(20.0), um(100.0)] {
            let product = w.to_si() * gaussian_divergence(w, wavelength);
            assert!((product / (wavelength.to_si() / PI) - 1.0).abs() < 1e-12);
        }
        // 633 nm through a 20 um waist opens at 10 mrad, which is 0.58 degrees.
        let theta = gaussian_divergence(um(20.0), wavelength);
        assert!((theta * 1e3 - 10.07).abs() < 0.05, "{} mrad", theta * 1e3);
    }

    /// Far from the waist the beam is a cone, and its measured opening angle is the
    /// closed-form divergence. Checked by propagating and differencing rather than by
    /// trusting the formula.
    #[test]
    fn the_far_field_opens_at_the_divergence_angle() {
        let wavelength = Length::nm(633.0);
        let waist = um(12.0);
        let zr = rayleigh_range(waist, wavelength);
        // Two Rayleigh ranges is 1.43 mm, so the grid needs a 2 um pitch to reach it.
        let beam = Grid::gaussian_waist(512, um(2.0), wavelength, waist);
        assert!(beam.max_exact_distance() > zr * 2.0);

        let near = Length::from_si(zr.to_si() * 1.0);
        let far = Length::from_si(zr.to_si() * 2.0);
        let r_near = beam.propagate(near).unwrap().beam_radius();
        let r_far = beam.propagate(far).unwrap().beam_radius();

        let measured_angle = (r_far - r_near).to_si() / (far - near).to_si();

        // Against the *exact secant* over the same interval, not against the asymptote.
        // Those are different numbers and the difference is not small: `w = w0√(1+u²)`
        // gives a secant of `(√5 − √2) w0/z_R` between one and two Rayleigh ranges,
        // which is 0.822 of the limit. Comparing a measurement taken two Rayleigh ranges
        // out against the asymptote would be an 18% disagreement that says nothing about
        // the propagator.
        let exact_secant = (gaussian_radius_at(waist, wavelength, far)
            - gaussian_radius_at(waist, wavelength, near))
        .to_si()
            / (far - near).to_si();
        assert!(
            (measured_angle / exact_secant - 1.0).abs() < 0.03,
            "measured {measured_angle:e} rad against an exact secant of {exact_secant:e}"
        );

        // And the secant really does tend to the divergence, which is a fact about the
        // closed form and needs no propagation to check. Ten Rayleigh ranges out it is
        // within a percent.
        let divergence = gaussian_divergence(waist, wavelength);
        assert!(
            (exact_secant / divergence - 0.822).abs() < 0.01,
            "the near secant should be 0.822 of the limit, got {}",
            exact_secant / divergence
        );
        let ten = Length::from_si(zr.to_si() * 10.0);
        let eleven = Length::from_si(zr.to_si() * 11.0);
        let far_secant = (gaussian_radius_at(waist, wavelength, eleven)
            - gaussian_radius_at(waist, wavelength, ten))
        .to_si()
            / (eleven - ten).to_si();
        assert!(
            (far_secant / divergence - 1.0).abs() < 0.01,
            "ten Rayleigh ranges out the beam is a cone, got {}",
            far_secant / divergence
        );
    }

    /// A propagation the grid cannot represent is refused, and the message says what
    /// the limit was. Returning a wrapped field would be the one failure a caller could
    /// not detect, because it looks like a picture.
    #[test]
    fn too_far_is_an_error_and_not_a_wrapped_field() {
        let beam = Grid::gaussian_waist(128, um(2.0), Length::nm(633.0), um(20.0));
        // N * pitch^2 / lambda = 128 * 4e-12 / 633e-9 = 0.809 mm.
        let limit = beam.max_exact_distance();
        assert!(
            (limit.in_mm() - 0.809).abs() < 0.005,
            "limit {} mm",
            limit.in_mm()
        );
        assert!(beam.propagate(limit * 0.99).is_ok());

        let err = beam
            .propagate(limit * 1.01)
            .expect_err("past the limit must be refused");
        match err {
            PropagationError::TooFar {
                requested,
                limit: l,
            } => {
                assert!(requested > l);
            }
        }
        assert!(err.to_string().contains("alias"), "{err}");
        // Backwards too: the limit is on the magnitude.
        assert!(beam.propagate(-limit * 1.01).is_err());
    }

    /// A finer grid reaches further, which is the trade the limit expresses: `z_max`
    /// goes as `N Δ²`, so doubling the pitch quadruples the reach at the cost of
    /// resolving less.
    #[test]
    fn the_reach_is_set_by_the_grid_and_the_wavelength() {
        let green = Length::nm(633.0);
        let a = Grid::gaussian_waist(128, um(2.0), green, um(20.0));
        let b = Grid::gaussian_waist(128, um(4.0), green, um(20.0));
        assert!(
            (b.max_exact_distance() / a.max_exact_distance() - 4.0).abs() < 1e-12,
            "doubling the pitch should quadruple the reach"
        );
        let wide = Grid::gaussian_waist(256, um(2.0), green, um(20.0));
        assert!(
            (wide.max_exact_distance() / a.max_exact_distance() - 2.0).abs() < 1e-12,
            "twice the samples, twice the reach"
        );
        // And a longer wavelength diffracts further, so it reaches less far.
        let red = Grid::gaussian_waist(128, um(2.0), Length::nm(1266.0), um(20.0));
        assert!((red.max_exact_distance() / a.max_exact_distance() - 0.5).abs() < 1e-12);
    }

    /// The thing a single-plane transform cannot express at all: the shadow of an
    /// aperture changes character with distance.
    ///
    /// Close up it is the aperture with fringed edges and a bright middle; far away it
    /// has become the Airy pattern of the far field. Between them the on-axis intensity
    /// oscillates, which is Fresnel diffraction and is exactly what
    /// `diffraction::fresnel_number` was there to identify without being able to
    /// compute.
    #[test]
    fn an_aperture_shadow_changes_with_distance() {
        use crate::diffraction::fresnel_number;

        let wavelength = Length::nm(633.0);
        let radius = um(40.0);
        let field = Grid::circular_aperture(512, um(1.0), wavelength, radius);
        let limit = field.max_exact_distance();

        // Near field: the Fresnel number is large and the shadow still looks like the
        // aperture — most of the light is inside the geometric edge.
        let near = limit * 0.1;
        assert!(
            fresnel_number(radius, near, wavelength) > 3.0,
            "the near case should be near field, F = {}",
            fresnel_number(radius, near, wavelength)
        );
        let near_field = field.propagate(near).unwrap();
        let inside_near = fraction_within(&near_field, radius);
        assert!(
            inside_near > 0.8,
            "close up the light is still inside the aperture's edge, got {inside_near:.3}"
        );

        // Further out it has spread past the geometric shadow.
        let far = limit * 0.95;
        let far_field = field.propagate(far).unwrap();
        let inside_far = fraction_within(&far_field, radius);
        assert!(
            inside_far < inside_near - 0.05,
            "further out the light should have spread: {inside_far:.3} against \
             {inside_near:.3}"
        );

        // On-axis intensity is not monotonic in distance -- the Fresnel oscillation --
        // which is the signature no single-plane transform reproduces.
        let axis = |z: Length| {
            let g = field.propagate(z).unwrap();
            let c = g.samples() / 2;
            g.intensity_at(c, c)
        };
        let samples: Vec<f64> = (1..=12).map(|i| axis(limit * (i as f64 / 13.0))).collect();
        let rises = samples.windows(2).filter(|w| w[1] > w[0]).count();
        assert!(
            rises > 0 && rises < samples.len() - 1,
            "the on-axis intensity should oscillate rather than march one way: \
             {samples:?}"
        );
    }

    /// Fraction of the power inside a radius.
    fn fraction_within(grid: &Grid, radius: Length) -> f64 {
        let intensity = grid.intensity();
        let total: f64 = intensity.iter().sum();
        if total <= 0.0 {
            return 0.0;
        }
        let centre = grid.samples() as f64 / 2.0;
        let r = radius.to_si();
        let mut inside = 0.0;
        for y in 0..grid.samples() {
            for x in 0..grid.samples() {
                let px = (x as f64 + 0.5 - centre) * grid.pitch;
                let py = (y as f64 + 0.5 - centre) * grid.pitch;
                if (px * px + py * py).sqrt() <= r {
                    inside += intensity[y * grid.samples() + x];
                }
            }
        }
        inside / total
    }

    /// Propagating in two hops is propagating once, which is what makes a bench of
    /// several elements meaningful.
    #[test]
    fn propagation_composes() {
        let beam = Grid::gaussian_waist(256, um(1.5), Length::nm(633.0), um(20.0));
        let one_hop = beam.propagate(um(600.0)).unwrap();
        let two_hops = beam
            .propagate(um(250.0))
            .unwrap()
            .propagate(um(350.0))
            .unwrap();
        let a = one_hop.intensity();
        let b = two_hops.intensity();
        let peak = a.iter().cloned().fold(0.0f64, f64::max);
        let worst = a
            .iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).abs() / peak)
            .fold(0.0f64, f64::max);
        assert!(
            worst < 1e-6,
            "one hop and two should agree, worst difference {worst:e}"
        );
    }

    /// The propagator is deterministic, like everything else here.
    #[test]
    fn propagation_is_bit_reproducible() {
        let build = || {
            Grid::gaussian_waist(128, um(2.0), Length::nm(633.0), um(15.0))
                .propagate(um(400.0))
                .unwrap()
        };
        let (a, b) = (build(), build());
        for i in 0..a.samples * a.samples {
            assert_eq!(a.re[i].to_bits(), b.re[i].to_bits());
            assert_eq!(a.im[i].to_bits(), b.im[i].to_bits());
        }
    }

    /// Grids that cannot work are refused at construction.
    #[test]
    #[should_panic(expected = "power of two")]
    fn a_non_power_of_two_grid_is_refused() {
        Grid::new(100, um(1.0), Length::nm(633.0));
    }

    #[test]
    #[should_panic(expected = "some distance")]
    fn a_zero_pitch_is_refused() {
        Grid::new(64, Length::ZERO, Length::nm(633.0));
    }
}
