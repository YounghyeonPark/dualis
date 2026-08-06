//! Aberrated point spread functions, from a pupil rather than from a formula.
//!
//! [`diffraction`](crate::diffraction) answers what a *perfect* system does: the Airy
//! pattern, the ideal MTF, the Rayleigh limit. Those are ceilings. A real instrument
//! sits below them, and how far below is the whole of what a lens is judged on — so
//! something has to compute the imperfect case, and no closed form covers it.
//!
//! The route is the standard one. A pupil carries a complex amplitude: how much light
//! gets through at each point, and how far out of step it arrives. Squaring the
//! Fourier transform of that gives the intensity in the image, which is the point
//! spread function. Aberrations enter as the phase, expanded in
//! [`Zernike`] polynomials because that is the basis an optical shop measures in.
//!
//! # Every result here is checked against the analytic module
//!
//! That is the point of the arrangement. Set the aberrations to zero and this module
//! must reproduce [`airy_intensity`](crate::diffraction::airy_intensity), which is
//! computed from Bessel functions and shares no code with it. Take the transform of
//! the PSF and it must reproduce [`mtf_ideal`](crate::diffraction::mtf_ideal), which
//! is the closed-form autocorrelation of a disc. Put in a small aberration and the
//! Strehl ratio must follow
//! [`strehl_from_wavefront_error`](crate::diffraction::strehl_from_wavefront_error),
//! the Maréchal approximation. Three independent agreements, none of them arranged.
//!
//! # Sampling, which is the thing that actually goes wrong
//!
//! A pupil of diameter `d` samples in an `N`-by-`N` grid produces a PSF sampled at
//! `N/d` pixels per `λ/D`. Two failure modes follow and both are silent:
//!
//! - **Too few pixels across the pupil** and its edge is a staircase, which puts
//!   energy in the wings that the optics did not.
//! - **Too little padding** (`N` close to `d`) and the PSF is sampled so coarsely
//!   that its core falls between pixels. The first Airy zero lands at `1.22 N/d`
//!   pixels, so `N = 4d` puts it at five — enough for a Strehl ratio, not enough for
//!   a profile.
//!
//! [`Pupil::circular`] takes both numbers rather than choosing for you, and
//! [`Psf::scale`] reports what they bought.

use std::f64::consts::{PI, TAU};

use crate::diffraction::FIRST_AIRY_ZERO;

/// One Zernike polynomial, by its radial order `n` and azimuthal frequency `m`.
///
/// Orthonormal over the unit disc: a coefficient given in waves **is** that mode's
/// RMS contribution to the wavefront error, so coefficients from different modes add
/// in quadrature and feed straight into a Strehl estimate. The other common
/// normalisation makes a coefficient the peak value instead, which is why two
/// wavefront specs quoting "0.1 waves of coma" can differ by a factor of `√8`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Zernike {
    pub n: u32,
    pub m: i32,
}

impl Zernike {
    /// Piston: a constant phase, which no image can see.
    pub const PISTON: Zernike = Zernike { n: 0, m: 0 };
    /// Tilt, which moves the image rather than blurring it.
    pub const TILT_X: Zernike = Zernike { n: 1, m: 1 };
    pub const TILT_Y: Zernike = Zernike { n: 1, m: -1 };
    /// Defocus. The one aberration a focus knob can remove.
    pub const DEFOCUS: Zernike = Zernike { n: 2, m: 0 };
    /// Astigmatism at 0 and 45 degrees.
    pub const ASTIGMATISM_0: Zernike = Zernike { n: 2, m: 2 };
    pub const ASTIGMATISM_45: Zernike = Zernike { n: 2, m: -2 };
    /// Coma, which makes a point look like a comet and does not average out.
    pub const COMA_X: Zernike = Zernike { n: 3, m: 1 };
    pub const COMA_Y: Zernike = Zernike { n: 3, m: -1 };
    pub const TREFOIL_0: Zernike = Zernike { n: 3, m: 3 };
    /// Spherical aberration: the one a single spherical surface cannot avoid.
    pub const SPHERICAL: Zernike = Zernike { n: 4, m: 0 };

    pub fn new(n: u32, m: i32) -> Zernike {
        Zernike { n, m }
    }

    /// Whether this is a valid mode: `|m| ≤ n` and `n - |m|` even.
    // Clippy on current stable suggests `u32::is_multiple_of`, which was stabilised in
    // 1.87 — later than the 1.78 this workspace declares and its CI verifies. A
    // declared MSRV is a promise to a consumer and a lint suggestion is not, so the
    // MSRV wins. Revisit if the floor ever moves past 1.87.
    #[allow(clippy::manual_is_multiple_of)]
    pub fn is_valid(&self) -> bool {
        let am = self.m.unsigned_abs();
        am <= self.n && (self.n - am) % 2 == 0
    }

    /// The polynomial at a point in polar coordinates, `rho` in 0..1.
    ///
    /// Zero outside the unit disc, where it is not defined and not orthonormal.
    pub fn value(&self, rho: f64, theta: f64) -> f64 {
        if !self.is_valid() || rho > 1.0 {
            return 0.0;
        }
        let am = self.m.unsigned_abs();
        let radial = radial_polynomial(self.n, am, rho);
        // sqrt(2(n+1)) for m != 0, sqrt(n+1) for m == 0: the factor that makes the
        // RMS over the disc exactly one.
        let norm = if am == 0 {
            ((self.n + 1) as f64).sqrt()
        } else {
            (2.0 * (self.n + 1) as f64).sqrt()
        };
        let angular = if self.m >= 0 {
            (am as f64 * theta).cos()
        } else {
            (am as f64 * theta).sin()
        };
        norm * radial * angular
    }

    /// A human name, for the modes that have one.
    pub fn name(&self) -> &'static str {
        match (self.n, self.m) {
            (0, 0) => "piston",
            (1, 1) => "tilt x",
            (1, -1) => "tilt y",
            (2, 0) => "defocus",
            (2, 2) => "astigmatism 0",
            (2, -2) => "astigmatism 45",
            (3, 1) => "coma x",
            (3, -1) => "coma y",
            (3, 3) => "trefoil 0",
            (3, -3) => "trefoil 30",
            (4, 0) => "spherical",
            _ => "higher order",
        }
    }
}

/// `R_n^m(rho)`, the radial part.
// See the note on `Zernike::is_valid` for why `is_multiple_of` is not used here.
#[allow(clippy::manual_is_multiple_of)]
fn radial_polynomial(n: u32, m: u32, rho: f64) -> f64 {
    if m > n || (n - m) % 2 != 0 {
        return 0.0;
    }
    let half = (n - m) / 2;
    let mut sum = 0.0;
    for k in 0..=half {
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        let numerator = factorial(n - k);
        let denominator = factorial(k) * factorial((n + m) / 2 - k) * factorial((n - m) / 2 - k);
        sum += sign * (numerator / denominator) * rho.powi((n - 2 * k) as i32);
    }
    sum
}

fn factorial(k: u32) -> f64 {
    (1..=k).map(f64::from).product::<f64>().max(1.0)
}

/// A complex amplitude sampled over a square grid.
///
/// Amplitude is the fraction of light getting through — 1 inside a clear aperture, 0
/// outside, and a partial value at the edge where a pixel is half covered. Phase is
/// the wavefront error in waves, so 0.25 is a quarter wave out of step.
#[derive(Clone, Debug)]
pub struct Pupil {
    samples: usize,
    diameter: f64,
    /// Transmitted amplitude per sample, including the antialiased rim.
    amplitude: Vec<f64>,
    /// Wavefront error in waves per sample.
    wavefront: Vec<f64>,
}

impl Pupil {
    /// A clear circular aperture of `diameter` samples inside a `samples`-square grid.
    ///
    /// `samples` must be a power of two. `diameter` may be fractional; the rim is
    /// antialiased by supersampling, which matters more than it sounds — a hard mask
    /// on a 64-pixel circle is a visible staircase, and a staircase diffracts.
    pub fn circular(samples: usize, diameter: f64) -> Pupil {
        Pupil::annular(samples, diameter, 0.0)
    }

    /// An annular aperture, with a central obscuration of `obscuration` times the
    /// outer radius.
    ///
    /// A reflecting objective or a Cassegrain telescope: the secondary mirror sits in
    /// the middle of the beam. The trade it makes is visible in the PSF — the core
    /// narrows, which improves resolution, while more energy moves into the rings,
    /// which destroys contrast on an extended scene.
    pub fn annular(samples: usize, diameter: f64, obscuration: f64) -> Pupil {
        assert!(
            samples.is_power_of_two() && samples >= 8,
            "the transform is radix two, so the grid must be a power of two"
        );
        assert!(
            diameter > 0.0 && diameter <= samples as f64,
            "the aperture must fit in the grid"
        );
        let obscuration = obscuration.clamp(0.0, 0.99);
        let radius = diameter / 2.0;
        let centre = samples as f64 / 2.0;
        // 4x4 supersampling of the rim. Enough to put the edge's contribution well
        // below the aberrations anyone would be modelling.
        const SUB: usize = 4;
        let mut amplitude = vec![0.0; samples * samples];
        for y in 0..samples {
            for x in 0..samples {
                let mut covered = 0;
                for sy in 0..SUB {
                    for sx in 0..SUB {
                        let px = x as f64 + (sx as f64 + 0.5) / SUB as f64 - centre;
                        let py = y as f64 + (sy as f64 + 0.5) / SUB as f64 - centre;
                        let r = (px * px + py * py).sqrt() / radius;
                        if r <= 1.0 && r >= obscuration {
                            covered += 1;
                        }
                    }
                }
                amplitude[y * samples + x] = covered as f64 / (SUB * SUB) as f64;
            }
        }
        Pupil {
            samples,
            diameter,
            amplitude,
            wavefront: vec![0.0; samples * samples],
        }
    }

    pub fn samples(&self) -> usize {
        self.samples
    }

    pub fn diameter(&self) -> f64 {
        self.diameter
    }

    /// Add a Zernike mode with the given RMS amplitude in waves.
    pub fn with_zernike(self, mode: Zernike, waves_rms: f64) -> Pupil {
        assert!(mode.is_valid(), "{mode:?} is not a Zernike mode");
        self.with_wavefront(|rho, theta| waves_rms * mode.value(rho, theta))
    }

    /// Add an arbitrary wavefront, as a function of `(rho, theta)` returning waves.
    ///
    /// Evaluated at every sample that transmits anything, with `rho` clamped to 1 at
    /// the rim. A pixel the aperture only half covers still carries light and so must
    /// carry a phase; leaving it at zero while its neighbours are aberrated would put
    /// a ring of artificial error around the pupil, and would make even a constant
    /// wavefront register as an aberration.
    pub fn with_wavefront(mut self, error: impl Fn(f64, f64) -> f64) -> Pupil {
        let radius = self.diameter / 2.0;
        let centre = self.samples as f64 / 2.0;
        for y in 0..self.samples {
            for x in 0..self.samples {
                let i = y * self.samples + x;
                if self.amplitude[i] <= 0.0 {
                    continue;
                }
                let px = x as f64 + 0.5 - centre;
                let py = y as f64 + 0.5 - centre;
                let rho = ((px * px + py * py).sqrt() / radius).min(1.0);
                self.wavefront[i] += error(rho, py.atan2(px));
            }
        }
        self
    }

    /// RMS wavefront error in waves, weighted by transmitted amplitude.
    ///
    /// The number a Strehl estimate is built from, and the number an interferometer
    /// reports. Piston is removed first, because a constant phase offset is not an
    /// aberration and including it would make a perfectly good wavefront look bad.
    pub fn rms_wavefront_error(&self) -> f64 {
        let total: f64 = self.amplitude.iter().sum();
        if total <= 0.0 {
            return 0.0;
        }
        let mean: f64 = self
            .amplitude
            .iter()
            .zip(self.wavefront.iter())
            .map(|(a, w)| a * w)
            .sum::<f64>()
            / total;
        let variance: f64 = self
            .amplitude
            .iter()
            .zip(self.wavefront.iter())
            .map(|(a, w)| a * (w - mean) * (w - mean))
            .sum::<f64>()
            / total;
        variance.max(0.0).sqrt()
    }

    /// The point spread function: the squared magnitude of the pupil's transform.
    ///
    /// Normalised so that an unaberrated pupil of the same shape peaks at exactly 1 —
    /// which makes the value at the centre the Strehl ratio, by definition rather
    /// than by a second calculation.
    pub fn psf(&self) -> Psf {
        let n = self.samples;
        let mut re = vec![0.0; n * n];
        let mut im = vec![0.0; n * n];
        for i in 0..n * n {
            let phase = TAU * self.wavefront[i];
            re[i] = self.amplitude[i] * phase.cos();
            im[i] = self.amplitude[i] * phase.sin();
        }
        fft2(&mut re, &mut im, n, false);

        let throughput: f64 = self.amplitude.iter().sum();
        let norm = if throughput > 0.0 {
            1.0 / (throughput * throughput)
        } else {
            0.0
        };
        let mut intensity = vec![0.0; n * n];
        for i in 0..n * n {
            intensity[i] = (re[i] * re[i] + im[i] * im[i]) * norm;
        }
        Psf {
            samples: n,
            intensity: fftshift(&intensity, n),
            scale: n as f64 / self.diameter,
        }
    }
}

/// A sampled point spread function.
#[derive(Clone, Debug)]
pub struct Psf {
    samples: usize,
    /// Intensity, DC at the centre, normalised to a perfect peak of 1.
    intensity: Vec<f64>,
    /// Pixels per `λ/D`.
    scale: f64,
}

impl Psf {
    pub fn samples(&self) -> usize {
        self.samples
    }

    /// Pixels per `λ/D`, which is `N/d`. The first Airy zero sits at `1.22` of these.
    pub fn scale(&self) -> f64 {
        self.scale
    }

    /// Intensity at an offset from the centre, in pixels.
    pub fn at(&self, dx: isize, dy: isize) -> f64 {
        let c = (self.samples / 2) as isize;
        let (x, y) = (c + dx, c + dy);
        if x < 0 || y < 0 || x >= self.samples as isize || y >= self.samples as isize {
            return 0.0;
        }
        self.intensity[y as usize * self.samples + x as usize]
    }

    /// The largest value anywhere.
    pub fn peak(&self) -> f64 {
        self.intensity.iter().cloned().fold(0.0, f64::max)
    }

    /// Strehl ratio: the on-axis intensity relative to an unaberrated pupil of the
    /// same shape.
    ///
    /// Read at the centre rather than at the peak, deliberately. An aberration with a
    /// tilt in it moves the brightest point off axis, and the Strehl ratio is about
    /// how much light lands where the image is supposed to be — [`Psf::peak`] is the
    /// other question.
    pub fn strehl(&self) -> f64 {
        self.at(0, 0)
    }

    /// Sum of the whole array. Not 1: it is `N²Σa²/(Σa)²` under this normalisation,
    /// which is what [`Psf`]'s Parseval check compares against.
    pub fn total_energy(&self) -> f64 {
        self.intensity.iter().sum()
    }

    /// Azimuthally averaged profile, as `(radius in λ/D, intensity)` pairs out to
    /// `max_radius`.
    ///
    /// Bins containing no pixels are omitted rather than reported as zero. Near the
    /// centre a fine binning asks about radii between one pixel and the next, and
    /// there is no measurement there — returning zero for it would be an answer the
    /// data does not contain.
    pub fn radial_profile(&self, max_radius: f64, bins: usize) -> Vec<(f64, f64)> {
        let bins = bins.max(1);
        let mut sums = vec![0.0; bins];
        let mut counts = vec![0usize; bins];
        let c = (self.samples / 2) as isize;
        let limit_px = max_radius * self.scale;
        for y in 0..self.samples as isize {
            for x in 0..self.samples as isize {
                let (dx, dy) = ((x - c) as f64, (y - c) as f64);
                let r = (dx * dx + dy * dy).sqrt();
                if r > limit_px {
                    continue;
                }
                let bin = ((r / limit_px) * bins as f64) as usize;
                let bin = bin.min(bins - 1);
                sums[bin] += self.intensity[y as usize * self.samples + x as usize];
                counts[bin] += 1;
            }
        }
        (0..bins)
            .filter(|&i| counts[i] > 0)
            .map(|i| {
                let r = (i as f64 + 0.5) / bins as f64 * max_radius;
                (r, sums[i] / counts[i] as f64)
            })
            .collect()
    }

    /// Fraction of the array's energy inside a radius, in `λ/D`.
    ///
    /// Against the whole array rather than against an infinite integral, so it is
    /// slightly high: the grid truncates the Airy wings, which fall off slowly enough
    /// that a percent or two of the light is outside any finite box.
    pub fn encircled_energy(&self, radius: f64) -> f64 {
        let total = self.total_energy();
        if total <= 0.0 {
            return 0.0;
        }
        let limit_px = radius * self.scale;
        let c = (self.samples / 2) as isize;
        let mut inside = 0.0;
        for y in 0..self.samples as isize {
            for x in 0..self.samples as isize {
                let (dx, dy) = ((x - c) as f64, (y - c) as f64);
                if (dx * dx + dy * dy).sqrt() <= limit_px {
                    inside += self.intensity[y as usize * self.samples + x as usize];
                }
            }
        }
        inside / total
    }

    /// The modulation transfer function, as the magnitude of the PSF's transform.
    ///
    /// By the Wiener-Khinchin relation this is the pupil's autocorrelation, which is
    /// why an unaberrated circular pupil reproduces
    /// [`mtf_ideal`](crate::diffraction::mtf_ideal) — the closed-form autocorrelation
    /// of a disc — without either calculation knowing about the other.
    pub fn mtf(&self) -> Mtf {
        let n = self.samples;
        // Undo the shift so the transform sees the PSF with its centre at the origin.
        let unshifted = fftshift(&self.intensity, n);
        let mut re = unshifted;
        let mut im = vec![0.0; n * n];
        fft2(&mut re, &mut im, n, false);
        let dc = (re[0] * re[0] + im[0] * im[0]).sqrt();
        let norm = if dc > 0.0 { 1.0 / dc } else { 0.0 };
        let mut magnitude = vec![0.0; n * n];
        for i in 0..n * n {
            magnitude[i] = (re[i] * re[i] + im[i] * im[i]).sqrt() * norm;
        }
        Mtf {
            samples: n,
            magnitude: fftshift(&magnitude, n),
            // The cutoff is one pupil diameter of separation in the autocorrelation,
            // which lands at `d` samples of the transform.
            cutoff_samples: n as f64 / self.scale,
        }
    }
}

/// A sampled modulation transfer function.
#[derive(Clone, Debug)]
pub struct Mtf {
    samples: usize,
    magnitude: Vec<f64>,
    /// Transform samples corresponding to the incoherent cutoff frequency.
    cutoff_samples: f64,
}

impl Mtf {
    /// Transfer along the x axis at a fraction of the cutoff frequency, interpolated.
    pub fn at_fraction(&self, s: f64) -> f64 {
        if s < 0.0 {
            return self.at_fraction(-s);
        }
        if s > 1.0 {
            return 0.0;
        }
        let px = s * self.cutoff_samples;
        let c = self.samples / 2;
        let lo = px.floor() as usize;
        let frac = px - lo as f64;
        let read = |i: usize| {
            if c + i < self.samples {
                self.magnitude[c * self.samples + c + i]
            } else {
                0.0
            }
        };
        read(lo) * (1.0 - frac) + read(lo + 1) * frac
    }

    /// Transfer averaged over direction, which is what a rotationally symmetric
    /// system should have and an astigmatic one should not.
    pub fn azimuthal_at_fraction(&self, s: f64) -> f64 {
        let px = (s * self.cutoff_samples).max(0.0);
        let c = (self.samples / 2) as f64;
        const STEPS: usize = 64;
        let mut sum = 0.0;
        for i in 0..STEPS {
            let a = TAU * i as f64 / STEPS as f64;
            let x = (c + px * a.cos()).round() as isize;
            let y = (c + px * a.sin()).round() as isize;
            if x >= 0 && y >= 0 && (x as usize) < self.samples && (y as usize) < self.samples {
                sum += self.magnitude[y as usize * self.samples + x as usize];
            }
        }
        sum / STEPS as f64
    }
}

// ---------------------------------------------------------------------------
// Transform
// ---------------------------------------------------------------------------

/// In-place radix-2 Cooley-Tukey, on split real and imaginary arrays.
///
/// The twiddle factors are computed from `cos` and `sin` at each step rather than
/// accumulated by repeated complex multiplication. Accumulating is faster and loses
/// several digits by the end of a long transform; this is a physics library, the
/// transforms here are small, and a PSF that is wrong in its fourth digit is not
/// worth the cycles saved. It also keeps the result independent of how the loop was
/// ordered, which matters for the same reason every other reduction in the workspace
/// is written in a fixed order.
fn fft_1d(re: &mut [f64], im: &mut [f64], inverse: bool) {
    let n = re.len();
    debug_assert!(n.is_power_of_two());

    // Bit-reversal permutation.
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }

    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2;
    while len <= n {
        let step = sign * TAU / len as f64;
        let mut base = 0;
        while base < n {
            for k in 0..len / 2 {
                let angle = step * k as f64;
                let (wr, wi) = (angle.cos(), angle.sin());
                let (a, b) = (base + k, base + k + len / 2);
                let (ur, ui) = (re[a], im[a]);
                let vr = re[b] * wr - im[b] * wi;
                let vi = re[b] * wi + im[b] * wr;
                re[a] = ur + vr;
                im[a] = ui + vi;
                re[b] = ur - vr;
                im[b] = ui - vi;
            }
            base += len;
        }
        len <<= 1;
    }

    if inverse {
        let scale = 1.0 / n as f64;
        for v in re.iter_mut() {
            *v *= scale;
        }
        for v in im.iter_mut() {
            *v *= scale;
        }
    }
}

/// Two-dimensional transform: rows, then columns.
///
/// Shared with [`propagation`](crate::propagation), which needs the same transform for
/// the angular spectrum. One implementation rather than two means the reversibility and
/// energy checks in either module cover both.
pub(crate) fn fft2(re: &mut [f64], im: &mut [f64], n: usize, inverse: bool) {
    let mut row_re = vec![0.0; n];
    let mut row_im = vec![0.0; n];
    for y in 0..n {
        row_re.copy_from_slice(&re[y * n..(y + 1) * n]);
        row_im.copy_from_slice(&im[y * n..(y + 1) * n]);
        fft_1d(&mut row_re, &mut row_im, inverse);
        re[y * n..(y + 1) * n].copy_from_slice(&row_re);
        im[y * n..(y + 1) * n].copy_from_slice(&row_im);
    }
    for x in 0..n {
        for y in 0..n {
            row_re[y] = re[y * n + x];
            row_im[y] = im[y * n + x];
        }
        fft_1d(&mut row_re, &mut row_im, inverse);
        for y in 0..n {
            re[y * n + x] = row_re[y];
            im[y * n + x] = row_im[y];
        }
    }
}

/// Swap quadrants so the zero frequency sits at the centre. Its own inverse for an
/// even-sized grid, which is why [`Psf::mtf`] uses it to undo itself.
pub(crate) fn fftshift(data: &[f64], n: usize) -> Vec<f64> {
    let half = n / 2;
    let mut out = vec![0.0; n * n];
    for y in 0..n {
        for x in 0..n {
            let sx = (x + half) % n;
            let sy = (y + half) % n;
            out[y * n + x] = data[sy * n + sx];
        }
    }
    out
}

/// Radius of the first Airy zero in units of `λ/D`: `3.8317/π`, or **1.2197**.
///
/// Not 0.61. That is the same zero measured in `λ/NA`, which is what
/// [`airy_radius`](crate::diffraction::airy_radius) returns, and the factor of two
/// between them is the numerical aperture: `NA = sin θ ≈ D/2f`. Mixing the two is the
/// most common way a resolution figure comes out wrong, and it is worth a named
/// constant on each side rather than a factor remembered at the call site.
pub const AIRY_ZERO_LAMBDA_OVER_D: f64 = FIRST_AIRY_ZERO / PI;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diffraction::{airy_intensity, encircled_energy, mtf_ideal};

    /// Orthonormality, which is the property the whole coefficient convention rests
    /// on: each mode has RMS one over the disc, and distinct modes are uncorrelated.
    ///
    /// Integrated numerically over a fine polar grid, and checked against the exact
    /// values 1 and 0 rather than against another implementation.
    #[test]
    fn the_zernike_modes_are_orthonormal_over_the_disc() {
        let modes = [
            Zernike::DEFOCUS,
            Zernike::ASTIGMATISM_0,
            Zernike::ASTIGMATISM_45,
            Zernike::COMA_X,
            Zernike::COMA_Y,
            Zernike::TREFOIL_0,
            Zernike::SPHERICAL,
        ];
        const NR: usize = 400;
        const NT: usize = 400;
        let inner = |a: Zernike, b: Zernike| {
            let mut sum = 0.0;
            let mut area = 0.0;
            for i in 0..NR {
                let rho = (i as f64 + 0.5) / NR as f64;
                for k in 0..NT {
                    let theta = TAU * (k as f64 + 0.5) / NT as f64;
                    let w = rho; // r dr dtheta
                    sum += a.value(rho, theta) * b.value(rho, theta) * w;
                    area += w;
                }
            }
            sum / area
        };
        for a in modes {
            let norm = inner(a, a);
            assert!(
                (norm - 1.0).abs() < 2e-3,
                "{} should have RMS one, got {norm}",
                a.name()
            );
            for b in modes {
                if a != b {
                    assert!(
                        inner(a, b).abs() < 2e-3,
                        "{} and {} should be orthogonal, got {}",
                        a.name(),
                        b.name(),
                        inner(a, b)
                    );
                }
            }
        }
    }

    /// Invalid modes are rejected rather than silently returning something.
    #[test]
    fn only_real_modes_exist() {
        assert!(Zernike::DEFOCUS.is_valid());
        assert!(Zernike::new(4, 2).is_valid());
        // |m| > n
        assert!(!Zernike::new(2, 3).is_valid());
        // n - |m| odd
        assert!(!Zernike::new(3, 0).is_valid());
        assert_eq!(Zernike::new(3, 0).value(0.5, 0.3), 0.0);
        // Outside the disc there is no polynomial.
        assert_eq!(Zernike::DEFOCUS.value(1.5, 0.0), 0.0);
    }

    /// A perfect pupil's PSF is the Airy pattern. Two entirely separate calculations —
    /// a Bessel series in `diffraction`, a fast Fourier transform here — and they
    /// agree on the core, the first zero and the first ring.
    ///
    /// This is the test that makes the rest of the module trustworthy.
    #[test]
    fn an_unaberrated_pupil_reproduces_the_airy_pattern() {
        let psf = Pupil::circular(512, 64.0).psf();
        // N/d = 8 pixels per lambda/D, so the first zero is at 4.9 px.
        assert!((psf.scale() - 8.0).abs() < 1e-12);
        assert!(
            (psf.strehl() - 1.0).abs() < 1e-9,
            "a perfect pupil peaks at 1"
        );

        // Both sides are averaged over the *same* pixels. Comparing a bin mean against
        // the closed form at the bin's centre would be measuring the binning of a
        // steep function, not the transform: at eight pixels per lambda/D the core
        // falls by half across one pixel, and that alone is worth a few percent.
        const BINS: usize = 30;
        let c = (psf.samples() / 2) as isize;
        let limit_px = 3.0 * psf.scale();
        let mut measured = vec![0.0; BINS];
        let mut analytic = vec![0.0; BINS];
        let mut counts = vec![0usize; BINS];
        for y in 0..psf.samples() as isize {
            for x in 0..psf.samples() as isize {
                let (dx, dy) = ((x - c) as f64, (y - c) as f64);
                let r_px = (dx * dx + dy * dy).sqrt();
                if r_px > limit_px {
                    continue;
                }
                let bin = (((r_px / limit_px) * BINS as f64) as usize).min(BINS - 1);
                measured[bin] += psf.at(x - c, y - c);
                // v = pi r, with r in lambda/D.
                analytic[bin] += airy_intensity(PI * r_px / psf.scale());
                counts[bin] += 1;
            }
        }
        let mut worst = 0.0f64;
        for bin in 0..BINS {
            if counts[bin] == 0 {
                continue;
            }
            let (m, a) = (
                measured[bin] / counts[bin] as f64,
                analytic[bin] / counts[bin] as f64,
            );
            let r = (bin as f64 + 0.5) / BINS as f64 * 3.0;
            assert!(
                (m - a).abs() < 5e-3,
                "at r = {r:.2} lambda/D the transform gives {m:.4} and the Bessel form \
                 gives {a:.4}"
            );
            worst = worst.max((m - a).abs());
        }
        // Well under a percent of the peak, everywhere out to three lambda/D. Two
        // calculations sharing no code.
        assert!(worst < 5e-3, "worst deviation {worst:e}");

        // And the first dark ring is dark, at the radius the closed form puts it. How
        // much energy sits inside it is checked separately, against `encircled_energy`.
        let ring_px = (AIRY_ZERO_LAMBDA_OVER_D * psf.scale()).round() as isize;
        assert!(
            psf.at(ring_px, 0) < 0.005,
            "the first zero at {ring_px} px should be dark, got {}",
            psf.at(ring_px, 0)
        );
    }

    /// Energy is conserved by the transform, which is Parseval's theorem and the
    /// cheapest possible check that the FFT is not quietly wrong.
    #[test]
    fn the_transform_conserves_energy() {
        let pupil = Pupil::circular(256, 32.0);
        let psf = pupil.psf();
        let n = pupil.samples() as f64;
        let sum_a: f64 = pupil.amplitude.iter().sum();
        let sum_a2: f64 = pupil.amplitude.iter().map(|a| a * a).sum();
        // Under this normalisation Parseval reads N^2 * sum(a^2) / sum(a)^2.
        let expected = n * n * sum_a2 / (sum_a * sum_a);
        assert!(
            (psf.total_energy() / expected - 1.0).abs() < 1e-9,
            "got {}, expected {expected}",
            psf.total_energy()
        );
        // And an aberration moves light around without creating any: the same pupil
        // with half a wave of coma has the same total.
        let aberrated = Pupil::circular(256, 32.0)
            .with_zernike(Zernike::COMA_X, 0.5)
            .psf();
        assert!(
            (aberrated.total_energy() / psf.total_energy() - 1.0).abs() < 1e-9,
            "an aberration is not a loss"
        );
    }

    /// 84% of a perfect image's light is inside the first dark ring, which
    /// `diffraction` computes from `1 - J0² - J1²`. The transform agrees, a little
    /// high because the grid truncates the wings.
    #[test]
    fn encircled_energy_agrees_with_the_closed_form() {
        let psf = Pupil::circular(512, 64.0).psf();
        let analytic = encircled_energy(FIRST_AIRY_ZERO);
        let measured = psf.encircled_energy(AIRY_ZERO_LAMBDA_OVER_D);
        assert!(
            (analytic - 0.8378).abs() < 1e-3,
            "the closed form should be 83.8%"
        );
        assert!(
            measured > analytic && measured - analytic < 0.05,
            "the grid truncates the wings, so {measured:.4} should sit a little above \
             {analytic:.4}"
        );
    }

    /// The MTF of a perfect circular pupil, from the transform of its PSF, against
    /// the closed-form autocorrelation of a disc. A third independent agreement, and
    /// the one that ties the image-plane calculation back to the frequency domain.
    #[test]
    fn the_mtf_reproduces_the_disc_autocorrelation() {
        let mtf = Pupil::circular(512, 64.0).psf().mtf();
        assert!(
            (mtf.at_fraction(0.0) - 1.0).abs() < 1e-9,
            "unity at zero frequency"
        );
        for s in [0.1, 0.25, 0.5, 0.75, 0.9] {
            let measured = mtf.at_fraction(s);
            let analytic = mtf_ideal(s);
            assert!(
                (measured - analytic).abs() < 0.02,
                "at {s} of cutoff the transform gives {measured:.4} and the closed \
                 form gives {analytic:.4}"
            );
        }
        // 0.391 at half cutoff, the figure every MTF chart is read against.
        assert!((mtf.at_fraction(0.5) - 0.391).abs() < 0.02);
        // Nothing past the cutoff.
        assert_eq!(mtf.at_fraction(1.2), 0.0);
    }

    /// Small aberrations follow the Maréchal approximation, which is what
    /// `diffraction::strehl_from_wavefront_error` computes. Checked across modes,
    /// because Maréchal depends only on the RMS and not on which aberration produced
    /// it — a claim worth testing rather than repeating.
    #[test]
    fn small_aberrations_follow_the_marechal_approximation() {
        use crate::diffraction::strehl_from_wavefront_error;
        for mode in [
            Zernike::DEFOCUS,
            Zernike::ASTIGMATISM_0,
            Zernike::COMA_X,
            Zernike::SPHERICAL,
            Zernike::TREFOIL_0,
        ] {
            for rms in [0.02, 0.05, 0.075] {
                let pupil = Pupil::circular(256, 64.0).with_zernike(mode, rms);
                // The pupil's own RMS should be the coefficient, since the modes are
                // orthonormal and piston is removed.
                let measured_rms = pupil.rms_wavefront_error();
                assert!(
                    (measured_rms / rms - 1.0).abs() < 0.03,
                    "{}: asked for {rms} waves RMS, pupil has {measured_rms:.4}",
                    mode.name()
                );

                let strehl = pupil.psf().strehl();
                let predicted = strehl_from_wavefront_error(rms);
                assert!(
                    (strehl - predicted).abs() < 0.02,
                    "{} at {rms} waves: transform gives Strehl {strehl:.4}, Maréchal \
                     predicts {predicted:.4}",
                    mode.name()
                );
            }
        }
    }

    /// The diffraction-limited threshold, arrived at from the pupil instead of from
    /// the approximation: a fourteenth of a wave RMS leaves a Strehl just past 0.8.
    #[test]
    fn a_fourteenth_of_a_wave_is_still_diffraction_limited() {
        let strehl = Pupil::circular(256, 64.0)
            .with_zernike(Zernike::COMA_X, 1.0 / 14.0)
            .psf()
            .strehl();
        assert!(
            strehl > 0.80 && strehl < 0.84,
            "lambda/14 RMS should sit just past the 0.8 threshold, got {strehl:.4}"
        );
        // A quarter wave RMS does not.
        let broken = Pupil::circular(256, 64.0)
            .with_zernike(Zernike::COMA_X, 0.25)
            .psf()
            .strehl();
        assert!(
            broken < 0.3,
            "a quarter wave RMS is badly aberrated, got {broken:.4}"
        );
    }

    /// What the geometric half of the crate could not express, and the reason this
    /// module exists: aberrations degrade an image in ways that differ from each
    /// other, and the difference is visible rather than a matter of degree.
    ///
    /// Defocus is symmetric and puts light into a ring. Coma is not symmetric and
    /// drags the brightest point off axis, which is why it cannot be focused out and
    /// why it ruins astrometry rather than merely softening a picture.
    #[test]
    fn coma_moves_the_image_and_defocus_does_not() {
        let defocused = Pupil::circular(512, 64.0)
            .with_zernike(Zernike::DEFOCUS, 0.2)
            .psf();
        let comatic = Pupil::circular(512, 64.0)
            .with_zernike(Zernike::COMA_X, 0.2)
            .psf();

        // Defocus stays symmetric: the same intensity either side of centre.
        let left = defocused.at(-6, 0);
        let right = defocused.at(6, 0);
        assert!(
            (left - right).abs() / left.max(right) < 1e-6,
            "defocus should stay symmetric: {left:.6} vs {right:.6}"
        );

        // Coma does not.
        let cl = comatic.at(-6, 0);
        let cr = comatic.at(6, 0);
        assert!(
            (cl - cr).abs() / cl.max(cr) > 0.2,
            "coma should be lopsided: {cl:.6} vs {cr:.6}"
        );
        // Both lose about the same amount on axis, since Maréchal only sees the RMS.
        assert!((defocused.strehl() - comatic.strehl()).abs() < 0.05);
        // But coma's brightest point is brighter than its on-axis value, because the
        // peak has moved off centre.
        assert!(
            comatic.peak() > comatic.strehl() * 1.001,
            "coma's peak {} should have left the axis (on-axis {})",
            comatic.peak(),
            comatic.strehl()
        );
        assert!(
            (defocused.peak() - defocused.strehl()).abs() < 1e-9,
            "defocus keeps its peak on axis"
        );
    }

    /// A central obscuration narrows the core and fills in the rings. Both halves of
    /// that trade are real, and a design that quotes only the first half is quoting
    /// half of it: a Cassegrain resolves finer detail and shows worse contrast on an
    /// extended scene, which is why one is used for stars and not for landscapes.
    #[test]
    fn an_obscuration_buys_resolution_with_contrast() {
        let filled = Pupil::circular(512, 64.0).psf();
        let annular = Pupil::annular(512, 64.0, 0.4).psf();

        // The core narrows: the annular PSF is darker at the filled pupil's half-power
        // radius because more of its light has left the core.
        let core_energy_filled = filled.encircled_energy(AIRY_ZERO_LAMBDA_OVER_D);
        let core_energy_annular = annular.encircled_energy(AIRY_ZERO_LAMBDA_OVER_D);
        assert!(
            core_energy_annular < core_energy_filled,
            "the obscuration should move energy out of the core: {core_energy_annular:.4} \
             vs {core_energy_filled:.4}"
        );

        // And the mid frequencies suffer while the high ones do not.
        let mtf_filled = filled.mtf();
        let mtf_annular = annular.mtf();
        assert!(
            mtf_annular.at_fraction(0.3) < mtf_filled.at_fraction(0.3),
            "mid-frequency contrast should fall"
        );
        assert!(
            mtf_annular.at_fraction(0.85) > mtf_filled.at_fraction(0.85),
            "high-frequency transfer should rise"
        );
        // Both still start at one and stop at the same cutoff, which the obscuration
        // does not move.
        assert!((mtf_annular.at_fraction(0.0) - 1.0).abs() < 1e-9);
        assert_eq!(mtf_annular.at_fraction(1.05), 0.0);
    }

    /// Astigmatism needs defocus to show which axis it prefers, and that is the
    /// physics rather than a limitation.
    ///
    /// At best focus the two axes are *equally* bad. Swapping x and y negates
    /// `ρ²cos2θ`, and negating a wavefront reflects its PSF through the origin, so the
    /// two operations cancel and the pattern comes out symmetric — a symmetric cross,
    /// not a line. Add defocus and the symmetry breaks, because defocus is
    /// rotationally symmetric and does not negate: one axis comes to a line focus
    /// while the other spreads. Reversing the defocus swaps them.
    ///
    /// That pair of line foci either side of best focus is the sagittal and tangential
    /// split a lens chart reports separately, and it is why refocusing an astigmatic
    /// system trades one blur for the other instead of removing it.
    #[test]
    fn astigmatism_needs_defocus_to_pick_an_axis() {
        let best_focus = Pupil::circular(512, 64.0)
            .with_zernike(Zernike::ASTIGMATISM_0, 0.25)
            .psf();
        let (x0, y0) = (best_focus.at(8, 0), best_focus.at(0, 8));
        assert!(
            (x0 - y0).abs() / x0.max(y0) < 1e-9,
            "at best focus astigmatism is equally bad both ways: {x0:.6} vs {y0:.6}"
        );

        let astigmatic = |defocus: f64| {
            Pupil::circular(512, 64.0)
                .with_zernike(Zernike::ASTIGMATISM_0, 0.25)
                .with_zernike(Zernike::DEFOCUS, defocus)
                .psf()
        };
        let near = astigmatic(0.25);
        let (nx, ny) = (near.at(8, 0), near.at(0, 8));
        assert!(
            (nx - ny).abs() / nx.max(ny) > 0.3,
            "defocused astigmatism should prefer an axis: {nx:.6} vs {ny:.6}"
        );

        // The other side of focus is the same pattern turned through ninety degrees.
        let far = astigmatic(-0.25);
        assert!(
            (far.at(0, 8) - nx).abs() / nx < 1e-9 && (far.at(8, 0) - ny).abs() / ny < 1e-9,
            "reversing the defocus should swap the axes"
        );

        // A rotationally symmetric aberration never prefers an axis, at any focus.
        let round = Pupil::circular(512, 64.0)
            .with_zernike(Zernike::SPHERICAL, 0.25)
            .with_zernike(Zernike::DEFOCUS, 0.25)
            .psf();
        assert!(
            (round.at(8, 0) - round.at(0, 8)).abs() / round.at(8, 0) < 1e-9,
            "spherical aberration is symmetric with or without defocus"
        );
    }

    /// Piston and tilt are not aberrations in the sense that matters, and the module
    /// treats them accordingly: piston is removed from the RMS entirely, and tilt
    /// moves the image without dimming it.
    #[test]
    fn piston_is_invisible_and_tilt_only_moves_things() {
        let piston = Pupil::circular(256, 64.0).with_zernike(Zernike::PISTON, 0.3);
        assert!(
            piston.rms_wavefront_error() < 1e-12,
            "a constant phase is not a wavefront error"
        );
        assert!((piston.psf().strehl() - 1.0).abs() < 1e-9);

        // Tilt keeps all its light, and keeps its peak, but not on axis.
        let tilted = Pupil::circular(256, 64.0)
            .with_zernike(Zernike::TILT_X, 0.3)
            .psf();
        let perfect = Pupil::circular(256, 64.0).psf();
        assert!(
            (tilted.peak() - perfect.peak()).abs() < 0.02,
            "tilt should not dim the image: {} vs {}",
            tilted.peak(),
            perfect.peak()
        );
        assert!(
            tilted.strehl() < tilted.peak() * 0.999,
            "but it should have moved off axis"
        );
    }

    /// The transform is deterministic, like everything else here. Same pupil, same
    /// bits, and no dependence on how the loops were ordered.
    #[test]
    fn the_transform_is_bit_reproducible() {
        let build = || {
            Pupil::circular(128, 32.0)
                .with_zernike(Zernike::COMA_X, 0.13)
                .with_zernike(Zernike::SPHERICAL, -0.07)
                .psf()
        };
        let a = build();
        let b = build();
        assert_eq!(a.intensity.len(), b.intensity.len());
        for (x, y) in a.intensity.iter().zip(b.intensity.iter()) {
            assert_eq!(x.to_bits(), y.to_bits(), "not bit-identical");
        }
        assert_eq!(a.strehl().to_bits(), b.strehl().to_bits());
    }

    /// The forward transform followed by the inverse returns what went in, which
    /// catches a sign or scaling error the magnitude-only tests would not.
    #[test]
    fn the_inverse_transform_undoes_the_forward_one() {
        const N: usize = 64;
        let mut re: Vec<f64> = (0..N * N)
            .map(|i| ((i * 37) % 101) as f64 / 101.0)
            .collect();
        let mut im: Vec<f64> = (0..N * N).map(|i| ((i * 53) % 97) as f64 / 97.0).collect();
        let (re0, im0) = (re.clone(), im.clone());
        fft2(&mut re, &mut im, N, false);
        fft2(&mut re, &mut im, N, true);
        for i in 0..N * N {
            assert!((re[i] - re0[i]).abs() < 1e-12, "real part at {i}");
            assert!((im[i] - im0[i]).abs() < 1e-12, "imaginary part at {i}");
        }
    }

    /// Grids that the transform cannot handle are refused at construction rather than
    /// producing a wrong answer.
    #[test]
    #[should_panic(expected = "power of two")]
    fn a_non_power_of_two_grid_is_refused() {
        Pupil::circular(100, 32.0);
    }

    #[test]
    #[should_panic(expected = "must fit")]
    fn an_aperture_larger_than_its_grid_is_refused() {
        Pupil::circular(64, 128.0);
    }
}
