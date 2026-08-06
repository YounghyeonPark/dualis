//! The discrete Fourier transform, as a kernel utility.
//!
//! Numerics rather than physics, which is why it belongs here alongside
//! [`integrator`](crate::integrator) and [`vector`](crate::vector) rather than in any
//! domain. It arrived in `dualis-optics` because that is where it was first needed —
//! twice, for pupil transforms and for angular-spectrum propagation — and it stayed
//! there one domain too long. An electrostatic solver needs the same transform for
//! Ewald summation, and a domain crate cannot reach into another one.
//!
//! # Accuracy over speed, deliberately
//!
//! The twiddle factors are computed from `cos` and `sin` at every butterfly rather than
//! accumulated by repeated complex multiplication. Accumulating is faster and loses
//! several digits by the end of a long transform. This is a physics library: a result
//! wrong in its fourth digit is not worth the cycles saved, and the transforms here are
//! small.
//!
//! It also keeps the answer independent of how the loops were ordered, which is the same
//! reason every other reduction in this workspace is written in a fixed sequence.
//!
//! # Radix two only
//!
//! Lengths must be powers of two, and that is enforced rather than worked around. A
//! mixed-radix implementation would be several times the code for a case nothing here
//! has needed, and silently padding a caller's array would change the frequency grid
//! underneath them.
//!
//! # Sign convention
//!
//! The forward transform carries `exp(-2πi jk/N)` and the inverse carries `exp(+2πi
//! jk/N)` with a `1/N` in front, so that `ifft(fft(x)) == x`. That is the convention
//! physics and signal processing agree on; some numerical libraries put the `1/N` on the
//! forward transform instead, and a field propagated with one and read back with the
//! other comes out scaled by `N²`.

/// In-place one-dimensional transform on split real and imaginary parts.
pub fn fft(re: &mut [f64], im: &mut [f64]) {
    transform_1d(re, im, false);
}

/// In-place one-dimensional inverse transform, including the `1/N`.
pub fn ifft(re: &mut [f64], im: &mut [f64]) {
    transform_1d(re, im, true);
}

/// In-place two-dimensional transform of an `n`-by-`n` array in row-major order.
pub fn fft2(re: &mut [f64], im: &mut [f64], n: usize) {
    transform_2d(re, im, n, false);
}

/// In-place two-dimensional inverse transform, including the `1/N²`.
pub fn ifft2(re: &mut [f64], im: &mut [f64], n: usize) {
    transform_2d(re, im, n, true);
}

/// Swap quadrants so the zero frequency sits at the centre of the array.
///
/// Its own inverse for an even-sized grid, which is why a caller can use it both to
/// centre a result for display and to un-centre one before transforming it again.
pub fn fftshift(data: &[f64], n: usize) -> Vec<f64> {
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

fn transform_1d(re: &mut [f64], im: &mut [f64], inverse: bool) {
    let n = re.len();
    assert_eq!(
        n,
        im.len(),
        "the real and imaginary parts must be the same length"
    );
    assert!(
        n.is_power_of_two(),
        "the transform is radix two, so the length must be a power of two, got {n}"
    );
    if n < 2 {
        return;
    }

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
        let step = sign * std::f64::consts::TAU / len as f64;
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

/// Rows, then columns. Separable, so the order does not matter mathematically — it is
/// fixed here anyway, because the floating-point result does depend on it.
fn transform_2d(re: &mut [f64], im: &mut [f64], n: usize, inverse: bool) {
    assert_eq!(re.len(), n * n, "the array must be n by n");
    assert_eq!(im.len(), n * n, "the array must be n by n");
    let mut row_re = vec![0.0; n];
    let mut row_im = vec![0.0; n];
    for y in 0..n {
        row_re.copy_from_slice(&re[y * n..(y + 1) * n]);
        row_im.copy_from_slice(&im[y * n..(y + 1) * n]);
        transform_1d(&mut row_re, &mut row_im, inverse);
        re[y * n..(y + 1) * n].copy_from_slice(&row_re);
        im[y * n..(y + 1) * n].copy_from_slice(&row_im);
    }
    for x in 0..n {
        for y in 0..n {
            row_re[y] = re[y * n + x];
            row_im[y] = im[y * n + x];
        }
        transform_1d(&mut row_re, &mut row_im, inverse);
        for y in 0..n {
            re[y * n + x] = row_re[y];
            im[y * n + x] = row_im[y];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    /// A single complex exponential lands in exactly one bin, and the rest are empty.
    ///
    /// The sharpest available check on a transform: it pins the sign convention, the
    /// normalisation and the bin ordering at once, and any error in the butterflies
    /// smears energy into neighbouring bins where there should be none.
    #[test]
    fn one_frequency_lands_in_one_bin() {
        const N: usize = 64;
        for m in [0usize, 1, 7, 31, 63] {
            let mut re: Vec<f64> = (0..N)
                .map(|k| (TAU * m as f64 * k as f64 / N as f64).cos())
                .collect();
            let mut im: Vec<f64> = (0..N)
                .map(|k| (TAU * m as f64 * k as f64 / N as f64).sin())
                .collect();
            fft(&mut re, &mut im);
            for j in 0..N {
                let magnitude = (re[j] * re[j] + im[j] * im[j]).sqrt();
                let expected = if j == m { N as f64 } else { 0.0 };
                assert!(
                    (magnitude - expected).abs() < 1e-10,
                    "frequency {m}: bin {j} has {magnitude}, expected {expected}"
                );
            }
        }
    }

    /// A delta is flat, which is the same statement read the other way round.
    #[test]
    fn a_delta_transforms_to_a_flat_spectrum() {
        const N: usize = 32;
        let mut re = vec![0.0; N];
        let mut im = vec![0.0; N];
        re[0] = 1.0;
        fft(&mut re, &mut im);
        for j in 0..N {
            assert!((re[j] - 1.0).abs() < 1e-12, "bin {j} real part {}", re[j]);
            assert!(im[j].abs() < 1e-12, "bin {j} imaginary part {}", im[j]);
        }
    }

    /// The inverse undoes the forward, in one dimension and in two.
    ///
    /// Catches a sign or scaling error that no magnitude-only test would: getting the
    /// `1/N` onto the wrong transform, or both, still gives plausible spectra.
    #[test]
    fn the_inverse_undoes_the_forward() {
        const N: usize = 64;
        let re0: Vec<f64> = (0..N).map(|i| ((i * 37) % 101) as f64 / 101.0).collect();
        let im0: Vec<f64> = (0..N).map(|i| ((i * 53) % 97) as f64 / 97.0).collect();

        let (mut re, mut im) = (re0.clone(), im0.clone());
        fft(&mut re, &mut im);
        ifft(&mut re, &mut im);
        for i in 0..N {
            assert!((re[i] - re0[i]).abs() < 1e-12, "1D real at {i}");
            assert!((im[i] - im0[i]).abs() < 1e-12, "1D imaginary at {i}");
        }

        const M: usize = 32;
        let re0: Vec<f64> = (0..M * M).map(|i| ((i * 29) % 89) as f64 / 89.0).collect();
        let im0: Vec<f64> = (0..M * M).map(|i| ((i * 41) % 83) as f64 / 83.0).collect();
        let (mut re, mut im) = (re0.clone(), im0.clone());
        fft2(&mut re, &mut im, M);
        ifft2(&mut re, &mut im, M);
        for i in 0..M * M {
            assert!((re[i] - re0[i]).abs() < 1e-12, "2D real at {i}");
            assert!((im[i] - im0[i]).abs() < 1e-12, "2D imaginary at {i}");
        }
    }

    /// Parseval: the transform moves energy around without creating or destroying any.
    /// `Σ|X|² = N Σ|x|²` in one dimension and `N² Σ|x|²` in two.
    #[test]
    fn the_transform_conserves_energy() {
        const N: usize = 128;
        let re0: Vec<f64> = (0..N)
            .map(|i| ((i * 17) % 71) as f64 / 71.0 - 0.5)
            .collect();
        let im0: Vec<f64> = (0..N)
            .map(|i| ((i * 23) % 59) as f64 / 59.0 - 0.5)
            .collect();
        let before: f64 = re0.iter().zip(im0.iter()).map(|(a, b)| a * a + b * b).sum();

        let (mut re, mut im) = (re0.clone(), im0.clone());
        fft(&mut re, &mut im);
        let after: f64 = re.iter().zip(im.iter()).map(|(a, b)| a * a + b * b).sum();
        assert!(
            (after / (N as f64 * before) - 1.0).abs() < 1e-12,
            "1D: {after} against {} ",
            N as f64 * before
        );

        const M: usize = 32;
        let re0: Vec<f64> = (0..M * M).map(|i| ((i * 13) % 61) as f64 / 61.0).collect();
        let im0 = vec![0.0; M * M];
        let before: f64 = re0.iter().map(|a| a * a).sum();
        let (mut re, mut im) = (re0, im0);
        fft2(&mut re, &mut im, M);
        let after: f64 = re.iter().zip(im.iter()).map(|(a, b)| a * a + b * b).sum();
        let n2 = (M * M) as f64;
        assert!((after / (n2 * before) - 1.0).abs() < 1e-12, "2D");
    }

    /// The two-dimensional transform is separable, and this checks that the
    /// rows-then-columns implementation really is it: the transform of an outer product
    /// is the outer product of the transforms.
    ///
    /// A strong test, because getting the row and column strides confused produces
    /// something that still looks like a spectrum.
    #[test]
    fn the_two_dimensional_transform_is_separable() {
        const N: usize = 16;
        let a: Vec<f64> = (0..N).map(|i| ((i * 7) % 13) as f64).collect();
        let b: Vec<f64> = (0..N).map(|i| ((i * 5) % 11) as f64).collect();

        // Transform the two one-dimensional signals on their own.
        let (mut ar, mut ai) = (a.clone(), vec![0.0; N]);
        let (mut br, mut bi) = (b.clone(), vec![0.0; N]);
        fft(&mut ar, &mut ai);
        fft(&mut br, &mut bi);

        // And the outer product as a two-dimensional array.
        let mut re = vec![0.0; N * N];
        let im0 = vec![0.0; N * N];
        for y in 0..N {
            for x in 0..N {
                re[y * N + x] = a[x] * b[y];
            }
        }
        let mut im = im0;
        fft2(&mut re, &mut im, N);

        for v in 0..N {
            for u in 0..N {
                // (Ar + i Ai)(Br + i Bi)
                let want_re = ar[u] * br[v] - ai[u] * bi[v];
                let want_im = ar[u] * bi[v] + ai[u] * br[v];
                let i = v * N + u;
                assert!(
                    (re[i] - want_re).abs() < 1e-9 && (im[i] - want_im).abs() < 1e-9,
                    "at ({u},{v}): got ({}, {}), expected ({want_re}, {want_im})",
                    re[i],
                    im[i]
                );
            }
        }
    }

    /// A real, even signal has a real spectrum. A symmetry the butterflies must preserve
    /// exactly, and one that a sign error in the twiddle factors breaks.
    #[test]
    fn a_real_even_signal_has_a_real_spectrum() {
        const N: usize = 64;
        let mut re: Vec<f64> = (0..N)
            .map(|i| {
                let k = if i <= N / 2 { i } else { N - i };
                (-(k as f64) * 0.3).exp()
            })
            .collect();
        let mut im = vec![0.0; N];
        fft(&mut re, &mut im);
        for (j, imaginary) in im.iter().enumerate() {
            assert!(
                imaginary.abs() < 1e-12,
                "bin {j} should be real, imaginary part {imaginary}"
            );
        }
    }

    /// Linearity, which is cheap to check and rules out an accidental nonlinearity in the
    /// scaling.
    #[test]
    fn the_transform_is_linear() {
        const N: usize = 32;
        let x: Vec<f64> = (0..N).map(|i| ((i * 11) % 17) as f64).collect();
        let y: Vec<f64> = (0..N).map(|i| ((i * 3) % 7) as f64).collect();

        let (mut xr, mut xi) = (x.clone(), vec![0.0; N]);
        let (mut yr, mut yi) = (y.clone(), vec![0.0; N]);
        fft(&mut xr, &mut xi);
        fft(&mut yr, &mut yi);

        let mut sr: Vec<f64> = x
            .iter()
            .zip(y.iter())
            .map(|(a, b)| 2.0 * a + 3.0 * b)
            .collect();
        let mut si = vec![0.0; N];
        fft(&mut sr, &mut si);

        for j in 0..N {
            let want = 2.0 * xr[j] + 3.0 * yr[j];
            assert!(
                (sr[j] - want).abs() < 1e-9,
                "bin {j}: {} against {want}",
                sr[j]
            );
            let want = 2.0 * xi[j] + 3.0 * yi[j];
            assert!((si[j] - want).abs() < 1e-9);
        }
    }

    /// Shifting is its own inverse for an even grid, which is what lets a caller centre a
    /// result and un-centre it with the same call.
    #[test]
    fn shifting_twice_returns_the_original() {
        const N: usize = 8;
        let data: Vec<f64> = (0..N * N).map(|i| i as f64).collect();
        let once = fftshift(&data, N);
        let twice = fftshift(&once, N);
        assert_eq!(data, twice);
        // And it really moved something: the corner becomes the centre.
        assert_eq!(once[(N / 2) * N + N / 2], data[0]);
    }

    /// Bit-reproducible, like everything else in this workspace.
    #[test]
    fn the_transform_is_bit_reproducible() {
        const N: usize = 64;
        let build = || {
            let mut re: Vec<f64> = (0..N * N).map(|i| ((i * 19) % 43) as f64).collect();
            let mut im: Vec<f64> = (0..N * N).map(|i| ((i * 31) % 37) as f64).collect();
            fft2(&mut re, &mut im, N);
            (re, im)
        };
        let (a_re, a_im) = build();
        let (b_re, b_im) = build();
        for i in 0..N * N {
            assert_eq!(a_re[i].to_bits(), b_re[i].to_bits());
            assert_eq!(a_im[i].to_bits(), b_im[i].to_bits());
        }
    }

    /// Lengths the transform cannot handle are refused rather than padded, which would
    /// change the frequency grid under the caller.
    #[test]
    #[should_panic(expected = "power of two")]
    fn a_non_power_of_two_length_is_refused() {
        let mut re = vec![0.0; 12];
        let mut im = vec![0.0; 12];
        fft(&mut re, &mut im);
    }

    #[test]
    #[should_panic(expected = "same length")]
    fn mismatched_halves_are_refused() {
        let mut re = vec![0.0; 8];
        let mut im = vec![0.0; 4];
        fft(&mut re, &mut im);
    }

    /// A degenerate length is a no-op rather than an error: the transform of one sample
    /// is that sample.
    #[test]
    fn a_single_sample_transforms_to_itself() {
        let mut re = vec![3.0];
        let mut im = vec![-1.0];
        fft(&mut re, &mut im);
        assert_eq!(re, vec![3.0]);
        assert_eq!(im, vec![-1.0]);
    }
}
