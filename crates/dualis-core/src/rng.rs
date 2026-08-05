//! A deterministic pseudo-random generator, and the sampling built on it.
//!
//! Every stochastic choice in a simulation has to be reproducible or the results
//! cannot be compared, so this is a plain seeded generator with no global state
//! and no entropy source. Two runs of the same scene draw the same numbers, on
//! every platform and in WebAssembly.
//!
//! # Determinism survives parallelism, but only through [`Rng::for_index`]
//!
//! A single [`Rng`] threaded through a trace is reproducible only while the draws
//! happen in one order. Hand that generator to several threads and the order
//! becomes whatever the scheduler decided this time, so the run stops being
//! repeatable — the invariant is lost precisely when the simulation gets big
//! enough to need it.
//!
//! [`Rng::for_index`] is the way out. It hashes a `(seed, index)` pair into an
//! independent stream, statelessly, so ray 10 000 can be drawn before ray 3 and
//! neither result changes. Seed by whatever identifies the work — pixel, sample,
//! bounce, cell, particle — and the answer no longer depends on the order the
//! work happened to be done in.

use glam::DVec3;

use crate::vector::basis_for;

/// Deterministic xorshift64* PRNG. Avoids a dependency and keeps scene
/// generation identical on every platform (important for WASM + tests).
///
/// `Clone` is deliberate: cloning captures the exact stream position, which is
/// how a speculative draw can be replayed or a substream forked at a known point.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Rng(u64);

/// SplitMix64's finaliser: an avalanche mix that decorrelates inputs differing by
/// a single bit. Used to turn structured indices — a pixel number, a bounce depth
/// — into seeds that behave like independent ones.
fn mix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

impl Rng {
    /// A generator from a seed.
    ///
    /// Every seed but zero is used as given; xorshift is stuck at zero forever,
    /// so that one value is replaced.
    pub fn new(seed: u64) -> Self {
        Rng(if seed == 0 {
            0x9E37_79B9_7F4A_7C15
        } else {
            seed
        })
    }

    /// An independent stream for one piece of work, identified by index.
    ///
    /// Stateless and order-free: the stream for index 10 000 is what it is
    /// whether or not index 3 was ever drawn. That is what lets a trace be
    /// parallel and bit-reproducible at the same time, and it is why every
    /// stochastic loop should seed per item rather than share one generator.
    ///
    /// Both arguments are mixed, so `(seed, index)` and `(index, seed)` differ and
    /// adjacent indices do not produce correlated streams.
    pub fn for_index(seed: u64, index: u64) -> Rng {
        Rng::new(mix64(seed ^ mix64(index)))
    }

    /// Fork a child stream and advance this one past it.
    ///
    /// For nesting that has no natural index — a recursive bounce that needs its
    /// own sampling without disturbing the caller's sequence. Where an index
    /// exists, prefer [`Rng::for_index`]: splitting is still order-dependent.
    pub fn split(&mut self) -> Rng {
        let drawn = self.next_u64();
        Rng::new(mix64(drawn))
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform in [0, 1).
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    /// Uniform point inside a disc of the given radius.
    pub fn in_disc(&mut self, radius: f64) -> (f64, f64) {
        let r = radius * self.unit().sqrt();
        let phi = std::f64::consts::TAU * self.unit();
        (r * phi.cos(), r * phi.sin())
    }

    /// Uniform direction on the unit sphere.
    pub fn on_sphere(&mut self) -> DVec3 {
        let z = self.range(-1.0, 1.0);
        let phi = std::f64::consts::TAU * self.unit();
        let r = (1.0 - z * z).max(0.0).sqrt();
        DVec3::new(r * phi.cos(), r * phi.sin(), z)
    }

    /// Uniform direction on the hemisphere about `normal`.
    pub fn on_hemisphere(&mut self, normal: DVec3) -> DVec3 {
        let d = self.on_sphere();
        if d.dot(normal) < 0.0 {
            -d
        } else {
            d
        }
    }

    /// Cosine-weighted direction about `normal` — the Lambertian scatter.
    ///
    /// A matte surface does not spray light evenly over the hemisphere: it sends
    /// it in proportion to the cosine of the angle from the normal, which is why
    /// it looks equally bright from every direction. Sampled by Malley's method
    /// (a uniform disc lifted onto the hemisphere), so the cosine weight is in the
    /// distribution and the estimator needs no correction factor.
    pub fn cosine_hemisphere(&mut self, normal: DVec3) -> DVec3 {
        let n = normal.normalize();
        let (x, y) = self.in_disc(1.0);
        let z = (1.0 - x * x - y * y).max(0.0).sqrt();
        let (u, v) = basis_for(n);
        (u * x + v * y + n * z).normalize()
    }

    /// A standard normal deviate, mean 0 and variance 1.
    ///
    /// Read noise, mechanical jitter, thermal fluctuation, Brownian motion: the
    /// noise in a simulation is Gaussian far more often than it is uniform. Plain
    /// Box-Muller rather than the polar form, because it draws exactly two numbers
    /// every time — a rejection loop would make stream consumption depend on the
    /// values drawn, and that is a needless dependency in something whose whole
    /// job is being predictable.
    pub fn gaussian(&mut self) -> f64 {
        // Guard the log against an exact zero, which `unit` can return.
        let u1 = self.unit().max(f64::MIN_POSITIVE);
        let u2 = self.unit();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// A normal deviate with the given mean and standard deviation.
    pub fn normal(&mut self, mean: f64, std_dev: f64) -> f64 {
        mean + std_dev * self.gaussian()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The property the whole determinism claim rests on: an indexed stream does
    /// not care when it was asked for. Drawing the indices backwards gives the
    /// same numbers as drawing them forwards, so a parallel trace that finishes
    /// its work in any order still produces one answer.
    #[test]
    fn indexed_streams_are_order_free() {
        let forward: Vec<f64> = (0..64)
            .map(|i| Rng::for_index(0xD0A1_15EE, i).unit())
            .collect();
        let backward: Vec<f64> = (0..64)
            .rev()
            .map(|i| Rng::for_index(0xD0A1_15EE, i).unit())
            .collect();
        let mut backward_reordered = backward;
        backward_reordered.reverse();
        assert_eq!(forward, backward_reordered);
    }

    /// Adjacent indices must not be correlated: a per-pixel seed of `y * w + x`
    /// is the common case, and if index i and i+1 shared most of their stream the
    /// image would show it as structure.
    #[test]
    fn adjacent_indices_are_decorrelated() {
        let draws: Vec<Vec<f64>> = (0..16)
            .map(|i| {
                let mut r = Rng::for_index(7, i);
                (0..8).map(|_| r.unit()).collect()
            })
            .collect();
        for i in 0..draws.len() {
            for j in (i + 1)..draws.len() {
                assert_ne!(draws[i], draws[j], "streams {i} and {j} collided");
            }
        }
        // Neighbouring streams should not even agree on their first value to a
        // few digits, which a weak mixer would allow.
        for i in 1..draws.len() {
            assert!(
                (draws[i][0] - draws[i - 1][0]).abs() > 1e-6,
                "streams {} and {i} start too close",
                i - 1
            );
        }
    }

    /// Swapping seed and index gives a different stream — otherwise seeding by
    /// (frame, pixel) and (pixel, frame) would alias.
    #[test]
    fn seed_and_index_are_not_interchangeable() {
        assert_ne!(
            Rng::for_index(3, 9).unit(),
            Rng::for_index(9, 3).unit(),
            "the pair must be ordered"
        );
    }

    /// Zero is the one seed xorshift cannot use — it would emit zero forever.
    #[test]
    fn the_degenerate_seed_is_handled() {
        let mut zero = Rng::new(0);
        let draws: Vec<f64> = (0..4).map(|_| zero.unit()).collect();
        assert!(draws.iter().all(|&v| v > 0.0 && v < 1.0), "{draws:?}");
        assert!(draws[0] != draws[1]);
    }

    /// A fork diverges from its parent rather than replaying it.
    #[test]
    fn a_split_stream_diverges_from_its_parent() {
        let mut parent = Rng::new(42);
        let mut child = parent.split();
        let p: Vec<f64> = (0..8).map(|_| parent.unit()).collect();
        let c: Vec<f64> = (0..8).map(|_| child.unit()).collect();
        assert_ne!(p, c);
        // Cloning captures the position exactly, which is what makes a draw
        // replayable.
        let mut a = Rng::new(42);
        let mut b = a.clone();
        assert_eq!(a.unit(), b.unit());
    }

    /// The generator's numbers are fixed, not merely reproducible-in-principle.
    /// This is the test that would fail if the algorithm, the seeding or the
    /// [0,1) conversion were ever changed — which is the point, since a changed
    /// stream silently invalidates every recorded result.
    #[test]
    fn the_stream_is_pinned() {
        let mut r = Rng::new(0x5A17_7E3D);
        // A cheap order-sensitive digest of 10 000 draws.
        let mut hash = 0u64;
        for _ in 0..10_000 {
            hash = hash.rotate_left(7).wrapping_mul(0x1000_0000_01B3)
                ^ (r.unit() * (1u64 << 53) as f64) as u64;
        }
        assert_eq!(hash, PINNED_DIGEST, "the generator's output has changed");
    }

    /// Changing this constant is never the fix. If this test fails, the stream
    /// moved, and every result recorded against the old one is now unreproducible.
    const PINNED_DIGEST: u64 = 6_777_642_030_472_145_829;

    /// Box-Muller has to actually be normal: zero mean, unit variance, and a
    /// tail that reaches past three sigma without running away.
    #[test]
    fn gaussians_are_standard_normal() {
        let mut r = Rng::new(1234);
        const N: usize = 100_000;
        let draws: Vec<f64> = (0..N).map(|_| r.gaussian()).collect();
        let mean = draws.iter().sum::<f64>() / N as f64;
        let variance = draws.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / N as f64;
        assert!(mean.abs() < 0.02, "mean {mean}");
        assert!((variance - 1.0).abs() < 0.03, "variance {variance}");
        // About 0.27% of a normal distribution lies beyond three sigma.
        let tail = draws.iter().filter(|d| d.abs() > 3.0).count() as f64 / N as f64;
        assert!((tail - 0.0027).abs() < 0.002, "three-sigma tail {tail}");
        assert!(draws.iter().all(|d| d.is_finite()));
    }

    /// A cosine-weighted sample stays in the hemisphere and leans towards the
    /// normal. The mean of the cosine over a Lambertian distribution is 2/3, which
    /// is a closed form and therefore worth testing against.
    #[test]
    fn cosine_sampling_leans_towards_the_normal() {
        let mut r = Rng::new(99);
        let n = DVec3::new(1.0, 2.0, -0.5).normalize();
        const N: usize = 50_000;
        let mut cos_sum = 0.0;
        for _ in 0..N {
            let d = r.cosine_hemisphere(n);
            let c = d.dot(n);
            assert!(c > -1e-9, "sample left the hemisphere: {c}");
            assert!((d.length() - 1.0).abs() < 1e-9);
            cos_sum += c;
        }
        let mean_cos = cos_sum / N as f64;
        assert!(
            (mean_cos - 2.0 / 3.0).abs() < 0.01,
            "Lambertian mean cosine should be 2/3, got {mean_cos}"
        );
    }

    /// Uniform hemisphere sampling has mean cosine 1/2, which is how it differs
    /// from the cosine-weighted one — and getting the two confused is a
    /// factor-of-4/3 error in every diffuse bounce.
    #[test]
    fn uniform_hemisphere_is_not_cosine_weighted() {
        let mut r = Rng::new(5);
        let n = DVec3::Z;
        const N: usize = 50_000;
        let mean_cos: f64 = (0..N).map(|_| r.on_hemisphere(n).dot(n)).sum::<f64>() / N as f64;
        assert!((mean_cos - 0.5).abs() < 0.01, "got {mean_cos}");
    }
}
