//! Many independent samples, run in parallel, with an answer that does not depend on how many
//! threads did the work.
//!
//! The other axis of parallelism in this workspace. `TreeNBody::with_threads` splits *one*
//! evaluation across cores; this splits *many* evaluations, which is the shape a Monte Carlo
//! study has and the shape a parameter sweep has. A hundred thousand detector realisations, a
//! thousand trajectories from perturbed initial conditions, a grid of designs each run to a
//! steady state — all the same problem: independent work items and a reduction at the end.
//!
//! # Why this can be parallel *and* still bit-for-bit
//!
//! Two decisions that were already made, meeting.
//!
//! [`Rng::for_index`] is stateless and addressed by index, so sample `i` draws the same numbers
//! whether it ran first, last, or on another core. Nothing is consumed from a shared stream and
//! there is no order to depend on.
//!
//! And results land in a slot chosen by index, never appended. Each worker owns a disjoint run
//! of the output and writes nothing else, exactly as `TreeNBody` does — so the vector that comes
//! back is a function of `(seed, count)` alone. Reduce it however you like; a fold over that
//! vector is in index order by construction.
//!
//! The failure this avoids is worth naming, because it is the usual one: a Monte Carlo that
//! draws from a shared generator gives a different answer on eight cores than on one, and the
//! difference looks like statistical noise. It is not noise, it is the result depending on the
//! scheduler, and no amount of averaging removes it.
//!
//! ```
//! use pantometry_core::{Ensemble, Rng};
//!
//! // A hundred thousand throws of a loaded die, in parallel.
//! let hits = Ensemble::new(20, 100_000)
//!     .with_threads(8)
//!     .run(|_, mut rng| u64::from(rng.unit() < 0.25));
//!
//! let heads: u64 = hits.iter().sum();
//! // Same seed, same count, same answer — on one thread or on eight.
//! assert_eq!(heads, Ensemble::new(20, 100_000).run(|_, mut rng| u64::from(rng.unit() < 0.25))
//!     .iter().sum::<u64>());
//! ```

use crate::Rng;

/// A set of independent samples to run.
///
/// Cheap to build and to copy; it holds a seed, a count and a thread count, and does the work in
/// [`Ensemble::run`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Ensemble {
    seed: u64,
    count: u64,
    threads: usize,
}

impl Ensemble {
    /// `count` samples, each drawing from `Rng::for_index(seed, i)`.
    ///
    /// Sequential until [`with_threads`](Ensemble::with_threads) says otherwise, which is the
    /// right default: threads are a performance decision and this crate does not make those for
    /// a caller who has not asked.
    pub fn new(seed: u64, count: u64) -> Ensemble {
        Ensemble {
            seed,
            count,
            threads: 1,
        }
    }

    /// How many threads to spread the samples over. 1 is sequential.
    ///
    /// **The answer does not change.** That is the whole point, and it is asserted rather than
    /// asserted-in-prose: a test runs the same ensemble at one, three and sixteen threads and
    /// compares the results bit for bit. If you find a thread count that changes an answer, the
    /// sample closure is reading something it does not own.
    ///
    /// Clamped to at least one, and never more than there are samples.
    pub fn with_threads(mut self, threads: usize) -> Ensemble {
        self.threads = threads.max(1);
        self
    }

    /// The seed every sample's generator is derived from.
    pub fn seed(&self) -> u64 {
        self.seed
    }

    /// How many samples.
    pub fn count(&self) -> u64 {
        self.count
    }

    /// Run every sample and collect the results **in index order**.
    ///
    /// The closure receives the sample's index and its own generator. Give it everything else it
    /// needs by capture; it must not mutate shared state, and `Fn` rather than `FnMut` is how
    /// that is enforced rather than requested.
    ///
    /// Index order matters more than it looks. A caller folding the result — a mean, a variance,
    /// a histogram — folds in that order whatever the thread count was, so the floating-point
    /// sum is the same sum. Collecting into a shared accumulator instead would make the answer
    /// depend on which thread finished first, in the last bits, invisibly.
    pub fn run<T, F>(&self, sample: F) -> Vec<T>
    where
        F: Fn(u64, Rng) -> T + Sync,
        T: Send + Default + Clone,
    {
        let n = self.count as usize;
        let mut out = vec![T::default(); n];
        let threads = self.threads.min(n.max(1));
        let seed = self.seed;

        if threads <= 1 {
            for (i, slot) in out.iter_mut().enumerate() {
                *slot = sample(i as u64, Rng::for_index(seed, i as u64));
            }
            return out;
        }

        #[cfg(not(target_family = "wasm"))]
        {
            let chunk = n.div_ceil(threads);
            let sample = &sample;
            // Disjoint slices, no reduction, nothing shared but the closure.
            std::thread::scope(|scope| {
                for (c, slice) in out.chunks_mut(chunk).enumerate() {
                    let base = (c * chunk) as u64;
                    scope.spawn(move || {
                        for (k, slot) in slice.iter_mut().enumerate() {
                            let i = base + k as u64;
                            *slot = sample(i, Rng::for_index(seed, i));
                        }
                    });
                }
            });
            out
        }

        // WebAssembly has no threads to spawn, so it takes the sequential path and gets the same
        // answer for a less interesting reason. `TreeNBody` resolves this the same way.
        #[cfg(target_family = "wasm")]
        {
            for (i, slot) in out.iter_mut().enumerate() {
                *slot = sample(i as u64, Rng::for_index(seed, i as u64));
            }
            out
        }
    }

    /// Run every sample and reduce to a mean and a standard error, folded in index order.
    ///
    /// The two numbers a Monte Carlo study is usually for: the estimate, and how much to trust
    /// it. The standard error is `s/√N` with the sample standard deviation, so it falls as
    /// `1/√N` — which is the rate this workspace asks a tolerance to be earned against, and the
    /// reason `samples` is reported beside it rather than left implicit.
    ///
    /// Returns `None` for fewer than two samples, because a variance over one is not a small
    /// number, it is not defined.
    pub fn estimate<F>(&self, sample: F) -> Option<Estimate>
    where
        F: Fn(u64, Rng) -> f64 + Sync,
    {
        if self.count < 2 {
            return None;
        }
        // Per-block partials rather than every sample, so a study is bounded by its *block*
        // count and not by its sample count. A hundred million f64 is 800 MB held for no reason;
        // this holds one `Partial` per BLOCK samples, which is 24 bytes per 4096.
        let blocks = self.blocks(
            |from, to, sample| Partial::of(from, to, self.seed, sample),
            &sample,
        );
        let total: Partial = blocks
            .iter()
            .copied()
            .reduce(Partial::merge)
            .expect("count >= 2 means at least one block");
        Some(total.finish())
    }

    /// Run the samples in fixed-size blocks and return one value per block.
    ///
    /// **The block size does not depend on the thread count**, and that is the whole reason for
    /// it. A reduction split per *thread* combines a different number of partial sums on four
    /// cores than on sixteen, and floating-point addition is not associative, so the answer
    /// moves — quietly, in the last bits, looking like nothing. Splitting per fixed block makes
    /// the association a function of `count` alone.
    ///
    /// A caller wanting a reduction this crate does not provide — a histogram, a maximum, a
    /// quantile — should build it here for the same reason.
    pub fn blocks<B, M, F>(&self, of_block: M, sample: &F) -> Vec<B>
    where
        M: Fn(u64, u64, &F) -> B + Sync,
        B: Send + Default + Clone,
        F: Sync,
    {
        let n_blocks = self.count.div_ceil(BLOCK) as usize;
        let mut out = vec![B::default(); n_blocks];
        let threads = self.threads.min(n_blocks.max(1));
        let (count, of_block) = (self.count, &of_block);

        let one = |slice: &mut [B], base: usize| {
            for (k, slot) in slice.iter_mut().enumerate() {
                let from = (base + k) as u64 * BLOCK;
                *slot = of_block(from, (from + BLOCK).min(count), sample);
            }
        };

        if threads <= 1 {
            one(&mut out, 0);
            return out;
        }

        #[cfg(not(target_family = "wasm"))]
        {
            let chunk = n_blocks.div_ceil(threads);
            std::thread::scope(|scope| {
                for (c, slice) in out.chunks_mut(chunk).enumerate() {
                    let base = c * chunk;
                    scope.spawn(move || one(slice, base));
                }
            });
            out
        }

        #[cfg(target_family = "wasm")]
        {
            one(&mut out, 0);
            out
        }
    }
}

/// Samples per reduction block. A power of two, and fixed: see [`Ensemble::blocks`].
///
/// 4096 doubles is 32 KB of intermediate per block, which stays in L1 while a block is folded,
/// and it keeps the number of partials small enough that combining them costs nothing.
const BLOCK: u64 = 4096;

/// One block's contribution to a mean and a variance.
///
/// Carries a count, a mean and the sum of squared deviations rather than raw power sums.
/// Merging two of these is Chan's parallel update, which is stable where `sum(x²) − n·mean²`
/// is not: that form subtracts two large nearly-equal numbers and loses every significant digit
/// exactly when a Monte Carlo has converged and the mean dwarfs the spread.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Partial {
    n: f64,
    mean: f64,
    m2: f64,
}

impl Partial {
    fn of<F: Fn(u64, Rng) -> f64>(from: u64, to: u64, seed: u64, sample: &F) -> Partial {
        let mut p = Partial::default();
        for i in from..to {
            let x = sample(i, Rng::for_index(seed, i));
            // Welford, in index order within the block.
            p.n += 1.0;
            let delta = x - p.mean;
            p.mean += delta / p.n;
            p.m2 += delta * (x - p.mean);
        }
        p
    }

    fn merge(a: Partial, b: Partial) -> Partial {
        if a.n == 0.0 {
            return b;
        }
        if b.n == 0.0 {
            return a;
        }
        let n = a.n + b.n;
        let delta = b.mean - a.mean;
        Partial {
            n,
            mean: a.mean + delta * (b.n / n),
            m2: a.m2 + b.m2 + delta * delta * (a.n * b.n / n),
        }
    }

    fn finish(self) -> Estimate {
        let variance = self.m2 / (self.n - 1.0);
        Estimate {
            mean: self.mean,
            standard_error: (variance / self.n).sqrt(),
            samples: self.n as u64,
        }
    }
}

/// What a Monte Carlo run came back with.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Estimate {
    /// The sample mean.
    pub mean: f64,
    /// `s/√N`: how far the mean is likely to be from the truth, not how spread the samples are.
    pub standard_error: f64,
    /// How many samples went into it. Quoted because a mean without one is not a measurement.
    pub samples: u64,
}

impl Estimate {
    /// The sample standard deviation — the spread of the samples themselves.
    ///
    /// Distinct from [`standard_error`](Estimate::standard_error), and confusing the two is the
    /// most common way to state a Monte Carlo result wrongly: the spread does not shrink with
    /// more samples and the error on the mean does.
    pub fn standard_deviation(&self) -> f64 {
        self.standard_error * (self.samples as f64).sqrt()
    }

    /// Whether a value sits within `k` standard errors of the mean.
    pub fn within(&self, k: f64, value: f64) -> bool {
        (value - self.mean).abs() <= k * self.standard_error
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The thread count does not change a single bit.**
    ///
    /// The claim the whole type exists to make. A Monte Carlo drawing from a shared generator
    /// gives a different answer on eight cores than on one, and the difference *looks like*
    /// statistical noise — so it is never investigated, and no amount of averaging removes it.
    ///
    /// Compared on `to_bits()` element by element rather than on a mean, because a mean can
    /// agree while the samples behind it are permuted.
    #[test]
    fn the_answer_does_not_depend_on_how_many_threads_produced_it() {
        let draw = |i: u64, mut rng: Rng| rng.gaussian() * (1.0 + i as f64 % 3.0);
        let one = Ensemble::new(4242, 5_000).run(draw);

        for threads in [2usize, 3, 7, 16, 64] {
            let many = Ensemble::new(4242, 5_000).with_threads(threads).run(draw);
            assert_eq!(one.len(), many.len());
            for (i, (a, b)) in one.iter().zip(&many).enumerate() {
                assert_eq!(
                    a.to_bits(),
                    b.to_bits(),
                    "sample {i} differed at {threads} threads: {a} against {b}"
                );
            }
        }

        // And the samples are not all the same value, which would make the above vacuous.
        let spread = one.iter().cloned().fold(f64::MIN, f64::max)
            - one.iter().cloned().fold(f64::MAX, f64::min);
        assert!(spread > 1.0, "the draws should vary; spread {spread}");
    }

    /// **A mean converges at `1/√N`, which is the rate that earns a tolerance.**
    ///
    /// The estimator is checked against a distribution whose mean is known exactly rather than
    /// against another run: a uniform draw on `[0, 1)` has mean `1/2`. Quadrupling the samples
    /// must halve the error, and that *rate* is the claim — a single count would only say the
    /// estimator is not wildly wrong.
    #[test]
    fn the_error_falls_as_one_over_root_n() {
        let err = |n: u64| {
            let e = Ensemble::new(7, n)
                .with_threads(4)
                .estimate(|_, mut rng| rng.unit())
                .expect("more than one sample");
            (e.mean - 0.5).abs()
        };
        let (coarse, fine) = (err(4_000), err(64_000));
        // Sixteen times the samples is four times the accuracy. Bounded loosely on purpose:
        // this is a random quantity and the band is wide enough that an honest run passes and
        // a broken estimator — one whose error does not fall at all — does not.
        let ratio = coarse / fine;
        assert!(
            (1.5..12.0).contains(&ratio),
            "16x the samples gave {ratio:.2}x the accuracy (expected about 4)"
        );

        // The reported standard error should bracket the truth most of the time. One estimate,
        // three sigma: this fails about three times in a thousand by chance, and the seed is
        // fixed, so it either passes forever or reports a real defect.
        let e = Ensemble::new(11, 40_000)
            .with_threads(4)
            .estimate(|_, mut rng| rng.unit())
            .unwrap();
        assert!(
            e.within(3.0, 0.5),
            "0.5 is {:.2} standard errors from {:.6}",
            (e.mean - 0.5).abs() / e.standard_error,
            e.mean
        );
        assert_eq!(e.samples, 40_000);

        // Spread and error-on-the-mean are different numbers, and confusing them is the usual
        // way to misreport a Monte Carlo. A uniform draw has standard deviation 1/√12.
        assert!(
            (e.standard_deviation() - (1.0f64 / 12.0).sqrt()).abs() < 0.01,
            "standard deviation {:.6}",
            e.standard_deviation()
        );
    }

    /// **Ten million samples, held in kilobytes rather than in eighty megabytes.**
    ///
    /// The reason `estimate` folds per block instead of collecting. `run` materialises every
    /// sample and is right when you want them; a study that only wants a mean should not pay
    /// 8 bytes times the sample count to get one, and at the sizes a Monte Carlo actually
    /// reaches — 1e8, 1e9 — paying it is not merely wasteful but impossible.
    ///
    /// Still thread-independent, which is the harder half: the block size is fixed, so the
    /// *association* of the additions is a function of the sample count alone and not of how
    /// many cores turned up.
    #[test]
    fn a_run_too_large_to_hold_still_agrees_across_threads() {
        let draw = |_: u64, mut rng: Rng| rng.unit();
        let n = 10_000_000;

        let one = Ensemble::new(5, n).estimate(draw).expect("plenty");
        for threads in [4usize, 16] {
            let many = Ensemble::new(5, n)
                .with_threads(threads)
                .estimate(draw)
                .expect("plenty");
            assert_eq!(
                one.mean.to_bits(),
                many.mean.to_bits(),
                "mean moved at {threads} threads: {} against {}",
                one.mean,
                many.mean
            );
            assert_eq!(one.standard_error.to_bits(), many.standard_error.to_bits());
            assert_eq!(one.samples, n);
        }

        // And it is right: a uniform draw has mean 1/2 and standard deviation 1/sqrt(12), and
        // ten million samples pin the mean to about a ten-thousandth.
        assert!(one.within(4.0, 0.5), "mean {:.8}", one.mean);
        assert!(
            (one.standard_deviation() - (1.0f64 / 12.0).sqrt()).abs() < 1e-3,
            "spread {:.6}",
            one.standard_deviation()
        );
    }

    /// **The blocked fold is more accurate than a flat sum, not merely cheaper.**
    ///
    /// A mean far from zero against a tiny spread is where the naive `sum(x²) − n·mean²` form
    /// loses every digit, and where a flat left-to-right sum of ten million values loses several
    /// to accumulation. Welford within a block and Chan's merge between blocks keeps both.
    ///
    /// Checked against a case whose answer is exact: `x = 1e9 + (i mod 2)` has mean
    /// `1e9 + 0.5` and variance exactly `0.25 · n/(n−1)`.
    #[test]
    fn the_estimator_survives_a_large_mean_and_a_small_spread() {
        let n = 1_000_000u64;
        let e = Ensemble::new(0, n)
            .with_threads(8)
            .estimate(|i, _| 1e9 + (i % 2) as f64)
            .expect("plenty");

        assert!(
            (e.mean - (1e9 + 0.5)).abs() < 1e-6,
            "mean {:.6} against 1000000000.5",
            e.mean
        );
        // Population variance 0.25, so the sample variance is 0.25·n/(n−1).
        let want = (0.25 * n as f64 / (n as f64 - 1.0)).sqrt();
        assert!(
            (e.standard_deviation() / want - 1.0).abs() < 1e-9,
            "spread {:.9} against {want:.9}",
            e.standard_deviation()
        );
    }

    /// Fewer than two samples has no variance, and says so rather than returning zero.
    #[test]
    fn one_sample_is_not_an_estimate() {
        assert!(Ensemble::new(1, 1).estimate(|_, mut r| r.unit()).is_none());
        assert!(Ensemble::new(1, 0).estimate(|_, mut r| r.unit()).is_none());
        assert!(Ensemble::new(1, 2).estimate(|_, mut r| r.unit()).is_some());
    }
}
