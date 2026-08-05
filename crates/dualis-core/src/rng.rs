//! A deterministic pseudo-random generator, and the sampling built on it.
//!
//! Every stochastic choice in a simulation has to be reproducible or the results
//! cannot be compared, so this is a plain seeded generator with no global state
//! and no entropy source. Two runs of the same scene draw the same numbers, on
//! every platform and in WebAssembly.

use glam::DVec3;
/// Deterministic xorshift64* PRNG. Avoids a dependency and keeps scene
/// generation identical on every platform (important for WASM + tests).
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed | 1)
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
}
