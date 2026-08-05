//! dualis-core: the physics underneath a simulated world.
//!
//! This is the layer that does not know what it is being used to simulate. It
//! knows that light has a wavelength, that a surface divides incident power
//! three ways, that glass bends light by an amount its refractive index decides,
//! that things move, and that a random choice has to be reproducible. What you
//! build on top of that — a microscope, a camera, a room — is somebody else's
//! problem.
//!
//! Units are **millimetres** for length, **nanometres** for wavelength and
//! **seconds** for time, everywhere and without exception.
//!
//! # The invariants
//!
//! Two rules hold throughout, and the tests exist to keep them holding:
//!
//! - **Energy is conserved at every surface.** Reflectance and transmittance are
//!   stored; absorptance is whatever is left. There is no way to write down a
//!   surface that returns more light than reached it.
//! - **Nothing is random.** Every stochastic choice comes from a seeded
//!   generator with no global state, so two runs of the same scene agree to the
//!   last bit, on every platform and in WebAssembly.
//!
//! # What is here
//!
//! | Module | |
//! | --- | --- |
//! | [`spectrum`] | Wavelength-dependent quantities: Planck's law, Gaussians, measured curves, filter bands |
//! | [`optics`] | What a surface does to light — Fresnel, reflectance, transmittance, absorptance |
//! | [`material`] | Refractive index against wavelength, and how much survives the glass |
//! | [`motion`] | Rigid motion and time gating: drift, oscillation, spin, strobe |
//! | [`geometry`] | Ray intersections, Snell's law, and the disc samplings that make ray bundles |
//! | [`rng`] | A deterministic generator, and sampling built on it |

pub mod geometry;
pub mod material;
pub mod motion;
pub mod optics;
pub mod rng;
pub mod spectrum;

pub use material::{Dispersion, Material};
pub use motion::{Motion, Strobe};
pub use optics::{fresnel_reflectance, SurfaceFinish, SurfaceOptics};
pub use rng::Rng;
pub use spectrum::Spectrum;
