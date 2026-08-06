//! dualis-optics: light, as a domain built on the `dualis-core` kernel.
//!
//! Spectral radiometry, surface optics, dispersion and ray geometry. This is the
//! first domain, and its job is partly to prove the kernel can hold one: it uses
//! the kernel's [`Violation`](dualis_core::Violation) for a surface that would
//! create light, the kernel's [`Rng`](dualis_core::Rng) for every stochastic
//! choice, and the kernel's units throughout — and the kernel knows nothing about
//! any of it.
//!
//! Units are SI in storage and dimensioned in the API. Wavelengths are
//! [`Length`](dualis_units::Length), so a path length cannot be passed as one; the
//! `_nm` fields inside [`Spectrum`] are the serialised form, where a readable
//! `"center_nm": 488` beats a correct `4.88e-7`.
//!
//! # The two invariants, in this domain
//!
//! **Energy is conserved at every surface.** [`SurfaceOptics`] stores reflectance
//! and transmittance; absorptance is the remainder. [`SurfaceOptics::validate`]
//! rejects a surface that would gain energy, and [`SurfaceOptics::split`]
//! renormalises rather than amplifying if one is traced anyway.
//!
//! **Nothing is random.** Every sampled scatter comes from a seeded generator, and
//! [`Rng::for_index`](dualis_core::Rng::for_index) makes that hold under
//! parallelism as well as in sequence.
//!
//! # What is in it
//!
//! | Module | |
//! | --- | --- |
//! | [`spectrum`] | Wavelength-dependent quantities: Planck's law, Gaussians, measured curves, filter bands with real edges and finite blocking |
//! | [`optics`] | What a surface does to light — Fresnel from the refractive indices, spectral R/T/A, coatings, polarisation split, scatter sampling |
//! | [`material`] | Refractive index against wavelength (Sellmeier and Cauchy, with a small glass catalogue), Abbe number, and how much light survives the glass |
//! | [`geometry`] | Rays and hits, intersections against caps, conics, planes, annuli and cylinders, Snell's law, and disc sampling |
//!
//! | [`diffraction`] | What a *perfect* system does: Airy patterns, encircled energy, the ideal MTF, Rayleigh and Abbe limits |
//! | [`wavefront`] | What an imperfect one does: Zernike aberrations, pupil transforms, aberrated PSFs and their MTFs |
//! | [`radiometry`] | Integrating a spectrum, and the difference between watts and photons |
//!
//! # What is deliberately not in it
//!
//! No scene graph, no meshes, no acceleration structure, no renderer.
//!
//! Wave optics is single-plane: a pupil transforms to an image and that is all. There
//! is no propagation between arbitrary planes, so a beam cannot be walked down a
//! bench, and there is no partial coherence — every PSF here is the incoherent one.

pub mod diffraction;
pub mod geometry;
pub mod material;
pub mod optics;
pub mod propagation;
pub mod radiometry;
pub mod spectrum;
pub mod wavefront;

pub use diffraction::{
    abbe_limit, airy_intensity, airy_radius, cutoff_frequency, depth_of_focus, encircled_energy,
    mtf_at, mtf_ideal, rayleigh_limit, strehl_from_wavefront_error,
};
pub use geometry::{
    annulus_intersect, cap_intersect, conic_intersect, conic_sag, cylinder_intersect,
    plane_intersect, refract, sag, Hit, Ray,
};
pub use material::{Dispersion, Material, C_LINE, D_LINE, F_LINE};
pub use optics::{
    brewster_angle, critical_angle, fresnel_reflectance, fresnel_split, Scatter, SurfaceFinish,
    SurfaceOptics,
};
pub use propagation::{
    gaussian_divergence, gaussian_radius_at, rayleigh_range, Grid, PropagationError,
};
pub use radiometry::SpectralPower;
pub use spectrum::{Spectrum, VISIBLE_RANGE};
pub use wavefront::{Mtf, Psf, Pupil, Zernike};

/// The kernel, re-exported so a consumer of this crate does not need to name it
/// separately to get at `Rng`, `Simulation` or the units.
pub use dualis_core as core;
pub use dualis_units as units;
