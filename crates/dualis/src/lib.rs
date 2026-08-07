//! dualis: physics for simulated worlds, in one dependency.
//!
//! A facade over the workspace. Nothing is implemented here — the point is that a
//! consumer writes `dualis = "0.1"` rather than naming seven crates, and that the
//! integration tests which need two domains at once have somewhere to live.
//!
//! ```
//! use dualis::prelude::*;
//!
//! // A lamp, a filter, and the question that needs both.
//! let lamp = SpectralPower::new(
//!     Spectrum::blackbody(3200.0),
//!     Power::w(1.0),
//!     VISIBLE_RANGE,
//! );
//! let green = Spectrum::bands(vec![[500.0, 560.0]], 0.95, 0.0);
//! let through = lamp.through(&green);
//! assert!(through < lamp.total());
//! ```
//!
//! # The dependency rule
//!
//! ```text
//! dualis-units      no dependencies but glam and serde
//! dualis-core       depends on units
//! dualis-optics     depends on core
//! dualis-thermal    depends on core
//! dualis-mechanics  depends on core
//! dualis-acoustic   depends on core
//! dualis-molecular  depends on core
//! dualis            depends on all of them
//! ```
//!
//! None of the five domains knows about any of the others. They meet on the kernel's
//! [`Exchange`](dualis_core::Exchange), and each one that arrived left the others
//! untouched — which is the claim the split was made to test, now held five times.

// Every public item carries a doc comment. Denied rather than warned: a public physics API
// whose `Length::mm` shows a blank summary in rustdoc is documented in the sense that a
// paragraph exists somewhere, and not in the sense a reader needs.
#![deny(missing_docs)]
pub use dualis_acoustic as acoustic;
pub use dualis_core as core;
pub use dualis_mechanics as mechanics;
pub use dualis_molecular as molecular;
pub use dualis_optics as optics;
pub use dualis_thermal as thermal;
pub use dualis_units as units;

/// Everything most simulations need, in one `use`.
pub mod prelude {
    pub use dualis_acoustic::{impedance, reflection_coefficient, End, Impedance, Tube};
    pub use dualis_core::conserved::quantity;
    pub use dualis_core::{
        audit, basis_for, oriented_against, reflect, velocity_verlet, Domain, Dynamics, Exchange,
        Flux, Integrator, Interface, Kind, Ledger, Motion, Newtonian, Rng, ScalarField, Schedule,
        Simulation, State, Strobe, Substance, VectorField, Violation,
    };
    pub use dualis_mechanics::{
        Body, ContactSystem, Coords, Ground, Inertia, NBody, RigidBody, TreeNBody, GRAVITATION,
    };
    pub use dualis_molecular::{Fluid, LennardJones, PeriodicBox, Thermostat};
    pub use dualis_optics::diffraction::{
        abbe_limit, airy_radius, cutoff_frequency, depth_of_focus, encircled_energy, mtf_at,
        rayleigh_limit, strehl_from_wavefront_error,
    };
    pub use dualis_optics::wavefront::AIRY_ZERO_LAMBDA_OVER_D;
    pub use dualis_optics::{
        fresnel_reflectance, fresnel_split, Hit, Material, Mtf, Psf, Pupil, Ray, Scatter,
        SpectralPower, Spectrum, SurfaceFinish, SurfaceOptics, Zernike, VISIBLE_RANGE,
    };
    pub use dualis_thermal::{Bar1D, Environment, LumpedMass, HEAT};
    pub use dualis_units::{
        AccelerationVec, Area, Damping, Density, Energy, Force, ForceVec, Frequency, HeatCapacity,
        Irradiance, Length, LengthVec, Mass, Momentum, MomentumVec, Power, Pressure, SpecificHeat,
        Stiffness, Temperature, Time, Velocity, VelocityVec, Volume, G0,
    };
}
