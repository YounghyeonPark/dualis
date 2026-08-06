//! dualis-core: the kernel a simulated world's physics is built on.
//!
//! This crate knows nothing about any particular physics. It knows that a
//! quantity can vary over space and time, that a process must answer for what it
//! conserves, that a system with no closed form has to be rolled forward, that
//! matter has properties several domains need at once, and that several domains
//! sharing a clock is a scheduling problem with real failure modes. What any of
//! that is *about* — light, heat, contact, sound — belongs to a domain crate.
//!
//! That separation is the point. `dualis-optics` depends on this crate; this crate
//! must never depend on it, or anything else that models a specific physics. If a
//! new domain needs the kernel changed, the kernel was wrong.
//!
//! # The two invariants
//!
//! Both survive the generalisation, and both are now enforced rather than
//! promised:
//!
//! - **Nothing is created or destroyed without being noticed.** A [`Ledger`] is
//!   what a process claims to hold and [`audit`] is the check; energy crossing
//!   between domains goes through [`Exchange`], which refuses to let a transfer
//!   silently lose some. This generalises what `SurfaceOptics` did for one
//!   quantity at one kind of boundary. Where a boundary is resolved into faces,
//!   the audit is per face — a redistribution that keeps the total but moves it
//!   to the wrong part of a mirror is the one bug a total-only check cannot see.
//! - **Nothing is random.** [`Rng::for_index`] gives every piece of work its own
//!   stateless stream, so a parallel simulation is still bit-reproducible — which
//!   is when the guarantee starts to matter, and when a single shared generator
//!   would have quietly lost it.
//!
//! # What is here
//!
//! | Module | |
//! | --- | --- |
//! | [`conserved`] | Conservation as an audit: ledgers, violations, tolerances |
//! | [`integrator`] | Fixed-step time evolution, and why symplectic beats accurate |
//! | [`sim`] | Several domains on one clock: quasi-static, multirate, iterative coupling |
//! | [`scene`] | Where two domains meet: shared boundaries, and flux that knows its place |
//! | [`field`] | Scalar and vector fields, with gradient, divergence, curl, Laplacian |
//! | [`substance`] | Thermal, mechanical and acoustic properties of matter |
//! | [`motion`] | Closed-form rigid motion and time gating |
//! | [`rng`] | A deterministic generator, and the sampling built on it |
//! | [`transform`] | The discrete Fourier transform, accurate rather than fast |
//! | [`vector`] | Basis construction and reflection — the vector maths no domain owns |
//!
//! Units come from `dualis-units` and are re-exported below, so a domain crate
//! needs one dependency rather than two.

pub mod conserved;
pub mod field;
pub mod integrator;
pub mod motion;
pub mod rng;
pub mod scene;
pub mod sim;
pub mod substance;
pub mod transform;
pub mod vector;

pub use conserved::{audit, Conserves, Ledger, Violation};
pub use field::{ScalarField, VectorField};
pub use integrator::{velocity_verlet, Dynamics, Integrator, Newtonian, State};
pub use motion::{Motion, Strobe};
pub use rng::Rng;
pub use scene::{Flux, Interface};
pub use sim::{Domain, Exchange, Kind, Schedule, Simulation};
pub use substance::Substance;
pub use transform::{fft, fft2, fftshift, ifft, ifft2};
pub use vector::{basis_for, oriented_against, reflect};

/// Everything from `dualis-units`, so that `use dualis_core::units::*` is enough
/// to write dimensioned physics.
pub mod units {
    pub use dualis_units::*;
}
