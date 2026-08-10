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
//! | [`bodies`] | The other shape a domain can be: a countable number of things at places |
//! | [`pose`] | Rigid motion — a rotation and a translation, and deliberately nothing more |
//! | [`substance`] | Thermal, mechanical and acoustic properties of matter |
//! | [`motion`] | Closed-form rigid motion and time gating |
//! | [`rng`] | A deterministic generator, and the sampling built on it |
//! | [`ensemble`] | Many independent samples in parallel, with an answer that does not depend on how many threads produced it |
//! | [`transform`] | The discrete Fourier transform, accurate rather than fast |
//! | [`vector`] | Basis construction and reflection — the vector maths no domain owns |
//!
//! [`scene`] here is **not** the `dualis-scene` crate, and the collision is worth naming. This
//! module is where two domains *meet* — an [`Interface`] cut into faces and a [`Flux`] that
//! knows which one it crossed — and it is kernel business because the audit runs on it. The
//! crate one layer up is where a domain *sits*, which is a statement about the world and not
//! about physics, and it is above every domain for that reason.
//!
//! # What a domain offers a layer above it
//!
//! Four optional accessors, none of which the kernel uses itself: [`ScalarField`] through
//! `as_field` for a continuum, [`Bodies`] through `as_bodies` for a countable set, [`Reading`]
//! through `readings` for the scalars a domain has when it has no picture, and `as_any` for a
//! caller that knows the concrete type.
//!
//! They exist so that a layer which must visit every domain never has to name one. **All four
//! are opt-in and default to nothing**, which is the hazard worth stating once: a domain that
//! forgets is silently absent from every table and every picture rather than failing to
//! compile. That has happened — four mechanics domains never opted into `as_any`, and an orbit
//! scene ran, conserved, and drew nothing at all.
//!
//! # A domain, and the audit that checks it
//!
//! ```
//! use dualis_core::conserved::quantity;
//! use dualis_core::{Domain, Exchange, Kind, Ledger, Schedule, Simulation, Violation};
//! use dualis_core::units::Time;
//!
//! /// A source that pays out a watt, and says so in its books.
//! struct Lamp { paid: f64 }
//! impl Domain for Lamp {
//!     fn name(&self) -> &str { "lamp" }
//!     fn kind(&self) -> Kind { Kind::QuasiStatic }
//!     fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
//!         let joules = 1.0 * dt.to_si();
//!         bus.publish(quantity::ENERGY, joules);
//!         self.paid += joules;
//!         Ok(())
//!     }
//!     // Negative: it is holding a debt, having handed the energy away.
//!     fn ledger(&self) -> Ledger { Ledger::new().with(quantity::ENERGY, -self.paid) }
//! }
//!
//! /// A sink that takes whatever is offered and keeps it.
//! struct Block { held: f64 }
//! impl Domain for Block {
//!     fn name(&self) -> &str { "block" }
//!     fn step(&mut self, _t: Time, _dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
//!         self.held += bus.take(quantity::ENERGY);
//!         Ok(())
//!     }
//!     fn ledger(&self) -> Ledger { Ledger::new().with(quantity::ENERGY, self.held) }
//!     // Opt in to being readable from outside. Without this, `domain_as` returns `None`:
//!     // the coupling never needs the concrete type, so it is not given away by default.
//!     fn as_any(&self) -> Option<&dyn std::any::Any> { Some(self) }
//!     fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> { Some(self) }
//! }
//!
//! let mut sim = Simulation::new(Schedule::Staggered)
//!     .with(Lamp { paid: 0.0 })
//!     .with(Block { held: 0.0 });
//! for _ in 0..10 {
//!     sim.advance(Time::ms(100.0)).expect("the books balance");
//! }
//!
//! // A joule went across, and the two ledgers cancel because nothing was lost.
//! let block: &Block = sim.domain_as("block").unwrap();
//! assert!((block.held - 1.0).abs() < 1e-12);
//! assert!(sim.ledger().get(quantity::ENERGY).unwrap().abs() < 1e-12);
//! ```
//!
//! Had `Block` consumed only half of what was published, `advance` would have returned a
//! [`Violation`] naming the channel rather than quietly losing the rest.
//!
//! Units come from `dualis-units` and are re-exported below, so a domain crate
//! needs one dependency rather than two.

// Every public item carries a doc comment. Denied rather than warned: a public physics API
// whose `Length::mm` shows a blank summary in rustdoc is documented in the sense that a
// paragraph exists somewhere, and not in the sense a reader needs.
#![deny(missing_docs)]
pub mod bodies;
pub mod conserved;
pub mod ensemble;
pub mod field;
pub mod integrator;
pub mod motion;
pub mod pose;
pub mod rng;
pub mod scene;
pub mod sim;
pub mod substance;
pub mod transform;
pub mod vector;

pub use bodies::Bodies;
pub use conserved::{audit, Conserves, Ledger, Violation};
pub use ensemble::{Ensemble, Estimate};
pub use field::{ScalarField, VectorField};
pub use integrator::{velocity_verlet, Dynamics, Integrator, Newtonian, State};
pub use motion::{Motion, Strobe};
pub use pose::Pose;
pub use rng::Rng;
pub use scene::{Flux, Interface};
pub use sim::{Domain, Exchange, Kind, Reading, Report, Schedule, Simulation};
pub use substance::Substance;
pub use transform::{fft, fft2, fftshift, ifft, ifft2};
pub use vector::{basis_for, oriented_against, reflect};

/// Everything from `dualis-units`, so that `use dualis_core::units::*` is enough
/// to write dimensioned physics.
pub mod units {
    pub use dualis_units::*;
}
