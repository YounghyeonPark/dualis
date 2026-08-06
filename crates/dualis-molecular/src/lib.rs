//! dualis-molecular: matter one atom at a time, as a domain on the `dualis-core` kernel.
//!
//! The fifth domain, and the first whose answers are *distributions* rather than values.
//! Optics gives a Fresnel coefficient, heat gives a temperature, sound gives a mode frequency;
//! a hundred atoms in a box give a trajectory that is chaotic by construction — perturb one
//! coordinate in its last bit and the paths separate completely inside a few hundred steps.
//!
//! That is the physics rather than a numerical failing, and it changes what a test can be. The
//! trajectory is not reproducible in any useful sense; what it *averages to* is exact, and
//! there is a great deal of it:
//!
//! | Claim | Exact form |
//! | --- | --- |
//! | Equipartition | `⟨KE⟩ = (3N − 3) k_BT / 2` |
//! | Ideal gas, at low density | `PV = N k_BT` |
//! | Virial pressure, at any density | `PV = N k_BT + ⟨Σ f·r⟩/3` |
//! | Lennard-Jones minimum | `−ε` exactly, at `2^(1/6) σ` |
//! | Momentum | Conserved to the last bit, by Newton's third law |
//! | Energy, unthermostatted | Bounded rather than drifting, because Verlet is symplectic |
//!
//! # What it added to the kernel
//!
//! Nothing, which is the fifth time that has held. It does lean on two kernel pieces harder
//! than any previous domain: [`Rng::for_index`](dualis_core::Rng::for_index), because a
//! Langevin thermostat is *made* of random numbers and a chaotic system offers no other route
//! back to a reproducible run; and the symplectic argument behind
//! [`velocity_verlet`](dualis_core::velocity_verlet), which the kernel proves on a harmonic
//! oscillator and which is relied on here for a many-body potential.
//!
//! # What is deliberately not here
//!
//! No bonds, angles or torsions, so no molecules — this is a monatomic fluid, and the name is
//! aspirational by one step. No electrostatics, which for a charged system is the hard part:
//! Coulomb falls off as `1/r` and cannot be cut off at all, so it needs Ewald summation or a
//! particle-mesh method, and that is a larger piece of work than everything here put together.
//! No constraints, no barostat, no free-energy machinery.
//!
//! What is here is the part with closed forms to check against, which is the same line every
//! other domain in this workspace is drawn along.

pub mod box_;
pub mod fluid;
pub mod potential;

pub use box_::{CellList, PeriodicBox};
pub use fluid::{
    reduced_density, reduced_temperature, sigma_of, temperature_from_reduced, unit_mass, Fluid,
    Thermostat,
};
pub use potential::{LennardJones, Pair};
