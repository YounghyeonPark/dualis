//! dualis: physics for simulated worlds, in one dependency.
//!
//! A facade over the workspace. Nothing is implemented here — the point is that a
//! consumer writes `dualis = "0.10"` rather than naming eleven crates, and that the
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
//! # Start here
//!
//! Three ideas carry the whole library.
//!
//! 1. **Units are types.** Dimensions live in the type, so `Length + Time` does not compile.
//!    One place may hold a factor of a thousand — a unit-bearing constructor — and `to_si()`
//!    is the only way back to a bare `f64`.
//! 2. **A domain is anything that steps.** [`Domain`](dualis_core::Domain) requires
//!    `name` and `step`; everything else has a default. Override `ledger` so the audit has
//!    something to check.
//! 3. **Domains never call each other.** They meet on
//!    [`Exchange`](dualis_core::Exchange), a bus of named channels carrying SI *amounts* —
//!    joules, not watts. A ledger says what you are holding, not what has passed through you.
//!
//! And the reason to pick this over a general-purpose engine: conservation is audited every
//! step, so a wrong model does not run quietly. `advance` returns a
//! [`Violation`](dualis_core::Violation) naming the quantity, the site, and the before and
//! after — a correctness signal you can act on without a human noticing first.
//!
//! ```text
//! energy destroyed at simulation: 5.000000e2 became 4.995000e2,
//! a relative change of 1.000e-3 against a tolerance of 1.000e-9
//! ```
//!
//! Be clear about the limit: the audit catches quantities appearing or vanishing, amounts
//! left unclaimed on the bus, and fluxes disagreeing face by face across a shared boundary.
//! It does *not* catch a model that is internally consistent and physically wrong — publish
//! a power where a joule was wanted and both sides agree perfectly about a number off by
//! `1/dt`. For that, check against something the code did not compute: a closed form, an
//! exact limit, or a convergence rate.
//!
//! `cargo run --example agents_quickstart` is all of the above as a running program,
//! including a deliberate leak so the failure is visible. `AGENTS.md` in the repository is
//! the one-page version.
//!
//! # The dependency rule
//!
//! ```text
//! dualis-units       no dependencies but glam and serde
//! dualis-core        depends on units          the kernel: what evolves, what it conserves
//! dualis-optics      depends on core   ┐
//! dualis-thermal     depends on core   │
//! dualis-mechanics   depends on core   ├ one crate per physics, and none knows another
//! dualis-acoustic    depends on core   │
//! dualis-molecular   depends on core   │
//! dualis-electrical  depends on core   ┘
//! dualis-scene       depends on core           where things are, and what a run looks like
//! dualis-view        depends on scene          how to draw that, chosen by the data's shape
//! dualis             depends on all of them
//! ```
//!
//! None of the eight domains knows about any of the others. They meet on the kernel's
//! [`Exchange`](dualis_core::Exchange), and each one that arrived left the others
//! untouched — which is the claim the split was made to test, now held six times.
//!
//! [`scene`] and [`view`] are layers up rather than domains, and they are bound by the same rule
//! from the other side: neither names a domain. A physics that arrives tomorrow is captured
//! without `scene` being edited, and drawn without `view` being edited, because the scene asks
//! each domain what it *offers* and the view dispatches on the shape of what came back.
//! `ARCHITECTURE.md` is the long version.

// Every public item carries a doc comment. Denied rather than warned: a public physics API
// whose `Length::mm` shows a blank summary in rustdoc is documented in the sense that a
// paragraph exists somewhere, and not in the sense a reader needs.
#![deny(missing_docs)]
pub use dualis_acoustic as acoustic;
pub use dualis_core as core;
pub use dualis_elastic as elastic;
pub use dualis_electrical as electrical;
pub use dualis_mechanics as mechanics;
pub use dualis_molecular as molecular;
pub use dualis_optics as optics;
pub use dualis_porous as porous;
pub use dualis_scene as scene;
pub use dualis_thermal as thermal;
pub use dualis_units as units;
pub use dualis_view as view;

/// Everything most simulations need, in one `use`.
pub mod prelude {
    pub use dualis_acoustic::{
        impedance, reflection_coefficient, End, Hall, Impedance, Room, Tube,
    };
    pub use dualis_core::conserved::quantity;
    pub use dualis_core::{
        audit, basis_for, oriented_against, reflect, velocity_verlet, Domain, Dynamics, Exchange,
        Flux, Integrator, Interface, Kind, Ledger, Motion, Newtonian, Report, Rng, ScalarField,
        Schedule, Simulation, State, Strobe, Substance, VectorField, Violation,
    };
    // `Substance` was in the prelude and could not be *built* from it: `Substance::bulk` leaves
    // `thermal: None`, which `LumpedMass` rightly refuses to step, and the three types needed to
    // supply one were not exported. A consumer wanting a material the catalogue does not carry
    // had to reach through two module paths to say so.
    pub use dualis_core::substance::{AcousticProps, MechanicalProps, ThermalProps};
    // `Reading`, `Tolerances` and `Bodies` were reachable only through `dualis::core`. Each is a
    // type a consumer meets while building a frame, setting an audit or writing a domain, which
    // is what a prelude is for.
    pub use dualis_core::{Bodies, Ensemble, Estimate, Pose, Reading, Tolerances};
    // `Block` rather than `Body`, because `dualis-mechanics` had that name first for a rigid
    // one. Two domains both wanted it, and the prelude is where that shows: a consumer who
    // imports everything gets one name for two quite different things, and the compiler catches
    // it only because both are exported here.
    pub use dualis_elastic::{Block, Elastic, Face};
    pub use dualis_electrical::{Conductor, Winding};
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
    pub use dualis_porous::{Basket, Bed, Grind, Liquid, Observable, Puck, Shot};
    pub use dualis_scene::{
        capture, sample_field, settle_framing, Extent, Frame, Panel, PanelData, Placement,
    };
    pub use dualis_thermal::{
        Bar1D, Environment, LumpedMass, Node, Solid3D, SteadyState, ThermalNetwork, HEAT,
    };
    pub use dualis_units::{
        AccelerationVec, Area, Conductance, Conductivity, Current, CurrentDensity, Damping,
        Density, ElectricField, Energy, Force, ForceVec, Frequency, HeatCapacity, Irradiance,
        Length, LengthVec, Mass, Momentum, MomentumVec, Power, Pressure, Resistance, Resistivity,
        SpecificHeat, Stiffness, Temperature, ThermalConductivity, ThermalExpansion, Time,
        Velocity, VelocityVec, Voltage, Volume, G0,
    };
}
