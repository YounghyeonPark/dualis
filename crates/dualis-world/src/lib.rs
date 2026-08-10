//! A world described as data, run against the dualis SDK, and drawn.
//!
//! This is the workspace's first consumer, and its first job is not to be a good application.
//! It is to be an *outside* user of the library — one that reaches for the API the way a
//! stranger would rather than the way its author remembers it — and to write down every place
//! that turns out to be awkward. A library with no consumers is a library whose ergonomics
//! nobody has measured.
//!
//! Findings are collected in `FRICTION.md` beside this crate. Sixteen of the twenty-two are
//! fixed — this crate is the record of what the API was like before, and the reason it changed.
//!
//! The layers it once carried are libraries now. `dualis-scene` owns capture, and this crate is
//! left with what an *application* actually is: a file format, the domain types that format names,
//! and one place saying how far each field extends. A consumer that wants to draw a run no longer
//! has to reach into a binary that is `publish = false` to do it.

#![deny(missing_docs)]

use dualis::prelude::*;
// `Reading` belongs below the layers, because it crosses one: a domain produces it and a view
// consumes it, and neither should have to depend on this application to name the type.
pub use dualis::core::Reading;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

// `Room` is in the prelude now (it was not, though `Tube` was — FRICTION 5). Still aliased,
// because this crate has a `DomainSpec::Room` variant of its own and that name is the right
// one on both sides. That collision is the app's problem and not the library's.
use dualis::molecular as dualis_molecular;
use dualis::prelude::Room as AcousticRoom;

pub mod beam;
pub mod heater;
pub mod light;

use beam::Beam;
use heater::Heater;
use light::Light;

/// What to simulate, in a form that can be written down rather than compiled in.
///
/// `deny_unknown_fields`, on this and on every type below it, and the reason is worth stating.
/// `serde` discards unknown keys by default, which is right for a wire protocol that must
/// tolerate a newer peer and wrong for a document somebody saved. A key this format does not
/// know is a typo or a version skew, and either way discarding it makes the file say something
/// the code will not do — silently, with the field falling back to its `Default`.
///
/// That happened. `main.rs`'s own built-in scene kept the pre-`release` spelling for two
/// commits after the format changed, and nothing failed: `mode` and `amplitude_pa` were
/// dropped, `Release::default()` filled in the same (1,1) mode at 1 Pa, and the output was
/// byte-identical. Editing those keys was a no-op that reported success.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Scene {
    /// Shown in the output; has no effect on the physics.
    pub title: String,
    /// How the domains interact. See [`ScheduleSpec`].
    #[serde(default)]
    pub schedule: ScheduleSpec,
    /// The domains, in declaration order — which is execution order for the staggered
    /// schedules, so it is part of the physics and not a formatting choice.
    pub domains: Vec<DomainSpec>,
    /// How long to run, in seconds.
    pub duration_s: f64,
    /// How many frames to capture over that duration.
    pub frames: usize,
    /// The relative conservation drift the run may accumulate before it is refused.
    ///
    /// Exposed because it is a property of the scene and not of the engine: a scene with a
    /// dissipative boundary legitimately drifts where a closed one does not.
    #[serde(default = "default_tolerance")]
    pub conservation_tolerance: f64,
}

fn default_tolerance() -> f64 {
    1e-6
}

/// Which coupling scheme to run the domains under.
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum ScheduleSpec {
    /// One pass, no feedback expected.
    OneWay,
    /// One pass in declaration order, each domain seeing the earlier ones' output.
    Staggered,
    /// Staggered, but every domain substeps to its own stability limit.
    ///
    /// The default, and it should be: an application picks a frame interval for reasons of
    /// its own — thirty a second, say — and that has nothing to do with any domain's CFL
    /// limit. Under [`ScheduleSpec::Staggered`] a frame interval larger than the limit is
    /// silently unstable. Under this one it is subcycled.
    #[default]
    Multirate,
}

/// One domain in a scene.
///
/// Still an enum rather than an open registry — a third party cannot add a variant — but the
/// kernel no longer forces that: `Simulation::with_boxed` takes a domain chosen at run time,
/// so `DomainSpec::build` is the only place that knows the types, and it hands back a
/// `Box<dyn Domain>` like anything else would.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case", deny_unknown_fields)]
pub enum DomainSpec {
    /// A two-dimensional box of air with rigid walls, released in a standing mode.
    Room {
        /// Domain name, and the handle the renderer uses to find it again.
        name: String,
        /// Across.
        width_m: f64,
        /// Up.
        height_m: f64,
        /// Grid resolution across the width.
        cells_across: usize,
        /// How the field starts. Defaults to the (1,1) mode at 1 Pa.
        #[serde(default)]
        release: Release,
    },
    /// A heat source with a finite tank, defined in this crate rather than the library.
    ///
    /// The publisher half of a coupled scene. See [`heater::Heater`] for why it is written
    /// here: a domain the library already provides tests the constructors and nothing else.
    Heater {
        /// Domain name.
        name: String,
        /// Element power.
        watts: f64,
        /// Joules it has to spend before it goes quiet.
        reserve_j: f64,
    },
    /// A beam that heats *where it lands*, over a shared boundary.
    ///
    /// The spatial publisher. `faces` has to equal the bar's `cells`: both sides build their
    /// own [`Interface`] and the kernel refuses a flux whose face count disagrees, naming
    /// both numbers. Stated twice because nothing derives one from the other — see
    /// `FRICTION.md`, finding 9.
    Beam {
        /// Domain name.
        name: String,
        /// The boundary to publish onto. Must match the bar's `exposes`.
        onto: String,
        /// Faces the boundary is cut into. Must equal the bar's `cells`.
        faces: usize,
        /// Area of one face.
        face_area_mm2: f64,
        /// Beam power.
        watts: f64,
        /// Joules it has to spend.
        reserve_j: f64,
        /// Gaussian waist, as a fraction of the boundary's span.
        waist_fraction: f64,
    },
    /// Bodies under their own gravity: a central mass with satellites on circular orbits.
    ///
    /// `dualis-mechanics`. Not a field — a countable number of things at places — so it is
    /// drawn as dots and `Domain::as_field` rightly declines to invent a continuum for it.
    Orbit {
        /// Domain name.
        name: String,
        /// The mass everything orbits, in kilograms.
        central_kg: f64,
        /// One satellite per radius, each started at the circular speed for that radius and
        /// spaced evenly in angle so they do not begin on top of each other.
        radii_m: Vec<f64>,
        /// Inclination of each orbit to the reference plane, in degrees. Repeats or pads
        /// with zero. Flat orbits make a flat picture and waste the third axis.
        #[serde(default)]
        inclinations_deg: Vec<f64>,
        /// Mass of each satellite. Small against the central one, or "circular" is a lie.
        satellite_kg: f64,
    },
    /// A ball bouncing on a floor through a penalty contact, losing energy to its dashpot.
    Bounce {
        /// Domain name.
        name: String,
        /// Drop height.
        drop_m: f64,
        /// Ball mass.
        mass_kg: f64,
        /// Contact stiffness, in N/m.
        stiffness: f64,
        /// Contact damping, in N·s/m. Zero bounces forever.
        damping: f64,
    },
    /// A Lennard-Jones fluid in a periodic box.
    ///
    /// `dualis-molecular`. Drawn as a slab: every atom projected onto the x-y plane, coloured
    /// by speed, which is what a molecular-dynamics snapshot conventionally shows.
    Atoms {
        /// Domain name.
        name: String,
        /// Unit cells per side; the count is `4·cells³`.
        cells: usize,
        /// Reduced number density, `ρ*`.
        density: f64,
        /// Reduced temperature to start at, `T*`.
        temperature: f64,
        /// If set, a Langevin bath holds it at this reduced temperature.
        #[serde(default)]
        thermostat_t: Option<f64>,
        /// Seed. Nothing here consults a clock.
        seed: u64,
    },
    /// A lumped thermal mass: one temperature, losing heat to still air.
    ///
    /// The consumer a dissipating domain needs. A `bounce` publishes its dashpot's losses
    /// onto the heat channel, and without something to take them the kernel refuses the step
    /// — correctly, because joules that left one domain and arrived nowhere are joules that
    /// went missing. It has no field, so it shows up in the numbers and not in the picture.
    Lump {
        /// Domain name.
        name: String,
        /// Volume of the thing being warmed.
        volume_cm3: f64,
        /// Conduction path length, for the Biot number.
        thickness_mm: f64,
        /// Starting temperature.
        initial_c: f64,
        /// Air temperature it loses to.
        ambient_c: f64,
        /// Surface area exposed to that air.
        area_cm2: f64,
    },
    /// A lamp on a coated surface: real spectra deciding how much becomes heat.
    ///
    /// `dualis-optics`. The absorbed fraction is the overlap of a blackbody at `colour_k`
    /// with the coating's absorptance across the visible range, so changing the colour
    /// temperature changes it — a cooler lamp puts more of its output where the coating is
    /// worse. It has no field, so it shows in the numbers and the bar shows in the picture.
    Light {
        /// Domain name.
        name: String,
        /// Lamp power over the visible range.
        watts: f64,
        /// Colour temperature of the blackbody, in kelvin. 3200 is tungsten.
        colour_k: f64,
        /// The surface. Only `aluminium` for now, whose reflectance falls off in the blue.
        finish: String,
        /// Joules it may spend before it goes dark.
        reserve_j: f64,
    },
    /// A one-dimensional conducting bar.
    Bar {
        /// Domain name, and the handle the renderer uses to find it again.
        name: String,
        /// Total length.
        length_mm: f64,
        /// How many cells to divide it into.
        cells: usize,
        /// Cross-sectional area.
        area_mm2: f64,
        /// Starting temperature, uniform, in celsius.
        initial_c: f64,
        /// If set, the bar exposes a boundary of this name that a beam can land on. One
        /// face per cell, which is the bar's own choice and the reason a beam has to be
        /// told the cell count separately.
        #[serde(default)]
        exposes: Option<Boundary>,
    },
    /// A room with a **ceiling**: the wave equation in three dimensions.
    ///
    /// `dualis-acoustic`'s `Hall`. `DomainSpec::Room` is a floor plan and does not have the
    /// vertical modes at all — not less accurately, at all — and a 2.4 m ceiling puts the first
    /// one at 71 Hz, well inside the range a room is judged on.
    ///
    /// Cells are cubes of the spacing the width sets, so height and depth are quantised. A zero
    /// depth is one node and reduces exactly to a `Room`.
    Hall {
        /// Domain name.
        name: String,
        /// Across.
        width_m: f64,
        /// Up.
        height_m: f64,
        /// And back. Quantised to whole cells like the height.
        depth_m: f64,
        /// Node count along the width, which sets the spacing for all three axes.
        nodes_across: usize,
        /// Which rigid-wall mode to release, as `[a, b, c]`.
        mode: [u32; 3],
        /// The amplitude to release it at.
        amplitude_pa: f64,
    },
    /// A three-dimensional conducting block.
    ///
    /// `dualis-thermal`'s `Solid3D`. The first domain in this format with a field that is
    /// genuinely a volume, and the reason the capture layer grew a third axis: everything before
    /// it was a line or a plane, so nothing had ever noticed that `Extent` could not describe a
    /// solid.
    ///
    /// Cells are cubes of one spacing, so the block's *shape* is the three counts and its size
    /// is `counts × cell_mm`. That is the domain's own restriction and the format does not
    /// paper over it: an anisotropic cell would make the stability limit anisotropic.
    Block {
        /// Domain name.
        name: String,
        /// Cells along x, y and z.
        cells: [usize; 3],
        /// The side of one cubic cell.
        cell_mm: f64,
        /// Starting temperature, uniform, in celsius.
        initial_c: f64,
        /// A cell to warm at the start, and by how much, so there is a hot spot to watch spread.
        ///
        /// **A statement about the initial state, not a delivery of heat.** It moves what the
        /// block holds and not what it has absorbed, so the audit's opening balance includes it
        /// — which is the honest bookkeeping and is why it is separate from any source.
        #[serde(default)]
        hot_spot: Option<HotSpot>,
    },
    /// A copper winding dissipating `I²R`, which is where a motor's heat actually comes from.
    ///
    /// `dualis-electrical`. Every other source in this format states its watts; this one
    /// *computes* them from a length of wire, a cross-section and a current, so getting the
    /// geometry wrong makes the number wrong in a way a closed form can catch. A stated number
    /// cannot be wrong, which is another way of saying it is not a model.
    ///
    /// Resistance rises with temperature — 0.393% per kelvin for copper — but the winding is not
    /// told what temperature it is by the simulation: a domain cannot read another's state
    /// inside the step loop. `at_c` is where you say it, and the feedback that causes thermal
    /// runaway is deliberately outside what this format expresses. See `dualis-electrical`'s
    /// module documentation.
    Winding {
        /// Domain name.
        name: String,
        /// Length of wire.
        length_m: f64,
        /// Cross-section, in square millimetres. 0.35 is roughly AWG 22.
        cross_section_mm2: f64,
        /// Current through it.
        amps: f64,
        /// The temperature its resistance is evaluated at.
        at_c: f64,
        /// Joules it may dissipate before it goes quiet. Not optional: without it the winding
        /// supplies energy from nowhere and the audit cannot see that, so the domain refuses.
        reserve_j: f64,
        /// If set, the winding's temperature is refreshed from a thermal network node between
        /// frames — `"network/node"`, e.g. `"motor/winding"`.
        ///
        /// This is the electro-thermal feedback, closed **here** because it cannot be closed
        /// inside the step loop: a domain has no way to read another's state, and that is the
        /// property the crate split exists to hold. The caller between frames can see both, so
        /// the loop lives at the application level or nowhere.
        ///
        /// It is not free. The temperature is refreshed once per *frame*, not once per substep,
        /// so the feedback lags by up to one frame interval. Scene 13 measures what that costs.
        #[serde(default)]
        tracks: Option<String>,
    },
    /// Several lumped bodies joined by conductances: junction, case, ambient.
    ///
    /// The one shape a `lump` cannot express. A `lump` reports the temperature of the whole
    /// thing, and the number a designer needs is the *winding*, which is hotter than the case by
    /// however much the joint between them resists. A network carries that drop explicitly.
    ///
    /// This is also the only spec whose parts refer to each other by name, so it is where the
    /// library's handle-based API meets a file. [`ThermalNetwork`] deliberately does not accept
    /// names for links — a link naming a node that does not exist balances the books perfectly
    /// while modelling something else — so the resolution happens here, once, where the error can
    /// quote the file's own vocabulary back.
    Network {
        /// Domain name.
        name: String,
        /// The bodies. Declared before the links, which refer to them by `name`.
        nodes: Vec<NetworkNode>,
        /// The conductances between them, in W/K.
        links: Vec<NetworkLink>,
        /// Which node heat arriving on the bus lands in — a winding, usually.
        absorbing: String,
    },
}

/// One body in a [`DomainSpec::Network`].
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NetworkNode {
    /// What the links call it. Unique within the network.
    pub name: String,
    /// One of `copper`, `aluminium`, `electrical_steel`, `fr4`, `pla`. Mixing them is the
    /// point: a motor is not a billet of its main metal.
    pub material: String,
    /// Volume of this body.
    pub volume_cm3: f64,
    /// Conduction path length within the body, for its Biot number.
    pub thickness_mm: f64,
    /// Starting temperature.
    pub initial_c: f64,
    /// If absent, the body is interior: it conducts to its neighbours and loses to nothing.
    #[serde(default)]
    pub loses_to: Option<NetworkAmbient>,
}

/// Still air a node loses heat to.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NetworkAmbient {
    /// Air temperature.
    pub ambient_c: f64,
    /// Surface area exposed to it.
    pub area_cm2: f64,
}

/// A conductance between two nodes of a [`DomainSpec::Network`].
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NetworkLink {
    /// Name of one node.
    pub from: String,
    /// Name of the other. Order does not matter; a conductance is symmetric.
    pub to: String,
    /// The conductance, in watts per kelvin. Repeats between the same pair add, as parallel
    /// paths do.
    pub w_per_k: f64,
}

/// How a room's field is set up before the clock starts.
///
/// A mode is the case with a closed form to check against, so it is what the tests use. A
/// pulse is the case worth looking at: it has no standing shape, so it travels, reflects off
/// the walls and interferes with itself, which is what a room actually does to a sound.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "as", rename_all = "kebab-case", deny_unknown_fields)]
pub enum Release {
    /// One standing mode, exactly. `cos(nx πx/Lx) cos(ny πy/Ly)`.
    Mode {
        /// Half-wavelengths across the width.
        nx: u32,
        /// Half-wavelengths up the height.
        ny: u32,
        /// Peak pressure, in pascals.
        amplitude_pa: f64,
    },
    /// A Gaussian bump, at rest.
    ///
    /// Released from rest it splits into outgoing waves in every direction, each carrying
    /// half the amplitude — worth knowing before reading a height off one of them.
    Pulse {
        /// Where, across.
        x_m: f64,
        /// Where, up.
        y_m: f64,
        /// Gaussian radius, in metres.
        radius_m: f64,
        /// Peak pressure, in pascals.
        amplitude_pa: f64,
    },
}

impl Default for Release {
    fn default() -> Release {
        Release::Mode {
            nx: 1,
            ny: 1,
            amplitude_pa: 1.0,
        }
    }
}

/// A boundary a bar offers for something else to publish onto.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Boundary {
    /// Both sides have to use the same name.
    pub name: String,
    /// Area of one face.
    pub face_area_mm2: f64,
}

impl DomainSpec {
    /// The name this domain will answer to.
    pub fn name(&self) -> &str {
        match self {
            DomainSpec::Room { name, .. }
            | DomainSpec::Bar { name, .. }
            | DomainSpec::Block { name, .. }
            | DomainSpec::Hall { name, .. }
            | DomainSpec::Heater { name, .. }
            | DomainSpec::Beam { name, .. }
            | DomainSpec::Orbit { name, .. }
            | DomainSpec::Bounce { name, .. }
            | DomainSpec::Atoms { name, .. }
            | DomainSpec::Lump { name, .. }
            | DomainSpec::Network { name, .. }
            | DomainSpec::Winding { name, .. }
            | DomainSpec::Light { name, .. } => name,
        }
    }

    /// Construct the domain this describes.
    ///
    /// Returns a box, which is what a builder chosen at run time can produce, and
    /// `Simulation::with_boxed` takes one. The names go straight through: they are `String`s
    /// out of a file and the constructors take `impl Into<String>`, so nothing is leaked and
    /// nothing is pretending to be `'static`.
    ///
    /// Fallible, and it has to be. An unrecognised `finish` used to return a lamp of zero watts
    /// with the mistake written into its *name* — which nothing printed, because a lamp has no
    /// panel. Worse, the early return skipped `with_reserve`, so the reserve stayed infinite,
    /// so `Light::ledger` reported nothing, so the audit had nothing to compare: the scene ran
    /// green at `conservation_tolerance(0.0)`, the strictest setting expressible, with the lamp
    /// doing nothing at all. One character — `aluminium` against `aluminum` — turned off the
    /// check this library exists for.
    pub fn build(&self) -> Result<Box<dyn Domain>, String> {
        let domain: Box<dyn Domain> = match self {
            DomainSpec::Room {
                name,
                width_m,
                height_m,
                cells_across,
                release,
            } => {
                let room = AcousticRoom::of_air(
                    name.clone(),
                    Length::m(*width_m),
                    Length::m(*height_m),
                    *cells_across,
                );
                Box::new(match release {
                    Release::Mode {
                        nx,
                        ny,
                        amplitude_pa,
                    } => room.released_in_mode(*nx, *ny, Pressure::from_si(*amplitude_pa)),
                    Release::Pulse {
                        x_m,
                        y_m,
                        radius_m,
                        amplitude_pa,
                    } => {
                        let (cx, cy, r, a) = (*x_m, *y_m, radius_m.max(1e-9), *amplitude_pa);
                        room.released_from(move |x, y| {
                            let (dx, dy) = (x.to_si() - cx, y.to_si() - cy);
                            Pressure::from_si(a * (-(dx * dx + dy * dy) / (r * r)).exp())
                        })
                    }
                })
            }
            DomainSpec::Heater {
                name,
                watts,
                reserve_j,
            } => Box::new(Heater::new(name.clone(), *watts, *reserve_j)),
            DomainSpec::Beam {
                name,
                onto,
                faces,
                face_area_mm2,
                watts,
                reserve_j,
                waist_fraction,
            } => Box::new(Beam::new(
                name.clone(),
                Interface::uniform(onto.clone(), *faces, Area::from_si(face_area_mm2 * 1e-6)),
                *watts,
                *reserve_j,
                *waist_fraction,
            )),
            DomainSpec::Orbit {
                name,
                central_kg,
                radii_m,
                inclinations_deg,
                satellite_kg,
            } => {
                let mut bodies = vec![Body::new(
                    Mass::kg(*central_kg),
                    LengthVec::m(0.0, 0.0, 0.0),
                    VelocityVec::ZERO,
                )];
                for (k, r) in radii_m.iter().enumerate() {
                    // Evenly spaced in angle, and each at the circular speed the library
                    // computes for its radius — so an ellipse in the picture is the
                    // integrator's opinion and not the initial condition's.
                    let a = k as f64 * std::f64::consts::TAU / radii_m.len().max(1) as f64;
                    let v = NBody::circular_speed(Mass::kg(*central_kg), Length::m(*r)).to_si();
                    // Tilt the orbit by rotating both position and velocity about the x axis,
                    // which keeps the speed exactly circular — inclining the position alone
                    // would put the body on an ellipse and call it an integrator error.
                    let inc = inclinations_deg.get(k).copied().unwrap_or(0.0).to_radians();
                    let (si, ci) = (inc.sin(), inc.cos());
                    let (px, py) = (r * a.cos(), r * a.sin());
                    let (vx, vy) = (-v * a.sin(), v * a.cos());
                    bodies.push(Body::new(
                        Mass::kg(*satellite_kg),
                        LengthVec::m(px, py * ci, py * si),
                        VelocityVec::from_si(glam::DVec3::new(vx, vy * ci, vy * si)),
                    ));
                }
                Box::new(NBody::new(name.clone(), &bodies))
            }
            DomainSpec::Bounce {
                name,
                drop_m,
                mass_kg,
                stiffness,
                damping,
            } => Box::new(ContactSystem::new(
                name.clone(),
                &[Body::new(
                    Mass::kg(*mass_kg),
                    LengthVec::m(0.0, 0.0, *drop_m),
                    VelocityVec::ZERO,
                )],
                AccelerationVec::from_si(-glam::DVec3::Z * G0.to_si()),
                Ground::floor(),
                Stiffness::from_si(*stiffness),
                Damping::from_si(*damping),
            )),
            DomainSpec::Atoms {
                name,
                cells,
                density,
                temperature,
                thermostat_t,
                seed,
            } => {
                let lj = LennardJones::reduced();
                let fluid = Fluid::lattice(
                    name.clone(),
                    lj,
                    dualis_molecular::unit_mass(),
                    *cells,
                    *density,
                )
                .thermalised(
                    dualis_molecular::temperature_from_reduced(*temperature, &lj),
                    *seed,
                );
                Box::new(match thermostat_t {
                    Some(t) => fluid.with_thermostat(Thermostat::Langevin {
                        target: dualis_molecular::temperature_from_reduced(*t, &lj),
                        damping: 2.0,
                    }),
                    None => fluid,
                })
            }
            DomainSpec::Light {
                name,
                watts,
                colour_k,
                finish,
                reserve_j,
            } => Box::new(
                Light::new(
                    name.clone(),
                    *watts,
                    *colour_k,
                    // One surface so far, and the scene names it anyway: a field that can
                    // hold only one value today is a field that says what would change.
                    match finish.as_str() {
                        "aluminium" => light::aluminium_mirror(),
                        other => {
                            return Err(format!(
                                "{name}: unknown finish {other:?}; known finishes are                                  \"aluminium\""
                            ))
                        }
                    },
                )
                .with_reserve(*reserve_j),
            ),
            DomainSpec::Lump {
                name,
                volume_cm3,
                thickness_mm,
                initial_c,
                ambient_c,
                area_cm2,
            } => Box::new(LumpedMass::new(
                name.clone(),
                Substance::aluminium_6061(),
                Volume::from_si(volume_cm3 * 1e-6),
                Length::mm(*thickness_mm),
                Temperature::celsius(*initial_c),
                Environment::still_air(
                    Temperature::celsius(*ambient_c),
                    Area::from_si(area_cm2 * 1e-4),
                ),
            )),
            DomainSpec::Winding {
                name,
                length_m,
                cross_section_mm2,
                amps,
                at_c,
                reserve_j,
                tracks: _,
            } => Box::new(
                dualis::electrical::Winding::of_copper(
                    name.clone(),
                    Length::m(*length_m),
                    cross_section_mm2 * 1e-6,
                    Temperature::celsius(*at_c),
                )
                .driven_at(Current::a(*amps))
                .with_reserve(*reserve_j),
            ),
            DomainSpec::Network {
                name,
                nodes,
                links,
                absorbing,
            } => {
                let mut network = ThermalNetwork::new(name.clone());
                for n in nodes {
                    let substance = match n.material.as_str() {
                        "copper" => Substance::copper(),
                        "aluminium" => Substance::aluminium_6061(),
                        "electrical_steel" => Substance::electrical_steel(),
                        "fr4" => Substance::fr4(),
                        "pla" => Substance::pla(),
                        other => {
                            return Err(format!(
                                "{name}/{}: unknown material {other:?}; known materials are \
                                 \"copper\", \"aluminium\", \"electrical_steel\", \"fr4\", \"pla\"",
                                n.name
                            ))
                        }
                    };
                    let volume = Volume::from_si(n.volume_cm3 * 1e-6);
                    let thickness = Length::mm(n.thickness_mm);
                    let initial = Temperature::celsius(n.initial_c);
                    match &n.loses_to {
                        Some(air) => network.node_losing_to(
                            n.name.clone(),
                            substance,
                            volume,
                            thickness,
                            initial,
                            Environment::still_air(
                                Temperature::celsius(air.ambient_c),
                                Area::from_si(air.area_cm2 * 1e-4),
                            ),
                        ),
                        None => {
                            network.node(n.name.clone(), substance, volume, thickness, initial)
                        }
                    };
                }

                // The resolver the library declines to provide, and the reason it declines: here
                // a name that is not a node is a *file* mistake, so the message can list what the
                // file did define. Inside the library the same lookup would have to fail silently
                // or invent an error vocabulary out of strings it did not choose.
                let find = |who: &str| -> Result<_, String> {
                    network.node_named(who).ok_or_else(|| {
                        let known: Vec<&str> = nodes.iter().map(|n| n.name.as_str()).collect();
                        format!(
                            "{name}: no node named {who:?}; this network defines {}",
                            known.join(", ")
                        )
                    })
                };
                let mut resolved = Vec::with_capacity(links.len());
                for l in links {
                    resolved.push((find(&l.from)?, find(&l.to)?, l.w_per_k));
                }
                let sink = find(absorbing)?;

                // Applied after resolution, so a bad name is reported before a half-built network
                // exists. `link` and `absorbing` both return `Violation`, which is not this
                // function's error type — a self-link, a negative conductance or a node from
                // another network. All three are file mistakes too, so they are worth the same
                // treatment rather than an `unwrap`.
                for (a, b, w) in resolved {
                    network
                        .link(a, b, Conductance::w_per_k(w))
                        .map_err(|e| format!("{name}: {e}"))?;
                }
                network
                    .absorbing(sink)
                    .map_err(|e| format!("{name}: {e}"))?;
                Box::new(network)
            }
            DomainSpec::Bar {
                name,
                length_mm,
                cells,
                area_mm2,
                initial_c,
                exposes,
            } => {
                let bar = Bar1D::new(
                    name.clone(),
                    Substance::aluminium_6061(),
                    *cells,
                    Length::mm(length_mm / *cells as f64),
                    Area::from_si(area_mm2 * 1e-6),
                    Temperature::celsius(*initial_c),
                );
                Box::new(match exposes {
                    Some(b) => bar.exposing(b.name.clone(), Area::from_si(b.face_area_mm2 * 1e-6)),
                    None => bar,
                })
            }
            DomainSpec::Hall {
                name,
                width_m,
                height_m,
                depth_m,
                nodes_across,
                mode,
                amplitude_pa,
            } => Box::new(
                dualis::acoustic::Hall::of_air(
                    name.clone(),
                    Length::m(*width_m),
                    Length::m(*height_m),
                    Length::m(*depth_m),
                    *nodes_across,
                )
                .released_in_mode(
                    (mode[0], mode[1], mode[2]),
                    dualis::units::Pressure::from_si(*amplitude_pa),
                ),
            ),
            DomainSpec::Block {
                name,
                cells,
                cell_mm,
                initial_c,
                hot_spot,
            } => {
                let mut block = dualis::thermal::Solid3D::new(
                    name.clone(),
                    Substance::aluminium_6061(),
                    (cells[0], cells[1], cells[2]),
                    Length::mm(*cell_mm),
                    Temperature::celsius(*initial_c),
                );
                if let Some(spot) = hot_spot {
                    block.set_temperature(
                        spot.at[0],
                        spot.at[1],
                        spot.at[2],
                        Temperature::celsius(initial_c + spot.above_k),
                    );
                }
                Box::new(block)
            }
        };
        Ok(domain)
    }
}

/// One cell of a [`DomainSpec::Block`], warmed at the start.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HotSpot {
    /// Which cell, as `[i, j, k]`. Out of range is ignored by the domain, which is what a caller
    /// writing a spot near a face means.
    pub at: [usize; 3],
    /// How much hotter than the block started, in kelvin.
    pub above_k: f64,
}

/// A scene that has been checked and turned into a runnable simulation.
pub struct World {
    scene: Scene,
    sim: Simulation,
}

impl World {
    /// Build a simulation from a scene.
    ///
    /// Fails only on a scene that cannot describe a simulation at all — no domains, a
    /// non-positive duration, no frames. Physical nonsense inside a domain is the domain's
    /// business and is reported by the audit at run time, which is where it belongs.
    pub fn build(scene: Scene) -> Result<World, String> {
        if scene.domains.is_empty() {
            return Err("a scene needs at least one domain".into());
        }
        // NaN spelled out rather than hidden in a negated comparison. A duration that is
        // not a number reaches `advance` as a step of NaN and poisons every field silently,
        // so it is worth refusing here where the message can say which field was wrong.
        if scene.duration_s <= 0.0 || scene.duration_s.is_nan() {
            return Err(format!(
                "duration must be positive, got {}",
                scene.duration_s
            ));
        }
        if scene.frames == 0 {
            return Err("a scene needs at least one frame".into());
        }

        // A `tracks` naming something that is not there would produce a winding whose
        // temperature never moves — which is indistinguishable from a correct run at a constant
        // temperature, and is therefore exactly the shape of failure this crate keeps finding:
        // not a wrong answer, an absent one that looks right. Refused here, where the message
        // can list what the scene does define.
        for spec in &scene.domains {
            let DomainSpec::Winding {
                name,
                tracks: Some(path),
                ..
            } = spec
            else {
                continue;
            };
            let Some((net_name, node_name)) = path.split_once('/') else {
                return Err(format!(
                    "{name}: tracks is {path:?}; it should be \"network/node\""
                ));
            };
            let target = scene.domains.iter().find_map(|d| match d {
                DomainSpec::Network { name, nodes, .. } if name == net_name => Some(nodes),
                _ => None,
            });
            let Some(nodes) = target else {
                return Err(format!(
                    "{name}: tracks names the network {net_name:?}, which this scene does not                      define"
                ));
            };
            if !nodes.iter().any(|n| n.name == node_name) {
                let known: Vec<&str> = nodes.iter().map(|n| n.name.as_str()).collect();
                return Err(format!(
                    "{name}: {net_name:?} has no node {node_name:?}; it has {}",
                    known.join(", ")
                ));
            }
        }

        let mut sim = Simulation::new(match scene.schedule {
            ScheduleSpec::OneWay => Schedule::OneWay,
            ScheduleSpec::Staggered => Schedule::Staggered,
            ScheduleSpec::Multirate => Schedule::Multirate,
        })
        .conservation_tolerance(scene.conservation_tolerance);

        // Two domains under one name is a scene that cannot describe a simulation, which is
        // what this function is for. `Simulation::domain` takes the *first* match, so without
        // this the second domain is never sampled and the first is drawn twice under the
        // second's geometry and label — a 500 C bar reported as 20 C, twice, with no warning.
        for (i, spec) in scene.domains.iter().enumerate() {
            if let Some(earlier) = scene.domains[..i].iter().find(|d| d.name() == spec.name()) {
                return Err(format!(
                    "two domains are both called {:?}; names are how they are looked up,                      so the second would be invisible",
                    earlier.name()
                ));
            }
        }

        for spec in &scene.domains {
            sim = sim.with_boxed(spec.build()?);
        }

        Ok(World { scene, sim })
    }

    /// The scene this was built from.
    pub fn scene(&self) -> &Scene {
        &self.scene
    }

    /// The simulation underneath, for a caller that wants the kernel's own accessors —
    /// `bus`, `ledger`, `domain_as`, `field`.
    pub fn simulation(&self) -> &Simulation {
        &self.sim
    }

    /// Where the clock is.
    pub fn time(&self) -> Time {
        self.sim.time()
    }

    /// Run to the end, capturing frames.
    ///
    /// Returns the frames, or the first [`Violation`] the audit raised — which stops the run,
    /// because a simulation that has stopped conserving is not producing frames worth
    /// drawing.
    pub fn run(&mut self) -> Result<Vec<Frame>, Violation> {
        let dt = Time::from_si(self.scene.duration_s / self.scene.frames as f64);
        let placed = self.placements();
        let mut frames = Vec::with_capacity(self.scene.frames + 1);
        frames.push(dualis::scene::capture(&self.sim, &placed));
        for _ in 0..self.scene.frames {
            self.sim.advance(dt)?;
            self.close_feedback();
            frames.push(dualis::scene::capture(&self.sim, &placed));
        }
        dualis::scene::settle_framing(&mut frames);
        Ok(frames)
    }

    /// Where each domain sits, keyed by the name the simulation knows it under.
    ///
    /// Built once per run rather than per frame: a placement is a property of the scene, and a
    /// scene does not move while it is being simulated. Rebuilding it 240 times would also have
    /// made a placement that *did* drift look like a working feature.
    pub fn placements(&self) -> BTreeMap<String, Placement> {
        self.scene
            .domains
            .iter()
            .map(|spec| (spec.name().to_string(), spec.placement()))
            .collect()
    }

    /// Refresh every tracking winding's temperature from the node it follows.
    ///
    /// **The one place in this application that couples two domains by hand**, and it is worth
    /// being precise about why that is allowed. Domains never read each other *inside* the step
    /// loop — they meet on the bus, which carries amounts and not state. This runs between
    /// frames, in the code that owns the simulation and could rebuild either domain from
    /// scratch, so nothing is being smuggled past the rule.
    ///
    /// It needed `Simulation::domain_as_mut`, which did not exist: a caller could read a domain
    /// and not write one, so this loop was unclosable from anywhere at all. `FRICTION.md` 18.
    ///
    /// Silent when a name does not resolve, and that is the wrong shape — a `tracks` pointing at
    /// a node that is not there produces a winding whose resistance never moves, which looks
    /// exactly like a correct run at a constant temperature. `World::build` validates the target
    /// up front so this cannot be reached with a bad name; the check there is the real one.
    fn close_feedback(&mut self) {
        let targets: Vec<(String, String, String)> = self
            .scene
            .domains
            .iter()
            .filter_map(|spec| match spec {
                DomainSpec::Winding {
                    name,
                    tracks: Some(path),
                    ..
                } => {
                    let (net, node) = path.split_once('/')?;
                    Some((name.clone(), net.to_string(), node.to_string()))
                }
                _ => None,
            })
            .collect();

        for (coil, net_name, node_name) in targets {
            let Some(net) = self.sim.domain_as::<ThermalNetwork>(&net_name) else {
                continue;
            };
            let Some(node) = net.node_named(&node_name) else {
                continue;
            };
            let t = net.temperature(node);
            if let Some(w) = self.sim.domain_as_mut::<dualis::electrical::Winding>(&coil) {
                w.at_temperature(t);
            }
        }
    }
}

/// Re-exported from [`dualis::scene`], which is where they live now.
///
/// They were defined here, because this crate was the only thing that had ever needed them and a
/// type nobody else can reach costs nothing to keep local. That stopped being true the moment
/// `publish = false` was the only reason a consumer could not draw a run: an application is not a
/// place to keep the shape of an answer.
pub use dualis::scene::{Extent, Frame, Panel, PanelData, Placement};

impl DomainSpec {
    /// Where this domain sits, and how much of it to sample.
    ///
    /// **The only place in this application that matches on a domain type to draw one**, and it
    /// is the right place: a scene format is a list of domain kinds already, so knowing them is
    /// its job rather than an admission. Every other such match is gone — the scalars a domain
    /// reports, the bodies it has, the unit its field is in, all of it is the domain's own.
    ///
    /// What is left here is the one thing a domain genuinely cannot answer: **how far a field
    /// extends**. `ScalarField` is a function of position and does not stop anywhere; a field
    /// that knew its own bounds would be a mesh. So the scene says, because the scene is what
    /// wrote the size down in the first place.
    ///
    /// Spelled out variant by variant rather than with a `_`, and that is deliberate. A wildcard
    /// here would give a newly added field domain no extent, and a domain with no extent is not
    /// drawn — the failure would be a missing panel in a report nobody was watching closely,
    /// which is exactly the shape of failure this crate keeps finding.
    pub fn placement(&self) -> Placement {
        match self {
            // A room is a rectangle, sampled at its own cell spacing so the picture has the
            // resolution the simulation does — no more, which would interpolate, and no less,
            // which would alias a mode into a smoother one.
            DomainSpec::Room {
                cells_across,
                width_m,
                height_m,
                ..
            } => {
                let nx = (*cells_across).max(3);
                let ny = ((height_m / width_m) * (nx - 1) as f64).round() as usize + 1;
                Placement::field(Extent::rectangle(
                    Length::m(*width_m),
                    Length::m(*height_m),
                    nx,
                    ny.max(3),
                ))
            }
            // A bar is a line. One row, and the report picks a profile rather than a heatmap
            // from that shape alone.
            DomainSpec::Bar {
                cells, length_mm, ..
            } => Placement::field(Extent::line(Length::mm(*length_mm), (*cells).max(2))),
            // A hall is a volume too, sampled at its own node count — which is the spacing the
            // width sets, applied to all three axes, so the picture has the resolution the
            // simulation does and no more.
            DomainSpec::Hall {
                width_m,
                height_m,
                depth_m,
                nodes_across,
                ..
            } => {
                let nx = (*nodes_across).max(2);
                let dx = width_m / (nx - 1) as f64;
                let nodes = |l: f64| ((l / dx).round().max(0.0) as usize) + 1;
                let (ny, nz) = (nodes(*height_m), nodes(*depth_m));
                Placement::field(Extent::volume(
                    LengthVec::m(
                        (nx - 1) as f64 * dx,
                        (ny - 1) as f64 * dx,
                        (nz - 1) as f64 * dx,
                    ),
                    nx,
                    ny,
                    nz,
                ))
            }
            // A block is a volume, and the extent says so. Sampled at the block's own cell
            // count, which is what `Extent::volume` is for — and which nothing in this format
            // could express until a domain with a 3D field arrived to need it.
            DomainSpec::Block { cells, cell_mm, .. } => Placement::field(Extent::volume(
                LengthVec::from_si(
                    glam::DVec3::new(cells[0] as f64, cells[1] as f64, cells[2] as f64)
                        * cell_mm
                        * 1e-3,
                ),
                cells[0],
                cells[1],
                cells[2],
            )),
            // Bodies, which carry their own positions and need no extent.
            DomainSpec::Orbit { .. } | DomainSpec::Bounce { .. } | DomainSpec::Atoms { .. } => {
                Placement::default()
            }
            // No picture at all: sources, sinks, a lumped mass, a graph of nodes. Their result
            // is a reading, and `Domain::readings` collects it without anybody placing them.
            DomainSpec::Heater { .. }
            | DomainSpec::Beam { .. }
            | DomainSpec::Lump { .. }
            | DomainSpec::Light { .. }
            | DomainSpec::Winding { .. }
            | DomainSpec::Network { .. } => Placement::default(),
        }
    }
}
