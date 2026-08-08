//! A world described as data, run against the dualis SDK, and drawn.
//!
//! This is the workspace's first consumer, and its first job is not to be a good application.
//! It is to be an *outside* user of the library — one that reaches for the API the way a
//! stranger would rather than the way its author remembers it — and to write down every place
//! that turns out to be awkward. A library with no consumers is a library whose ergonomics
//! nobody has measured.
//!
//! Findings are collected in `FRICTION.md` beside this crate. Five of the six are fixed —
//! this crate is the record of what the API was like before, and the reason it changed.

#![deny(missing_docs)]

use dualis::prelude::*;
use serde::{Deserialize, Serialize};

// `Room` is in the prelude now (it was not, though `Tube` was — FRICTION 5). Still aliased,
// because this crate has a `DomainSpec::Room` variant of its own and that name is the right
// one on both sides. That collision is the app's problem and not the library's.
use dualis::molecular as dualis_molecular;
use dualis::prelude::Room as AcousticRoom;

pub mod beam;
pub mod heater;
pub mod light;
pub mod render;

use beam::Beam;
use heater::Heater;
use light::Light;

/// What to simulate, in a form that can be written down rather than compiled in.
#[derive(Clone, Debug, Serialize, Deserialize)]
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
#[serde(rename_all = "kebab-case")]
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
#[serde(tag = "kind", rename_all = "kebab-case")]
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
}

/// How a room's field is set up before the clock starts.
///
/// A mode is the case with a closed form to check against, so it is what the tests use. A
/// pulse is the case worth looking at: it has no standing shape, so it travels, reflects off
/// the walls and interferes with itself, which is what a room actually does to a sound.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "as", rename_all = "kebab-case")]
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
            | DomainSpec::Heater { name, .. }
            | DomainSpec::Beam { name, .. }
            | DomainSpec::Orbit { name, .. }
            | DomainSpec::Bounce { name, .. }
            | DomainSpec::Atoms { name, .. }
            | DomainSpec::Lump { name, .. }
            | DomainSpec::Light { name, .. } => name,
        }
    }

    /// Construct the domain this describes.
    ///
    /// Returns a box, which is what a builder chosen at run time can produce, and
    /// `Simulation::with_boxed` takes one. The names go straight through: they are `String`s
    /// out of a file and the constructors take `impl Into<String>`, so nothing is leaked and
    /// nothing is pretending to be `'static`.
    pub fn build(&self) -> Box<dyn Domain> {
        match self {
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
                    // hold only one value today is a field that says what would change, and
                    // the alternative is a scene format that has to grow a key later.
                    match finish.as_str() {
                        "aluminium" => light::aluminium_mirror(),
                        other => {
                            return Box::new(Light::new(
                                format!("{name} (unknown finish {other:?})"),
                                0.0,
                                *colour_k,
                                light::aluminium_mirror(),
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
        }
    }
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

        let mut sim = Simulation::new(match scene.schedule {
            ScheduleSpec::OneWay => Schedule::OneWay,
            ScheduleSpec::Staggered => Schedule::Staggered,
            ScheduleSpec::Multirate => Schedule::Multirate,
        })
        .conservation_tolerance(scene.conservation_tolerance);

        for spec in &scene.domains {
            sim = sim.with_boxed(spec.build());
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
        let mut frames = Vec::with_capacity(self.scene.frames + 1);
        frames.push(self.capture());
        for _ in 0..self.scene.frames {
            self.sim.advance(dt)?;
            frames.push(self.capture());
        }
        Ok(frames)
    }

    /// Sample every drawable domain onto a grid at the current time.
    ///
    /// Note what is *not* here: any mention of `Room` or `Bar1D`. A domain that has a field
    /// hands one over through `Domain::as_field`, and this samples it. Adding a domain type
    /// to the scene format no longer means editing the renderer — which is what
    /// `ScalarField` was for, and was not reachable until the trait had a way to offer one.
    fn capture(&self) -> Frame {
        let t = self.sim.time();
        let panels = self
            .scene
            .domains
            .iter()
            .filter_map(|spec| match self.sim.field(spec.name()) {
                Some(field) => Some(sample(spec, field, t)),
                // FRICTION 11: `Domain::as_field` covers the domains that *are* fields, and
                // there is no counterpart for the ones that are a countable number of bodies.
                // So a renderer wanting to draw an orbit or a box of atoms is back to
                // downcasting and knowing the type — which is finding 3 again, for the other
                // half of the domains.
                None => self.bodies(spec),
            })
            .collect();
        Frame {
            time_s: t.to_si(),
            panels,
        }
    }

    /// The domains that are bodies rather than fields, as dots.
    fn bodies(&self, spec: &DomainSpec) -> Option<Panel> {
        let (positions, values, bounds, boxed, unit) = match spec {
            DomainSpec::Orbit { name, radii_m, .. } => {
                let n = self.sim.domain_as::<NBody>(name)?;
                let r = radii_m.iter().cloned().fold(1.0f64, f64::max) * 1.2;
                let (mut p, mut v) = (Vec::new(), Vec::new());
                for i in 0..n.count() {
                    let b = n.body(i);
                    let x = b.position.to_si();
                    p.push([x.x, x.y, x.z]);
                    v.push(b.velocity.to_si().length());
                }
                (p, v, [-r, -r, -r, r, r, r], false, "m/s")
            }
            DomainSpec::Bounce { name, drop_m, .. } => {
                let c = self.sim.domain_as::<ContactSystem>(name)?;
                let (mut p, mut v) = (Vec::new(), Vec::new());
                for i in 0..c.count() {
                    let b = c.body(i);
                    let x = b.position.to_si();
                    p.push([x.x, x.y, x.z]);
                    v.push(b.velocity.to_si().length());
                }
                let half = drop_m * 0.45;
                (
                    p,
                    v,
                    [-half, -half, 0.0, half, half, *drop_m * 1.05],
                    true,
                    "m/s",
                )
            }
            DomainSpec::Atoms { name, .. } => {
                let f = self.sim.domain_as::<Fluid>(name)?;
                let l = f.bounds().length;
                let (mut p, mut v) = (Vec::new(), Vec::new());
                for i in 0..f.count() {
                    let x = f.position(i);
                    p.push([x.x, x.y, x.z]);
                    v.push(f.velocity(i).length());
                }
                // The periodic cell is a real boundary, so it gets drawn.
                (p, v, [0.0, 0.0, 0.0, l, l, l], true, "m/s")
            }
            _ => return None,
        };
        Some(Panel {
            name: spec.name().to_string(),
            unit,
            data: PanelData::Points {
                positions,
                values,
                bounds,
                boxed,
            },
        })
    }
}

/// Sample a field over the extent the scene says it occupies.
///
/// The extent has to come from the scene rather than from the field, because `ScalarField`
/// is a function of position and says nothing about where it stops. That is the right
/// division — a field that knew its own bounds would be a mesh — but it does mean the caller
/// supplies them, and here the scene already has them.
fn sample(spec: &DomainSpec, field: &dyn ScalarField, t: Time) -> Panel {
    let (nx, ny, w, h, unit, offset) = match spec {
        DomainSpec::Room {
            cells_across,
            width_m,
            height_m,
            ..
        } => {
            let nx = (*cells_across).max(3);
            let ny = ((height_m / width_m) * (nx - 1) as f64).round() as usize + 1;
            (nx, ny.max(3), *width_m, *height_m, "Pa", 0.0)
        }
        // A heater has no field, so `Domain::as_field` returns `None` and `capture` never
        // reaches here with one. Matching on it anyway rather than an `unreachable!`, because
        // a panic reachable only through a future edit is worse than a panel nobody draws.
        // Sources have no field, and bodies are handled by `World::bodies`, so neither
        // reaches here. Matched rather than left to an `unreachable!`, because a panic a
        // future edit could reach is worse than a panel nobody draws.
        DomainSpec::Heater { .. }
        | DomainSpec::Beam { .. }
        | DomainSpec::Orbit { .. }
        | DomainSpec::Bounce { .. }
        | DomainSpec::Atoms { .. }
        | DomainSpec::Lump { .. }
        | DomainSpec::Light { .. } => (1, 1, 0.0, 0.0, "J", 0.0),
        DomainSpec::Bar {
            cells, length_mm, ..
        } => (
            (*cells).max(2),
            1,
            length_mm * 1e-3,
            0.0,
            "C",
            -273.15, // the bar reports kelvin; a picture of a room wants celsius
        ),
    };

    let mut values = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let x = if nx > 1 {
                w * i as f64 / (nx - 1) as f64
            } else {
                0.0
            };
            let y = if ny > 1 {
                h * j as f64 / (ny - 1) as f64
            } else {
                0.0
            };
            values.push(field.at(LengthVec::m(x, y, 0.0), t) + offset);
        }
    }
    Panel {
        name: spec.name().to_string(),
        unit,
        data: PanelData::Field { nx, ny, values },
    }
}

/// One captured instant: every drawable domain, sampled.
#[derive(Clone, Debug)]
pub struct Frame {
    /// Simulation time, in seconds.
    pub time_s: f64,
    /// One per drawable domain, in scene order.
    pub panels: Vec<Panel>,
}

/// One domain, captured for drawing.
#[derive(Clone, Debug)]
pub struct Panel {
    /// Which domain this came from.
    pub name: String,
    /// What the values mean, for the legend.
    pub unit: &'static str,
    /// The shape of what was captured.
    pub data: PanelData,
}

/// A continuum sampled on a grid, or a finite number of bodies at positions.
///
/// Two shapes rather than one because the domains genuinely are two kinds of thing. A room
/// and a bar are fields: defined everywhere, and a picture of one is a raster. An orbit, a
/// bouncing ball and a box of atoms are not: they are a countable number of things at places,
/// and rasterising them would be inventing a continuum they do not have.
///
/// This is the distinction `ScalarField` cannot express, and it is why `Domain::as_field`
/// returns `None` for three of the five domains rather than something contrived.
#[derive(Clone, Debug)]
pub enum PanelData {
    /// A field, sampled onto a grid.
    Field {
        /// Samples across.
        nx: usize,
        /// Samples up. One for a bar, which has no second dimension.
        ny: usize,
        /// Row-major, `nx * ny` values.
        values: Vec<f64>,
    },
    /// Bodies at positions in space, each carrying a value to colour it by.
    ///
    /// Three dimensions, because the physics is three-dimensional and always was: `NBody`,
    /// `ContactSystem` and `Fluid` all carry `DVec3`. Flattening to a plane was the
    /// *renderer's* simplification, and it threw away a whole axis of every orbit and every
    /// box of atoms. The projection is the picture's business, not the simulation's.
    Points {
        /// Where each body is, in metres.
        positions: Vec<[f64; 3]>,
        /// One per body — a speed, a height, whatever the scene is about.
        values: Vec<f64>,
        /// `[x0, y0, z0, x1, y1, z1]`, the region to draw. Fixed for the whole run so a body
        /// moving is a body moving and not the frame rescaling underneath it.
        bounds: [f64; 6],
        /// Whether to draw the bounding box as a wireframe. True for a periodic cell, which
        /// is a real wall; false for an orbit, whose box is only the extent of the drawing.
        boxed: bool,
    },
}

impl Panel {
    /// The scalar values, whichever shape this is.
    pub fn values(&self) -> &[f64] {
        match &self.data {
            PanelData::Field { values, .. } | PanelData::Points { values, .. } => values,
        }
    }

    /// Grid size, for a field. `None` for points, which have no grid.
    pub fn grid(&self) -> Option<(usize, usize)> {
        match self.data {
            PanelData::Field { nx, ny, .. } => Some((nx, ny)),
            PanelData::Points { .. } => None,
        }
    }
}
