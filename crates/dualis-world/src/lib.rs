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
use dualis::prelude::Room as AcousticRoom;

pub mod beam;
pub mod heater;
pub mod render;

use beam::Beam;
use heater::Heater;

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
        /// Which standing mode to start in, as (nx, ny).
        mode: [u32; 2],
        /// How loud, in pascals.
        amplitude_pa: f64,
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
            | DomainSpec::Beam { name, .. } => name,
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
                mode,
                amplitude_pa,
            } => Box::new(
                AcousticRoom::of_air(
                    name.clone(),
                    Length::m(*width_m),
                    Length::m(*height_m),
                    *cells_across,
                )
                .released_in_mode(
                    mode[0],
                    mode[1],
                    Pressure::from_si(*amplitude_pa),
                ),
            ),
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
            .filter_map(|spec| {
                let field = self.sim.field(spec.name())?;
                Some(sample(spec, field, t))
            })
            .collect();
        Frame {
            time_s: t.to_si(),
            panels,
        }
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
        DomainSpec::Heater { .. } | DomainSpec::Beam { .. } => (1, 1, 0.0, 0.0, "J", 0.0),
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
        nx,
        ny,
        values,
        unit,
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

/// One domain's field, sampled onto a grid for drawing.
#[derive(Clone, Debug)]
pub struct Panel {
    /// Which domain this came from.
    pub name: String,
    /// Samples across.
    pub nx: usize,
    /// Samples up. One for a bar, which has no second dimension.
    pub ny: usize,
    /// Row-major, `nx * ny` values in whatever unit the domain reports.
    pub values: Vec<f64>,
    /// What the values mean, for the legend.
    pub unit: &'static str,
}
