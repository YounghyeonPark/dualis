//! A world described as data, run against the dualis SDK, and drawn.
//!
//! This is the workspace's first consumer, and its first job is not to be a good application.
//! It is to be an *outside* user of the library — one that reaches for the API the way a
//! stranger would rather than the way its author remembers it — and to write down every place
//! that turns out to be awkward. A library with no consumers is a library whose ergonomics
//! nobody has measured.
//!
//! Findings so far are collected in `FRICTION.md` beside this crate, and the ones that show
//! up in the code are marked `FRICTION:` where they bite.

#![deny(missing_docs)]

use dualis::prelude::*;
use serde::{Deserialize, Serialize};

// FRICTION 5: `Room` is not in the prelude, though `Tube` is and they are the same crate's
// two headline types. Reached for through the module instead. Aliased because this crate has
// a `DomainSpec::Room` variant of its own, which is the natural name on both sides.
use dualis::acoustic::Room as AcousticRoom;

pub mod render;

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
/// An enum rather than a registry of constructors, because `Simulation::with` takes a
/// concrete type. See `FRICTION.md`, finding 1.
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
    },
}

impl DomainSpec {
    /// The name this domain will answer to.
    pub fn name(&self) -> &str {
        match self {
            DomainSpec::Room { name, .. } | DomainSpec::Bar { name, .. } => name,
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
            // FRICTION 1: this has to be a match with one arm per domain type, and the arms
            // cannot be factored into a `Box<dyn Domain>` because `Simulation::with` takes
            // `impl Domain + 'static` by value. A scene format that a third party could
            // extend would need the kernel to accept a boxed domain.
            //
            // FRICTION 2: `leak`. Every domain constructor takes `name: &'static str`, and a
            // name read out of a file is a `String`. Leaking is bounded by the number of
            // domains in a scene, so it is survivable, but it is the API telling a consumer
            // that names are compile-time things when for an application they are data.
            sim = match spec {
                DomainSpec::Room {
                    name,
                    width_m,
                    height_m,
                    cells_across,
                    mode,
                    amplitude_pa,
                } => sim.with(
                    AcousticRoom::of_air(
                        leak(name),
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
                DomainSpec::Bar {
                    name,
                    length_mm,
                    cells,
                    area_mm2,
                    initial_c,
                } => sim.with(Bar1D::new(
                    leak(name),
                    Substance::aluminium_6061(),
                    *cells,
                    Length::mm(length_mm / *cells as f64),
                    Area::from_si(area_mm2 * 1e-6),
                    Temperature::celsius(*initial_c),
                )),
            };
        }

        Ok(World { scene, sim })
    }

    /// The scene this was built from.
    pub fn scene(&self) -> &Scene {
        &self.scene
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
    fn capture(&self) -> Frame {
        let t = self.sim.time();
        let mut panels = Vec::new();
        for spec in &self.scene.domains {
            // FRICTION 3: reading a domain's state back needs `domain_as::<T>`, which needs
            // the concrete type — so the renderer knows every domain type too, exactly like
            // `build` does. `ScalarField` is the abstraction that would avoid this, but a
            // `&dyn ScalarField` cannot be obtained from a `&dyn Domain`.
            match spec {
                DomainSpec::Room { name, .. } => {
                    if let Some(room) = self.sim.domain_as::<AcousticRoom>(name) {
                        panels.push(sample_room(name, room, t));
                    }
                }
                DomainSpec::Bar { name, .. } => {
                    if let Some(bar) = self.sim.domain_as::<Bar1D>(name) {
                        panels.push(sample_bar(name, bar));
                    }
                }
            }
        }
        Frame {
            time_s: t.to_si(),
            panels,
        }
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

fn sample_room(name: &str, room: &AcousticRoom, t: Time) -> Panel {
    let (nx, ny) = room.cells();
    let (w, h) = (room.width().to_si(), room.height().to_si());
    let probe = Length::from_si(w / (nx - 1) as f64);
    let mut values = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let p = LengthVec::m(
                w * i as f64 / (nx - 1) as f64,
                h * j as f64 / (ny - 1) as f64,
                0.0,
            );
            let _ = probe;
            values.push(room.at(p, t));
        }
    }
    Panel {
        name: name.to_string(),
        nx,
        ny,
        values,
        unit: "Pa",
    }
}

fn sample_bar(name: &str, bar: &Bar1D) -> Panel {
    let n = bar.cell_count();
    let values = (0..n)
        .map(|i| bar.temperature_at(i).to_si() - 273.15)
        .collect();
    Panel {
        name: name.to_string(),
        nx: n,
        ny: 1,
        values,
        unit: "C",
    }
}

/// Turn a scene's `String` name into the `&'static str` every constructor demands.
///
/// See `FRICTION.md`, finding 2. Bounded by the number of domains in a scene, which is why
/// this is survivable rather than merely wrong, but it is a leak and it is the API's fault.
fn leak(name: &str) -> &'static str {
    Box::leak(name.to_string().into_boxed_str())
}
