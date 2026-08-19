//! A world described as data, run against the dualis SDK, and drawn.
//!
//! This is the workspace's first consumer, and its first job is not to be a good application.
//! It is to be an *outside* user of the library — one that reaches for the API the way a
//! stranger would rather than the way its author remembers it — and to write down every place
//! that turns out to be awkward. A library with no consumers is a library whose ergonomics
//! nobody has measured.
//!
//! Findings are collected in `FRICTION.md` beside this crate. Eighteen of the twenty-four are
//! fixed — this crate is the record of what the API was like before, and the reason it changed.
//! The count in that file is under test; this line is not, and had been stale for two releases.
//!
//! The layers it once carried are libraries now. `dualis-scene` owns capture, and this crate is
//! left with what an *application* actually is: a file format, the domain types that format names,
//! and one place saying how far each field extends. A consumer that wants to draw a run no longer
//! has to reach into a binary that is `publish = false` to do it.

#![deny(missing_docs)]

use dualis::core::mixture::Mix;
use dualis::prelude::*;
use dualis::units::ThermalConductivity;
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
pub mod fit;
pub mod heater;
pub mod light;
pub mod verify;

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
    /// Which revision of this format the file is written in.
    ///
    /// **Absent means 1**, which is what every scene written before this field existed is. That
    /// default is the whole design: a reader can tell an old file from a new one, and does not
    /// have to guess from which keys happen to be present.
    ///
    /// A version this build does not know is [refused](Scene::check_version) rather than read
    /// hopefully. `serde`'s `deny_unknown_fields` already catches a key that was added; it cannot
    /// catch a key whose *meaning* changed, and that is what a version number is for.
    ///
    /// The promise attached to it is narrow and worth stating exactly: **within one version, a
    /// file that loads today loads tomorrow.** A change that would alter what an existing file
    /// means bumps this. A change that only adds an optional key does not.
    #[serde(default = "default_format")]
    pub format: u32,
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
    /// Per-quantity overrides, by channel name — `"energy"`, `"momentum"`, `"mass"`, `"charge"`,
    /// `"photons"`.
    ///
    /// The default above applies to everything not named here. A scene mixing schemes of
    /// different achievable accuracy needs this: a Barnes-Hut tree gives up exact momentum by
    /// construction while energy in a rigid room is exact to `1e-15`, and under one number either
    /// the momentum check refuses a correct run or the energy check stops seeing anything.
    ///
    /// A name this format does not know is **refused**, not ignored. A typo here would silently
    /// leave a quantity at the default — which is the shape of failure that turns an audit off,
    /// and this format has been bitten by exactly that once already with `aluminum`.
    #[serde(default)]
    pub tolerance_for: std::collections::BTreeMap<String, f64>,
    /// Substances this scene defines, by the name its domains will use.
    ///
    /// The catalogue holds nine materials and the world holds hundreds of thousands, so this is the
    /// door that matters: a `Substance` is data, so anything with a datasheet can be written here and
    /// used exactly as a catalogue entry is. Nothing downstream can tell the difference, and the
    /// [Stefan front of a declared gallium](https://docs.rs/dualis) is checked against Neumann's exact
    /// solution to the same standard as ice.
    ///
    /// ```json
    /// "materials": {
    ///   "gallium": {
    ///     "name": "gallium (99.99%)",
    ///     "density": 5904.0,
    ///     "thermal": { "conductivity": 40.6, "specific_heat": 371.0,
    ///                  "expansion": 1.8e-5, "emissivity": 0.1 },
    ///     "fusion": { "melting_point": 302.91, "latent_heat": 80160.0 }
    ///   }
    /// }
    /// ```
    ///
    /// A `BTreeMap` and not a `HashMap`, because the names come back out in error messages and a
    /// message whose word order changes between runs of the same file is not a message anybody can
    /// diff. Determinism is a promise this workspace makes about every byte it emits.
    ///
    /// Three things are refused rather than accepted quietly, and [`Palette`] says why each: an
    /// impossible substance, a name that shadows the catalogue, and — the one worth reading —
    /// a declaration that nothing goes on to use.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub materials: BTreeMap<String, Substance>,
    /// Composites this scene defines, by the name its domains will use.
    ///
    /// The other half of [`Scene::materials`], and the one that needed `Mix` to exist. A motor is copper,
    /// steel, magnets and air; a board is FR-4, copper and solder; a buffer is wax in a metal matrix. Each
    /// wants to be **one** material a coarse grid can hold, and its properties are not the properties of
    /// its main constituent.
    ///
    /// ```json
    /// "composites": {
    ///   "buffer": {
    ///     "parts": [
    ///       { "material": "octadecane", "volume_fraction": 0.8 },
    ///       { "material": "aluminium",  "volume_fraction": 0.2 }
    ///     ],
    ///     "conductivity_w_per_m_k": 12.5,
    ///     "emissivity": 0.9
    ///   }
    /// }
    /// ```
    ///
    /// **`volume_fraction`, spelled out, because volume and mass are the trap.** Volumetric heat capacity
    /// is volume-additive and specific heat is mass-weighted, and confusing them is worth 46% on a copper
    /// and FR-4 board — see [`Mix::specific_heat`](dualis::core::mixture::Mix::specific_heat). Wax filling
    /// 80% of a volume is 54.7% of the mass, so its latent heat dilutes to 133 kJ/kg and not 195.
    ///
    /// **The conductivity is the caller's and is checked.** No single value exists for a composite without
    /// knowing its microstructure, so this format cannot compute one — but it can refuse an impossible
    /// one, and it does: a value outside the Voigt and Reuss bounds is rejected with both bounds and the
    /// tighter Hashin–Shtrikman pair in the message, because those are the numbers somebody choosing needs
    /// and no scene file has anywhere else to get them.
    ///
    /// **The emissivity is the caller's and cannot be checked**, because it is a property of the surface
    /// and a mixture has no surface. A half-copper board is 0.05 bare and 0.9 under solder mask.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub composites: BTreeMap<String, CompositeSpec>,
    /// Where each domain sits, by the name it answers to. See [`PoseSpec`].
    ///
    /// A name-keyed map rather than a field on every [`DomainSpec`] variant, and the reason is
    /// not style: serde's `flatten` — the only way to add one optional field to fifteen
    /// variants without writing it fifteen times — **silently disables
    /// `deny_unknown_fields`**, which is this format's whole defence against a typo running
    /// something other than what the file says. A map beside `materials` and `composites`
    /// keeps every type strict and reads the way those do.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub poses: BTreeMap<String, PoseSpec>,
    /// The stage the experiment stands on. See [`EnvironmentSpec`].
    ///
    /// **Absent means what every scene written before this key existed means**: each domain's
    /// own assumption, unchanged — a bounce falls at standard gravity, an orbit falls only
    /// toward its own bodies. That default is what lets this key be added without a format
    /// bump: no existing file changes meaning.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub environment: Option<EnvironmentSpec>,
}

/// Where one domain sits in the world: a translation, and optionally a turn.
///
/// # What placing a part does today, and what it does not
///
/// It moves where the domain is **captured and drawn** — its field's samples, its bodies, its
/// box in the viewport — and it is applied through [`Pose`], the kernel's rigid motion, so a
/// rotation preserves every distance and angle exactly. What it does **not** do is make two
/// placed parts interact: domains meet on the bus and over an [`Interface`], and an interface
/// is matched by face count rather than by geometry, so two blocks placed against each other
/// are two blocks that are drawn against each other. `ARCHITECTURE.md` names that gap — "two
/// parts have no way to touch" — as the one that arrives first in practice, and placing them
/// is the half of it that can be done without answering the other half.
///
/// Saying so here rather than letting somebody discover it: a scene that places two parts in
/// contact and runs happily, with no heat crossing between them, is the shape of failure this
/// format exists to make impossible.
///
/// ```json
/// "poses": {
///   "armature": { "at_m": [0.12, 0.0, 0.0] },
///   "housing":  { "at_m": [0.0, 0.0, 0.05],
///                 "turn": { "axis": [0.0, 0.0, 1.0], "degrees": 45.0 } }
/// }
/// ```
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PoseSpec {
    /// Where the domain's origin lands in the world, in metres.
    #[serde(default)]
    pub at_m: [f64; 3],
    /// A turn about an axis through that origin. Absent is no rotation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub turn: Option<TurnSpec>,
}

/// A rotation, as an axis and an angle.
///
/// Axis-and-angle rather than a quaternion or three Euler angles, because it is the form a
/// person can write down and check: a quaternion's four numbers have a constraint between them
/// that a file cannot express, and Euler angles need a convention that every writer remembers
/// differently. An axis is a direction and an angle is an angle.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TurnSpec {
    /// The axis, in the domain's own coordinates. Normalised on the way in; a zero vector is
    /// refused rather than normalised into a NaN.
    pub axis: [f64; 3],
    /// How far to turn about it, in degrees. Degrees because a file is written by a person.
    pub degrees: f64,
}

impl PoseSpec {
    /// The kernel's [`Pose`] this describes, or why it is not one.
    pub fn to_pose(&self, site: &str) -> Result<Pose, String> {
        if !self.at_m.iter().all(|v| v.is_finite()) {
            return Err(format!(
                "poses.{site}: at_m is {:?}, and a position has to be a number",
                self.at_m
            ));
        }
        let translation = LengthVec::m(self.at_m[0], self.at_m[1], self.at_m[2]);
        let Some(turn) = &self.turn else {
            return Ok(Pose::at(translation));
        };
        if !turn.axis.iter().all(|v| v.is_finite()) || !turn.degrees.is_finite() {
            return Err(format!(
                "poses.{site}: the turn is not a number — axis {:?}, {} degrees",
                turn.axis, turn.degrees
            ));
        }
        let axis = glam::DVec3::new(turn.axis[0], turn.axis[1], turn.axis[2]);
        if axis.length() < 1e-12 {
            return Err(format!(
                "poses.{site}: the turn's axis is {:?}, which has no direction — a rotation needs one, and normalising a zero vector gives a NaN rather than an error",
                turn.axis
            ));
        }
        Ok(Pose::new(
            translation,
            glam::DQuat::from_axis_angle(axis.normalize(), turn.degrees.to_radians()),
        ))
    }
}

/// The stage an experiment stands on — the conditions that are the *scene's* and not any one
/// domain's.
///
/// This block exists because the stage used to be an assumption. Gravity was hardcoded into
/// the bounce's constructor, so a scene on the Moon could not be written, a scene in free
/// fall could not be written, and nothing in any file said `9.80665` anywhere a reader could
/// see it. A layer boundary turns assumptions into statements; this is that boundary for the
/// conditions an experiment is performed under.
///
/// # Every domain must answer a stated environment, and there are three answers
///
/// - **Consume it.** A bounce under stated gravity falls at that gravity — the Moon at 1.62,
///   a drop tower at 0.
/// - **Refuse it.** Bodies under their own mutual gravity cannot also stand in a stated
///   uniform field — the domain has no way to consume one, and quietly not consuming it
///   would run a different experiment than the file describes. The build refuses, naming the
///   domain.
/// - **Dismiss it, with the measurement that earns the dismissal.** A molecular fluid under
///   stated gravity ignores it *correctly*: for argon, `m·g·σ` — the gravitational energy
///   across one molecular diameter — is `2.2e-34 J` against a well depth of `1.65e-21 J`, a
///   ratio of `1.3e-13`. The dismissal is reported as a build note (printed by `--check` and
///   shown in the editor) rather than left silent, because a stated condition that is
///   ignored without a word is this format's oldest failure shape.
///
/// The fourth answer — silence — is not available.
///
/// # What is deliberately not here yet
///
/// No ground plane: the bounce owns its own floor today, and a stage-level ground earns its
/// place when a *second* domain needs to stand on the same one — the rule `World::advance`
/// was made public under. No ambient temperature: every thermal domain already states its own
/// `ambient_c` explicitly, which is more honest than a default it might not mean. Each moves
/// here when something needs them shared, not before.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnvironmentSpec {
    /// Uniform gravity, in m/s², pulling along −z. Zero is free fall; standard Earth is
    /// 9.80665, which is also what an *absent* environment means for the domains that fall.
    ///
    /// A magnitude, not a vector: −z is the format's down, and a scene that wants sideways
    /// gravity wants a rotated scene. Negative is refused rather than read as "up".
    pub gravity_m_per_s2: f64,
}

impl EnvironmentSpec {
    /// Refuse a stage that does not describe one.
    fn check(&self) -> Result<(), String> {
        if !self.gravity_m_per_s2.is_finite() || self.gravity_m_per_s2 < 0.0 {
            return Err(format!(
                "environment: gravity_m_per_s2 is {}; it is a magnitude along -z, so it must \
                 be finite and non-negative — a scene wanting gravity in another direction \
                 wants a rotated scene",
                self.gravity_m_per_s2
            ));
        }
        Ok(())
    }
}

/// One composite in [`Scene::composites`].
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CompositeSpec {
    /// What it is made of, and how much of the volume each part is.
    pub parts: Vec<CompositePart>,
    /// The conductivity to use, in W/m·K. Refused if no microstructure of these parts could have it.
    pub conductivity_w_per_m_k: f64,
    /// Emissivity of whatever surface it presents, 0 to 1.
    pub emissivity: f64,
}

/// One constituent of a [`CompositeSpec`].
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CompositePart {
    /// A catalogue name or one of [`Scene::materials`].
    ///
    /// **Not another composite.** A composite of composites needs a declaration order, and a `BTreeMap`
    /// has none by design — but the real reason is that nesting changes what the bounds mean: the
    /// Hashin–Shtrikman pair is a *two-phase* result and a mixture of three things has none, which
    /// `a_mixture.rs` records rather than papers over. Flatten it: list all the parts with their fractions
    /// of the whole.
    pub material: String,
    /// This part's share of the **volume**, not of the mass. The fractions must sum to one.
    pub volume_fraction: f64,
}

fn default_tolerance() -> f64 {
    1e-6
}

/// The revision this build writes, and the highest it can read.
///
/// One, still. It goes to two the first time a change would make an existing file mean something
/// different — not when a key is added, which an old file simply does not have.
pub const FORMAT: u32 = 1;

fn default_format() -> u32 {
    1
}

impl Scene {
    /// Where each domain sits, by name — the domain's own extent, moved by its [`PoseSpec`].
    ///
    /// **The one place poses are applied.** `World` uses it to capture and the editors use it
    /// to draw; a shell that composed the pose itself would be a second implementation of the
    /// same composition, and the two would drift the first time one of them learned something.
    ///
    /// A malformed pose is skipped here rather than reported, because this runs on a scene that
    /// may not have been built yet — `World::build` is where a bad pose is refused, and it is
    /// called before anything is run.
    pub fn placements(&self) -> BTreeMap<String, Placement> {
        self.domains
            .iter()
            .map(|spec| {
                let mut placement = spec.placement();
                if let Some(pose) = self.poses.get(spec.name()) {
                    if let Ok(p) = pose.to_pose(spec.name()) {
                        placement.pose = p;
                    }
                }
                (spec.name().to_string(), placement)
            })
            .collect()
    }

    /// Refuse a file from a future this build does not know.
    ///
    /// Called by [`World::build`], so nothing can run a scene it half-understands. The failure
    /// being prevented is the one this format has already had once in a smaller form: a key that
    /// is not read leaves a field at its default and the run proceeds, quietly doing something
    /// other than what the file says.
    ///
    /// Refusing forward rather than attempting a downgrade is deliberate. A newer file may use a
    /// key this build does not have, or the same key differently; reading it on a best-effort
    /// basis produces a run that is plausible and not the one that was written down.
    pub fn check_version(&self) -> Result<(), String> {
        if self.format == 0 {
            return Err(format!(
                "format 0 is not a version this format ever had; scenes are {FORMAT} or, if the \
                 key is absent, 1"
            ));
        }
        if self.format > FORMAT {
            return Err(format!(
                "this scene is format {} and this build reads up to {FORMAT}. It was written by a \
                 newer dualis; upgrade rather than run it, because a key this build does not know \
                 would be left at a default and the run would not be the one in the file.",
                self.format
            ));
        }
        Ok(())
    }
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
    ///
    /// The samples sit on the walls, so `cells_across` samples are `cells_across − 1`
    /// intervals, and the **height is quantised** to whole intervals of the spacing the width
    /// sets — exactly as [`DomainSpec::Hall`] documents for three dimensions. A 3.1 m height
    /// at 61 samples across 4.4 m is really 3.08 m, and `verify`'s resolution sweep measured
    /// the consequence before this sentence existed: the quantised height converges on the
    /// stated one at first order, and every mode frequency rides along with it.
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
    /// A block of conductor with two electrodes, **solved** for its potential.
    ///
    /// `dualis-electrical`'s `Conductor`. Every other source in this format states its watts and
    /// `winding` computes them from `I²R`; this one computes them from a *shape*. Nobody states a
    /// resistance — it comes out of the solve, and for a uniform block it comes out as `ρL/A`
    /// exactly, which is what makes a notch or a via checkable at all.
    Conductor {
        /// Domain name.
        name: String,
        /// Cells along x, y and z. Current runs along x, between the two end faces.
        cells: [usize; 3],
        /// The side of one cubic cell.
        cell_mm: f64,
        /// Resistivity of the bulk material, in ohm-metres. Copper is 1.724e-8.
        resistivity_ohm_m: f64,
        /// The potential difference across the end faces.
        volts: f64,
        /// Cells to make insulating, as `[i, j, k]`. A notch, a via wall, a crack.
        ///
        /// This is the whole reason to solve rather than to state: a block with a notch has a
        /// resistance that is a property of its geometry and no formula gives it.
        #[serde(default)]
        blocked: Vec<[usize; 3]>,
    },
    /// A basket of packed grounds with liquid driven through it.
    ///
    /// `dualis-porous`'s `Puck`. The flow is not stated anywhere: Darcy's law is solved on the
    /// permeability the grind and the packing give, and everything else — how long the shot takes,
    /// where the liquid goes, what comes out — follows from it.
    ///
    /// `channel_porosity` is how a fault is put into a scene. A basket is otherwise uniform,
    /// because an even tamp gives an even bed and channelling is a *defect* rather than something
    /// the physics produces on its own.
    Puck {
        /// Domain name.
        name: String,
        /// Cells along x, y and z. **`y` is the flow axis**, so the middle one is the bed's depth.
        cells: [usize; 3],
        /// The side of one cubic cell.
        cell_mm: f64,
        /// The basket's inside radius. Leave a few cells of box outside it for the metal.
        radius_mm: f64,
        /// Sieve diameter of the grind. 250 µm is a conventional espresso setting.
        grind_um: f64,
        /// Inter-particle void fraction after tamping. About 0.45.
        porosity: f64,
        /// The pressure held across the bed.
        bar: f64,
        /// Brew temperature, and the basket's own starting temperature unless `wall_c` says
        /// otherwise.
        brew_c: f64,
        /// The basket's starting temperature, for the cold-portafilter case.
        #[serde(default)]
        wall_c: Option<f64>,
        /// Porosity of the ring against the basket wall — a puck that shrank away from it.
        ///
        /// `None` is an even bed. Give a number above `porosity` and the flow takes the ring.
        #[serde(default)]
        channel_porosity: Option<f64>,
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
        /// What it is made of, one of [`MATERIALS`]. Absent is `aluminium`, which is what every
        /// block written before this key existed is.
        #[serde(default)]
        material: Option<String>,
        /// Boxes of cells made of something else — a layer, a coating, a joint, an inclusion.
        ///
        /// Applied in order, so a later region overwrites an earlier one where they overlap. That
        /// is worth stating because it is how a coating *on* a layer is written and there is no
        /// other way to mean it.
        #[serde(default)]
        regions: Vec<Region>,
        /// A cell to warm at the start, and by how much, so there is a hot spot to watch spread.
        ///
        /// **A statement about the initial state, not a delivery of heat.** It moves what the
        /// block holds and not what it has absorbed, so the audit's opening balance includes it
        /// — which is the honest bookkeeping and is why it is separate from any source.
        #[serde(default)]
        hot_spot: Option<HotSpot>,
        /// Designed parts to fill this block from, each an STL and a material. See [`PartSpec`].
        ///
        /// **Stating any part makes every cell outside them void** — nothing, rather than the
        /// block's bulk material. That is what an assembly in air *is*, and the alternative is
        /// what this format did before there was a void at all: two parts conducting through
        /// the gap between them as though the gap were metal.
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        parts: Vec<PartSpec>,
        /// Faces that lose heat to something, and to what. See [`CoolingSpec`].
        ///
        /// **Absent is a block insulated on all six faces**, which is what every scene written
        /// before this key existed is — and, until the domain could shed heat at all, what every
        /// three-dimensional thermal scene in this repository was. Such a scene warms for as long
        /// as it runs and never settles anywhere, which answers a pulse and not the question a
        /// designer asks.
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        cooling: Vec<CoolingSpec>,
        /// Boxes that **generate** heat, and how many watts. See [`DissipationSpec`].
        ///
        /// The key without which no real application can be written. Every other source in this
        /// format hands watts to the bus, and the bus carries an amount and no location by
        /// design — heat arriving there spreads to a uniform rise over everything that can hold
        /// it, which is the only choice that adds no information and the wrong answer for a die,
        /// a winding, a brake disc or a laser absorber. All of them dissipate *somewhere*, and
        /// the gradient between there and the heatsink is the question.
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        dissipation: Vec<DissipationSpec>,
    },
    /// A solid body under load, solved for its displacement — `dualis-elastic`'s `Block`.
    ///
    /// The first domain in this format whose unknown is a **vector**, and the first that answers
    /// what a shape *does* rather than what it holds. Elliptic: there is no time in it, so every
    /// frame is the answer for that frame's loads and the run is a sequence of static solves.
    ///
    /// # What makes it worth having beside a thermal block
    ///
    /// `follows` names a `block` whose temperature field becomes this body's **stress-free
    /// strain**, `α(T − T_ref)` element by element. That is the coupling a digital twin is for:
    /// the platform could already compute a power module's temperature field to four figures and
    /// could do nothing with it, and the failure mode of a real module is not its temperature but
    /// the solder fatigue that comes of silicon at 2.6e-6 per kelvin sitting on solder at 2.15e-5.
    ///
    /// The two grids must match, and a mismatch is refused rather than interpolated: an element
    /// and a cell are the same box or the coupling is a guess about which cell a corner belongs to.
    Structure {
        /// Domain name.
        name: String,
        /// Elements along x, y and z.
        cells: [usize; 3],
        /// The side of one cubic element.
        cell_mm: f64,
        /// What it is made of, one of [`MATERIALS`]. Absent is `aluminium`.
        #[serde(default)]
        material: Option<String>,
        /// Boxes of elements made of something else — a layer, a joint, an inclusion. Same rules
        /// as a block's [`Region`], and `initial_c` is refused because a structure has no
        /// temperature of its own.
        #[serde(default)]
        regions: Vec<Region>,
        /// Faces that are held, and how. See [`HeldSpec`].
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        held: Vec<HeldSpec>,
        /// Faces under pressure. See [`PressedSpec`].
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        pressed: Vec<PressedSpec>,
        /// The name of a `block` whose temperature drives this body's expansion.
        ///
        /// Absent is a body with no thermal strain at all, which is every structure loaded only
        /// by what `pressed` states.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        follows: Option<String>,
        /// The temperature at which the body is the size it was drawn — its stress-free state.
        ///
        /// **Not the same as the block's starting temperature**, and conflating them is the
        /// commonest way to get a thermal-stress answer that is confidently wrong: a power module
        /// is assembled at its solder's reflow temperature and *starts* at whatever the room is,
        /// so it is already stressed before it is switched on. Required when `follows` is stated,
        /// because there is no sensible default for it.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        reference_c: Option<f64>,
    },
    /// Incompressible flow on a staggered grid — `dualis-fluid`'s `Channel`.
    ///
    /// The domain this workspace's own documentation calls the hardest to trust, and it says why:
    /// "it looks like a fluid" is the easiest wrong answer in computational physics to accept,
    /// because a scheme with the wrong viscosity still makes plausible vortices and one that
    /// quietly loses momentum still makes a pretty picture. So a scene using it should be written
    /// around one of the three exact solutions that exist — Poiseuille, Couette or Taylor–Green —
    /// and scene `26` is.
    ///
    /// Nothing to draw: velocity is a vector on a staggered grid and this format's field panel is
    /// a scalar at cell centres, so a channel reports its scalars and exports no geometry. That is
    /// the same shape as the other twelve fieldless scenes and is stated rather than discovered.
    Channel {
        /// Domain name.
        name: String,
        /// Cells along x, y and z. `y` is across the channel, which is where the walls are.
        cells: [usize; 3],
        /// The side of one cubic cell.
        cell_mm: f64,
        /// What is flowing. See [`FluidSpec`].
        fluid: FluidSpec,
        /// The walls, or absent for a periodic box. See [`WallsSpec`].
        #[serde(default, skip_serializing_if = "Option::is_none")]
        walls: Option<WallsSpec>,
        /// A body force per unit mass, m/s² — a pressure gradient written as what it does.
        ///
        /// `[g, 0, 0]` down a walled channel is Poiseuille flow, whose steady profile is
        /// `u(y) = (g/2ν)·y(h−y)` and whose mean is `g·h²/(12ν)`. Absent is no drive, which for a
        /// channel started from rest is a fluid that stays at rest.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        drive_m_per_s2: Option<[f64; 3]>,
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

/// Where a scene's `parts` find their bytes.
///
/// A scene names a part with a string, and on a machine with a filesystem that string is a path.
/// In a browser it is not: there is no filesystem, the page already **has** the bytes because
/// somebody dropped a file on it, and the string is a label. Both are the same question — *give me
/// the bytes called this* — so the scene format does not have to know which world it is in, and
/// the same file runs in both.
///
/// The default is [`OnDisk`], which is what [`World::build`] uses and what the CLI has always
/// done. [`Uploaded`] is the other one, and it is not a browser-only convenience: it is how a test
/// builds an assembly without writing a fixture to a temporary directory, which is how the tests
/// for this trait are written.
pub trait Parts {
    /// The bytes for `name`, or why not — the message reaches the user with the scene's own site
    /// prefix in front of it, so it should say what was looked for and where.
    fn bytes(&self, name: &str) -> Result<Vec<u8>, String>;
}

/// [`Parts`] backed by the filesystem: the name is a path.
#[derive(Debug, Clone, Copy, Default)]
pub struct OnDisk;

impl Parts for OnDisk {
    fn bytes(&self, name: &str) -> Result<Vec<u8>, String> {
        std::fs::read(name).map_err(|e| format!("{name}: {e}"))
    }
}

/// [`Parts`] held in memory: the name is a label somebody chose.
///
/// What a browser has after a drop, and what a test has instead of a temporary directory. The
/// refusal names what *is* here, because "no such part" with nothing else in it is the least
/// useful thing this can say to somebody who has just uploaded three files and misspelled one.
#[derive(Debug, Clone, Default)]
pub struct Uploaded {
    files: BTreeMap<String, Vec<u8>>,
}

impl Uploaded {
    /// An empty set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add or replace one file. Returns `self` so a caller can chain.
    pub fn with(mut self, name: impl Into<String>, bytes: impl Into<Vec<u8>>) -> Self {
        self.files.insert(name.into(), bytes.into());
        self
    }

    /// Add or replace one file in place.
    pub fn insert(&mut self, name: impl Into<String>, bytes: impl Into<Vec<u8>>) {
        self.files.insert(name.into(), bytes.into());
    }

    /// Forget one file, and say whether there was one.
    pub fn remove(&mut self, name: &str) -> bool {
        self.files.remove(name).is_some()
    }

    /// The names held, in order.
    pub fn names(&self) -> Vec<&str> {
        self.files.keys().map(String::as_str).collect()
    }

    /// How many bytes are held in total, which is the number a page showing an upload list wants.
    pub fn total_bytes(&self) -> usize {
        self.files.values().map(Vec::len).sum()
    }
}

impl Parts for Uploaded {
    fn bytes(&self, name: &str) -> Result<Vec<u8>, String> {
        self.files.get(name).cloned().ok_or_else(|| {
            if self.files.is_empty() {
                format!("{name}: no file by that name, and nothing has been uploaded")
            } else {
                format!(
                    "{name}: no file by that name; uploaded are {}",
                    self.names().join(", ")
                )
            }
        })
    }
}

/// The one material name this format reserves: a box or a part of **nothing**.
///
/// Not a substance with a very low conductivity, which is what a scene had to write before there
/// was a void at all — and which still conducts, still stores heat and still sets a stability
/// limit. Measured, on three identical bars differing only in one cell: the far end warms 50 K
/// through copper, still warms through the catalogue's poorest insulator, and moves only by what
/// radiation carries across nothing.
pub const VOID: &str = "void";

/// The materials a scene may name without declaring them: [`Substance::CATALOGUE`], verbatim.
///
/// It is an alias now and used to be a hand-written copy — eight names beside the catalogue's nine,
/// which is why `water` could not be named in any of eleven releases. A format's spelling of a
/// material belongs beside the material.
pub const MATERIALS: [&str; 9] = Substance::CATALOGUE;

/// The substances a scene can resolve a name to: the catalogue, plus whatever the scene declared.
///
/// # Why this is a type and not a function
///
/// It records what it was asked for. A declaration nothing asks for is [refused](Palette::unused),
/// and "asked for" has exactly one honest definition: it went through [`Palette::get`]. The
/// alternative is a second list of the places a material name can appear in this format, which would
/// be a list that agrees with the builder until somebody adds a field — the defect [`MATERIALS`] was
/// just cured of, reintroduced one level up.
pub struct Palette {
    declared: BTreeMap<String, Substance>,
    asked: std::collections::BTreeSet<String>,
}

impl Palette {
    /// Validate a scene's declarations and prepare to resolve names against them.
    ///
    /// Three refusals here, and each is a wrong answer this format would otherwise produce:
    ///
    /// - **An impossible substance.** [`Substance::check`], with the declared name as the site. A
    ///   negative conductivity or a zero density reaches a sweep as `NaN` and poisons a field with
    ///   nothing saying which number it came from.
    /// - **A name that shadows the catalogue.** Two files that both say `"copper"` have to mean the
    ///   same copper, or a run's material is a property of the file it was launched from and no
    ///   comparison between two runs means anything.
    /// - **A declaration nothing names.** See [`Palette::unused`] — this is the one that is not
    ///   obvious and it is the important one.
    pub fn new(declared: &BTreeMap<String, Substance>) -> Result<Palette, String> {
        Palette::with_composites(declared, &BTreeMap::new())
    }

    /// The same, plus the composites the scene declared.
    ///
    /// Two passes, and the order is the whole design: substances first, because a composite's parts are
    /// resolved against them, and composites second. A composite may not name another composite — see
    /// [`CompositePart::material`] — which is what keeps this two passes rather than a dependency graph.
    ///
    /// A part that resolves to a catalogue entry or a declared substance counts as **used**, so a
    /// substance declared only to be mixed is not reported as dead weight. That is the interaction between
    /// this and [`Palette::unused`] that a first draft got wrong: declaring a wax purely to put it in a
    /// composite was refused as unused, which is the opposite of the rule's purpose.
    pub fn with_composites(
        declared: &BTreeMap<String, Substance>,
        composites: &BTreeMap<String, CompositeSpec>,
    ) -> Result<Palette, String> {
        for (key, substance) in declared {
            if key.is_empty() {
                return Err("materials: a declared material needs a name".to_string());
            }
            if key == VOID {
                return Err(format!(
                    "materials.{key}: {VOID:?} is this format's word for **nothing** and a scene \
                     may not redefine it as a substance — a region or a part naming it would \
                     then be solid in one file and empty in another"
                ));
            }
            if Substance::from_name(key).is_some() {
                return Err(format!(
                    "materials.{key}: {key:?} is already a catalogue material and a scene may not \
                     redefine it — two files saying {key:?} have to mean the same substance, or no \
                     comparison between two runs means anything. Pick another name",
                ));
            }
            substance
                .check()
                .map_err(|e| format!("materials.{key}: {e}"))?;
        }
        let mut palette = Palette {
            declared: declared.clone(),
            asked: std::collections::BTreeSet::new(),
        };

        for (key, spec) in composites {
            if key.is_empty() {
                return Err("composites: a declared composite needs a name".to_string());
            }
            if Substance::from_name(key).is_some() {
                return Err(format!(
                    "composites.{key}: {key:?} is already a catalogue material and a scene may not \
                     redefine it"
                ));
            }
            if declared.contains_key(key) {
                return Err(format!(
                    "composites.{key}: {key:?} is also declared as a material, and one name cannot mean \
                     two things"
                ));
            }
            if composites.contains_key(key)
                && spec
                    .parts
                    .iter()
                    .any(|p| composites.contains_key(&p.material))
            {
                return Err(format!(
                    "composites.{key}: a part names another composite, which this format does not \
                     support — flatten it and list every part's fraction of the whole. See the docs on \
                     `CompositePart::material` for why nesting is not just an ordering problem"
                ));
            }
            let mut parts = Vec::with_capacity(spec.parts.len());
            for (n, part) in spec.parts.iter().enumerate() {
                let site = format!("composites.{key}.parts[{n}]");
                let substance = palette.get(&site, &part.material)?;
                parts.push((substance, part.volume_fraction));
            }
            let mix = Mix::of(&parts).map_err(|e| format!("composites.{key}: {e}"))?;
            let conductivity = ThermalConductivity::w_per_m_k(spec.conductivity_w_per_m_k);
            let substance = mix
                .as_substance(key, conductivity, spec.emissivity)
                .map_err(|e| {
                    // The refusal `Mix` gives names the outer bounds. A scene has nowhere else to learn
                    // the tighter pair, and that is usually the number somebody wants, so it is added
                    // here rather than left for them to compute.
                    match mix.hashin_shtrikman() {
                        Some((lo, hi)) => format!(
                            "composites.{e}. For an isotropic microstructure the Hashin–Shtrikman pair \
                             is {} to {} W/m/K, which is the narrower range and usually the one to \
                             choose inside",
                            lo.to_si(),
                            hi.to_si()
                        ),
                        None => format!("composites.{e}"),
                    }
                })?;
            palette.declared.insert(key.clone(), substance);
        }
        // The parts asked for above stay recorded, which is what makes a substance mixed into a composite
        // count as used. The composites themselves are *inserted* rather than asked for, so each is still
        // subject to the unused check until a domain names it — the behaviour wanted, and it falls out of
        // the distinction between `get` and `declared` rather than needing a rule.
        //
        // A first draft cleared this set here, reasoning that naming a part is not using the composite.
        // True and irrelevant: it also erased the parts, so a wax declared purely to be mixed was refused
        // as dead weight — the rule firing on the case it exists to protect. Three tests failed on one line.
        Ok(palette)
    }

    /// Resolve the name a scene wrote, recording that it was asked for.
    ///
    /// The catalogue first, then the declarations — an order that cannot matter, because
    /// [`Palette::new`] has already refused any declaration that could collide.
    ///
    /// The error lists every name that would have worked, because a caller who guessed wrong has no
    /// other way to find out: this is a JSON file with no completion and no type checker behind it.
    pub fn get(&mut self, site: &str, material: &str) -> Result<Substance, String> {
        self.asked.insert(material.to_string());
        if let Some(s) = Substance::from_name(material) {
            return Ok(s);
        }
        if let Some(s) = self.declared.get(material) {
            return Ok(s.clone());
        }
        let mut known: Vec<String> = MATERIALS.iter().map(|m| format!("{m:?}")).collect();
        known.extend(self.declared.keys().map(|k| format!("{k:?} (declared)")));
        Err(format!(
            "{site}: unknown material {material:?}; known materials are {}",
            known.join(", ")
        ))
    }

    /// Every declared name that no domain asked for.
    ///
    /// # Why this is an error and not a warning
    ///
    /// Because `material` on a block is **optional and defaults to aluminium**. A scene that declares
    /// a substance and then does not name it — a field left off, a name misspelled at the *use* site
    /// where the declaration is spelled right — runs as a block of aluminium. It runs, it audits, it
    /// renders, and it answers a question about the wrong material with nothing anywhere saying so.
    ///
    /// That is the same failure a region selecting no cells has, one level up, and it gets the same
    /// treatment for the same reason. The cost of being wrong about this is a line deleted from a
    /// file; the cost of the other way is a number somebody believes.
    pub fn unused(&self) -> Vec<&str> {
        self.declared
            .keys()
            .filter(|k| !self.asked.contains(*k))
            .map(|k| k.as_str())
            .collect()
    }
}

/// One body in a [`DomainSpec::Network`].
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NetworkNode {
    /// What the links call it. Unique within the network.
    pub name: String,
    /// One of [`MATERIALS`]. Mixing them is the point: a motor is not a billet of its main metal.
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
            | DomainSpec::Structure { name, .. }
            | DomainSpec::Channel { name, .. }
            | DomainSpec::Hall { name, .. }
            | DomainSpec::Conductor { name, .. }
            | DomainSpec::Puck { name, .. }
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
    ///
    /// Takes the [`Palette`] by mutable reference because resolving a material name is what *records*
    /// that the name was used, and a declared substance nobody used is refused. See
    /// [`Palette::unused`].
    ///
    /// Takes the scene's [`EnvironmentSpec`] because the stage is stated once and consumed
    /// here, where the domains are constructed — `None` is a scene written before the stage
    /// could be stated, and means each domain's own assumption.
    pub fn build(
        &self,
        palette: &mut Palette,
        environment: Option<&EnvironmentSpec>,
        files: &dyn Parts,
        notes: &mut Vec<String>,
    ) -> Result<Box<dyn Domain>, String> {
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
                // Refused, not ignored: these bodies fall only toward each other, and the
                // domain has no way to consume a uniform field. Building it anyway would run
                // a different experiment than the file states — in free fall while the file
                // says Earth — which is the failure this format exists to make impossible.
                if let Some(env) = environment {
                    if env.gravity_m_per_s2 != 0.0 {
                        return Err(format!(
                            "{name}: an orbit cannot stand in the stated uniform gravity of \
                             {} m/s^2 — its bodies fall only toward each other. State \
                             \"gravity_m_per_s2\": 0.0 for a free-fall stage, or remove the \
                             environment and let each domain keep its own assumption",
                            env.gravity_m_per_s2
                        ));
                    }
                }
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
            } => {
                // The stage's gravity when one is stated, the old hardcoded assumption when
                // not — which keeps every scene written before the stage existed meaning
                // exactly what it meant. Zero is a legitimate stage: a drop tower's ball
                // floats on its dashpot.
                let g = environment.map_or(G0.to_si(), |e| e.gravity_m_per_s2);
                Box::new(ContactSystem::new(
                    name.clone(),
                    &[Body::new(
                        Mass::kg(*mass_kg),
                        LengthVec::m(0.0, 0.0, *drop_m),
                        VelocityVec::ZERO,
                    )],
                    AccelerationVec::from_si(-glam::DVec3::Z * g),
                    Ground::floor(),
                    Stiffness::from_si(*stiffness),
                    Damping::from_si(*damping),
                ))
            }
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
                                "{name}: unknown finish {other:?}; known finishes are \"aluminium\""
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
                    let substance =
                        palette.get(&format!("{name}/{}", n.name), n.material.as_str())?;
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
                        None => network.node(n.name.clone(), substance, volume, thickness, initial),
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
            DomainSpec::Conductor {
                name,
                cells,
                cell_mm,
                resistivity_ohm_m,
                volts,
                blocked,
            } => {
                let mut c = dualis::electrical::Conductor::new(
                    name.clone(),
                    (cells[0], cells[1], cells[2]),
                    Length::mm(*cell_mm),
                    dualis::units::Resistivity::ohm_m(*resistivity_ohm_m),
                    dualis::units::Voltage::v(*volts),
                );
                // A practical insulator rather than a literal zero: a conductance of exactly zero
                // makes the system singular for any cell it isolates, and the solver would be
                // right to refuse. Twelve orders of magnitude is a crack, not a wire.
                let insulator = dualis::units::Resistivity::ohm_m(resistivity_ohm_m * 1e12);
                for at in blocked {
                    c.set_resistivity(at[0], at[1], at[2], insulator);
                }
                // Re-solved here rather than left to the first step, because `blocked` changed
                // the problem after the constructor solved the unnotched one — and the first
                // frame is captured before anything steps.
                c.solve(1e-12);
                Box::new(c)
            }
            DomainSpec::Puck {
                name,
                cells,
                cell_mm,
                radius_mm,
                grind_um,
                porosity,
                bar,
                brew_c,
                wall_c,
                channel_porosity,
            } => {
                let mut p = dualis::porous::Puck::new(
                    name.clone(),
                    dualis::porous::Basket {
                        counts: (cells[0], cells[1], cells[2]),
                        cell: Length::mm(*cell_mm),
                        radius: Length::mm(*radius_mm),
                        grind: dualis::porous::Grind::sieved(Length::from_si(grind_um * 1e-6)),
                        porosity: *porosity,
                        pressure: dualis::units::Pressure::from_si(bar * 1e5),
                        temperature: dualis::units::Temperature::celsius(*brew_c),
                        ..dualis::porous::Basket::espresso()
                    },
                );
                if let Some(loose) = channel_porosity {
                    // The ring of packed cells against the wall, computed before the repack so
                    // that widening one cell does not change what counts as the ring.
                    let (nx, ny, nz) = (cells[0], cells[1], cells[2]);
                    let mut ring = vec![false; nx * ny * nz];
                    for k in 0..nz {
                        for j in 0..ny {
                            for i in 0..nx {
                                if !p.is_packed(i, j, k) {
                                    continue;
                                }
                                ring[i + nx * (j + ny * k)] = i == 0
                                    || k == 0
                                    || i + 1 == nx
                                    || k + 1 == nz
                                    || !p.is_packed(i - 1, j, k)
                                    || !p.is_packed(i + 1, j, k)
                                    || !p.is_packed(i, j, k - 1)
                                    || !p.is_packed(i, j, k + 1);
                            }
                        }
                    }
                    p.repack(*loose, |i, j, k| ring[i + nx * (j + ny * k)]);
                }
                if let Some(c) = wall_c {
                    p.set_wall_temperature(dualis::units::Temperature::celsius(*c));
                    p.set_inlet_temperature(dualis::units::Temperature::celsius(*brew_c));
                }
                Box::new(p)
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
            DomainSpec::Channel {
                name,
                cells,
                cell_mm,
                fluid,
                walls,
                drive_m_per_s2,
            } => {
                if cells.contains(&0) {
                    return Err(format!("{name}: a channel needs at least one cell"));
                }
                if cell_mm.is_nan() || *cell_mm <= 0.0 {
                    return Err(format!("{name}: cell_mm must be positive, got {cell_mm}"));
                }
                let stuff = match fluid {
                    FluidSpec::Named(n) => match n.as_str() {
                        "water" => dualis::fluid::Fluid::water(),
                        "air" => dualis::fluid::Fluid::air(),
                        other => {
                            return Err(format!(
                                "{name}: {other:?} is not a fluid this format knows. It knows \
                                 \"water\" and \"air\", or state `density` and \
                                 `kinematic_viscosity` outright"
                            ));
                        }
                    },
                    FluidSpec::Stated {
                        density,
                        kinematic_viscosity,
                    } => {
                        if !(density.is_finite() && *density > 0.0) {
                            return Err(format!("{name}: density must be positive, got {density}"));
                        }
                        if !(kinematic_viscosity.is_finite() && *kinematic_viscosity > 0.0) {
                            return Err(format!(
                                "{name}: kinematic_viscosity must be positive, got \
                                 {kinematic_viscosity}. It is `mu/rho` in m^2/s, which for water \
                                 is about 1e-6 — a thousand times smaller than the dynamic \
                                 viscosity tables quote beside it"
                            ));
                        }
                        dualis::fluid::Fluid {
                            density: dualis::units::Density::from_si(*density),
                            kinematic_viscosity: dualis::units::Diffusivity::from_si(
                                *kinematic_viscosity,
                            ),
                        }
                    }
                };
                let boundary = match walls {
                    Some(w) => dualis::fluid::Walls::Sliding {
                        low: w.lower_m_per_s,
                        high: w.upper_m_per_s,
                    },
                    None => dualis::fluid::Walls::None,
                };
                let mut channel = dualis::fluid::Channel::new(
                    name.clone(),
                    (cells[0], cells[1], cells[2]),
                    Length::mm(*cell_mm),
                    stuff,
                    boundary,
                );
                if let Some(g) = drive_m_per_s2 {
                    if g.iter().any(|c| !c.is_finite()) {
                        return Err(format!("{name}: {g:?} is not a body force"));
                    }
                    channel.drive(glam::DVec3::new(g[0], g[1], g[2]));
                }
                Box::new(channel)
            }
            DomainSpec::Structure {
                name,
                cells,
                cell_mm,
                material,
                regions,
                held,
                pressed,
                follows,
                reference_c,
            } => {
                if cells.contains(&0) {
                    return Err(format!("{name}: a structure needs at least one element"));
                }
                // NaN spelled out rather than folded into a negated comparison, which clippy is
                // right to object to and which this format has been caught by before.
                if cell_mm.is_nan() || *cell_mm <= 0.0 {
                    return Err(format!("{name}: cell_mm must be positive, got {cell_mm}"));
                }
                let base = palette.get(name, material.as_deref().unwrap_or("aluminium"))?;
                let elastic_of = |substance: &Substance, site: &str| {
                    substance
                        .mechanical
                        .map(|m| dualis::elastic::Elastic {
                            youngs_modulus: m.youngs_modulus,
                            poisson_ratio: m.poisson_ratio,
                            density: substance.density,
                        })
                        .ok_or_else(|| {
                            format!(
                        "{site}: {} has no mechanical properties, so it cannot be a structure — \
                         state `youngs_modulus`, `poisson_ratio` and `yield_strength` in its \
                         `mechanical` block",
                        substance.name
                    )
                        })
                };

                let counts = (cells[0], cells[1], cells[2]);
                let mut body = dualis::elastic::Block::new(
                    name.clone(),
                    counts,
                    Length::mm(*cell_mm),
                    elastic_of(&base, name)?,
                );
                // Expansion coefficients per element, kept beside the body because `Elastic` has
                // no place for one — an elastic domain that carried a thermal coefficient would be
                // a domain that knew about temperature.
                let mut alpha = vec![
                    base.thermal.map_or(0.0, |t| t.expansion.to_si());
                    counts.0 * counts.1 * counts.2
                ];

                for (n, r) in regions.iter().enumerate() {
                    let site = format!("{name}/regions[{n}]");
                    if r.initial_c.is_some() {
                        return Err(format!(
                            "{site}: a structure has no temperature of its own, so a region \
                             cannot start at one — its expansion comes from the block it `follows`"
                        ));
                    }
                    if r.material == VOID {
                        return Err(format!(
                            "{site}: a structure has no void — an element with nothing in it has \
                             no stiffness, and an unsupported node is a rigid-body mode the solve \
                             cannot resolve"
                        ));
                    }
                    let empty = (0..3).any(|a| r.to[a] <= r.from[a]);
                    let outside = (0..3).any(|a| r.to[a] > cells[a]);
                    if empty || outside {
                        return Err(format!(
                            "{site}: {:?}..{:?} selects no elements of a {}x{}x{} body; `to` is \
                             one past the last, so a single element is `to = from + 1`",
                            r.from, r.to, cells[0], cells[1], cells[2]
                        ));
                    }
                    let substance = palette.get(&site, &r.material)?;
                    let inside = |i: usize, j: usize, k: usize| {
                        let p = [i, j, k];
                        (0..3).all(|a| p[a] >= r.from[a] && p[a] < r.to[a])
                    };
                    body.fill(elastic_of(&substance, &site)?, inside);
                    let expansion = substance.thermal.map_or(0.0, |t| t.expansion.to_si());
                    for k in 0..counts.2 {
                        for j in 0..counts.1 {
                            for i in 0..counts.0 {
                                if inside(i, j, k) {
                                    alpha[i + counts.0 * (j + counts.1 * k)] = expansion;
                                }
                            }
                        }
                    }
                }

                for h in held {
                    let face = match h.face {
                        FaceSpec::XMin => dualis::elastic::Face::XLow,
                        FaceSpec::XMax => dualis::elastic::Face::XHigh,
                        FaceSpec::YMin => dualis::elastic::Face::YLow,
                        FaceSpec::YMax => dualis::elastic::Face::YHigh,
                        FaceSpec::ZMin => dualis::elastic::Face::ZLow,
                        FaceSpec::ZMax => dualis::elastic::Face::ZHigh,
                    };
                    match h.how {
                        Hold::Clamp => body.clamp(face),
                        Hold::Roller => body.roller(face),
                    }
                }
                if body.free_dofs()
                    == body.node_counts().0 * body.node_counts().1 * body.node_counts().2 * 3
                {
                    return Err(format!(
                        "{name}: nothing is held, so this body has six rigid-body modes and no \
                         unique displacement. Hold at least three faces — rollers on three that \
                         meet at a corner leave it free to change size in every direction"
                    ));
                }

                for pr in pressed {
                    let face = match pr.face {
                        FaceSpec::XMin => dualis::elastic::Face::XLow,
                        FaceSpec::XMax => dualis::elastic::Face::XHigh,
                        FaceSpec::YMin => dualis::elastic::Face::YLow,
                        FaceSpec::YMax => dualis::elastic::Face::YHigh,
                        FaceSpec::ZMin => dualis::elastic::Face::ZLow,
                        FaceSpec::ZMax => dualis::elastic::Face::ZHigh,
                    };
                    if !pr.mpa.is_finite() {
                        return Err(format!("{name}: {} is not a pressure", pr.mpa));
                    }
                    body.press(face, Pressure::from_si(pr.mpa * 1e6));
                }

                match (follows, reference_c) {
                    (Some(_), None) => {
                        return Err(format!(
                            "{name}: a structure that follows a block needs `reference_c`, the \
                             temperature at which it is the size it was drawn. There is no \
                             default — a module assembled at its solder's reflow temperature is \
                             already stressed at room temperature, and guessing which is meant \
                             would change the answer without saying so"
                        ));
                    }
                    (None, Some(_)) => {
                        return Err(format!(
                            "{name}: `reference_c` says when this body is stress-free, which only \
                             means something if it `follows` a block whose temperature moves it"
                        ));
                    }
                    _ => {}
                }
                if follows.is_some() && alpha.iter().all(|a| *a == 0.0) {
                    notes.push(format!(
                        "{name}: follows a block, and nothing it is made of has an expansion \
                         coefficient — the temperature will move it by exactly nothing"
                    ));
                }

                Box::new(body)
            }
            DomainSpec::Block {
                name,
                cells,
                cell_mm,
                initial_c,
                material,
                regions,
                hot_spot,
                parts,
                cooling,
                dissipation,
            } => {
                let bulk = palette.get(name, material.as_deref().unwrap_or("aluminium"))?;
                let mut block = dualis::thermal::Solid3D::new(
                    name.clone(),
                    bulk,
                    (cells[0], cells[1], cells[2]),
                    Length::mm(*cell_mm),
                    Temperature::celsius(*initial_c),
                );
                // Designed parts first, and everything outside them is nothing. Before
                // `regions`, because a region is a box stated against the block's own grid and
                // is the way to say "this corner of that part is something else" — the same
                // last-writer-wins order the format already documents.
                if !parts.is_empty() {
                    let origin = LengthVec::m(0.0, 0.0, 0.0);
                    let counts = (cells[0], cells[1], cells[2]);
                    let mut occupied = vec![false; cells[0] * cells[1] * cells[2]];
                    for (n, part) in parts.iter().enumerate() {
                        let site = format!("{name}/parts[{n}]");
                        let bytes = files.bytes(&part.stl).map_err(|e| format!("{site}: {e}"))?;
                        let mesh = dualis::shape::Mesh::from_stl(&bytes)
                            .map_err(|e| format!("{site}: {}: {e}", part.stl))?;
                        let voxels = dualis::shape::Voxels::onto(
                            &mesh,
                            origin,
                            counts,
                            Length::mm(*cell_mm),
                        )
                        .map_err(|e| format!("{site}: {e}"))?;
                        let substance = palette.get(&site, &part.material)?;
                        // What voxelising cost, reported rather than left in the crate. The
                        // number with no symptom is the volume error: a rib finer than the grid
                        // does not fail, it disappears, and the run is perfectly well behaved
                        // about a different object.
                        let loss = voxels.loss();
                        notes.push(format!(
                            "{site}: {} filled {} cells at {cell_mm} mm — volume {:+.2}%, {:.0}% of it \
                             in boundary cells, {} thin run(s), {} triangle(s) under a cell{}",
                            part.stl,
                            voxels.filled(),
                            loss.volume_error * 100.0,
                            loss.boundary_fraction * 100.0,
                            loss.thin_runs,
                            loss.small_triangles,
                            if loss.ambiguous_rows > 0 {
                                format!(", {} ROW(S) THE RASTERISER COULD NOT DECIDE", loss.ambiguous_rows)
                            } else {
                                String::new()
                            }
                        ));
                        block.fill(substance, |i, j, k| voxels.contains(i, j, k));
                        for k in 0..counts.2 {
                            for j in 0..counts.1 {
                                for i in 0..counts.0 {
                                    if voxels.contains(i, j, k) {
                                        occupied[i + counts.0 * (j + counts.1 * k)] = true;
                                    }
                                }
                            }
                        }
                    }
                    // A part that voxelised to nothing at all is refused. It would leave a
                    // block of pure void that runs, conserves and answers about an empty box —
                    // the same shape as a region that selects no cells, and refused for the
                    // same reason.
                    if !occupied.iter().any(|o| *o) {
                        return Err(format!(
                            "{name}: every part voxelised to no cells at {cell_mm} mm. The \
                             meshes may be smaller than one cell, or outside the block"
                        ));
                    }
                    block = block.empty(|i, j, k| !occupied[i + counts.0 * (j + counts.1 * k)]);
                }
                for (n, r) in regions.iter().enumerate() {
                    let site = format!("{name}/regions[{n}]");
                    let empty = (0..3).any(|a| r.to[a] <= r.from[a]);
                    let outside = (0..3).any(|a| r.to[a] > cells[a]);
                    if empty || outside {
                        return Err(format!(
                            "{site}: {:?}..{:?} selects no cells of a {}x{}x{} block; `to` is one \
                             past the last cell, so a single cell is `to = from + 1`",
                            r.from, r.to, cells[0], cells[1], cells[2]
                        ));
                    }
                    let inside = |i: usize, j: usize, k: usize| {
                        let p = [i, j, k];
                        (0..3).all(|a| p[a] >= r.from[a] && p[a] < r.to[a])
                    };
                    // A box of nothing, which is how a scene states an air gap without carrying
                    // a mesh. Applied in the same order as any other region — later wins where
                    // they overlap — so a gap cut into a part reads the way a coating on a layer
                    // does.
                    if r.material == VOID {
                        if r.initial_c.is_some() {
                            return Err(format!(
                                "{site}: a {VOID:?} region cannot start at a temperature, because there is nothing there to have one"
                            ));
                        }
                        block = block.empty(inside);
                        continue;
                    }
                    let substance = palette.get(&site, &r.material)?;
                    block.fill(substance, inside);
                    if let Some(celsius) = r.initial_c {
                        if !celsius.is_finite() {
                            return Err(format!("{site}: initial_c is {celsius}"));
                        }
                        for k in 0..cells[2] {
                            for j in 0..cells[1] {
                                for i in 0..cells[0] {
                                    if inside(i, j, k) {
                                        block.set_temperature(
                                            i,
                                            j,
                                            k,
                                            Temperature::celsius(celsius),
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
                // Exposed faces, before the hot spot so a refusal names the file's own
                // mistake rather than arriving after the block has been half set up.
                let mut seen = std::collections::BTreeSet::new();
                for (n, cool) in cooling.iter().enumerate() {
                    let site = format!("{name}/cooling[{n}]");
                    cool.check(&site)?;
                    if !seen.insert(cool.face) {
                        return Err(format!(
                            "{site}: {:?} is already cooled by an earlier entry, and a face \
                             cannot lose heat to two different airs — add the areas, or state \
                             the one film that is really there",
                            cool.face
                        ));
                    }
                    block = block.losing_from(
                        cool.face.to_face(),
                        Environment {
                            ambient: Temperature::celsius(cool.ambient_c),
                            convection_w_per_m2_k: cool.convection_w_per_m2_k,
                            area: Area::from_si(cool.area_cm2 * 1e-4),
                        },
                    );
                }
                if let Some(spot) = hot_spot {
                    block.set_temperature(
                        spot.at[0],
                        spot.at[1],
                        spot.at[2],
                        Temperature::celsius(initial_c + spot.above_k),
                    );
                }
                // **How much of a clearance's answer is a boundary condition.** A gap is charged
                // as two infinite parallel plates, which is exact when the sides of the gap are
                // mirrors — and the block's own outer faces are, because an insulated boundary is
                // implemented as one. A user who meant the gap to be *open* — two parts floating
                // in vacuum rather than inside a housing — would get the plain view factor, and
                // the two agree only when the gap is narrow compared with the surfaces. Reported
                // rather than chosen: which reading is right is a statement about the geometry
                // that the grid cannot make, and a number is worth more than a caveat.
                for (n, d) in dissipation.iter().enumerate() {
                    let site = format!("{name}/dissipation[{n}]");
                    let empty = (0..3).any(|a| d.to[a] <= d.from[a]);
                    let outside = (0..3).any(|a| d.to[a] > cells[a]);
                    if empty || outside {
                        return Err(format!(
                            "{site}: {:?}..{:?} selects no cells of a {}x{}x{} block; `to` is one \
                             past the last cell, so a single cell is `to = from + 1`",
                            d.from, d.to, cells[0], cells[1], cells[2]
                        ));
                    }
                    if !d.watts.is_finite() {
                        return Err(format!("{site}: {} is not a number of watts", d.watts));
                    }
                    let inside = |i: usize, j: usize, k: usize| {
                        let p = [i, j, k];
                        (0..3).all(|a| p[a] >= d.from[a] && p[a] < d.to[a])
                    };
                    // Everything in the box is nothing, so the watts have nowhere to go. The
                    // library allows that state because a caller building an assembly cell by cell
                    // passes through it; a scene may not, because a file saying 45 W and meaning
                    // none is the quiet kind of wrong.
                    let solid = block.cells_on_where(&inside);
                    if solid == 0 {
                        return Err(format!(
                            "{site}: every cell in {:?}..{:?} is empty, so the {} W it states \
                             would be generated nowhere",
                            d.from, d.to, d.watts
                        ));
                    }
                    block = block.dissipating(d.watts, inside);
                }

                for patch in block.gap_patches() {
                    let f = patch.view_factor();
                    if f >= 0.95 {
                        continue;
                    }
                    notes.push(format!(
                        "{name}: a clearance of {:.1} mm across {:.1} x {:.1} mm carries the infinite-plate \
                         exchange, exact for the mirrored sides this block has. Open to space the \
                         same pair would see {:.3} of each other, so that reading is worth {:.1}x here{}",
                        patch.distance * 1e3,
                        patch.span.0 * 1e3,
                        patch.span.1 * 1e3,
                        f,
                        1.0 / f.max(1e-12),
                        if patch.rectangular {
                            ""
                        } else {
                            ", and the sheet is not a rectangle so that factor is an upper bound"
                        }
                    ));
                }
                Box::new(block)
            }
        };
        Ok(domain)
    }
}

/// One designed part filled into a [`DomainSpec::Block`] from a mesh.
///
/// # An STL is read as millimetres, and this is the way an assembly goes wrong
///
/// The format is unitless and every CAD tool writes millimetres, so that is what
/// `Mesh::from_stl` reads. A file exported in metres arrives a **thousand times too small** and
/// voxelises to no cells at all, which this format refuses — but the refusal says the meshes
/// may be smaller than a cell, and the reason it is usually right is this one.
///
/// # Where a part sits is what its own file says
///
/// An STL carries absolute coordinates, so two parts exported from one assembly land where the
/// assembly put them, in the block's own axes. There is no offset here to get wrong and no pose:
/// `poses` places a **domain**, and every part of an assembly is inside one domain.
///
/// A part whose mesh reaches outside the block's cells is **refused with both boxes named**,
/// rather than cropped. A part with its corner cut off runs, audits, renders and answers about a
/// different shape.
///
/// # What the cell size costs, and why the report is printed
///
/// Voxelising is where a designed shape meets a grid, and the thing with no symptom is *how much
/// of the shape survives*: a 0.5 mm rib at 2 mm cells is not a thin rib, it is gone, and the
/// simulation runs perfectly well without it. `--check` prints the [`dualis::shape::Loss`]
/// for every part — volume error, the share of the volume in boundary cells, runs one or two
/// cells thick, and rows the rasteriser could not decide — because a number nobody sees is a
/// number nobody acts on.
///
/// # The file is read from disk, and a browser has no disk
///
/// `stl` is a path, resolved when the scene is built. That works for the CLI and the native
/// editor and **fails in the browser**, where there is no filesystem and the error says which
/// file it could not find. Stated here rather than discovered: a page that wants geometry has to
/// carry it, and that is a format decision this key does not make.
///
/// ```json
/// "parts": [
///   { "stl": "bracket.stl", "material": "aluminium" },
///   { "stl": "insert.stl",  "material": "copper" }
/// ]
/// ```
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PartSpec {
    /// Path to a binary or ASCII STL, relative to wherever the scene is run from.
    pub stl: String,
    /// What the part is made of: a catalogue name or one of [`Scene::materials`].
    pub material: String,
}

/// How a face of a structure is held.
///
/// The distinction is not decoration and it changes the answer by a factor: a **clamp** holds all
/// three components, a **roller** holds only the one normal to the face. A bar between two rollers
/// cannot lengthen and is free to fatten, which is what "held along one axis" means and gives
/// `σ = −Eε₀`; the same bar between two clamps also has its ends held sideways, and comes out
/// **1.82 times** stiffer — a three-dimensional answer to a different question.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Hold {
    /// All three components held: a built-in support.
    Clamp,
    /// Only the normal component held: a symmetry plane, or a surface it can slide along.
    Roller,
}

/// A face of a structure that is held, and how.
///
/// ```json
/// "held": [
///   { "face": "z-min", "as": "clamp" },
///   { "face": "x-min", "as": "roller" }
/// ]
/// ```
///
/// **A structure with nothing held cannot be solved**, and is refused rather than returning a
/// number: an elliptic system with no support has six rigid-body modes, the operator is singular
/// along them, and a solver asked for a displacement would answer with whatever the iteration
/// happened to drift into. Three rollers on three faces meeting at a corner is the minimum, and
/// leaves a body free to change size in every direction.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct HeldSpec {
    /// Which face. Same spelling as [`FaceSpec`].
    pub face: FaceSpec,
    /// Clamp or roller.
    #[serde(rename = "as")]
    pub how: Hold,
}

/// A face of a structure under pressure, in megapascals.
///
/// ```json
/// "pressed": [ { "face": "z-max", "mpa": 1.5 } ]
/// ```
///
/// Positive presses **into** the body and negative pulls on it, which is the sign convention of
/// pressure everywhere else and the opposite of the stress convention — stated because getting it
/// backwards produces a perfectly well-behaved answer to the other question.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct PressedSpec {
    /// Which face.
    pub face: FaceSpec,
    /// Pressure in megapascals, positive inward.
    pub mpa: f64,
}

/// What is flowing: a name, or the two numbers that define one.
///
/// ```json
/// "fluid": "water"
/// "fluid": { "density": 998.0, "kinematic_viscosity": 1.004e-6 }
/// ```
///
/// **Kinematic** viscosity, `μ/ρ`, in m²/s — the one that appears in the equations and in every
/// closed form. Tables quote the dynamic viscosity as often, and for water the two differ by a
/// factor of a thousand: the sort of error dimensions cannot catch and a name can, which is why
/// the key spells it out rather than saying `viscosity`.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, PartialEq)]
#[serde(untagged)]
pub enum FluidSpec {
    /// One of `water` or `air`, at twenty degrees.
    Named(String),
    /// Stated outright.
    Stated {
        /// Density, kg/m³.
        density: f64,
        /// Kinematic viscosity, m²/s.
        kinematic_viscosity: f64,
    },
}

/// The walls of a channel, and how fast they are sliding.
///
/// ```json
/// "walls": { "lower_m_per_s": 0.0, "upper_m_per_s": 0.0 }
/// ```
///
/// No-slip at `y = 0` and `y = h`. Both zero is a **channel** and Poiseuille flow when it is
/// driven; one moving and no drive is **Couette**. Absent is a periodic box in every direction,
/// which is what Taylor–Green lives in and is a different problem rather than a simpler one.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct WallsSpec {
    /// Speed of the `y = 0` wall along x, m/s.
    #[serde(default)]
    pub lower_m_per_s: f64,
    /// Speed of the `y = h` wall along x, m/s.
    #[serde(default)]
    pub upper_m_per_s: f64,
}

/// A box of cells that generates heat, in watts.
///
/// ```json
/// "dissipation": [
///   { "watts": 45.0, "from": [3, 3, 6], "to": [6, 6, 7] }
/// ]
/// ```
///
/// The watts are the **total** for the box, not a figure per cell: a die dissipating 45 W
/// dissipates 45 W whether the grid gives it eight cells or eight thousand. That is what makes a
/// dissipating scene survive `verify`'s resolution sweep, and it is the opposite of the choice a
/// per-cell figure would force — the answer would then grow with the mesh and there would be no
/// convergence to measure.
///
/// Void cells inside the box generate nothing and take no share, so a box drawn around a part and
/// its clearance heats the part at the full rate. A box that selects only void, or no cell at all,
/// is **refused**: a scene that says 45 W and means none is a mistake worth naming.
///
/// Boxes may overlap. Later wins, the same last-writer rule `regions` and `parts` follow, because
/// two files saying the same thing have to mean the same thing.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct DissipationSpec {
    /// Total watts generated in this box.
    pub watts: f64,
    /// The first cell of the box, inclusive.
    pub from: [usize; 3],
    /// One past the last cell, so a single cell is `to = from + 1`.
    pub to: [usize; 3],
}

/// One face of a [`DomainSpec::Block`] losing heat to still or moving air.
///
/// # What the numbers mean
///
/// `area_cm2` is the **whole face's** area, not a cell's, and the domain charges each boundary
/// cell its share. That is what keeps the number a scene writes independent of the grid it
/// chose: refining `cells` does not change what the file says the part exposes, which is the
/// property `verify`'s resolution sweep needs to compare two grids of one problem rather than
/// two problems.
///
/// `convection_w_per_m2_k` has no default and there is no right one — still air is about 5 to
/// 10, forced air 25 to 100, water hundreds. The **radiative** half is not stated here at all:
/// it comes from the material's emissivity, because that is a property of the surface and not
/// of the air beside it. A black surface at room temperature radiates about 6 W·m⁻²·K⁻¹, the
/// same order as still air, so a scene that named only convection would be wrong by about a
/// factor of two and would look reasonable while it was.
///
/// ```json
/// "cooling": [
///   { "face": "z-max", "ambient_c": 20.0, "convection_w_per_m2_k": 25.0, "area_cm2": 4.0 }
/// ]
/// ```
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CoolingSpec {
    /// Which face.
    pub face: FaceSpec,
    /// The temperature of whatever it loses to.
    pub ambient_c: f64,
    /// Convective coefficient, W·m⁻²·K⁻¹.
    pub convection_w_per_m2_k: f64,
    /// The whole face's area, in square centimetres.
    pub area_cm2: f64,
}

/// Which outer face of a block a [`CoolingSpec`] is about.
///
/// The block's own axes, which a `poses` entry is what places in the world. Spelled out rather
/// than indexed because a file is read by a person and `faces[3]` is not.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub enum FaceSpec {
    /// The `x = 0` face.
    XMin,
    /// The far face along `x`.
    XMax,
    /// The `y = 0` face.
    YMin,
    /// The far face along `y`.
    YMax,
    /// The `z = 0` face.
    ZMin,
    /// The far face along `z`.
    ZMax,
}

impl FaceSpec {
    /// The domain's own face.
    fn to_face(self) -> dualis::thermal::Face {
        use dualis::thermal::Face;
        match self {
            FaceSpec::XMin => Face::XMin,
            FaceSpec::XMax => Face::XMax,
            FaceSpec::YMin => Face::YMin,
            FaceSpec::YMax => Face::YMax,
            FaceSpec::ZMin => Face::ZMin,
            FaceSpec::ZMax => Face::ZMax,
        }
    }
}

impl CoolingSpec {
    /// Refuse a face that does not describe a surface losing heat.
    fn check(&self, site: &str) -> Result<(), String> {
        //  first, so the comparison that follows is between two numbers and the
        // NaN case is not being caught by a negated inequality — which is what clippy asks
        // for and is also the clearer reading.
        if !self.area_cm2.is_finite() || self.area_cm2 <= 0.0 {
            return Err(format!(
                "{site}: area_cm2 is {}, and a face with no area loses nothing — which would \
                 run as an insulated block and answer a different question in silence",
                self.area_cm2
            ));
        }
        if !self.convection_w_per_m2_k.is_finite() || self.convection_w_per_m2_k < 0.0 {
            return Err(format!(
                "{site}: convection_w_per_m2_k is {}, and a negative film would carry heat \
                 uphill",
                self.convection_w_per_m2_k
            ));
        }
        if !self.ambient_c.is_finite() {
            return Err(format!("{site}: ambient_c is {}", self.ambient_c));
        }
        Ok(())
    }
}

/// A box of cells in a [`DomainSpec::Block`] made of something other than the block.
///
/// # Half-open, and empty is an error
///
/// `to` is one past the last cell, so a single cell is `to = from + 1` and a twelve-cell layer along
/// z is `from: [0,0,12], to: [1,1,24]`. That is the convention every index in this workspace uses,
/// and mixing it with an inclusive one somewhere is worse than either.
///
/// A region that selects **no cells** is refused rather than ignored. It is the whole silent failure
/// this format can have: a mistyped bound gives a block of one material that runs, audits, renders
/// and answers the wrong question, with nothing anywhere saying the coating was not applied.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Region {
    /// What this box is made of, one of [`MATERIALS`] — or [`VOID`], for a box of nothing.
    ///
    /// `"void"` is the one reserved spelling in this format. It is how a scene says *air gap*
    /// without carrying a mesh: the cells hold no substance, conduct to nothing, take no share
    /// of the bus and have no temperature. Two parts either side of one exchange **radiation**
    /// and not conduction — see `Solid3D::empty` for what that models and what it leaves out.
    ///
    /// Reserved rather than resolved, so a scene cannot declare a material called `void` and
    /// quietly mean a solid: [`Palette`] refuses the name for the same reason it refuses a
    /// catalogue spelling, which is that two files saying the same word have to mean the same
    /// thing.
    pub material: String,
    /// The first cell, as `[i, j, k]`.
    pub from: [usize; 3],
    /// One past the last cell, as `[i, j, k]`.
    pub to: [usize; 3],
    /// What this box starts at, in celsius. Absent is the block's own `initial_c`.
    ///
    /// **A region could say what a box was made of and not what it held**, and that gap was
    /// found the way this crate is meant to find things: by trying to write a scene of a hot
    /// part beside a cool lid and discovering the format could say it in neither direction.
    /// Heat off the bus is spread over the whole block by design — the plain channel carries an
    /// amount and no location — so a source cannot warm one part of an assembly either.
    ///
    /// Like [`HotSpot`], this is **a statement about the initial state and not a delivery of
    /// heat**: it moves what the block holds, not what it has absorbed, so the audit's opening
    /// balance includes it.
    ///
    /// Refused on a `void` region, because nothing has no temperature.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial_c: Option<f64>,
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
    // `pub(crate)` for the `verify` module, whose instrumented loop reads the ledger and the
    // stability limits between the advances `World::run` makes without measuring.
    pub(crate) sim: Simulation,
    /// What the build dismissed and why — see [`World::notes`].
    notes: Vec<String>,
    /// Linear expansion per kelvin, per element, for each structure that follows a block.
    ///
    /// Kept here rather than in the elastic domain because a domain carrying a thermal coefficient
    /// would be a domain that knew about temperature, and rule 4 says it may not. The coupling is
    /// the application's to make, and this is the application.
    expansion: BTreeMap<String, Vec<f64>>,
}

impl World {
    /// Build a simulation from a scene, reading any `parts` from the **filesystem**.
    ///
    /// Fails only on a scene that cannot describe a simulation at all — no domains, a
    /// non-positive duration, no frames. Physical nonsense inside a domain is the domain's
    /// business and is reported by the audit at run time, which is where it belongs.
    ///
    /// See [`World::build_with`] for a scene whose parts are not on a disk, which is every scene
    /// in a browser.
    pub fn build(scene: Scene) -> Result<World, String> {
        Self::build_with(scene, &OnDisk)
    }

    /// Build a simulation from a scene, taking `parts` from `files`.
    ///
    /// The same builder, the same audit, the same messages — the only thing that changes is where
    /// a part's bytes come from. That is deliberate: a scene that runs from a page and a scene that
    /// runs from a terminal have to be the *same scene*, or the browser is a demo rather than the
    /// product.
    pub fn build_with(scene: Scene, files: &dyn Parts) -> Result<World, String> {
        // Before anything else, and before any field is read: a file from a newer build must not
        // be half-run.
        scene.check_version()?;
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
                    "{name}: tracks names the network {net_name:?}, which this scene does not define"
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
        for (name, tol) in &scene.tolerance_for {
            // Matched against the kernel's constants rather than passed through, because
            // `Tolerances::with` takes `&'static str` on purpose: two spellings of one channel are
            // two channels, and a scene file is exactly where a second spelling comes from.
            let channel = match name.as_str() {
                "energy" => quantity::ENERGY,
                "momentum" => quantity::MOMENTUM,
                "mass" => quantity::MASS,
                "charge" => quantity::CHARGE,
                "photons" => quantity::PHOTONS,
                other => {
                    return Err(format!(
                        "tolerance_for names an unknown quantity {other:?}; known: energy, momentum, mass, charge, photons"
                    ))
                }
            };
            sim = sim.conservation_tolerance_for(channel, *tol);
        }

        // Two domains under one name is a scene that cannot describe a simulation, which is
        // what this function is for. `Simulation::domain` takes the *first* match, so without
        // this the second domain is never sampled and the first is drawn twice under the
        // second's geometry and label — a 500 C bar reported as 20 C, twice, with no warning.
        for (i, spec) in scene.domains.iter().enumerate() {
            if let Some(earlier) = scene.domains[..i].iter().find(|d| d.name() == spec.name()) {
                return Err(format!(
                    "two domains are both called {:?}; names are how they are looked up, so the second would be invisible",
                    earlier.name()
                ));
            }
        }

        // Validated before any domain is built, so an impossible substance is reported as the
        // declaration it is rather than as a `NaN` somewhere downstream.
        if let Some(env) = &scene.environment {
            env.check()?;
        }
        // A pose naming a domain the scene does not define would place nothing, silently — the
        // same shape as `tracks` pointing at a node that is not there, and refused the same way,
        // with the file's own vocabulary quoted back.
        for (name, pose) in &scene.poses {
            if !scene.domains.iter().any(|d| d.name() == name) {
                let known: Vec<&str> = scene.domains.iter().map(|d| d.name()).collect();
                return Err(format!(
                    "poses.{name}: this scene defines no domain called {name:?}; it defines {}",
                    known.join(", ")
                ));
            }
            pose.to_pose(name)?;
        }
        let mut palette = Palette::with_composites(&scene.materials, &scene.composites)?;
        let mut notes = Vec::new();
        let mut expansion: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        for spec in &scene.domains {
            sim = sim.with_boxed(spec.build(
                &mut palette,
                scene.environment.as_ref(),
                files,
                &mut notes,
            )?);
            // The dismissals — a stated condition a domain ignores for a measured reason.
            // Collected here, in the composition root, because the reason is about the pair
            // (this stage, this domain) and neither owns it alone. Reported by `--check` and
            // the editor, because a dismissal nobody can see is indistinguishable from the
            // silence it exists to replace.
            if let (DomainSpec::Atoms { name, .. }, Some(env)) = (spec, &scene.environment) {
                if env.gravity_m_per_s2 != 0.0 {
                    notes.push(format!(
                        "{name}: the stated gravity is dismissed at this scale — the \
                         gravitational energy across one molecular diameter, m·g·σ, is about \
                         1.3e-13 of the Lennard-Jones well depth for argon, and the fluid's \
                         physics would not move by a pixel"
                    ));
                }
            }
        }
        // **What a structure that follows a block needs, checked here rather than at the first
        // step.** A build-time refusal names the file; a step-time one names a number.
        for spec in &scene.domains {
            let DomainSpec::Structure {
                name,
                cells,
                cell_mm,
                material,
                regions,
                follows: Some(block),
                ..
            } = spec
            else {
                continue;
            };
            let source = scene
                .domains
                .iter()
                .find(|d| d.name() == block)
                .ok_or_else(|| {
                    format!(
                        "{name}: follows {block:?}, which is not a domain in this scene. The \
                         scene has {}",
                        scene
                            .domains
                            .iter()
                            .map(DomainSpec::name)
                            .collect::<Vec<_>>()
                            .join(", ")
                    )
                })?;
            let DomainSpec::Block {
                cells: bc,
                cell_mm: bmm,
                ..
            } = source
            else {
                return Err(format!(
                    "{name}: follows {block:?}, which is not a block — only a block has a \
                     temperature field for a structure to expand with"
                ));
            };
            // **The grids must be the same grid.** An element and a cell are the same box, or the
            // coupling is a guess about which cell a corner belongs to — and an interpolation
            // invented in a builder would be a physics decision nobody stated.
            if bc != cells || (bmm - cell_mm).abs() > 1e-12 {
                return Err(format!(
                    "{name}: follows {block:?} but their grids differ — {cells:?} at {cell_mm} mm \
                     against {bc:?} at {bmm} mm. An element and a cell have to be the same box"
                ));
            }
            // The expansion coefficients, resolved the way the body's materials were.
            let base = palette.get(name, material.as_deref().unwrap_or("aluminium"))?;
            let counts = (cells[0], cells[1], cells[2]);
            let mut alpha = vec![
                base.thermal.map_or(0.0, |t| t.expansion.to_si());
                counts.0 * counts.1 * counts.2
            ];
            for (n, r) in regions.iter().enumerate() {
                let site = format!("{name}/regions[{n}]");
                let substance = palette.get(&site, &r.material)?;
                let e = substance.thermal.map_or(0.0, |t| t.expansion.to_si());
                for k in 0..counts.2 {
                    for j in 0..counts.1 {
                        for i in 0..counts.0 {
                            let p = [i, j, k];
                            if (0..3).all(|a| p[a] >= r.from[a] && p[a] < r.to[a]) {
                                alpha[i + counts.0 * (j + counts.1 * k)] = e;
                            }
                        }
                    }
                }
            }
            expansion.insert(name.clone(), alpha);
        }

        // And after, because "used" means a domain asked for it, which is only known once they all
        // have. A block's `material` is optional and defaults to aluminium, so a declaration nothing
        // named is a run on the wrong substance that reports nothing.
        let unused = palette.unused();
        if !unused.is_empty() {
            return Err(format!(
                "materials: {} declared and used by nothing — a block with no `material` runs as aluminium, so this scene would answer about a material it did not mean. Name it or delete it",
                unused
                    .iter()
                    .map(|u| format!("{u:?}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }

        // **The coupling must be reachable, and a silent failure here is invisible.** Every
        // structure that follows a block is written to between steps through
        // `Simulation::domain_as_mut`, which returns `None` for a domain that has not implemented
        // `as_any_mut` — an opt-in trait method that defaults to nothing. Measured before this
        // check existed: `dualis-elastic::Block` implemented `as_any` and not `as_any_mut`, so the
        // whole coupling did nothing and the scene reported **zero strain, zero stress and zero
        // strain energy** — which reads as *no thermal stress* rather than as *not connected*, and
        // is the more believable of the two.
        for (structure, block) in scene.domains.iter().filter_map(|d| match d {
            DomainSpec::Structure {
                name,
                follows: Some(b),
                ..
            } => Some((name.clone(), b.clone())),
            _ => None,
        }) {
            if sim.domain_as::<dualis::thermal::Solid3D>(&block).is_none() {
                return Err(format!(
                    "{structure}: {block:?} cannot be read back as a thermal block, so its                      temperature cannot drive anything"
                ));
            }
            if sim
                .domain_as_mut::<dualis::elastic::Block>(&structure)
                .is_none()
            {
                return Err(format!(
                    "{structure}: this body cannot be written to between steps, so the                      temperature would reach it and change nothing. A domain has to implement                      `Domain::as_any_mut` to be coupled into"
                ));
            }
        }

        let mut world = World {
            scene,
            sim,
            notes,
            expansion,
        };
        // **Once before anything runs**, so the first captured frame already carries the strain the
        // starting temperature implies. Not a formality: a power module is assembled at its
        // solder's reflow temperature and sits at room temperature before it is switched on, so it
        // is **already stressed** at `t = 0` and a frame showing zero there would be the one
        // untrue frame in the run. It also fixes the reading set from the start, which is what
        // keeps a conditional column from appearing halfway through a table.
        world.close_expansion();
        world.resolve_structures();
        Ok(world)
    }

    /// The scene this was built from.
    pub fn scene(&self) -> &Scene {
        &self.scene
    }

    /// What the build wants a reader to know and the run will not say.
    ///
    /// Two kinds so far, and both are things that otherwise happen in silence. A **dismissal**:
    /// a stated condition a domain correctly ignored, with the measurement that earns it — a
    /// molecular fluid under stated gravity. And a **rasterisation report**: what a designed
    /// part cost when it met the grid, because a rib finer than a cell does not fail, it
    /// disappears, and the run is perfectly well behaved about a different object.
    pub fn notes(&self) -> &[String] {
        &self.notes
    }

    /// The simulation underneath, for a caller that wants the kernel's own accessors —
    /// `bus`, `ledger`, `domain_as`, `field`.
    pub fn simulation(&self) -> &Simulation {
        &self.sim
    }

    /// The simulation underneath, mutably, for a caller that wants to write to a domain
    /// **between** steps.
    ///
    /// The same allowance `Domain::as_any_mut` documents and for the same reason: this is the
    /// owner of the simulation, outside the step loop, holding `&mut World` already — it could
    /// drop the world and rebuild it, so denying it a write was never protecting anything. What
    /// it does *not* do is let one domain read another inside `step`, which is the rule the
    /// crate split exists to hold.
    ///
    /// Setting an initial condition a scene file cannot state is the use that asked for it: a
    /// whole part warmed rather than one cell, which `hot_spot` cannot say.
    pub fn simulation_mut(&mut self) -> &mut Simulation {
        &mut self.sim
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
            self.advance(dt)?;
            frames.push(dualis::scene::capture(&self.sim, &placed));
        }
        dualis::scene::settle_framing(&mut frames);
        Ok(frames)
    }

    /// Advance the clock by `dt` and close the between-frames feedback — exactly one
    /// iteration of [`World::run`]'s loop, made public.
    ///
    /// This is how a caller runs a world **frame by frame**: capture between advances with
    /// [`dualis::scene::capture`], stream a long run to a screen as it computes, or stop one
    /// early. The `verify` battery needed this loop first and reached it through crate
    /// privacy; the editor needed it second, from outside, and a loop two consumers need is
    /// an API rather than an implementation detail. The feedback is closed inside rather than
    /// left to the caller, because a caller who forgets it gets a winding whose resistance
    /// never moves — a run that is wrong in exactly the way nothing reports.
    pub fn advance(&mut self, dt: Time) -> Result<dualis::core::Report, Violation> {
        let report = self.sim.advance(dt)?;
        self.close_feedback();
        // **After the advance, and then solved again.** The capture that follows shows the block's
        // new temperature, so the body beside it has to answer *that* temperature and not the one
        // before. Setting the strain first instead left the structure one window stale — invisible
        // near a steady state and a 0.5% error on a 60 s run, which is exactly the size of lag that
        // gets read as a physics result.
        self.close_expansion();
        self.resolve_structures();
        Ok(report)
    }

    /// Hand each structure the stress-free strain its block's temperature implies.
    ///
    /// `α(T − T_ref)` element by element. This is the coupling a digital twin is for, and it lives
    /// here rather than in either domain because rule 4 forbids a domain knowing about another:
    /// `dualis-elastic` takes an eigenstrain and never learns it came from a temperature, and
    /// `dualis-thermal` never learns that anybody read it.
    ///
    /// A void cell has no temperature, and an element over one gets **no strain** rather than a
    /// `NaN` — which would otherwise reach the load vector and bring the whole solve back
    /// not-a-number, reported as a failure to converge rather than as what it was.
    /// Bring every structure back to equilibrium after its strain changed.
    ///
    /// Called once at build so the first captured frame is a solved body rather than one handed a
    /// strain and not yet asked — which reported the strain energy of a fully constrained body,
    /// strains of zero and an infinite residual: three readings disagreeing about one instant.
    /// And after every advance, so a frame's body is the answer to that frame's temperature.
    pub(crate) fn resolve_structures(&mut self) {
        let names: Vec<String> = self
            .scene
            .domains
            .iter()
            .filter_map(|d| match d {
                DomainSpec::Structure { name, .. } => Some(name.clone()),
                _ => None,
            })
            .collect();
        for name in names {
            if let Some(body) = self.sim.domain_as_mut::<dualis::elastic::Block>(&name) {
                body.solve(1e-10);
            }
        }
    }

    pub(crate) fn close_expansion(&mut self) {
        let wanted: Vec<(String, String, f64)> = self
            .scene
            .domains
            .iter()
            .filter_map(|spec| match spec {
                DomainSpec::Structure {
                    name,
                    follows: Some(block),
                    reference_c: Some(reference),
                    ..
                } => Some((name.clone(), block.clone(), *reference)),
                _ => None,
            })
            .collect();

        for (structure, block, reference_c) in wanted {
            let Some(alpha) = self.expansion.get(&structure).cloned() else {
                continue;
            };
            let reference = Temperature::celsius(reference_c).to_si();
            let Some(source) = self.sim.domain_as::<dualis::thermal::Solid3D>(&block) else {
                continue;
            };
            let (nx, ny, nz) = source.counts();
            let mut strain = vec![0.0; alpha.len()];
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        let at = i + nx * (j + ny * k);
                        if at >= strain.len() {
                            continue;
                        }
                        let t = source.temperature_at(i, j, k).to_si();
                        if t.is_finite() {
                            strain[at] = alpha[at] * (t - reference);
                        }
                    }
                }
            }
            if let Some(body) = self.sim.domain_as_mut::<dualis::elastic::Block>(&structure) {
                body.stress_free_strain(|i, j, k| {
                    strain.get(i + nx * (j + ny * k)).copied().unwrap_or(0.0)
                });
            }
        }
    }

    /// Where each domain sits, keyed by the name the simulation knows it under.
    ///
    /// Built once per run rather than per frame: a placement is a property of the scene, and a
    /// scene does not move while it is being simulated. Rebuilding it 240 times would also have
    /// made a placement that *did* drift look like a working feature.
    pub fn placements(&self) -> BTreeMap<String, Placement> {
        self.scene.placements()
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
    pub(crate) fn close_feedback(&mut self) {
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
            // A conductor is a volume, sampled at its own cells: the potential is what the
            // solve produces, and a picture of it is where the current is going.
            DomainSpec::Conductor { cells, cell_mm, .. } => Placement::field(Extent::volume(
                LengthVec::from_si(
                    glam::DVec3::new(cells[0] as f64, cells[1] as f64, cells[2] as f64)
                        * cell_mm
                        * 1e-3,
                ),
                cells[0],
                cells[1],
                cells[2],
            )),
            // A basket is a volume, sampled at its own cells. `y` is the flow axis, so a slice
            // at fixed `k` is a vertical cut through the bed — which is the picture.
            DomainSpec::Puck { cells, cell_mm, .. } => Placement::field(Extent::volume(
                LengthVec::from_si(
                    glam::DVec3::new(cells[0] as f64, cells[1] as f64, cells[2] as f64)
                        * cell_mm
                        * 1e-3,
                ),
                cells[0],
                cells[1],
                cells[2],
            )),
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
            // A structure's field is its **nodes**, not its elements: the displacement lives at
            // the corners and there is one more of them along every axis. Getting that off by one
            // would draw a body one element short and nothing would say so.
            DomainSpec::Structure { cells, cell_mm, .. } => Placement::field(Extent::volume(
                LengthVec::from_si(
                    glam::DVec3::new(cells[0] as f64, cells[1] as f64, cells[2] as f64)
                        * cell_mm
                        * 1e-3,
                ),
                cells[0] + 1,
                cells[1] + 1,
                cells[2] + 1,
            )),
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
            // A channel's field is its **speed**, sampled at its own cell centres. The caveat is
            // in the unit and in `Channel::as_field`'s doc rather than in a refusal to draw: a
            // magnitude shows where the fluid is moving fast and not which way it is going. Drawn
            // at all because this crate's own docs say "it looks like a fluid" is the easiest
            // wrong answer to accept, and the answer to that is closed forms **and** a picture.
            DomainSpec::Channel { cells, cell_mm, .. } => Placement::field(Extent::volume(
                LengthVec::from_si(
                    glam::DVec3::new(cells[0] as f64, cells[1] as f64, cells[2] as f64)
                        * cell_mm
                        * 1e-3,
                ),
                cells[0],
                cells[1],
                cells[2],
            )),
        }
    }
}
