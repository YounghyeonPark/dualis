//! Where simulated things are, and what a run of them looks like.
//!
//! The middle layer. [`dualis-core`](https://docs.rs/dualis-core) says what evolves and what it
//! conserves; a view says how to draw it; this says **where things sit and what shape their
//! output is**.
//!
//! # It knows no physics, and that is the point
//!
//! Nothing here names a domain. It asks each one what it offers —
//! [`as_field`](dualis_core::Domain::as_field) for a continuum,
//! [`as_bodies`](dualis_core::Domain::as_bodies) for a countable set,
//! [`readings`](dualis_core::Domain::readings) for scalars — and builds a [`Frame`] from the
//! answers. A physics that arrives tomorrow gets captured without this crate being edited, which
//! is the property the whole workspace is arranged around.
//!
//! That was not free. Pulling this layer out of the application found three places where it had
//! been matching on domain types instead: one for scalars, one for bodies, one for a field's
//! extent. The first two became trait methods on `Domain`. The third became [`Placement`].
//!
//! # Two kinds of placement, and why they are one type with two fields
//!
//! A [`Pose`] changes what the physics computes: two solids in contact, a grid rotated against
//! its neighbour. It lives in the kernel because physics needs it.
//!
//! A [`Placement::marker`] is a position handed to something that *has* no geometry, purely so a
//! viewer can put it somewhere. A thermal network node has a capacity and not a position, and a
//! conductance is not a distance — giving one a coordinate is a statement about a diagram.
//!
//! They are separated by **which crate they live in**, not by a naming convention. `Pose` is
//! below; `marker` is here, above every domain, where no physics can reach it. If they shared a
//! home a drawing coordinate would eventually arrive in a conductance and nothing would fail
//! loudly.

#![deny(missing_docs)]
#![forbid(unsafe_code)]

use dualis_core::{Pose, Reading, ScalarField, Simulation};
use dualis_units::{Length, LengthVec, Time};
use std::collections::BTreeMap;

/// The region a field occupies, in its own coordinates, and how finely to sample it.
///
/// A [`ScalarField`] is a function of position and **does not know
/// where it stops** — that is the right division, since a field that knew its own bounds would be
/// a mesh, but it means somebody has to say. This is where.
///
/// `samples` is a request rather than a property: the same field can be captured coarsely for a
/// thumbnail and finely for a paper, and neither is more true than the other.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Extent {
    /// The low corner, in the domain's own coordinates.
    pub min: LengthVec,
    /// The high corner.
    pub max: LengthVec,
    /// How many samples across and up. One row means a line rather than a plane.
    pub samples: (usize, usize),
}

impl Extent {
    /// A box from two corners, sampled `nx` by `ny`.
    pub fn new(min: LengthVec, max: LengthVec, nx: usize, ny: usize) -> Extent {
        Extent {
            min,
            max,
            samples: (nx.max(1), ny.max(1)),
        }
    }

    /// A line along x, for a domain with one dimension.
    pub fn line(length: Length, cells: usize) -> Extent {
        Extent::new(
            LengthVec::ZERO,
            LengthVec::from_si(glam::DVec3::new(length.to_si(), 0.0, 0.0)),
            cells.max(2),
            1,
        )
    }

    /// A rectangle in the x–y plane.
    pub fn rectangle(width: Length, height: Length, nx: usize, ny: usize) -> Extent {
        Extent::new(
            LengthVec::ZERO,
            LengthVec::from_si(glam::DVec3::new(width.to_si(), height.to_si(), 0.0)),
            nx,
            ny,
        )
    }
}

/// Where one domain sits in the world, and how big it is.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Placement {
    /// The rigid motion from the domain's coordinates into the world's. **Physical**: it changes
    /// what the physics computes once anything reads across it.
    pub pose: Pose,
    /// The region to sample, for a domain that is a field. `None` for everything else.
    pub extent: Option<Extent>,
    /// A position for a domain that has no geometry at all, so a viewer can place it on a
    /// diagram. **Presentational**, and unreachable from any physics — see the module docs.
    pub marker: Option<LengthVec>,
}

impl Placement {
    /// Placed by a pose, with no field and no marker.
    pub fn at(pose: Pose) -> Placement {
        Placement {
            pose,
            ..Placement::default()
        }
    }

    /// A field of this extent, at the origin.
    pub fn field(extent: Extent) -> Placement {
        Placement {
            extent: Some(extent),
            ..Placement::default()
        }
    }

    /// Somewhere to draw a domain that has nowhere to be.
    pub fn marked(marker: LengthVec) -> Placement {
        Placement {
            marker: Some(marker),
            ..Placement::default()
        }
    }

    /// The same placement, moved.
    pub fn with_pose(mut self, pose: Pose) -> Placement {
        self.pose = pose;
        self
    }
}

/// One instant of a run: every domain's output, in whatever shape it has.
#[derive(Clone, Debug)]
pub struct Frame {
    /// Simulation time, in seconds.
    pub time_s: f64,
    /// One per domain that had something to draw, in the order they were placed.
    pub panels: Vec<Panel>,
    /// Named scalars from every domain, drawable or not.
    pub readings: Vec<Reading>,
}

/// One domain, captured.
#[derive(Clone, Debug)]
pub struct Panel {
    /// Which domain this came from.
    pub name: String,
    /// What the values mean, for a legend.
    pub unit: &'static str,
    /// The shape of what was captured.
    pub data: PanelData,
}

/// A continuum sampled on a grid, or a finite number of bodies at positions.
///
/// Two shapes because domains genuinely are two kinds of thing, and collapsing them would mean
/// inventing a continuum for the bodies or a body count for the field.
#[derive(Clone, Debug)]
pub enum PanelData {
    /// A field, sampled onto a grid.
    Field {
        /// Samples across.
        nx: usize,
        /// Samples up. One for a line.
        ny: usize,
        /// Row-major, `nx * ny` values.
        values: Vec<f64>,
    },
    /// Bodies at positions **in world coordinates**, each with a value to colour it by.
    Points {
        /// Where each body is.
        positions: Vec<[f64; 3]>,
        /// One per body.
        values: Vec<f64>,
        /// `[x0, y0, z0, x1, y1, z1]` — the region to draw, fixed for the whole run by
        /// [`settle_framing`] so a body moving is a body moving and not the picture rescaling.
        bounds: [f64; 6],
        /// Whether that box is a **real wall** — a periodic cell — rather than a drawing margin.
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

    /// The grid shape, for a field.
    pub fn grid(&self) -> Option<(usize, usize)> {
        match self.data {
            PanelData::Field { nx, ny, .. } => Some((nx, ny)),
            PanelData::Points { .. } => None,
        }
    }
}

/// Capture every placed domain at the simulation's current time.
///
/// Asks each domain what it is rather than being told: a field is sampled over its
/// [`Extent`], a body set is read through [`Bodies`](dualis_core::Bodies), and every domain
/// contributes its [`readings`](dualis_core::Domain::readings) whether or not it drew anything.
///
/// A domain with no placement is still read for scalars. A domain that is a field and has no
/// extent is **not** drawn, because nobody said how big it is — and a guess would be a picture of
/// a region the caller never chose.
pub fn capture(sim: &Simulation, placed: &BTreeMap<String, Placement>) -> Frame {
    let t = sim.time();
    let mut panels = Vec::new();
    let mut readings = Vec::new();

    for domain in sim.domains() {
        let name = domain.name().to_string();
        readings.extend(domain.readings());
        let placement = placed.get(&name).copied().unwrap_or_default();

        if let (Some(field), Some(extent)) = (domain.as_field(), placement.extent) {
            panels.push(sample(&name, field, extent, placement.pose, t));
        } else if let Some(bodies) = domain.as_bodies() {
            panels.push(points(&name, bodies, placement.pose));
        }
    }
    Frame {
        time_s: t.to_si(),
        panels,
        readings,
    }
}

/// Sample a field over the extent it was placed with.
fn sample(name: &str, field: &dyn ScalarField, extent: Extent, pose: Pose, t: Time) -> Panel {
    let (nx, ny) = extent.samples;
    let (lo, hi) = (extent.min.to_si(), extent.max.to_si());
    let mut values = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let u = if nx > 1 {
                i as f64 / (nx - 1) as f64
            } else {
                0.5
            };
            let v = if ny > 1 {
                j as f64 / (ny - 1) as f64
            } else {
                0.5
            };
            let local = LengthVec::from_si(lo + (hi - lo) * glam::DVec3::new(u, v, 0.0));
            // Sampled in the domain's own coordinates. The pose is where the *result* goes, not
            // where the question is asked — a field does not know it has been placed.
            let _ = pose;
            values.push(field.at(local, t));
        }
    }
    Panel {
        name: name.to_string(),
        unit: field.unit(),
        data: PanelData::Field { nx, ny, values },
    }
}

/// Read a body set, in world coordinates.
fn points(name: &str, bodies: &dyn dualis_core::Bodies, pose: Pose) -> Panel {
    let n = bodies.count();
    let mut positions = Vec::with_capacity(n);
    let mut values = Vec::with_capacity(n);
    for i in 0..n {
        let p = pose.point_to_world(bodies.position(i)).to_si();
        positions.push([p.x, p.y, p.z]);
        values.push(bodies.value(i));
    }
    let (bounds, boxed) = match bodies.cell() {
        Some((lo, hi)) => {
            let (lo, hi) = (
                pose.point_to_world(lo).to_si(),
                pose.point_to_world(hi).to_si(),
            );
            (
                [
                    lo.x.min(hi.x),
                    lo.y.min(hi.y),
                    lo.z.min(hi.z),
                    lo.x.max(hi.x),
                    lo.y.max(hi.y),
                    lo.z.max(hi.z),
                ],
                true,
            )
        }
        None => {
            // Measured, and widened over the run later. Nothing physical sits at this edge.
            let r = positions
                .iter()
                .flat_map(|p| p.iter())
                .fold(0.0f64, |m, v| m.max(v.abs()))
                * 1.2;
            let r = if r > 0.0 { r } else { 1.0 };
            ([-r, -r, -r, r, r, r], false)
        }
    };
    Panel {
        name: name.to_string(),
        unit: bodies.value_unit(),
        data: PanelData::Points {
            positions,
            values,
            bounds,
            boxed,
        },
    }
}

/// Give every body panel one framing for the whole run.
///
/// [`capture`] sees a frame at a time and cannot see the future, so a panel without a real wall
/// comes back framed to *that* frame — and a body crossing the picture would look still while the
/// picture moved. Call this once on a finished run.
///
/// Panels with a real wall are left alone: a periodic cell is a boundary condition and does not
/// grow because a run is longer.
pub fn settle_framing(frames: &mut [Frame]) {
    let names: Vec<String> = frames
        .first()
        .map(|f| {
            f.panels
                .iter()
                .filter(|p| matches!(p.data, PanelData::Points { boxed: false, .. }))
                .map(|p| p.name.clone())
                .collect()
        })
        .unwrap_or_default();

    for name in names {
        let mut widest = [0.0f64; 6];
        for frame in frames.iter() {
            if let Some(Panel {
                data: PanelData::Points { bounds, .. },
                ..
            }) = frame.panels.iter().find(|p| p.name == name)
            {
                for k in 0..3 {
                    widest[k] = widest[k].min(bounds[k]);
                    widest[k + 3] = widest[k + 3].max(bounds[k + 3]);
                }
            }
        }
        for frame in frames.iter_mut() {
            if let Some(Panel {
                data: PanelData::Points { bounds, .. },
                ..
            }) = frame.panels.iter_mut().find(|p| p.name == name)
            {
                *bounds = widest;
            }
        }
    }
}
