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
//! Knowing no physics does not make a layer complete, and this one was not. [`Extent`] described
//! a *plane* until a domain with a genuinely three-dimensional field arrived, and the gap did not
//! announce itself: the sampler built its position as `(u, v, 0)`, so a solid would have been
//! captured as its `z = 0` face and drawn as a perfectly plausible picture of a block. Nothing
//! here could have caught that, because every field this crate had ever been handed was flat. The
//! lesson is not about `Extent` — it is that a layer's assumptions are only visible from below.
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
    /// How many samples along x, y and z. A count of one collapses that axis: `(n, 1, 1)` is a
    /// line, `(n, m, 1)` a plane, and all three above one a volume.
    ///
    /// **Three, not two.** It was a pair until a domain with a genuinely three-dimensional field
    /// arrived, and a pair does not fail when handed one — it samples the `z = min` slice and
    /// returns a perfectly plausible picture of a solid. That is the shape of failure this
    /// workspace keeps finding: not a wrong answer, a *narrower* one, with nothing to say so.
    pub samples: (usize, usize, usize),
}

impl Extent {
    /// A box from two corners, sampled `nx` by `ny` by `nz`.
    pub fn new(min: LengthVec, max: LengthVec, nx: usize, ny: usize, nz: usize) -> Extent {
        Extent {
            min,
            max,
            samples: (nx.max(1), ny.max(1), nz.max(1)),
        }
    }

    /// How many samples this asks for in total.
    pub fn count(&self) -> usize {
        self.samples.0 * self.samples.1 * self.samples.2
    }

    /// How many of the three axes actually vary — 1, 2 or 3.
    ///
    /// What a view should dispatch on, rather than on the counts directly: a `(60, 1, 1)` line
    /// and a `(1, 60, 1)` line are the same kind of thing to draw and differ only in which way
    /// they point.
    pub fn dimensions(&self) -> usize {
        let (nx, ny, nz) = self.samples;
        [nx, ny, nz].iter().filter(|&&n| n > 1).count().max(1)
    }

    /// A line along x, for a domain with one dimension.
    pub fn line(length: Length, cells: usize) -> Extent {
        Extent::new(
            LengthVec::ZERO,
            LengthVec::from_si(glam::DVec3::new(length.to_si(), 0.0, 0.0)),
            cells.max(2),
            1,
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
            1,
        )
    }

    /// A box from the origin, sampled on a regular grid.
    ///
    /// The counts are given separately from the size because they are different questions: a
    /// block that is long and thin is still worth sampling evenly *in space*, which means more
    /// samples along the long axis and not larger ones.
    pub fn volume(size: LengthVec, nx: usize, ny: usize, nz: usize) -> Extent {
        Extent::new(LengthVec::ZERO, size, nx, ny, nz)
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
    /// A field, sampled on a grid of one, two or three dimensions.
    Field {
        /// Samples along x.
        nx: usize,
        /// Samples along y. One for a line.
        ny: usize,
        /// Samples along z. One for a line or a plane.
        ///
        /// A view that ignores this draws the `z = 0` slice of a solid and calls it the solid.
        /// It is a separate field rather than folded into `ny` for exactly that reason: an
        /// `nx * (ny*nz)` grid would still *render*, as a plane with the slices stacked into a
        /// stripe, and would look like a picture rather than like a mistake.
        nz: usize,
        /// `nx * ny * nz` values, x fastest, then y, then z.
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
    pub fn grid(&self) -> Option<(usize, usize, usize)> {
        match self.data {
            PanelData::Field { nx, ny, nz, .. } => Some((nx, ny, nz)),
            PanelData::Points { .. } => None,
        }
    }

    /// One z-slice of a field, as an `nx * ny` plane, or `None` if this is not a field or the
    /// slice is out of range.
    ///
    /// What a two-dimensional view should call rather than reading `values` directly. A view that
    /// takes the first `nx * ny` entries gets slice zero and no indication that there were
    /// others; this makes the choice explicit and countable.
    pub fn slice(&self, k: usize) -> Option<&[f64]> {
        match &self.data {
            PanelData::Field { nx, ny, nz, values } => {
                (k < *nz).then(|| &values[k * nx * ny..(k + 1) * nx * ny])
            }
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
    let (nx, ny, nz) = extent.samples;
    let (lo, hi) = (extent.min.to_si(), extent.max.to_si());
    // A single sample along an axis is taken at the **middle** of it rather than at `min`. For a
    // flat extent the two are the same point; for an extent with real thickness that was asked
    // for at one sample, the middle is the honest representative and the low face is a corner.
    let along = |i: usize, n: usize| {
        if n > 1 {
            i as f64 / (n - 1) as f64
        } else {
            0.5
        }
    };
    let mut values = Vec::with_capacity(extent.count());
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let f = glam::DVec3::new(along(i, nx), along(j, ny), along(k, nz));
                let local = LengthVec::from_si(lo + (hi - lo) * f);
                // Sampled in the domain's own coordinates. The pose is where the *result* goes,
                // not where the question is asked — a field does not know it has been placed.
                let _ = pose;
                values.push(field.at(local, t));
            }
        }
    }
    Panel {
        name: name.to_string(),
        unit: field.unit(),
        data: PanelData::Field { nx, ny, nz, values },
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
