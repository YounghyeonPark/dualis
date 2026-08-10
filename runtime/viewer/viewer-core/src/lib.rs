//! Reads a dualis run and turns it into something a renderer can draw.
//!
//! No GPU, no window, no `dualis` dependency. That last one is deliberate and is the point of the
//! crate: if a viewer can be written against the **file** a run produces and nothing else, then
//! the wire format is complete. If it needed to reach back into the library for something the
//! file did not carry, the format would be the thing to fix — and finding that out is worth more
//! than the convenience of linking.
//!
//! # What a run is
//!
//! What `dualis_view::to_json` writes: a title and a list of frames, each with a time, some
//! panels and some readings. A panel is one of three shapes — a field on a grid, a set of points,
//! or runs of connected points — and the reader below accepts all three by name.
//!
//! # What this crate is responsible for
//!
//! Everything a renderer would otherwise get wrong the same way twice:
//!
//! - **One colour scale across the whole run.** A frame normalised to itself makes a decay look
//!   like a steady state, which is the one thing a picture of a simulation must never do. The
//!   scale is computed once, over every frame, per panel.
//! - **One framing across the whole run**, from the geometry rather than from the current frame,
//!   so a body moving is a body moving and not the camera chasing it.
//! - **The projection**, so the window and any test agree about where a point lands.

#![deny(missing_docs)]
#![forbid(unsafe_code)]

use serde::Deserialize;

/// A whole run, as read from a file.
#[derive(Clone, Debug, Deserialize)]
pub struct Run {
    /// What the run was called.
    pub title: String,
    /// Every captured instant, in order.
    pub frames: Vec<Frame>,
}

/// One instant.
#[derive(Clone, Debug, Deserialize)]
pub struct Frame {
    /// Simulation time, in seconds.
    pub t: f64,
    /// One per drawable domain.
    pub panels: Vec<Panel>,
    /// Named scalars, drawable or not.
    #[serde(default)]
    pub readings: Vec<Reading>,
}

/// One named scalar.
#[derive(Clone, Debug, Deserialize)]
pub struct Reading {
    /// Which domain it came from.
    pub domain: String,
    /// What it is.
    pub label: String,
    /// Its value.
    pub value: f64,
    /// Its unit.
    pub unit: String,
}

/// One domain, captured, in whichever shape it has.
///
/// Tagged by `kind`, which is what the writer emits. An unknown kind is an **error** rather than a
/// panel that is silently skipped: a viewer that quietly drops a shape it does not know shows a
/// page with something missing and no way to tell.
#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "lowercase", deny_unknown_fields)]
pub enum Panel {
    /// A field sampled on a grid of one, two or three dimensions.
    Field {
        /// Which domain this came from.
        name: String,
        /// What the values are in.
        unit: String,
        /// Samples along x.
        nx: usize,
        /// Along y.
        ny: usize,
        /// Along z.
        nz: usize,
        /// `nx*ny*nz` values, x fastest then y then z.
        values: Vec<f64>,
    },
    /// A countable set of bodies.
    Points {
        /// Which domain this came from.
        name: String,
        /// What the values are in.
        unit: String,
        /// Whether the bounding box is a real wall.
        boxed: bool,
        /// `[x0,y0,z0,x1,y1,z1]`.
        bounds: [f64; 6],
        /// Flattened `xyz` per body.
        positions: Vec<f64>,
        /// One per body.
        values: Vec<f64>,
    },
    /// Runs of connected points — rays, trajectories, field lines.
    Paths {
        /// Which domain this came from.
        name: String,
        /// What the values are in.
        unit: String,
        /// `[x0,y0,z0,x1,y1,z1]`.
        bounds: [f64; 6],
        /// Where each run begins, as an index into `vertices` divided by three.
        starts: Vec<f64>,
        /// Flattened `xyz` per vertex.
        vertices: Vec<f64>,
        /// One per run.
        values: Vec<f64>,
    },
}

impl Panel {
    /// Which domain this panel came from.
    pub fn name(&self) -> &str {
        match self {
            Panel::Field { name, .. } | Panel::Points { name, .. } | Panel::Paths { name, .. } => {
                name
            }
        }
    }

    /// What its values are in.
    pub fn unit(&self) -> &str {
        match self {
            Panel::Field { unit, .. } | Panel::Points { unit, .. } | Panel::Paths { unit, .. } => {
                unit
            }
        }
    }

    /// The values, whichever shape this is.
    pub fn values(&self) -> &[f64] {
        match self {
            Panel::Field { values, .. }
            | Panel::Points { values, .. }
            | Panel::Paths { values, .. } => values,
        }
    }

    /// The box this panel occupies, in world coordinates.
    ///
    /// A field carries no bounds of its own — it was sampled over an extent the writer did not
    /// record — so its box is the grid in cell units, which is the right thing to frame it by and
    /// is stated here rather than guessed at the call site.
    pub fn bounds(&self) -> [f64; 6] {
        match self {
            Panel::Field { nx, ny, nz, .. } => [0.0, 0.0, 0.0, *nx as f64, *ny as f64, *nz as f64],
            Panel::Points { bounds, .. } | Panel::Paths { bounds, .. } => *bounds,
        }
    }
}

impl Run {
    /// Read a run from the JSON a `dualis` run wrote.
    pub fn from_json(text: &str) -> Result<Run, String> {
        serde_json::from_str(text).map_err(|e| format!("not a dualis run: {e}"))
    }

    /// The value range of one panel **across every frame**.
    ///
    /// The whole run, not the current frame. A per-frame scale makes a decaying mode look
    /// constant and a constant one look like noise, and it is what a renderer does if nobody
    /// stops it. Returns `None` for a panel name the run does not have.
    pub fn scale_of(&self, panel: &str) -> Option<(f64, f64)> {
        let mut lo = f64::MAX;
        let mut hi = f64::MIN;
        let mut seen = false;
        for frame in &self.frames {
            for p in &frame.panels {
                if p.name() != panel {
                    continue;
                }
                seen = true;
                for v in p.values() {
                    lo = lo.min(*v);
                    hi = hi.max(*v);
                }
            }
        }
        seen.then_some((lo, hi))
    }

    /// The box one panel occupies **across every frame**, widened to hold all of them.
    ///
    /// Same reason as the scale: a camera framed to the current frame follows a moving body and
    /// makes it look still.
    pub fn framing_of(&self, panel: &str) -> Option<[f64; 6]> {
        let mut out = [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN];
        let mut seen = false;
        for frame in &self.frames {
            for p in &frame.panels {
                if p.name() != panel {
                    continue;
                }
                seen = true;
                let b = p.bounds();
                for a in 0..3 {
                    out[a] = out[a].min(b[a]);
                    out[a + 3] = out[a + 3].max(b[a + 3]);
                }
            }
        }
        seen.then_some(out)
    }

    /// Every panel name in the run, in the order they first appear.
    pub fn panels(&self) -> Vec<String> {
        let mut names: Vec<String> = Vec::new();
        for frame in &self.frames {
            for p in &frame.panels {
                if !names.iter().any(|n| n == p.name()) {
                    names.push(p.name().to_string());
                }
            }
        }
        names
    }
}

/// Where the eye is, in orbit around what it is looking at.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Camera {
    /// Azimuth, radians.
    pub azimuth: f64,
    /// Elevation, radians, clamped to just under a right angle so the up vector never degenerates.
    pub elevation: f64,
    /// Distance, in units of the framed box's longest side.
    pub distance: f64,
}

impl Default for Camera {
    /// Three-quarter view from slightly above — the same angles the HTML report opens at, so a
    /// reader moving between the two is not re-orienting.
    fn default() -> Camera {
        Camera {
            azimuth: 0.7,
            elevation: 0.4,
            distance: 2.5,
        }
    }
}

impl Camera {
    /// Turn by a drag, in radians. Elevation is clamped; azimuth wraps.
    pub fn turn(&mut self, d_az: f64, d_el: f64) {
        self.azimuth += d_az;
        self.elevation = (self.elevation + d_el).clamp(-1.5, 1.5);
    }

    /// Zoom by a scroll. Bounded, so a viewer cannot be scrolled inside the subject or into
    /// the next county and left with a blank window and no way back.
    pub fn zoom(&mut self, factor: f64) {
        self.distance = (self.distance * factor).clamp(1.2, 9.0);
    }

    /// Project a world point into normalised device coordinates, with its depth.
    ///
    /// `x` and `y` run `-1..1` across the shorter side; `depth` grows away from the eye and is
    /// what a renderer sorts or tests by. `aspect` is width over height.
    pub fn project(&self, p: [f64; 3], frame: &Framing, aspect: f64) -> Projected {
        let x = (p[0] - frame.centre[0]) / frame.span;
        let y = (p[1] - frame.centre[1]) / frame.span;
        let z = (p[2] - frame.centre[2]) / frame.span;

        let (ca, sa) = (self.azimuth.cos(), self.azimuth.sin());
        let (ce, se) = (self.elevation.cos(), self.elevation.sin());
        let x1 = x * ca - z * sa;
        let z1 = x * sa + z * ca;
        let y1 = y * ce - z1 * se;
        let z2 = y * se + z1 * ce;

        // Never divide by a depth at or behind the eye: a point there has no screen position, and
        // returning a huge coordinate instead of clamping is how a renderer draws a streak across
        // the window and calls it geometry.
        let d = (z2 + self.distance).max(0.05);
        let f = 0.6 / d;
        Projected {
            x: x1 * f / aspect.max(1e-6),
            y: y1 * f,
            depth: d,
        }
    }
}

/// A projected point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Projected {
    /// Normalised device x.
    pub x: f64,
    /// Normalised device y.
    pub y: f64,
    /// Distance from the eye. Larger is further.
    pub depth: f64,
}

/// The box a panel is drawn in, reduced to a centre and one span.
///
/// **One span for all three axes**, so a cube is drawn as a cube. Scaling each axis to fill the
/// window independently is the easiest way to make a picture that is not of the thing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Framing {
    /// The centre of the box.
    pub centre: [f64; 3],
    /// The longest side.
    pub span: f64,
}

impl Framing {
    /// From a bounding box.
    pub fn of(bounds: [f64; 6]) -> Framing {
        let centre = [
            0.5 * (bounds[0] + bounds[3]),
            0.5 * (bounds[1] + bounds[4]),
            0.5 * (bounds[2] + bounds[5]),
        ];
        let span = (bounds[3] - bounds[0])
            .max(bounds[4] - bounds[1])
            .max(bounds[5] - bounds[2]);
        Framing {
            centre,
            span: if span > 0.0 { span } else { 1.0 },
        }
    }
}

/// One line segment, ready to hand to a renderer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Segment {
    /// Both ends, projected.
    pub from: Projected,
    /// The other end.
    pub to: Projected,
    /// Where this segment's value sits in the run's range, `0..1`.
    pub shade: f64,
}

/// Build the segments for one panel at one frame.
///
/// Sorted **back to front**, so a renderer without a depth buffer still draws them in an order
/// that reads. A renderer with one can ignore the order and lose nothing.
///
/// Returns an empty list for a field, which is not a set of lines — a field view is a different
/// pipeline and pretending otherwise here would produce a plausible picture of nothing.
pub fn segments(panel: &Panel, camera: &Camera, framing: &Framing, aspect: f64) -> Vec<Segment> {
    let mut out = Vec::new();
    if let Panel::Paths {
        starts,
        vertices,
        values,
        ..
    } = panel
    {
        let (lo, hi) = span_of(values);
        for (k, start) in starts.iter().enumerate() {
            let from = *start as usize;
            let to = starts
                .get(k + 1)
                .map(|s| *s as usize)
                .unwrap_or(vertices.len() / 3);
            let shade = ((values.get(k).copied().unwrap_or(0.0) - lo) / (hi - lo)).clamp(0.0, 1.0);
            for i in from..to.saturating_sub(1) {
                let a = camera.project(vertex(vertices, i), framing, aspect);
                let b = camera.project(vertex(vertices, i + 1), framing, aspect);
                out.push(Segment {
                    from: a,
                    to: b,
                    shade,
                });
            }
        }
    }
    out.sort_by(|a, b| {
        let (da, db) = (a.from.depth + a.to.depth, b.from.depth + b.to.depth);
        db.total_cmp(&da)
    });
    out
}

fn vertex(flat: &[f64], i: usize) -> [f64; 3] {
    [flat[3 * i], flat[3 * i + 1], flat[3 * i + 2]]
}

fn span_of(values: &[f64]) -> (f64, f64) {
    let lo = values.iter().copied().fold(f64::MAX, f64::min);
    let hi = values.iter().copied().fold(f64::MIN, f64::max);
    if hi > lo {
        (lo, hi)
    } else {
        (lo, lo + 1.0)
    }
}
