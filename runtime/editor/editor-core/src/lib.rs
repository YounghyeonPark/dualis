//! The GUI-free half of the scene editor, for the same reason `viewer-core` is the GUI-free
//! half of the viewer: everything that could be got wrong the same way twice lives where a
//! test can reach it, and the shell is left holding a text box, a canvas and an event loop.
//!
//! What lives here:
//!
//! - **Checking** — the same two steps `dualis-world --check` runs, parse then build, with the
//!   parse error carried as `line:column` because that is what an editor puts a squiggle
//!   under.
//! - **Placement geometry** — every placed extent as eight posed corners, ready to wireframe,
//!   with the union bounds a camera fits to. The corners go through [`Pose::point_to_world`]
//!   even though no scene can state a pose yet, so the day the format grows one this crate
//!   does not need to learn about it.
//! - **Running and verifying** — thin passes over [`World::run`] and
//!   [`dualis_world::verify::verify`], returning the run's JSON (which `viewer-core` reads)
//!   and the battery's rendered report.
//!
//! # The two halves, and which one this is
//!
//! ARCHITECTURE.md's platform rules split an editor into an authoring half, which is the
//! composition root and may name domains, and an inspection half, which must dispatch on the
//! shape of the data so that an eleventh physics costs no editor edit. This crate is the
//! authoring half's machinery: it consumes `Scene` and `DomainSpec` through `dualis-world`'s
//! public API — `DomainSpec::placement()` is where domain knowledge already legitimately
//! lives — and hands the shell *shapes*: boxes, points, paths, readings. The shell's painting
//! code never sees a domain name except to label things, which is what keeps the viewport
//! open to domains that do not exist yet.

#![deny(missing_docs)]

use dualis::units::LengthVec;
use dualis_world::{Scene, World};

/// One placed extent, as the eight corners of its box in world coordinates, in metres.
///
/// Corner order is the binary one — bit 0 is x, bit 1 is y, bit 2 is z, low corner for a
/// clear bit — which is the same order `Camera::fit` walks and the order [`EDGES`] indexes.
#[derive(Clone, Debug, PartialEq)]
pub struct PlacedBox {
    /// The domain this box belongs to, for a label beside the wireframe.
    pub name: String,
    /// Eight corners, `[x, y, z]` each, in metres.
    pub corners: [[f64; 3]; 8],
}

/// The twelve edges of a box, as index pairs into [`PlacedBox::corners`].
///
/// Written out rather than generated because twelve constant pairs cannot be wrong quietly,
/// and a loop that generates them can — an edge between corners differing in two bits draws a
/// diagonal, and a wireframe with a diagonal in it still looks like geometry.
pub const EDGES: [(usize, usize); 12] = [
    (0, 1),
    (2, 3),
    (4, 5),
    (6, 7), // along x
    (0, 2),
    (1, 3),
    (4, 6),
    (5, 7), // along y
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7), // along z
];

/// What checking the editor's text produced. Every field is present on both success and
/// failure, because an editor redraws from whatever state it has: an error should not blank
/// the viewport that was drawn from the last good parse.
#[derive(Clone, Debug, Default)]
pub struct Checked {
    /// The parse or build error, `None` when the scene is runnable. A parse error leads with
    /// `line:column`, which is what an editor puts a squiggle under.
    pub error: Option<String>,
    /// One line of what the scene is — title, domains, duration, frames — for the header.
    pub summary: Option<String>,
    /// What the build dismissed and why — a stated condition a domain correctly ignored,
    /// with the measurement that earns it. Shown under the summary, because a dismissal
    /// nobody can see is the silence it exists to replace.
    pub notes: Vec<String>,
    /// Every placed extent, posed into world coordinates.
    pub boxes: Vec<PlacedBox>,
    /// The union of every box, `[x0, y0, z0, x1, y1, z1]`, for the camera to fit. `None` when
    /// nothing in the scene has geometry — which the shell must say rather than draw an empty
    /// viewport that reads as "framed on nothing".
    pub bounds: Option<[f64; 6]>,
}

/// Parse and build the text, and lay out its geometry.
///
/// The same two steps `dualis-world --check` runs, in the same order, so the editor and the
/// CLI cannot disagree about what a valid scene is. Geometry is laid out from the *parsed*
/// scene even when the build fails, because "the beam's face count disagrees with the bar's"
/// is exactly when a person wants to be looking at the boxes.
pub fn check(text: &str) -> Checked {
    let scene: Scene = match serde_json::from_str(text) {
        Ok(s) => s,
        Err(e) => {
            return Checked {
                error: Some(format!("{}:{}: {e}", e.line(), e.column())),
                ..Checked::default()
            }
        }
    };

    let (error, notes) = match World::build(scene.clone()) {
        Ok(world) => (None, world.notes().to_vec()),
        Err(e) => (Some(e), Vec::new()),
    };
    let mut out = Checked {
        error,
        notes,
        summary: Some(format!(
            "{}: {} domain(s), {:.3} s in {} frames",
            scene.title,
            scene.domains.len(),
            scene.duration_s,
            scene.frames
        )),
        ..Checked::default()
    };

    // The scene's own placements, so a pose the file states moves the box on screen — and so
    // this shell and `World` cannot disagree about where anything is.
    let placed = scene.placements();
    for spec in &scene.domains {
        let placement = placed
            .get(spec.name())
            .copied()
            .unwrap_or_else(|| spec.placement());
        let Some(extent) = placement.extent else {
            continue;
        };
        let (lo, hi) = (extent.min.to_si(), extent.max.to_si());
        let mut corners = [[0.0; 3]; 8];
        for (i, corner) in corners.iter_mut().enumerate() {
            let local = LengthVec::m(
                if i & 1 == 0 { lo.x } else { hi.x },
                if i & 2 == 0 { lo.y } else { hi.y },
                if i & 4 == 0 { lo.z } else { hi.z },
            );
            let world = placement.pose.point_to_world(local).to_si();
            *corner = [world.x, world.y, world.z];
        }
        out.boxes.push(PlacedBox {
            name: spec.name().to_string(),
            corners,
        });
    }

    let mut bounds = [
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    for b in &out.boxes {
        for c in &b.corners {
            for a in 0..3 {
                bounds[a] = bounds[a].min(c[a]);
                bounds[a + 3] = bounds[a + 3].max(c[a]);
            }
        }
    }
    if bounds[0].is_finite() {
        out.bounds = Some(bounds);
    }
    out
}

/// Run the scene and return the run as JSON — the same bytes `dualis-world scene.json out.json`
/// writes, which is the format `viewer-core` reads. A violation is the error, worded by the
/// kernel.
pub fn run(text: &str) -> Result<String, String> {
    let scene: Scene =
        serde_json::from_str(text).map_err(|e| format!("{}:{}: {e}", e.line(), e.column()))?;
    let title = scene.title.clone();
    let mut world = World::build(scene)?;
    let frames = world.run().map_err(|v| {
        format!(
            "the audit stopped the run at t = {:.4} s: {v}",
            world.time().to_si()
        )
    })?;
    Ok(dualis::view::to_json(&title, &frames))
}

/// How a streamed run ended, when it did not fail.
///
/// Its own type rather than a bare `Ok(())`, because a stopped run and a finished one must
/// not be confusable: a partial run that reads as complete is a picture of something that did
/// not happen, which is this workspace's oldest failure shape wearing a scrub slider.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RunEnd {
    /// Every frame the scene asked for was captured.
    Finished,
    /// The stop flag was raised between frames; what was emitted is a prefix, and the caller
    /// must say so wherever the frames are shown.
    Stopped,
}

/// Run the scene, emitting the run-so-far as JSON after every captured frame.
///
/// This is what makes a run **watchable while it happens**: each `emit` payload is a complete,
/// readable run — `viewer-core` parses every one — containing the frames captured so far, and
/// the last payload is byte-identical to what [`run`] returns for the same text, which the
/// tests pin. Intermediate payloads are the run *unsettled*: `settle_framing` runs once at the
/// end, exactly as [`World::run`] does, so a shell scrubbing mid-run sees each frame's own
/// framing and the picture settles when the run does — the honest rendering of a run that is
/// not finished yet.
///
/// `stop` is read between frames. A violation still emits nothing extra: the frames already
/// emitted stand, and the error carries the kernel's own words, so a shell can leave the
/// partial run on screen *beside* the reason it ended — which is precisely the view somebody
/// debugging a violation wants.
pub fn run_streaming(
    text: &str,
    stop: &std::sync::atomic::AtomicBool,
    mut emit: impl FnMut(String),
) -> Result<RunEnd, String> {
    use std::sync::atomic::Ordering;

    let scene: Scene =
        serde_json::from_str(text).map_err(|e| format!("{}:{}: {e}", e.line(), e.column()))?;
    let title = scene.title.clone();
    let mut world = World::build(scene.clone())?;
    let dt = dualis::units::Time::from_si(scene.duration_s / scene.frames as f64);
    let placed = world.placements();

    let mut frames = vec![dualis::scene::capture(world.simulation(), &placed)];
    emit(dualis::view::to_json(&title, &frames));
    for _ in 0..scene.frames {
        if stop.load(Ordering::Relaxed) {
            return Ok(RunEnd::Stopped);
        }
        world.advance(dt).map_err(|v| {
            format!(
                "the audit stopped the run at t = {:.4} s: {v}",
                world.time().to_si()
            )
        })?;
        frames.push(dualis::scene::capture(world.simulation(), &placed));
        emit(dualis::view::to_json(&title, &frames));
    }
    dualis::scene::settle_framing(&mut frames);
    emit(dualis::view::to_json(&title, &frames));
    Ok(RunEnd::Finished)
}

/// Run the verification battery and return its rendered report, with the findings count so
/// the shell can be as loud as the CLI's exit code.
pub fn verify(text: &str, deep: bool) -> Result<(String, usize), String> {
    let scene: Scene =
        serde_json::from_str(text).map_err(|e| format!("{}:{}: {e}", e.line(), e.column()))?;
    let battery = dualis_world::verify::verify(&scene, deep)?;
    Ok((battery.render(), battery.findings.len()))
}

#[cfg(test)]
mod tests {
    use super::*;

    const ROOM: &str = r#"{
  "title": "a room",
  "duration_s": 0.005,
  "frames": 2,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.0, "height_m": 2.0,
      "cells_across": 21 }
  ]
}"#;

    /// The check mirrors `--check`: a valid scene has no error, a summary, and geometry.
    #[test]
    fn a_valid_scene_checks_clean_and_lays_out() {
        let c = check(ROOM);
        assert!(c.error.is_none(), "{:?}", c.error);
        assert_eq!(c.boxes.len(), 1);
        let b = c.bounds.expect("a room has geometry");
        assert_eq!(b[3] - b[0], 4.0, "the box spans the stated width");
        assert_eq!(b[4] - b[1], 2.0);
    }

    /// A parse error carries `line:column`, which is what an editor puts a squiggle under —
    /// and the failure must not blank the rest of the state, it IS the state.
    #[test]
    fn a_parse_error_names_its_line_and_column() {
        let c = check("{ \"title\": ");
        let e = c.error.expect("truncated JSON is an error");
        assert!(
            e.starts_with("1:"),
            "no position an editor could squiggle: {e}"
        );
        assert!(c.boxes.is_empty() && c.bounds.is_none());
    }

    /// A build error keeps the parsed geometry, because a cross-domain disagreement is
    /// exactly when a person wants to be looking at the boxes.
    #[test]
    fn a_build_error_keeps_the_geometry() {
        let two_names = ROOM.replace(
            "\"domains\": [",
            r#""domains": [
    { "kind": "room", "name": "room", "width_m": 1.0, "height_m": 1.0, "cells_across": 5 },"#,
        );
        let c = check(&two_names);
        assert!(
            c.error
                .as_deref()
                .is_some_and(|e| e.contains("both called")),
            "{:?}",
            c.error
        );
        assert_eq!(c.boxes.len(), 2, "the boxes survive the refusal");
    }

    /// Every edge joins corners differing in exactly one bit — an edge, not a diagonal.
    #[test]
    fn the_edges_are_edges_and_cover_the_box() {
        for (a, b) in EDGES {
            assert_eq!((a ^ b).count_ones(), 1, "{a}-{b} is a diagonal");
        }
        // Each corner is touched by exactly three edges, as a cube's corners are.
        for corner in 0..8 {
            let touching = EDGES
                .iter()
                .filter(|(a, b)| *a == corner || *b == corner)
                .count();
            assert_eq!(touching, 3);
        }
    }

    /// The run's JSON is the wire format: `viewer-core` — which never links `dualis` — reads
    /// it back whole. This is the editor standing on the same contract the viewer proved.
    #[test]
    fn a_run_round_trips_through_the_wire_format() {
        let json = run(ROOM).expect("the room runs");
        let run = viewer_core::Run::from_json(&json).expect("the viewer's reader accepts it");
        assert_eq!(run.frames.len(), 3, "two frames plus the initial capture");
        assert!(run.frames[0].panels.iter().any(|p| p.name() == "room"));
    }

    /// The verify pass returns the same report the CLI prints, findings counted.
    #[test]
    fn verify_reports_and_counts_findings() {
        let (report, findings) = verify(ROOM, false).expect("the battery runs");
        assert_eq!(findings, 0, "{report}");
        assert!(report.contains("determinism     two runs, identical bytes"));
    }

    /// Streaming emits one readable run per capture, and its last payload is byte-identical
    /// to [`run`]'s — the stream lands exactly where the batch does, so watching a run happen
    /// and reading it afterwards are the same run.
    #[test]
    fn a_streamed_run_is_watchable_and_lands_on_the_batch_answer() {
        let stop = std::sync::atomic::AtomicBool::new(false);
        let mut payloads = Vec::new();
        let end = run_streaming(ROOM, &stop, |j| payloads.push(j)).expect("the room runs");
        assert_eq!(end, RunEnd::Finished);
        // The initial capture, one per frame, and the settled final.
        assert_eq!(payloads.len(), 4);
        for p in &payloads {
            viewer_core::Run::from_json(p).expect("every payload is a whole, readable run");
        }
        assert_eq!(
            payloads.last().unwrap(),
            &run(ROOM).unwrap(),
            "the stream must land on the batch answer to the byte"
        );
    }

    /// The stop flag ends a run between frames, and the ending says Stopped — a prefix must
    /// never be confusable with the run.
    #[test]
    fn a_stopped_run_says_so() {
        use std::sync::atomic::{AtomicBool, Ordering};
        let stop = AtomicBool::new(false);
        let mut emitted = 0;
        let end = run_streaming(ROOM, &stop, |_| {
            emitted += 1;
            stop.store(true, Ordering::Relaxed);
        })
        .expect("stopping is not a violation");
        assert_eq!(end, RunEnd::Stopped);
        assert_eq!(emitted, 1, "stopped after the first capture");
    }
}

/// One drawn cell of a field: where it is in the world, and what colour it is.
///
/// World space and no camera, deliberately. Two shells draw these now — the native window and
/// the browser — and the projection is `viewer-core`'s in both, so what lives here is the half
/// that could be got wrong the same way twice: **where the samples sit inside the placed box,
/// and what colour a value is.** The camera was already shared for that reason; this is the
/// same argument one layer along.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Splat {
    /// Where, in world metres.
    pub at: [f64; 3],
    /// Colour and coverage, straight sRGB with an alpha.
    pub rgba: [u8; 4],
}

/// What a field's splats are, and the sentence a reader needs beside them.
#[derive(Clone, Debug)]
pub struct Splatted {
    /// The cells to draw, in grid order. A shell projects, sorts by depth and paints.
    pub splats: Vec<Splat>,
    /// Cells skipped per axis. One means every cell is drawn.
    pub stride: usize,
    /// Whether the colours are Planck's or the conventional ramp — see [`field_splats`].
    pub physical: bool,
    /// The note the canvas must carry, saying which of the two a reader is looking at.
    pub note: &'static str,
}

/// The most cells one field may draw in a frame.
///
/// A 100³ field is a million splats and a painter that stops painting. Past this the field is
/// subsampled with a stride and the note says so — a silently decimated picture is a picture of
/// a coarser simulation than the one that ran.
pub const MAX_SPLATS: usize = 8000;

/// Turn a field panel into cells to draw, coloured by physics where physics gives a colour.
///
/// # Two colourings, and the physics decides which
///
/// A temperature field is drawn in the colour a body at that temperature **actually is** —
/// Planck's law through the CIE matching functions, from [`dualis::view::colour`] — whenever
/// anything in it is hot enough to emit visible light. That is not a palette: a melting block
/// glows the orange a melting block glows, and nothing here picked it.
///
/// Below that, physics gives no colour. A body at 300 K emits nothing visible and this
/// workspace holds no visible reflectance for it to have instead, so the field falls back to a
/// conventional ramp — which says *more* and *less* and does not pretend to say *looks like* —
/// and [`Splatted::note`] states which, because a false colour mistaken for a real one is a
/// wrong answer that looks right.
///
/// The dispatch is on the panel's **unit**, which is data rather than a domain name, so a
/// physics added next year is coloured correctly with no edit here.
pub fn field_splats(
    corners: &[[f64; 3]; 8],
    counts: (usize, usize, usize),
    values: &[f64],
    unit: &str,
    scale: Option<(f64, f64)>,
) -> Splatted {
    let (nx, ny, nz) = counts;
    if nx == 0 || ny == 0 || nz == 0 || values.len() < nx * ny * nz {
        return Splatted {
            splats: Vec::new(),
            stride: 1,
            physical: false,
            note: "field: the panel's values do not fill its grid — not drawn",
        };
    }

    // The box's axes. Corner 0 is the low one and bits 0, 1, 2 step one axis each, which is the
    // order [`EDGES`] is written against.
    let o = corners[0];
    let axis = |c: [f64; 3]| [c[0] - o[0], c[1] - o[1], c[2] - o[2]];
    let (ax, ay, az) = (axis(corners[1]), axis(corners[2]), axis(corners[4]));

    let to_kelvin = |v: f64| match unit {
        "K" => Some(v),
        "C" => Some(v + 273.15),
        _ => None,
    };
    let hottest = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let peak_glow = to_kelvin(hottest).map_or(0.0, dualis::view::glow_fraction);
    let physical = peak_glow > 1e-6;

    let total = nx * ny * nz;
    let stride = if total > MAX_SPLATS {
        ((total as f64 / MAX_SPLATS as f64).cbrt().ceil() as usize).max(1)
    } else {
        1
    };
    let frac = |i: usize, n: usize| {
        if n > 1 {
            i as f64 / (n - 1) as f64
        } else {
            0.5
        }
    };

    let mut splats = Vec::new();
    for k in (0..nz).step_by(stride) {
        for j in (0..ny).step_by(stride) {
            for i in (0..nx).step_by(stride) {
                let v = values[i + nx * (j + ny * k)];
                if !v.is_finite() {
                    continue;
                }
                let (u, w, t) = (frac(i, nx), frac(j, ny), frac(k, nz));
                let at = [
                    o[0] + ax[0] * u + ay[0] * w + az[0] * t,
                    o[1] + ax[1] * u + ay[1] * w + az[1] * t,
                    o[2] + ax[2] * u + ay[2] * w + az[2] * t,
                ];
                // One scale across the whole run, never per frame, for the reason
                // `viewer-core` states.
                let s = match scale {
                    Some((lo, hi)) if hi > lo => ((v - lo) / (hi - lo)).clamp(0.0, 1.0),
                    _ => 0.5,
                };
                let rgba = if physical {
                    let kelvin = to_kelvin(v).unwrap_or(0.0);
                    let [r, g, b] = dualis::view::blackbody_srgb(kelvin);
                    // Brightness is the glow relative to this field's own hottest cell, so a
                    // cool corner of a glowing block is dark rather than merely bluer — which
                    // is what a photograph of it looks like.
                    let rel = (dualis::view::glow_fraction(kelvin) / peak_glow).clamp(0.0, 1.0);
                    [r, g, b, ((rel.sqrt() * 235.0) as u8).max(6)]
                } else {
                    let (r, g, b) = ramp(s);
                    [r, g, b, (30.0 + 200.0 * s * s) as u8]
                };
                splats.push(Splat { at, rgba });
            }
        }
    }

    let note = match (physical, stride) {
        (true, 1) => "field: colour is Planck's, not a palette",
        (true, _) => "field: Planck colour, subsampled — see the report for every cell",
        (false, 1) => "field: false colour — nothing here is hot enough to glow",
        (false, _) => "field: false colour, subsampled — see the report for every cell",
    };
    Splatted {
        splats,
        stride,
        physical,
        note,
    }
}

/// The conventional ramp, for a quantity physics gives no colour to: cool blue through to warm.
///
/// A designer's choice and stated as one. It says *more* and *less*, which for a pressure
/// swinging about zero is the whole truth and for a temperature is a stand-in.
fn ramp(t: f64) -> (u8, u8, u8) {
    let t = t.clamp(0.0, 1.0);
    (
        (60.0 + 195.0 * t) as u8,
        (90.0 + 40.0 * (1.0 - (2.0 * t - 1.0).abs())) as u8,
        (230.0 - 170.0 * t) as u8,
    )
}
