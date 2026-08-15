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

    let mut out = Checked {
        error: World::build(scene.clone()).err(),
        summary: Some(format!(
            "{}: {} domain(s), {:.3} s in {} frames",
            scene.title,
            scene.domains.len(),
            scene.duration_s,
            scene.frames
        )),
        ..Checked::default()
    };

    for spec in &scene.domains {
        let placement = spec.placement();
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
}
