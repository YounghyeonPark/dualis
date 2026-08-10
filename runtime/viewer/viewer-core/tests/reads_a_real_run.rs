//! The viewer reads what a run actually wrote, with no `dualis` dependency at all.
//!
//! That last clause is the claim under test. The fixtures in `tests/runs/` are genuine output —
//! `dualis-world` on two shipped scenes and the `optical_bench` example — trimmed to a couple of
//! frames each so a hundred kilobytes of generated JSON does not enter the tree. Nothing here
//! links the library, so if a shape could not be drawn from the file alone the fix would have to
//! be in the wire format, which is where it belongs.

use viewer_core::{segments, Camera, Framing, Panel, Projected, Run};

fn load(name: &str) -> Run {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/runs")
        .join(name);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    Run::from_json(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// **All three panel shapes parse out of real files.**
///
/// A field, a set of bodies and a set of paths, each from a different producer. A reader that
/// handles two of three is a reader that shows a blank card for the third.
#[test]
fn every_shape_a_run_can_write_is_read_back() {
    let block = load("block.json");
    let field = block.frames[0]
        .panels
        .iter()
        .find(|p| matches!(p, Panel::Field { .. }))
        .expect("a hot spot in a block is a field");
    let Panel::Field {
        nx, ny, nz, values, ..
    } = field
    else {
        unreachable!()
    };
    assert_eq!((*nx, *ny, *nz), (9, 9, 9));
    assert_eq!(values.len(), 9 * 9 * 9, "every sample, not one slice");
    assert_eq!(field.unit(), "K", "the field carries the unit it holds");

    let orbits = load("orbits.json");
    let points = orbits.frames[0]
        .panels
        .iter()
        .find(|p| matches!(p, Panel::Points { .. }))
        .expect("an orbit is bodies");
    let Panel::Points {
        positions,
        values,
        boxed,
        ..
    } = points
    else {
        unreachable!()
    };
    assert_eq!(positions.len(), values.len() * 3, "three coordinates each");
    assert!(!boxed, "an orbit has no wall to draw");

    let bench = load("bench.json");
    let paths = bench.frames[0]
        .panels
        .iter()
        .find(|p| matches!(p, Panel::Paths { .. }))
        .expect("a traced bench is paths");
    let Panel::Paths {
        starts,
        vertices,
        values,
        ..
    } = paths
    else {
        unreachable!()
    };
    assert_eq!(starts.len(), values.len(), "one value per run");
    assert!(
        starts.len() > 100,
        "61 rays at three fields: {}",
        starts.len()
    );
    assert_eq!(vertices.len() % 3, 0);
}

/// **The colour scale is taken over the whole run, not the current frame.**
///
/// The failure this prevents is the one every renderer makes for free: normalise each frame to
/// itself and a decaying mode looks constant, while a constant one looks like noise. Checked by
/// taking a run whose values genuinely move between frames and showing the run-wide range is
/// wider than either frame's.
#[test]
fn the_scale_spans_the_run_rather_than_a_frame() {
    let run = load("block.json");
    let (lo, hi) = run.scale_of("block").expect("the block is in the run");

    let frame_span = |k: usize| {
        let v = run.frames[k]
            .panels
            .iter()
            .find(|p| p.name() == "block")
            .expect("present in every frame")
            .values();
        (
            v.iter().copied().fold(f64::MAX, f64::min),
            v.iter().copied().fold(f64::MIN, f64::max),
        )
    };
    let (a_lo, a_hi) = frame_span(0);
    let (b_lo, b_hi) = frame_span(1);

    assert!(a_hi > b_hi, "the spot must cool between the two frames");
    assert_eq!(lo, a_lo.min(b_lo));
    assert_eq!(hi, a_hi.max(b_hi));
    assert!(
        hi - lo > b_hi - b_lo,
        "the run's range must be wider than the later frame's: {} against {}",
        hi - lo,
        b_hi - b_lo
    );
    assert!(run.scale_of("not a panel").is_none());
}

/// **The framing is taken over the whole run too, so a moving body does not move the camera.**
#[test]
fn the_framing_holds_still_while_the_bodies_move() {
    let run = load("orbits.json");
    let box_of = run.framing_of("system").expect("the orbit is in the run");
    let framing = Framing::of(box_of);

    for frame in &run.frames {
        for p in frame.panels.iter().filter(|p| p.name() == "system") {
            let b = p.bounds();
            for a in 0..3 {
                assert!(
                    b[a] >= box_of[a] - 1e-9,
                    "a frame reaches outside the run's box"
                );
                assert!(b[a + 3] <= box_of[a + 3] + 1e-9);
            }
        }
    }
    assert!(framing.span > 0.0);
    // One span for all three axes, so a cube is drawn as a cube.
    let sides = [
        box_of[3] - box_of[0],
        box_of[4] - box_of[1],
        box_of[5] - box_of[2],
    ];
    assert_eq!(framing.span, sides.iter().copied().fold(0.0, f64::max));
}

/// **The projection is a projection: on the axis, symmetric, and never divides by a depth at the
/// eye.**
///
/// Three properties rather than a golden number, because a golden number pins the arithmetic and
/// says nothing about whether it is a projection. The last one is the guard that matters: a point
/// at or behind the eye has no screen position, and returning a huge coordinate instead of
/// clamping is how a renderer draws a streak across the window and calls it geometry.
#[test]
fn the_projection_behaves_like_one() {
    let framing = Framing {
        centre: [0.0; 3],
        span: 1.0,
    };
    let camera = Camera {
        scale: 1.0,
        azimuth: 0.0,
        elevation: 0.0,
        distance: 2.5,
    };

    // The centre projects to the centre.
    let o = camera.project([0.0, 0.0, 0.0], &framing, 1.0);
    assert!(o.x.abs() < 1e-12 && o.y.abs() < 1e-12, "{o:?}");

    // Symmetric about it.
    let up = camera.project([0.0, 0.4, 0.0], &framing, 1.0);
    let down = camera.project([0.0, -0.4, 0.0], &framing, 1.0);
    assert!((up.y + down.y).abs() < 1e-12, "{up:?} {down:?}");

    // Nearer is bigger, which is what perspective is.
    let near = camera.project([0.0, 0.4, -1.0], &framing, 1.0);
    let far = camera.project([0.0, 0.4, 1.0], &framing, 1.0);
    assert!(near.y > far.y, "near {near:?} far {far:?}");
    assert!(near.depth < far.depth);

    // Aspect squeezes x and leaves y, so a wide window does not stretch the subject.
    let wide = camera.project([0.4, 0.4, 0.0], &framing, 2.0);
    let square = camera.project([0.4, 0.4, 0.0], &framing, 1.0);
    assert!((wide.x * 2.0 - square.x).abs() < 1e-12);
    assert!((wide.y - square.y).abs() < 1e-12);

    // A point behind the eye is clamped rather than turned into a huge coordinate.
    let behind = camera.project([0.0, 0.0, 100.0], &framing, 1.0);
    assert!(behind.depth >= 0.05, "{behind:?}");
    assert!(behind.x.is_finite() && behind.y.is_finite());
    assert!(
        behind.y.abs() < 100.0,
        "no streak across the window: {behind:?}"
    );
}

/// **Segments come out of a real bench, back to front, and only from paths.**
#[test]
fn a_traced_bench_becomes_line_segments() {
    let run = load("bench.json");
    let panel = run.frames[0]
        .panels
        .iter()
        .find(|p| matches!(p, Panel::Paths { .. }))
        .expect("paths");
    let framing = Framing::of(run.framing_of(panel.name()).unwrap());
    let camera = Camera::default();
    let span = run.scale_of(panel.name()).expect("the bench is in the run");
    let lines = segments(panel, &camera, &framing, 1.5, span);

    let Panel::Paths {
        starts, vertices, ..
    } = panel
    else {
        unreachable!()
    };
    // Each run of `n` vertices gives `n-1` segments, so the total is vertices minus runs.
    assert_eq!(lines.len(), vertices.len() / 3 - starts.len());
    assert!(lines.len() > 500, "only {} segments", lines.len());

    // Sorted back to front, so a renderer with no depth buffer still draws in an order that reads.
    for pair in lines.windows(2) {
        let (a, b) = (
            pair[0].from.depth + pair[0].to.depth,
            pair[1].from.depth + pair[1].to.depth,
        );
        assert!(a >= b - 1e-12, "not depth sorted: {a} then {b}");
    }
    // The shades span the run's own range rather than sitting on one value.
    let lo = lines.iter().map(|s| s.shade).fold(f64::MAX, f64::min);
    let hi = lines.iter().map(|s| s.shade).fold(f64::MIN, f64::max);
    assert!(
        (lo - 0.0).abs() < 1e-12 && (hi - 1.0).abs() < 1e-12,
        "{lo} to {hi}"
    );

    // A field is not a set of lines, and saying so beats drawing a plausible picture of nothing.
    let block = load("block.json");
    let field = block.frames[0]
        .panels
        .iter()
        .find(|p| matches!(p, Panel::Field { .. }))
        .unwrap();
    assert!(segments(field, &camera, &framing, 1.5, (0.0, 1.0)).is_empty());
}

/// **A kind the reader does not know is refused, not skipped.**
///
/// The failure being prevented is a viewer that silently drops a panel it cannot parse: the
/// window opens, something is missing, and nothing anywhere says what. A run written by a newer
/// library than the viewer is exactly when that happens.
#[test]
fn an_unknown_panel_kind_is_an_error() {
    let text = r#"{"title":"t","frames":[{"t":0.0,"panels":[
        {"name":"x","unit":"","kind":"isosurface","values":[1.0]}
    ],"readings":[]}]}"#;
    let err = Run::from_json(text).expect_err("an unknown kind must not parse");
    assert!(
        err.contains("isosurface"),
        "the message should name it: {err}"
    );

    // And a known kind missing a field it needs is refused too, rather than defaulted.
    let text = r#"{"title":"t","frames":[{"t":0.0,"panels":[
        {"name":"x","unit":"K","kind":"field","nx":2,"ny":1,"values":[1.0,2.0]}
    ],"readings":[]}]}"#;
    assert!(Run::from_json(text).is_err(), "nz is missing and matters");
}

/// **The camera's controls stay inside the range that keeps a picture on screen.**
#[test]
fn the_camera_cannot_be_driven_off_the_end() {
    let mut c = Camera::default();
    for _ in 0..200 {
        c.turn(0.1, 0.1);
    }
    assert!(c.elevation <= 1.5, "{c:?}");
    for _ in 0..200 {
        c.turn(0.0, -0.1);
    }
    assert!(c.elevation >= -1.5, "{c:?}");

    for _ in 0..200 {
        c.zoom(1.2);
    }
    assert!(c.distance <= 9.0, "{c:?}");
    for _ in 0..200 {
        c.zoom(0.8);
    }
    assert!(c.distance >= 1.2, "{c:?}");

    // And the projection still works at both ends, which is the reason for the clamps.
    let framing = Framing {
        centre: [0.0; 3],
        span: 1.0,
    };
    let p: Projected = c.project([0.5, 0.5, 0.5], &framing, 1.6);
    assert!(p.x.is_finite() && p.y.is_finite() && p.depth > 0.0);
}

/// **The shading is measured against the run, not against the frame in hand.**
///
/// The failure this prevents is the one `report.rs` documents from the other side: a scale that
/// re-fits every frame makes a quantity look constant while it changes by orders of magnitude.
/// A parcel of water that leaves a shower screen clean and reaches the spout loaded would render
/// mid-ramp the whole way, because at every instant it sits halfway between that instant's
/// lightest and darkest.
///
/// `Run::scale_of` existed for this and was tested, and the renderer did not call it. So this
/// checks the consumption rather than the accessor: **the same value in two frames whose local
/// ranges differ must come out the same shade.**
#[test]
fn the_shading_uses_the_runs_scale_and_not_the_frames() {
    let run = load("bench.json");
    let name = run.frames[0]
        .panels
        .iter()
        .find(|p| matches!(p, Panel::Paths { .. }))
        .expect("paths")
        .name()
        .to_string();
    let span = run.scale_of(&name).expect("the panel is in the run");
    let camera = Camera::default();
    let framing = Framing::of(run.framing_of(&name).expect("bounds"));

    // A value at the bottom of the run's range shades to 0 wherever it appears.
    let panel = &run.frames[0]
        .panels
        .iter()
        .find(|p| p.name() == name)
        .expect("panel");
    let wide = segments(panel, &camera, &framing, 1.5, span);
    assert!(!wide.is_empty(), "the bench traces to line segments");

    // The same panel against a span ten times wider must shade everything *lower*, which is only
    // true if the span is doing the work rather than the panel's own values.
    let stretched = (span.0, span.0 + (span.1 - span.0) * 10.0);
    let narrow = segments(panel, &camera, &framing, 1.5, stretched);
    let (a, b) = (
        wide.iter().map(|s| s.shade).fold(0.0f64, f64::max),
        narrow.iter().map(|s| s.shade).fold(0.0f64, f64::max),
    );
    println!("  brightest against the run {a:.4}, against a 10x span {b:.4}");
    assert!(
        b < a * 0.2 + 1e-9,
        "a ten-times-wider span must darken everything: {b:.4} against {a:.4}"
    );

    // A degenerate span is allowed and puts everything at the bottom rather than dividing by zero.
    let flat = segments(panel, &camera, &framing, 1.5, (1.0, 1.0));
    assert!(
        flat.iter().all(|s| s.shade.is_finite()),
        "a run with one value must not produce NaN shades"
    );
}

/// **A tall thin run and a cubic one both fill the frame.**
///
/// `Framing` normalises by the longest side, so a fixed camera distance is a distance chosen for
/// a cube: a portafilter at 67 mm across and 118 mm tall rendered at 15% of the frame height. In
/// a window that is a scroll away from fixed. In `--snapshot` it is not, and a snapshot of a
/// subject at a tenth of the frame is a weak check as much as a poor picture.
///
/// So the check is on the *fit*, at three quite different shapes, and it is two-sided: too small
/// wastes the frame and too large runs off it.
#[test]
fn the_camera_frames_whatever_shape_the_run_is() {
    for (label, bounds) in [
        ("a cube", [-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]),
        (
            "a portafilter",
            [-0.0335, -0.115, -0.0335, 0.0335, 0.003, 0.0335],
        ),
        ("a plate", [-1.0, -0.02, -1.0, 1.0, 0.02, 1.0]),
    ] {
        let framing = Framing::of(bounds);
        let mut camera = Camera::default();
        camera.fit(bounds, &framing, 16.0 / 9.0, 0.85);

        let mut worst: f64 = 0.0;
        for i in 0..8 {
            let c = [
                if i & 1 == 0 { bounds[0] } else { bounds[3] },
                if i & 2 == 0 { bounds[1] } else { bounds[4] },
                if i & 4 == 0 { bounds[2] } else { bounds[5] },
            ];
            let q = camera.project(c, &framing, 16.0 / 9.0);
            worst = worst.max(q.x.abs()).max(q.y.abs());
            assert!(
                q.depth > 0.05,
                "{label}: the camera must stay outside the subject, depth {:.4}",
                q.depth
            );
        }
        println!(
            "  {label:<14} focal {:.3}, furthest corner at {worst:.3}",
            camera.scale
        );
        assert!(
            (0.75..0.95).contains(&worst),
            "{label}: the box should nearly fill the frame, furthest corner at {worst:.3}"
        );
    }
}
