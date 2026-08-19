//! Every view this crate has, driven by frames that no physics produced.
//!
//! The frames here are written out by hand. That is the point rather than a shortcut: it is the
//! only way to test the claim this layer actually makes, which is that **the view is chosen by
//! the shape of the data and by nothing else**. A test that ran a real simulation to get its
//! frames could not distinguish "the report drew a heatmap because the data is a 2D grid" from
//! "the report drew a heatmap because that domain was a room".
//!
//! It also fixes the values, so the assertions can be about numbers rather than about the file
//! being non-empty — a renderer that silently drew nothing would pass every `len() > 0` check
//! ever written for it.

use pantometry_scene::{Frame, Panel, PanelData};
use pantometry_view::{html, readings_csv, svg, to_json};

/// A 3D field, a 2D field, a 1D field, some bodies, and two readings — one of each shape.
fn frames() -> Vec<Frame> {
    (0..4)
        .map(|k| {
            let t = k as f64 * 0.25;
            Frame {
                time_s: t,
                panels: vec![
                    Panel {
                        name: "sheet".into(),
                        unit: "Pa",
                        data: PanelData::Field {
                            nx: 3,
                            ny: 2,
                            // Row-major, and deliberately not symmetric: a renderer that
                            // transposed nx and ny would still produce six cells.
                            nz: 1,
                            values: vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0 + t],
                        },
                    },
                    Panel {
                        name: "lump".into(),
                        unit: "K",
                        data: PanelData::Field {
                            nx: 2,
                            ny: 2,
                            nz: 3,
                            // Twelve values, x fastest then y then z, and every slice different
                            // — so a view that drew slice 0 three times, or that read the array
                            // as one 2×6 plane, produces something this test can tell apart.
                            values: vec![
                                300.0,
                                301.0,
                                302.0,
                                303.0,
                                310.0,
                                311.0,
                                312.0,
                                313.0,
                                320.0,
                                321.0,
                                322.0,
                                323.0 + t,
                            ],
                        },
                    },
                    Panel {
                        name: "wire".into(),
                        unit: "K",
                        data: PanelData::Field {
                            nx: 4,
                            ny: 1,
                            nz: 1,
                            values: vec![300.0, 310.0, 305.0 + t, 300.0],
                        },
                    },
                    Panel {
                        name: "rays".into(),
                        unit: "nm",
                        // Two paths of three vertices each, going different ways, so a view that
                        // drew only the first or joined them into one is distinguishable.
                        data: PanelData::paths(
                            vec![
                                vec![[0.0, 0.0, 0.0], [1.0, 0.5, 0.0], [2.0, 0.0, 1.0 + t]],
                                vec![[0.0, 0.0, 0.0], [1.0, -0.5, 0.0], [2.0, -1.0, -1.0]],
                            ],
                            vec![486.1, 656.3],
                        ),
                    },
                    Panel {
                        name: "specks".into(),
                        unit: "m/s",
                        data: PanelData::Points {
                            positions: vec![[t, 0.0, 0.0], [0.0, 1.0, -1.0]],
                            // Climbs past the static body, which matters for
                            // `the_scale_does_not_move_between_frames`: while the moving value
                            // stayed under 2.0 the run-wide range *equalled* the first frame's
                            // and that test agreed with itself for the wrong reason.
                            values: vec![t * 12.0, 2.0],
                            bounds: [-1.0, -1.0, -1.0, 1.0, 1.0, 1.0],
                            boxed: true,
                        },
                    },
                ],
                readings: vec![
                    pantometry_core::Reading::new("box", "reserve", 10.0 - t, "J"),
                    pantometry_core::Reading::new("box", "temperature", 20.0 + t, "C"),
                ],
            }
        })
        .collect()
}

/// **Four shapes, four views, and the crate never sees a domain.**
///
/// Nothing in `frames()` came from a simulation. If the report can still pick a profile for the
/// one-row field, a heatmap for the two-row one, a scene for the bodies and a chart for the
/// readings, then it is dispatching on shape — which is the whole claim.
#[test]
fn the_report_picks_a_view_per_shape() {
    let page = html("hand-built", &frames());

    for kind in [
        "profile", "heatmap", "volume", "slices", "layout", "scene", "series",
    ] {
        assert!(
            page.contains(&format!("data-kind=\"{kind}\"")),
            "no {kind} view"
        );
    }
    // Five panels plus one card for the readings — and the 3D field gets **two**, a render and
    // a montage, because neither answers the other's question.
    assert_eq!(page.matches("class=\"card\"").count(), 7);

    // **The volume is not drawn as a plane.** `lump` is 2×2×3, and a view that ignored `nz`
    // would render its first four values and call it done — a perfectly plausible heatmap of a
    // solid, which is why the shape has to be checked rather than the picture. The dispatch
    // above puts it on `slices`, and the wire format carries the third count.
    assert!(
        page.contains("\"nz\":3"),
        "the third axis did not reach the page"
    );
    assert!(
        page.contains("3.230") || page.contains("3.23"),
        "the last slice's values are missing from the page"
    );

    // Self-contained, which is the promise that makes it useful to someone with no toolchain:
    // no network, no library, nothing to install. Checked because it is the kind of thing one
    // convenient `<script src=...>` undoes forever.
    for forbidden in ["http://", "https://", "src=", "@import"] {
        assert!(
            !page.contains(forbidden),
            "reaches outside itself: {forbidden}"
        );
    }
    assert!(page.starts_with("<!doctype html>"));
    assert!(page.trim_end().ends_with("</html>"));
}

/// **The units travel with the values.**
///
/// A legend that is separate from the data is a legend that can be wrong about it. `wire` is
/// kelvin and `sheet` is pascals in the same report, which is exactly the case where getting
/// this wrong produces a plausible picture.
#[test]
fn every_asset_carries_its_units() {
    let f = frames();
    let page = html("units", &f);
    for unit in ["Pa", "K", "m/s"] {
        assert!(page.contains(unit), "the report lost {unit}");
    }

    // The CSV puts them in the header, one column per `domain.label`.
    let csv = readings_csv(&f);
    let header = csv.lines().next().expect("a header");
    assert_eq!(header, "t_s,box.reserve [J],box.temperature [C]");
    assert_eq!(csv.lines().count(), f.len() + 1);

    // And the rows are the values, in order, so a plot of column 2 is a plot of the reserve.
    let last: Vec<f64> = csv
        .lines()
        .last()
        .unwrap()
        .split(',')
        .map(|s| s.parse().expect("a number"))
        .collect();
    assert_eq!(last.len(), 3);
    assert!((last[0] - 0.75).abs() < 1e-12, "t was {}", last[0]);
    assert!((last[1] - 9.25).abs() < 1e-12, "reserve was {}", last[1]);
}

/// **A run with nothing drawable still produces a report worth opening.**
///
/// The case that motivated the readings channel: a winding and a thermal network have no picture
/// at all, and for them the scalar *is* the result. An empty page here would be the layer failing
/// exactly where it is most needed.
#[test]
fn readings_alone_are_enough() {
    let bare: Vec<Frame> = frames()
        .into_iter()
        .map(|f| Frame {
            panels: Vec::new(),
            ..f
        })
        .collect();

    let page = html("nothing to draw", &bare);
    assert!(page.contains("data-kind=\"series\""), "no series view");
    assert_eq!(page.matches("class=\"card\"").count(), 1);
    assert!(page.contains("\"reserve\""), "the labels travel with it");

    // The filmstrip is empty rather than an SVG of a blank canvas, so a caller can tell the
    // difference between "nothing to draw" and "the renderer broke" and say something true.
    assert!(svg("nothing to draw", &bare, 4).is_empty());

    // The table is not empty, because this is precisely the run the table exists for.
    assert_eq!(readings_csv(&bare).lines().count(), bare.len() + 1);
}

/// **The filmstrip holds one colour scale across the whole run.**
///
/// A picture that renormalises per frame makes a decay look like a steady state. `specks` has a
/// body climbing from 0 to 9 while another sits at 2.0 throughout, so the first frame's own range
/// is `0..2` and the run's is `0..9`. Under per-frame normalisation the static body would be
/// drawn at the top of the scale in frame 0 and at a fifth of it in frame 3 — without having
/// moved.
#[test]
fn the_scale_does_not_move_between_frames() {
    let f = frames();
    let all = svg("a run", &f, 4);
    let first_only = svg("one frame", &f[..1], 1);

    // Same domain, same first frame, different runs: if the scale were per-frame or per-call the
    // first frame's markup would differ between these two.
    let first_cell = |s: &str| {
        let i = s.find("<circle").expect("bodies are drawn");
        s[i..i + 120].to_string()
    };
    assert_ne!(
        first_cell(&all),
        first_cell(&first_only),
        "a one-frame run and a four-frame run happen to agree, so this test proves nothing"
    );
}

/// **The JSON says what shape each panel is, and does not lose an axis.**
///
/// `positions` is flattened for the wire, which is the kind of place a z becomes a y. Two bodies
/// at six coordinates, and the second one's z is -1.
#[test]
fn the_json_keeps_all_three_axes() {
    let text = to_json("wire format", &frames());
    assert!(text.contains("\"kind\": \"points\""));
    assert!(text.contains("\"kind\": \"field\""));
    assert!(text.contains("\"boxed\": true"));

    let i = text.find("\"positions\":").expect("positions");
    let list = &text[i..text[i..].find(']').unwrap() + i + 1];
    assert_eq!(list.matches(',').count(), 5, "six coordinates: {list}");
    assert!(list.contains("-1.000000e0"), "the z was lost: {list}");
}

/// **Every view kind a card declares is dispatched, and every card has somewhere to write.**
///
/// The failure this prevents is adding a view and wiring it half way. A canvas with a `data-kind`
/// nothing draws is a blank rectangle; a card whose caption id nothing matches is a caption that
/// silently goes nowhere. Both look like a rendering that has not finished loading.
///
/// Checked structurally rather than by running the page, and the limit is worth naming: **nothing
/// in CI executes the viewer's JavaScript.** This catches a missing branch and a missing element,
/// not a branch that draws nothing. The volume view was verified once by running the script under
/// node against three real scenes and counting lit pixels — 4.8% for an acoustic mode, 2.3% for a
/// busbar, 0.1% for a single hot cell, which is why that last one now says so in its caption.
#[test]
fn every_declared_view_is_wired_to_something_that_draws_it() {
    let page = html("wiring", &frames());

    // Each `data-kind` on a canvas must have a branch in `drawAll`.
    let kinds: Vec<&str> = page
        .match_indices("data-kind=\"")
        .map(|(i, m)| {
            let rest = &page[i + m.len()..];
            &rest[..rest.find('"').expect("a closing quote")]
        })
        .collect();
    assert!(kinds.len() >= 7, "only {} canvases: {kinds:?}", kinds.len());

    // **Inside `drawAll`, not anywhere on the page.** The first version of this searched the whole
    // document, and passed with the dispatch deliberately broken — because the *rotation* filter
    // also names the kind, and `||v.kind==="volume"` matched there. A wiring check that any other
    // mention satisfies is not a wiring check.
    let body = {
        let start = page.find("function drawAll(){").expect("drawAll exists");
        let rest = &page[start..];
        &rest[..rest
            .find(
                "
}",
            )
            .expect("and it ends")]
    };
    for kind in &kinds {
        // `series` is the fall-through arm and names no draw call of its own.
        let dispatched = *kind == "series" || body.contains(&format!("v.kind===\"{kind}\") draw"));
        assert!(
            dispatched,
            "{kind} is declared on a canvas and drawAll never draws it"
        );
    }

    // Each canvas's slot must have a caption element, and the ids must be distinct — a volume has
    // two views of one panel, and keying the caption by panel made the second overwrite the first.
    let slots: Vec<&str> = page
        .match_indices("data-slot=\"")
        .map(|(i, m)| {
            let rest = &page[i + m.len()..];
            &rest[..rest.find('"').expect("a closing quote")]
        })
        .collect();
    assert_eq!(slots.len(), kinds.len(), "every canvas needs a slot");
    for slot in &slots {
        assert!(
            page.contains(&format!("id=\"cap-{slot}\"")),
            "{slot} has no caption element to write to"
        );
    }
    let mut unique = slots.clone();
    unique.sort_unstable();
    unique.dedup();
    assert_eq!(
        unique.len(),
        slots.len(),
        "two canvases share a slot: {slots:?}"
    );

    // And the volume and the montage really are the same panel seen twice.
    assert!(slots.contains(&"lump-volume") && slots.contains(&"lump-slices"));
}

/// **A path is a run, not a bag of points.**
///
/// `PanelData::paths` flattens runs into one array with an index, and the thing that can go wrong
/// is losing which vertex belongs to which path — at which point a traced ray becomes a scatter
/// and the one property that makes it a ray is gone.
#[test]
fn paths_keep_their_runs() {
    let f = frames();
    let panel = f[0].panels.iter().find(|p| p.name == "rays").unwrap();

    assert_eq!(panel.values(), &[486.1, 656.3]);
    assert_eq!(panel.path(0).unwrap().len(), 3);
    assert_eq!(panel.path(1).unwrap().len(), 3);
    assert_eq!(panel.path(0).unwrap()[1], [1.0, 0.5, 0.0]);
    assert_eq!(panel.path(1).unwrap()[1], [1.0, -0.5, 0.0]);
    assert!(panel.path(2).is_none(), "there is no third path");
    assert!(panel.grid().is_none(), "a path panel is not a grid");

    // A run of fewer than two points is dropped rather than kept as a degenerate line, and the
    // values follow the runs that survived rather than the ones that were offered.
    let thin = PanelData::paths(
        vec![vec![[0.0; 3]], vec![[0.0; 3], [1.0, 0.0, 0.0]]],
        vec![1.0, 2.0],
    );
    let panel = Panel {
        name: "thin".into(),
        unit: "",
        data: thin,
    };
    assert_eq!(
        panel.values(),
        &[2.0],
        "the surviving run keeps its own value"
    );
    assert_eq!(panel.path(0).unwrap().len(), 2);

    // And the box is measured from the geometry rather than left at a default.
    let PanelData::Paths { bounds, .. } = panel.data else {
        panic!("expected paths");
    };
    assert_eq!(bounds, [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]);
}
