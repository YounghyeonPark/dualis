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

use dualis_scene::{Frame, Panel, PanelData};
use dualis_view::{html, readings_csv, svg, to_json};

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
                    dualis_core::Reading::new("box", "reserve", 10.0 - t, "J"),
                    dualis_core::Reading::new("box", "temperature", 20.0 + t, "C"),
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

    for kind in ["profile", "heatmap", "slices", "scene", "series"] {
        assert!(
            page.contains(&format!("data-kind=\"{kind}\"")),
            "no {kind} view"
        );
    }
    // Four panels plus one card for the readings.
    assert_eq!(page.matches("class=\"card\"").count(), 5);

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
