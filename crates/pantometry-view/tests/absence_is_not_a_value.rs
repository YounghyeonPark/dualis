//! What each writer does with a sample that is not a number.
//!
//! A field gained cells with nothing in them, and every writer here had an answer for them that
//! looked like an answer for something. Measured on a real scene before the fix: the glTF carried
//! 252 points for a grid with 120 solid cells, so two thirds of what went to Blender was material
//! that is not there; the SVG drew the empty ones in the **middle** colour of the ramp, because
//! `NaN as i32` saturates to zero and zero is the middle bucket, so an empty cell was the same
//! colour as a cell at the average; and the JSON writers would have emitted the literal `NaN`,
//! which is not a JSON token at all.
//!
//! One rule, four writers: **absence is not a value.** These are on the writers. The sampler that
//! decides what is absent has its own, in `pantometry-thermal`.

use pantometry_scene::{Frame, Panel, PanelData};

/// A 3x1x3 field with the middle *column* empty, so the hole survives any projection.
///
/// The values are `300 + 10k`, distinct per slice, so an assertion can say where in the array a
/// hole is and not merely how many there are.
fn holed() -> Frame {
    let mut values = Vec::new();
    for k in 0..3 {
        for i in 0..3 {
            values.push(if i == 1 {
                f64::NAN
            } else {
                300.0 + 10.0 * k as f64
            });
        }
    }
    Frame {
        time_s: 0.0,
        panels: vec![Panel {
            name: "block".into(),
            unit: "K",
            data: PanelData::Field {
                nx: 3,
                ny: 1,
                nz: 3,
                values,
            },
        }],
        readings: Vec::new(),
    }
}

/// **The glTF carries a point for what is there and none for what is not.**
///
/// This is the one output that leaves the workspace — Blender, three.js, a USD tool — so a
/// clearance exported as vertices becomes a wrong model in somebody else's software, where nothing
/// about it says where it came from. Nine cells, three empty, six points.
#[test]
fn the_gltf_leaves_the_holes_out() {
    let doc = pantometry_view::gltf("holed", &holed()).document;
    assert!(
        doc.contains("\"count\":6"),
        "six of the nine cells have something in them, and only those are geometry: {doc}"
    );
    assert!(
        !doc.contains("\"count\":9"),
        "and the empty three are not quietly along for the ride: {doc}"
    );
}

/// **The JSON writers say `null`, in the right places, and never say `NaN`.**
///
/// `NaN` is not a token this format has. A document containing one is refused outright by a strict
/// parser — which turns a hole in a field into a file nobody can open — and accepted by a lenient
/// one, which is worse, because it comes back as a float and gets plotted.
///
/// Positional rather than a count: the array has to be `300, nothing, 300` and not three nulls
/// swept to the end, because a writer that dropped the empties would shift every later cell.
#[test]
fn the_json_writers_spell_absence_null_and_leave_it_where_it_was() {
    for (what, text) in [
        ("to_json", pantometry_view::to_json("holed", &[holed()])),
        ("html", pantometry_view::html("holed", &[holed()])),
    ] {
        // Spelled as it would appear *in an array*, because the report's own script names the
        // JavaScript constant legitimately — a bare search for `NaN` matches that and says the
        // writer is broken when it is not.
        assert!(
            !text.contains(",NaN") && !text.contains("[NaN") && !text.contains("NaN]"),
            "{what}: `NaN` is not a JSON token and must not reach a values array"
        );
        assert!(
            text.contains("3.000000e2,null,3.000000e2")
                && text.contains("3.100000e2,null,3.100000e2")
                && text.contains("3.200000e2,null,3.200000e2"),
            "{what}: the hole stays in the slot it was in"
        );
    }
}

/// **The SVG draws no square where there is nothing**, and the background shows through.
///
/// The failure was silent in the exact way that matters: nothing about a mid-ramp square says it is
/// not a measurement. Counted by the squares the raster emits — the filmstrip draws the middle
/// z-slice of a volume, which here is three cells with one of them empty.
#[test]
fn the_svg_draws_no_square_for_an_empty_cell() {
    let svg = pantometry_view::svg("holed", &[holed()], 1);
    assert_eq!(
        svg.matches("h1v1h-1z").count(),
        2,
        "one square per cell that has a value, and none for the one that does not"
    );
}
