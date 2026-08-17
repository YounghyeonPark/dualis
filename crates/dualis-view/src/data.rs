//! A run as data: a table for a plot, and a document for a viewer this crate does not contain.
//!
//! The two assets that are not pictures. Both exist because a picture is not always the answer —
//! a researcher with a plotting stack they already trust wants the numbers, and one building
//! their own viewer wants the frames.
//!
//! Written by hand rather than derived from `serde`, and that is deliberate for the JSON: the
//! moment anything reads it, it is a wire format, and a wire format should look like a decision
//! somebody made rather than like whatever the field names happened to be.
//!
//! It is also the **only** asset here that carries a three-dimensional field whole. The two
//! pictures show one slice of one, because a flat canvas cannot do otherwise; this carries every
//! sample, and is what a caller with a volume renderer should take.

use dualis_scene::{Frame, PanelData};

/// Every domain's scalars, one row per frame.
///
/// **The asset for the domains that have no picture at all** — a heater, a lamp, a winding, a
/// thermal network — and for several of those the scalar *is* the result. It is also the shape
/// a plot wants and the shape a spreadsheet wants.
///
/// Columns are `domain.label` so two networks with a node called `winding` do not collide, and
/// the unit is in the header rather than in a separate legend nobody reads.
pub fn readings_csv(frames: &[Frame]) -> String {
    let Some(first) = frames.first() else {
        return String::from("t_s\n");
    };
    let mut out = String::from("t_s");
    for r in &first.readings {
        out.push_str(&format!(",{}.{} [{}]", r.domain, r.label, r.unit));
    }
    out.push('\n');
    for frame in frames {
        out.push_str(&format!("{:.9}", frame.time_s));
        for r in &frame.readings {
            out.push_str(&format!(",{:.9}", r.value));
        }
        out.push('\n');
    }
    out
}

/// The frames as JSON, for a viewer this crate does not contain.
///
/// Fields as grids, bodies as positions in space, and the readings beside them. Written by hand
/// rather than derived, so the shape is chosen here and stays where a reader can see it — this
/// is a wire format the moment anything consumes it, and it should look deliberate.
pub fn to_json(title: &str, frames: &[Frame]) -> String {
    let mut out = format!("{{\n  \"title\": {},\n  \"frames\": [\n", quote(title));
    for (fi, frame) in frames.iter().enumerate() {
        out.push_str(&format!("    {{ \"t\": {:.6}, \"panels\": [", frame.time_s));
        for (pi, panel) in frame.panels.iter().enumerate() {
            out.push_str(&format!(
                "\n      {{ \"name\": {}, \"unit\": {}, ",
                quote(&panel.name),
                quote(panel.unit)
            ));
            match &panel.data {
                PanelData::Field { nx, ny, nz, values } => out.push_str(&format!(
                    "\"kind\": \"field\", \"nx\": {nx}, \"ny\": {ny}, \"nz\": {nz}, \
                     \"values\": {}",
                    numbers(values)
                )),
                PanelData::Paths {
                    vertices,
                    starts,
                    values,
                    bounds,
                } => {
                    let flat: Vec<f64> = vertices.iter().flatten().copied().collect();
                    let heads: Vec<f64> = starts.iter().map(|k| *k as f64).collect();
                    out.push_str(&format!(
                        "\"kind\": \"paths\", \"bounds\": {}, \"starts\": {}, \
                         \"vertices\": {}, \"values\": {}",
                        numbers(bounds),
                        numbers(&heads),
                        numbers(&flat),
                        numbers(values)
                    ));
                }
                PanelData::Points {
                    positions,
                    values,
                    bounds,
                    boxed,
                } => {
                    let flat: Vec<f64> = positions.iter().flatten().copied().collect();
                    out.push_str(&format!(
                        "\"kind\": \"points\", \"boxed\": {boxed}, \"bounds\": {}, \
                         \"positions\": {}, \"values\": {}",
                        numbers(bounds),
                        numbers(&flat),
                        numbers(values)
                    ));
                }
            }
            out.push_str(if pi + 1 == frame.panels.len() {
                " }"
            } else {
                " },"
            });
        }
        out.push_str("\n    ], \"readings\": [");
        for (ri, r) in frame.readings.iter().enumerate() {
            out.push_str(&format!(
                "\n      {{ \"domain\": {}, \"label\": {}, \"unit\": {}, \"value\": {:.6e} }}{}",
                quote(&r.domain),
                quote(&r.label),
                quote(r.unit),
                r.value,
                if ri + 1 == frame.readings.len() {
                    ""
                } else {
                    ","
                }
            ));
        }
        out.push_str(if fi + 1 == frames.len() {
            "\n    ] }\n"
        } else {
            "\n    ] },\n"
        });
    }
    out.push_str("  ]\n}\n");
    out
}

fn numbers(v: &[f64]) -> String {
    let mut s = String::from("[");
    for (i, x) in v.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        // `NaN` is not a JSON literal and would make the whole document unreadable to a strict
        // parser; `null` is the spelling for a sample that is not there, and numpy, pandas and
        // `JSON.parse` all take it.
        if x.is_finite() {
            s.push_str(&format!("{x:.6e}"));
        } else {
            s.push_str("null");
        }
    }
    s.push(']');
    s
}

fn quote(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if c.is_control() => out.push(' '),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}
