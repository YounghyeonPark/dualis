//! A filmstrip of heatmaps, as SVG, with no dependency.
//!
//! SVG because it is text: a `format!` and a file write, no encoder and no fonts. The same
//! reasoning as the library's examples — and the same code, very nearly, which is
//! `FRICTION.md` finding 4: the plotting the examples use lives under `examples/common/` and
//! cannot be reached from another crate.

use crate::Frame;

/// Draw every frame side by side, one row per domain.
pub fn filmstrip(title: &str, frames: &[Frame], columns: usize) -> String {
    let panels = frames.first().map(|f| f.panels.len()).unwrap_or(0);
    if panels == 0 {
        return String::new();
    }
    let (cell, pad, top) = (150.0f64, 8.0f64, 44.0f64);
    let columns = columns.max(1).min(frames.len().max(1));
    let rows = frames.len().div_ceil(columns);
    let w = columns as f64 * (cell + pad) + pad;
    let h = top + rows as f64 * panels as f64 * (cell + pad + 14.0) + pad;

    let mut s = format!(
        "<svg xmlns='http://www.w3.org/2000/svg' width='{w:.0}' height='{h:.0}' \
         viewBox='0 0 {w:.0} {h:.0}'>\n\
         <rect width='100%' height='100%' fill='#faf8f2'/>\n\
         <text x='{pad}' y='26' font-family='sans-serif' font-size='17' fill='#222'>{}</text>\n",
        escape(title)
    );

    // One scale for the whole strip, so frames are comparable with each other. Rescaling per
    // frame would make a decaying wave look like it never decays.
    let extent = frames
        .iter()
        .flat_map(|f| f.panels.iter())
        .flat_map(|p| p.values.iter())
        .fold(0.0f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE);

    for (k, frame) in frames.iter().enumerate() {
        let (col, row) = (k % columns, k / columns);
        let x0 = pad + col as f64 * (cell + pad);
        for (pi, panel) in frame.panels.iter().enumerate() {
            let y0 = top + (row * panels + pi) as f64 * (cell + pad + 14.0);
            s.push_str(&raster(panel, x0, y0, cell, extent));
            s.push_str(&format!(
                "<text x='{x0:.1}' y='{:.1}' font-family='sans-serif' font-size='9' \
                 fill='#555'>{} t={:.4}s</text>\n",
                y0 + cell + 11.0,
                escape(&panel.name),
                frame.time_s
            ));
        }
    }
    s.push_str("</svg>\n");
    s
}

/// The most cells to draw along either axis of a thumbnail.
///
/// A panel is 150 px, so beyond this the cells are sub-pixel and the detail is not being
/// seen. The *physics* is still sampled at full grid resolution in `World::capture` — this
/// coarsens the picture, not the measurement, and `tests/scene.rs` reads the full values.
const MAX_DRAWN: usize = 48;

/// How many colour steps the diverging ramp is quantised to, per side.
///
/// The quantisation is what makes the output small: cells of the same colour are collected
/// into one `<path>` instead of getting a `<rect>` each. One rect per cell put a 61x43 room
/// at 2.2 MB across twelve frames — about 70 bytes a cell — where a path subpath is
/// `M47 42h1v1h-1z`, fourteen. Forty-eight steps a side is finer than the eye separates on a
/// 3 px cell, so nothing visible is lost buying that.
const LEVELS: i32 = 48;

fn raster(p: &crate::Panel, x0: f64, y0: f64, size: f64, extent: f64) -> String {
    let (sx, sy) = (p.nx.max(1), p.ny.max(1));
    // Nearest neighbour, so an extremum survives if it lands on a chosen sample. Averaging
    // would be smoother and would hide exactly the overshoot worth seeing.
    let nx = sx.min(MAX_DRAWN);
    let ny = sy.min(MAX_DRAWN);

    // One bucket per quantised level, addressed by level + LEVELS so negatives fit.
    let mut buckets: Vec<String> = vec![String::new(); (2 * LEVELS + 1) as usize];
    for j in 0..ny {
        for i in 0..nx {
            let v = p.values[(j * sy / ny) * sx + (i * sx / nx)] / extent;
            let level = (v.clamp(-1.0, 1.0) * LEVELS as f64).round() as i32;
            // Rows are drawn top-down and the field's y runs up, so flip.
            buckets[(level + LEVELS) as usize].push_str(&format!("M{i} {}h1v1h-1z", ny - 1 - j));
        }
    }

    // Integer cell coordinates, scaled into place by the group, so no coordinate in the
    // path data needs a decimal point.
    let mut s = format!(
        "<g transform='translate({x0:.2} {y0:.2}) scale({:.4} {:.4})' \
         shape-rendering='crispEdges'>\n",
        size / nx as f64,
        size / ny as f64
    );
    for (k, d) in buckets.iter().enumerate() {
        if d.is_empty() {
            continue;
        }
        let t = (k as i32 - LEVELS) as f64 / LEVELS as f64;
        s.push_str(&format!("<path fill='{}' d='{d}'/>\n", diverging(t)));
    }
    s.push_str("</g>\n");
    s.push_str(&format!(
        "<rect x='{x0:.1}' y='{y0:.1}' width='{size:.1}' height='{size:.1}' fill='none' \
         stroke='#bbb' stroke-width='0.7'/>\n"
    ));
    s
}

/// Blue for negative, red for positive, near-white at zero — so the sign of a pressure is
/// visible, which a single-ended ramp hides.
fn diverging(t: f64) -> String {
    let t = t.clamp(-1.0, 1.0);
    let (r, g, b) = if t >= 0.0 {
        (
            255.0,
            255.0 - 175.0 * t.powf(0.75),
            250.0 - 220.0 * t.powf(0.75),
        )
    } else {
        let a = (-t).powf(0.75);
        (250.0 - 220.0 * a, 255.0 - 150.0 * a, 255.0)
    };
    format!("#{:02x}{:02x}{:02x}", r as u8, g as u8, b as u8)
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;")
}
