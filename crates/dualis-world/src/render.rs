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
        .flat_map(|p| p.values().iter())
        .fold(0.0f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE);

    for (k, frame) in frames.iter().enumerate() {
        let (col, row) = (k % columns, k / columns);
        let x0 = pad + col as f64 * (cell + pad);
        for (pi, panel) in frame.panels.iter().enumerate() {
            let y0 = top + (row * panels + pi) as f64 * (cell + pad + 14.0);
            s.push_str(&draw(panel, x0, y0, cell, extent));
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

/// One panel, in whichever shape it came in.
fn draw(p: &crate::Panel, x0: f64, y0: f64, size: f64, extent: f64) -> String {
    let body = match &p.data {
        crate::PanelData::Field { nx, ny, values } => {
            raster(*nx, *ny, values, x0, y0, size, extent)
        }
        crate::PanelData::Points {
            positions,
            values,
            bounds,
            boxed,
        } => scatter(positions, values, bounds, *boxed, x0, y0, size, extent),
    };
    body + &format!(
        "<rect x='{x0:.1}' y='{y0:.1}' width='{size:.1}' height='{size:.1}' fill='none' \
         stroke='#bbb' stroke-width='0.7'/>\n"
    )
}

/// Azimuth and elevation of the view, in radians.
///
/// Not isometric. A true isometric view puts the three axes at 120 degrees and makes a cube
/// ambiguous — the near-top and far-bottom corners land on the same point, which is a famous
/// optical illusion and a bad way to read a periodic box. These angles are close enough to
/// read as a cube and far enough off to be unambiguous.
const AZIMUTH: f64 = 0.61; // ~35 degrees
const ELEVATION: f64 = 0.49; // ~28 degrees

/// Project a point to the drawing plane, and report how far away it is.
///
/// Returns `(screen x, screen y, depth)`, with depth increasing away from the viewer.
fn project(p: [f64; 3]) -> (f64, f64, f64) {
    let (sa, ca) = (AZIMUTH.sin(), AZIMUTH.cos());
    let (se, ce) = (ELEVATION.sin(), ELEVATION.cos());
    let across = p[0] * ca - p[1] * sa;
    let depth = p[0] * sa + p[1] * ca;
    let up = p[2] * ce - depth * se;
    (across, up, depth)
}

/// The eight corners of a box.
fn corners(b: &[f64; 6]) -> [[f64; 3]; 8] {
    let (x0, y0, z0, x1, y1, z1) = (b[0], b[1], b[2], b[3], b[4], b[5]);
    [
        [x0, y0, z0],
        [x1, y0, z0],
        [x1, y1, z0],
        [x0, y1, z0],
        [x0, y0, z1],
        [x1, y0, z1],
        [x1, y1, z1],
        [x0, y1, z1],
    ]
}

/// Bodies as dots in an axonometric view, in a frame fixed for the whole run.
///
/// Three things make the depth readable, and the picture is flat without any of them.
/// Bodies are drawn **back to front**, so a near one covers a far one rather than whichever
/// happened to be last in the array. Radius grows toward the viewer. And the colour is mixed
/// toward the plate for distance, the way an aerial perspective works.
///
/// One `<circle>` each rather than the quantised paths a raster gets, because there are tens
/// or hundreds of these and not thousands, and because a dot's *position* is the information.
#[allow(clippy::too_many_arguments)]
fn scatter(
    positions: &[[f64; 3]],
    values: &[f64],
    bounds: &[f64; 6],
    boxed: bool,
    x0: f64,
    y0: f64,
    size: f64,
    extent: f64,
) -> String {
    // Frame the projected box, not the box itself: a rotated cube is wider than its side.
    let projected: Vec<(f64, f64, f64)> = corners(bounds).iter().map(|c| project(*c)).collect();
    let (mut ax0, mut ay0, mut ax1, mut ay1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    let (mut d0, mut d1) = (f64::MAX, f64::MIN);
    for (a, u, d) in &projected {
        ax0 = ax0.min(*a);
        ax1 = ax1.max(*a);
        ay0 = ay0.min(*u);
        ay1 = ay1.max(*u);
        d0 = d0.min(*d);
        d1 = d1.max(*d);
    }
    let span = (ax1 - ax0).max(ay1 - ay0).max(1e-30);
    // One scale for both axes, so a cube is not drawn as a cuboid.
    let pad = 0.06 * size;
    let inner = size - 2.0 * pad;
    let to_screen = |a: f64, u: f64| {
        (
            x0 + pad + (a - ax0) / span * inner,
            // SVG y runs down and the world's up runs up.
            y0 + pad + inner - (u - ay0) / span * inner,
        )
    };
    let depth_of = |d: f64| ((d - d0) / (d1 - d0).max(1e-30)).clamp(0.0, 1.0);

    let mut s = String::new();
    if boxed {
        // Twelve edges: four on the bottom, four on the top, four uprights.
        const EDGES: [(usize, usize); 12] = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        for (a, b) in EDGES {
            let (pa, pb) = (projected[a], projected[b]);
            let (xa, ya) = to_screen(pa.0, pa.1);
            let (xb, yb) = to_screen(pb.0, pb.1);
            // Far edges thinner and paler, which is enough to say which face is behind.
            let far = depth_of((pa.2 + pb.2) * 0.5);
            s.push_str(&format!(
                "<line x1='{xa:.2}' y1='{ya:.2}' x2='{xb:.2}' y2='{yb:.2}' \
                 stroke='#8d9099' stroke-width='{:.2}' opacity='{:.2}'/>\n",
                0.9 - 0.35 * far,
                0.55 - 0.3 * far
            ));
        }
    }

    // Painter's algorithm: farthest first.
    let mut order: Vec<usize> = (0..positions.len()).collect();
    order.sort_by(|a, b| {
        project(positions[*b])
            .2
            .partial_cmp(&project(positions[*a]).2)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let base = (inner / (positions.len().max(1) as f64).sqrt() / 2.6).clamp(0.9, 5.0);
    for k in order {
        let (a, u, d) = project(positions[k]);
        let (x, y) = to_screen(a, u);
        let far = depth_of(d);
        let r = base * (1.32 - 0.5 * far);
        let v = values.get(k).copied().unwrap_or(0.0) / extent;
        s.push_str(&format!(
            "<circle cx='{x:.2}' cy='{y:.2}' r='{r:.2}' fill='{}'/>\n",
            faded(v, far)
        ));
    }
    s
}

/// The diverging colour, mixed toward the plate by distance.
///
/// Aerial perspective: what is far away loses contrast against the background. Without it a
/// hundred atoms are an even wash and the box has no depth at all.
fn faded(t: f64, far: f64) -> String {
    let (r, g, b) = ramp(t);
    let m = 0.62 * far.clamp(0.0, 1.0);
    let mix = |c: f64, p: f64| c + (p - c) * m;
    format!(
        "#{:02x}{:02x}{:02x}",
        mix(r, 250.0) as u8,
        mix(g, 248.0) as u8,
        mix(b, 242.0) as u8
    )
}

#[allow(clippy::too_many_arguments)]
fn raster(
    sx: usize,
    sy: usize,
    values: &[f64],
    x0: f64,
    y0: f64,
    size: f64,
    extent: f64,
) -> String {
    let (sx, sy) = (sx.max(1), sy.max(1));
    // Nearest neighbour, so an extremum survives if it lands on a chosen sample. Averaging
    // would be smoother and would hide exactly the overshoot worth seeing.
    let nx = sx.min(MAX_DRAWN);
    let ny = sy.min(MAX_DRAWN);

    // One bucket per quantised level, addressed by level + LEVELS so negatives fit.
    let mut buckets: Vec<String> = vec![String::new(); (2 * LEVELS + 1) as usize];
    for j in 0..ny {
        for i in 0..nx {
            let v = values[(j * sy / ny) * sx + (i * sx / nx)] / extent;
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
    s
}

/// Blue for negative, red for positive, near-white at zero — so the sign of a pressure is
/// visible, which a single-ended ramp hides.
fn diverging(t: f64) -> String {
    let (r, g, b) = ramp(t);
    format!("#{:02x}{:02x}{:02x}", r as u8, g as u8, b as u8)
}

/// The ramp itself, as raw channels, so distance can be mixed in before it is written out.
fn ramp(t: f64) -> (f64, f64, f64) {
    let t = t.clamp(-1.0, 1.0);
    if t >= 0.0 {
        (
            255.0,
            255.0 - 175.0 * t.powf(0.75),
            250.0 - 220.0 * t.powf(0.75),
        )
    } else {
        let a = (-t).powf(0.75);
        (250.0 - 220.0 * a, 255.0 - 150.0 * a, 255.0)
    }
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;")
}
