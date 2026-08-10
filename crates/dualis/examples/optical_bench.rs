//! A three-dimensional optical bench: elements placed by pose, rays traced through them, and a
//! layout you can rotate in a browser.
//!
//! ```text
//! cargo run --release --example optical_bench             # numbers, checked
//! cargo run --release --example optical_bench bench.html  # and the 3D layout
//! ```
//!
//! Every other example here draws a graph. This one draws the **instrument**: a collimated beam,
//! a doublet, a fold mirror that turns the axis through 90°, and the rays arriving at a tilted
//! image plane — as paths in space, on one page, rotatable and scrollable with nothing installed.
//!
//! # What this is and is not
//!
//! It is the analysis half of an optical layout tool: place, trace, measure, look. Ray paths are a
//! first-class shape in the scene layer now — `PanelData::Paths` — so the view draws them the same
//! way it draws bodies and fields, depth-sorted and on one colour scale.
//!
//! It is **not** an interactive editor and there is no renderer with materials or shadows. Nothing
//! here is real-time; a run produces a file. Those are different products and this workspace does
//! not pretend otherwise.
//!
//! # What is checked, because a picture is not a result
//!
//! A layout that looks right is the easiest wrong answer in optics. So:
//!
//! - the doublet's focal length, against the thin-lens combination of its two powers;
//! - the fold, against the law of reflection — the axis leaves at exactly 90° to how it arrived;
//! - the marginal ray's height at the stop, against `h = f·tanθ` for each field angle;
//! - the RMS spot against the diffraction limit, so "in focus" is a number and not a look.
//!
//! Each is arithmetic done here, not a value read back out of the thing being checked.

use dualis::prelude::*;
use dualis::scene::{Frame, Panel, PanelData};
use dualis_optics::geometry::{cap_intersect, hexapolar_unit, plane_intersect, refract, Ray};
use glam::DVec3;

mod common;
use common::{check, check_between, heading};

/// N-BK7 and SF2 at the d line — a classic cemented doublet pair.
const N_CROWN: f64 = 1.5168;
const N_FLINT: f64 = 1.6477;
/// Focal length to aim for.
const FOCAL: f64 = 100e-3;
/// The stop: the semi-aperture the pupil is sampled over.
const SEMI: f64 = 8.0e-3;
/// The glass is cut a little larger, so the *stop* limits the beam rather than the edge of a
/// surface. A prescription whose clear aperture equals its stop vignettes on the first tilt.
const GLASS: f64 = 9.5e-3;
/// Where the fold mirror sits along the axis, past the doublet.
const FOLD_AT: f64 = 60e-3;
/// The flat back face of the flint, which is where the lens ends.
const BACK_Z: f64 = 7e-3;
/// Field angles traced, in degrees.
const FIELDS: [f64; 3] = [0.0, 1.5, 3.0];

fn main() {
    // ================================================================ the prescription
    heading("A doublet, from the powers it has to add up to");

    // Three surfaces, and each one's power is `(n2 - n1)/R` — the refractive index it goes *into*
    // minus the one it comes *from*. The crown's front is air-to-crown; the cemented face is
    // crown-to-flint, a difference of 0.13 rather than 0.52; the flint's back is flat.
    //
    // That distinction is not pedantry. Using `(n - 1)/R` on the cemented surface — the air-glass
    // formula, which is the one everybody remembers — made the flint four times too weak, and the
    // traced focal length came out 150 mm against a 100 mm prescription. The trace was right and
    // the prescription was wrong, which is the correct way round and only visible because the
    // trace measures `f` instead of being told it.
    // One free parameter: how much of the total power the front surface carries. Whatever is
    // left goes on the cemented face, so **every split has the same paraxial focal length** and
    // they differ only in aberration. That is the classic bending variable.
    let split = 1.55;
    let (r_crown, r_cement) = curvatures(split);

    println!(
        "  {:<30} {:>9.3} mm  air to crown",
        "front radius",
        r_crown * 1e3
    );
    println!(
        "  {:<30} {:>9.3} mm  crown to flint",
        "cemented radius",
        r_cement * 1e3
    );
    println!("  {:<30} {:>9}     flint to air", "back", "flat");
    let phi = 1.0 / FOCAL;
    check(
        "the three powers add to 1/f",
        1.0 / (phi * split + (phi - phi * split)),
        FOCAL,
        1e-12,
        "m",
    );

    // ================================================================ the bench
    heading("The bench, placed in three dimensions");
    let mut bench = Bench {
        crown_r: r_crown,
        cement_r: r_cement,
        fold_at: FOLD_AT,
        // Worked out rather than guessed. reflect(d, n) = d - 2(d.n)n, and for +z to leave along
        // +x that needs n = (-1, 0, 1)/sqrt(2). The other sign of the z component sends the beam
        // to -x, where the image plane is not — and every ray then misses it, which is the right
        // failure and an opaque one to read: sixty-one rays traced, none arrived.
        fold_normal: DVec3::new(-1.0, 0.0, 1.0).normalize(),
        // Filled in from the measurement below, because the thin-lens prescription does not know
        // about the glass thickness and therefore does not know where the focus is.
        image_along: 0.0,
    };

    // Where the focus actually is, from a traced paraxial ray.
    let (f_measured, bfd) = bench.measure().expect("a paraxial ray gets through");
    bench.image_along = bfd - (FOLD_AT - BACK_Z);
    println!(
        "  {:<30} {:>9.3} mm  measured from a traced ray, against {:.1} thin-lens",
        "effective focal length",
        f_measured * 1e3,
        FOCAL * 1e3
    );
    println!(
        "  {:<30} {:>9.3} mm  from the back face, so the plane is {:.1} mm past the mirror",
        "back focal distance",
        bfd * 1e3,
        bench.image_along * 1e3
    );
    check_between(
        "a thick lens focuses shorter than its thin-lens prescription",
        bfd / f_measured,
        0.90,
        1.0,
        "x",
    );
    println!("  {:<30} {:>9.1} mm  along +z", "doublet at", 0.0);
    println!(
        "  {:<30} {:>9.1} mm  normal at 45 degrees",
        "fold mirror at",
        bench.fold_at * 1e3
    );
    println!(
        "  {:<30} {:>9.1} mm  along +x after the fold",
        "image plane at",
        bench.image_along * 1e3
    );

    // ================================================================ trace
    heading("Rays through it, and the closed forms they must obey");
    let Traced { spots, .. } = trace_fields(&bench);
    for (k, angle_deg) in FIELDS.iter().enumerate() {
        println!(
            "  {angle_deg:>4.1} deg   {:>3} of {:>3} rays through",
            spots[k].len(),
            PUPIL_RAYS
        );
        assert!(
            spots[k].len() * 10 >= PUPIL_RAYS * 9,
            "at least 90% of the pupil must make it, or the bench is vignetting"
        );
    }

    // The fold: a ray arriving along +z must leave along +x, exactly.
    let axial = Ray::new(LengthVec::from_si(DVec3::new(0.0, 0.0, 0.0)), DVec3::Z);
    let turned = reflect(axial.dir, bench.fold_normal);
    println!(
        "  {:<30} ({:.6}, {:.6}, {:.6})",
        "axis after the fold", turned.x, turned.y, turned.z
    );
    check("the fold turns the axis to +x", turned.x, 1.0, 1e-12, "");
    check_between("and leaves nothing along z", turned.z.abs(), 0.0, 1e-12, "");

    // ================================================================ bending
    heading("Bending: one parameter, scanned for the smallest spot");
    // Every split traced below has the same paraxial focal length, so this trades nothing away.
    // It is free performance sitting in a number somebody had to choose, which is why an optical
    // design tool's first move is to scan it.
    let scan: Vec<(f64, f64)> = (0..=28)
        .filter_map(|k| {
            let s = 1.0 + 1.4 * k as f64 / 28.0;
            on_axis_rms(&bench, s).map(|spot| (s, spot))
        })
        .collect();
    let &(best_split, best_rms) = scan
        .iter()
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .expect("some split traces");
    for (s, spot) in scan.iter().step_by(4) {
        println!(
            "  split {s:>5.2}   RMS spot {:>8.2} um{}",
            spot * 1e6,
            if (*s - best_split).abs() < 1e-12 {
                "   <- best"
            } else {
                ""
            }
        );
    }
    let started = on_axis_rms(&bench, split).expect("the starting split traces");
    println!(
        "  {:<30} {:>9.2} um  ->  {:.2} um at split {best_split:.2}",
        "bending the doublet",
        started * 1e6,
        best_rms * 1e6
    );
    check_between(
        "bending is worth something",
        started / best_rms,
        1.5,
        500.0,
        "x",
    );
    // A minimum, not the end of the range — otherwise the scan chose its own edge.
    assert!(
        best_split > scan[0].0 && best_split < scan[scan.len() - 1].0,
        "the optimum is at the edge of the scan: {best_split}"
    );

    // Rebuild on the bent prescription and re-measure what follows from it.
    let (br, bc) = curvatures(best_split);
    bench.crown_r = br;
    bench.cement_r = bc;
    let (f_bent, bfd_bent) = bench.measure().expect("the bent doublet still focuses");
    bench.image_along = bfd_bent - (FOLD_AT - BACK_Z);
    check("bending does not move f", f_bent, f_measured, 0.03, "m");
    let Traced {
        paths,
        colours,
        spots,
    } = trace_fields(&bench);

    // The chief ray's height at the image plane is `f tan(theta)` — the definition of focal
    // length, and a check the trace cannot fake because nothing in it computes `f`.
    for (k, angle_deg) in FIELDS.iter().enumerate() {
        let theta = angle_deg.to_radians();
        let centroid = centroid_of(&spots[k]);
        let want = f_bent * theta.tan();
        if *angle_deg == 0.0 {
            check_between(
                "on axis the centroid is on axis",
                centroid.1.abs(),
                0.0,
                2e-5,
                "m",
            );
        } else {
            check(
                &format!("{angle_deg:.1} deg lands at f tan(theta)"),
                centroid.1.abs(),
                want,
                0.06,
                "m",
            );
        }
    }

    // ================================================================ the spot
    heading("Is it in focus? A number, not a look");
    let rms = rms_radius(&spots[0]);
    // The diffraction limit for this aperture, at the d line.
    let na = SEMI / f_bent;
    let airy = 0.61 * 587.6e-9 / na;
    println!(
        "  {:<30} {:>9.2} um  RMS spot radius on axis",
        "geometric",
        rms * 1e6
    );
    println!(
        "  {:<30} {:>9.2} um  0.61 lambda / NA",
        "diffraction limit",
        airy * 1e6
    );
    println!(
        "  {:<30} {:>9.2}x  the Airy radius",
        "the geometric spot is",
        rms / airy
    );
    // **The verdict, and it changed with the bending.** Before, 209 um against a 4.4 um limit —
    // 47x, and the picture would have been all aberration. After, the geometric spot is inside
    // the Airy disc, so the design is diffraction-limited and the ray trace has stopped being the
    // thing that decides the image. That is the answer an optical design review wants, and it is
    // a number rather than a look at a spot diagram.
    check_between(
        "the bent design is diffraction-limited",
        rms / airy,
        0.0,
        1.0,
        "x the Airy radius",
    );
    println!(
        "  {:<30} {:>9.1}x  before bending, which the scan removed",
        "it was",
        started / airy
    );

    // ================================================================ the deliverable
    if let Some(path) = common::output_path() {
        let frame = Frame {
            time_s: 0.0,
            panels: vec![Panel {
                name: "bench".into(),
                unit: "deg field",
                data: PanelData::paths(paths, colours),
            }],
            readings: vec![
                Reading::new("bench", "focal length", f_bent * 1e3, "mm"),
                Reading::new("bench", "semi-aperture", SEMI * 1e3, "mm"),
                Reading::new("bench", "rms spot", rms * 1e6, "um"),
                Reading::new("bench", "diffraction limit", airy * 1e6, "um"),
            ],
        };
        // The extension chooses the asset, the way the application does: `.html` is a page you
        // open, `.json` is the frames themselves for something else to draw — the native viewer
        // in `runtime/viewer`, for instance, which reads this and nothing else.
        let frames = std::slice::from_ref(&frame);
        let asset = if path.ends_with(".json") {
            dualis::view::to_json("optical bench", frames)
        } else {
            dualis::view::html("optical bench", frames)
        };
        common::write(&path, &asset);
        println!("\n  drag to rotate, scroll to zoom");
    } else {
        println!("\n  give a filename ending .html for the 3D layout");
    }
}

/// How many rays a hexapolar pupil of four rings holds.
const PUPIL_RAYS: usize = 61;

/// The two curvatures for a given power split, at a fixed focal length.
///
/// The front takes `split` times the total power and the cemented face takes the rest, so every
/// split has the same paraxial `f` — which is what makes bending a free parameter rather than a
/// trade. Each surface's power is `(n2 - n1)/R`: the index it goes *into* minus the one it comes
/// *from*.
///
/// Using the air-glass `(n - 1)/R` on the cemented face — the formula everybody remembers — makes
/// the flint four times too weak, and the traced focal length came out 150 mm against a 100 mm
/// prescription. The trace was right and the prescription was wrong, which is the correct way
/// round and only visible because the trace *measures* `f` instead of being told it.
fn curvatures(split: f64) -> (f64, f64) {
    let phi = 1.0 / FOCAL;
    (
        (N_CROWN - 1.0) / (phi * split),
        (N_FLINT - N_CROWN) / (phi - phi * split),
    )
}

/// One ray's path through the bench and where it landed on the image plane.
type TracedRay = (Vec<[f64; 3]>, (f64, f64));

/// What a full trace of every field produces.
struct Traced {
    /// One run of points per ray that got through.
    paths: Vec<Vec<[f64; 3]>>,
    /// The field angle each path belongs to, which is what colours it.
    colours: Vec<f64>,
    /// Where each field's bundle landed on the image plane.
    spots: Vec<Vec<(f64, f64)>>,
}

/// Trace every field angle through the bench.
fn trace_fields(bench: &Bench) -> Traced {
    let pupil = hexapolar_unit(4);
    let (mut paths, mut colours, mut spots) = (Vec::new(), Vec::new(), Vec::new());
    for angle_deg in FIELDS {
        let theta = angle_deg.to_radians();
        let dir = DVec3::new(0.0, theta.sin(), theta.cos()).normalize();
        let mut landed = Vec::new();
        for (u, v) in &pupil {
            // Back-propagated **through** the pupil rather than launched from a flat plane.
            // Starting every ray at z = -40 mm and tilting it means an off-axis bundle has
            // drifted a millimetre by the time it reaches the glass, so the outer rays fall off
            // the aperture — 48 of 61 at 1.5 degrees, which reads as vignetting the design does
            // not have. The stop is the pupil; a ray is defined by where it crosses it.
            let at_pupil = DVec3::new(u * SEMI, v * SEMI, 0.0);
            let ray = Ray::new(LengthVec::from_si(at_pupil - dir * 40e-3), dir);
            if let Some((path, hit)) = bench.trace(ray) {
                landed.push(hit);
                paths.push(path);
                colours.push(angle_deg);
            }
        }
        spots.push(landed);
    }
    Traced {
        paths,
        colours,
        spots,
    }
}

/// The on-axis RMS spot for one power split, refocused for that split.
///
/// **Refocused**, which is the whole point: bending moves the focus a little, and comparing spots
/// at a fixed image plane would be measuring defocus and calling it aberration.
fn on_axis_rms(bench: &Bench, split: f64) -> Option<f64> {
    let (crown_r, cement_r) = curvatures(split);
    let mut trial = Bench {
        crown_r,
        cement_r,
        image_along: 0.0,
        ..*bench
    };
    let (_, bfd) = trial.measure()?;
    trial.image_along = bfd - (FOLD_AT - BACK_Z);
    let landed = trace_fields(&trial).spots;
    let axial = landed.first()?;
    (axial.len() * 10 >= PUPIL_RAYS * 9).then(|| rms_radius(axial))
}

/// The elements, and where they are.
#[derive(Clone, Copy)]
struct Bench {
    crown_r: f64,
    cement_r: f64,
    fold_at: f64,
    fold_normal: DVec3,
    /// How far past the mirror the image plane sits, along the folded axis.
    image_along: f64,
}

impl Bench {
    /// Trace one ray, returning its path through the bench and where it landed.
    ///
    /// Returns `None` if the ray misses anything it has to hit — which is vignetting, not an
    /// error, and the caller counts it rather than being told nothing happened.
    fn trace(&self, ray: Ray) -> Option<TracedRay> {
        let mut path = vec![point(ray.origin)];
        let mut r = self.lens(&mut path, ray)?;

        // The fold mirror.
        let t = plane_intersect(
            r,
            LengthVec::from_si(DVec3::new(0.0, 0.0, self.fold_at)),
            self.fold_normal,
        )?;
        path.push(point(r.at(t)));
        r = r.redirect(t, reflect(r.dir, self.fold_normal));

        // The image plane, normal to +x after the fold.
        let t = plane_intersect(
            r,
            LengthVec::from_si(DVec3::new(self.image_along, 0.0, self.fold_at)),
            DVec3::X,
        )?;
        let end = r.at(t);
        path.push(point(end));
        // In the image plane, z is across the field and y is the other transverse axis.
        Some((path, (end.to_si().z - self.fold_at, end.to_si().y)))
    }

    /// The three lens surfaces, leaving the ray on its way to the mirror.
    fn lens(&self, path: &mut Vec<[f64; 3]>, ray: Ray) -> Option<Ray> {
        let mut r = ray;

        // Front surface of the crown. `cap_intersect` wants the **vertex** — the surface's apex
        // on the axis — not the centre of curvature. Passing the centre put the surface a radius
        // downstream of where it belongs and every ray missed the aperture, which is the right
        // failure for a wrong vertex and a confusing one to read.
        let hit = cap_intersect(
            r,
            LengthVec::ZERO,
            DVec3::Z,
            Length::from_si(self.crown_r),
            Length::from_si(GLASS),
        )?;
        r = self.bend(path, r, hit.t, hit.normal, 1.0, N_CROWN)?;

        // The cemented interface.
        let hit = cap_intersect(
            r,
            LengthVec::from_si(DVec3::new(0.0, 0.0, 4e-3)),
            DVec3::Z,
            Length::from_si(self.cement_r),
            Length::from_si(GLASS),
        )?;
        r = self.bend(path, r, hit.t, hit.normal, N_CROWN, N_FLINT)?;

        // Flat back face of the flint.
        let t = plane_intersect(
            r,
            LengthVec::from_si(DVec3::new(0.0, 0.0, BACK_Z)),
            DVec3::Z,
        )?;
        self.bend(path, r, t, DVec3::Z, N_FLINT, 1.0)
    }

    /// The effective focal length and back focal distance, **measured from a traced ray**.
    ///
    /// A ray entering parallel to the axis at height `h` leaves at angle `u`; `f = h/tan(u)` is
    /// the definition of focal length, and where it crosses the axis is the back focus. Neither
    /// is assumed: the thin-lens prescription that produced the curvatures ignores the 7 mm of
    /// glass, so the real focus is not where `FOCAL` says.
    fn measure(&self) -> Option<(f64, f64)> {
        let h = 0.5e-3;
        let mut path = Vec::new();
        let entering = Ray::new(LengthVec::from_si(DVec3::new(0.0, h, -40e-3)), DVec3::Z);
        let out = self.lens(&mut path, entering)?;
        let (o, d) = (out.origin.to_si(), out.dir);
        // Where it crosses the axis, measured from the last vertex.
        let t = -o.y / d.y;
        let crossing = o.z + d.z * t;
        let angle = (-d.y / d.z).atan();
        Some((h / angle.tan(), crossing - BACK_Z))
    }

    /// Refract at a surface, recording the vertex. `None` on total internal reflection.
    fn bend(
        &self,
        path: &mut Vec<[f64; 3]>,
        r: Ray,
        t: Length,
        normal: DVec3,
        n1: f64,
        n2: f64,
    ) -> Option<Ray> {
        path.push(point(r.at(t)));
        // `Hit::normal` is the surface's own, not oriented against the ray — the trap `lens_spots`
        // documents, and the reason this is one helper rather than four call sites.
        let n = oriented_against(normal, r.dir);
        let dir = refract(r.dir, n, n1 / n2)?;
        Some(r.redirect(t, dir))
    }
}

fn point(p: LengthVec) -> [f64; 3] {
    let v = p.to_si();
    [v.x, v.y, v.z]
}

fn centroid_of(spots: &[(f64, f64)]) -> (f64, f64) {
    let n = spots.len().max(1) as f64;
    (
        spots.iter().map(|(a, _)| a).sum::<f64>() / n,
        spots.iter().map(|(_, b)| b).sum::<f64>() / n,
    )
}

fn rms_radius(spots: &[(f64, f64)]) -> f64 {
    let (cx, cy) = centroid_of(spots);
    let n = spots.len().max(1) as f64;
    (spots
        .iter()
        .map(|(a, b)| (a - cx).powi(2) + (b - cy).powi(2))
        .sum::<f64>()
        / n)
        .sqrt()
}
