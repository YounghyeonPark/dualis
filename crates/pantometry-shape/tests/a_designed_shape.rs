//! A mesh read, measured and rasterised, against what the geometry says it should be.
//!
//! Four closed forms, and they check different halves:
//!
//! - **A box's mesh volume is `w·h·d` exactly.** The divergence-theorem sum has no sampling and no
//!   tolerance in it, so this is machine precision and it checks the *reader*: a vertex read into the
//!   wrong slot, a millimetre conversion applied twice, a winding taken backwards all fail it.
//! - **A cell-aligned box rasterises exactly.** Zero volume error at every resolution, which checks the
//!   rasteriser against the one case where the answer is not approximate. An `O(dx)` scheme and an exact
//!   one look identical at a single resolution otherwise.
//! - **A box's boundary layer is the box less its interior**, `1 − (nx−2)(ny−2)(nz−2)/(nx·ny·nz)`, with
//!   no tolerance at all. The sharpest check available on the layer's thickness.
//! - **A sphere's boundary layer is `A·dx/V`**: first order, with a coefficient bracketed in
//!   `[0.25, 1.5]` by counting column ends. Against the *mesh's* own volume rather than `4πr³/3`, because
//!   a tessellated sphere is not a sphere and the analytic value would fold the tessellation's error into
//!   the grid's.
//! - **A cube's cell count is exact at every resolution.** A regression test rather than a new idea, and
//!   the shape a box cannot stand in for: a cube's face diagonal runs through cell centres and a box's
//!   generally does not.
//!
//! The sphere's *volume* error is deliberately not among them. It partly cancels, and what is left is a
//! lattice-point count that a refinement can make worse — asserted rather than assumed. It does **not**
//! cancel to zero: the rasteriser over-fills a convex body, measurably, and that bias is recorded rather
//! than asserted away.
//!
//! Two things that are checked here because nothing else can check them. **Where the rows were sampled**,
//! by asking that no cell centre outside the sphere is filled — every aggregate above survives a sampler
//! moved a quarter cell off centre, and that one does not. And **that the degenerate-ray retry fires at
//! all**, by counting the rows that needed it, because `ambiguous_rows` — the case where every retry
//! fails — has never been produced by any mesh here and cannot serve as the witness.
//!
//! And then the half that is not a closed form: **what the rasterisation lost**. A thin plate, meshed and
//! then voxelised at a cell thicker than the plate, is the failure this crate exists to make visible — it
//! rasterises to nothing, or to a sheet one cell thick, and neither reports itself as an error anywhere
//! downstream.

use glam::DVec3;
use pantometry_shape::{Mesh, Triangle, Voxels};
use pantometry_units::{Length, Volume};

/// A closed box from `low` to `high`, in metres, wound outward.
///
/// Written out rather than loaded, so the tests do not depend on a file — and because a box is the one
/// shape whose every measurement here has an exact answer.
fn box_mesh(low: DVec3, high: DVec3) -> Mesh {
    let v = |x: f64, y: f64, z: f64| DVec3::new(x, y, z);
    let (l, h) = (low, high);
    let corner = [
        v(l.x, l.y, l.z),
        v(h.x, l.y, l.z),
        v(h.x, h.y, l.z),
        v(l.x, h.y, l.z),
        v(l.x, l.y, h.z),
        v(h.x, l.y, h.z),
        v(h.x, h.y, h.z),
        v(l.x, h.y, h.z),
    ];
    // Each face as two triangles, wound counter-clockwise seen from outside.
    let face = [
        [0, 3, 2, 1], // z low, seen from below
        [4, 5, 6, 7], // z high
        [0, 1, 5, 4], // y low
        [3, 7, 6, 2], // y high
        [0, 4, 7, 3], // x low
        [1, 2, 6, 5], // x high
    ];
    let mut triangles = Vec::with_capacity(12);
    for f in face {
        triangles.push(Triangle {
            a: corner[f[0]],
            b: corner[f[1]],
            c: corner[f[2]],
        });
        triangles.push(Triangle {
            a: corner[f[0]],
            b: corner[f[2]],
            c: corner[f[3]],
        });
    }
    Mesh::new(triangles)
}

/// A sphere of radius `r` at the origin, tessellated into `bands × sectors` quads split into triangles.
///
/// Closed by construction: the poles are fans and every other edge is shared by exactly two triangles,
/// and the vertices are computed once per `(band, sector)` so shared ones are bit-identical — which is
/// what `Mesh::is_closed` requires and what an exporter does not always give.
///
/// The poles are written out rather than computed, and that is not tidiness. `sin(PI)` is `1.22e-16` and
/// not zero, so the obvious formula puts the south "pole" on a ring of `sectors` distinct points a tenth
/// of a femtometre across — a surface with a puncture in it, geometrically as well as by bit pattern. The
/// north pole is fine by the same formula and the south one is not, which is exactly the kind of asymmetry
/// that survives an eyeball check of the picture.
fn sphere_mesh(r: f64, bands: usize, sectors: usize) -> Mesh {
    let point = |b: usize, s: usize| {
        if b == 0 {
            return DVec3::new(0.0, 0.0, r);
        }
        if b == bands {
            return DVec3::new(0.0, 0.0, -r);
        }
        let phi = std::f64::consts::PI * b as f64 / bands as f64;
        let theta = 2.0 * std::f64::consts::PI * (s % sectors) as f64 / sectors as f64;
        DVec3::new(
            r * phi.sin() * theta.cos(),
            r * phi.sin() * theta.sin(),
            r * phi.cos(),
        )
    };
    let mut triangles = Vec::new();
    for b in 0..bands {
        for s in 0..sectors {
            let (a, bb, c, d) = (
                point(b, s),
                point(b, s + 1),
                point(b + 1, s + 1),
                point(b + 1, s),
            );
            // Wound so the normal points outward, which is what makes `Mesh::volume` positive. Reversing
            // these leaves every edge shared by exactly two triangles, so `is_closed` still passes and
            // only the sign of the volume says which way round it is — see the assertion below.
            if b == 0 {
                triangles.push(Triangle { a, b: d, c });
            } else if b + 1 == bands {
                triangles.push(Triangle { a, b: c, c: bb });
            } else {
                triangles.push(Triangle { a, b: c, c: bb });
                triangles.push(Triangle { a, b: d, c });
            }
        }
    }
    Mesh::new(triangles)
}

/// **A box's mesh volume is exactly `w·h·d`, and the reader round-trips through STL.**
///
/// The divergence-theorem sum has no sampling in it, so the only thing between it and the exact answer is
/// floating point. That makes it the check on the *reader*: a coordinate read into the wrong slot, the
/// millimetre conversion applied twice or not at all, or a face wound inward all fail it and none of them
/// would fail a tolerance-based check on a rasterisation.
///
/// Both STL flavours, because the binary and ASCII paths share no code and the format's own
/// binary-or-text test is famously done wrong — on the word `solid`, which a binary header may contain.
#[test]
fn a_box_measures_exactly_and_survives_both_stl_flavours() {
    // Not 30 x 20 x 10. Those are exact in `f32`, so the round-trip below came out bit-perfect and the
    // tolerance guarding it measured **zero** — a check on a configuration with none of the effect it
    // names in it. A tenth of a millimetre is not a dyadic rational and the f32 error is real: 4.3e-8.
    let (w, h, d) = (0.0301, 0.0203, 0.0107);
    let mesh = box_mesh(DVec3::ZERO, DVec3::new(w, h, d));
    let want = w * h * d;

    assert!(mesh.is_closed(), "a box of twelve triangles is closed");
    assert!(
        (mesh.volume().to_si() / want - 1.0).abs() < 1e-15,
        "volume {} against {want}",
        mesh.volume().to_si()
    );
    let faces = 2.0 * (w * h + w * d + h * d);
    assert!(
        (mesh.area().to_si() / faces - 1.0).abs() < 1e-15,
        "area {} against {faces}",
        mesh.area().to_si()
    );

    // Binary: the 84-byte header and count, then 50 bytes a facet. Millimetres, as every CAD tool
    // writes and as `from_stl` documents.
    let mut binary = vec![0u8; 84];
    binary[80..84].copy_from_slice(&(mesh.triangles().len() as u32).to_le_bytes());
    for t in mesh.triangles() {
        binary.extend_from_slice(&[0u8; 12]); // the normal, which the reader ignores
        for v in [t.a, t.b, t.c] {
            for c in [v.x, v.y, v.z] {
                binary.extend_from_slice(&((c * 1e3) as f32).to_le_bytes());
            }
        }
        binary.extend_from_slice(&[0u8; 2]);
    }
    let read = Mesh::from_stl(&binary).expect("binary STL parses");
    assert_eq!(read.triangles().len(), 12);
    // `f32` in the file, so this is the format's precision and not the reader's.
    let f32_error = (read.volume().to_si() / want - 1.0).abs();
    println!("  the binary round-trip costs {f32_error:.3e}, which is f32 and not the reader");
    assert!(
        f32_error > 1e-9 && f32_error < 1e-6,
        "the file stores f32, so this should be near 1e-7 -- and *nonzero*, or the coordinates chosen are dyadic and the check has no content: {f32_error:e}"
    );

    let mut ascii = String::from("solid box\n");
    for t in mesh.triangles() {
        ascii.push_str("  facet normal 0 0 0\n    outer loop\n");
        for v in [t.a, t.b, t.c] {
            ascii.push_str(&format!(
                "      vertex {} {} {}\n",
                v.x * 1e3,
                v.y * 1e3,
                v.z * 1e3
            ));
        }
        ascii.push_str("    endloop\n  endfacet\n");
    }
    ascii.push_str("endsolid box\n");
    let read = Mesh::from_stl(ascii.as_bytes()).expect("ASCII STL parses");
    assert_eq!(read.triangles().len(), 12);
    // The ASCII writer uses `{}`, which round-trips an f64 exactly, so this one legitimately is zero.
    assert!(
        (read.volume().to_si() / want - 1.0).abs() < 1e-15,
        "ascii volume {} against {want}",
        read.volume().to_si()
    );
}

/// **A cell-aligned box rasterises exactly, at every resolution.**
///
/// The one case where the grid can hold the shape with nothing left over, so the volume error is zero
/// rather than small — and that distinction is the point. An `O(dx)` rasteriser and an exact one look the
/// same at a single resolution on a curved shape; on this one they do not.
///
/// It also pins the grid's bookkeeping. The box is 30×20×10 mm at 2 mm, so 15×10×5 = 750 cells should be
/// inside, and that is asserted as a count rather than as a volume: a volume can come out right with the
/// cells in the wrong place.
#[test]
fn a_cell_aligned_box_rasterises_with_nothing_lost() {
    let mesh = box_mesh(DVec3::ZERO, DVec3::new(0.030, 0.020, 0.010));
    for (cell_mm, cells) in [(2.0, (15, 10, 5)), (1.0, (30, 20, 10)), (0.5, (60, 40, 20))] {
        let want_cells = cells.0 * cells.1 * cells.2;
        let voxels = Voxels::of(&mesh, Length::mm(cell_mm)).expect("a box is closed");
        let loss = voxels.loss();
        println!(
            "  {cell_mm} mm: {} cells of {want_cells}, volume error {:.2e}, thin runs {}, ambiguous {}",
            voxels.filled(),
            loss.volume_error,
            loss.thin_runs,
            loss.ambiguous_rows
        );
        assert_eq!(
            voxels.filled(),
            want_cells,
            "{cell_mm} mm: an aligned box has exactly this many cells in it"
        );
        assert!(
            loss.volume_error.abs() < 1e-12,
            "{cell_mm} mm: nothing should be lost, volume error is {:.3e}",
            loss.volume_error
        );
        assert_eq!(loss.ambiguous_rows, 0, "{cell_mm} mm: no degenerate rows");
        assert!(loss.is_clean(), "{cell_mm} mm: {loss:?}");

        // And the boundary layer is countable here rather than estimated: the cells with no face on
        // the outside are the box shrunk by one cell each way, so the fraction is exact and there is
        // no tolerance to earn. It is the sharpest check on `boundary_fraction` available -- a layer
        // one cell too thick or one cell too thin fails it outright.
        let (nx, ny, nz) = (cells.0 as f64, cells.1 as f64, cells.2 as f64);
        let interior = (nx - 2.0) * (ny - 2.0) * (nz - 2.0);
        let want = 1.0 - interior / (nx * ny * nz);
        println!(
            "    boundary {:.6} against an exact {want:.6}",
            loss.boundary_fraction
        );
        assert!(
            (loss.boundary_fraction - want).abs() < 1e-12,
            "{cell_mm} mm: the exposed cells are the box less its interior, {:.6} against {want:.6}",
            loss.boundary_fraction
        );
    }
}

/// **A sphere's boundary layer is first order with a bracketed coefficient, and its volume error is not
/// a convergent quantity at all.**
///
/// The second half is the finding, and it is the opposite of what the obvious test would have asserted.
/// The signed volume error cancels — cells the surface bulges out of against cells it cuts into — and what
/// is left is a lattice-point count, erratic by nature. **A refinement can make it worse**, and the sweep
/// below asserts that it does rather than only saying so: on powers of two the sequence looks tidy, and it
/// is between them that it is caught. Three points would fit any convergence order you liked, which is why
/// none is claimed.
///
/// Sliding the mesh a third of a cell off the grid changes none of it, so it is not an alignment artefact
/// either.
///
/// What does behave is the **boundary layer**: the volume in cells with a face on the outside. That is a
/// surface area rather than a cancellation, so it is `A·dx/V` — first order, with a coefficient. Both are
/// asserted, and the coefficient is asserted as a *bracket* that column-counting proves rather than as
/// the 0.82 that is measured, because 0.82 is a property of the staircase that nothing here derives.
///
/// The bracket, for a convex body: a filled cell is exposed only if it is the last filled cell in some
/// column along ±x, ±y or ±z. Columns along `z` that meet the body number `πR²/dx²`, the sphere's shadow,
/// and each has two ends — so the `z` ends alone give `2πR²/dx²` cells and all three axes give at most
/// `6πR²/dx²`. Against `A = 4πR²` that is
///
/// ```text
///   0.5 ≤ boundary volume / (A · dx) ≤ 1.5
/// ```
///
/// Everything is measured against the **mesh's** own volume, not `4πR³/3`. A tessellated sphere is a
/// polyhedron holding 0.10% less than the sphere it approximates; comparing to the analytic value would
/// fold the tessellation's error into the grid's and report a rasterisation error that is partly neither.
#[test]
fn a_sphere_has_a_first_order_boundary_and_an_erratic_volume_error() {
    let r = 0.010;
    let mesh = sphere_mesh(r, 64, 128);
    assert!(mesh.is_closed(), "the tessellation is closed");
    let meshed = mesh.volume().to_si();
    let analytic = 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
    println!(
        "  the tessellation holds {:.4}% less than the sphere it approximates",
        (1.0 - meshed / analytic) * 100.0
    );
    assert!(
        meshed < analytic && meshed / analytic > 0.99,
        "a polyhedron inscribed in a sphere holds a little less: {meshed:e} against {analytic:e}"
    );

    let mut boundary = Vec::new();
    for cell_mm in [2.0, 1.0, 0.5, 0.25] {
        let voxels = Voxels::of(&mesh, Length::mm(cell_mm)).expect("closed");
        let loss = voxels.loss();
        // boundary *volume* over A·dx, which is the bracketed quantity. The fraction is of V, so
        // multiplying by V/(A·dx) is the same as dividing by the sphere's own 3·dx/R.
        let coefficient = loss.boundary_fraction / (3.0 * cell_mm * 1e-3 / r);
        println!(
            "  {cell_mm} mm: volume error {:+.4}%, boundary {:.4} of the volume, \
             coefficient {coefficient:.4}",
            loss.volume_error * 100.0,
            loss.boundary_fraction
        );

        // The bracket column counting proves, and it is [0.25, 1.5] rather than the [0.5, 1.5] an
        // earlier version asserted. The lower half was wrong: a column with exactly *one* filled cell
        // has one end, not two, so the floor is the column count and not twice it.
        assert!(
            (0.25..=1.5).contains(&coefficient),
            "{cell_mm} mm: column counting brackets this in [0.25, 1.5], is {coefficient:.4}"
        );
        // And the measured value, which is what actually holds the connectivity choice down. See
        // below for why the derived bracket cannot.
        if cell_mm <= 0.5 {
            assert!(
                (0.75..=0.90).contains(&coefficient),
                "{cell_mm} mm: six-connectivity settles near 0.82 and twenty-six near 1.49; this is \
                 {coefficient:.4}, which is neither"
            );
        }
        boundary.push((cell_mm, loss.boundary_fraction, coefficient));
    }

    // **The bracket above cannot see the choice this crate documents most carefully.** `boundary_share`
    // counts a cell as exposed if one of its **six** face neighbours is outside, not one of its
    // twenty-six — a deliberate call, because a cell touching the outside only at a corner exchanges
    // nothing through a seven-point stencil. Swap six for twenty-six and the coefficient goes to
    // 1.11 / 1.31 / 1.41 / 1.45, all of which sit inside [0.25, 1.5]. The derived bracket spans a factor
    // of six and the two candidate definitions differ by 1.8, so only the measured value separates them,
    // which is why the band above is there and is labelled as measured rather than derived.
    //
    // The box test above cannot see it either, and for a reason worth writing down: on a solid
    // rectangular block the six- and twenty-six-exposed sets are *identical*. The shape with the exact
    // answer is the shape on which the question does not arise.

    // First order. The rate is exactly `2·c(dx)/c(dx/2)` by the definition of `c`, so the band is a band
    // on the coefficient's drift and nothing else — and the drift is measured at no more than 6% per
    // halving over this range, giving 2 × [0.94, 1.06]. Widened to [1.85, 2.15] for the margin.
    // No pair is skipped: all three are inside, and the coarsest (1.897) is nearer 2 than the next.
    for (n, pair) in boundary.windows(2).enumerate() {
        let rate = pair[0].1 / pair[1].1;
        let drift = pair[1].2 / pair[0].2;
        println!(
            "    {} mm to {} mm: layer thinned {rate:.3}x, coefficient drifted {:+.2}%",
            pair[0].0,
            pair[1].0,
            (drift - 1.0) * 100.0
        );
        assert!(
            (1.85..=2.15).contains(&rate),
            "pair {n}: a one-cell layer over a fixed area halves when the cell halves, this went {rate:.3}x"
        );
        assert!(
            (drift - 1.0).abs() < 0.06,
            "pair {n}: the band above is derived from this drift being under 6%, and it is {:.2}%",
            (drift - 1.0) * 100.0
        );
    }
    assert!(
        boundary[0].1 > 0.4,
        "2 mm on a sphere of radius 10 mm leaves 43% of the volume in boundary cells, and a caller \
         choosing that cell size should be told so: {:.4}",
        boundary[0].1
    );

    // **Where the rows were sampled, checked directly.** Everything above is an aggregate, and an
    // aggregate is blind to the sampler being off centre: moving it a quarter cell changes the volume
    // error at 2 mm from +5.5% to +1.7% and leaves every assertion in this file standing. A cell is
    // filled if and only if its *centre* is inside, so on a sphere the filled set is exactly the cells
    // whose centres are within `R` — and the tessellation only ever pulls the surface *in*, never out.
    let voxels = Voxels::of(&mesh, Length::mm(0.5)).expect("closed");
    let (nx, ny, nz) = voxels.counts();
    let cell = voxels.cell().to_si();
    let origin = voxels.origin().to_si();
    let mut furthest_inside: f64 = 0.0;
    let mut nearest_outside = f64::INFINITY;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let centre = origin
                    + glam::DVec3::new(i as f64 + 0.5, j as f64 + 0.5, k as f64 + 0.5) * cell;
                let d = centre.length();
                if voxels.contains(i, j, k) {
                    furthest_inside = furthest_inside.max(d);
                } else {
                    nearest_outside = nearest_outside.min(d);
                }
            }
        }
    }
    // The inscribed radius of a 64-band tessellation: the facets chord across, so the surface dips to
    // `R·cos(π/64)` between vertices and a centre outside that may legitimately be excluded.
    let inscribed = r * (std::f64::consts::PI / 64.0).cos();
    println!(
        "  the furthest filled centre is at {:.5} mm and the nearest empty one at {:.5} mm, \
         between the inscribed {:.5} and R = {:.5}",
        furthest_inside * 1e3,
        nearest_outside * 1e3,
        inscribed * 1e3,
        r * 1e3
    );
    assert!(
        furthest_inside <= r,
        "no cell centre outside the sphere may be filled, and one at {:.6} mm is",
        furthest_inside * 1e3
    );
    assert!(
        nearest_outside >= inscribed,
        "no cell centre inside the tessellation's inscribed radius may be empty, and one at {:.6} mm is",
        nearest_outside * 1e3
    );

    // **`small_triangles`, which nothing else looks at.** It is not evidence that anything was lost —
    // a large flat face tessellated into slivers loses nothing — but a mesh whose facets are mostly
    // below the cell's own area is a mesh with detail the grid cannot hold, and that is worth being
    // told before the run rather than after. Here it is the whole tessellation at a coarse cell and a
    // shrinking fraction at a fine one, which is the behaviour that makes it readable.
    let (coarse_count, fine_count) = (
        Voxels::of(&mesh, Length::mm(2.0))
            .expect("closed")
            .loss()
            .small_triangles,
        Voxels::of(&mesh, Length::mm(0.2))
            .expect("closed")
            .loss()
            .small_triangles,
    );
    println!(
        "  of {} facets, {coarse_count} are under a 2 mm cell face and {fine_count} under a 0.2 mm one",
        mesh.triangles().len()
    );
    assert_eq!(
        coarse_count,
        mesh.triangles().len(),
        "every facet of a 64-band sphere of radius 10 mm is smaller than a 2 mm square"
    );
    assert!(
        fine_count > 0 && fine_count < mesh.triangles().len(),
        "and at 0.2 mm only some of them are, which is what says the count tracks the cell: {fine_count}"
    );

    // The erraticism, swept where it shows. A coarser tessellation, because this is a property of the
    // lattice and not of the triangles -- the boundary fractions come out identical to four decimals on
    // both meshes, which is the check that says so.
    println!("  and the volume error over a sweep that is not powers of two:");
    let coarse = sphere_mesh(r, 32, 64);
    let mut errors = Vec::new();
    for cell_mm in [3.0, 2.5, 2.0, 1.5, 1.25, 1.0, 0.8, 0.625, 0.5] {
        let loss = Voxels::of(&coarse, Length::mm(cell_mm))
            .expect("closed")
            .loss();
        println!(
            "    {cell_mm:5} mm: {:+.4}%, boundary {:.5}",
            loss.volume_error * 100.0,
            loss.boundary_fraction
        );
        errors.push((cell_mm, loss.volume_error));
    }

    let worsened: Vec<_> = errors
        .windows(2)
        .filter(|p| p[1].1.abs() > p[0].1.abs())
        .map(|p| (p[0].0, p[1].0))
        .collect();
    println!("    refining made it worse at {worsened:?}");
    assert!(
        !worsened.is_empty(),
        "the whole point of this test is that the signed volume error is not monotone under \
         refinement; if it has become monotone, either the rasteriser changed or this sweep no longer \
         samples where it is not, and the claim in the documentation above needs re-earning"
    );

    // **The cancellation is real but it is not symmetric, and the number says so.** The obvious next
    // assertion is that the signed error has no bias — mean small against its own scatter — and it is
    // false: measured over these nine, the mean is `+2.56%` against an RMS of `4.19%`, a ratio of 0.61,
    // and over the finer half alone it is 0.89. This rasteriser **over-fills a convex body**, and a bias
    // check written with any threshold under that would be a tolerance chosen to make a test pass.
    //
    // So it is measured and printed and not asserted. The thing a bias check would have been aimed at —
    // a surface placed a fraction of a cell too far out — is caught directly and sharply by the
    // furthest-filled-centre check above, which needs no statistics.
    let signed: Vec<f64> = errors.iter().map(|e| e.1).collect();
    let mean = signed.iter().sum::<f64>() / signed.len() as f64;
    let rms = (signed.iter().map(|e| e * e).sum::<f64>() / signed.len() as f64).sqrt();
    println!(
        "    mean {:+.4}% against an RMS of {:.4}%, a ratio of {:.3} — biased, not centred",
        mean * 100.0,
        rms * 100.0,
        mean.abs() / rms
    );
    assert!(
        mean > 0.0 && mean.abs() < rms,
        "the over-fill is a positive bias smaller than the scatter it sits in: mean {:+.4}% against \
         RMS {:.4}%",
        mean * 100.0,
        rms * 100.0
    );
    assert!(
        signed.iter().any(|e| *e < 0.0),
        "and at least one sample still overshoots the other way, so the bias does not swamp the \
         cancellation — this rests on a single point of the nine (1.5 mm) and would retire quietly if \
         the sweep were changed"
    );
}

/// **A cube rasterises exactly, and it is the shape that says the edge case is handled.**
///
/// This is a regression test and the defect it holds down was silent in the worst way. A ray along `+x`
/// that passes through the diagonal edge shared by a face's two triangles hits **both** of them, so the
/// crossings come out `[0, 0, L, L]`: an even count, which parity approves, paired as `(0, 0)` and
/// `(L, L)`, filling nothing. The row comes back empty, `ambiguous_rows` stays at zero, and the loss
/// report says the rasterisation was clean.
///
/// A cube is where that happens every time, because its face diagonal is `z = y` and every cell centre on
/// the diagonal plane is on it — an eight-cell cube loses 64 of its 512 cells, a plane straight through
/// the middle. The 30×20×10 box above never sees it: its diagonal is `z = y/2`, which needs an even `y`
/// centre and a cell-centred grid has none. **The first version of this crate passed every box test with
/// the bug in it**, and what caught it was a cube in a different crate's integration test.
///
/// So: cubes, at cell sizes that divide them, asserted on the cell *count*. Anything less than the whole
/// count means a row was dropped.
#[test]
fn a_cube_keeps_the_rows_its_face_diagonal_runs_through() {
    for (side_mm, cell_mm, n) in [
        (8.0, 1.0, 8),
        (8.0, 2.0, 4),
        (10.0, 0.5, 20),
        (12.0, 3.0, 4),
        (5.0, 1.0, 5),
    ] {
        let side = side_mm * 1e-3;
        let mesh = box_mesh(DVec3::ZERO, DVec3::splat(side));
        let voxels = Voxels::of(&mesh, Length::mm(cell_mm)).expect("a cube is closed");
        let loss = voxels.loss();
        println!(
            "  a {side_mm} mm cube at {cell_mm} mm: {} of {} cells, {} rows retried, {} ambiguous, error {:.2e}",
            voxels.filled(),
            n * n * n,
            loss.retried_rows,
            loss.ambiguous_rows,
            loss.volume_error
        );
        assert_eq!(
            voxels.filled(),
            n * n * n,
            "{side_mm} mm at {cell_mm} mm: rows on the face diagonal were dropped"
        );
        assert_eq!(
            loss.ambiguous_rows, 0,
            "{side_mm} mm at {cell_mm} mm: the fixed perturbations should recover every one of them"
        );
        // And they were *needed*. Without this the retry path is only ever exercised by inference from
        // a cell count, so a retry that silently did nothing would look like a mesh that never asked
        // for one. `ambiguous_rows` cannot serve here: no mesh in this suite has ever produced one.
        assert!(
            loss.retried_rows > 0,
            "{side_mm} mm at {cell_mm} mm: a cube on cell boundaries sends rows through its face diagonals, so some row must have needed a moved ray"
        );
        assert!(
            loss.volume_error.abs() < 1e-12,
            "{side_mm} mm at {cell_mm} mm: {loss:?}"
        );
    }

    // And the plane that was lost is solid. The count above would also pass if the diagonal plane were
    // dropped and an equal number of cells were wrongly added somewhere else, which is not a failure
    // anyone would construct on purpose but is exactly what a count alone cannot exclude.
    let voxels =
        Voxels::of(&box_mesh(DVec3::ZERO, DVec3::splat(0.008)), Length::mm(1.0)).expect("closed");
    let (nx, _, nz) = voxels.counts();
    let mut on_diagonal = 0;
    for k in 0..nz {
        for i in 0..nx {
            if voxels.contains(i, k, k) {
                on_diagonal += 1;
            }
        }
    }
    println!("  and the diagonal plane holds {on_diagonal} cells, not 0");
    assert_eq!(
        on_diagonal, 64,
        "the plane where j equals k is the one the diagonal edge runs through, and it is solid metal"
    );
}

/// **A plate thinner than a cell is reported, not silently dropped.**
///
/// The failure this crate exists for. A 0.4 mm plate voxelised at 2 mm has nothing the grid can hold: the
/// cell centres either miss it entirely, giving a shape that is *not there*, or catch it and give a sheet
/// one cell thick that is five times too thick. Neither makes a solver fail — it runs, it audits clean,
/// and it answers a question about a different object.
///
/// So both outcomes are asserted, and both have to be *visible* in [`Loss`]: an enormous volume error, or
/// a thin run, or both. `is_clean` is what a caller checks and it must be false.
#[test]
fn a_plate_thinner_than_a_cell_is_reported_rather_than_lost() {
    // 40 x 40 x 0.4 mm — a shim, and the sort of thing a real assembly is full of.
    let plate = box_mesh(DVec3::ZERO, DVec3::new(0.040, 0.040, 0.0004));
    let meshed = plate.volume().to_si();

    for cell_mm in [2.0, 1.0] {
        let voxels = Voxels::of(&plate, Length::mm(cell_mm)).expect("closed");
        let loss = voxels.loss();
        // Across every row, not down `j = 0` — that one is the half cell of margin and is always empty.
        let (nx, ny, nz) = voxels.counts();
        let thickness = (0..nz)
            .filter(|k| (0..ny).any(|j| (0..nx).any(|i| voxels.contains(i, j, *k))))
            .count();
        println!(
            "  {cell_mm} mm on a 0.4 mm plate: {} cells, {thickness} thick, volume error {:.1}%, \
             thin runs {}",
            voxels.filled(),
            loss.volume_error * 100.0,
            loss.thin_runs
        );
        assert!(
            !loss.is_clean(),
            "{cell_mm} mm cannot hold a 0.4 mm plate and must say so: {loss:?}"
        );
        // Either it vanished or it is one cell thick. Both are wrong and both are reported; what must
        // not happen is a plausible-looking answer with a clean loss report.
        assert!(
            voxels.filled() == 0 || loss.thin_runs > 0,
            "a plate at {cell_mm} mm is either gone or one cell thick"
        );
    }

    // **Two cells thick, which is the other half of `thin_runs` and was never reached.** The plate above
    // only ever rasterises one cell thick, so dropping `|| run == 2` from the counter changed nothing.
    // A 0.8 mm plate at 0.4 mm is exactly two cells — resolved in the sense that the volume is right,
    // and still not resolved in the sense that matters, because a seven-point stencil has no interior
    // in it and a trilinear element has one element.
    let two = box_mesh(DVec3::ZERO, DVec3::new(0.040, 0.040, 0.0008));
    let voxels = Voxels::of(&two, Length::mm(0.4)).expect("closed");
    println!(
        "  a 0.8 mm plate at 0.4 mm: volume error {:.2e}, thin runs {}",
        voxels.loss().volume_error,
        voxels.loss().thin_runs
    );
    assert!(
        voxels.loss().volume_error.abs() < 1e-12,
        "two cells hold the plate's volume exactly"
    );
    assert!(
        voxels.loss().thin_runs > 0,
        "and two cells is still thin, which is the case the counter's second clause is for"
    );
    assert!(
        !voxels.loss().is_clean(),
        "so the report is not clean even though the volume is exact — which is the whole point of having a thin-run count beside a volume error: {:?}",
        voxels.loss()
    );

    // And at a cell that can hold it, the report comes back clean and the volume is right.
    let voxels = Voxels::of(&plate, Length::mm(0.1)).expect("closed");
    println!(
        "  0.1 mm: volume error {:.3}%, thin runs {}",
        voxels.loss().volume_error * 100.0,
        voxels.loss().thin_runs
    );
    assert!(
        voxels.loss().volume_error.abs() < 1e-9,
        "an aligned plate at a cell that divides it is exact"
    );
    assert!(
        (voxels.volume().to_si() / meshed - 1.0).abs() < 1e-9,
        "and holds the plate's own volume"
    );
}

/// **An open mesh is refused, and a mesh with a one-bit gap is open.**
///
/// Parity has no meaning through a surface with a hole in it: a ray entering and never leaving leaves the
/// row's crossings odd, and filling by any rule at all would be inventing a shape. So `Voxels::of` refuses
/// rather than producing something that renders.
///
/// The second half is the one worth having. STL stores **no topology**, so two triangles share an edge
/// only if their vertices were written identically, and `is_closed` matches the bit patterns. Moving one
/// vertex by a single unit in the last place opens the mesh — and that is the true answer rather than a
/// strict one, because a ray passes through a gap of one bit as readily as through a gap of a millimetre.
#[test]
fn an_open_mesh_is_refused_and_one_bit_opens_it() {
    let mut triangles = box_mesh(DVec3::ZERO, DVec3::new(0.01, 0.01, 0.01))
        .triangles()
        .to_vec();
    let whole = Mesh::new(triangles.clone());
    assert!(whole.is_closed());
    assert!(Voxels::of(&whole, Length::mm(1.0)).is_ok());

    // A missing facet: the obvious hole.
    triangles.pop();
    let holed = Mesh::new(triangles.clone());
    assert!(!holed.is_closed(), "eleven triangles cannot close a box");
    let error = Voxels::of(&holed, Length::mm(1.0)).expect_err("refused");
    assert!(
        error.contains("not closed") && error.contains("parity"),
        "the message should say what cannot be done and why: {error}"
    );

    // And the one-bit version, which is what an exporter actually produces.
    let mut triangles = whole.triangles().to_vec();
    let nudged = f64::from_bits(triangles[0].a.x.to_bits() + 1);
    triangles[0].a.x = nudged;
    assert!(
        !Mesh::new(triangles).is_closed(),
        "a vertex moved by one unit in the last place leaves an edge unshared, and a ray goes through \
         it as readily as through a visible gap"
    );
}

/// **The same mesh and cell size give the same voxels, bit for bit.**
///
/// The workspace promises reproducibility across platforms, optimisation levels and thread counts, and a
/// rasteriser that perturbs a degenerate ray is exactly where that promise could be lost. The offsets are
/// a fixed list tried in order rather than anything sampled, so the answer is a function of the inputs.
///
/// Checked by rasterising twice and comparing every cell, and on a shape whose rays *do* hit vertices —
/// a box aligned to the grid, where every scanline through a face lands on an edge.
#[test]
fn rasterising_twice_gives_the_same_cells() {
    let mesh = box_mesh(DVec3::ZERO, DVec3::new(0.010, 0.010, 0.010));
    for cell_mm in [1.0, 0.7, 0.3] {
        let a = Voxels::of(&mesh, Length::mm(cell_mm)).expect("closed");
        let b = Voxels::of(&mesh, Length::mm(cell_mm)).expect("closed");
        let (nx, ny, nz) = a.counts();
        assert_eq!(a.counts(), b.counts());
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    assert_eq!(
                        a.contains(i, j, k),
                        b.contains(i, j, k),
                        "{cell_mm} mm: cell ({i}, {j}, {k}) differs between two runs"
                    );
                }
            }
        }
        assert_eq!(a.volume(), b.volume());
        assert_eq!(a.loss(), b.loss());
    }
}

/// **A cell size that is not a number, and an empty mesh, are refused.**
#[test]
fn the_mistakes_a_caller_makes_are_refused() {
    let mesh = box_mesh(DVec3::ZERO, DVec3::new(0.01, 0.01, 0.01));
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            Voxels::of(&mesh, Length::from_si(bad)).is_err(),
            "a cell size of {bad} was accepted"
        );
    }
    let empty = Mesh::new(Vec::new());
    assert!(!empty.is_closed(), "nothing is not a closed surface");
    assert_eq!(empty.volume(), Volume::from_si(0.0));
    assert!(
        empty.bounds().is_none(),
        "the bounds of nothing are not a box"
    );
    assert!(Voxels::of(&empty, Length::mm(1.0)).is_err());

    assert!(Mesh::from_stl(b"not an stl at all").is_err());
    assert!(
        Mesh::from_stl(b"solid x\n  facet normal 0 0 0\n    outer loop\n      vertex 1 2 3\n")
            .is_err(),
        "a facet with one vertex is incomplete and should say so"
    );
}
