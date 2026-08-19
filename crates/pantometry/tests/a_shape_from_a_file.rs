//! A designed part becomes a simulation, and the answer is the one the geometry predicts.
//!
//! This is the claim `pantometry-shape` was built to make good: a mesh a person exported can drive the
//! physics with **no change to any domain**, because a domain's `fill` already takes
//! `Fn(usize, usize, usize) -> bool` and a rasterisation already is one. Nothing here reaches inside
//! `Solid3D`; it hands over a closure and asks the same questions it would ask of a hand-built grid.
//!
//! It also lives where it has to. `pantometry-shape` depends on `pantometry-units` and nothing else in the
//! workspace, and no domain depends on it — the facade is the only crate that can see both, so a test
//! joining them belongs here and in neither of them.
//!
//! # The cells outside the part are not empty
//!
//! `Solid3D` has no void. A grid is a box of cells and every one of them is some substance, so the cells
//! the part does not occupy are filled with another material rather than with nothing. That is a real
//! limitation and these tests are written around it: the closed forms below count both substances by
//! volume, which is honest, and is not the same as having simulated a part surrounded by air.

use glam::DVec3;
use pantometry::prelude::*;
use pantometry::shape::{Mesh, Triangle, Voxels};
use pantometry_core::{Domain, Exchange, Schedule};

/// A closed box from `low` to `high`, in metres, wound outward.
fn box_mesh(low: DVec3, high: DVec3) -> Mesh {
    let (l, h) = (low, high);
    let corner = [
        DVec3::new(l.x, l.y, l.z),
        DVec3::new(h.x, l.y, l.z),
        DVec3::new(h.x, h.y, l.z),
        DVec3::new(l.x, h.y, l.z),
        DVec3::new(l.x, l.y, h.z),
        DVec3::new(h.x, l.y, h.z),
        DVec3::new(h.x, h.y, h.z),
        DVec3::new(l.x, h.y, h.z),
    ];
    let face = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [3, 7, 6, 2],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
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

fn named(name: &str) -> Substance {
    Substance::from_name(name).expect("in the catalogue")
}

/// A block whose cells are `part` exactly where `voxels` says and `around` elsewhere.
fn block_from(voxels: &Voxels, part: &Substance, around: &Substance) -> Solid3D {
    let mut block = Solid3D::new(
        "part",
        around.clone(),
        voxels.counts(),
        voxels.cell(),
        Temperature::celsius(20.0),
    );
    // The whole coupling between geometry and physics, in one line, and the same call the
    // hand-built tests make.
    block.fill(part.clone(), |i, j, k| voxels.contains(i, j, k));
    block
}

/// **A rasterised part holds the heat capacity its geometry says it does.**
///
/// The closed form is `ρ_p·c_p·V + ρ_a·c_a·(V_block − V)` with `V` the volume the *mesh* encloses —
/// measured by the divergence theorem, which involves no grid at all — against the capacity the thermal
/// domain reports for the cells the rasterisation actually filled. Nothing is fitted: a voxeliser that
/// filled the wrong region, the wrong number of cells, or filled them with the wrong substance makes
/// these disagree.
///
/// It is a cell-aligned box, so the agreement is **exact** rather than close, and that is the point. On a
/// curved shape the two differ by the boundary layer and the comparison would need a tolerance — which
/// would hide a rasteriser that was systematically one cell out, the error most likely to be there.
#[test]
fn a_rasterised_part_carries_its_own_heat_capacity() {
    // 24 x 16 x 8 mm at 2 mm: 12 x 8 x 4 cells and nothing left over.
    let mesh = box_mesh(DVec3::ZERO, DVec3::new(0.024, 0.016, 0.008));
    let voxels = Voxels::of(&mesh, Length::mm(2.0)).expect("a box is closed");
    assert!(
        voxels.loss().volume_error.abs() < 1e-12,
        "an aligned box loses nothing: {:?}",
        voxels.loss()
    );
    assert_eq!(voxels.filled(), 12 * 8 * 4, "and fills the cells it should");

    let (aluminium, glass) = (named("aluminium"), named("borosilicate"));
    let block = block_from(&voxels, &aluminium, &glass);

    let meshed = mesh.volume();
    let around = Volume::from_si(block.volume().to_si() - meshed.to_si());
    let want = aluminium
        .heat_capacity(meshed)
        .expect("has a specific heat")
        .to_si()
        + glass
            .heat_capacity(around)
            .expect("has a specific heat")
            .to_si();
    let got = block.heat_capacity().to_si();
    println!(
        "  {} cells of aluminium in {} of glass: {got:.6} J/K against {want:.6} J/K",
        voxels.filled(),
        block.volume().to_si() / voxels.cell().to_si().powi(3),
    );
    assert!(
        (got / want - 1.0).abs() < 1e-12,
        "the cells hold what the mesh's own volume says: {got} against {want}"
    );

    // And the cells are in the right *place*, not merely the right number. A part shifted by one cell
    // has the same capacity and is a different object.
    let (nx, ny, nz) = voxels.counts();
    for (i, j, k, want_part) in [
        (nx / 2, ny / 2, nz / 2, true),
        (0, 0, 0, false),
        (nx - 1, ny - 1, nz - 1, false),
    ] {
        let is_part = *block.substance_at(i, j, k) == aluminium;
        assert_eq!(
            is_part,
            want_part,
            "cell ({i}, {j}, {k}) should {} be part of the object",
            if want_part { "" } else { "not" }
        );
    }
}

/// **A hole through the part removes exactly the hole's own volume.**
///
/// The test the previous one cannot do. A solid box and a box with a shaft through it rasterise to
/// different cell sets, and the difference in heat capacity has a closed form — `(ρ_p c_p − ρ_a c_a)`
/// times the shaft's volume — that involves no simulation. A voxeliser that ignored the inner surface,
/// or filled the hole because parity went wrong on a two-shelled mesh, gives back the solid box's number
/// and this fires.
///
/// Two shells in one mesh is the case worth having, because it is where parity earns its keep: a ray
/// through the shaft crosses four surfaces, and filling between sorted *pairs* leaves the middle empty
/// without anything having to know that the inner shell is a hole rather than a second part. The inner
/// shell is wound **inward**, which is the STL convention for a cavity and is also what makes
/// `Mesh::volume` subtract it without being told.
#[test]
fn a_hole_through_the_part_removes_exactly_its_own_volume() {
    let (w, h, d) = (0.024, 0.024, 0.008);
    let outer = box_mesh(DVec3::ZERO, DVec3::new(w, h, d));

    // An 8 x 8 mm shaft through the middle, on cell boundaries, spanning the full thickness so it is a
    // hole and not a cavity. Reversed winding makes it a void.
    let (lo, hi) = (0.008, 0.016);
    let shaft: Vec<Triangle> = box_mesh(DVec3::new(lo, lo, 0.0), DVec3::new(hi, hi, d))
        .triangles()
        .iter()
        .map(|t| Triangle {
            a: t.a,
            b: t.c,
            c: t.b,
        })
        .collect();
    let mut both = outer.triangles().to_vec();
    both.extend(shaft);
    let holed = Mesh::new(both);

    let shaft_volume = (hi - lo) * (hi - lo) * d;
    let want_volume = w * h * d - shaft_volume;
    println!(
        "  the mesh encloses {:.4} mm3; the box less the shaft is {:.4} mm3",
        holed.volume().to_si() * 1e9,
        want_volume * 1e9
    );
    assert!(
        (holed.volume().to_si() / want_volume - 1.0).abs() < 1e-12,
        "an inward-wound inner shell subtracts: {:e} against {want_volume:e}",
        holed.volume().to_si()
    );

    let pierced = Voxels::of(&holed, Length::mm(2.0)).expect("both shells are closed");
    let solid = Voxels::of(&outer, Length::mm(2.0)).expect("closed");
    println!(
        "  {} cells with the shaft against {} without, volume error {:.3e}, ambiguous rows {}",
        pierced.filled(),
        solid.filled(),
        pierced.loss().volume_error,
        pierced.loss().ambiguous_rows
    );
    assert_eq!(
        pierced.counts(),
        solid.counts(),
        "the same bounding box either way, so the two blocks are comparable"
    );
    assert!(
        pierced.loss().volume_error.abs() < 1e-12,
        "cell-aligned inside and out, so the hole costs nothing extra: {:?}",
        pierced.loss()
    );

    // The middle of the shaft must be empty. Without this the volume could come out right with the
    // wrong cells — a wall one cell too thick against a hole one cell too small, say.
    let (nx, ny, nz) = pierced.counts();
    assert!(
        !pierced.contains(nx / 2, ny / 2, nz / 2),
        "the middle of the shaft is a hole, not metal"
    );
    assert!(
        solid.contains(nx / 2, ny / 2, nz / 2),
        "and it is metal when the shaft is not there, which is what makes the pair a measurement"
    );

    let (aluminium, glass) = (named("aluminium"), named("borosilicate"));
    let removed = block_from(&solid, &aluminium, &glass)
        .heat_capacity()
        .to_si()
        - block_from(&pierced, &aluminium, &glass)
            .heat_capacity()
            .to_si();
    let shaft_volume = Volume::from_si(shaft_volume);
    let want = aluminium
        .heat_capacity(shaft_volume)
        .expect("has one")
        .to_si()
        - glass.heat_capacity(shaft_volume).expect("has one").to_si();
    println!("  the shaft removed {removed:.6} J/K; the closed form says {want:.6} J/K");
    assert!(
        (removed / want - 1.0).abs() < 1e-12,
        "the hole costs exactly its own volume: {removed} against {want}"
    );
}

/// **The part that came out of a file steps, audits, and settles where its geometry says.**
///
/// The two above check the geometry. This one checks that what came out is a *simulation*: it runs under
/// the same conservation audit as everything else, and the equilibrium it reaches is one the mesh
/// determines rather than the grid.
///
/// The closed form is the capacity-weighted mean. An insulated block holds its enthalpy, so whatever the
/// heat does on the way, it ends uniform at
///
/// ```text
///   T = (C_part·T_part + C_around·T_around) / (C_part + C_around)
/// ```
///
/// and `C_part` is `ρc` times the **mesh's** volume. A rasterisation one cell out moves that answer, and
/// moves it by more than the tolerance here — the second half of the test measures by how much, so the
/// tolerance is not taken on trust.
#[test]
fn a_part_from_a_mesh_runs_under_the_same_audit() {
    // An 8 mm cube at 1 mm, so equilibrium is reached in a second rather than an hour.
    let mesh = box_mesh(DVec3::ZERO, DVec3::new(0.008, 0.008, 0.008));
    let voxels = Voxels::of(&mesh, Length::mm(1.0)).expect("closed");
    assert_eq!(voxels.filled(), 8 * 8 * 8);
    let (copper, aluminium) = (named("copper"), named("aluminium"));

    let (hot, cold) = (Temperature::celsius(120.0), Temperature::celsius(20.0));
    let mut block = Solid3D::new(
        "cube",
        aluminium.clone(),
        voxels.counts(),
        voxels.cell(),
        cold,
    );
    block.fill(copper.clone(), |i, j, k| voxels.contains(i, j, k));
    let (nx, ny, nz) = voxels.counts();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                if voxels.contains(i, j, k) {
                    block.set_temperature(i, j, k, hot);
                }
            }
        }
    }

    let meshed = mesh.volume();
    let around = Volume::from_si(block.volume().to_si() - meshed.to_si());
    let c_part = copper.heat_capacity(meshed).expect("has one").to_si();
    let c_around = aluminium.heat_capacity(around).expect("has one").to_si();
    let want = (c_part * hot.to_si() + c_around * cold.to_si()) / (c_part + c_around);

    // The tolerance the assertions below are sized against, named once so the teeth check at the end
    // is tied to it rather than to a second hard-coded number.
    const SETTLED_TOLERANCE_K: f64 = 1e-3;

    let dt = Time::s(1e-4);
    assert!(
        block.stability_ratio(dt) < 1.0,
        "the step has to be inside the limit before any of this means anything: {}",
        block.stability_ratio(dt)
    );

    // Through `Simulation`, so the conservation audit is the one that runs. A part read from a file
    // gets no exemption from it.
    let mut world = Simulation::new(Schedule::Multirate).with(block);
    for n in 0..20_000 {
        world.advance(dt).unwrap_or_else(|v| {
            panic!("step {n}: a part read from a mesh audits like any other: {v}")
        });
    }
    let block = world
        .domain_as::<Solid3D>("cube")
        .expect("still there, and still a Solid3D");

    let spread = block.peak_temperature().to_si() - block.coldest_temperature().to_si();
    let settled = block.mean_temperature().to_si();
    println!(
        "  settled at {settled:.5} K against a weighted mean of {want:.5} K, spread {spread:.2e} K"
    );
    assert!(
        spread < 1e-3,
        "two seconds is long enough for an 11 mm block of metal to become uniform: {spread:.3e} K"
    );
    assert!(
        (settled - want).abs() < SETTLED_TOLERANCE_K,
        "an insulated block ends at its capacity-weighted mean: {settled} against {want}"
    );
    // Absolute, and the scale it is absolute against is worth stating rather than leaving implied:
    // the copper started 100 K above the aluminium, so about 71 J actually crossed between them on the
    // way here and the block holds a little over a kilojoule. A nanojoule is 1.4e-11 of the transfer.
    // `absorbed_energy` is measured from the starting state and an isolated block takes in nothing, so
    // the true value is exactly zero and this is a floating-point margin, not a physical one.
    let crossed = c_part * (hot.to_si() - want);
    println!(
        "  {crossed:.2} J crossed inside the block; it took in {:.3e} J",
        block.absorbed_energy().to_si()
    );
    assert!(
        block.absorbed_energy().to_si().abs() < 1e-9,
        "and it took nothing in, being insulated: {} J against {crossed:.2} J that moved inside it",
        block.absorbed_energy().to_si()
    );

    // What the tolerance is worth. One cell of copper is this many kelvin on the answer, so a
    // rasterisation off by a single cell out of 512 fails the assertion above by a wide margin.
    let one_cell = Volume::from_si(voxels.cell().to_si().powi(3));
    let lost = copper.heat_capacity(one_cell).expect("has one").to_si();
    let gained = aluminium.heat_capacity(one_cell).expect("has one").to_si();
    let shifted = ((c_part - lost) * hot.to_si() + (c_around + gained) * cold.to_si())
        / (c_part - lost + c_around + gained);
    println!(
        "  one cell of {} moves the answer {:.4} K, against a tolerance of 0.001 K",
        voxels.filled(),
        (shifted - want).abs()
    );
    assert!(
        (shifted - want).abs() > 50.0 * SETTLED_TOLERANCE_K,
        "the tolerance has to be small compared to the smallest error it must catch, and one cell \
         moves the answer only {:.4} K",
        (shifted - want).abs()
    );
}

/// **Nothing about a domain changed to make this work.**
///
/// The architectural claim, checked rather than asserted. `Solid3D::fill` takes
/// `impl Fn(usize, usize, usize) -> bool`; a hand-written closure and a rasterisation are the same
/// argument to it, and a block built either way is indistinguishable afterwards.
///
/// The comparison is cell by cell and not by capacity, because two different fills can hold the same
/// heat — which is exactly the mistake this is meant to exclude.
#[test]
fn a_rasterised_fill_and_a_handwritten_one_are_the_same_call() {
    let mesh = box_mesh(
        DVec3::new(0.002, 0.002, 0.002),
        DVec3::new(0.014, 0.010, 0.008),
    );
    let voxels = Voxels::of(&mesh, Length::mm(1.0)).expect("closed");
    let (nx, ny, nz) = voxels.counts();
    let (steel, glass) = (named("stainless_304"), named("borosilicate"));

    let from_file = block_from(&voxels, &steel, &glass);

    // The same region written out by hand, in the grid's own indices. The mesh sits from 2 to 14 mm
    // in x and the grid's origin is half a cell below 2 mm, so the part starts at cell 1.
    let origin = voxels.origin().to_si();
    let cell = voxels.cell().to_si();
    let by_hand_start = |low: f64, axis: f64| ((low - axis) / cell).round() as usize;
    let (i0, j0, k0) = (
        by_hand_start(0.002, origin.x),
        by_hand_start(0.002, origin.y),
        by_hand_start(0.002, origin.z),
    );
    let mut by_hand = Solid3D::new(
        "part",
        glass.clone(),
        voxels.counts(),
        voxels.cell(),
        Temperature::celsius(20.0),
    );
    by_hand.fill(steel.clone(), |i, j, k| {
        (i0..i0 + 12).contains(&i) && (j0..j0 + 8).contains(&j) && (k0..k0 + 6).contains(&k)
    });

    let mut same = 0;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                assert_eq!(
                    from_file.substance_at(i, j, k),
                    by_hand.substance_at(i, j, k),
                    "cell ({i}, {j}, {k}) differs between the file and the hand-written fill"
                );
                same += 1;
            }
        }
    }
    println!(
        "  {same} cells agree, and {} of them are the part",
        voxels.filled()
    );
    assert_eq!(
        voxels.filled(),
        12 * 8 * 6,
        "the part is what the mesh says"
    );

    // And they step identically, which is the statement that nothing downstream can tell them apart.
    let dt = Time::s(1e-4);
    let (mut a, mut b) = (from_file, by_hand);
    let mut bus = Exchange::new();
    a.set_temperature(nx / 2, ny / 2, nz / 2, Temperature::celsius(200.0));
    b.set_temperature(nx / 2, ny / 2, nz / 2, Temperature::celsius(200.0));
    let mut t = Time::from_si(0.0);
    for _ in 0..500 {
        a.step(t, dt, &mut bus).expect("stable");
        b.step(t, dt, &mut bus).expect("stable");
        t += dt;
    }
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                assert_eq!(
                    a.temperature_at(i, j, k).to_si().to_bits(),
                    b.temperature_at(i, j, k).to_si().to_bits(),
                    "after 500 steps, cell ({i}, {j}, {k}) differs bit for bit"
                );
            }
        }
    }
}
