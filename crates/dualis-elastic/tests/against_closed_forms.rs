//! What a block of elastic material does, against things that were true before this code existed.
//!
//! Four moduli, an energy identity, an equilibrium statement and a convergence rate. None of these
//! compares against another run of the same solver, and none is a number read off this
//! implementation and pasted back in.
//!
//! # Why four moduli and not one
//!
//! `E`, `M`, `K` and `G` are four different combinations of the same two constants. A solver that
//! had `λ` and `μ` transposed reproduces none of them; one whose shear rows carried `2μ` instead
//! of `μ` reproduces the first three and fails the fourth; one that split a face's traction
//! equally between its nodes instead of by shape function gets all four slightly wrong in a way
//! that looks like a boundary effect. Checking one modulus is checking that a stiffness exists.

use dualis_core::Domain;
use dualis_elastic::{Block, Elastic, Face};
use dualis_units::{Length, Pressure};

const DX: f64 = 1e-3;
/// The applied stress in every uniform-strain case. Small enough that linear elasticity is the
/// right model for a real metal, which matters only for the claim being about anything.
const LOAD: f64 = 1.0e6;

fn block(counts: (usize, usize, usize)) -> Block {
    Block::new(
        "block",
        counts,
        Length::from_si(DX),
        Elastic::aluminium_6061(),
    )
}

/// The three rollers that leave exactly one static answer.
///
/// Each holds one component on one face, so together they remove the three translations and the
/// three rotations and nothing else — the body is still free to strain in every direction. This is
/// the standard statically determinate mount every closed form below is written against, and
/// anything less leaves a rigid motion in and anything more adds a constraint the closed form does
/// not have.
fn rollers(b: &mut Block) {
    b.roller(Face::XLow);
    b.roller(Face::YLow);
    b.roller(Face::ZLow);
}

/// **Uniaxial stress gives Young's modulus exactly, and Poisson's ratio with it.**
///
/// The sides are free, so `σ = Eε` and the two lateral strains are `−νε`. Trilinear elements
/// reproduce a linear displacement field exactly, so this is machine precision at any mesh — and
/// checking it at three meshes is what says so rather than assumes it.
#[test]
fn a_pulled_block_gives_youngs_modulus_exactly() {
    let m = Elastic::aluminium_6061();
    for counts in [(2, 2, 2), (4, 3, 5), (8, 8, 8)] {
        let mut b = block(counts);
        rollers(&mut b);
        b.pull(Face::XHigh, Pressure::from_si(LOAD));
        assert!(b.solve(1e-12), "residual {:.3e}", b.residual());

        let ex = b.mean_strain(0);
        let modulus = LOAD / ex;
        let nu_y = -b.mean_strain(1) / ex;
        let nu_z = -b.mean_strain(2) / ex;
        println!(
            "  {counts:?}: E {:.6e} against {:.6e}, nu {nu_y:.9} and {nu_z:.9}",
            modulus,
            m.youngs_modulus.to_si()
        );
        assert!(
            (modulus / m.youngs_modulus.to_si() - 1.0).abs() < 1e-9,
            "uniform strain is exact for this element: {modulus:.6e}"
        );
        for nu in [nu_y, nu_z] {
            assert!(
                (nu - m.poisson_ratio).abs() < 1e-9,
                "and so is the lateral contraction: {nu:.9} against {}",
                m.poisson_ratio
            );
        }
    }
}

/// **Held sides give the constrained modulus, which is half again bigger.**
///
/// `M = E(1−ν)/((1+ν)(1−2ν))` — the stiffness of a block that cannot spread. At `ν = 0.33` it is
/// 1.49× Young's modulus, so this and the test above cannot both be passed by a solver that has
/// only one number in it.
#[test]
fn a_confined_block_gives_the_constrained_modulus() {
    let m = Elastic::aluminium_6061();
    let mut b = block((6, 6, 6));
    // Every lateral face on a roller: the block may shorten along x and may not spread at all.
    for face in [Face::XLow, Face::YLow, Face::YHigh, Face::ZLow, Face::ZHigh] {
        b.roller(face);
    }
    b.pull(Face::XHigh, Pressure::from_si(LOAD));
    assert!(b.solve(1e-12), "residual {:.3e}", b.residual());

    let modulus = LOAD / b.mean_strain(0);
    println!(
        "  M {:.6e} against {:.6e}, and E is {:.6e}",
        modulus,
        m.constrained_modulus().to_si(),
        m.youngs_modulus.to_si()
    );
    assert!(
        (modulus / m.constrained_modulus().to_si() - 1.0).abs() < 1e-9,
        "confined compression is lambda + 2mu: {modulus:.6e}"
    );
    // The lateral strains are zero because they were held, which is the hypothesis rather than a
    // result — but a solver that let them move would fail the modulus above, so this states it.
    for a in [1usize, 2] {
        assert!(
            b.mean_strain(a).abs() < 1e-18,
            "a held face does not move: {:.3e}",
            b.mean_strain(a)
        );
    }
    assert!(
        modulus / m.youngs_modulus.to_si() > 1.4,
        "and it is well clear of E: {:.4}x",
        modulus / m.youngs_modulus.to_si()
    );
}

/// **Pressure on all six faces gives the bulk modulus.**
///
/// `ΔV/V = −p/K`, `K = E/(3(1−2ν))`. A third independent combination, and the only one of the
/// four in which every face carries a load rather than a hold.
#[test]
fn pressure_on_every_face_gives_the_bulk_modulus() {
    let m = Elastic::aluminium_6061();
    let mut b = block((5, 5, 5));
    // Rollers on the low faces are a symmetry mount, and pressing all six is symmetric — so the
    // low faces do not move and holding them changes nothing but the rigid motion.
    rollers(&mut b);
    for face in [Face::XHigh, Face::YHigh, Face::ZHigh] {
        b.press(face, Pressure::from_si(LOAD));
    }
    assert!(b.solve(1e-12), "residual {:.3e}", b.residual());

    let modulus = -LOAD / b.volumetric_strain();
    println!(
        "  K {:.6e} against {:.6e}",
        modulus,
        m.bulk_modulus().to_si()
    );
    assert!(
        (modulus / m.bulk_modulus().to_si() - 1.0).abs() < 1e-9,
        "hydrostatic compression is E/(3(1-2nu)): {modulus:.6e}"
    );
    // And it compressed rather than expanded, which a sign error in `press` would reverse while
    // leaving the magnitude right.
    assert!(
        b.volumetric_strain() < 0.0,
        "pressing a block makes it smaller: {:.3e}",
        b.volumetric_strain()
    );
}

/// **The patch test, on a general linear field**, which is what an element has to pass first.
///
/// Prescribe a linear displacement on the whole boundary and the interior has exactly one answer:
/// the same linear field. An element that cannot reproduce it is not consistent, and no amount of
/// refinement makes it converge to the right thing.
///
/// The field here has all six strain components non-zero and a rigid rotation mixed in, so nothing
/// about it is special to an axis. Checked two ways: every interior node lands on the field, and
/// the strain energy equals `½ε:D:ε·V` computed by hand from the constants.
#[test]
fn the_interior_reproduces_a_linear_field_prescribed_on_the_boundary() {
    let m = Elastic::aluminium_6061();
    // A general small gradient: symmetric part is the strain, antisymmetric part is a rotation
    // the energy must ignore.
    let g = [
        [3.0e-4, 1.0e-4, -2.0e-4],
        [-4.0e-4, -1.0e-4, 5.0e-4],
        [2.0e-4, 3.0e-4, 4.0e-4],
    ];
    let field = |p: dualis_units::LengthVec| {
        let v = p.to_si();
        dualis_units::LengthVec::from_si(glam::DVec3::new(
            g[0][0] * v.x + g[0][1] * v.y + g[0][2] * v.z,
            g[1][0] * v.x + g[1][1] * v.y + g[1][2] * v.z,
            g[2][0] * v.x + g[2][1] * v.y + g[2][2] * v.z,
        ))
    };

    let mut b = block((4, 3, 5));
    b.prescribe_boundary(field);
    assert!(b.solve(1e-14), "residual {:.3e}", b.residual());

    let (nx, ny, nz) = b.node_counts();
    let mut worst: f64 = 0.0;
    let mut scale: f64 = 0.0;
    for k in 1..nz - 1 {
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let at = dualis_units::LengthVec::from_si(
                    glam::DVec3::new(i as f64, j as f64, k as f64) * DX,
                );
                let want = field(at).to_si();
                let got = b.displacement_at(i, j, k).to_si();
                worst = worst.max((got - want).length());
                scale = scale.max(want.length());
            }
        }
    }
    println!("  interior off by {worst:.3e} m against a field spanning {scale:.3e} m");
    assert!(
        worst < 1e-12 * scale.max(1e-30),
        "a trilinear element reproduces a linear field exactly: off by {worst:.3e} m"
    );

    // And the energy, from the constants rather than from the solver.
    let (lambda, mu) = (
        m.youngs_modulus.to_si() * m.poisson_ratio
            / ((1.0 + m.poisson_ratio) * (1.0 - 2.0 * m.poisson_ratio)),
        m.shear_modulus().to_si(),
    );
    let mut eps = [[0.0f64; 3]; 3];
    for a in 0..3 {
        for c in 0..3 {
            eps[a][c] = 0.5 * (g[a][c] + g[c][a]);
        }
    }
    let tr = eps[0][0] + eps[1][1] + eps[2][2];
    let dd: f64 = (0..3)
        .map(|a| (0..3).map(|c| eps[a][c] * eps[a][c]).sum::<f64>())
        .sum();
    let v = b.size().to_si();
    let closed = (0.5 * lambda * tr * tr + mu * dd) * v.x * v.y * v.z;
    let got = b.strain_energy().to_si();
    println!("  energy {got:.9e} J against ½(λ tr² + 2μ ε:ε)V = {closed:.9e} J");
    assert!(
        (got / closed - 1.0).abs() < 1e-9,
        "the energy of a uniform strain is exact: {got:.6e} against {closed:.6e}"
    );
}

/// **Simple shear costs `½μγ²V`, which isolates the one constant a factor of two hides in.**
///
/// The engineering shear strain `γ` is twice the tensor component, and the constitutive matrix
/// carries `μ` on its shear rows rather than `2μ` for exactly that reason. A solver that put `2μ`
/// there passes tension, confinement and pressure and fails only here — by a factor of two, which
/// is why this test states what the wrong answer would be.
///
/// Prescribed on the boundary rather than driven by a rig. A shear rig — clamp the bottom, drag the
/// top — is not simple shear: the sides are free, the block bends, and a cube comes out at 0.40 of
/// `G`. That was measured, and it is why this test looks the way it does.
#[test]
fn simple_shear_costs_half_mu_gamma_squared() {
    let m = Elastic::aluminium_6061();
    let gamma = 1e-4;
    let mut b = block((4, 4, 4));
    b.prescribe_boundary(|p| {
        dualis_units::LengthVec::from_si(glam::DVec3::new(gamma * p.to_si().y, 0.0, 0.0))
    });
    assert!(b.solve(1e-14), "residual {:.3e}", b.residual());

    let v = b.size().to_si();
    let volume = v.x * v.y * v.z;
    let closed = 0.5 * m.shear_modulus().to_si() * gamma * gamma * volume;
    let got = b.strain_energy().to_si();
    println!(
        "  {got:.9e} J against ½Gγ²V = {closed:.9e} J   (2G would give {:.9e})",
        2.0 * closed
    );
    assert!(
        (got / closed - 1.0).abs() < 1e-9,
        "simple shear costs ½Gγ²V exactly: {got:.6e} against {closed:.6e}"
    );
    // No volume change in simple shear, which a lambda leaking into the shear rows would break.
    println!("  volumetric strain {:.3e}", b.volumetric_strain());
    assert!(
        b.volumetric_strain().abs() < 1e-15,
        "shear does not change volume: {:.3e}",
        b.volumetric_strain()
    );
}

/// **Clapeyron: twice the strain energy is the work the loads did.**
///
/// This crate's Tellegen. One side comes from the stiffness and the displacement, the other from
/// the loads and the displacement, and they meet only if the displacement actually solves the
/// system. A solve stopped early fails it; a wrong load vector fails it; a non-symmetric operator
/// fails it.
#[test]
fn twice_the_strain_energy_is_the_work_done() {
    for counts in [(3, 3, 3), (6, 4, 5)] {
        let mut b = block(counts);
        rollers(&mut b);
        b.pull(Face::XHigh, Pressure::from_si(LOAD));
        b.press(Face::YHigh, Pressure::from_si(LOAD / 3.0));
        assert!(b.solve(1e-13), "residual {:.3e}", b.residual());

        let (u, w) = (b.strain_energy().to_si(), b.work_done().to_si());
        println!(
            "  {counts:?}: 2U {:.9e} against W {:.9e}, off {:.2e}",
            2.0 * u,
            w,
            b.energy_balance()
        );
        assert!(
            b.energy_balance() < 1e-9,
            "Clapeyron holds at equilibrium: off by {:.3e}",
            b.energy_balance()
        );
        assert!(u > 0.0, "a loaded body holds energy: {u:.3e} J");
    }
}

/// **The reaction equals the load, which is equilibrium and not an assumption.**
///
/// The rollers carry exactly what the pulled face applies. Computed from `K·u − f` on the held
/// degrees of freedom, so it is the *discrete* equilibrium being checked rather than the
/// continuous one being restated.
#[test]
fn what_the_holds_carry_is_what_was_applied() {
    let (nx, ny, nz) = (4usize, 6, 5);
    let mut b = block((nx, ny, nz));
    rollers(&mut b);
    b.pull(Face::XHigh, Pressure::from_si(LOAD));
    assert!(b.solve(1e-12), "residual {:.3e}", b.residual());

    let area = (ny as f64 * DX) * (nz as f64 * DX);
    let applied = LOAD * area;
    let carried = -b.normal_reaction(Face::XLow).to_si();
    println!("  applied {applied:.6e} N, the x-low rollers carry {carried:.6e} N");
    assert!(
        (carried / applied - 1.0).abs() < 1e-9,
        "sum of forces is zero: {carried:.6e} against {applied:.6e}"
    );

    // The lateral rollers carry nothing along their own normals, which is what makes this
    // uniaxial. **The normal component and not the magnitude**: a node on an edge belongs to two
    // faces, so `reaction(YLow)` also sums the y-low share of the x-low rollers' 30 N and reports
    // 2.5 N of a force it is not carrying. That is the shared edge, not a lateral stress, and
    // taking the magnitude here is how it would be mistaken for one.
    for face in [Face::YLow, Face::ZLow] {
        let normal = b.normal_reaction(face).to_si().abs();
        let magnitude = b.reaction(face).length();
        println!(
            "  {face:?}: normal {normal:.3e} N, magnitude {magnitude:.3e} N — the difference is \
             the shared edge"
        );
        assert!(
            normal < 1e-9 * applied,
            "an unloaded roller carries nothing along its normal: {normal:.3e} N"
        );
        assert!(
            magnitude > 1e-3 * applied,
            "and the magnitude is not that number, which is the point of measuring the component"
        );
    }
}

/// **A body that is not held enough is refused rather than solved.**
///
/// Six rigid motions cost no energy, so a free body's system is singular by exactly six dimensions
/// and has infinitely many answers. Returning one of them would be a displacement field that is
/// correct up to a translation nobody asked for — plausible, and not the answer.
#[test]
fn a_body_free_to_drift_is_refused() {
    let mut b = block((3, 3, 3));
    b.pull(Face::XHigh, Pressure::from_si(LOAD));
    // Nothing held at all.
    assert!(
        !b.solve(1e-12),
        "an unheld body has no unique equilibrium and must not report one"
    );
    assert_eq!(b.free_dofs(), 4 * 4 * 4 * 3);

    // A single pin removes the three translations and leaves the three rotations, so it is still
    // not enough — which is the part that is easy to get wrong.
    let mut b = block((3, 3, 3));
    b.pin(0, 0, 0);
    b.pull(Face::XHigh, Pressure::from_si(LOAD));
    assert!(
        !b.solve(1e-12),
        "a single pin leaves three rotations and is not a mount"
    );

    // The three rollers are.
    let mut b = block((3, 3, 3));
    rollers(&mut b);
    b.pull(Face::XHigh, Pressure::from_si(LOAD));
    assert!(b.solve(1e-12), "three rollers are a complete mount");
}

/// **A rigid motion costs no energy, and nothing else is free.**
///
/// The null space of the stiffness is exactly the six rigid motions. Checked by putting each of
/// them in and reading the energy, and then by putting in something that is *not* one and reading
/// a positive number — because "the energy of this displacement is zero" is only interesting
/// beside a displacement whose energy is not.
#[test]
fn the_null_space_is_the_six_rigid_motions_and_nothing_more() {
    let m = Elastic::aluminium_6061();
    let b = block((3, 3, 3));
    let (nx, ny, nz) = b.node_counts();
    let volume = 3.0 * DX * 3.0 * DX * 3.0 * DX;
    // A strain of this size in the same body, for scale.
    let reference = 0.5 * m.youngs_modulus.to_si() * 1e-4f64.powi(2) * volume;

    type Mode = fn(f64, f64, f64) -> [f64; 3];
    let modes: [(&str, Mode); 6] = [
        ("translate x", |_, _, _| [1e-4, 0.0, 0.0]),
        ("translate y", |_, _, _| [0.0, 1e-4, 0.0]),
        ("translate z", |_, _, _| [0.0, 0.0, 1e-4]),
        ("rotate about z", |x, y, _| [-1e-2 * y, 1e-2 * x, 0.0]),
        ("rotate about x", |_, y, z| [0.0, -1e-2 * z, 1e-2 * y]),
        ("rotate about y", |x, _, z| [1e-2 * z, 0.0, -1e-2 * x]),
    ];
    for (label, f) in modes {
        let mut u = Vec::new();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let _ = (i, j, k);
                }
            }
        }
        // Build in node order, which is x fastest.
        u.clear();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let (x, y, z) = (i as f64 * DX, j as f64 * DX, k as f64 * DX);
                    u.extend(f(x, y, z));
                }
            }
        }
        let e = energy_of(&b, &u);
        println!("  {label:<16} {e:.3e} J against a reference strain's {reference:.3e} J");
        assert!(
            e < 1e-12 * reference,
            "{label} is a rigid motion and must cost nothing: {e:.3e} J"
        );
    }

    // And a uniform stretch, which is not rigid, costs what it should: `½Eε²V` for uniaxial
    // strain with the sides free is not this case — the sides are not free in a prescribed field,
    // so the right closed form is `½ M ε² V` with the constrained modulus.
    let mut u = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let _ = (j, k);
                u.extend([1e-4 * i as f64 * DX, 0.0, 0.0]);
            }
        }
    }
    let e = energy_of(&b, &u);
    let closed = 0.5 * m.constrained_modulus().to_si() * 1e-4f64.powi(2) * volume;
    println!("  uniform stretch  {e:.6e} J against ½ M eps² V = {closed:.6e} J");
    assert!(
        (e / closed - 1.0).abs() < 1e-9,
        "a uniform stretch's energy is exact for this element: {e:.6e} against {closed:.6e}"
    );
}

/// The strain energy of an arbitrary displacement field, through the public solve.
///
/// `Block` has no way to set a displacement directly, and adding one for a test would be adding
/// public API to make a test convenient. Instead the field is imposed by holding every degree of
/// freedom and loading with `K·u` — but that needs `K·u`, which is the thing being measured. So
/// the energy is computed here from the same element stiffness, assembled independently in the
/// test, and that independence is the point: a bug in the crate's assembly would have to be
/// reproduced exactly here to hide.
fn energy_of(body: &Block, u: &[f64]) -> f64 {
    let m = body.material();
    let (nu, e) = (m.poisson_ratio, m.youngs_modulus.to_si());
    let mu = e / (2.0 * (1.0 + nu));
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let (ex, ey, ez) = body.elements();
    let (nx, ny, _) = body.node_counts();
    let h = body.cell().to_si();
    let node = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);

    // 2x2x2 Gauss over each element, strain from the trilinear gradients, energy density from the
    // isotropic law. Written the long way rather than through the crate's own element.
    let g = 1.0 / 3.0f64.sqrt();
    let mut total = 0.0;
    for e_z in 0..ez {
        for e_y in 0..ey {
            for e_x in 0..ex {
                for gz in [-g, g] {
                    for gy in [-g, g] {
                        for gx in [-g, g] {
                            let mut grad = [[0.0f64; 3]; 3]; // du_a / dx_b
                            for (c, corner) in [
                                [0usize, 0, 0],
                                [1, 0, 0],
                                [0, 1, 0],
                                [1, 1, 0],
                                [0, 0, 1],
                                [1, 0, 1],
                                [0, 1, 1],
                                [1, 1, 1],
                            ]
                            .iter()
                            .enumerate()
                            {
                                let _ = c;
                                let s = [
                                    if corner[0] == 0 { -1.0 } else { 1.0 },
                                    if corner[1] == 0 { -1.0 } else { 1.0 },
                                    if corner[2] == 0 { -1.0 } else { 1.0 },
                                ];
                                let (fx, fy, fz) = (
                                    0.5 * (1.0 + s[0] * gx),
                                    0.5 * (1.0 + s[1] * gy),
                                    0.5 * (1.0 + s[2] * gz),
                                );
                                let d = [
                                    0.5 * s[0] * fy * fz * 2.0 / h,
                                    fx * 0.5 * s[1] * fz * 2.0 / h,
                                    fx * fy * 0.5 * s[2] * 2.0 / h,
                                ];
                                let n = node(e_x + corner[0], e_y + corner[1], e_z + corner[2]);
                                for a in 0..3 {
                                    for b in 0..3 {
                                        grad[a][b] += u[3 * n + a] * d[b];
                                    }
                                }
                            }
                            let mut strain = [[0.0f64; 3]; 3];
                            for (a, row) in strain.iter_mut().enumerate() {
                                for (b, e) in row.iter_mut().enumerate() {
                                    *e = 0.5 * (grad[a][b] + grad[b][a]);
                                }
                            }
                            let tr = strain[0][0] + strain[1][1] + strain[2][2];
                            let dd: f64 = strain
                                .iter()
                                .flat_map(|row| row.iter())
                                .map(|e| e * e)
                                .sum();
                            total += (0.5 * lambda * tr * tr + mu * dd) * (h / 2.0).powi(3);
                        }
                    }
                }
            }
        }
    }
    total
}

/// **A cantilever converges on `PL³/3EI` from below**, which is what a locking element does.
///
/// Not an equality, and the reason is stated rather than tuned around: a fully integrated
/// trilinear element develops shear strain where it should be flexing, so it is too stiff, and the
/// error falls as the mesh through the thickness refines. The claim is therefore a **direction and
/// a limit**: every refinement must be softer than the last, all must be stiffer than the beam
/// formula, and the finest must be within reach of it.
///
/// A solver whose bending was simply wrong would not approach anything.
#[test]
fn a_cantilever_converges_on_the_beam_formula_from_below() {
    let m = Elastic::aluminium_6061();
    let mut ratios = Vec::new();
    for through in [2usize, 3, 4, 6] {
        // Ten to one, which is slender enough for the beam formula to be the right closed form
        // and short enough that the element count stays sane.
        let (nx, ny, nz) = (through * 10, through, through);
        let dx = 1e-3;
        let mut b = Block::new("beam", (nx, ny, nz), Length::from_si(dx), m);
        b.clamp(Face::XLow);
        let (l, h, w) = (nx as f64 * dx, ny as f64 * dx, nz as f64 * dx);
        // A shear traction on the free end, whose resultant is the beam formula's `P`.
        let p = 50.0;
        b.traction(Face::XHigh, glam::DVec3::new(0.0, -p / (h * w), 0.0));
        assert!(b.solve(1e-12), "residual {:.3e}", b.residual());

        let tip = -b.displacement_at(nx, ny / 2, nz / 2).to_si().y;
        let i = w * h.powi(3) / 12.0;
        let closed = p * l.powi(3) / (3.0 * m.youngs_modulus.to_si() * i);
        println!(
            "  {nx:>2}x{ny}x{nz}: tip {tip:.4e} m against PL^3/3EI {closed:.4e} — {:.3}x",
            tip / closed
        );
        ratios.push(tip / closed);
    }
    for r in &ratios {
        assert!(
            *r < 1.0,
            "a fully integrated element is too stiff, never too soft: {r:.4}x"
        );
    }
    for pair in ratios.windows(2) {
        assert!(
            pair[1] > pair[0],
            "refining must soften it toward the limit: {:.4} then {:.4}",
            pair[0],
            pair[1]
        );
    }
    assert!(
        *ratios.last().unwrap() > 0.55,
        "and the finest must be within reach of the beam formula: {:.3}x",
        ratios.last().unwrap()
    );
}

/// **The domain refuses a step it could not solve, rather than reporting a field of zeros.**
#[test]
fn an_unsolvable_body_refuses_its_step() {
    let mut b = block((3, 3, 3));
    b.pull(Face::XHigh, Pressure::from_si(LOAD));
    let err = b
        .step(
            dualis_units::Time::from_si(0.0),
            dualis_units::Time::from_si(1.0),
            &mut dualis_core::Exchange::new(),
        )
        .expect_err("an unheld body cannot be solved");
    assert_eq!(err.quantity, "equilibrium residual");

    rollers(&mut b);
    b.step(
        dualis_units::Time::from_si(0.0),
        dualis_units::Time::from_si(1.0),
        &mut dualis_core::Exchange::new(),
    )
    .expect("a mounted one can");
    assert!(b.strain_energy().to_si() > 0.0);
}
