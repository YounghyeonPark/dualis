//! The trilinear hexahedron, and why the stiffness is assembled from an energy.
//!
//! # Why not finite differences on the Navier–Cauchy equation
//!
//! `(λ+μ)∇(∇·u) + μ∇²u = 0` differenced directly is the obvious thing and it is a trap. The
//! operator that comes out is not symmetric at a boundary unless every one-sided difference is
//! chosen to make it so, and conjugate gradients needs symmetry — not as an optimisation but for
//! the method to mean anything. Worse, the `∇(∇·u)` term on a collocated grid admits a
//! checkerboard displacement with zero energy, so the solve is free to add one and nothing in the
//! residual objects.
//!
//! Assembling from the strain energy avoids both by construction. The operator is the Hessian of
//!
//! ```text
//!   U = ∫ ½ ( λ (tr ε)² + 2μ ε:ε ) dV
//! ```
//!
//! which is symmetric because a Hessian is, and positive semi-definite because `U ≥ 0` with
//! equality only on rigid motions. Those six rigid modes are the *only* null space, which is a
//! statement the tests check rather than assume.
//!
//! # What this element is good at, and what it is not
//!
//! Trilinear shape functions reproduce any **linear** displacement field exactly, so uniform
//! strain — tension, compression, hydrostatic pressure, simple shear — comes out to machine
//! precision at any mesh size. That is what makes four separate moduli checkable as equalities
//! rather than as convergences.
//!
//! It is poor at **bending**. Fully integrated, it develops spurious shear strain when it should
//! be flexing, and a slender beam comes out too stiff — shear locking. That is not a defect being
//! hidden: it is why the cantilever test here is a convergence toward `PL³/3EI` from below rather
//! than an equality, and why anything asking about bending wants more elements through the
//! thickness than it would need for tension.
//!
//! Full 2×2×2 integration rather than reduced, deliberately. Reduced integration cures the locking
//! and buys hourglass modes — zero-energy patterns that are not rigid motions — which is trading a
//! stiffness error for a *singularity*, and a singular operator with a spurious null space is the
//! failure this file's first paragraph is about.

/// Degrees of freedom in one element: eight nodes, three each.
pub(crate) const DOF: usize = 24;

/// The local node order: `a + 2b + 4c` for the corner at `(a, b, c)`.
///
/// Stated rather than implied, because every index in the assembly depends on it and a
/// transposition here produces a stiffness that is still symmetric, still positive definite, and
/// wrong in a way no invariant catches.
pub(crate) const CORNERS: [[usize; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
];

/// The 24×24 stiffness of one cubic element of side `h`, for Lamé constants `lambda` and `mu`.
///
/// Row-major, `DOF*DOF`. Computed once per body: every element is the same cube of the same
/// material, so assembling it per element would be the same arithmetic `n³` times.
pub(crate) fn stiffness(h: f64, lambda: f64, mu: f64) -> Vec<f64> {
    let mut k = vec![0.0; DOF * DOF];
    // 2×2×2 Gauss. The points are at ±1/√3 in the reference cube and each carries unit weight.
    let g = 1.0 / 3.0f64.sqrt();
    let det_j = (h / 2.0).powi(3);
    // The reference cube spans ±1, so a derivative in it is `2/h` of one in space.
    let scale = 2.0 / h;

    for gz in [-g, g] {
        for gy in [-g, g] {
            for gx in [-g, g] {
                // Shape-function gradients in space, one row of three per node.
                let mut grad = [[0.0f64; 3]; 8];
                for (n, c) in CORNERS.iter().enumerate() {
                    let s = [
                        if c[0] == 0 { -1.0 } else { 1.0 },
                        if c[1] == 0 { -1.0 } else { 1.0 },
                        if c[2] == 0 { -1.0 } else { 1.0 },
                    ];
                    let (fx, fy, fz) = (
                        0.5 * (1.0 + s[0] * gx),
                        0.5 * (1.0 + s[1] * gy),
                        0.5 * (1.0 + s[2] * gz),
                    );
                    grad[n] = [
                        0.5 * s[0] * fy * fz * scale,
                        fx * 0.5 * s[1] * fz * scale,
                        fx * fy * 0.5 * s[2] * scale,
                    ];
                }

                // B, six strain rows by 24 columns, in Voigt order
                // `[εxx, εyy, εzz, γyz, γzx, γxy]`.
                let mut b = [[0.0f64; DOF]; 6];
                for (n, g) in grad.iter().enumerate() {
                    let (dx, dy, dz) = (g[0], g[1], g[2]);
                    let (cx, cy, cz) = (3 * n, 3 * n + 1, 3 * n + 2);
                    b[0][cx] = dx;
                    b[1][cy] = dy;
                    b[2][cz] = dz;
                    b[3][cy] = dz;
                    b[3][cz] = dy;
                    b[4][cx] = dz;
                    b[4][cz] = dx;
                    b[5][cx] = dy;
                    b[5][cy] = dx;
                }

                // D·B, then Bᵀ(D·B). `D` is block diagonal: a 3×3 of `λ` with `2μ` on its
                // diagonal, and `μ` on each shear row. Written out rather than as a matrix,
                // because the shear rows carry `μ` and not `2μ` and that factor is the classic
                // place to lose a two.
                let mut db = [[0.0f64; DOF]; 6];
                for j in 0..DOF {
                    let tr = b[0][j] + b[1][j] + b[2][j];
                    db[0][j] = lambda * tr + 2.0 * mu * b[0][j];
                    db[1][j] = lambda * tr + 2.0 * mu * b[1][j];
                    db[2][j] = lambda * tr + 2.0 * mu * b[2][j];
                    db[3][j] = mu * b[3][j];
                    db[4][j] = mu * b[4][j];
                    db[5][j] = mu * b[5][j];
                }
                for i in 0..DOF {
                    for j in 0..DOF {
                        let mut acc = 0.0;
                        for r in 0..6 {
                            acc += b[r][i] * db[r][j];
                        }
                        k[i * DOF + j] += acc * det_j;
                    }
                }
            }
        }
    }
    k
}

/// `∫ ∇N dV` over one element, twenty-four numbers in the same node-then-component order as
/// [`stiffness`].
///
/// The geometric half of an **eigenstrain** load. A body with a stress-free strain `ε₀` carries a
/// nodal load `∫ Bᵀ D ε₀ dV`, and for an isotropic `ε₀ = e·I` the product `D ε₀` is
/// `(3λ+2μ)·e` on each of the three normal rows and nothing on the shear rows — so the whole
/// element load is `(3λ+2μ)·e` times this vector, and this vector depends on the geometry alone.
///
/// It is exact and needs no quadrature. `∫ ∂N_n/∂x dV = s_x·h²/4`, where `s_x` is `−1` or `+1`
/// depending on which side of the element node `n` sits: a trilinear shape function's derivative
/// integrates to the flux it carries through the face it faces. The eight entries of each
/// component sum to zero, which is the statement that a stress-free strain applies **no net
/// force** — a body expanding into nothing does not push itself across the room.
pub(crate) fn eigen_load(h: f64) -> [f64; DOF] {
    let mut f = [0.0; DOF];
    let quarter = h * h / 4.0;
    for (n, c) in CORNERS.iter().enumerate() {
        for axis in 0..3 {
            let sign = if c[axis] == 0 { -1.0 } else { 1.0 };
            f[3 * n + axis] = sign * quarter;
        }
    }
    f
}

/// Lamé constants from the engineering pair.
///
/// `ν → 0.5` is incompressibility and `λ` diverges; the caller is expected to have refused that
/// already, and [`crate::Elastic::new`] does.
pub(crate) fn lame(youngs: f64, poisson: f64) -> (f64, f64) {
    let mu = youngs / (2.0 * (1.0 + poisson));
    let lambda = youngs * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    (lambda, mu)
}
