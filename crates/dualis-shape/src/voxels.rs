//! A mesh rasterised into cells, and an account of what the cells could not hold.

use crate::mesh::Mesh;
use dualis_units::{Length, LengthVec, Volume};
use glam::DVec3;

/// What a rasterisation lost, reported rather than left to be discovered.
///
/// A missing feature has no symptom. It does not make a solver fail, produce a `NaN` or trip the
/// conservation audit — it produces a smooth, plausible answer about a different object, which is the
/// failure this workspace is organised around not having. So every [`Voxels`] carries one of these and
/// [`Loss::is_clean`] says in one call whether anything here needs reading.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Loss {
    /// `voxel volume / mesh volume − 1`.
    ///
    /// The aggregate of what the grid kept and did not, and it is **signed**, because the cells the
    /// surface bulges out of and the ones it cuts into partly cancel.
    ///
    /// That cancellation is why this is *not* a discretisation error with an order, and it is worth
    /// saying plainly because the expectation is so natural. Rasterising a sphere of radius 10 mm — a
    /// 32×64 tessellation, and the numbers move a few tenths of a point with that — at 2.5, 2.0 and
    /// 1.5 mm gives `+4.9%, +5.8%, −2.3%`: the first refinement makes it **worse** and the second changes
    /// its sign. Sliding the mesh relative to the grid changes none of it. What is left after the
    /// cancellation is a lattice-point count, and its error is erratic by nature.
    ///
    /// So use it as a *result* and not as a bound on the next resolution.
    pub volume_error: f64,
    /// The share of the voxel volume in cells that have a face on the outside.
    ///
    /// **The rasterisation's uncertainty rather than its error**: the volume that sat close enough to the
    /// surface for the answer to have gone either way.
    ///
    /// Unlike the volume error this is clean, because it is a surface area rather than a cancellation:
    /// the layer is one cell thick over an area `A`, so the fraction is `A·dx/V` — **first order, with a
    /// coefficient**. For a sphere that is `c·3dx/R`, and the tests measure `c ≈ 0.82`.
    ///
    /// It is the number to put in front of someone choosing a cell size. *Forty-three percent of your
    /// object's volume is in cells that could have gone either way* says what 2 mm on a 20 mm ball means
    /// in a way that a step count does not.
    ///
    /// # It does **not** bound the volume error, and an earlier version of this sentence said it did
    ///
    /// The argument is seductive and wrong: "only cells the surface passes through can be
    /// misclassified, and those are these". They are not. This counts exposed cells that came out
    /// **filled**; a cell the surface passes through whose centre landed outside is misclassified too and
    /// appears in neither the numerator nor the denominator. The two are also fractions of different
    /// volumes — this of the voxel volume, the error of the mesh's.
    ///
    /// The 0.4 mm plate in `a_designed_shape.rs` is the counterexample and it was in the suite the whole
    /// time: at a 2 mm cell it reports a volume error of `+4.0` against a boundary fraction of `1.0`.
    pub boundary_fraction: f64,
    /// Solid runs one or two cells thick, counted along all three axes.
    ///
    /// A feature two cells across is not resolved by any scheme here — a seven-point stencil has no
    /// interior in it, a trilinear element has one element — and a feature *one* cell across exists only
    /// by luck of where the surface fell relative to the cell centres. Move the mesh half a cell and it
    /// may not be there at all.
    pub thin_runs: usize,
    /// Triangles smaller than one face of a cell. See [`Mesh::triangles_below`].
    pub small_triangles: usize,
    /// Scanlines whose first ray was degenerate and which a perturbed one decided.
    ///
    /// Not a loss — these rows came out right — but the count says the mesh has geometry lined up with
    /// the grid, which is the common case rather than the exotic one: a cube on cell boundaries sends
    /// **every** row through a face diagonal, and this reads in the hundreds.
    ///
    /// It is here because the alternative is a mechanism nothing can see working. [`Loss::ambiguous_rows`]
    /// counts only the rows where *all* the perturbations failed, and no mesh in the test suite has ever
    /// produced one — so without this field the whole retry path would be exercised only by inference
    /// from a cell count, and a broken retry would look exactly like a mesh that never needed it.
    pub retried_rows: usize,
    /// Scanlines no ray could decide, after every perturbation was tried.
    ///
    /// A ray through a closed surface must cross it an even number of times, and must not pass through
    /// an edge. A row that fails both tests on all four rays cannot be filled by parity. **These rows are
    /// left empty**, which is visible as a slot missing from the shape rather than as a subtly wrong
    /// fill — and the count is here so it is not only visible in a picture.
    ///
    /// **No mesh in the test suite has produced one**, which is stated rather than left to be assumed:
    /// the branch is written and reasoned about and is not covered by a measurement, so a caller who sees
    /// a nonzero value here is in territory this crate has not walked. [`Loss::retried_rows`] is the
    /// nearby thing that *is* exercised.
    pub ambiguous_rows: usize,
}

impl Loss {
    /// Whether anything was lost that a caller should look at.
    ///
    /// The thresholds are deliberately loose and deliberately stated: **2% of volume**, and *any* thin run
    /// or ambiguous row at all. Volume error is a smooth thing that a caller trades against cost, so it
    /// gets a number; a feature one cell thick and a row that could not be filled are not trade-offs, they
    /// are things that either happened or did not.
    ///
    /// Three fields are deliberately **not** here, for two different reasons.
    ///
    /// [`Loss::boundary_fraction`] and [`Loss::small_triangles`] are left out because no threshold on
    /// either is right for more than one physics: a diffusion problem run to steady state barely notices
    /// a boundary layer that would move a stress concentration by a factor of two, and a large flat face
    /// tessellated into a thousand slivers loses nothing at all. Putting a number on those would be this
    /// library guessing, which is the one thing it does not do — so they are reported and the judgement is
    /// the caller's.
    ///
    /// [`Loss::retried_rows`] is left out because it is not a loss. Those rows came out right.
    pub fn is_clean(&self) -> bool {
        self.volume_error.abs() < 0.02 && self.thin_runs == 0 && self.ambiguous_rows == 0
    }
}

/// A mesh rasterised onto a grid of cubes.
///
/// The grid is the mesh's own bounding box, grown to whole cells, with the shape sitting inside it. Cells
/// are **inside or outside** — there is no partial cell, because the domains this feeds have no partial
/// cell either.
#[derive(Clone, Debug)]
pub struct Voxels {
    counts: (usize, usize, usize),
    dx: f64,
    origin: DVec3,
    inside: Vec<bool>,
    loss: Loss,
}

impl Voxels {
    /// Rasterise `mesh` onto cubes of side `cell`.
    ///
    /// # By scanline, and why the crossings are counted a row at a time
    ///
    /// A point-in-mesh test per cell would cast one ray per cell and cost `cells × triangles`. Casting one
    /// ray along `x` for each `(j, k)` row and filling between sorted pairs of crossings costs
    /// `rows × triangles` — the same answer for a factor of `nx` less work, and it is not only faster: a
    /// whole row decided from one sorted list of crossings **cannot** disagree with itself about parity,
    /// where per-cell tests can and do near a surface.
    ///
    /// # Degenerate rays are detected, perturbed deterministically, and then admitted
    ///
    /// Detected, and that word is doing the work. **Parity is not enough to find them.** A ray through
    /// the edge two triangles share hits both, so the count stays even and the two crossings coincide;
    /// the row pairs them against each other, fills nothing, and reports success. A cube loses a whole
    /// diagonal plane of rows to this at every resolution — 64 cells of 512 on an eight-cell cube — with
    /// no error anywhere. So a hit on an edge is a case of its own, distinct from both a miss and a
    /// clean crossing, and a row that produces one is retried whatever its parity says.
    ///
    /// The retry moves the ray, and the offsets are **fixed** — a short list of small multiples of the
    /// cell — so the same mesh and cell size give bit-for-bit the same voxels on every platform and at
    /// every optimisation level, which is a promise the whole workspace makes.
    ///
    /// A row still degenerate or still odd after all of them is left **empty** and counted in
    /// [`Loss::ambiguous_rows`]. Empty is the honest failure: a hole is visible, and a row filled by
    /// guessing which crossing to drop is a wrong shape that looks right.
    ///
    /// Refuses an open mesh, because parity has no meaning through a surface with a hole in it.
    pub fn of(mesh: &Mesh, cell: Length) -> Result<Voxels, String> {
        let dx = cell.to_si();
        if !(dx.is_finite() && dx > 0.0) {
            return Err(format!("a cell size must be finite and positive, is {dx}"));
        }
        if !mesh.is_closed() {
            return Err(
                "the mesh is not closed: some edge is not shared by exactly two triangles, so a ray \
                 can pass through the surface and parity cannot say what is inside. STL stores no \
                 topology, so this is matched on the vertices as written — see `Mesh::is_closed`"
                    .to_string(),
            );
        }
        let (low, high) = mesh
            .bounds()
            .ok_or_else(|| "an empty mesh has no bounds to rasterise".to_string())?;
        let (low, high) = (low.to_si(), high.to_si());

        // Whole cells, with the shape centred in them. The half cell of margin means the surface never
        // lands exactly on the outer boundary, where a crossing is hardest to count.
        let span = high - low;
        let counts = (
            ((span.x / dx).ceil() as usize + 2).max(1),
            ((span.y / dx).ceil() as usize + 2).max(1),
            ((span.z / dx).ceil() as usize + 2).max(1),
        );
        let grid = DVec3::new(
            counts.0 as f64 * dx,
            counts.1 as f64 * dx,
            counts.2 as f64 * dx,
        );
        let origin = low - (grid - span) * 0.5;
        Ok(Voxels::rasterise(mesh, origin, counts, dx))
    }

    /// Rasterise `mesh` onto a grid somebody else chose: a stated origin, cell count and cell.
    ///
    /// **This is what assembly needs.** [`Voxels::of`] gives every mesh its own box and its own
    /// origin, so two parts come back on two grids that share no cell and cannot be adjacent to
    /// each other — `ARCHITECTURE.md` names that as the gap that arrives first in practice.
    /// Rasterised onto one grid, two parts occupy neighbouring cells of the same array, and a
    /// domain that fills from both gets a conducting interface between them for free, because
    /// its stencil already crosses a face between cells of different materials.
    ///
    /// The grid is the caller's and the mesh's coordinates are read as they are written: an STL
    /// carries absolute positions, so where two parts sit relative to each other is what their
    /// files already say. No pose is applied here — `Pose` places a *domain*, and both parts are
    /// inside one domain now.
    ///
    /// **A mesh that does not fit is refused, with both boxes named.** Cropping it silently is
    /// the failure this workspace keeps finding: a part with its corner cut off runs, audits,
    /// renders and answers a question about a different shape.
    pub fn onto(
        mesh: &Mesh,
        origin: LengthVec,
        counts: (usize, usize, usize),
        cell: Length,
    ) -> Result<Voxels, String> {
        let dx = cell.to_si();
        if !(dx.is_finite() && dx > 0.0) {
            return Err(format!("a cell size must be finite and positive, is {dx}"));
        }
        if counts.0 == 0 || counts.1 == 0 || counts.2 == 0 {
            return Err(format!("a grid of {counts:?} cells holds nothing"));
        }
        if !mesh.is_closed() {
            return Err(
                "the mesh is not closed: some edge is not shared by exactly two triangles, so a \
                 ray can pass through the surface and parity cannot say what is inside"
                    .to_string(),
            );
        }
        let (low, high) = mesh
            .bounds()
            .ok_or_else(|| "an empty mesh has no bounds to rasterise".to_string())?;
        let (low, high) = (low.to_si(), high.to_si());
        let o = origin.to_si();
        let far = o + DVec3::new(
            counts.0 as f64 * dx,
            counts.1 as f64 * dx,
            counts.2 as f64 * dx,
        );
        let fits = low.x >= o.x
            && low.y >= o.y
            && low.z >= o.z
            && high.x <= far.x
            && high.y <= far.y
            && high.z <= far.z;
        if !fits {
            return Err(format!(
                "the mesh spans ({:.4}, {:.4}, {:.4}) to ({:.4}, {:.4}, {:.4}) m and the grid \
                 covers ({:.4}, {:.4}, {:.4}) to ({:.4}, {:.4}, {:.4}) m, so part of it would be \
                 cut off — a part with its corner missing runs and audits and answers about a \
                 different shape, so it is refused rather than cropped",
                low.x, low.y, low.z, high.x, high.y, high.z, o.x, o.y, o.z, far.x, far.y, far.z
            ));
        }
        Ok(Voxels::rasterise(mesh, o, counts, dx))
    }

    /// The scanline itself, on whatever grid it is handed. Shared by [`Voxels::of`] and
    /// [`Voxels::onto`] so there is one rasteriser and not two that agree until they do not.
    fn rasterise(mesh: &Mesh, origin: DVec3, counts: (usize, usize, usize), dx: f64) -> Voxels {
        let mut inside = vec![false; counts.0 * counts.1 * counts.2];
        let mut ambiguous_rows = 0;
        let mut retried_rows = 0;
        // Fixed offsets, tried in order. Irregular multiples so a second attempt does not land on the
        // same symmetry the first one did, and constant so the result is reproducible.
        const NUDGE: [(f64, f64); 4] = [(0.0, 0.0), (0.19, 0.07), (-0.11, 0.23), (0.31, -0.29)];

        let mut crossings: Vec<f64> = Vec::new();
        for k in 0..counts.2 {
            for j in 0..counts.1 {
                let mut filled = false;
                for (attempt, (dy, dz)) in NUDGE.into_iter().enumerate() {
                    let y = origin.y + (j as f64 + 0.5 + dy) * dx;
                    let z = origin.z + (k as f64 + 0.5 + dz) * dx;
                    crossings.clear();
                    let mut degenerate = false;
                    for t in mesh.triangles() {
                        match hit_x(t.a, t.b, t.c, y, z) {
                            Hit::Miss => {}
                            Hit::At(x) => crossings.push(x),
                            // One is enough to spoil the row, and the rest of the triangles cannot
                            // un-spoil it.
                            Hit::Degenerate => {
                                degenerate = true;
                                break;
                            }
                        }
                    }
                    // Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87 —
                    // later than the 1.78 this workspace declares and its CI verifies. A declared MSRV
                    // is a promise to a consumer and a lint suggestion is not.
                    #[allow(clippy::manual_is_multiple_of)]
                    let odd = crossings.len() % 2 != 0;
                    if degenerate || odd {
                        continue;
                    }
                    crossings.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
                    for pair in crossings.chunks_exact(2) {
                        for i in 0..counts.0 {
                            let x = origin.x + (i as f64 + 0.5) * dx;
                            if x > pair[0] && x < pair[1] {
                                inside[i + counts.0 * (j + counts.1 * k)] = true;
                            }
                        }
                    }
                    if attempt > 0 {
                        retried_rows += 1;
                    }
                    filled = true;
                    break;
                }
                if !filled {
                    ambiguous_rows += 1;
                }
            }
        }

        let mut voxels = Voxels {
            counts,
            dx,
            origin,
            inside,
            loss: Loss {
                volume_error: 0.0,
                boundary_fraction: 0.0,
                thin_runs: 0,
                small_triangles: mesh.triangles_below(Length::from_si(dx)),
                retried_rows,
                ambiguous_rows,
            },
        };
        let meshed = mesh.volume().to_si();
        voxels.loss.volume_error = if meshed != 0.0 {
            voxels.volume().to_si() / meshed - 1.0
        } else {
            f64::NAN
        };
        voxels.loss.boundary_fraction = voxels.boundary_share();
        voxels.loss.thin_runs = voxels.count_thin_runs();
        voxels
    }

    /// Cells along each axis.
    pub fn counts(&self) -> (usize, usize, usize) {
        self.counts
    }

    /// The cell side.
    pub fn cell(&self) -> Length {
        Length::from_si(self.dx)
    }

    /// The low corner of cell `(0, 0, 0)`.
    pub fn origin(&self) -> LengthVec {
        LengthVec::from_si(self.origin)
    }

    /// Whether a cell is inside the surface. Out-of-range indices are outside.
    ///
    /// This is the predicate a domain's `fill` takes, and the whole coupling between this crate and the
    /// physics:
    ///
    /// ```no_run
    /// # use dualis_shape::Voxels;
    /// # fn go(voxels: &Voxels, block: &mut impl FnMut(&dyn Fn(usize, usize, usize) -> bool)) {
    /// block(&|i, j, k| voxels.contains(i, j, k));
    /// # }
    /// ```
    pub fn contains(&self, i: usize, j: usize, k: usize) -> bool {
        if i >= self.counts.0 || j >= self.counts.1 || k >= self.counts.2 {
            return false;
        }
        self.inside[i + self.counts.0 * (j + self.counts.1 * k)]
    }

    /// How many cells are inside.
    pub fn filled(&self) -> usize {
        self.inside.iter().filter(|b| **b).count()
    }

    /// The volume those cells occupy.
    pub fn volume(&self) -> Volume {
        Volume::from_si(self.filled() as f64 * self.dx.powi(3))
    }

    /// What the rasterisation lost. See [`Loss`].
    pub fn loss(&self) -> Loss {
        self.loss
    }

    /// The share of filled cells with a face on the outside.
    ///
    /// Six neighbours, not twenty-six: a cell touching the outside only at an edge or a corner is not one
    /// the seven-point stencils here exchange anything through, and counting it would make the layer
    /// thicker than the physics sees it.
    fn boundary_share(&self) -> f64 {
        let filled = self.filled();
        if filled == 0 {
            return 0.0;
        }
        let (nx, ny, nz) = self.counts;
        let mut on_surface = 0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if !self.contains(i, j, k) {
                        continue;
                    }
                    let exposed = i == 0
                        || j == 0
                        || k == 0
                        || !self.contains(i - 1, j, k)
                        || !self.contains(i + 1, j, k)
                        || !self.contains(i, j - 1, k)
                        || !self.contains(i, j + 1, k)
                        || !self.contains(i, j, k - 1)
                        || !self.contains(i, j, k + 1);
                    if exposed {
                        on_surface += 1;
                    }
                }
            }
        }
        on_surface as f64 / filled as f64
    }

    /// Solid runs one or two cells long, along all three axes.
    fn count_thin_runs(&self) -> usize {
        let (nx, ny, nz) = self.counts;
        let mut thin = 0;
        let mut tally = |run: usize| {
            if run == 1 || run == 2 {
                thin += 1;
            }
        };
        for k in 0..nz {
            for j in 0..ny {
                let mut run = 0;
                for i in 0..nx {
                    if self.contains(i, j, k) {
                        run += 1;
                    } else {
                        tally(run);
                        run = 0;
                    }
                }
                tally(run);
            }
        }
        for k in 0..nz {
            for i in 0..nx {
                let mut run = 0;
                for j in 0..ny {
                    if self.contains(i, j, k) {
                        run += 1;
                    } else {
                        tally(run);
                        run = 0;
                    }
                }
                tally(run);
            }
        }
        for j in 0..ny {
            for i in 0..nx {
                let mut run = 0;
                for k in 0..nz {
                    if self.contains(i, j, k) {
                        run += 1;
                    } else {
                        tally(run);
                        run = 0;
                    }
                }
                tally(run);
            }
        }
        thin
    }
}

/// What a ray along `+x` did to one triangle.
enum Hit {
    /// It missed, or ran parallel to the triangle's plane. Either way the triangle is not on this row.
    Miss,
    /// It crossed cleanly, at this `x`.
    At(f64),
    /// It went through an edge, or grazed the plane, so this row's parity cannot be trusted.
    ///
    /// **Parity does not detect this, which is why it has its own variant.** A ray through the edge two
    /// triangles share hits *both*, and the count stays even — so a box whose face diagonal passes
    /// through a cell centre yields crossings `[0, 0, 24, 24]`, pairs them as `(0, 0)` and `(24, 24)`,
    /// fills nothing, and reports no error. A cube at any resolution loses a whole diagonal plane of
    /// rows that way, silently. The row has to be retried on a moved ray instead.
    Degenerate,
}

/// Where a ray along `+x` through `(y, z)` meets a triangle.
///
/// Möller–Trumbore reduced to a fixed direction, with the degenerate cases separated out rather than
/// rounded into an answer. Both thresholds are on **scale-free** quantities so a mesh in metres and the
/// same mesh in millimetres are judged alike:
///
/// - `u`, `v` and `1 − u − v` are barycentric, and near zero means the hit is on an edge;
/// - `det / |e1 × e2|` is exactly `−n·x̂`, the triangle normal's component along the ray, so near zero
///   means edge-on. A triangle *at* zero contributes nothing and is a miss; one merely close gives a
///   crossing position divided by that small number, which is a position not worth having.
fn hit_x(a: DVec3, b: DVec3, c: DVec3, y: f64, z: f64) -> Hit {
    // Both are six orders above the rounding and four below anything a real mesh does on purpose, and
    // that gap is the whole justification — neither is a measurement.
    //
    // `u` and `v` come out of two divisions and a handful of products, so a well-conditioned triangle
    // carries a few ulps and a sliver a few thousand; `1e-9` leaves six orders over the worst of that. In
    // the other direction, a *deliberate* feature 1e-9 of a facet across is below any manufacturing
    // tolerance and below `f32`, which is what the file it came from stores. So the band between "this is
    // rounding" and "this is geometry" is wide, and where in it the threshold sits does not matter.
    //
    // What is *not* claimed: the tests do not distinguish `1e-9` from any value between about `1e-16` and
    // `1e-4`, because the case they exercise — a cube's face diagonal — is degenerate exactly rather than
    // nearly. Finding a mesh that lands in between is possible and has not been done.
    const ON_EDGE: f64 = 1e-9;
    const GRAZING: f64 = 1e-9;

    let origin = DVec3::new(0.0, y, z);
    let direction = DVec3::X;
    let (e1, e2) = (b - a, c - a);
    let twice_area = e1.cross(e2).length();
    if twice_area == 0.0 {
        // A degenerate triangle: three collinear vertices, which exporters do produce. It has no
        // inside, so it is on no row.
        return Hit::Miss;
    }
    let h = direction.cross(e2);
    let det = e1.dot(h);
    let along = det / twice_area;
    if along == 0.0 {
        return Hit::Miss;
    }
    if along.abs() < GRAZING {
        return Hit::Degenerate;
    }
    let inv = 1.0 / det;
    let s = origin - a;
    let u = inv * s.dot(h);
    let q = s.cross(e1);
    let v = inv * direction.dot(q);
    let w = 1.0 - u - v;
    if u < -ON_EDGE || v < -ON_EDGE || w < -ON_EDGE {
        return Hit::Miss;
    }
    if u < ON_EDGE || v < ON_EDGE || w < ON_EDGE {
        return Hit::Degenerate;
    }
    Hit::At(inv * e2.dot(q))
}
