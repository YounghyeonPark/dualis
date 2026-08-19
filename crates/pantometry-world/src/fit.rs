//! Choosing a grid for an assembly, by **measuring** each choice rather than predicting it.
//!
//! `ARCHITECTURE.md` names this as the third of the three gaps that block a real assembly, and it
//! is the only one it left open on purpose: "picking a cell that resolves the smallest feature and
//! reporting what was traded is a layer above this one, and would be the first place in the
//! workspace that *guesses* — so it has to guess visibly."
//!
//! It does not have to guess at all. [`Voxels::loss`] already measures what a cell size cost, so a
//! proposal can rasterise at every candidate and report what happened. The only guess left is the
//! *rule* for which row to recommend — which is one sentence on [`Fit::recommended`], and can be
//! disagreed with by reading the table it came from.
//!
//! **Prediction would be wrong here in a specific, documented way.** `Loss::volume_error` is a
//! lattice-point count after the bulges and the cuts have cancelled, and its own doc says so: a
//! sphere at 2.5, 2.0 and 1.5 mm gives `+4.9%, +5.8%, −2.3%` — the first refinement makes it worse
//! and the second changes its sign. A tool that extrapolated from one resolution to the next would
//! be confidently wrong on the commonest shape there is.

use pantometry::shape::{Mesh, Voxels};
use pantometry::units::{Length, LengthVec};

/// What one part costs at one cell size.
#[derive(Debug, Clone, PartialEq)]
pub struct PartAt {
    /// The part's name, as the scene will spell it.
    pub name: String,
    /// Cells this part filled. **Zero is the failure that matters**: a part finer than the grid
    /// does not error, it disappears, and the run is well behaved about a different object.
    pub filled: usize,
    /// `voxel volume / mesh volume − 1`, as a result and not as a bound on the next row.
    pub volume_error: f64,
    /// The share of the voxel volume in cells with a face on the outside — the rasterisation's
    /// *uncertainty*, which is the quantity that actually falls as the grid refines.
    pub boundary_fraction: f64,
    /// Runs one cell wide, where the answer had to be all or nothing.
    pub thin_runs: usize,
    /// Triangles smaller than one cell face: detail the grid cannot hold, counted before anything
    /// is rasterised.
    pub features_below: usize,
    /// Rows the rasteriser could not decide even after retrying. **Any is a red flag**, not a cost.
    pub ambiguous_rows: usize,
}

/// One candidate grid, with what every part cost on it.
#[derive(Debug, Clone, PartialEq)]
pub struct Candidate {
    /// The cell side, in metres.
    pub cell_m: f64,
    /// The grid this implies, from the union of the parts' bounds.
    pub counts: (usize, usize, usize),
    /// `nx · ny · nz`, which is what a run actually costs.
    pub total_cells: usize,
    /// Per part, in the order they were given.
    pub parts: Vec<PartAt>,
}

impl Candidate {
    /// Whether this grid is usable at all: every part is there, and nothing was undecidable.
    ///
    /// Separate from *good*, which is what [`Fit::recommended`] is about. A grid on which a part
    /// vanishes is not a coarse answer, it is an answer about a different assembly.
    pub fn sound(&self) -> bool {
        self.parts
            .iter()
            .all(|p| p.filled > 0 && p.ambiguous_rows == 0)
    }

    /// The worst boundary fraction across the parts — the grid's uncertainty, set by whichever
    /// part it resolves least.
    pub fn worst_boundary(&self) -> f64 {
        self.parts
            .iter()
            .map(|p| p.boundary_fraction)
            .fold(0.0, f64::max)
    }
}

/// A ladder of candidate grids for one assembly, measured.
#[derive(Debug, Clone, PartialEq)]
pub struct Fit {
    /// The corner the grid starts from, which every part is placed against.
    pub origin_m: [f64; 3],
    /// The assembly's extent, in metres — the union of the parts' bounds.
    pub extent_m: [f64; 3],
    /// Coarse first, so the cheapest usable row is the first sound one.
    pub candidates: Vec<Candidate>,
    /// Why the ladder stopped, for a caller that wanted a finer row than it got.
    pub stopped: String,
}

impl Fit {
    /// The row to recommend, and **the rule is a sentence rather than a score**: the coarsest grid
    /// on which every part is present, no row was undecidable, and the worst boundary fraction is
    /// at or under `uncertainty`.
    ///
    /// Boundary fraction rather than volume error, because volume error is erratic by construction
    /// and cannot be steered by — see this module's own note. Coarsest rather than finest, because
    /// the cost of a grid is cubic in its refinement and the caller can always read further down
    /// the table.
    ///
    /// `None` when no row qualifies, which is a real answer rather than a failure: either the parts
    /// are finer than the budget allows, or one of them is not watertight.
    pub fn recommended(&self, uncertainty: f64) -> Option<&Candidate> {
        self.candidates
            .iter()
            .find(|c| c.sound() && c.worst_boundary() <= uncertainty)
    }

    /// The table, as the CLI and the page both show it.
    pub fn render(&self) -> String {
        let mm = |v: f64| v * 1e3;
        let mut out = format!(
            "assembly {:.1} x {:.1} x {:.1} mm from ({:.1}, {:.1}, {:.1}) mm\n",
            mm(self.extent_m[0]),
            mm(self.extent_m[1]),
            mm(self.extent_m[2]),
            mm(self.origin_m[0]),
            mm(self.origin_m[1]),
            mm(self.origin_m[2]),
        );
        out.push_str(
            "\n  cell_mm          grid     cells   filled  volume%  bndry%  thin  fine  ?  part\n",
        );
        for c in &self.candidates {
            for (n, p) in c.parts.iter().enumerate() {
                let head = if n == 0 {
                    format!(
                        "{:7.3}  {:>4}x{:<4}x{:<4} {:>8}",
                        mm(c.cell_m),
                        c.counts.0,
                        c.counts.1,
                        c.counts.2,
                        c.total_cells
                    )
                } else {
                    " ".repeat(32)
                };
                out.push_str(&format!(
                    "{head} {:>8} {:+8.2} {:7.1} {:5} {:5} {:>2}  {}\n",
                    p.filled,
                    p.volume_error * 100.0,
                    p.boundary_fraction * 100.0,
                    p.thin_runs,
                    p.features_below,
                    p.ambiguous_rows,
                    p.name,
                ));
            }
        }
        out.push_str(&format!("\n  {}\n", self.stopped));
        out
    }

    /// The `cells`, `cell_mm` and `parts` of a scene built on `candidate`, ready to paste.
    ///
    /// The origin is not in it, and that is not an omission: `Voxels::onto` places every part
    /// against the grid's own corner, so an assembly whose parts share a CAD origin keeps its
    /// relative placement with no further statement. A scene that needs the parts moved says so
    /// with `poses`.
    pub fn scene_fragment(&self, candidate: &Candidate, material: &str) -> String {
        let parts: Vec<String> = candidate
            .parts
            .iter()
            .map(|p| {
                format!(
                    "      {{ \"stl\": {:?}, \"material\": {material:?} }}",
                    p.name
                )
            })
            .collect();
        format!(
            "    \"cells\": [{}, {}, {}],\n    \"cell_mm\": {:.4},\n    \"parts\": [\n{}\n    ]",
            candidate.counts.0,
            candidate.counts.1,
            candidate.counts.2,
            candidate.cell_m * 1e3,
            parts.join(",\n"),
        )
    }
}

/// Measure a ladder of grids for these parts.
///
/// The ladder comes from the geometry rather than from round numbers: the **thinnest dimension of
/// any part** is what a grid has to resolve, so the candidates give it 1, 2, 4, 8 … cells and stop
/// when the whole grid would exceed `budget_cells`. That makes every row a statement about the
/// assembly — "the thinnest wall gets four cells" — instead of about millimetres, which are a unit
/// and not a feature.
///
/// Fails only when there is nothing to measure: no parts, or a part with no triangles, or an
/// assembly that is flat in every direction.
pub fn propose(parts: &[(String, Mesh)], budget_cells: usize) -> Result<Fit, String> {
    if parts.is_empty() {
        return Err("a grid is for something: give at least one part".to_string());
    }
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    let mut thinnest = f64::INFINITY;
    for (name, mesh) in parts {
        let (a, b) = mesh
            .bounds()
            .ok_or_else(|| format!("{name}: the mesh has no triangles, so it has no extent"))?;
        let (a, b) = (a.to_si(), b.to_si());
        for axis in 0..3 {
            lo[axis] = lo[axis].min(a[axis]);
            hi[axis] = hi[axis].max(b[axis]);
            let side = b[axis] - a[axis];
            if side > 0.0 {
                thinnest = thinnest.min(side);
            }
        }
    }
    let extent = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
    if !thinnest.is_finite() || thinnest <= 0.0 {
        return Err(
            "every part is flat in every direction, so there is nothing to fill".to_string(),
        );
    }

    let origin = LengthVec::from_si(glam::DVec3::new(lo[0], lo[1], lo[2]));
    let mut candidates = Vec::new();
    let mut stopped = String::new();
    for step in 0..12u32 {
        // Cells across the thinnest feature: 1, 2, 4, 8 …
        let across = 1u32 << step;
        let cell = thinnest / across as f64;
        let counts = (
            ((extent[0] / cell).ceil() as usize).max(1),
            ((extent[1] / cell).ceil() as usize).max(1),
            ((extent[2] / cell).ceil() as usize).max(1),
        );
        let total = counts.0.saturating_mul(counts.1).saturating_mul(counts.2);
        if total > budget_cells {
            stopped = format!(
                "stopped at {across} cell(s) across the thinnest feature: the next grid is \
                 {total} cells, past the budget of {budget_cells}"
            );
            break;
        }
        let mut measured = Vec::new();
        for (name, mesh) in parts {
            let voxels = Voxels::onto(mesh, origin, counts, Length::from_si(cell))
                .map_err(|e| format!("{name}: {e}"))?;
            let loss = voxels.loss();
            measured.push(PartAt {
                name: name.clone(),
                filled: voxels.filled(),
                volume_error: loss.volume_error,
                boundary_fraction: loss.boundary_fraction,
                thin_runs: loss.thin_runs,
                features_below: mesh.triangles_below(Length::from_si(cell)),
                ambiguous_rows: loss.ambiguous_rows,
            });
        }
        candidates.push(Candidate {
            cell_m: cell,
            counts,
            total_cells: total,
            parts: measured,
        });
    }
    if stopped.is_empty() {
        stopped = "the ladder ran to twelve rows before the budget did".to_string();
    }
    Ok(Fit {
        origin_m: lo,
        extent_m: extent,
        candidates,
        stopped,
    })
}
