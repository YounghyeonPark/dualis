//! Current as a field, so `I²R` is a consequence of a shape rather than a number somebody typed.
//!
//! [`Winding`](crate::Winding) states its resistance from `ρL/A` and dissipates `I²R`. That is
//! exactly right for a wire, which *is* a uniform bar, and it is the model a motor designer
//! actually uses. It is also the whole of what a lumped electrical model can say: given `R`, it
//! returns `I²R`, and the interesting question — where the heat is, and why `R` is what it is —
//! is assumed rather than answered.
//!
//! A [`Conductor`] is the other half. It is a block of material with a conductivity per cell and
//! two electrodes, and it **solves** for the potential:
//!
//! ```text
//!   ∇·(σ ∇φ) = 0        inside
//!   φ = 0, φ = V        on the two electrodes
//!   J·n = 0             everywhere else
//! ```
//!
//! From that, `J = −σ∇φ` and the dissipation is `∫ σ|∇φ|² dV`. Nobody states a resistance; it
//! comes out, and for a uniform bar it comes out as `ρL/A` **exactly** — which is what makes this
//! checkable rather than merely plausible.
//!
//! # What the field formulation buys
//!
//! Three things a lumped resistor cannot say:
//!
//! - **A shape that is not a bar has no `ρL/A`.** A constriction, a via, a busbar with a corner,
//!   a contact patch — each has a resistance that is a property of its geometry, and the closed
//!   forms that exist for them (spreading resistance is `ρ/4a` for a circular contact into a
//!   half-space) are limits rather than formulas you can apply to a shape.
//! - **Where the heat is.** `I²R` gives a total. A current crowding into a corner dissipates in
//!   that corner, and a joint fails at the hot spot rather than at the average. This domain hands
//!   the density over as a field, so a thermal domain on the other side of the bus can take it
//!   *where it landed* rather than as a lump.
//! - **Series and parallel are consequences, not cases.** Two materials in series add their
//!   resistances and two side by side add their conductances, and neither is coded — both fall
//!   out of the same solve. That is what the tests check, because a formulation that got one of
//!   them wrong would still look like electricity.
//!
//! # Quasi-static, and why the solve is not a march
//!
//! This is [`Kind::QuasiStatic`]. Charge relaxes in a metal in about `ε/σ` — 1.5×10⁻¹⁹ s for
//! copper — so on any timescale a simulation cares about the current distribution is the
//! *solution* of an elliptic problem and not the state of a marched one. Nothing here has a
//! stability limit, and the step is a solve rather than an advance.
//!
//! # The solve, and what it refuses to do quietly
//!
//! Conjugate gradients on the symmetric positive-definite system the finite-volume
//! discretisation gives. Deterministic: a fixed iteration order, a fixed starting vector, no
//! threads, no clock. [`Conductor::residual`] reports what it achieved and
//! [`Conductor::converged`] whether it met the tolerance — and a step that did **not** converge
//! returns a [`Violation`] rather than a plausible-looking potential field.
//!
//! That last one is deliberate and is the failure this workspace keeps finding. An iterative
//! solver that stops at its iteration cap returns something shaped exactly like an answer: smooth,
//! bounded, roughly right in the middle and wrong at the edges. Nothing downstream can tell.

use glam::DVec3;
use pantometry_core::conserved::quantity;
use pantometry_core::{Domain, Exchange, Kind, Ledger, Reading, ScalarField, Violation};
use pantometry_units::{
    Conductivity, Current, CurrentDensity, Energy, Length, LengthVec, Power, Resistance,
    Resistivity, Time, Voltage,
};

use crate::HEAT;

/// How hard the solver tries before it gives up, as a multiple of the cell count.
///
/// Conjugate gradients converges in at most `n` iterations in exact arithmetic; in floating point
/// it usually needs far fewer and occasionally a few more. Four times the cell count is generous
/// enough that hitting it means the system is ill-conditioned rather than merely large — an
/// insulating island, a conductivity ratio of 10¹⁵ — and that is worth reporting rather than
/// grinding at.
const ITERATION_BUDGET: usize = 4;

/// A block of conducting material with two electrodes, solved for its potential.
///
/// Cells are cubes. The electrodes are the whole `x = 0` and `x = L` faces, held at fixed
/// potentials; every other surface is insulating, so no current leaves through it. That is the
/// four-terminal arrangement a resistance is *defined* by, and it is what makes `ρL/A` the exact
/// answer for a uniform block rather than an approximation to it.
#[derive(Clone, Debug)]
pub struct Conductor {
    name: String,
    counts: (usize, usize, usize),
    dx: f64,
    /// Conductivity per cell, S/m.
    sigma: Vec<f64>,
    /// Potential per cell, V. The solve's output.
    phi: Vec<f64>,
    drive: f64,
    residual: f64,
    converged: bool,
    dissipated: f64,
    tolerance: f64,
    /// Iterations `step` will spend before it refuses. `None` is the default budget.
    max_iterations: Option<usize>,
}

impl Conductor {
    /// A uniform block of `counts` cubic cells of side `dx`, driven by a potential difference
    /// across its x faces.
    pub fn new(
        name: impl Into<String>,
        counts: (usize, usize, usize),
        dx: Length,
        material: Resistivity,
        drive: Voltage,
    ) -> Conductor {
        let counts = (counts.0.max(1), counts.1.max(1), counts.2.max(1));
        let cells = counts.0 * counts.1 * counts.2;
        let sigma = material.conductivity().to_si();
        let mut built = Conductor {
            name: name.into(),
            counts,
            dx: dx.to_si(),
            sigma: vec![sigma; cells],
            phi: vec![0.0; cells],
            drive: drive.to_si(),
            residual: f64::INFINITY,
            converged: false,
            dissipated: 0.0,
            tolerance: 1e-12,
            max_iterations: None,
        };
        // **Solved at construction**, because a quasi-static domain has no state before its
        // solve. Leaving it unsolved meant the first captured frame reported a resistance
        // computed from a potential field of zeros — 24x below the value `ρL/A` puts a floor
        // under, beside a residual of `inf` that nothing was reading. A number that wrong is
        // easy to spot; the point is that nothing *stopped* it being reported.
        //
        // The mutators below invalidate it. They set `converged` false and `residual` infinite,
        // so a caller who changes the material and reads the answer without re-solving is told.
        built.solve(built.tolerance);
        built
    }

    /// What [`Domain::step`] asks of the solver, and what it refuses below.
    ///
    /// Exposed because a caller with a hard time budget may prefer a bounded solve to an
    /// unbounded one — and because the refusal path needs to be reachable from a test. A domain
    /// whose only failure mode cannot be provoked is a domain whose failure mode is untested.
    pub fn with_solver(mut self, tolerance: f64, max_iterations: usize) -> Conductor {
        self.tolerance = tolerance;
        self.max_iterations = Some(max_iterations);
        self
    }

    /// Cell counts along x, y and z.
    pub fn counts(&self) -> (usize, usize, usize) {
        self.counts
    }

    /// The cell side.
    pub fn spacing(&self) -> Length {
        Length::from_si(self.dx)
    }

    /// The block's extent — cells times spacing.
    pub fn size(&self) -> LengthVec {
        let (nx, ny, nz) = self.counts;
        LengthVec::from_si(DVec3::new(nx as f64, ny as f64, nz as f64) * self.dx)
    }

    /// The cross-section the current passes through.
    pub fn section(&self) -> pantometry_units::Area {
        let (_, ny, nz) = self.counts;
        pantometry_units::Area::from_si((ny * nz) as f64 * self.dx * self.dx)
    }

    /// The potential difference across the electrodes.
    pub fn drive(&self) -> Voltage {
        Voltage::from_si(self.drive)
    }

    /// Give one cell a different material.
    ///
    /// The point of the whole domain: a block that is not one material has no `ρL/A`, and this is
    /// how it stops being one. Out of range is ignored.
    pub fn set_resistivity(&mut self, i: usize, j: usize, k: usize, material: Resistivity) {
        if let Some(idx) = self.index(i, j, k) {
            self.sigma[idx] = material.conductivity().to_si();
            self.converged = false;
            self.residual = f64::INFINITY;
        }
    }

    /// Give a whole slab of cells one material, by a predicate on cell indices.
    ///
    /// The readable way to build a series or parallel arrangement, and the way the closed-form
    /// tests do it.
    pub fn set_region(
        &mut self,
        mut which: impl FnMut(usize, usize, usize) -> bool,
        material: Resistivity,
    ) {
        let (nx, ny, nz) = self.counts;
        let sigma = material.conductivity().to_si();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if which(i, j, k) {
                        self.sigma[i + nx * (j + ny * k)] = sigma;
                    }
                }
            }
        }
        self.converged = false;
        self.residual = f64::INFINITY;
    }

    /// The flat index of a cell, or `None` out of range.
    pub fn index(&self, i: usize, j: usize, k: usize) -> Option<usize> {
        let (nx, ny, nz) = self.counts;
        (i < nx && j < ny && k < nz).then(|| i + nx * (j + ny * k))
    }

    /// The potential at one cell centre.
    pub fn potential_at(&self, i: usize, j: usize, k: usize) -> Voltage {
        let (nx, ny, nz) = self.counts;
        let idx = self
            .index(i.min(nx - 1), j.min(ny - 1), k.min(nz - 1))
            .expect("clamped indices are in range");
        Voltage::from_si(self.phi[idx])
    }

    /// The current density at one cell, by central differences on the potential.
    ///
    /// `J = −σ∇φ`. At a cell against an insulating face the normal component is zero by
    /// construction, because there is no neighbour to differ from.
    pub fn current_density_at(&self, i: usize, j: usize, k: usize) -> DVec3 {
        let (nx, ny, nz) = self.counts;
        let (i, j, k) = (i.min(nx - 1), j.min(ny - 1), k.min(nz - 1));
        let here = self.index(i, j, k).expect("clamped");
        let axis = |lo: Option<usize>, hi: Option<usize>| -> f64 {
            match (lo, hi) {
                (Some(a), Some(b)) => (self.phi[b] - self.phi[a]) / (2.0 * self.dx),
                // Against a face: one-sided, over one spacing rather than two.
                (None, Some(b)) => (self.phi[b] - self.phi[here]) / self.dx,
                (Some(a), None) => (self.phi[here] - self.phi[a]) / self.dx,
                (None, None) => 0.0,
            }
        };
        let grad = DVec3::new(
            axis(
                i.checked_sub(1).and_then(|a| self.index(a, j, k)),
                self.index(i + 1, j, k),
            ),
            axis(
                j.checked_sub(1).and_then(|b| self.index(i, b, k)),
                self.index(i, j + 1, k),
            ),
            axis(
                k.checked_sub(1).and_then(|c| self.index(i, j, c)),
                self.index(i, j, k + 1),
            ),
        );
        -self.sigma[here] * grad
    }

    /// The total current through the block, measured at the driven electrode.
    ///
    /// Measured rather than derived: the current is the sum of what actually crosses the
    /// electrode faces in the solved field. A solve that had not converged would report a current
    /// that disagreed with the one measured at the other electrode, which is exactly the check
    /// [`Conductor::current_balance`] makes.
    pub fn current(&self) -> Current {
        Current::from_si(self.electrode_current(true))
    }

    /// How much the current in disagrees with the current out, relative to the current itself.
    ///
    /// Zero for a converged solve, because charge does not accumulate. This is the number that
    /// says whether the answer is an answer — and it is measured from the two electrodes
    /// independently rather than being a residual the solver reports about itself.
    pub fn current_balance(&self) -> f64 {
        // In at the driven electrode against out at the grounded one, which is the negative of
        // what flows *into* the block there.
        let (a, b) = (self.electrode_current(true), -self.electrode_current(false));
        let scale = a.abs().max(b.abs());
        if scale <= 0.0 {
            0.0
        } else {
            (a - b).abs() / scale
        }
    }

    /// The resistance the geometry has, `V/I`.
    ///
    /// **Not stated anywhere.** For a uniform block this comes out as `ρL/A` to machine
    /// precision; for anything else it comes out as whatever the shape gives, which is the reason
    /// the domain exists.
    pub fn resistance(&self) -> Resistance {
        let i = self.electrode_current(true);
        if i.abs() <= 0.0 {
            return Resistance::from_si(f64::INFINITY);
        }
        Resistance::from_si(self.drive / i)
    }

    /// The power dissipated, `∫σ|∇φ|²dV`, summed over the faces where the gradient actually is.
    ///
    /// Computed from the field rather than as `V·I`, so that the two agreeing is a check and not
    /// a tautology. They agree to machine precision for a converged solve — that is Tellegen's
    /// theorem, and it is the sharpest single statement about whether the discretisation is
    /// self-consistent.
    pub fn dissipation(&self) -> Power {
        let mut total = 0.0;
        for (a, b, g) in self.faces() {
            let dphi = self.phi_of(b) - self.phi_of(a);
            total += g * dphi * dphi;
        }
        Power::from_si(total)
    }

    /// Energy dissipated over the run.
    pub fn dissipated_energy(&self) -> Energy {
        Energy::from_si(self.dissipated)
    }

    /// Whether the last solve met its tolerance.
    pub fn converged(&self) -> bool {
        self.converged
    }

    /// The relative residual the last solve reached.
    pub fn residual(&self) -> f64 {
        self.residual
    }

    /// Solve for the potential, to a relative residual of `tolerance`.
    ///
    /// Conjugate gradients, which is exact for this system in `n` steps in exact arithmetic and
    /// is symmetric positive definite because the conductances are positive and the coupling is
    /// symmetric. Deterministic: fixed order, fixed start, no threads.
    ///
    /// Returns whether it converged. A caller that ignores the answer gets a field that looks
    /// like a field, which is why [`Domain::step`] refuses instead.
    pub fn solve(&mut self, tolerance: f64) -> bool {
        let budget = self
            .max_iterations
            .unwrap_or(ITERATION_BUDGET * self.phi.len() + 32);
        self.solve_within(tolerance, budget)
    }

    /// Solve, spending at most `max_iterations`.
    ///
    /// Returns whether the tolerance was met. **A `false` here is the whole reason the method
    /// returns anything**: the potential field left behind is smooth, bounded and shaped exactly
    /// like an answer, and nothing downstream can tell it from one.
    pub fn solve_within(&mut self, tolerance: f64, max_iterations: usize) -> bool {
        let n = self.phi.len();
        let budget = max_iterations;

        // b holds what the electrodes inject; A is the interior coupling plus the electrode
        // conductances on the diagonal.
        let b = self.source();
        let mut x = std::mem::take(&mut self.phi);
        if x.len() != n {
            x = vec![0.0; n];
        }
        let mut r = b.clone();
        let ax = self.apply(&x);
        for (ri, axi) in r.iter_mut().zip(&ax) {
            *ri -= axi;
        }
        let mut p = r.clone();
        let mut rr: f64 = r.iter().map(|v| v * v).sum();
        let scale: f64 = b
            .iter()
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt()
            .max(f64::MIN_POSITIVE);

        let mut iterations = 0;
        while rr.sqrt() / scale > tolerance && iterations < budget {
            let ap = self.apply(&p);
            let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
            if pap <= 0.0 {
                // Not positive definite, which for this discretisation means every conductance
                // vanished. Stopping is right; pretending is not.
                break;
            }
            let alpha = rr / pap;
            for (xi, pi) in x.iter_mut().zip(&p) {
                *xi += alpha * pi;
            }
            for (ri, api) in r.iter_mut().zip(&ap) {
                *ri -= alpha * api;
            }
            let rr_next: f64 = r.iter().map(|v| v * v).sum();
            let beta = rr_next / rr;
            for (pi, ri) in p.iter_mut().zip(&r) {
                *pi = ri + beta * *pi;
            }
            rr = rr_next;
            iterations += 1;
        }

        self.phi = x;
        self.residual = rr.sqrt() / scale;
        self.converged = self.residual <= tolerance;
        self.converged
    }

    /// Every interior face, as `(cell a, cell b, conductance)`.
    ///
    /// The conductance across a face between two cells of different material is the **harmonic**
    /// mean, not the arithmetic one, because two half-cells in series add their resistances. An
    /// arithmetic mean is the classic mistake here and it is invisible for a uniform block —
    /// which is why the series test uses two materials four orders of magnitude apart.
    fn faces(&self) -> Vec<(Side, Side, f64)> {
        let (nx, ny, nz) = self.counts;
        let area = self.dx * self.dx;
        let mut out = Vec::new();
        // Electrode faces: cell centre to the electrode is half a spacing.
        for k in 0..nz {
            for j in 0..ny {
                let low = i_index(0, j, k, nx, ny);
                out.push((
                    Side::Electrode(false),
                    Side::Cell(low),
                    self.sigma[low] * area / (0.5 * self.dx),
                ));
                let high = i_index(nx - 1, j, k, nx, ny);
                out.push((
                    Side::Cell(high),
                    Side::Electrode(true),
                    self.sigma[high] * area / (0.5 * self.dx),
                ));
            }
        }
        let mut interior = |a: usize, b: usize| {
            let (sa, sb) = (self.sigma[a], self.sigma[b]);
            let g = if sa <= 0.0 || sb <= 0.0 {
                0.0
            } else {
                // Two half-cells in series: 1/G = dx/2/(sa*A) + dx/2/(sb*A).
                area / (0.5 * self.dx / sa + 0.5 * self.dx / sb)
            };
            out.push((Side::Cell(a), Side::Cell(b), g));
        };
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx - 1 {
                    interior(i_index(i, j, k, nx, ny), i_index(i + 1, j, k, nx, ny));
                }
            }
        }
        for k in 0..nz {
            for j in 0..ny - 1 {
                for i in 0..nx {
                    interior(i_index(i, j, k, nx, ny), i_index(i, j + 1, k, nx, ny));
                }
            }
        }
        for k in 0..nz - 1 {
            for j in 0..ny {
                for i in 0..nx {
                    interior(i_index(i, j, k, nx, ny), i_index(i, j, k + 1, nx, ny));
                }
            }
        }
        out
    }

    fn phi_of(&self, s: Side) -> f64 {
        match s {
            Side::Cell(i) => self.phi[i],
            Side::Electrode(high) => {
                if high {
                    self.drive
                } else {
                    0.0
                }
            }
        }
    }

    /// `A·x` for the finite-volume operator: the sum over faces of `G·(x_a − x_b)`.
    fn apply(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; x.len()];
        for (a, b, g) in self.faces() {
            match (a, b) {
                (Side::Cell(i), Side::Cell(j)) => {
                    let d = g * (x[i] - x[j]);
                    y[i] += d;
                    y[j] -= d;
                }
                (Side::Cell(i), Side::Electrode(_)) | (Side::Electrode(_), Side::Cell(i)) => {
                    y[i] += g * x[i];
                }
                (Side::Electrode(_), Side::Electrode(_)) => {}
            }
        }
        y
    }

    /// What the electrodes inject: `G·φ_electrode` for each face touching one.
    fn source(&self) -> Vec<f64> {
        let mut b = vec![0.0; self.phi.len()];
        for (a, c, g) in self.faces() {
            match (a, c) {
                (Side::Electrode(high), Side::Cell(i)) | (Side::Cell(i), Side::Electrode(high)) => {
                    b[i] += g * if high { self.drive } else { 0.0 };
                }
                _ => {}
            }
        }
        b
    }

    /// The current flowing **from** one electrode **into** the block.
    ///
    /// One convention, stated once. Positive at the driven electrode, which sources; negative at
    /// the grounded one, which sinks. The first version of this carried a per-face sign *and* a
    /// negation at the end, and the two cancelled into an answer of exactly the right magnitude
    /// and the wrong sign — which every closed-form test caught immediately, because a resistance
    /// cannot be negative.
    fn electrode_current(&self, high: bool) -> f64 {
        let phi_e = if high { self.drive } else { 0.0 };
        let mut total = 0.0;
        for (a, b, g) in self.faces() {
            let cell = match (a, b) {
                (Side::Electrode(h), Side::Cell(i)) | (Side::Cell(i), Side::Electrode(h))
                    if h == high =>
                {
                    Some(i)
                }
                _ => None,
            };
            if let Some(i) = cell {
                total += g * (phi_e - self.phi[i]);
            }
        }
        total
    }
}

/// One end of a face: an interior cell, or one of the two electrodes.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Side {
    Cell(usize),
    /// `true` for the driven electrode at `x = L`, `false` for the grounded one at `x = 0`.
    Electrode(bool),
}

fn i_index(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
    i + nx * (j + ny * k)
}

impl Domain for Conductor {
    fn books_balance(&self) -> bool {
        true
    }

    fn name(&self) -> &str {
        &self.name
    }

    /// Quasi-static: charge relaxes in `ε/σ`, which for copper is 1.5×10⁻¹⁹ s. On any timescale a
    /// simulation cares about, the current distribution is a solution and not a state.
    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let want = self.tolerance;
        if !self.solve(want) {
            return Err(Violation {
                quantity: "solver residual".to_string(),
                site: format!("{} (conjugate gradients did not converge)", self.name),
                before: want,
                after: self.residual,
                scale: 1.0,
                tolerance: want,
            });
        }
        let joules = self.dissipation().to_si() * dt.to_si();
        self.dissipated += joules;
        bus.publish(HEAT, joules);
        Ok(())
    }

    /// What it has left to give, which for a source held at a fixed voltage is a negative number
    /// that keeps getting more negative.
    ///
    /// The same bookkeeping [`Winding`](crate::Winding) uses: a domain paying joules out has to
    /// say so, or the audit sees energy appear from nowhere. It is not a reserve — an ideal
    /// voltage source has none — so what is reported is the debt.
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, -self.dissipated)
    }

    /// The resistance, the current, the dissipation and how well the solve converged.
    ///
    /// **The residual is a reading**, not merely an internal number, and that is the point of
    /// having it here: an iterative solve that quietly stopped early produces a field shaped like
    /// an answer, and the only thing that would ever say otherwise is a column somebody can look
    /// at.
    fn readings(&self) -> Vec<Reading> {
        vec![
            Reading::new(&self.name, "resistance", self.resistance().to_si(), "ohm"),
            Reading::new(&self.name, "current", self.current().to_si(), "A"),
            Reading::new(&self.name, "dissipating", self.dissipation().to_si(), "W"),
            Reading::new(&self.name, "spent", self.dissipated, "J"),
            Reading::new(&self.name, "residual", self.residual, ""),
        ]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }

    /// The **potential**, in volts, as a field.
    ///
    /// Not the current density, which is a vector and has no `ScalarField`. A view wanting to see
    /// where the current crowds should take `|J|` from
    /// [`current_density_at`](Conductor::current_density_at); the potential is what the solve
    /// produces and what a contour plot of an electrical problem conventionally shows.
    fn as_field(&self) -> Option<&dyn ScalarField> {
        Some(self)
    }
}

impl ScalarField for Conductor {
    fn unit(&self) -> &'static str {
        "V"
    }

    /// Trilinear between cell centres, clamped at the faces.
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let (nx, ny, nz) = self.counts;
        let q = p.to_si() / self.dx - DVec3::splat(0.5);
        if q.is_nan() {
            return self.phi[0];
        }
        let axis = |v: f64, n: usize| -> (usize, f64) {
            let last = n.saturating_sub(1);
            if v <= 0.0 {
                (0, 0.0)
            } else if v >= last as f64 {
                (last, 0.0)
            } else {
                let i = v.floor();
                (i as usize, v - i)
            }
        };
        let (i, fx) = axis(q.x, nx);
        let (j, fy) = axis(q.y, ny);
        let (k, fz) = axis(q.z, nz);
        let (i1, j1, k1) = (
            (i + 1).min(nx - 1),
            (j + 1).min(ny - 1),
            (k + 1).min(nz - 1),
        );
        let g = |a: usize, b: usize, c: usize| self.phi[i_index(a, b, c, nx, ny)];
        let lerp = |lo: f64, hi: f64, t: f64| lo * (1.0 - t) + hi * t;
        let z0 = lerp(
            lerp(g(i, j, k), g(i1, j, k), fx),
            lerp(g(i, j1, k), g(i1, j1, k), fx),
            fy,
        );
        let z1 = lerp(
            lerp(g(i, j, k1), g(i1, j, k1), fx),
            lerp(g(i, j1, k1), g(i1, j1, k1), fx),
            fy,
        );
        lerp(z0, z1, fz)
    }

    /// `∇φ`, whose negative times `σ` is the current density.
    fn gradient(&self, p: LengthVec, t: Time, h: Length) -> DVec3 {
        let d = h.to_si().max(self.dx);
        let sample = |o: DVec3| self.at(LengthVec::from_si(p.to_si() + o), t);
        DVec3::new(
            (sample(DVec3::X * d) - sample(-DVec3::X * d)) / (2.0 * d),
            (sample(DVec3::Y * d) - sample(-DVec3::Y * d)) / (2.0 * d),
            (sample(DVec3::Z * d) - sample(-DVec3::Z * d)) / (2.0 * d),
        )
    }

    /// Zero. A quasi-static potential does not evolve; it is re-solved when something changes.
    fn rate(&self, _p: LengthVec, _t: Time, _dt: Time) -> f64 {
        0.0
    }
}

/// The current density as a vector, for a caller that wants `J` rather than `φ`.
impl Conductor {
    /// `|J|` at a cell, which is what a picture of current crowding wants.
    pub fn current_density_magnitude(&self, i: usize, j: usize, k: usize) -> CurrentDensity {
        CurrentDensity::from_si(self.current_density_at(i, j, k).length())
    }

    /// The conductivity of one cell.
    pub fn conductivity_at(&self, i: usize, j: usize, k: usize) -> Conductivity {
        let (nx, ny, nz) = self.counts;
        let idx = self
            .index(i.min(nx - 1), j.min(ny - 1), k.min(nz - 1))
            .expect("clamped");
        Conductivity::from_si(self.sigma[idx])
    }
}
