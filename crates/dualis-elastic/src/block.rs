//! A block of elastic material, held and loaded, solved for its displacement.

use dualis_core::conserved::quantity;
use dualis_core::{Domain, Exchange, Kind, Ledger, Reading, Violation};
use dualis_units::{Energy, Force, Length, LengthVec, Pressure, Time};
use glam::DVec3;

use crate::element::{stiffness, CORNERS, DOF};
use crate::Elastic;

/// Conjugate gradients gets this many iterations per degree of freedom before it gives up.
const ITERATION_BUDGET: usize = 3;

/// One of the six faces of the block.
///
/// Named by axis and end rather than by "top" or "left", because which way is up is a property of
/// how somebody drew it and not of the body.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Face {
    /// `x = 0`.
    XLow,
    /// `x = L`.
    XHigh,
    /// `y = 0`.
    YLow,
    /// `y = L`.
    YHigh,
    /// `z = 0`.
    ZLow,
    /// `z = L`.
    ZHigh,
}

impl Face {
    /// Which axis the face's normal is along: 0, 1 or 2.
    pub fn axis(self) -> usize {
        match self {
            Face::XLow | Face::XHigh => 0,
            Face::YLow | Face::YHigh => 1,
            Face::ZLow | Face::ZHigh => 2,
        }
    }

    /// Whether it is the high end of that axis.
    pub fn is_high(self) -> bool {
        matches!(self, Face::XHigh | Face::YHigh | Face::ZHigh)
    }

    /// The outward normal.
    pub fn normal(self) -> DVec3 {
        let mut n = DVec3::ZERO;
        n[self.axis()] = if self.is_high() { 1.0 } else { -1.0 };
        n
    }
}

/// A rectangular block of one material, on a grid of cubic elements.
///
/// Displacements live at the **nodes** — the corners — so a grid of `(nx, ny, nz)` elements has
/// `(nx+1, ny+1, nz+1)` of them. That is where a trilinear element's unknowns are, and putting
/// them at cell centres instead is the first step toward the checkerboard the element module
/// describes.
///
/// # Nothing is solved until it is held
///
/// A free body has six rigid motions that cost no energy, so the system is singular by exactly six
/// dimensions and there is no unique answer. [`Block::solve`] refuses rather than returning one of
/// the infinitely many, and the refusal names the problem — which is more useful than a
/// displacement field that is correct up to a translation nobody asked for.
#[derive(Clone, Debug)]
pub struct Block {
    name: String,
    /// Elements along each axis.
    counts: (usize, usize, usize),
    /// Nodes along each axis, one more than the elements.
    nodes: (usize, usize, usize),
    dx: f64,
    material: Elastic,
    /// The 24×24 element stiffness, once.
    ke: Vec<f64>,
    /// Displacement, three per node.
    u: Vec<f64>,
    /// Whether each degree of freedom is held.
    held: Vec<bool>,
    /// What it is held *at*. Zero unless something prescribed a motion.
    value: Vec<f64>,
    /// Applied nodal force, three per node.
    load: Vec<f64>,
    residual: f64,
    converged: bool,
    saved: Option<Box<Saved>>,
}

/// The state a [`Domain::checkpoint`] puts aside: what the body is doing and what is being done
/// to it.
#[derive(Clone, Debug)]
struct Saved {
    u: Vec<f64>,
    held: Vec<bool>,
    value: Vec<f64>,
    load: Vec<f64>,
}

impl Block {
    /// A block of `counts` cubic elements of side `cell`.
    ///
    /// Nothing is held and nothing is loaded, so it is not solved yet — see the type's docs.
    pub fn new(
        name: impl Into<String>,
        counts: (usize, usize, usize),
        cell: Length,
        material: Elastic,
    ) -> Block {
        let counts = (counts.0.max(1), counts.1.max(1), counts.2.max(1));
        let nodes = (counts.0 + 1, counts.1 + 1, counts.2 + 1);
        let n = nodes.0 * nodes.1 * nodes.2 * 3;
        let dx = cell.to_si();
        let (lambda, mu) = material.lame();
        Block {
            name: name.into(),
            counts,
            nodes,
            dx,
            material,
            ke: stiffness(dx, lambda, mu),
            u: vec![0.0; n],
            held: vec![false; n],
            value: vec![0.0; n],
            load: vec![0.0; n],
            residual: f64::INFINITY,
            converged: false,
            saved: None,
        }
    }

    /// Elements along each axis.
    pub fn elements(&self) -> (usize, usize, usize) {
        self.counts
    }

    /// Nodes along each axis.
    pub fn node_counts(&self) -> (usize, usize, usize) {
        self.nodes
    }

    /// The element side.
    pub fn cell(&self) -> Length {
        Length::from_si(self.dx)
    }

    /// The block's outside dimensions.
    pub fn size(&self) -> LengthVec {
        LengthVec::from_si(
            DVec3::new(
                self.counts.0 as f64,
                self.counts.1 as f64,
                self.counts.2 as f64,
            ) * self.dx,
        )
    }

    /// What it is made of.
    pub fn material(&self) -> Elastic {
        self.material
    }

    fn node(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.nodes.0 * (j + self.nodes.1 * k)
    }

    /// The nodes on a face, as `(i, j, k)`.
    fn face_nodes(&self, face: Face) -> Vec<(usize, usize, usize)> {
        let (nx, ny, nz) = self.nodes;
        let fixed = if face.is_high() {
            [nx - 1, ny - 1, nz - 1][face.axis()]
        } else {
            0
        };
        let mut out = Vec::new();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let along = [i, j, k][face.axis()];
                    if along == fixed {
                        out.push((i, j, k));
                    }
                }
            }
        }
        out
    }

    /// Hold every component of every node on a face. A bond, not a bearing.
    pub fn clamp(&mut self, face: Face) {
        for (i, j, k) in self.face_nodes(face) {
            let n = self.node(i, j, k);
            for a in 0..3 {
                self.held[3 * n + a] = true;
                self.value[3 * n + a] = 0.0;
            }
        }
        self.invalidate();
    }

    /// Hold only the component normal to a face, leaving it free to slide in the other two.
    ///
    /// A roller, and the boundary condition every closed-form case in this crate is stated
    /// against. A symmetry plane is the same thing: material on the other side would stop the
    /// normal motion and permit the rest.
    pub fn roller(&mut self, face: Face) {
        let axis = face.axis();
        for (i, j, k) in self.face_nodes(face) {
            let n = self.node(i, j, k);
            self.held[3 * n + axis] = true;
            self.value[3 * n + axis] = 0.0;
        }
        self.invalidate();
    }

    /// Hold one node completely, to remove a rigid motion a roller set does not.
    pub fn pin(&mut self, i: usize, j: usize, k: usize) {
        let (nx, ny, nz) = self.nodes;
        let n = self.node(i.min(nx - 1), j.min(ny - 1), k.min(nz - 1));
        for a in 0..3 {
            self.held[3 * n + a] = true;
            self.value[3 * n + a] = 0.0;
        }
        self.invalidate();
    }

    /// Hold every node on the outer surface at a displacement of your choosing.
    ///
    /// # This is the patch test
    ///
    /// Prescribe a linear field on the whole boundary and the interior has exactly one answer:
    /// the same linear field. An element that cannot reproduce it is not consistent and will not
    /// converge to the right thing no matter how fine the mesh — which is why this is the test a
    /// finite element is expected to pass before anything else is believed about it.
    ///
    /// It is also the only way to isolate a single modulus. A shear *rig* — clamp the bottom, drag
    /// the top — is not simple shear: the sides are free, the block bends, and what comes out for a
    /// cube is 0.40 of `G` rather than `G`. Prescribing the field removes the rig from the answer.
    ///
    /// The field is given the node's position and returns its displacement.
    pub fn prescribe_boundary(&mut self, field: impl Fn(LengthVec) -> LengthVec) {
        let (nx, ny, nz) = self.nodes;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let surface =
                        i == 0 || j == 0 || k == 0 || i + 1 == nx || j + 1 == ny || k + 1 == nz;
                    if !surface {
                        continue;
                    }
                    let at = LengthVec::from_si(DVec3::new(i as f64, j as f64, k as f64) * self.dx);
                    let u = field(at).to_si();
                    let n = self.node(i, j, k);
                    for a in 0..3 {
                        self.held[3 * n + a] = true;
                        self.value[3 * n + a] = u[a];
                    }
                }
            }
        }
        self.invalidate();
    }

    /// Apply a uniform traction over a face, as a vector in world axes.
    ///
    /// # Consistent nodal loads, not equal shares
    ///
    /// A uniform pressure on a face of `n×m` elements does **not** put the same force on every
    /// node. A trilinear element's shape functions integrate to a quarter of the face area each,
    /// so a node shared by four elements gets four quarters and a corner node gets one — the
    /// interior nodes carry four times what the corners do.
    ///
    /// Splitting the total equally is the obvious mistake and it is nearly invisible: the
    /// resultant is right, equilibrium holds, and the answer is wrong only near the edges, where
    /// it looks like a boundary effect somebody would expect anyway. It also breaks the exact
    /// moduli, which is how this crate would notice.
    pub fn traction(&mut self, face: Face, traction: DVec3) {
        let axis = face.axis();
        let (ea, eb) = match axis {
            0 => (self.counts.1, self.counts.2),
            1 => (self.counts.0, self.counts.2),
            _ => (self.counts.0, self.counts.1),
        };
        let quarter = self.dx * self.dx / 4.0;
        // Walk the face's elements and hand each of its four corner nodes a quarter of its area.
        for b in 0..eb {
            for a in 0..ea {
                for (da, db) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                    let (i, j, k) = match axis {
                        0 => (
                            if face.is_high() { self.counts.0 } else { 0 },
                            a + da,
                            b + db,
                        ),
                        1 => (
                            a + da,
                            if face.is_high() { self.counts.1 } else { 0 },
                            b + db,
                        ),
                        _ => (
                            a + da,
                            b + db,
                            if face.is_high() { self.counts.2 } else { 0 },
                        ),
                    };
                    let n = self.node(i, j, k);
                    for c in 0..3 {
                        self.load[3 * n + c] += traction[c] * quarter;
                    }
                }
            }
        }
        self.invalidate();
    }

    /// Press a face inward with a uniform pressure. Positive presses.
    pub fn press(&mut self, face: Face, pressure: Pressure) {
        self.traction(face, -face.normal() * pressure.to_si());
    }

    /// Pull a face outward along its own normal. Positive pulls.
    pub fn pull(&mut self, face: Face, stress: Pressure) {
        self.traction(face, face.normal() * stress.to_si());
    }

    /// Clear every load, leaving the holds.
    pub fn unload(&mut self) {
        self.load.fill(0.0);
        self.invalidate();
    }

    fn invalidate(&mut self) {
        self.converged = false;
        self.residual = f64::INFINITY;
    }

    /// Displacement at one node.
    pub fn displacement_at(&self, i: usize, j: usize, k: usize) -> LengthVec {
        let (nx, ny, nz) = self.nodes;
        let n = self.node(i.min(nx - 1), j.min(ny - 1), k.min(nz - 1));
        LengthVec::from_si(DVec3::new(
            self.u[3 * n],
            self.u[3 * n + 1],
            self.u[3 * n + 2],
        ))
    }

    /// The mean normal strain along an axis, from the two end faces.
    ///
    /// Measured from the displacement of the faces rather than from a gradient inside, because
    /// that is what a strain gauge on the outside of a part would read and because it is the
    /// quantity every closed form below is written in.
    pub fn mean_strain(&self, axis: usize) -> f64 {
        let axis = axis.min(2);
        let (lo, hi) = match axis {
            0 => (Face::XLow, Face::XHigh),
            1 => (Face::YLow, Face::YHigh),
            _ => (Face::ZLow, Face::ZHigh),
        };
        let mean_of = |face: Face| {
            let nodes = self.face_nodes(face);
            let sum: f64 = nodes
                .iter()
                .map(|(i, j, k)| self.u[3 * self.node(*i, *j, *k) + axis])
                .sum();
            sum / nodes.len().max(1) as f64
        };
        let length = [self.counts.0, self.counts.1, self.counts.2][axis] as f64 * self.dx;
        (mean_of(hi) - mean_of(lo)) / length
    }

    /// The fractional change in volume, to first order: the sum of the three normal strains.
    pub fn volumetric_strain(&self) -> f64 {
        (0..3).map(|a| self.mean_strain(a)).sum()
    }

    /// The net force the holds are carrying on one face, as a vector.
    ///
    /// `K·u − f` summed over the held degrees of freedom of that face's nodes. It is what a load
    /// cell under the part would read, and it is an independent route to a stress: for a uniaxial
    /// test the reaction over the area must equal the applied traction, which is equilibrium
    /// rather than an assumption.
    ///
    /// # Read the component, not the magnitude
    ///
    /// A node on an edge belongs to **two** faces, so a hold placed by one face is summed by both.
    /// With rollers on `x`, `y` and `z` low and a pull along `x`, the `y`-low face reports 2.5 N
    /// out of the 30 N the `x`-low rollers are carrying — its share of the shared edge — even
    /// though it is carrying nothing itself.
    ///
    /// The component along the face's own normal is unambiguous, because only that face's roller
    /// holds it. The magnitude is not, and taking one is the mistake this paragraph exists to stop.
    pub fn reaction(&self, face: Face) -> DVec3 {
        let ku = self.apply(&self.u);
        let mut out = DVec3::ZERO;
        for (i, j, k) in self.face_nodes(face) {
            let n = self.node(i, j, k);
            for a in 0..3 {
                if self.held[3 * n + a] {
                    out[a] += ku[3 * n + a] - self.load[3 * n + a];
                }
            }
        }
        out
    }

    /// The reaction normal to a face — the unambiguous component. See [`Block::reaction`].
    pub fn normal_reaction(&self, face: Face) -> Force {
        Force::from_si(self.reaction(face)[face.axis()])
    }

    /// Strain energy, `½ uᵀKu`.
    pub fn strain_energy(&self) -> Energy {
        let ku = self.apply(&self.u);
        Energy::from_si(0.5 * self.u.iter().zip(&ku).map(|(a, b)| a * b).sum::<f64>())
    }

    /// Work done by the applied loads, `Σ f·u`.
    ///
    /// Independent of [`Block::strain_energy`]: one comes from the stiffness, the other from the
    /// loads. At equilibrium they are in the ratio Clapeyron's theorem gives, and that they are is
    /// a check rather than a restatement.
    pub fn work_done(&self) -> Energy {
        Energy::from_si(
            self.load
                .iter()
                .zip(&self.u)
                .map(|(f, u)| f * u)
                .sum::<f64>(),
        )
    }

    /// How far `2U = Σf·u` is from holding, relative to the work.
    pub fn energy_balance(&self) -> f64 {
        let (u, w) = (self.strain_energy().to_si(), self.work_done().to_si());
        if w.abs() <= 0.0 {
            0.0
        } else {
            (2.0 * u - w).abs() / w.abs()
        }
    }

    /// Whether the last solve met its tolerance.
    pub fn converged(&self) -> bool {
        self.converged
    }

    /// The relative residual the last solve reached.
    pub fn residual(&self) -> f64 {
        self.residual
    }

    /// How many degrees of freedom are free to move.
    pub fn free_dofs(&self) -> usize {
        self.held.iter().filter(|h| !**h).count()
    }

    /// Solve for the displacement.
    ///
    /// Conjugate gradients on the held-out system: a degree of freedom that is held is removed
    /// from the iteration entirely rather than penalised with a large diagonal, so the condition
    /// number is the physics' and not a number somebody picked.
    ///
    /// Returns whether it converged. `false` means the answer is shaped like an answer and is not
    /// one, which is why [`Domain::step`] refuses on it.
    pub fn solve(&mut self, tolerance: f64) -> bool {
        let budget = ITERATION_BUDGET * self.u.len() + 64;
        self.solve_within(tolerance, budget)
    }

    /// Solve, spending at most `max_iterations`.
    pub fn solve_within(&mut self, tolerance: f64, max_iterations: usize) -> bool {
        // **The iterate starts at the prescribed values, not at zero.** Everything below keeps
        // the search direction zero on a held degree of freedom, so whatever `x` holds there at
        // the start it holds at the end — which is how a prescribed motion enters the system
        // without a penalty stiffness or a row elimination.
        let mut x = self.value.clone();
        let b = self.load.clone();

        let ax = self.apply(&x);
        let mut r = vec![0.0; x.len()];
        for (i, ri) in r.iter_mut().enumerate() {
            *ri = if self.held[i] { 0.0 } else { b[i] - ax[i] };
        }
        let mut p = r.clone();
        let mut rr: f64 = r.iter().map(|v| v * v).sum();
        // Measured against whichever is bigger: the loads, or the first residual. Under pure
        // displacement control there are no loads at all, and dividing by their norm would report
        // a relative residual of `1e300` for a solve that is converging perfectly well.
        let scale: f64 = b
            .iter()
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt()
            .max(rr.sqrt())
            .max(f64::MIN_POSITIVE);

        let mut iterations = 0;
        while rr.sqrt() / scale > tolerance && iterations < max_iterations {
            let mut ap = self.apply(&p);
            for (i, held) in self.held.iter().enumerate() {
                if *held {
                    ap[i] = 0.0;
                }
            }
            let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
            if pap <= 0.0 {
                // Not positive definite on the free set, which for this operator means the holds
                // left a rigid motion in. Stopping is right; returning one of the infinitely many
                // answers is not.
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

        self.u = x;
        self.residual = rr.sqrt() / scale;
        self.converged = self.residual <= tolerance;
        self.converged
    }

    /// `K·x`, assembled element by element.
    fn apply(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; x.len()];
        let (ex, ey, ez) = self.counts;
        let mut map = [0usize; 8];
        for e_z in 0..ez {
            for e_y in 0..ey {
                for e_x in 0..ex {
                    for (c, corner) in CORNERS.iter().enumerate() {
                        map[c] = self.node(e_x + corner[0], e_y + corner[1], e_z + corner[2]);
                    }
                    for a in 0..DOF {
                        let row = 3 * map[a / 3] + a % 3;
                        let mut acc = 0.0;
                        for b in 0..DOF {
                            acc += self.ke[a * DOF + b] * x[3 * map[b / 3] + b % 3];
                        }
                        y[row] += acc;
                    }
                }
            }
        }
        y
    }
}

impl Domain for Block {
    fn name(&self) -> &str {
        &self.name
    }

    /// [`Kind::QuasiStatic`]: there is no state to roll forward and no stability limit. A static
    /// equilibrium is the solution of an elliptic problem, not the result of marching one.
    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }

    fn step(&mut self, _t: Time, _dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        if self.solve(1e-10) {
            return Ok(());
        }
        Err(Violation::at(
            self.name.clone(),
            "equilibrium residual",
            self.residual,
        ))
    }

    /// The strain energy the body is holding.
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.strain_energy().to_si())
    }

    fn readings(&self) -> Vec<Reading> {
        vec![
            Reading::new(
                &self.name,
                "strain energy",
                self.strain_energy().to_si(),
                "J",
            ),
            Reading::new(&self.name, "strain x", self.mean_strain(0), ""),
            Reading::new(&self.name, "strain y", self.mean_strain(1), ""),
            Reading::new(&self.name, "strain z", self.mean_strain(2), ""),
            Reading::new(&self.name, "volume change", self.volumetric_strain(), ""),
            Reading::new(&self.name, "residual", self.residual, ""),
        ]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn checkpoint(&mut self) {
        self.saved = Some(Box::new(Saved {
            u: self.u.clone(),
            held: self.held.clone(),
            value: self.value.clone(),
            load: self.load.clone(),
        }));
    }

    fn restore(&mut self) {
        if let Some(s) = self.saved.take() {
            self.u = s.u;
            self.held = s.held;
            self.value = s.value;
            self.load = s.load;
        }
    }
}
