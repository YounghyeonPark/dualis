//! Barnes-Hut gravity: `O(n log n)`, genuinely parallel, and bit-reproducible.
//!
//! [`NBody`](crate::NBody) sums every pair. That is exact, and it conserves momentum
//! to the last bit because each pair's force is computed once and applied with
//! opposite signs — but it costs `n²` and the trick that makes it exact is precisely
//! what makes it awkward to parallelise, since two threads handling different pairs
//! both need to write to the same body.
//!
//! A tree gives that up deliberately. Distant groups of bodies are replaced by their
//! centre of mass, so each body's force is computed **independently** from the tree,
//! reading shared state and writing only its own slot.
//!
//! # Which is why this is the one that parallelises safely
//!
//! There is no reduction. Each thread owns a disjoint range of the output, so there is
//! nothing to merge and no floating-point summation whose order could vary. The
//! headline claim of this workspace — that determinism survives parallelism — is not
//! argued here, it is *executed*: [`TreeNBody::with_threads`] changes how many threads
//! run and
//! `tree::tests::parallel_and_sequential_agree_bit_for_bit` asserts the results are
//! identical to the bit across one, two, four and eight of them.
//!
//! Compare with what the exact solver would need. Splitting `i < j` pairs across
//! threads means several threads accumulating into the same body, and floating-point
//! addition is not associative, so the answer would depend on which thread finished
//! first. Making that deterministic means either locking or a fixed-order merge, and
//! both cost more than they save.
//!
//! # What it costs: momentum stops being exact
//!
//! Body `i` sees the tree's approximation of the rest, and body `j` sees a *different*
//! approximation of the rest, so the two forces between them are no longer equal and
//! opposite. Their sum is no longer zero, and total momentum drifts.
//!
//! That is not a bug to be fixed; it is the price of the approximation, and it shrinks
//! with the opening angle. [`TreeNBody::with_theta`] controls it: at `θ = 0` the tree
//! degenerates to direct summation and the drift vanishes, at `θ = 0.5` it is a few
//! parts in `10⁵` per step, and at `θ = 1.0` it is worse and the force itself is only
//! good to a percent.
//!
//! So a simulation using this cannot audit momentum at
//! [`NBody`](crate::NBody)'s tolerance, and the right response is to know why rather
//! than to loosen the number until it passes.
//!
//! # Threads on WebAssembly
//!
//! `wasm32` targets have no threads to spawn, so the parallel path compiles to the
//! sequential one there. The bit-identical test therefore passes trivially on wasm
//! rather than proving anything — it is the native runs in CI that establish the
//! claim, and the wasm run that establishes the *answer* is the same one.

use dualis_core::{velocity_verlet, Domain, Exchange, Kind, Ledger, Newtonian, Violation};
use dualis_units::{Energy, Length, Mass, Time};
use glam::DVec3;

use crate::{conserved, Body, Coords, GRAVITATION};

/// One cell of the octree.
#[derive(Clone, Debug)]
struct Cell {
    /// Geometric centre of the cube this cell covers.
    centre: DVec3,
    /// Half the cube's side.
    half: f64,
    /// Total mass below here.
    mass: f64,
    /// Mass-weighted position sum; divided through once the tree is built.
    moment: DVec3,
    /// Second moments `Σ m rᵢrⱼ` about the origin, upper triangle in the order
    /// `xx, xy, xz, yy, yz, zz`.
    ///
    /// Accumulated about the origin rather than about the centre of mass, because the
    /// centre of mass is not known until every body has been inserted. Shifting to the
    /// centre of mass is one subtraction, done once when the tree is finished.
    second: [f64; 6],
    /// The one body in this cell, if it is a leaf holding exactly one.
    body: Option<usize>,
    /// Indices of the eight children, or none for a leaf.
    children: Option<[usize; 8]>,
}

impl Cell {
    fn new(centre: DVec3, half: f64) -> Cell {
        Cell {
            centre,
            half,
            mass: 0.0,
            moment: DVec3::ZERO,
            second: [0.0; 6],
            body: None,
            children: None,
        }
    }

    /// Centre of mass, once the tree is complete.
    fn com(&self) -> DVec3 {
        if self.mass > 0.0 {
            self.moment / self.mass
        } else {
            self.centre
        }
    }

    /// The traceless quadrupole tensor about this cell's centre of mass,
    /// `Qᵢⱼ = Σ m (3 rᵢrⱼ - r² δᵢⱼ)`, as the upper triangle.
    ///
    /// Traceless by construction, which is what makes it the *correction* to the
    /// monopole rather than a second, redundant description of the same mass.
    fn quadrupole(&self) -> [f64; 6] {
        if self.mass <= 0.0 {
            return [0.0; 6];
        }
        // Shift the second moments from the origin to the centre of mass:
        // Σ m (r - c)i (r - c)j = Σ m ri rj - M ci cj.
        let c = self.com();
        let m = self.mass;
        let s = [
            self.second[0] - m * c.x * c.x,
            self.second[1] - m * c.x * c.y,
            self.second[2] - m * c.x * c.z,
            self.second[3] - m * c.y * c.y,
            self.second[4] - m * c.y * c.z,
            self.second[5] - m * c.z * c.z,
        ];
        let trace = s[0] + s[3] + s[5];
        [
            3.0 * s[0] - trace,
            3.0 * s[1],
            3.0 * s[2],
            3.0 * s[3] - trace,
            3.0 * s[4],
            3.0 * s[5] - trace,
        ]
    }
}

/// `Q · d`, and `dᵀ Q d`, from the packed upper triangle.
fn contract(q: &[f64; 6], d: DVec3) -> (DVec3, f64) {
    let qd = DVec3::new(
        q[0] * d.x + q[1] * d.y + q[2] * d.z,
        q[1] * d.x + q[3] * d.y + q[4] * d.z,
        q[2] * d.x + q[4] * d.y + q[5] * d.z,
    );
    (qd, qd.dot(d))
}

/// An octree over a set of point masses.
///
/// Built by inserting bodies in index order, so the structure is a function of the
/// positions alone and not of anything the machine was doing at the time.
#[derive(Clone, Debug)]
struct Octree {
    cells: Vec<Cell>,
    softening: f64,
    /// Whether to add the quadrupole correction when a cell is accepted.
    quadrupole: bool,
}

/// Deepest subdivision before coincident bodies are simply pooled.
///
/// Two bodies at the same position would otherwise subdivide until the cube
/// underflowed. Sixty-four levels is far past any real configuration, and pooling
/// rather than diverging is the honest failure: bodies that close together are not
/// resolved by this method anyway, which is what softening is for.
const MAX_DEPTH: u32 = 64;

impl Octree {
    fn build(positions: &[DVec3], masses: &[f64], softening: f64, quadrupole: bool) -> Octree {
        let mut tree = Octree {
            cells: Vec::with_capacity(positions.len() * 2),
            softening,
            quadrupole,
        };
        if positions.is_empty() {
            tree.cells.push(Cell::new(DVec3::ZERO, 1.0));
            return tree;
        }

        // A cube that holds everything. Padded so a body cannot sit exactly on a face,
        // where the octant it belongs to would depend on a rounding.
        let mut lo = positions[0];
        let mut hi = positions[0];
        for p in positions.iter() {
            lo = lo.min(*p);
            hi = hi.max(*p);
        }
        let centre = (lo + hi) * 0.5;
        let extent = (hi - lo).max_element();
        let half = if extent > 0.0 {
            extent * 0.5 * 1.01
        } else {
            1.0
        };

        tree.cells.push(Cell::new(centre, half));
        for i in 0..positions.len() {
            tree.insert(0, i, positions, masses, 0);
        }
        tree
    }

    /// Which child a point belongs in: one bit per axis, so the order is fixed.
    fn octant(centre: DVec3, p: DVec3) -> usize {
        usize::from(p.x >= centre.x)
            | (usize::from(p.y >= centre.y) << 1)
            | (usize::from(p.z >= centre.z) << 2)
    }

    fn subdivide(&mut self, cell: usize) -> [usize; 8] {
        let (centre, half) = (self.cells[cell].centre, self.cells[cell].half);
        let quarter = half * 0.5;
        let mut children = [0usize; 8];
        for (octant, slot) in children.iter_mut().enumerate() {
            let offset = DVec3::new(
                if octant & 1 != 0 { quarter } else { -quarter },
                if octant & 2 != 0 { quarter } else { -quarter },
                if octant & 4 != 0 { quarter } else { -quarter },
            );
            *slot = self.cells.len();
            self.cells.push(Cell::new(centre + offset, quarter));
        }
        self.cells[cell].children = Some(children);
        children
    }

    fn insert(
        &mut self,
        cell: usize,
        body: usize,
        positions: &[DVec3],
        masses: &[f64],
        depth: u32,
    ) {
        let (m, p) = (masses[body], positions[body]);
        self.cells[cell].mass += m;
        self.cells[cell].moment += p * m;
        let s = &mut self.cells[cell].second;
        s[0] += m * p.x * p.x;
        s[1] += m * p.x * p.y;
        s[2] += m * p.x * p.z;
        s[3] += m * p.y * p.y;
        s[4] += m * p.y * p.z;
        s[5] += m * p.z * p.z;

        if let Some(children) = self.cells[cell].children {
            let octant = Self::octant(self.cells[cell].centre, positions[body]);
            self.insert(children[octant], body, positions, masses, depth + 1);
            return;
        }

        match self.cells[cell].body {
            None => {
                self.cells[cell].body = Some(body);
            }
            Some(existing) => {
                if depth >= MAX_DEPTH {
                    // Coincident, or close enough that subdividing would not separate
                    // them. Leave both pooled in this cell's monopole; `body` keeps the
                    // first so the cell still has a direct term.
                    return;
                }
                self.cells[cell].body = None;
                let children = self.subdivide(cell);
                for moved in [existing, body] {
                    let octant = Self::octant(self.cells[cell].centre, positions[moved]);
                    self.insert(children[octant], moved, positions, masses, depth + 1);
                }
            }
        }
    }

    /// Gravitational acceleration on body `target`, by traversal.
    ///
    /// The traversal order is fixed — an explicit stack, children pushed in octant
    /// order — so the sum accumulates the same way every time. That is what makes the
    /// answer independent of which thread computed it.
    fn acceleration_on(
        &self,
        target: usize,
        positions: &[DVec3],
        masses: &[f64],
        theta: f64,
    ) -> DVec3 {
        let g = GRAVITATION.to_si();
        let soft2 = self.softening * self.softening;
        let p = positions[target];
        let mut acc = DVec3::ZERO;
        let mut stack = vec![0usize];

        while let Some(index) = stack.pop() {
            let cell = &self.cells[index];
            if cell.mass <= 0.0 {
                continue;
            }
            let com = cell.com();
            let d = com - p;
            let r2 = d.length_squared() + soft2;

            let is_leaf = cell.children.is_none();
            // Barnes-Hut's opening criterion: a cell may be treated as one mass when
            // its width subtends less than `theta` from here. Squared, to keep the
            // square root out of the inner loop.
            let opens = if is_leaf {
                false
            } else {
                let width = cell.half * 2.0;
                width * width > theta * theta * r2
            };

            if opens {
                if let Some(children) = cell.children {
                    for child in children.iter().rev() {
                        stack.push(*child);
                    }
                }
                continue;
            }

            if is_leaf {
                // A leaf's own body must not pull on itself.
                if cell.body == Some(target) && (cell.mass - masses[target]).abs() < 1e-30 {
                    continue;
                }
            }
            if r2 <= soft2 && soft2 == 0.0 {
                continue;
            }
            let inv_r = r2.sqrt().recip();
            let inv_r3 = inv_r * inv_r * inv_r;
            acc += d * (g * cell.mass * inv_r3);

            // The next term of the multipole expansion. A cell is not really a point,
            // and the quadrupole is the leading correction for its being spread out —
            // one order in `s/d` better than the monopole, so the same accuracy comes at
            // a wider opening angle and fewer cells opened.
            if self.quadrupole && !is_leaf {
                let q = cell.quadrupole();
                if q.iter().any(|v| *v != 0.0) {
                    let (qd, dqd) = contract(&q, d);
                    let inv_r5 = inv_r3 * inv_r * inv_r;
                    let inv_r7 = inv_r5 * inv_r * inv_r;
                    // From the gradient of `Q_ij d_i d_j / (2 r^5)` with `d = c - p`:
                    //
                    //   a = G ( -Q·d / r^5  +  5/2 (dᵀQd) d / r^7 )
                    //
                    // Worth writing out, because the sign is not guessable and getting
                    // it backwards makes the correction *increase* the error — which is
                    // exactly what the accuracy test caught the first time. Sanity: for
                    // mass elongated along the line of sight, `Q` along that line is
                    // positive and this pulls harder than the monopole, as it should,
                    // since the near end of the distribution is closer than its centre.
                    acc += (d * (2.5 * dqd * inv_r7) - qd * inv_r5) * g;
                }
            }
        }
        acc
    }
}

/// Point masses under gravity, evaluated through an octree.
///
/// The approximate, parallel counterpart to [`NBody`](crate::NBody). Read the module
/// documentation before using it for anything whose momentum matters.
pub struct TreeNBody {
    name: String,
    masses: Vec<f64>,
    positions: Coords,
    velocities: Coords,
    softening: f64,
    theta: f64,
    quadrupole: bool,
    threads: usize,
    saved: Option<(Coords, Coords)>,
}

impl TreeNBody {
    /// A Barnes-Hut tree over these bodies: `O(n log n)`, parallel, and approximate.
    ///
    /// Approximate in a way worth knowing before using it — each body sees its own truncated
    /// expansion of the rest, so their mutual forces no longer cancel and momentum drifts.
    /// The drift is set by [`TreeNBody::with_theta`] and vanishes at zero.
    pub fn new(name: impl Into<String>, bodies: &[Body]) -> TreeNBody {
        TreeNBody {
            name: name.into(),
            masses: bodies.iter().map(|b| b.mass.to_si()).collect(),
            positions: Coords(bodies.iter().map(|b| b.position.to_si()).collect()),
            velocities: Coords(bodies.iter().map(|b| b.velocity.to_si()).collect()),
            softening: 0.0,
            theta: 0.5,
            quadrupole: false,
            threads: 1,
            saved: None,
        }
    }

    /// Include the quadrupole correction when a cell is accepted as a single mass.
    ///
    /// A cell is not a point, and the quadrupole is the leading correction for its
    /// being spread out. It buys one more order in `s/d`, so the same force accuracy
    /// comes at a wider opening angle — and since the cost of the tree is dominated by
    /// how many cells get opened, a wider angle is a real saving that more than pays
    /// for the six extra numbers per cell.
    ///
    /// It does not repair the momentum drift. That comes from body `i` and body `j`
    /// seeing different approximations of each other, and a better approximation is
    /// still a different one.
    pub fn with_quadrupole(mut self, quadrupole: bool) -> TreeNBody {
        self.quadrupole = quadrupole;
        self
    }

    /// Opening angle. Smaller is more accurate and slower; `0` is exact and `O(n²)`.
    ///
    /// 0.5 is the conventional working value. Above about 1 the monopole approximation
    /// is being asked to stand in for groups that are not far enough away for it, and
    /// the force error goes past a percent.
    pub fn with_theta(mut self, theta: f64) -> TreeNBody {
        self.theta = theta.max(0.0);
        self
    }

    /// Soften the singularity, as [`NBody::with_softening`](crate::NBody::with_softening).
    pub fn with_softening(mut self, softening: Length) -> TreeNBody {
        self.softening = softening.to_si().abs();
        self
    }

    /// How many threads to evaluate forces on. 1 is sequential.
    ///
    /// The result does not depend on this number. That is the entire point, and it is
    /// asserted rather than asserted-about: see the module docs.
    pub fn with_threads(mut self, threads: usize) -> TreeNBody {
        self.threads = threads.max(1);
        self
    }

    /// How many bodies.
    pub fn count(&self) -> usize {
        self.masses.len()
    }

    /// The opening angle in use. Smaller is more accurate and slower; zero is exact and
    /// `O(n²)`.
    pub fn theta(&self) -> f64 {
        self.theta
    }

    /// One body.
    pub fn body(&self, index: usize) -> Body {
        Body {
            mass: Mass::from_si(self.masses[index]),
            position: dualis_units::LengthVec::from_si(self.positions.0[index]),
            velocity: dualis_units::VelocityVec::from_si(self.velocities.0[index]),
        }
    }

    /// Total linear momentum. Not conserved here, unlike in `NBody` — see
    /// [`TreeNBody::new`].
    pub fn momentum(&self) -> dualis_units::MomentumVec {
        let p = self
            .masses
            .iter()
            .zip(self.velocities.0.iter())
            .map(|(m, v)| *v * *m)
            .fold(DVec3::ZERO, |a, b| a + b);
        dualis_units::MomentumVec::from_si(p)
    }

    /// Summed `½mv²` over every body.
    pub fn kinetic_energy(&self) -> Energy {
        let k: f64 = self
            .masses
            .iter()
            .zip(self.velocities.0.iter())
            .map(|(m, v)| 0.5 * m * v.length_squared())
            .sum();
        Energy::from_si(k)
    }

    /// Accelerations of every body, evaluated through a freshly built tree.
    ///
    /// Public because it is the thing worth comparing against
    /// [`NBody`](crate::NBody)'s exact answer, and the thing whose independence from
    /// the thread count is the claim.
    pub fn accelerations(&self, positions: &[DVec3]) -> Vec<DVec3> {
        let tree = Octree::build(positions, &self.masses, self.softening, self.quadrupole);
        let n = positions.len();
        let mut out = vec![DVec3::ZERO; n];
        let threads = self.threads.min(n.max(1));

        if threads <= 1 {
            for (i, a) in out.iter_mut().enumerate() {
                *a = tree.acceleration_on(i, positions, &self.masses, self.theta);
            }
            return out;
        }

        #[cfg(not(target_family = "wasm"))]
        {
            let chunk = n.div_ceil(threads);
            let masses = &self.masses;
            let theta = self.theta;
            // Each thread owns a disjoint slice of the output and writes nothing else,
            // so there is no reduction and therefore no order to depend on.
            std::thread::scope(|scope| {
                for (c, slice) in out.chunks_mut(chunk).enumerate() {
                    let tree = &tree;
                    let base = c * chunk;
                    scope.spawn(move || {
                        for (k, a) in slice.iter_mut().enumerate() {
                            *a = tree.acceleration_on(base + k, positions, masses, theta);
                        }
                    });
                }
            });
            out
        }

        // No threads to spawn here, so the answer comes from the same code path the
        // sequential branch uses. Identical results, for a less interesting reason.
        #[cfg(target_family = "wasm")]
        {
            for (i, a) in out.iter_mut().enumerate() {
                *a = tree.acceleration_on(i, positions, &self.masses, self.theta);
            }
            out
        }
    }
}

impl Newtonian for TreeNBody {
    type Coords = Coords;

    fn acceleration(&self, x: &Coords, _t: Time) -> Coords {
        Coords(self.accelerations(&x.0))
    }
}

impl Domain for TreeNBody {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    fn max_stable_dt(&self, _now: Time) -> Time {
        let g = GRAVITATION.to_si();
        let mut fastest = 0.0f64;
        for i in 0..self.masses.len() {
            for j in (i + 1)..self.masses.len() {
                let r = (self.positions.0[j] - self.positions.0[i]).length();
                let r = (r * r + self.softening * self.softening).sqrt();
                if r <= 0.0 {
                    continue;
                }
                let m = self.masses[i] + self.masses[j];
                fastest = fastest.max((g * m / (r * r * r)).sqrt());
            }
        }
        if fastest <= 0.0 {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(1.0 / (20.0 * fastest))
    }

    fn step(&mut self, t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let mut x = std::mem::replace(&mut self.positions, Coords(Vec::new()));
        let mut v = std::mem::replace(&mut self.velocities, Coords(Vec::new()));
        velocity_verlet(self, &mut x, &mut v, t, dt);
        self.positions = x;
        self.velocities = v;
        Ok(())
    }

    /// Momentum per axis, exactly as [`NBody`](crate::NBody) reports it — but it will
    /// not hold to the same tolerance, and the module docs say why.
    /// Momentum per axis, accumulated **one body at a time**.
    ///
    /// The per-body `add` is the whole point and not a style choice. `Ledger` records the
    /// largest single contribution as an entry's `scale`, and [`audit`](dualis_core::audit)
    /// judges a change against that rather than against the total — because a correct
    /// system's total momentum is usually exactly zero, and a relative tolerance on zero
    /// means nothing.
    ///
    /// Handing the pre-summed total to one `with` call sets the scale to `|total|`, which for
    /// a symmetric system is `0.0`. `audit` then skips the entry entirely at its
    /// `scale < 1e-300` guard, and the momentum audit silently does not run. That is what
    /// this did. See `NBody::ledger` for how it was found.
    fn ledger(&self) -> Ledger {
        let mut ledger = Ledger::new();
        for (m, v) in self.masses.iter().zip(self.velocities.0.iter()) {
            let p = *v * *m;
            ledger.add(conserved::MOMENTUM_X, p.x);
            ledger.add(conserved::MOMENTUM_Y, p.y);
            ledger.add(conserved::MOMENTUM_Z, p.z);
        }
        ledger
    }

    fn checkpoint(&mut self) {
        self.saved = Some((self.positions.clone(), self.velocities.clone()));
    }

    fn restore(&mut self) {
        if let Some((x, v)) = self.saved.clone() {
            self.positions = x;
            self.velocities = v;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NBody;
    use dualis_core::Rng;
    use dualis_units::{LengthVec, VelocityVec};

    /// A reproducible cloud of bodies, from the kernel's deterministic generator.
    fn cloud(count: usize, seed: u64) -> Vec<Body> {
        (0..count)
            .map(|i| {
                let mut rng = Rng::for_index(seed, i as u64);
                Body::new(
                    Mass::kg(1e9 + rng.range(0.0, 9e9)),
                    LengthVec::m(
                        rng.range(-500.0, 500.0),
                        rng.range(-500.0, 500.0),
                        rng.range(-500.0, 500.0),
                    ),
                    VelocityVec::m_per_s(
                        rng.range(-0.1, 0.1),
                        rng.range(-0.1, 0.1),
                        rng.range(-0.1, 0.1),
                    ),
                )
            })
            .collect()
    }

    /// **The test this module exists for.** One thread, two, four, eight — the same
    /// numbers, to the bit.
    ///
    /// The workspace has claimed since its first commit that determinism survives
    /// parallelism, and until now no code in it actually ran in parallel. This does,
    /// and the guarantee holds because each thread writes only its own slice of the
    /// output: there is no reduction, so there is no summation order to vary.
    #[test]
    fn parallel_and_sequential_agree_bit_for_bit() {
        let bodies = cloud(400, 0xBADC0FFEE);
        let positions: Vec<DVec3> = bodies.iter().map(|b| b.position.to_si()).collect();

        let reference = TreeNBody::new("t", &bodies)
            .with_theta(0.5)
            .with_threads(1)
            .accelerations(&positions);

        for threads in [2usize, 3, 4, 8, 16] {
            let parallel = TreeNBody::new("t", &bodies)
                .with_theta(0.5)
                .with_threads(threads)
                .accelerations(&positions);
            assert_eq!(parallel.len(), reference.len());
            for (i, (a, b)) in parallel.iter().zip(reference.iter()).enumerate() {
                assert_eq!(
                    (a.x.to_bits(), a.y.to_bits(), a.z.to_bits()),
                    (b.x.to_bits(), b.y.to_bits(), b.z.to_bits()),
                    "body {i} differs on {threads} threads"
                );
            }
        }
    }

    /// And it holds through a whole integration, not only one force evaluation — so an
    /// entire run is reproducible whatever the thread count.
    #[test]
    fn a_threaded_run_reproduces_a_sequential_one() {
        let bodies = cloud(200, 7);
        let run = |threads: usize| {
            let mut system = TreeNBody::new("t", &bodies)
                .with_theta(0.5)
                .with_softening(Length::m(1.0))
                .with_threads(threads);
            let mut bus = Exchange::new();
            let dt = Time::s(50.0);
            let mut t = Time::ZERO;
            for _ in 0..20 {
                system.step(t, dt, &mut bus).unwrap();
                t += dt;
            }
            (0..system.count())
                .map(|i| system.body(i).position.to_si())
                .collect::<Vec<_>>()
        };
        let sequential = run(1);
        for threads in [2usize, 4, 8] {
            for (i, (a, b)) in run(threads).iter().zip(sequential.iter()).enumerate() {
                assert_eq!(
                    (a.x.to_bits(), a.y.to_bits(), a.z.to_bits()),
                    (b.x.to_bits(), b.y.to_bits(), b.z.to_bits()),
                    "body {i} after 20 steps on {threads} threads"
                );
            }
        }
    }

    /// The tree converges to the exact answer as the opening angle closes, and at
    /// `theta = 0` it *is* the exact answer.
    ///
    /// Measured against [`NBody`]'s direct summation, which shares no traversal code
    /// with the tree — the same arrangement as the wavefront module against the
    /// analytic one.
    #[test]
    fn the_tree_converges_to_direct_summation() {
        let bodies = cloud(150, 11);
        let positions: Vec<DVec3> = bodies.iter().map(|b| b.position.to_si()).collect();
        let coords = Coords(positions.clone());

        let exact = NBody::new("exact", &bodies)
            .with_softening(Length::m(1.0))
            .acceleration(&coords, Time::ZERO);
        let scale = exact.0.iter().map(|a| a.length()).fold(0.0f64, f64::max);
        assert!(scale > 0.0);

        let worst_for = |theta: f64| {
            let approx = TreeNBody::new("tree", &bodies)
                .with_theta(theta)
                .with_softening(Length::m(1.0))
                .accelerations(&positions);
            approx
                .iter()
                .zip(exact.0.iter())
                .map(|(a, b)| (*a - *b).length() / scale)
                .fold(0.0f64, f64::max)
        };

        // theta = 0 never opens a cell as a monopole, so it reduces to every pair.
        assert!(
            worst_for(0.0) < 1e-12,
            "a closed opening angle should be exact, got {:e}",
            worst_for(0.0)
        );
        // And the error falls monotonically as the angle closes.
        let (wide, working, tight) = (worst_for(1.0), worst_for(0.5), worst_for(0.2));
        assert!(
            wide > working && working > tight,
            "error should fall with theta: {wide:e} > {working:e} > {tight:e}"
        );
        assert!(
            working < 0.02,
            "the conventional theta = 0.5 should be good to a couple of percent, got \
             {working:e}"
        );
    }

    /// What the tree gives up, stated as a measurement rather than a caveat.
    ///
    /// Body `i` sees the tree's approximation of everything else, and body `j` sees a
    /// different one, so their mutual forces no longer cancel and the total is no
    /// longer zero. [`NBody`] holds momentum to 1e-13 of itself over thousands of
    /// steps; this does not, and the gap closes as the opening angle does.
    #[test]
    fn barnes_hut_gives_up_exact_momentum() {
        let bodies = cloud(150, 13);
        let positions: Vec<DVec3> = bodies.iter().map(|b| b.position.to_si()).collect();

        // Sum of |m a| as the scale the residual should be judged against — the same
        // reasoning the kernel's Ledger uses for its own tolerance.
        let residual_for = |theta: f64| {
            let system = TreeNBody::new("tree", &bodies)
                .with_theta(theta)
                .with_softening(Length::m(1.0));
            let acc = system.accelerations(&positions);
            let mut net = DVec3::ZERO;
            let mut scale = 0.0;
            for (a, m) in acc.iter().zip(system.masses.iter()) {
                net += *a * *m;
                scale += (*a * *m).length();
            }
            net.length() / scale
        };

        // Exact pairing at theta = 0: the forces cancel as they should.
        assert!(
            residual_for(0.0) < 1e-14,
            "at theta = 0 the tree is exact and momentum must cancel, got {:e}",
            residual_for(0.0)
        );
        // Open it up and they do not.
        let working = residual_for(0.5);
        assert!(
            working > 1e-6,
            "theta = 0.5 should visibly break the cancellation, got {working:e}"
        );
        assert!(
            working < 1e-2,
            "but it should still be a small fraction of the forces, got {working:e}"
        );
        // And the breakage shrinks with the angle, which is what makes it a knob
        // rather than a defect.
        assert!(
            residual_for(0.2) < working,
            "closing the angle should restore the cancellation"
        );
    }

    /// The quadrupole buys an order in the opening angle, which is the whole reason to
    /// carry six more numbers per cell.
    ///
    /// A cell is not a point. The monopole pretends it is and the error goes as
    /// `(s/d)²`; adding the next term of the expansion pushes that to `(s/d)³`, so the
    /// same accuracy is reached at a wider angle and fewer cells have to be opened.
    /// Measured against direct summation, which knows nothing about either.
    #[test]
    fn the_quadrupole_buys_accuracy_at_the_same_angle() {
        let bodies = cloud(200, 29);
        let positions: Vec<DVec3> = bodies.iter().map(|b| b.position.to_si()).collect();
        let exact = NBody::new("exact", &bodies)
            .with_softening(Length::m(1.0))
            .acceleration(&Coords(positions.clone()), Time::ZERO);
        let scale = exact.0.iter().map(|a| a.length()).fold(0.0f64, f64::max);

        let worst = |theta: f64, quadrupole: bool| {
            let approx = TreeNBody::new("tree", &bodies)
                .with_theta(theta)
                .with_quadrupole(quadrupole)
                .with_softening(Length::m(1.0))
                .accelerations(&positions);
            approx
                .iter()
                .zip(exact.0.iter())
                .map(|(a, b)| (*a - *b).length() / scale)
                .fold(0.0f64, f64::max)
        };

        for theta in [0.3f64, 0.5, 0.8] {
            let (mono, quad) = (worst(theta, false), worst(theta, true));
            assert!(
                quad < mono,
                "at theta = {theta} the quadrupole should beat the monopole: {quad:e} \
                 against {mono:e}"
            );
        }

        // The gain is one power of `d/s`, so it grows as the angle *closes* — which is
        // the opposite of what one might guess. A wide angle accepts cells that are
        // relatively large, and for those the term after the quadrupole is not small
        // either, so the correction fixes a smaller share of what is wrong. Measured:
        // about 6.5 times better at theta = 0.3 and only 2.3 at 0.8.
        let gain = |theta: f64| worst(theta, false) / worst(theta, true);
        assert!(
            gain(0.3) > gain(0.8),
            "the gain should grow as the angle closes: {} at 0.3 against {} at 0.8",
            gain(0.3),
            gain(0.8)
        );
        assert!(
            gain(0.5) > 1.5,
            "a working angle should gain appreciably, got {}",
            gain(0.5)
        );

        // And theta = 0 is still exact either way, since no cell is ever accepted.
        assert!(worst(0.0, true) < 1e-12);
    }

    /// What the quadrupole does *not* fix. The momentum drift comes from two bodies
    /// seeing different approximations of each other, and a better approximation is
    /// still a different one.
    #[test]
    fn the_quadrupole_does_not_restore_exact_momentum() {
        let bodies = cloud(150, 31);
        let positions: Vec<DVec3> = bodies.iter().map(|b| b.position.to_si()).collect();
        let residual = |quadrupole: bool| {
            let system = TreeNBody::new("tree", &bodies)
                .with_theta(0.5)
                .with_quadrupole(quadrupole)
                .with_softening(Length::m(20.0));
            let acc = system.accelerations(&positions);
            let mut net = DVec3::ZERO;
            let mut scale = 0.0;
            for (a, m) in acc.iter().zip(system.masses.iter()) {
                net += *a * *m;
                scale += (*a * *m).length();
            }
            net.length() / scale
        };
        // Smaller, because the forces themselves are better -- but still nothing like
        // the exact solver's cancellation.
        let with = residual(true);
        assert!(with > 1e-9, "still no cancellation, got {with:e}");
        assert!(
            with < residual(false),
            "a better force is still a smaller residual"
        );
    }

    /// The quadrupole is still bit-reproducible across thread counts, which the extra
    /// state and the extra arithmetic could easily have broken.
    #[test]
    fn the_quadrupole_is_still_thread_independent() {
        let bodies = cloud(300, 37);
        let positions: Vec<DVec3> = bodies.iter().map(|b| b.position.to_si()).collect();
        let build = |threads: usize| {
            TreeNBody::new("t", &bodies)
                .with_theta(0.6)
                .with_quadrupole(true)
                .with_threads(threads)
                .accelerations(&positions)
        };
        let reference = build(1);
        for threads in [2usize, 4, 8] {
            for (i, (a, b)) in build(threads).iter().zip(reference.iter()).enumerate() {
                assert_eq!(
                    (a.x.to_bits(), a.y.to_bits(), a.z.to_bits()),
                    (b.x.to_bits(), b.y.to_bits(), b.z.to_bits()),
                    "body {i} differs on {threads} threads"
                );
            }
        }
    }

    /// Real physics still comes out: a two-body orbit through the tree keeps Kepler's
    /// period, because two bodies never trigger the approximation at all.
    #[test]
    fn a_two_body_orbit_is_unaffected_by_the_tree() {
        let central = Mass::from_si(5.9722e24);
        let r = Length::from_si(6_771e3);
        let v = NBody::circular_speed(central, r);
        let period = NBody::orbital_period(central, r);
        let bodies = [
            Body::new(central, LengthVec::ZERO, VelocityVec::ZERO),
            Body::new(
                Mass::kg(1000.0),
                LengthVec::from_si(DVec3::X * r.to_si()),
                VelocityVec::from_si(DVec3::Y * v.to_si()),
            ),
        ];
        let mut system = TreeNBody::new("orbit", &bodies)
            .with_theta(0.5)
            .with_threads(4);
        let steps = 2000;
        let dt = period / steps as f64;
        let mut bus = Exchange::new();
        let mut worst = 0.0f64;
        let mut t = Time::ZERO;
        for _ in 0..steps {
            system.step(t, dt, &mut bus).unwrap();
            t += dt;
            let radius = (system.body(1).position - system.body(0).position).length();
            worst = worst.max((radius / r - 1.0).abs());
        }
        assert!(
            worst < 1e-5,
            "two bodies are always a direct pair, so the orbit should stay round: {worst:e}"
        );
    }

    /// Bodies at the same place do not subdivide forever. The failure mode a naive
    /// octree has is a stack overflow, which is not a failure a physics library should
    /// be able to reach from data.
    #[test]
    fn coincident_bodies_do_not_subdivide_forever() {
        let same = LengthVec::m(10.0, -3.0, 4.0);
        let bodies = [
            Body::new(Mass::kg(1e10), same, VelocityVec::ZERO),
            Body::new(Mass::kg(2e10), same, VelocityVec::ZERO),
            Body::new(Mass::kg(3e10), same, VelocityVec::ZERO),
            Body::new(
                Mass::kg(1e10),
                LengthVec::m(500.0, 0.0, 0.0),
                VelocityVec::ZERO,
            ),
        ];
        let system = TreeNBody::new("pile", &bodies)
            .with_softening(Length::m(1.0))
            .with_threads(2);
        let positions: Vec<DVec3> = bodies.iter().map(|b| b.position.to_si()).collect();
        let acc = system.accelerations(&positions);
        assert_eq!(acc.len(), 4);
        for (i, a) in acc.iter().enumerate() {
            assert!(a.is_finite(), "body {i} got {a}");
        }
        // The distant body is pulled towards the pile, which is the only thing that can
        // be asserted about a configuration this degenerate.
        assert!(acc[3].x < 0.0, "the outlier should be pulled inwards");
    }

    /// An empty or single-body system is not a special case anyone should have to
    /// think about.
    #[test]
    fn degenerate_populations_are_handled() {
        let empty = TreeNBody::new("none", &[]);
        assert!(empty.accelerations(&[]).is_empty());
        assert_eq!(empty.momentum().length().to_si(), 0.0);

        let alone = [Body::new(
            Mass::kg(1e12),
            LengthVec::m(1.0, 2.0, 3.0),
            VelocityVec::m_per_s(0.0, 1.0, 0.0),
        )];
        let one = TreeNBody::new("one", &alone).with_threads(4);
        let acc = one.accelerations(&[DVec3::new(1.0, 2.0, 3.0)]);
        assert_eq!(acc.len(), 1);
        assert!(
            acc[0].length() < 1e-30,
            "nothing pulls on a lone body, got {}",
            acc[0]
        );
    }

    /// Under the kernel's scheduler, with the audit tolerance the approximation earns
    /// rather than the one the exact solver earns — and with the opening angle shown to
    /// be the knob that sets it.
    ///
    /// [`NBody`] is audited at 1e-11 and holds. This cannot be, for the reason in the
    /// module docs, so the honest response is to find out what it *does* hold to and
    /// why that number moves.
    #[test]
    fn the_scheduler_audits_the_tree_at_the_tolerance_theta_earns() {
        use dualis_core::{Schedule, Simulation};

        let bodies = cloud(120, 17);
        // One window, from the same starting configuration every time, so the only
        // thing varying is the opening angle. Softened generously: a random cloud
        // throws up close pairs, and a close pair's tree error swamps the effect being
        // measured. Over a long run the drift also compounds as the cloud collapses
        // into more of them — that is real, and it is a statement about the
        // configuration rather than about theta.
        let drift_at = |theta: f64, tolerance: f64| -> Result<f64, Violation> {
            let mut sim = Simulation::new(Schedule::Multirate)
                .conservation_tolerance(tolerance)
                .with(
                    TreeNBody::new("tree", &bodies)
                        .with_theta(theta)
                        .with_softening(Length::m(20.0))
                        .with_threads(4),
                );
            let before = sim.ledger().get(conserved::MOMENTUM_X).unwrap();
            sim.advance(Time::s(200.0))?;
            let after = sim.ledger().get(conserved::MOMENTUM_X).unwrap();
            Ok((after - before).abs() / before.abs())
        };

        // What the conventional angle actually holds to.
        let working = drift_at(0.5, 1e-2).expect("1e-2 is what theta = 0.5 earns");
        assert!(
            working > 1e-6 && working < 1e-2,
            "theta = 0.5 should drift measurably but not much, got {working:e}"
        );

        // The audit has teeth: ask for the exact solver's standard and it refuses,
        // naming the axis.
        let err = drift_at(0.5, 1e-9).expect_err("1e-9 is not something a tree can hold");
        assert_eq!(err.quantity, "momentum_x");

        // Close the angle and the drift closes with it, which is what makes this a knob
        // rather than a defect to be tolerated.
        let tight = drift_at(0.1, 1e-2).expect("a tighter angle drifts less, not more");
        assert!(
            tight < working / 5.0,
            "closing theta from 0.5 to 0.1 should cut the drift well back: {tight:e} \
             against {working:e}"
        );
        // And at a closed angle the tree is the exact solver, so the drift is only the
        // force-to-acceleration round trip again.
        let exact = drift_at(0.0, 1e-2).expect("theta = 0 is exact");
        assert!(
            exact < 1e-12,
            "a closed angle should restore exact cancellation, got {exact:e}"
        );
    }

    /// The tree is a function of the positions and nothing else: building it twice
    /// gives the same structure, and so the same forces.
    #[test]
    fn the_tree_is_reproducible_across_builds() {
        let bodies = cloud(80, 23);
        let positions: Vec<DVec3> = bodies.iter().map(|b| b.position.to_si()).collect();
        let once = TreeNBody::new("t", &bodies).accelerations(&positions);
        let twice = TreeNBody::new("t", &bodies).accelerations(&positions);
        for (a, b) in once.iter().zip(twice.iter()) {
            assert_eq!(a.x.to_bits(), b.x.to_bits());
            assert_eq!(a.y.to_bits(), b.y.to_bits());
            assert_eq!(a.z.to_bits(), b.z.to_bits());
        }
    }
}
