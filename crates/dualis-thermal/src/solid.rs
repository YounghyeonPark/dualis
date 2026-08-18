//! Conduction through a block, in three dimensions.
//!
//! [`Bar1D`](crate::Bar1D) resolves a gradient along one axis, and that is the right model for a
//! rod, a wall or a beam landing on a strip. A heat sink is not any of those. Heat spreading
//! sideways out of a hot spot is exactly what a fin, a spreader plate and a mounting boss are
//! *for*, and a one-dimensional model cannot show it: it has nowhere for the heat to go but
//! along.
//!
//! # What the third dimension costs
//!
//! Twice, and the second time is the one that hurts.
//!
//! A block of `n` cells on a side is `n³` cells rather than `n`. And the explicit stability
//! limit tightens with each axis, because the limit is on the *sum* of what the three
//! directions do in one step:
//!
//! ```text
//!   1D      α·dt/dx² ≤ 1/2
//!   2D      α·dt/dx² ≤ 1/4
//!   3D      α·dt/dx² ≤ 1/6
//! ```
//!
//! So a 3D block at the same spacing takes three times as many steps as a bar, each of them n²
//! times more work. `Room` in `dualis-acoustic` records the same trade for the wave equation and
//! reaches a factor of √3 rather than 3, because a wave's limit is on the wave speed and a
//! diffusion limit is on the diffusivity — one is linear in the sum and the other in its square
//! root.
//!
//! That is why [`LumpedMass`](crate::LumpedMass) and [`ThermalNetwork`](crate::ThermalNetwork)
//! are not going anywhere. In aluminium at a millimetre the step is **2.41 ms**, so a motor
//! housing over its two-thousand-second thermal time constant is **828,000 steps** — each of
//! them a sweep over however many cells the housing is, which at that spacing is around a
//! million. A graph of four nodes answers the same question immediately. Use the cheapest model
//! whose reduction still holds, and `LumpedMass::biot_number` is how you find out whether it
//! does.
//!
//! What a block is *for* is the case where the reduction does not hold: a hot spot, a spreader,
//! a gradient across a joint. Those are questions a lumped model cannot answer at any price.
//!
//! # What it is checked against
//!
//! A separable cosine mode is an **exact eigenvector of the discrete operator**, not merely an
//! approximate solution of the continuum one. On a cell-centred grid with insulated faces the
//! mode `cos(aπ(i+½)/nx)·cos(bπ(j+½)/ny)·cos(cπ(k+½)/nz)` decays by exactly the same factor
//! every step, and that factor is known in closed form. So the test is an equality at machine
//! precision rather than a tolerance, and it is sensitive to a swapped axis, a dropped term or a
//! wrong spacing in a way a smooth decaying blob would not be.
//!
//! The *continuum* rate is then a second test: the discrete eigenvalue approaches
//! `−α·π²·(1/Lx² + 1/Ly² + 1/Lz²)` at second order, and refining the grid quarters the error.
//! Rate rather than value, because a scheme that is first order where it claims to be second is
//! the defect this workspace has already shipped once.
//!
//! # A workaround that is gone
//!
//! This domain briefly offered `Domain::as_bodies` alongside `as_field`, so a viewer could get
//! its cells as a point cloud. That was not a design choice; it was cover for `dualis-scene`'s
//! `Extent` being two-dimensional, which would have captured a block as its `z = 0` face.
//!
//! `Extent` is three-dimensional now, so the cover is unnecessary — and it was never free. A
//! domain that is two shapes at once makes the picture depend on whether somebody remembered to
//! set an extent, which is a mode nothing announces. It is a field, and only a field.

use dualis_core::conserved::quantity;
use dualis_core::Reading;
use dualis_core::{
    units::{
        Area, Conductance, Energy, HeatCapacity, Length, LengthVec, Temperature, Time, Volume,
    },
    Domain, Exchange, Ledger, ScalarField, Substance, Violation,
};
use glam::DVec3;

use crate::{Environment, HEAT};
use dualis_units::STEFAN_BOLTZMANN;
use std::collections::BTreeMap;

/// The largest Fourier number an explicit three-dimensional sweep is stable at.
///
/// `1/(2d)` for `d` dimensions, from requiring the amplification factor of the worst-resolved
/// mode to stay inside the unit circle. Public because a caller sizing a grid needs it before
/// there is anything to ask.
pub const STABLE_FOURIER_3D: f64 = 1.0 / 6.0;

/// A rectangular block, conducting in three dimensions, of one material or of several.
///
/// Cells are **cubes** of a single spacing rather than boxes of three. That is a deliberate
/// restriction and the same one `Room::of_air` makes: an anisotropic cell makes the stability
/// limit anisotropic and the truncation error different along each axis, so the grid would be
/// resolving one direction better than another for a reason that had nothing to do with the
/// physics. A block that is longer than it is thick is more cells along, not longer cells.
///
/// Faces are **insulated**. Every boundary cell exchanges with the neighbours it has and no
/// others, which is what makes the total exactly conserved rather than conserved to a tolerance.
/// A block that should lose heat gets that from a domain on the other side of the bus.
///
/// [`fill`](Solid3D::fill) gives cells a different substance from the constructor's, which is what
/// makes a coating, a joint or a layered wall expressible. Read that method for the one number the
/// scheme turns on — the conductivity **on a face** is the harmonic mean of the two cells', not
/// A surface conductance put in series with the half cell of solid behind it.
///
/// `1/(1/G_film + dx/(2kA))`. The film is charged against the **surface**, and a finite-volume
/// cell knows only its centre, so the half cell between them is part of the path. Leaving it
/// out sheds too much heat and — measured — turns a second-order boundary into a first-order
/// one, which is the defect this workspace has twice found in an acoustic wall.
///
/// An infinite conductivity (a substance that does not say) leaves the film alone, which is the
/// right limit rather than a special case.
fn series_with_half_cell(film: f64, conductivity: f64, area: f64, dx: f64) -> f64 {
    if film <= 0.0 {
        return 0.0;
    }
    let half = 2.0 * conductivity * area / dx;
    if !half.is_finite() {
        return film;
    }
    if half <= 0.0 {
        return 0.0;
    }
    1.0 / (1.0 / film + 1.0 / half)
}

/// One outer face of a block.
///
/// Named rather than indexed because a caller says which side of a part is exposed, and
/// `faces[3]` is a thing nobody can read back. The axis pairs are low and high along `x`, `y`
/// and `z` in the block's own coordinates — a [`Pose`](dualis_core::Pose) is what puts those in
/// the world.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Face {
    /// The `x = 0` face.
    XMin,
    /// The far face along `x`.
    XMax,
    /// The `y = 0` face.
    YMin,
    /// The far face along `y`.
    YMax,
    /// The `z = 0` face.
    ZMin,
    /// The far face along `z`.
    ZMax,
}

impl Face {
    /// All six, in a fixed order — for a block exposed on every side.
    pub const ALL: [Face; 6] = [
        Face::XMin,
        Face::XMax,
        Face::YMin,
        Face::YMax,
        Face::ZMin,
        Face::ZMax,
    ];

    /// Whether the cell at `(i, j, k)` of a `counts`-shaped block lies on this face.
    fn holds(&self, (i, j, k): (usize, usize, usize), counts: (usize, usize, usize)) -> bool {
        let (nx, ny, nz) = counts;
        match self {
            Face::XMin => i == 0,
            Face::XMax => i + 1 == nx,
            Face::YMin => j == 0,
            Face::YMax => j + 1 == ny,
            Face::ZMin => k == 0,
            Face::ZMax => k + 1 == nz,
        }
    }
}

/// the arithmetic one, and the difference is not a refinement away.
#[derive(Clone, Debug)]
pub struct Solid3D {
    name: String,
    /// Every substance in the block. Index 0 is the constructor's, and a block nobody has
    /// [`fill`](Solid3D::fill)ed holds only that one.
    materials: Vec<Substance>,
    /// Which of `materials` each cell is, indexed as `cells` is.
    which: Vec<u32>,
    /// Cell-centre temperatures, indexed `x + nx*(y + ny*z)`.
    cells: Vec<f64>,
    saved: Vec<f64>,
    counts: (usize, usize, usize),
    dx: Length,
    absorbed: f64,
    /// What [`stored_heat`](Solid3D::stored_heat) is measured from — see [`Bar1D`](crate::Bar1D)
    /// for why an enthalpy reference is chosen for precision rather than for physics.
    reference: f64,
    /// Face conductivities in W/m/K, rebuilt by [`resolve`](Solid3D::resolve): `kx` is indexed
    /// `i + (nx+1)*(j + ny*k)` and holds the face between cells `i-1` and `i`, so `kx[0]` and
    /// `kx[nx]` are the outer faces and are **zero** — which is the insulated boundary, stated in
    /// the conductance rather than in the stencil's index arithmetic.
    kx: Vec<f64>,
    ky: Vec<f64>,
    kz: Vec<f64>,
    /// Per-cell `dx / C_i`, which is what multiplies `Σ k_f ΔT` to give a rate of temperature.
    mobility: Vec<f64>,
    /// Per-cell heat capacity `ρ_i c_i dx³`, in J/K. `NaN` where a substance does not say.
    capacity: Vec<f64>,
    /// Per-cell melting point in kelvin, or `+∞` for a substance that does not melt — so the phase
    /// branch is a comparison that is simply never true rather than an `Option` to unwrap per cell.
    melt_point: Vec<f64>,
    /// Per-cell latent heat in **joules**: `ρ L dx³`. Zero where there is no phase change, which is
    /// what the sweep tests.
    ///
    /// # It was kelvin, and two phases broke that
    ///
    /// `L/c_p` — 163 K for ice — kept the whole update in one unit and the numbers `O(1)` rather than
    /// `O(3e5)`, which is the reference-point argument `runtime/gpu` paid 1660× for. It is correct
    /// while a cell has **one** specific heat.
    ///
    /// With two it is not, and it fails quietly. The kelvin figure is normalised by the *solid's*
    /// `c_p` and then multiplied by the cell's *current* capacity, which for a liquid cell is the
    /// liquid's — so freezing a cell of water cost `4182/2050` = **2.04×** what it should, and the
    /// front came out **27% short** with nothing pointing at the latent heat. Joules do not have a
    /// normalisation to get wrong.
    latent: Vec<f64>,
    /// Per-cell capacity of the **solid** phase, `ρ c_s dx³`, in J/K.
    ///
    /// The enthalpy map needs both capacities by name: `C_s` below the melting point and `C_l` above
    /// it. `capacity` is the *mixed* one and is what the sweep divides a flux by, which is right —
    /// a fully solid or fully liquid cell has its own, and a mushy cell does not change temperature
    /// at all, so what it would have divided by never matters.
    cap_solid: Vec<f64>,
    /// Per-cell capacity of the **liquid** phase, `ρ c_l dx³`. Equal to `cap_solid` under the
    /// one-phase model.
    cap_liquid: Vec<f64>,
    /// How much of each cell has melted, `0` solid and `1` liquid. Derived from the temperature by
    /// [`resolve`](Solid3D::resolve) and then carried by the sweep.
    melted: Vec<f64>,
    /// The phase at the last [`checkpoint`](Domain::checkpoint), alongside `saved`.
    saved_melted: Vec<f64>,
    /// `max_i dx·Σ_f k_f / C_i`, the reciprocal of the largest stable step. See
    /// [`max_stable_dt`](Domain::max_stable_dt) for why it is a maximum over cells and not a
    /// function of the fastest material.
    worst_rate: f64,
    /// Cells that are **not part of the block** — see [`Solid3D::empty`].
    void: Vec<bool>,
    /// Pairs of solid cells that face each other across a run of void, with the radiative
    /// exchange coefficient between them. Rebuilt whenever the void or the materials change.
    ///
    /// `(a, b, coefficient)` where the coefficient is `σA/(1/ε₁ + 1/ε₂ − 1)`, so the exchange is
    /// `coefficient · (Tₐ⁴ − T_b⁴)` — see [`Solid3D::empty`] for what this models and what it
    /// does not.
    gaps: Vec<(usize, usize, f64)>,
    /// Each cell's conduction sum, `Σ_f k_f`, kept so the stability limit can be rebuilt at
    /// the block's **current** temperature — see [`Solid3D::worst_rate_now`].
    face_sum: Vec<f64>,
    /// Which faces lose heat, and to what. Empty is the adiabatic block every scene had until
    /// now: six insulated faces and no steady state to reach.
    exposed: BTreeMap<Face, Environment>,
    /// Joules given up to those environments over the run, and the reason the books still
    /// balance: the ledger is `stored + lost`, so what leaves the cells arrives in this
    /// counter and the total moves only by what crossed the bus. `LumpedMass` keeps the same
    /// pair for the same reason.
    lost: f64,
    /// The saved counterpart of [`Solid3D::lost`].
    ///
    /// Saved because `LumpedMass` learned this the expensive way: rewinding the cells and not
    /// the losses makes an iterative sweep grow its books by one sweep of shed heat per
    /// iteration, and the audit reports energy created from nothing. Measured there at 1567.6 J
    /// becoming 1600.9 J over forty advances.
    saved_lost: f64,
    /// Whether any material names a **liquid** phase whose properties differ from its solid.
    ///
    /// The one thing that decides whether the sweep has to rebuild its operator each step. A block of
    /// ice does; a block of aluminium, or of ice under the one-phase model, does not.
    two_phase: bool,
}

/// One patch of a clearance: two facing surfaces and how far apart they are.
///
/// A gap is found pair by pair, along one axis, one cell face against the one directly across
/// from it. A *patch* is what those pairs add up to — a connected sheet of facing area at one
/// separation, which is the shape a view factor is a statement about.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GapPatch {
    /// How many facing cell pairs are in it.
    pub pairs: usize,
    /// The patch's extent across the gap, in metres — the bounding box of the facing cells.
    pub span: (f64, f64),
    /// How far the two surfaces are apart, in metres.
    pub distance: f64,
    /// Whether the facing cells fill their bounding box.
    ///
    /// [`view_factor`](GapPatch::view_factor) is exact for a rectangle and an **upper bound** for
    /// anything else, because a patch with a bite out of it sees less of itself than its bounding
    /// box does. Reported rather than corrected: an L-shaped clearance is not a shape the closed
    /// form covers, and a number invented for it would be worth less than knowing it is a bound.
    pub rectangular: bool,
}

impl GapPatch {
    /// The view factor from one of the facing surfaces to the other, `F₁₂`.
    ///
    /// The exact closed form for two equal, parallel, directly-opposed rectangles — which is the
    /// geometry the pairing produces by construction, since a pair is a cell and the cell across
    /// from it. With `X = a/c` and `Y = b/c`:
    ///
    /// ```text
    /// F = 2/(πXY) · [ ln √((1+X²)(1+Y²)/(1+X²+Y²))
    ///                 + X√(1+Y²)·atan(X/√(1+Y²)) + Y√(1+X²)·atan(Y/√(1+X²))
    ///                 − X·atan X − Y·atan Y ]
    /// ```
    ///
    /// **This is not what the exchange is charged**, and the difference is the point. What the
    /// block charges is `F̄ = 1`, which is exact when the sides of the gap are mirrors — and the
    /// block's own outer faces are exactly that, because an insulated boundary is implemented as
    /// a mirror and a mirror extends the two surfaces to infinity. `F₁₂` is what the same pair
    /// would exchange with **nothing** at the sides, open to space. The two agree when the gap is
    /// narrow compared with the surfaces and diverge without limit when it is not, so this is the
    /// number that says how much of the answer is resting on reading the boundary one way.
    pub fn view_factor(&self) -> f64 {
        let (x, y) = (self.span.0 / self.distance, self.span.1 / self.distance);
        // A patch with no extent sees nothing; two surfaces with nothing between them see each
        // other entirely. Both are limits of the form rather than sentinels, and both are here
        // because this is public: a caller can build a `GapPatch` this block would never produce.
        if !self.distance.is_finite() || self.span.0.is_nan() || self.span.1.is_nan() {
            return 0.0;
        }
        if self.distance <= 0.0 {
            return if self.span.0 > 0.0 && self.span.1 > 0.0 {
                1.0
            } else {
                0.0
            };
        }
        if x <= 0.0 || y <= 0.0 {
            return 0.0;
        }
        let (x2, y2) = (x * x, y * y);
        // **The closed form cancels itself away for a distant patch.** The bracket below is a
        // small difference of terms of order `X`, and the answer is of order `X²` — so at `X` of
        // 1e-3 four digits are already gone and the ratio to the exact `A/πc²` limit *turns
        // around*: measured 0.999993 at X = 1/300, 1.00009 at 1/1000 and 1.007 at 1/3000, growing
        // where it should be converging. A function that is quietly worse the further apart two
        // surfaces are is the shape of defect this repository looks for, so the small-argument
        // regime gets the series instead, `F = XY/π · (1 − (X²+Y²)/3 + O(X⁴))`, which is exact to
        // 1e-8 either side of the crossover.
        if x < 0.01 && y < 0.01 {
            return (x * y / std::f64::consts::PI * (1.0 - (x2 + y2) / 3.0)).clamp(0.0, 1.0);
        }
        let (rx, ry) = ((1.0 + x2).sqrt(), (1.0 + y2).sqrt());
        let t = ((1.0 + x2) * (1.0 + y2) / (1.0 + x2 + y2)).sqrt().ln()
            + x * ry * (x / ry).atan()
            + y * rx * (y / rx).atan()
            - x * x.atan()
            - y * y.atan();
        (2.0 / (std::f64::consts::PI * x * y) * t).clamp(0.0, 1.0)
    }
}

impl Solid3D {
    /// A block of `counts` cubic cells of side `dx`, all starting at `initial`.
    ///
    /// Each count is forced to at least one. A block one cell thick in two directions is a
    /// legitimate thing to ask for and reduces exactly to a bar, which is how the closed-form
    /// tests check the three axes against each other.
    pub fn new(
        name: impl Into<String>,
        substance: Substance,
        counts: (usize, usize, usize),
        dx: Length,
        initial: Temperature,
    ) -> Solid3D {
        let counts = (counts.0.max(1), counts.1.max(1), counts.2.max(1));
        let cells = vec![initial.to_si(); counts.0 * counts.1 * counts.2];
        let mut block = Solid3D {
            name: name.into(),
            materials: vec![substance],
            which: vec![0; cells.len()],
            saved: cells.clone(),
            cells,
            counts,
            dx,
            absorbed: 0.0,
            reference: initial.to_si(),
            kx: Vec::new(),
            ky: Vec::new(),
            kz: Vec::new(),
            mobility: Vec::new(),
            capacity: Vec::new(),
            melt_point: Vec::new(),
            latent: Vec::new(),
            melted: Vec::new(),
            saved_melted: Vec::new(),
            cap_solid: Vec::new(),
            cap_liquid: Vec::new(),
            worst_rate: 0.0,
            void: vec![false; counts.0 * counts.1 * counts.2],
            gaps: Vec::new(),
            face_sum: Vec::new(),
            exposed: BTreeMap::new(),
            lost: 0.0,
            saved_lost: 0.0,
            two_phase: false,
        };
        block.resolve();
        block
    }

    /// How much of one cell has melted: `0` entirely solid, `1` entirely liquid.
    ///
    /// Always zero for a substance with no [`FusionProps`](dualis_core::substance::FusionProps), which is every
    /// entry in the catalogue but [`Substance::ice`].
    pub fn melted_fraction_at(&self, i: usize, j: usize, k: usize) -> f64 {
        let (nx, ny, nz) = self.counts;
        self.melted[i.min(nx - 1) + nx * (j.min(ny - 1) + ny * k.min(nz - 1))]
    }

    /// Declare how much of a cell has melted, for an initial condition a temperature cannot express.
    ///
    /// The counterpart to [`set_temperature`](Solid3D::set_temperature), and it exists because **a
    /// temperature is not a state** for a substance that melts: 0 °C is ice, water, or any mixture,
    /// and Stefan's problem starts from liquid at exactly the melting point. Without this the only
    /// way to say that is to start a hair above it, which puts sensible heat in the initial condition
    /// that the closed form does not have.
    ///
    /// Ignored for a substance that does not melt, and out of range is ignored, matching
    /// `set_temperature`.
    ///
    /// # Supercooling is not representable, and this keeps it that way
    ///
    /// The state is one monotone number, so a fraction and a temperature cannot disagree: asking for
    /// liquid raises the temperature to at least the melting point and asking for solid lowers it to
    /// at most. Liquid below freezing and solid above it are real states of real matter and this
    /// model does not have them — the sharp-interface problem assumes the interface is *at* the
    /// melting point, which is what makes Neumann's solution its solution.
    pub fn set_melted_fraction(&mut self, i: usize, j: usize, k: usize, fraction: f64) {
        let Some(idx) = self.index(i, j, k) else {
            return;
        };
        if self.latent[idx] <= 0.0 {
            return;
        }
        let phi = fraction.clamp(0.0, 1.0);
        self.melted[idx] = phi;
        let point = self.melt_point[idx];
        if phi >= 1.0 {
            self.cells[idx] = self.cells[idx].max(point);
        } else if phi <= 0.0 {
            self.cells[idx] = self.cells[idx].min(point);
        } else {
            self.cells[idx] = point;
        }
    }

    /// The volume that has melted, summed over every cell.
    ///
    /// The integral quantity, and the one worth reading rather than a front *position*: a front is
    /// only a position if the problem is one-dimensional, and this is the same number in three. For
    /// a column of cells it gives the position anyway — `melted_volume / area` is how far in the
    /// front has reached, including the partial cell it is currently inside, which is what makes a
    /// front measurable to better than one cell.
    pub fn melted_volume(&self) -> Volume {
        Volume::from_si(self.melted.iter().sum::<f64>() * self.cell_volume())
    }

    /// Put a different substance in every cell the predicate accepts.
    ///
    /// Composable: call it once per layer, per coating, per inclusion. Cells nobody claims keep the
    /// constructor's substance, so a block is never partly undefined.
    ///
    /// # The harmonic mean, and why it is not a detail
    ///
    /// Heat crosses a face, and a face has two materials touching it. The conductance of the half
    /// cell either side is in **series**, so the conductivity that governs the face is
    ///
    /// ```text
    ///   k_face = 2 k_L k_R / (k_L + k_R)
    /// ```
    ///
    /// the harmonic mean, and the arithmetic mean `(k_L + k_R)/2` is not an alternative convention
    /// — it is wrong. With aluminium against borosilicate the two differ by **38×** (2.21 against
    /// 84.1 W/m/K), and the arithmetic one short-circuits the interface: a wall it models has 4.2%
    /// less resistance than its own layers add up to at 24 cells, and reaching the harmonic answer
    /// to 0.1% would take about a thousand.
    ///
    /// The harmonic mean is not merely better. With the material interface on a cell face it makes
    /// the discrete series resistance **exactly** `Σ Lᵢ/(kᵢA)` at every resolution, which is why
    /// `a_layered_wall.rs` is an equality and not a tolerance.
    ///
    /// "On a cell face" is the whole condition, and it is one this method cannot break: cells are
    /// whole, so an interface is always on a face. A scheme that placed a layer boundary partway
    /// through a cell would be first order there whichever mean it used, because the cell would be a
    /// mixture and no single conductivity describes one.
    ///
    /// # What it costs
    ///
    /// Five `f64` and a `u32` per cell beyond the temperatures: three face conductivities, a
    /// capacity, its reciprocal scaled by `dx`, and which substance the cell is. All rebuilt here
    /// and none of it in the sweep, which is the trade — the alternative is six harmonic means and
    /// six divisions per cell per step.
    ///
    /// A uniform block pays the same. That is deliberate: a fast path for one material would be a
    /// second implementation of the same physics, exercised by the tests that came before this
    /// method and by nothing after it.
    pub fn fill(
        &mut self,
        substance: Substance,
        which: impl Fn(usize, usize, usize) -> bool,
    ) -> &mut Solid3D {
        let id = match self.materials.iter().position(|s| *s == substance) {
            Some(at) => at as u32,
            None => {
                self.materials.push(substance);
                (self.materials.len() - 1) as u32
            }
        };
        let (nx, ny, nz) = self.counts;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if which(i, j, k) {
                        self.which[i + nx * (j + ny * k)] = id;
                    }
                }
            }
        }
        self.resolve();
        self
    }

    /// The substance in one cell. Out of range reads the nearest one in range.
    pub fn substance_at(&self, i: usize, j: usize, k: usize) -> &Substance {
        let (nx, ny, nz) = self.counts;
        let idx = i.min(nx - 1) + nx * (j.min(ny - 1) + ny * k.min(nz - 1));
        &self.materials[self.which[idx] as usize]
    }

    /// How many distinct substances are in the block. One until something has [`fill`]ed it.
    ///
    /// [`fill`]: Solid3D::fill
    pub fn substances(&self) -> usize {
        self.materials.len()
    }

    /// The conductance of the face between two cells, or `None` if they are not face neighbours.
    ///
    /// `k_face · A / dx`, which for cubic cells is `k_face · dx`. This is the number the sweep
    /// actually uses, so a caller checking a joint against `Σ Lᵢ/(kᵢA)` by hand is checking the
    /// same arithmetic the march does rather than a restatement of it.
    pub fn face_conductance(
        &self,
        a: (usize, usize, usize),
        b: (usize, usize, usize),
    ) -> Option<Conductance> {
        let (nx, ny, _) = self.counts;
        self.index(a.0, a.1, a.2)?;
        self.index(b.0, b.1, b.2)?;
        let step = |p: usize, q: usize| (p as isize - q as isize).abs();
        let (di, dj, dk) = (step(a.0, b.0), step(a.1, b.1), step(a.2, b.2));
        let k = match (di, dj, dk) {
            (1, 0, 0) => self.kx[a.0.max(b.0) + (nx + 1) * (a.1 + ny * a.2)],
            (0, 1, 0) => self.ky[a.0 + nx * (a.1.max(b.1) + (ny + 1) * a.2)],
            (0, 0, 1) => self.kz[a.0 + nx * (a.1 + ny * a.2.max(b.2))],
            _ => return None,
        };
        Some(Conductance::from_si(k * self.dx.to_si()))
    }

    /// Rebuild everything derived from the materials: face conductivities, mobilities, the limit.
    ///
    /// Called by the constructor and by every [`fill`](Solid3D::fill), on the rule this workspace
    /// learned from `Puck::repack` — a mutator that leaves a cached solve stale reports a number
    /// that was true about the previous object, and it looks exactly like a number.
    fn resolve(&mut self) {
        let (nx, ny, nz) = self.counts;
        let dx = self.dx.to_si();
        let volume = Volume::from_si(dx * dx * dx);

        // Per material, the **solid** pair and the **liquid** pair. A material whose `FusionProps`
        // names no liquid uses the solid pair for both, which is the one-phase model — exact whenever
        // the liquid sits at the melting point, because a face with no temperature difference across
        // it carries no heat whatever its conductivity.
        //
        // `NaN` for a substance that does not say, which `step` refuses on rather than stepping with
        // a plausible default.
        let props: Vec<((f64, f64), (f64, f64))> = self
            .materials
            .iter()
            .map(|s| {
                let solid = (
                    s.thermal.map_or(f64::NAN, |t| t.conductivity.to_si()),
                    s.heat_capacity(volume).map_or(f64::NAN, |c| c.to_si()),
                );
                let liquid = match s.fusion.and_then(|f| f.liquid) {
                    Some(t) => (
                        t.conductivity.to_si(),
                        s.mass_of(volume).to_si() * t.specific_heat.to_si(),
                    ),
                    None => solid,
                };
                (solid, liquid)
            })
            .collect();

        // A liquid phase that differs from the solid is what makes the operator move with the front.
        // A material whose `liquid` is absent, or identical to its solid, needs no per-step rebuild.
        self.two_phase =
            self.materials
                .iter()
                .any(|s| match (s.fusion.and_then(|f| f.liquid), s.thermal) {
                    (Some(l), Some(t)) => {
                        l.conductivity != t.conductivity || l.specific_heat != t.specific_heat
                    }
                    _ => false,
                });

        // The phase state first, because everything below depends on it.
        //
        // Rebuilt from temperature only when there is not one yet. Once the sweep is carrying a
        // fraction, re-deriving it would throw away every partially melted cell — and this runs every
        // step for a two-phase block, because a cell's conductivity moves with its fraction.
        let fresh = self.melted.len() != self.cells.len();
        let carried = std::mem::take(&mut self.melted);
        self.melt_point = Vec::with_capacity(self.cells.len());
        self.latent = Vec::with_capacity(self.cells.len());
        self.melted = Vec::with_capacity(self.cells.len());
        for (c, prior) in carried
            .iter()
            .copied()
            .chain(std::iter::repeat(0.0))
            .take(self.cells.len())
            .enumerate()
        {
            let s = &self.materials[self.which[c] as usize];
            let (point, latent) = match (s.fusion, s.thermal) {
                (Some(f), Some(_)) => (
                    f.melting_point.to_si(),
                    s.latent_energy(volume).map_or(0.0, |e| e.to_si()),
                ),
                // A substance with fusion but no thermal properties cannot be stepped at all: its
                // capacity is `NaN` and the sweep refuses rather than guessing.
                _ => (f64::INFINITY, 0.0),
            };
            self.melt_point.push(point);
            self.latent.push(latent);
            self.melted.push(if fresh {
                if self.cells[c] > point {
                    1.0
                } else {
                    0.0
                }
            } else {
                prior
            });
        }
        // Only on a fresh build. Doing it every time would overwrite a checkpoint, and `resolve` is
        // called from the sweep now.
        if fresh {
            self.saved_melted.clone_from(&self.melted);
        }

        // **Mixed by melt fraction, arithmetically, and that choice is the one to argue about.**
        //
        // A half-melted cell is half of each, and across a cell the two phases sit side by side
        // rather than in series — heat crossing the mush passes through both at once, so a parallel
        // mean is what matches. The *face* between two cells stays the harmonic mean of whatever the
        // two cells came out as, because a face **is** series. Mixing one way and joining the other
        // is not an inconsistency; they are different geometries.
        //
        // It applies only inside the mush, one or two cells wide, and that is where the scheme's
        // remaining first-order error lives.
        let mixed: Vec<(f64, f64)> = (0..self.cells.len())
            .map(|c| {
                // Nothing conducts nothing and holds nothing. Zeroing here rather than at every
                // reader is what makes void fall out of the rest of the arithmetic: the face
                // mean is already the harmonic one and guarded at zero, so a face touching void
                // carries zero without knowing what void is.
                if self.void[c] {
                    return (0.0, 0.0);
                }
                let (solid, liquid) = props[self.which[c] as usize];
                let phi = self.melted[c];
                (
                    solid.0 + phi * (liquid.0 - solid.0),
                    solid.1 + phi * (liquid.1 - solid.1),
                )
            })
            .collect();

        // The series mean of the two half cells. Written as `2ab/(a+b)` and guarded at zero, because
        // a perfect insulator is a legitimate fill and `0/0` is not a boundary condition.
        let series = |a: f64, b: f64| {
            let sum = a + b;
            if sum > 0.0 {
                2.0 * a * b / sum
            } else {
                0.0
            }
        };
        let k_of = |cell: usize| mixed[cell].0;

        self.kx = vec![0.0; (nx + 1) * ny * nz];
        self.ky = vec![0.0; nx * (ny + 1) * nz];
        self.kz = vec![0.0; nx * ny * (nz + 1)];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let c = i + nx * (j + ny * k);
                    if i > 0 {
                        self.kx[i + (nx + 1) * (j + ny * k)] = series(k_of(c - 1), k_of(c));
                    }
                    if j > 0 {
                        self.ky[i + nx * (j + (ny + 1) * k)] = series(k_of(c - nx), k_of(c));
                    }
                    if k > 0 {
                        self.kz[i + nx * (j + ny * k)] = series(k_of(c - nx * ny), k_of(c));
                    }
                }
            }
        }

        self.cap_solid = (0..self.cells.len())
            .map(|c| props[self.which[c] as usize].0 .1)
            .collect();
        self.cap_liquid = (0..self.cells.len())
            .map(|c| props[self.which[c] as usize].1 .1)
            .collect();
        self.capacity = vec![0.0; self.cells.len()];
        self.mobility = vec![0.0; self.cells.len()];
        self.face_sum = vec![0.0; self.cells.len()];
        self.worst_rate = 0.0;
        let mut nowhere = false;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let c = i + nx * (j + ny * k);
                    let cap = mixed[c].1;
                    self.capacity[c] = cap;
                    // A void cell has no capacity, and `dx/0` is an infinity that would poison
                    // the limit and the sweep alike. Zero mobility is the honest value: nothing
                    // moves in a cell that holds nothing.
                    self.mobility[c] = if self.void[c] { 0.0 } else { dx / cap };
                    let sum = self.kx[i + (nx + 1) * (j + ny * k)]
                        + self.kx[i + 1 + (nx + 1) * (j + ny * k)]
                        + self.ky[i + nx * (j + (ny + 1) * k)]
                        + self.ky[i + nx * (j + 1 + (ny + 1) * k)]
                        + self.kz[i + nx * (j + ny * k)]
                        + self.kz[i + nx * (j + ny * (k + 1))];
                    // NaN loses every comparison, so `max` would step straight past a substance that
                    // does not say what it conducts. Tested for, and the sweep is refused.
                    //
                    // An exposed cell also loses to air, and that conductance leaves the cell
                    // exactly as a face conductance does — so it belongs in the same sum. It is
                    // divided by `dx` because `mobility` is `dx/capacity`: the face terms are
                    // conductivities and this one is already a conductance.
                    self.face_sum[c] = sum;
                    let rate = self.mobility[c] * sum;
                    if rate.is_nan() {
                        nowhere = true;
                    } else {
                        self.worst_rate = self.worst_rate.max(rate);
                    }
                }
            }
        }
        if nowhere {
            self.worst_rate = f64::NAN;
        }
        self.find_gaps();
    }

    /// Which solid cells face each other across void, and how strongly they radiate.
    ///
    /// Along each grid line, a run of void with solid at both ends is a gap, and the two cells
    /// at its ends see each other. Diagonals are not paired: a face is what radiates, and the
    /// grid's faces are axis-aligned.
    fn find_gaps(&mut self) {
        self.gaps.clear();
        if !self.void.iter().any(|v| *v) {
            return;
        }
        let (nx, ny, nz) = self.counts;
        let dx = self.dx.to_si();
        let area = dx * dx;
        let emissivity = |c: usize| {
            self.materials[self.which[c] as usize]
                .thermal
                .map_or(0.0, |t| t.emissivity)
        };

        // One pass per axis. `line` walks the indices of a single row, column or pillar.
        let mut walk = |line: Vec<usize>| {
            let mut last_solid: Option<usize> = None;
            let mut void_between = false;
            for c in line {
                if self.void[c] {
                    void_between = last_solid.is_some();
                    continue;
                }
                if let (Some(a), true) = (last_solid, void_between) {
                    let (ea, eb) = (emissivity(a), emissivity(c));
                    // The parallel-plate series: `1/ε₁ + 1/ε₂ − 1`. A surface that does not
                    // radiate at all makes the pair carry nothing, which is the right limit and
                    // also keeps the reciprocal finite.
                    if ea > 0.0 && eb > 0.0 {
                        let resistance = 1.0 / ea + 1.0 / eb - 1.0;
                        if resistance > 0.0 {
                            self.gaps
                                .push((a, c, STEFAN_BOLTZMANN.to_si() * area / resistance));
                        }
                    }
                }
                last_solid = Some(c);
                void_between = false;
            }
        };

        for k in 0..nz {
            for j in 0..ny {
                walk((0..nx).map(|i| i + nx * (j + ny * k)).collect());
            }
        }
        for k in 0..nz {
            for i in 0..nx {
                walk((0..ny).map(|j| i + nx * (j + ny * k)).collect());
            }
        }
        for j in 0..ny {
            for i in 0..nx {
                walk((0..nz).map(|k| i + nx * (j + ny * k)).collect());
            }
        }
    }

    /// The clearances in this block, as sheets of facing area rather than as cell pairs.
    ///
    /// A pair is what the exchange is computed on; a patch is what a **view factor** is a
    /// statement about, and the two are not the same object. Grouped by axis and separation, then
    /// by connectivity across the gap, so two unrelated clearances at the same width do not
    /// average into one patch that describes neither.
    ///
    /// Empty for a block with no void, which is every block that existed before there was one.
    pub fn gap_patches(&self) -> Vec<GapPatch> {
        let (nx, ny, _) = self.counts;
        let dx = self.dx.to_si();
        // A pair's axis and separation are recoverable from the two indices: the stride between
        // them says which way it ran, and the multiple of that stride says how far.
        let strides = [(0usize, 1usize), (1, nx), (2, nx * ny)];
        let mut sheets: BTreeMap<(usize, usize), Vec<(usize, usize)>> = BTreeMap::new();
        for (a, b, _) in &self.gaps {
            let delta = b - a;
            let Some(&(axis, stride)) = strides.iter().find(|(ax, st)| {
                delta % st == 0
                    && delta / st > 1
                    && match ax {
                        0 => delta < nx,
                        1 => delta < nx * ny,
                        _ => true,
                    }
            }) else {
                continue;
            };
            // The **void** between them, not the distance between their centres. A view factor is
            // a statement about two surfaces, and the surfaces are the faces bounding the gap: a
            // one-cell clearance puts them one cell apart while their centres are two.
            let cells = delta / stride - 1;
            let (i, j, k) = (a % nx, (a / nx) % ny, a / (nx * ny));
            // The two coordinates *across* the gap, which is what a patch is measured in.
            let across = match axis {
                0 => (j, k),
                1 => (i, k),
                _ => (i, j),
            };
            sheets.entry((axis, cells)).or_default().push(across);
        }

        let mut out = Vec::new();
        for ((_, cells), mut face) in sheets {
            face.sort_unstable();
            let mut seen = vec![false; face.len()];
            let index: BTreeMap<(usize, usize), usize> =
                face.iter().enumerate().map(|(n, p)| (*p, n)).collect();
            for start in 0..face.len() {
                if seen[start] {
                    continue;
                }
                // Four-connected flood fill over the facing cells, so a patch is a sheet somebody
                // could point at rather than every pair that happens to share a width.
                let mut stack = vec![start];
                seen[start] = true;
                let mut members = Vec::new();
                while let Some(n) = stack.pop() {
                    let (u, v) = face[n];
                    members.push((u, v));
                    for (du, dv) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
                        let (nu, nv) = (u as i64 + du, v as i64 + dv);
                        if nu < 0 || nv < 0 {
                            continue;
                        }
                        if let Some(&m) = index.get(&(nu as usize, nv as usize)) {
                            if !seen[m] {
                                seen[m] = true;
                                stack.push(m);
                            }
                        }
                    }
                }
                let (u0, u1) = (
                    members.iter().map(|m| m.0).min().unwrap_or(0),
                    members.iter().map(|m| m.0).max().unwrap_or(0),
                );
                let (v0, v1) = (
                    members.iter().map(|m| m.1).min().unwrap_or(0),
                    members.iter().map(|m| m.1).max().unwrap_or(0),
                );
                let (wu, wv) = (u1 + 1 - u0, v1 + 1 - v0);
                out.push(GapPatch {
                    pairs: members.len(),
                    span: (wu as f64 * dx, wv as f64 * dx),
                    distance: cells as f64 * dx,
                    rectangular: members.len() == wu * wv,
                });
            }
        }
        out
    }

    /// Mark cells as **nothing** — not a substance, not part of the block.
    ///
    /// A grid had no void until now, so the cells a part did not occupy were some other material,
    /// and an assembly of two parts in air was two parts buried in whatever the block was made
    /// of. `ARCHITECTURE.md` names it: "insulating it is a substance with a low conductivity,
    /// which is not the same thing". A low conductivity still conducts, still stores heat, and
    /// still sets a stability limit; nothing does none of those.
    ///
    /// A void cell holds no heat, conducts to nothing — every face it touches carries zero, which
    /// the harmonic mean already gives for a zero conductivity — takes no share of what arrives
    /// on the bus, and is left out of every average. Its temperature is **not a number**, because
    /// there is nothing there to have one, and `temperature_at` says so rather than returning a
    /// zero somebody would plot.
    ///
    /// # What crosses a gap, and what does not
    ///
    /// **Radiation does.** Two solid cells facing each other along a grid line across a run of
    /// void exchange `σA(T₁⁴ − T₂⁴)/(1/ε₁ + 1/ε₂ − 1)`, the parallel-plate series, with an
    /// exchange factor of **one**.
    ///
    /// That factor was written down here as a known approximation — a wide gap has a view factor
    /// well under one, so charging it as one was said to couple a wide gap too hard. Measuring it
    /// says otherwise, and the correction is worth more than the caveat was. The sides of a gap in
    /// this model are the **block's own outer faces**, and an insulated boundary is implemented
    /// as a mirror; a mirror puts an image of each surface beyond it and the images tile the
    /// plane, so the pair *is* two infinite parallel plates and one is exact at every width. What
    /// the old note described was a different geometry from the one the model has.
    ///
    /// The geometry it described is a real one, though — two parts floating in vacuum, open to
    /// space, where most of what leaves one surface does miss the other. Which of the two a scene
    /// means is a statement about its boundary that a grid cannot infer, so it is not chosen here.
    /// It is *reported*: [`gap_patches`](Solid3D::gap_patches) groups a clearance into the sheets
    /// a view factor is a statement about and [`GapPatch::view_factor`] gives the open-gap number,
    /// so the difference between the two readings is a factor somebody can see. For a 32 mm square
    /// part 16 mm under a lid it is **2.4x**.
    ///
    /// Still not here: a radiative-exchange solver. Side walls that are real material at a real
    /// temperature do not radiate into the gap at all, and that is a domain rather than a boundary
    /// condition.
    ///
    /// **Convection does not.** A gap full of air carries heat by moving that air, which needs a
    /// Rayleigh number and a correlation, and a correlation is not a closed form. So a gap in
    /// air is coupled *less* here than it really is, by however much the convection would have
    /// carried — and for a millimetre-scale gap at modest temperatures that is the same order as
    /// the radiation, so the answer is a lower bound rather than an estimate.
    ///
    /// **Conduction does not**, and that is right: there is nothing there to conduct through.
    ///
    /// It is still not a fluid. A part in air also loses heat to the room, and that is
    /// [`Solid3D::losing_from`], which is about the block's **outer** faces.
    pub fn empty(mut self, which: impl Fn(usize, usize, usize) -> bool) -> Solid3D {
        let (nx, ny, nz) = self.counts;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if which(i, j, k) {
                        self.void[i + nx * (j + ny * k)] = true;
                    }
                }
            }
        }
        self.resolve();
        self
    }

    /// Whether this cell is void — nothing rather than a substance.
    pub fn is_void(&self, i: usize, j: usize, k: usize) -> bool {
        self.index(i, j, k).map(|c| self.void[c]).unwrap_or(false)
    }

    /// How many cells hold nothing.
    pub fn void_cells(&self) -> usize {
        self.void.iter().filter(|v| **v).count()
    }

    /// Expose a face to an environment, so the block can lose heat through it.
    ///
    /// **Until this existed a `Solid3D` was adiabatic on all six faces**, which means no
    /// three-dimensional thermal scene could reach a steady state: every one of them warmed
    /// for as long as it ran. That is honest for a pulse and useless for the question a
    /// designer actually asks — *what temperature does this run at* — which is the question a
    /// chip package, a magnet busbar, a motor and a factory cell are all asking.
    ///
    /// The loss is [`Environment::loss_from`]'s, convective **and** radiative, applied per
    /// boundary cell at that cell's own temperature and its own emissivity. Per cell rather
    /// than to a mean, because that is what makes a gradient: the middle of a face runs hotter
    /// than its edge, and a lumped loss says it does not.
    ///
    /// The environment's `area` is the **whole face's** area; each boundary cell is charged its
    /// share, `area / cells_on_the_face`. Stating the face's area rather than a cell's keeps
    /// the number a caller writes independent of the grid they chose, which is what lets the
    /// same scene be refined without becoming a different problem — the property
    /// `dualis-world`'s `verify` sweep depends on.
    ///
    /// Exposing a face **tightens the stability limit**, because a cell that can also lose to
    /// air has more conductance leaving it. [`Solid3D::max_stable_dt`] accounts for it; a
    /// caller stepping by hand past the returned limit is refused as ever.
    pub fn losing_from(mut self, face: Face, environment: Environment) -> Solid3D {
        self.exposed.insert(face, environment);
        self.resolve();
        self
    }

    /// Heat given up to the exposed faces' environments over the run.
    pub fn lost_energy(&self) -> Energy {
        Energy::from_si(self.lost)
    }

    /// Which faces are exposed, and to what.
    pub fn exposed_faces(&self) -> impl Iterator<Item = (Face, &Environment)> + '_ {
        self.exposed.iter().map(|(f, e)| (*f, e))
    }

    /// The stability rate at the block's **present** state: conduction, plus each exposed
    /// cell's loss conductance at that cell's own temperature.
    ///
    /// State-dependent on purpose, which is what [`Domain::max_stable_dt`] means by "the
    /// largest step this domain can take **from `now`**". A block cooling from 1000 °C is
    /// stiffer at the start than at the end, and a limit that ignored that would be correct
    /// only at the end.
    ///
    /// Costs a pass over the boundary cells and nothing at all for an unexposed block, which is
    /// every scene that existed before this.
    fn worst_rate_now(&self) -> f64 {
        // Gaps as well as films: a pair of plates facing each other across void may have no
        // conducting face at all, so the radiative exchange is the *only* thing setting their
        // step. Returning early on `exposed` alone handed such a pair an infinite limit, and a
        // march at infinity is a NaN block reported as a substance with no diffusivity.
        if (self.exposed.is_empty() && self.gaps.is_empty()) || self.worst_rate.is_nan() {
            return self.worst_rate;
        }
        let (nx, ny, nz) = self.counts;
        let dx = self.dx.to_si();
        let mut worst = self.worst_rate;
        // A gap's radiative conductance, linearised at the **hotter** of the pair, which is the
        // conservative end: `dq/dT = 4σA T³/(1/ε₁+1/ε₂−1)` grows with temperature, so the hotter
        // cell's tangent bounds the pair's. Charged to both cells, because either could be the
        // one the step is too long for.
        for (a, b, coefficient) in &self.gaps {
            let hotter = self.cells[*a].max(self.cells[*b]);
            let g = 4.0 * coefficient * hotter.powi(3);
            for c in [*a, *b] {
                if self.capacity[c] > 0.0 {
                    worst = worst.max(g / self.capacity[c]);
                }
            }
        }
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let c = i + nx * (j + ny * k);
                    let loss = self.loss_conductance_at(self.temperature_at(i, j, k), (i, j, k));
                    if loss == 0.0 {
                        continue;
                    }
                    worst = worst.max(self.mobility[c] * (self.face_sum[c] + loss / dx));
                }
            }
        }
        worst
    }

    /// The film's flux out of a cell, in the sweep's own units — a conductivity times a
    /// temperature difference, so that `mobility · this` is a rate in kelvin per second.
    ///
    /// The **secant** conductance of the bare surface, which is the exact one for this flux,
    /// put in series with the half cell of solid behind it. Positive means heat leaving.
    fn film_flux(&self, old: &[f64], cell: (usize, usize, usize), c: usize) -> f64 {
        let thermal = self.substance_at(cell.0, cell.1, cell.2).thermal;
        let emissivity = thermal.map_or(0.0, |t| t.emissivity);
        let conductivity = thermal.map_or(f64::INFINITY, |t| t.conductivity.to_si());
        let dx = self.dx.to_si();
        let here = Temperature::from_si(old[c]);
        let mut out = 0.0;
        for (face, env) in &self.exposed {
            if !face.holds(cell, self.counts) {
                continue;
            }
            let share = env.area.to_si() / self.cells_on(*face).max(1) as f64;
            let gap = here.to_si() - env.ambient.to_si();
            if gap == 0.0 {
                continue;
            }
            let bare = Environment {
                ambient: env.ambient,
                convection_w_per_m2_k: env.convection_w_per_m2_k,
                area: Area::from_si(share),
            }
            .loss_from(here, emissivity)
            .to_si();
            let g = series_with_half_cell(bare / gap, conductivity, share, dx);
            out += g * gap / dx;
        }
        out
    }

    /// How many cells lie on a face.
    fn cells_on(&self, face: Face) -> usize {
        let (nx, ny, nz) = self.counts;
        match face {
            Face::XMin | Face::XMax => ny * nz,
            Face::YMin | Face::YMax => nx * nz,
            Face::ZMin | Face::ZMax => nx * ny,
        }
    }

    /// The loss conductance an exposed cell carries, in W/K, summed over the faces it lies on.
    ///
    /// Evaluated at the temperature handed in, and the callers hand in **the cell's own current
    /// temperature** rather than the ambient. That is not a refinement, it is the difference
    /// between a limit that holds and one that does not: the radiative conductance is `4εσT³`
    /// and a part at 1000 °C carries **26 times** what the same surface carries at room
    /// temperature — measured, `ε = 0.9`, 136 against 5.14 W/m²·K. A limit linearised about
    /// ambient would hand a hot block a step an order too long, and an explicit boundary past
    /// its limit does not diverge loudly: it oscillates about the air it is losing to while the
    /// conservation audit stays perfectly happy, because the heat really did leave.
    ///
    /// Two things make it conservative rather than merely plausible, and each was measured
    /// wrong first.
    ///
    /// The **derivative** `4εσT³` rather than the secant `q/(T−T∞)`, because `T⁴` is convex and
    /// the derivative is the larger of the two — four times it in the hot limit.
    ///
    /// And at `max(T_cell, T_ambient)`, because the derivative only dominates the secant for
    /// `T ≥ T∞`. Below ambient it reverses without bound, and the first version of this used
    /// the cell alone: a 20 mm cell at 20 °C inside a 700 °C radiant enclosure was handed a
    /// step whose true ratio was **13**, and one accepted step took it to 8847 °C. `4T∞³`
    /// bounds `(T+T∞)(T²+T∞²)` for every `T ≤ T∞`, so the larger of the two points is safe on
    /// both sides.
    ///
    /// The film is put **in series with the half cell of solid between the cell centre and the
    /// surface**, `2k·A/dx`. A finite-volume cell knows only its centre and the film acts at the
    /// surface, so the half cell between them is part of the path; charging the whole film
    /// against the centre sheds too much, by the cell Biot number `h·dx/2k`.
    ///
    /// **And it took two corrections to reach second order, not one.** The series form fixes
    /// the *space*: without it a review measured 2.08, 2.04, 2.02 per grid doubling against the
    /// interior's 4.0 — the signature this workspace has twice caught in an acoustic wall. With
    /// it the boundary was still first order, at 1.27, 1.71, 1.87, because the film was applied
    /// as a pass *after* the conduction sweep. That split's error carries a coefficient growing
    /// as `1/dx` while the step falls as `dx²`, so their product falls as `dx`. Applying the
    /// film inside the same explicit update — see `film_flux`, called from the sweep — makes
    /// them one operator, and `tests/the_cooled_boundary_order.rs` measures four per doubling at
    /// `Bi = 1.7` and at `Bi = 17`.
    fn loss_conductance_at(&self, at: Temperature, cell: (usize, usize, usize)) -> f64 {
        self.exposed
            .iter()
            .filter(|(face, _)| face.holds(cell, self.counts))
            .map(|(face, env)| {
                let share = env.area.to_si() / self.cells_on(*face).max(1) as f64;
                let thermal = self.substance_at(cell.0, cell.1, cell.2).thermal;
                let emissivity = thermal.map_or(0.0, |t| t.emissivity);
                let hot = at.to_si().max(env.ambient.to_si());
                let h = env.convection_w_per_m2_k
                    + 4.0 * emissivity * STEFAN_BOLTZMANN.to_si() * hot.powi(3);
                series_with_half_cell(
                    h * share,
                    thermal.map_or(f64::INFINITY, |t| t.conductivity.to_si()),
                    share,
                    self.dx.to_si(),
                )
            })
            .sum()
    }

    /// How many cells along each axis.
    pub fn counts(&self) -> (usize, usize, usize) {
        self.counts
    }

    /// The cell side.
    pub fn spacing(&self) -> Length {
        self.dx
    }

    /// The block's extent, which is `counts × dx` — the outer faces, not the cell centres.
    pub fn size(&self) -> LengthVec {
        let (nx, ny, nz) = self.counts;
        LengthVec::from_si(DVec3::new(nx as f64, ny as f64, nz as f64) * self.dx.to_si())
    }

    /// The flat index of a cell, or `None` if any component is out of range.
    ///
    /// Returned rather than panicking because the natural way to write a stencil is to ask for a
    /// neighbour that may not exist, and a boundary is exactly where that happens.
    pub fn index(&self, i: usize, j: usize, k: usize) -> Option<usize> {
        let (nx, ny, nz) = self.counts;
        (i < nx && j < ny && k < nz).then(|| i + nx * (j + ny * k))
    }

    /// Where the centre of cell `(i, j, k)` is, in the block's own coordinates.
    pub fn centre_of(&self, i: usize, j: usize, k: usize) -> LengthVec {
        LengthVec::from_si(
            DVec3::new(i as f64 + 0.5, j as f64 + 0.5, k as f64 + 0.5) * self.dx.to_si(),
        )
    }

    /// The temperature of one cell. Out of range reads the nearest one in range.
    pub fn temperature_at(&self, i: usize, j: usize, k: usize) -> Temperature {
        let (nx, ny, nz) = self.counts;
        let idx = self
            .index(i.min(nx - 1), j.min(ny - 1), k.min(nz - 1))
            .expect("clamped indices are in range");
        // **Not a number where there is nothing.** A void cell has no temperature, and a zero
        // or an ambient here is a value somebody would plot, average or believe. Every reader in
        // this workspace that draws a field already skips a non-finite sample.
        if self.void[idx] {
            return Temperature::from_si(f64::NAN);
        }
        Temperature::from_si(self.cells[idx])
    }

    /// Set one cell, for an initial condition a constructor cannot express.
    ///
    /// **This does not change what the block has absorbed.** It is a statement about the initial
    /// state, not a delivery of heat, so `stored_heat` moves and `absorbed_energy` does not — and
    /// a simulation started this way and then audited will show the difference as its opening
    /// balance rather than as a leak. Use [`deposit`](Solid3D::deposit) for heat that arrived.
    ///
    /// Out of range is ignored rather than a panic: a caller writing a hot spot in a loop over a
    /// radius is the expected use, and clipping at the edge is what they mean.
    pub fn set_temperature(&mut self, i: usize, j: usize, k: usize, t: Temperature) {
        if let Some(idx) = self.index(i, j, k) {
            self.cells[idx] = t.to_si();
            // A temperature is not a state for a substance that melts — 0 °C is ice, water, or any
            // mixture of the two. Naming a temperature therefore names a phase as well, and the
            // choice is the unmelted one at the point itself: a caller clamping a face below
            // freezing means ice, and a caller clamping it above means liquid.
            if self.latent[idx] > 0.0 {
                self.melted[idx] = if t.to_si() > self.melt_point[idx] {
                    1.0
                } else {
                    0.0
                };
            }
        }
    }

    /// Put joules into one cell, as heat that arrived there.
    ///
    /// Counts toward [`absorbed_energy`](Solid3D::absorbed_energy), so the books balance. Out of
    /// range is ignored, which would silently lose energy — so it does not: the joules are
    /// refused, and nothing is added to either total.
    pub fn deposit(&mut self, i: usize, j: usize, k: usize, joules: Energy) {
        let Some(idx) = self.index(i, j, k) else {
            return;
        };
        let rise = joules.to_si() / self.capacity[idx];
        self.add_kelvin(idx, rise);
        self.absorbed += joules.to_si();
    }

    /// Mean over every cell. Every cell has the same volume, so this is the volume average.
    pub fn mean_temperature(&self) -> Temperature {
        // Over what is there. Averaging in the void's placeholder would drag the number toward
        // a value nothing holds, and the more of the box is empty the further it would drag.
        let (sum, n) = (0..self.cells.len())
            .filter(|c| !self.void[*c])
            .fold((0.0, 0usize), |(s, n), c| (s + self.cells[c], n + 1));
        Temperature::from_si(if n == 0 { f64::NAN } else { sum / n as f64 })
    }

    /// The hottest cell — the number a hot spot exists to produce, and the one a lumped model
    /// reports as the mean.
    pub fn peak_temperature(&self) -> Temperature {
        Temperature::from_si(
            (0..self.cells.len())
                .filter(|c| !self.void[*c])
                .map(|c| self.cells[c])
                .fold(f64::MIN, f64::max),
        )
    }

    /// The coldest cell.
    pub fn coldest_temperature(&self) -> Temperature {
        Temperature::from_si(
            (0..self.cells.len())
                .filter(|c| !self.void[*c])
                .map(|c| self.cells[c])
                .fold(f64::MAX, f64::min),
        )
    }

    /// Heat taken from the bus over the run.
    pub fn absorbed_energy(&self) -> Energy {
        Energy::from_si(self.absorbed)
    }

    /// `α·dt/dx²` for the constructor's substance — the classic Fourier number.
    ///
    /// This is the quantity [`mode_amplification`](Solid3D::mode_amplification) is written in terms
    /// of, so it stays the textbook one and does not become shape- or fill-aware. For a block of
    /// one material with at least two cells on every axis it *is* the stability number, and
    /// `fourier_number(max_stable_dt) == 1/6` exactly.
    ///
    /// It is **not** the stability number for a thin block or a filled one, and those are the two
    /// cases where a caller sizing a step by hand needs the other one — see
    /// [`stability_ratio`](Solid3D::stability_ratio), which is what `step` actually enforces.
    pub fn fourier_number(&self, dt: Time) -> f64 {
        let Some(alpha) = self.materials[0].diffusivity() else {
            return f64::INFINITY;
        };
        alpha.to_si() * dt.to_si() / (self.dx.to_si() * self.dx.to_si())
    }

    /// `dt` as a fraction of the largest step this block is stable at — exactly one at the limit.
    ///
    /// The number `step` refuses on, and the honest one for a block that is thin or filled. It is
    /// `dt · maxᵢ (Σ_f k_f / Cᵢ) · dx`: a maximum over **cells**, because stability is a statement
    /// about a row of the update matrix and each cell has its own.
    pub fn stability_ratio(&self, dt: Time) -> f64 {
        dt.to_si() * self.worst_rate_now()
    }

    /// The heat capacity of the whole block, which for a filled one is a sum and not a product.
    pub fn heat_capacity(&self) -> HeatCapacity {
        HeatCapacity::from_si(self.capacity.iter().sum())
    }

    /// The exact per-step amplification of one separable cosine mode.
    ///
    /// `(a, b, c)` are half-wave counts along the three axes: `(1, 0, 0)` is the longest mode
    /// along x with the other two flat. Returns the factor the mode's amplitude is multiplied by
    /// in one step of `dt`, which is `1 + F·Σ(−4 sin²(mπ/2n))` — **exact**, not a linearisation,
    /// because that mode is an eigenvector of the discrete operator this domain steps with.
    ///
    /// Public because it is what makes this domain checkable without a reference implementation:
    /// a caller can predict an amplitude arbitrarily far ahead and compare. It is also what a
    /// grid designer wants, since a factor outside `(−1, 1]` is the instability itself.
    pub fn mode_amplification(&self, mode: (usize, usize, usize), dt: Time) -> f64 {
        let f = self.fourier_number(dt);
        1.0 + f * self.mode_eigenvalue_dx2(mode)
    }

    /// The discrete Laplacian eigenvalue of a mode, times `dx²` — dimensionless, in `[-12, 0]`.
    fn mode_eigenvalue_dx2(&self, mode: (usize, usize, usize)) -> f64 {
        let (nx, ny, nz) = self.counts;
        let term = |m: usize, n: usize| {
            if n <= 1 {
                // One cell across an axis is a mirror against itself: no gradient is
                // representable, so that direction contributes nothing at all.
                return 0.0;
            }
            let s = (m as f64 * std::f64::consts::PI / (2.0 * n as f64)).sin();
            -4.0 * s * s
        };
        term(mode.0, nx) + term(mode.1, ny) + term(mode.2, nz)
    }

    /// Fill the block with one separable cosine mode about a mean.
    ///
    /// The initial condition the closed-form tests use, and a genuinely useful one for anybody
    /// checking a grid: it is the only shape whose future this domain can state exactly.
    pub fn release_mode(&mut self, mode: (usize, usize, usize), mean: Temperature, amplitude: f64) {
        let (nx, ny, nz) = self.counts;
        let phase = |m: usize, n: usize, i: usize| {
            if n <= 1 {
                1.0
            } else {
                (m as f64 * std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos()
            }
        };
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = i + nx * (j + ny * k);
                    self.cells[idx] = mean.to_si()
                        + amplitude
                            * phase(mode.0, nx, i)
                            * phase(mode.1, ny, j)
                            * phase(mode.2, nz, k);
                }
            }
        }
        self.saved.clone_from(&self.cells);
    }

    /// The amplitude of one mode currently present, by projection.
    ///
    /// The counterpart to [`release_mode`](Solid3D::release_mode), and what makes a decay
    /// measurable rather than merely visible. Cosine modes on this grid are orthogonal, so this
    /// is exact for a block holding one of them and is the correct coefficient for a block
    /// holding several.
    pub fn mode_amplitude(&self, mode: (usize, usize, usize)) -> f64 {
        let (nx, ny, nz) = self.counts;
        let phase = |m: usize, n: usize, i: usize| {
            if n <= 1 {
                1.0
            } else {
                (m as f64 * std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos()
            }
        };
        // ⟨cos²⟩ is 1/2 per axis that actually varies, and 1 for an axis that cannot.
        let norm = [(mode.0, nx), (mode.1, ny), (mode.2, nz)]
            .iter()
            .map(|&(m, n)| if n <= 1 || m == 0 { 1.0 } else { 0.5 })
            .product::<f64>();
        let mut sum = 0.0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = i + nx * (j + ny * k);
                    sum += self.cells[idx]
                        * phase(mode.0, nx, i)
                        * phase(mode.1, ny, j)
                        * phase(mode.2, nz, k);
                }
            }
        }
        sum / (self.cells.len() as f64 * norm)
    }

    /// The volume of the whole block.
    pub fn volume(&self) -> Volume {
        let (nx, ny, nz) = self.counts;
        Volume::from_si((nx * ny * nz) as f64 * self.cell_volume())
    }

    fn cell_volume(&self) -> f64 {
        let dx = self.dx.to_si();
        dx * dx * dx
    }

    /// Heat held, measured from the state the block started in.
    ///
    /// Weighted cell by cell, so a filled block's books balance for the same reason a uniform
    /// one's do: the sweep moves `G_f·ΔT` across a face and takes it off one side and puts it on
    /// the other, and this is the sum that is therefore constant.
    ///
    /// **In enthalpy, not in temperature**, so a melting front is on the books. A cell holding at its
    /// melting point while it absorbs 306 mJ per cubic millimetre has taken that heat in and its
    /// temperature says nothing about it; an audit reading temperature alone would call it a leak.
    fn stored_heat(&self) -> f64 {
        (0..self.cells.len())
            .map(|c| self.enthalpy_joules(c) - self.reference_enthalpy_joules(c))
            .sum()
    }

    /// The cell's state as one number, in kelvin above its melting point.
    ///
    /// Monotone and invertible, which is the whole reason the phase change needs no iteration and
    /// conserves exactly:
    ///
    /// ```text
    ///   solid    e = T − T_m ≤ 0
    ///   mush     e = φ·ℓ ∈ [0, ℓ]      at T = T_m
    ///   liquid   e = ℓ + T − T_m ≥ ℓ
    /// ```
    ///
    /// This is the enthalpy method, bookkept as a temperature and a fraction rather than as a single
    /// enthalpy field — identical arithmetic, and it keeps `cells` holding kelvin so everything that
    /// reads a temperature still can.
    ///
    /// It is **not** the apparent-heat-capacity method, and the difference is the failure that method
    /// has: smearing `L` over a narrow temperature interval lets a cell cross the whole interval in
    /// one step and skip the latent heat, which runs the front fast and conserves nothing. Here a
    /// step that overshoots the mush deposits the remainder as sensible heat on the far side, because
    /// the inverse map says where the energy goes rather than a branch guessing.
    fn enthalpy_joules(&self, c: usize) -> f64 {
        let (t, point, latent) = (self.cells[c], self.melt_point[c], self.latent[c]);
        if latent <= 0.0 {
            // No phase change. Measured from the block's own reference rather than from a melting
            // point at infinity, and that choice is precision: subtracting a nearby number keeps
            // digits that subtracting 273.15 from 293.15 does not.
            return self.capacity[c] * (t - self.reference);
        }
        if self.melted[c] >= 1.0 {
            latent + self.cap_liquid[c] * (t - point)
        } else if self.melted[c] <= 0.0 {
            self.cap_solid[c] * (t - point)
        } else {
            self.melted[c] * latent
        }
    }

    /// What [`enthalpy_kelvin`](Solid3D::enthalpy_kelvin) was when the block started, so the ledger
    /// reports what arrived rather than what is there.
    fn reference_enthalpy_joules(&self, c: usize) -> f64 {
        let (point, latent) = (self.melt_point[c], self.latent[c]);
        if latent <= 0.0 {
            return 0.0;
        }
        if self.reference > point {
            latent + self.cap_liquid[c] * (self.reference - point)
        } else {
            self.cap_solid[c] * (self.reference - point)
        }
    }

    /// Add `rise` kelvin of *enthalpy* to a cell, which is not the same as adding kelvin to it.
    ///
    /// The one place heat becomes state, so that a phase change cannot be forgotten at one of the
    /// several doors heat comes in through — the sweep, [`deposit`](Solid3D::deposit), and the plain
    /// channel all pass through here. A cell at its melting point takes the whole of it as melting
    /// and does not warm at all.
    fn add_kelvin(&mut self, c: usize, rise: f64) {
        if self.latent[c] > 0.0 {
            // `rise` is what the cell *would* have warmed by; the joules it stands for are that times
            // whichever capacity the cell has now, which is the mixed one and is exactly right for a
            // cell that is wholly one phase. A mushy cell does not change temperature, so the flux
            // that reaches it has to be converted through the capacity it had when the flux was
            // computed — the same one.
            let joules = rise * self.capacity[c];
            let e = self.enthalpy_joules(c) + joules;
            let (t, phi) = self.state_from_enthalpy(c, e);
            self.cells[c] = t;
            self.melted[c] = phi;
        } else {
            self.cells[c] += rise;
        }
    }

    /// The inverse: a state in kelvin back to a temperature and a melted fraction.
    fn state_from_enthalpy(&self, c: usize, e: f64) -> (f64, f64) {
        let (point, latent) = (self.melt_point[c], self.latent[c]);
        if e <= 0.0 {
            (point + e / self.cap_solid[c], 0.0)
        } else if e >= latent {
            (point + (e - latent) / self.cap_liquid[c], 1.0)
        } else {
            (point, e / latent)
        }
    }
}

impl Domain for Solid3D {
    fn books_balance(&self) -> bool {
        true
    }

    fn name(&self) -> &str {
        &self.name
    }

    /// `minᵢ Cᵢ / (dx · Σ_f k_f)` — the tightest row of the update matrix, and `dx²/(6α)` when the
    /// block is one material with at least two cells on every axis.
    ///
    /// # Why a maximum over cells, and what summing the faces buys
    ///
    /// Stability is Gershgorin's condition on one row: the update is `I + dt·D⁻¹L`, row `i` has
    /// diagonal `−θᵢ` and off-diagonals summing to `+θᵢ` for `θᵢ = dt Σ_f G_f / Cᵢ`, so every
    /// eigenvalue is in `[1 − 2θ_max, 1]` and `θ_max ≤ 1` is the limit.
    ///
    /// On a uniform grid that disc is **tight** rather than cautious — the sharpest representable
    /// mode saturates it — and this is where `1/2`, `1/4` and `1/6` come from. They are not three
    /// dimensions; they are however many axes have more than one cell, and the block is stepped at
    /// whichever applies. This used to report `dx²/(6α)` for every shape, which cost a bar-shaped
    /// block three times the steps for nothing.
    ///
    /// On a filled grid the value is the opposite one: the limit is usually far **looser** than
    /// `dx²/(6·α_max)`, because `k_f ≤ 2·min(k_L, k_R)` means a cell cannot be heated through a
    /// face faster than its worse side allows. Heat does not reach a fast material at that
    /// material's own rate; it reaches it at the rate the neighbour delivers. Measured on one
    /// aluminium cell embedded in borosilicate, the honest limit is **75×** the one aluminium's
    /// diffusivity would name, and that factor is wall-clock.
    ///
    /// It can in principle go the other way — up to a factor of two, when a cell's neighbours
    /// conduct better *and* store more — but that needs volumetric heat capacity to vary as widely
    /// as conductivity, and across solids it varies by one order of magnitude where conductivity
    /// varies by four. So the tightening is a bound this catalogue cannot reach, and the loosening
    /// is a saving any coating or inclusion gets.
    ///
    /// # A third of the bar's, for a cube
    ///
    /// A third of what [`Bar1D`](crate::Bar1D) reports for the same spacing and material. That
    /// factor is the reason `Schedule::Multirate` exists: a block and a lumped mass in one scene
    /// differ by five orders of magnitude in the step they can take, and a single global step
    /// would make the cheap domain pay the expensive one's bill.
    ///
    /// # Stable is not accurate, and this is only the first
    ///
    /// At **exactly** this step the sharpest mode the grid can hold has an amplification factor of
    /// `−1`. It flips sign every step and never decays. That is what marginal stability means and
    /// it is not a defect — the scheme does not diverge there, which is the whole of what a
    /// stability limit claims.
    ///
    /// It does mean sharp initial data is carried badly. A point source excites that mode as hard
    /// as anything can, and the peak comes out **1.96×** the exact answer while the conservation
    /// audit stays exact to the last bit. At half this step it is 1.005×.
    ///
    /// So take a fraction of it when the initial condition is sharp. `Schedule::Multirate` divides
    /// by `ceil(dt / limit)` and so usually lands comfortably inside, but a caller stepping by
    /// hand can sit exactly on it. `cargo run --example heat_in_three_dimensions` is the
    /// measurement.
    fn max_stable_dt(&self, _now: Time) -> Time {
        let rate = self.worst_rate_now();
        if rate.is_nan() {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(1.0 / rate)
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let ratio = self.stability_ratio(dt);
        if ratio.is_nan() {
            return Err(Violation::at(
                &self.name,
                "substance has no diffusivity",
                f64::INFINITY,
            ));
        }
        // Reported in Fourier-number units, which is `dt/(6·max_stable_dt)` — the classic
        // `α·dt/dx²` for a block of one material with two cells on every axis, and the number the
        // limit is expressed in for every other block.
        if ratio > 1.0 + 6e-12 {
            return Err(Violation {
                quantity: "Fourier number".to_string(),
                site: format!("{} (explicit 3D conduction)", self.name),
                before: STABLE_FOURIER_3D,
                after: ratio * STABLE_FOURIER_3D,
                scale: STABLE_FOURIER_3D,
                tolerance: 1e-12,
            });
        }

        // Heat off the plain channel, which carries an amount and no location.
        //
        // **Spread evenly**, and this is the one place where the 3D domain must not copy the 1D
        // one. `Bar1D` puts placeless heat in its first cell, and that is defensible for a bar,
        // which has an end that a surface absorbing light would plausibly be. A block has six
        // faces and no distinguished cell, so choosing one would invent a location the bus never
        // carried — and a hot spot that came from a tie-break is worse than no hot spot, because
        // it looks like physics. Even spreading is the unique choice that adds no information.
        // Heat that *does* have a place arrives through `deposit` or over an `Interface`.
        // Spread to a uniform **rise**, which for a filled block means in proportion to each cell's
        // capacity rather than in equal joules. Equal joules would warm the low-capacity material
        // more and so would say where the heat landed, which the bus never carried.
        let gained = bus.take_share(HEAT, dt);
        if gained != 0.0 {
            self.absorbed += gained;
            let rise = gained / self.capacity.iter().sum::<f64>();
            for c in 0..self.cells.len() {
                // Void takes no share: it has no capacity to contribute to the sum above and no
                // temperature to raise. Warming it would put joules where there is nothing to
                // hold them, and the ledger would then disagree with the block.
                if self.void[c] {
                    continue;
                }
                self.add_kelvin(c, rise);
            }
        }

        // The seven-point stencil in conductance form: `Cᵢ ΔTᵢ = dt Σ_f G_f (T_f − Tᵢ)`.
        //
        // It is the **face flux** that is computed, and each face is read once from each side with
        // opposite sign, so `Σ Cᵢ Tᵢ` is conserved to the last bit rather than to a tolerance —
        // and it stays that way when the capacities differ, which a stencil written as
        // `f·(Σ T_n − 6T)` cannot do because there is no per-cell `f` in it.
        let (nx, ny, nz) = self.counts;
        let old = self.cells.clone();
        let dts = dt.to_si();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let c = i + nx * (j + ny * k);
                    let t = old[c];
                    let (xr, yr) = ((nx + 1) * (j + ny * k), nx * (j + (ny + 1) * k));
                    let mut flux = 0.0;
                    if i > 0 {
                        flux += self.kx[i + xr] * (old[c - 1] - t);
                    }
                    if i + 1 < nx {
                        flux += self.kx[i + 1 + xr] * (old[c + 1] - t);
                    }
                    if j > 0 {
                        flux += self.ky[i + yr] * (old[c - nx] - t);
                    }
                    if j + 1 < ny {
                        flux += self.ky[i + nx + yr] * (old[c + nx] - t);
                    }
                    if k > 0 {
                        flux += self.kz[c] * (old[c - nx * ny] - t);
                    }
                    if k + 1 < nz {
                        flux += self.kz[c + nx * ny] * (old[c + nx * ny] - t);
                    }
                    debug_assert_eq!(t, self.cells[c], "the sweep reads `old` and writes `cells`");

                    // **The film is part of this flux, not a pass after it.** Applying it
                    // separately is Lie splitting, and the split's error carries a coefficient
                    // that grows as `1/dx` while the step falls as `dx²` — so the product falls
                    // as `dx`, and the boundary came out **first order** while the interior was
                    // second. Measured before this line moved here: ratios 1.27, 1.71, 1.87 per
                    // grid doubling, approaching two rather than four. In the same update they
                    // are one operator and the order is the interior's.
                    //
                    // Read off `old` like every other term, so the sweep stays a function of the
                    // state it began with.
                    let shed = if self.exposed.is_empty() {
                        0.0
                    } else {
                        self.film_flux(&old, (i, j, k), c)
                    };
                    self.add_kelvin(c, dts * self.mobility[c] * (flux - shed));
                    // What that removed, in joules, counted in the same statement that removes
                    // it so `stored + lost` stays exact.
                    if shed != 0.0 {
                        self.lost += dts * self.mobility[c] * shed * self.capacity[c];
                    }
                }
            }
        }

        // **What the gaps carry.** A pair of solid cells facing each other across void exchange
        // radiation, and it is applied here — in the same explicit pass as the conduction sweep,
        // read off the same `old` state — for the reason the surface film moved into the sweep:
        // a separate pass is Lie splitting, whose error carries a coefficient growing as `1/dx`
        // against a step falling as `dx²`, and the boundary comes out an order worse than the
        // interior.
        //
        // **Antisymmetric**, so the pair conserves exactly: what leaves one arrives in the
        // other, in the same statement, and `Σ Cᵢ Tᵢ` is unchanged to the last bit as it is for
        // a conduction face.
        if !self.gaps.is_empty() {
            let dts = dt.to_si();
            for (a, b, coefficient) in self.gaps.clone() {
                let (ta, tb) = (old[a], old[b]);
                let watts = coefficient * (ta.powi(4) - tb.powi(4));
                if watts == 0.0 {
                    continue;
                }
                let joules = watts * dts;
                self.add_kelvin(a, -joules / self.capacity[a]);
                self.add_kelvin(b, joules / self.capacity[b]);
            }
        }

        // **A two-phase block re-derives its conductivities, because they moved.**
        //
        // A cell's `k` and `c` depend on how much of it has melted, so a front that advanced changed
        // the operator. `resolve` carries the fractions the sweep just produced rather than rebuilding
        // them from temperature, and recomputes the faces, the capacities and the limit from them.
        //
        // Only when there is a liquid phase to mix toward: `latent > 0` alone is the one-phase model,
        // whose properties do not depend on the fraction at all, and paying for a rebuild there would
        // slow every non-melting block for nothing.
        //
        // Measured, so the cost is on the record rather than assumed: a `resolve` is 1.6x a `step` at
        // 40 cells and 4.2x at 4096, so a two-phase sweep runs 2.6x to 5.2x a one-phase one. That is
        // the price of the simple version. An incremental update touching only the mush — one or two
        // cells wide — is the optimisation available if a problem ever needs it, and none does yet.
        // The film is applied inside the sweep above, in the same explicit update as
        // conduction — see . It used to run here, as a pass afterwards, and that
        // split cost the boundary an order.
        if self.two_phase {
            self.resolve();
        }
        Ok(())
    }

    /// Heat gained since the start. The faces are insulated, so this is exactly what came in.
    fn ledger(&self) -> Ledger {
        // `stored + lost`, so the total moves only by what crossed the bus and the claim in
        // `books_balance` survives a block that sheds heat to air. `LumpedMass` carries the
        // same pair for the same reason; what differs is that this one keeps the claim,
        // because every joule it sheds is counted here in the same step it leaves a cell.
        Ledger::new().with(quantity::ENERGY, self.stored_heat() + self.lost)
    }

    /// The temperatures **and** the phase, because a temperature alone is not a state.
    ///
    /// A cell at 0 °C is ice, water or any mixture, so saving `cells` and not `melted` would restore
    /// a block that was half melted as one that was entirely solid at the same temperature — losing
    /// 306 mJ per cubic millimetre with nothing to say it had gone. `Schedule::Iterative` and the
    /// audit's retry both restore, so this is a live path and not a precaution.
    fn checkpoint(&mut self) {
        self.saved.clone_from(&self.cells);
        self.saved_melted.clone_from(&self.melted);
        // The losses too, and `LumpedMass` records what it costs to forget: an iterative sweep
        // that rewinds the cells and not the counter grows its books by one sweep of shed heat
        // per iteration, and the audit reports energy created from nothing.
        self.saved_lost = self.lost;
    }

    fn restore(&mut self) {
        self.cells.clone_from(&self.saved);
        self.melted.clone_from(&self.saved_melted);
        self.lost = self.saved_lost;
    }

    fn supports_restore(&self) -> bool {
        true
    }

    /// Peak, mean and coldest, in celsius, and what it has absorbed.
    ///
    /// All three ends of the distribution, because the whole reason to pay for a 3D grid is that
    /// they differ. A block reported by its mean alone is a `LumpedMass` that cost `n³` times as
    /// much, and the gap between peak and mean is the number that says whether the reduction
    /// would have been honest.
    fn readings(&self) -> Vec<Reading> {
        let mut out = vec![
            Reading::new(
                &self.name,
                "peak",
                self.peak_temperature().to_si() - 273.15,
                "C",
            ),
            Reading::new(
                &self.name,
                "mean",
                self.mean_temperature().to_si() - 273.15,
                "C",
            ),
            Reading::new(
                &self.name,
                "coldest",
                self.coldest_temperature().to_si() - 273.15,
                "C",
            ),
            Reading::new(&self.name, "absorbed", self.absorbed_energy().to_si(), "J"),
        ];
        // The melted volume, and **only** for a block that can melt. A column of zeros in every
        // report tells a reader nothing and costs a line in all of them; the condition is a property
        // of the block's materials fixed at construction, so it is not a mode that can surprise
        // anybody mid-run. Without it a phase change is the one thing this domain does that a report
        // could not see, and the temperature would say a cell was holding still.
        if self.latent.iter().any(|l| *l > 0.0) {
            out.push(Reading::new(
                &self.name,
                "melted",
                self.melted_volume().to_si() * 1e9,
                "mm3",
            ));
        }
        out
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }

    /// A temperature field, so nothing above has to know this is a block.
    fn as_field(&self) -> Option<&dyn ScalarField> {
        Some(self)
    }
}

impl ScalarField for Solid3D {
    /// **Kelvin**, because that is what the cells hold. See [`Bar1D`](crate::Bar1D).
    fn unit(&self) -> &'static str {
        "K"
    }

    /// Trilinear between cell centres, clamped at the faces, and **masked over void**.
    ///
    /// Clamped rather than extrapolated: outside an insulated face the temperature is not
    /// defined, and continuing the gradient would draw a block hotter than any cell in it.
    ///
    /// Masked because `self.cells` still holds whatever an emptied cell held when it was emptied,
    /// and reading it raw is how a gap came out of every exporter as a piece of the block sitting
    /// at ambient forever — a plausible, unchanging number, which is worse than no number. The
    /// weights are renormalised over the solid corners, so a sample at a cell centre is that cell
    /// exactly, a sample inside a clearance is `NaN`, and a sample straddling the two is the
    /// material's own value rather than a blend with something that is not there.
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let (nx, ny, nz) = self.counts;
        let q = p.to_si() / self.dx.to_si() - DVec3::splat(0.5);
        // NaN spelled out rather than folded into a comparison: a visualiser can hand one over,
        // and it must not reach the cast below. Answered with a `NaN` rather than with cell zero —
        // a question about nowhere has no answer, and cell zero may itself be empty, in which case
        // the old fallback returned a frozen number for a point that was never asked about.
        if q.is_nan() {
            return f64::NAN;
        }
        let axis = |v: f64, n: usize| -> (usize, f64) {
            let last = n.saturating_sub(1);
            if v <= 0.0 {
                return (0, 0.0);
            }
            if v >= last as f64 {
                return (last, 0.0);
            }
            let i = v.floor();
            (i as usize, v - i)
        };
        let (i, fx) = axis(q.x, nx);
        let (j, fy) = axis(q.y, ny);
        let (k, fz) = axis(q.z, nz);

        // **A sample belongs to the cell it is in**, and if that cell is empty there is nothing
        // there to sample. Without this the masked weights below still answer — with the value of
        // whichever solid neighbour the sample leans towards — and a panel that samples across the
        // extent rather than on cell centres lands most of a clearance's first layer inside the
        // material's half. Measured on a two-layer gap: 36 of the 72 empty cells came out solid.
        let nearest = |v: f64, n: usize| (v.round().max(0.0) as usize).min(n.saturating_sub(1));
        if self.void[nearest(q.x, nx) + nx * (nearest(q.y, ny) + ny * nearest(q.z, nz))] {
            return f64::NAN;
        }
        let (i1, j1, k1) = (
            (i + 1).min(nx - 1),
            (j + 1).min(ny - 1),
            (k + 1).min(nz - 1),
        );

        let mut sum = 0.0;
        let mut weight = 0.0;
        for (a, wa) in [(i, 1.0 - fx), (i1, fx)] {
            for (b, wb) in [(j, 1.0 - fy), (j1, fy)] {
                for (c, wc) in [(k, 1.0 - fz), (k1, fz)] {
                    let w = wa * wb * wc;
                    let at = a + nx * (b + ny * c);
                    // A corner with no weight is not a corner, so a sample sitting on a cell
                    // centre never consults the neighbour it does not use — which is what makes
                    // an exactly-sampled grid exact rather than nearly so.
                    if w > 0.0 && !self.void[at] {
                        sum += w * self.cells[at];
                        weight += w;
                    }
                }
            }
        }
        if weight > 0.0 {
            sum / weight
        } else {
            f64::NAN
        }
    }

    /// Central differences on the cell grid, mirrored at the faces **and at void**.
    ///
    /// A clearance is a boundary, and the block already knows what to do at one: mirror. Walking
    /// off the outer edge and walking into nothing are the same situation — there is no sample
    /// that way — so both return the asking cell and the difference becomes one-sided, which is
    /// the insulated condition rather than a slope towards a number that is not there.
    ///
    /// `NaN` inside a clearance, for the same reason [`temperature_at`](Solid3D::temperature_at)
    /// gives one: a gradient there is a value somebody would plot.
    fn gradient(&self, p: LengthVec, _t: Time, _h: Length) -> DVec3 {
        let (i, j, k) = self.nearest_cell(p);
        if self.absent(i, j, k) {
            return DVec3::splat(f64::NAN);
        }
        let d = 2.0 * self.dx.to_si();
        let m = |a: isize, b: isize, c: isize| self.mirrored_into(i, j, k, a, b, c);
        DVec3::new(
            (m(i + 1, j, k) - m(i - 1, j, k)) / d,
            (m(i, j + 1, k) - m(i, j - 1, k)) / d,
            (m(i, j, k + 1) - m(i, j, k - 1)) / d,
        )
    }

    /// `∇²T` on the seven-point stencil.
    ///
    /// The Laplacian of the *temperature*, which is what the trait asks for and is a statement
    /// about the field rather than about the material. It is the operator the sweep uses only when
    /// the block is one material; for a filled one the sweep uses `∇·(k∇T)`, and
    /// [`rate`](ScalarField::rate) is where that appears.
    fn laplacian(&self, p: LengthVec, _t: Time, _h: Length) -> f64 {
        let (i, j, k) = self.nearest_cell(p);
        if self.absent(i, j, k) {
            return f64::NAN;
        }
        let dx = self.dx.to_si();
        let centre = self.mirrored(i, j, k);
        // Each neighbour that is not there mirrors back to the centre and contributes nothing to
        // `sum - 6·centre`, which is exactly what an insulated face does — and exactly what a face
        // touching a clearance carries, since the harmonic face mean is already zero there.
        let m = |a: isize, b: isize, c: isize| self.mirrored_into(i, j, k, a, b, c);
        let sum = m(i - 1, j, k)
            + m(i + 1, j, k)
            + m(i, j - 1, k)
            + m(i, j + 1, k)
            + m(i, j, k - 1)
            + m(i, j, k + 1);
        (sum - 6.0 * centre) / (dx * dx)
    }

    /// `∂T/∂t = (1/Cᵢ)·Σ_f G_f (T_f − Tᵢ)`, which is `α∇²T` where the block is one material.
    ///
    /// The conductance form rather than `α·laplacian`, so that it is the sweep's own operator at a
    /// point in a filled block too — read off the same face conductivities, not reconstructed from
    /// a diffusivity the block may not have a single value of.
    ///
    /// Conduction only: heat arriving over the bus is a source this cannot see.
    fn rate(&self, p: LengthVec, _t: Time, _dt: Time) -> f64 {
        let (nx, ny, nz) = self.counts;
        let (i, j, k) = self.nearest_cell(p);
        let clamp = |v: isize, n: usize| v.clamp(0, n as isize - 1) as usize;
        let (i, j, k) = (clamp(i, nx), clamp(j, ny), clamp(k, nz));
        let c = i + nx * (j + ny * k);
        // Nothing does not warm at a rate, and the `is_finite` guard below would otherwise report
        // that it warms at **zero** — `mobility` is `1/C`, a cell with no capacity has none, and
        // `inf * 0.0` is the `NaN` that guard was written to catch coming from somewhere else.
        // A confident zero here is indistinguishable from a solid cell in equilibrium.
        if self.void[c] {
            return f64::NAN;
        }
        let t = self.cells[c];
        let (xr, yr) = ((nx + 1) * (j + ny * k), nx * (j + (ny + 1) * k));
        let mut flux = 0.0;
        if i > 0 {
            flux += self.kx[i + xr] * (self.cells[c - 1] - t);
        }
        if i + 1 < nx {
            flux += self.kx[i + 1 + xr] * (self.cells[c + 1] - t);
        }
        if j > 0 {
            flux += self.ky[i + yr] * (self.cells[c - nx] - t);
        }
        if j + 1 < ny {
            flux += self.ky[i + nx + yr] * (self.cells[c + nx] - t);
        }
        if k > 0 {
            flux += self.kz[c] * (self.cells[c - nx * ny] - t);
        }
        if k + 1 < nz {
            flux += self.kz[c + nx * ny] * (self.cells[c + nx * ny] - t);
        }
        let rate = self.mobility[c] * flux;
        if rate.is_finite() {
            rate
        } else {
            0.0
        }
    }
}

impl Solid3D {
    /// Whether the cell a stencil is asking about is not there — off the grid, or empty.
    ///
    /// The two cases are one case for a derivative: neither is a sample. Signed, so a stencil can
    /// ask about a neighbour it has already walked off the edge to reach.
    fn absent(&self, i: isize, j: isize, k: isize) -> bool {
        let (nx, ny, nz) = self.counts;
        let outside = |v: isize, n: usize| v < 0 || v >= n as isize;
        if outside(i, nx) || outside(j, ny) || outside(k, nz) {
            return true;
        }
        self.void[i as usize + nx * (j as usize + ny * k as usize)]
    }

    /// [`mirrored`](Solid3D::mirrored), but a neighbour that is **not there** mirrors back to the
    /// cell that asked rather than to the edge of the grid.
    ///
    /// The distinction only shows up inside a block: `mirrored` clamps, so a stencil at cell 3
    /// reaching into a clearance at cell 4 would be handed cell 4's frozen value. This hands back
    /// cell 3, which is the one-sided difference a boundary deserves.
    fn mirrored_into(&self, ci: isize, cj: isize, ck: isize, i: isize, j: isize, k: isize) -> f64 {
        if self.absent(i, j, k) {
            self.mirrored(ci, cj, ck)
        } else {
            self.mirrored(i, j, k)
        }
    }

    /// The value at a cell, with out-of-range indices mirrored back — an insulated face.
    ///
    /// A mirror rather than a zero, and the distinction is the whole boundary condition: a zero
    /// neighbour is a face held at absolute zero and would drain the block, where a mirror is a
    /// face with no gradient across it and so no flow through it.
    fn mirrored(&self, i: isize, j: isize, k: isize) -> f64 {
        let (nx, ny, nz) = self.counts;
        let clamp = |v: isize, n: usize| v.clamp(0, n as isize - 1) as usize;
        let idx = clamp(i, nx) + nx * (clamp(j, ny) + ny * clamp(k, nz));
        self.cells[idx]
    }

    /// The cell a point falls in, as signed indices so a stencil can walk off the edge.
    fn nearest_cell(&self, p: LengthVec) -> (isize, isize, isize) {
        let q = p.to_si() / self.dx.to_si();
        let one = |v: f64| {
            if v.is_nan() {
                0
            } else {
                v.floor().clamp(-1.0, 1e9) as isize
            }
        };
        (one(q.x), one(q.y), one(q.z))
    }
}
