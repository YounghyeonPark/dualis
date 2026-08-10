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
    units::{Energy, Length, LengthVec, Temperature, Time, Volume},
    Domain, Exchange, Ledger, ScalarField, Substance, Violation,
};
use glam::DVec3;

use crate::HEAT;

/// The largest Fourier number an explicit three-dimensional sweep is stable at.
///
/// `1/(2d)` for `d` dimensions, from requiring the amplification factor of the worst-resolved
/// mode to stay inside the unit circle. Public because a caller sizing a grid needs it before
/// there is anything to ask.
pub const STABLE_FOURIER_3D: f64 = 1.0 / 6.0;

/// A rectangular block of one material, conducting in three dimensions.
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
#[derive(Clone, Debug)]
pub struct Solid3D {
    name: String,
    substance: Substance,
    /// Cell-centre temperatures, indexed `x + nx*(y + ny*z)`.
    cells: Vec<f64>,
    saved: Vec<f64>,
    counts: (usize, usize, usize),
    dx: Length,
    absorbed: f64,
    /// What [`stored_heat`](Solid3D::stored_heat) is measured from — see [`Bar1D`](crate::Bar1D)
    /// for why an enthalpy reference is chosen for precision rather than for physics.
    reference: f64,
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
        Solid3D {
            name: name.into(),
            substance,
            saved: cells.clone(),
            cells,
            counts,
            dx,
            absorbed: 0.0,
            reference: initial.to_si(),
        }
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
        let capacity = self.cell_capacity();
        self.cells[idx] += joules.to_si() / capacity;
        self.absorbed += joules.to_si();
    }

    /// Mean over every cell. Every cell has the same volume, so this is the volume average.
    pub fn mean_temperature(&self) -> Temperature {
        Temperature::from_si(self.cells.iter().sum::<f64>() / self.cells.len() as f64)
    }

    /// The hottest cell — the number a hot spot exists to produce, and the one a lumped model
    /// reports as the mean.
    pub fn peak_temperature(&self) -> Temperature {
        Temperature::from_si(self.cells.iter().copied().fold(f64::MIN, f64::max))
    }

    /// The coldest cell.
    pub fn coldest_temperature(&self) -> Temperature {
        Temperature::from_si(self.cells.iter().copied().fold(f64::MAX, f64::min))
    }

    /// Heat taken from the bus over the run.
    pub fn absorbed_energy(&self) -> Energy {
        Energy::from_si(self.absorbed)
    }

    /// `α·dt/dx²`. Must stay at or under [`STABLE_FOURIER_3D`].
    pub fn fourier_number(&self, dt: Time) -> f64 {
        let Some(alpha) = self.substance.diffusivity() else {
            return f64::INFINITY;
        };
        alpha.to_si() * dt.to_si() / (self.dx.to_si() * self.dx.to_si())
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

    fn cell_capacity(&self) -> f64 {
        self.substance
            .heat_capacity(Volume::from_si(self.cell_volume()))
            .map(|c| c.to_si())
            .unwrap_or(f64::INFINITY)
    }

    /// Heat held, measured from the temperature the block started at.
    fn stored_heat(&self) -> f64 {
        self.cell_capacity() * self.cells.iter().map(|t| t - self.reference).sum::<f64>()
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
}

impl Domain for Solid3D {
    fn name(&self) -> &str {
        &self.name
    }

    /// `dx²/(6α)` — the explicit limit with all three axes counted.
    ///
    /// A third of what [`Bar1D`](crate::Bar1D) reports for the same spacing and material. That
    /// factor is the reason `Schedule::Multirate` exists: a block and a lumped mass in one scene
    /// differ by five orders of magnitude in the step they can take, and a single global step
    /// would make the cheap domain pay the expensive one's bill.
    fn max_stable_dt(&self, _now: Time) -> Time {
        let Some(alpha) = self.substance.diffusivity() else {
            return Time::from_si(f64::INFINITY);
        };
        let dx = self.dx.to_si();
        Time::from_si(STABLE_FOURIER_3D * dx * dx / alpha.to_si())
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let f = self.fourier_number(dt);
        if !f.is_finite() {
            return Err(Violation::at(&self.name, "substance has no diffusivity", f));
        }
        if f > STABLE_FOURIER_3D + 1e-12 {
            return Err(Violation {
                quantity: "Fourier number".to_string(),
                site: format!("{} (explicit 3D conduction)", self.name),
                before: STABLE_FOURIER_3D,
                after: f,
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
        let gained = bus.take_share(HEAT, dt);
        if gained != 0.0 {
            self.absorbed += gained;
            let per_cell = gained / (self.cells.len() as f64 * self.cell_capacity());
            for cell in self.cells.iter_mut() {
                *cell += per_cell;
            }
        }

        // The seven-point stencil, with insulated faces by mirroring.
        let (nx, ny, nz) = self.counts;
        let previous = self.cells.clone();
        let old = Solid3D {
            cells: previous,
            ..self.clone()
        };
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let (i_, j_, k_) = (i as isize, j as isize, k as isize);
                    let centre = old.mirrored(i_, j_, k_);
                    let sum = old.mirrored(i_ - 1, j_, k_)
                        + old.mirrored(i_ + 1, j_, k_)
                        + old.mirrored(i_, j_ - 1, k_)
                        + old.mirrored(i_, j_ + 1, k_)
                        + old.mirrored(i_, j_, k_ - 1)
                        + old.mirrored(i_, j_, k_ + 1);
                    self.cells[i + nx * (j + ny * k)] = centre + f * (sum - 6.0 * centre);
                }
            }
        }
        Ok(())
    }

    /// Heat gained since the start. The faces are insulated, so this is exactly what came in.
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.stored_heat())
    }

    fn checkpoint(&mut self) {
        self.saved.clone_from(&self.cells);
    }

    fn restore(&mut self) {
        self.cells.clone_from(&self.saved);
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
        vec![
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
        ]
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

    /// Trilinear between cell centres, clamped at the faces.
    ///
    /// Clamped rather than extrapolated: outside an insulated face the temperature is not
    /// defined, and continuing the gradient would draw a block hotter than any cell in it.
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let (nx, ny, nz) = self.counts;
        let q = p.to_si() / self.dx.to_si() - DVec3::splat(0.5);
        // NaN spelled out rather than folded into a comparison: a visualiser can hand one over,
        // and it must not reach the cast below.
        if q.is_nan() {
            return self.cells[0];
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
        let (i1, j1, k1) = (
            (i + 1).min(nx - 1),
            (j + 1).min(ny - 1),
            (k + 1).min(nz - 1),
        );

        let get = |a: usize, b: usize, c: usize| self.cells[a + nx * (b + ny * c)];
        let lerp = |lo: f64, hi: f64, t: f64| lo * (1.0 - t) + hi * t;
        let z0 = lerp(
            lerp(get(i, j, k), get(i1, j, k), fx),
            lerp(get(i, j1, k), get(i1, j1, k), fx),
            fy,
        );
        let z1 = lerp(
            lerp(get(i, j, k1), get(i1, j, k1), fx),
            lerp(get(i, j1, k1), get(i1, j1, k1), fx),
            fy,
        );
        lerp(z0, z1, fz)
    }

    /// Central differences on the cell grid, mirrored at the faces.
    fn gradient(&self, p: LengthVec, _t: Time, _h: Length) -> DVec3 {
        let (i, j, k) = self.nearest_cell(p);
        let d = 2.0 * self.dx.to_si();
        DVec3::new(
            (self.mirrored(i + 1, j, k) - self.mirrored(i - 1, j, k)) / d,
            (self.mirrored(i, j + 1, k) - self.mirrored(i, j - 1, k)) / d,
            (self.mirrored(i, j, k + 1) - self.mirrored(i, j, k - 1)) / d,
        )
    }

    /// The seven-point stencil — the same one [`Domain::step`] uses, which is the point.
    fn laplacian(&self, p: LengthVec, _t: Time, _h: Length) -> f64 {
        let (i, j, k) = self.nearest_cell(p);
        let dx = self.dx.to_si();
        let centre = self.mirrored(i, j, k);
        let sum = self.mirrored(i - 1, j, k)
            + self.mirrored(i + 1, j, k)
            + self.mirrored(i, j - 1, k)
            + self.mirrored(i, j + 1, k)
            + self.mirrored(i, j, k - 1)
            + self.mirrored(i, j, k + 1);
        (sum - 6.0 * centre) / (dx * dx)
    }

    /// `∂T/∂t = α∇²T`. Conduction only — heat arriving over the bus is a source this cannot see.
    fn rate(&self, p: LengthVec, t: Time, _dt: Time) -> f64 {
        let Some(alpha) = self.substance.diffusivity() else {
            return 0.0;
        };
        alpha.to_si() * self.laplacian(p, t, self.dx)
    }
}

impl Solid3D {
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
