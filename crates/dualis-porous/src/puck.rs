//! The bed itself: a cylinder of packed solid in a basket, with liquid driven through it.

use dualis_core::conserved::quantity;
use dualis_core::{Domain, Exchange, Kind, Ledger, Reading, ScalarField, Substance, Violation};
use dualis_units::{
    Density, Energy, Length, LengthVec, Mass, MassFlow, Power, Pressure, Temperature,
    ThermalConductivity, Time, Velocity,
};
use glam::DVec3;

use crate::{Bed, Grind, Liquid};

/// Conjugate gradients gets this many iterations per cell before it gives up.
const ITERATION_BUDGET: usize = 4;

/// The reference the enthalpy ledger is measured from, K.
///
/// Not absolute zero. The bed sits near 360 K and the changes worth seeing are of order 1 K, so a
/// ledger of absolute enthalpies would be subtracting two numbers that agree to five digits every
/// step — the arithmetic that cost this workspace three orders of magnitude on the GPU port. The
/// stencil and the advection are both linear, so shifting the origin commutes with them exactly.
const T_REF: f64 = 273.15;

/// Which of the bed's several simultaneous fields to look at.
///
/// A packed bed under flow has a temperature, a pressure, a speed and an extraction state at every
/// point, and all four are true at once. [`Domain::as_field`] can nominate one;
/// [`Puck::field`] hands over any of them, and `dualis_scene::sample_field` turns one into a panel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Observable {
    /// Temperature, K.
    Temperature,
    /// Pressure, Pa. Zero at the outlet.
    Pressure,
    /// Speed of the liquid in the pores, m/s. **Not** the Darcy velocity — see
    /// [`Puck::pore_velocity_at`].
    Speed,
    /// How much of what could dissolve already has, as a fraction from 0 to 1.
    ///
    /// This is the field that shows channelling: an even puck extracts evenly, and a channel is
    /// visible as a column that reached 1 while the bed beside it is still at 0.4.
    Extraction,
    /// Solute concentration in the pore liquid, kg/m³.
    Concentration,
}

/// A cylinder of packed, wetted solid in a basket, driven by a pressure difference.
///
/// The grid is a box and the basket is a cylinder in it, of a radius you give. Cells outside that
/// radius are **the basket wall** — a real material with a real heat capacity, no flow through it
/// and nothing to extract. That is not decoration: a portafilter is a couple of hundred grams of
/// metal against eighteen grams of coffee, and it is the reason a shot pulled into a cold one is a
/// different drink.
///
/// Leave the box a few cells wider than the basket. A cylinder that fills the grid puts the metal
/// in the corners, where it has the right heat capacity and the wrong shape — and a cut through
/// the axis then crosses none of it, which is visible the first time anybody draws one.
///
/// # Axes
///
/// **`y` is the flow axis.** `j = 0` is the inlet — the shower screen, held at the pump pressure —
/// and `j = ny−1` is the outlet, the basket's holes, at atmosphere. `x` and `z` are the radial
/// plane. So a slice at fixed `k` is a **vertical cross-section** through the puck, and the one at
/// `k = nz/2` cuts through the axis.
#[derive(Clone, Debug)]
pub struct Puck {
    name: String,
    counts: (usize, usize, usize),
    dx: f64,
    bed: Bed,
    liquid: Liquid,
    wall: Substance,
    /// Cell is inside the basket: packed solid rather than wall.
    packed: Vec<bool>,
    porosity: Vec<f64>,
    grind: Vec<Grind>,
    /// Pressure, Pa. The solve's output.
    pressure: Vec<f64>,
    /// Temperature, K.
    temperature: Vec<f64>,
    /// Extractable solute still inside the particles, kg.
    solids: Vec<f64>,
    /// What each cell started with, kg. The reference the equilibrium isotherm is written against.
    initial_solids: Vec<f64>,
    /// Heat capacity per cell, J/K. Cached: it moves only when the packing does.
    capacity: Vec<f64>,
    /// Thermal conductivity per cell, W/m/K. Cached for the same reason.
    lambda: Vec<f64>,
    /// The face list, rebuilt at the start of every solve. Rebuilding it inside the CG loop cost
    /// more than the solve did.
    face_cache: Vec<(Side, Side, f64)>,
    /// Solute dissolved in this cell's pore liquid, kg.
    dissolved: Vec<f64>,
    drive: f64,
    inlet_temperature: f64,
    tolerance: f64,
    max_iterations: Option<usize>,
    residual: f64,
    converged: bool,
    /// Totals over the run.
    delivered_volume: f64,
    delivered_solute: f64,
    delivered_enthalpy: f64,
    admitted_enthalpy: f64,
    /// The dry mass the basket was dosed with, kg. Fixed at construction.
    dose: f64,
    /// What was extractable at the start, kg.
    extractable: f64,
    /// Instantaneous flow, m³/s, from the last solve.
    flow: f64,
    saved: Option<Box<Saved>>,
}

#[derive(Clone, Debug)]
struct Saved {
    pressure: Vec<f64>,
    temperature: Vec<f64>,
    solids: Vec<f64>,
    dissolved: Vec<f64>,
    delivered_volume: f64,
    delivered_solute: f64,
    delivered_enthalpy: f64,
    admitted_enthalpy: f64,
}

/// One end of a face: a cell, or one of the two boundaries the flow is driven between.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Side {
    Cell(usize),
    Inlet,
    Outlet,
}

/// Everything a basket is, in one place.
///
/// # Why a struct and not ten arguments
///
/// It was ten arguments, and `clippy::too_many_arguments` was silenced to allow it. Two of them
/// were `Length` and two were temperatures, so transposing a pair compiled and ran. This is the
/// shape the workspace's own consumer keeps arriving at from the other direction: a caller who
/// wants to change one thing should not have to restate nine.
///
/// [`Basket::espresso`] is a conventional double basket. Change the field you are asking about.
#[derive(Clone, Debug)]
pub struct Basket {
    /// Cells along `(x, y, z)`. **`y` is the flow axis**, so `ny` sets the bed's depth.
    pub counts: (usize, usize, usize),
    /// Cell side. Cells are cubes.
    pub cell: Length,
    /// The basket's inside radius. Cells outside it are wall.
    ///
    /// Explicit rather than "the largest circle that fits", which was the first version and was
    /// wrong in a way that only showed in a picture: a circle inscribed in the grid touches the
    /// box at the four axes, so the metal lives in the *corners* and a cut through the axis
    /// crosses no metal at all. A basket has a jacket all the way round. Leave a few cells for it.
    pub radius: Length,
    /// What the grounds are.
    pub bed: Bed,
    /// How finely they are ground.
    pub grind: Grind,
    /// Inter-particle void fraction after tamping. About 0.45 for espresso.
    pub porosity: f64,
    /// What is being pushed through.
    pub liquid: Liquid,
    /// What the basket is made of.
    pub wall: Substance,
    /// The pressure held across the bed.
    pub pressure: Pressure,
    /// The starting temperature of everything — bed, liquid and wall alike. Move the wall on its
    /// own with [`Puck::set_wall_temperature`].
    pub temperature: Temperature,
}

impl Basket {
    /// A conventional 58 mm double basket, 20 mm deep, on a 2 mm grid with a 4 mm jacket.
    ///
    /// 17.6 g of 250 µm grind at `ε = 0.45`, 9 bar, 93 °C — the settings the crate's two
    /// calibrated numbers were fitted to, and the ones a shot is compared against.
    pub fn espresso() -> Basket {
        Basket {
            counts: (33, 10, 33),
            cell: Length::from_si(2e-3),
            radius: Length::from_si(29e-3),
            bed: Bed::roasted_coffee(),
            grind: Grind::espresso(),
            porosity: 0.45,
            liquid: Liquid::water(),
            wall: Substance::stainless_304(),
            pressure: Pressure::from_si(9.0e5),
            temperature: Temperature::celsius(93.0),
        }
    }
}

impl Puck {
    /// A basket, dosed and wetted and ready to run.
    ///
    /// The dose is **not** a parameter: it is `(1−ε)·ρ_particle·V_bed`, which is what that volume
    /// of that grind packed to that porosity actually weighs. Setting a dose *and* a porosity
    /// *and* a depth independently would let a caller specify a bed that does not exist.
    pub fn new(name: impl Into<String>, basket: Basket) -> Puck {
        let Basket {
            counts,
            cell,
            radius,
            bed,
            grind,
            porosity,
            liquid,
            wall,
            pressure,
            temperature,
        } = basket;
        let counts = (counts.0.max(1), counts.1.max(1), counts.2.max(1));
        let (nx, ny, nz) = counts;
        let cells = nx * ny * nz;
        let dx = cell.to_si();
        let porosity = porosity.clamp(0.01, 0.99);

        // The basket, tested at cell centres, centred on the grid.
        let radius = radius.to_si().min(0.5 * (nx.min(nz) as f64) * dx);
        let (cx, cz) = (0.5 * nx as f64 * dx, 0.5 * nz as f64 * dx);
        let mut packed = vec![false; cells];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let (x, z) = ((i as f64 + 0.5) * dx - cx, (k as f64 + 0.5) * dx - cz);
                    packed[i + nx * (j + ny * k)] = x * x + z * z <= radius * radius;
                }
            }
        }

        let volume = dx * dx * dx;
        let solid_per_cell = (1.0 - porosity) * bed.solid_density.to_si() * volume;
        let dose = solid_per_cell * packed.iter().filter(|p| **p).count() as f64;
        let extractable = dose * bed.soluble_fraction;
        let solids: Vec<f64> = packed
            .iter()
            .map(|p| {
                if *p {
                    solid_per_cell * bed.soluble_fraction
                } else {
                    0.0
                }
            })
            .collect();

        let mut built = Puck {
            name: name.into(),
            counts,
            dx,
            bed,
            liquid,
            wall,
            porosity: packed
                .iter()
                .map(|p| if *p { porosity } else { 0.0 })
                .collect(),
            grind: vec![grind; cells],
            packed,
            pressure: vec![0.0; cells],
            temperature: vec![temperature.to_si(); cells],
            initial_solids: solids.clone(),
            solids,
            capacity: vec![0.0; cells],
            lambda: vec![0.0; cells],
            face_cache: Vec::new(),
            dissolved: vec![0.0; cells],
            drive: pressure.to_si(),
            inlet_temperature: temperature.to_si(),
            tolerance: 1e-13,
            max_iterations: None,
            residual: f64::INFINITY,
            converged: false,
            delivered_volume: 0.0,
            delivered_solute: 0.0,
            delivered_enthalpy: 0.0,
            admitted_enthalpy: 0.0,
            dose,
            extractable,
            flow: 0.0,
            saved: None,
        };
        built.refresh_properties();
        // Solved at construction, for the reason `Conductor` records: a quasi-static field read
        // before its solve is a field of zeros wearing the shape of an answer.
        built.solve(built.tolerance);
        built
    }

    /// The temperature of the liquid arriving at the inlet.
    ///
    /// Separate from the bed's own temperature, because they are separate in the machine: the
    /// group head delivers water at whatever the boiler and the path to it produce, and the basket
    /// is at whatever it was left at.
    pub fn set_inlet_temperature(&mut self, t: Temperature) {
        self.inlet_temperature = t.to_si();
    }

    /// The pressure across the bed.
    pub fn set_drive(&mut self, p: Pressure) {
        self.drive = p.to_si();
        self.resolve();
    }

    /// Set the temperature of every cell, bed and wall alike.
    pub fn set_temperature(&mut self, t: Temperature) {
        self.temperature.fill(t.to_si());
        self.resolve();
    }

    /// Set the temperature of the basket wall only, leaving the bed where it is.
    ///
    /// The cold-portafilter case, which is the one worth being able to state separately.
    pub fn set_wall_temperature(&mut self, t: Temperature) {
        for (idx, packed) in self.packed.iter().enumerate() {
            if !packed {
                self.temperature[idx] = t.to_si();
            }
        }
        self.resolve();
    }

    /// Change the porosity of the cells a predicate selects.
    ///
    /// This is how a defect is put into a puck: a loosely packed column against the wall, a crack
    /// from a knocked basket, a dry patch. Nothing else in this crate makes a puck uneven, because
    /// an even puck is what an even tamp gives and channelling is a *fault*, not a phenomenon the
    /// physics produces on its own.
    pub fn repack(&mut self, porosity: f64, which: impl Fn(usize, usize, usize) -> bool) {
        let (nx, ny, nz) = self.counts;
        let volume = self.dx.powi(3);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = i + nx * (j + ny * k);
                    if self.packed[idx] && which(i, j, k) {
                        self.porosity[idx] = porosity.clamp(0.01, 0.99);
                        let solid =
                            (1.0 - self.porosity[idx]) * self.bed.solid_density.to_si() * volume;
                        self.solids[idx] = solid * self.bed.soluble_fraction;
                    }
                }
            }
        }
        self.recount_dose();
        self.refresh_properties();
        self.resolve();
    }

    /// Change the grind of the cells a predicate selects.
    pub fn regrind(&mut self, grind: Grind, which: impl Fn(usize, usize, usize) -> bool) {
        let (nx, ny, nz) = self.counts;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = i + nx * (j + ny * k);
                    if self.packed[idx] && which(i, j, k) {
                        self.grind[idx] = grind;
                    }
                }
            }
        }
        self.resolve();
    }

    fn recount_dose(&mut self) {
        let volume = self.dx.powi(3);
        self.dose = (0..self.porosity.len())
            .filter(|i| self.packed[*i])
            .map(|i| (1.0 - self.porosity[i]) * self.bed.solid_density.to_si() * volume)
            .sum();
        self.extractable = self.dose * self.bed.soluble_fraction;
    }

    /// Re-solve after something that changed the flow problem.
    ///
    /// # Every mutator ends solved, and this is why
    ///
    /// The constructor solves, for the reason `Conductor` writes down: a quasi-static field read
    /// before its solve is a field of zeros wearing the shape of an answer. Marking the solve
    /// stale and leaving it to the next [`Domain::step`] is worse than that, not better — the
    /// field left behind is not zeros but *the previous answer*, which is smooth, bounded and
    /// exactly the right order of magnitude.
    ///
    /// It was measured. `repack` widening the ring to 0.60 left `flow_rate` reporting the even
    /// bed's flow to the last digit — a 42% error, and 1.0000 where the closed form says 1.4188.
    /// Nothing in the suite saw it, because every test stepped the puck before reading anything
    /// from it and a step re-solves.
    ///
    /// Re-solving costs one warm-started conjugate-gradient run, which from the previous field is
    /// a handful of iterations. None of these methods is in a hot loop.
    fn resolve(&mut self) {
        self.converged = false;
        self.residual = f64::INFINITY;
        self.solve(self.tolerance);
    }

    /// Recompute the per-cell heat capacity and conductivity.
    ///
    /// Cached rather than computed in the inner loop: both are read once per face per step, and
    /// the wall's come through an `Option` and a `to_si`. Neither depends on temperature, so a
    /// cache is exact rather than a lag.
    fn refresh_properties(&mut self) {
        let volume = self.dx.powi(3);
        let wall_c = self
            .wall
            .thermal
            .as_ref()
            .map(|t| t.specific_heat.to_si())
            .unwrap_or(500.0);
        let wall_k = self
            .wall
            .thermal
            .as_ref()
            .map(|t| t.conductivity.to_si())
            .unwrap_or(15.0);
        for idx in 0..self.capacity.len() {
            if self.packed[idx] {
                let e = self.porosity[idx];
                self.capacity[idx] =
                    volume * (e * self.liquid.rho_c() + (1.0 - e) * self.bed.rho_c());
                self.lambda[idx] = Self::bed_conductivity(
                    self.liquid.conductivity.to_si(),
                    self.bed.conductivity.to_si(),
                    e,
                );
            } else {
                self.capacity[idx] = volume * self.wall.density.to_si() * wall_c;
                self.lambda[idx] = wall_k;
            }
        }
    }

    /// The effective thermal conductivity of a saturated bed, by **Maxwell–Eucken with the liquid as
    /// the continuous phase**.
    ///
    /// ```text
    ///   k = k_l · (2k_l + k_s − 2(1−ε)(k_l − k_s)) / (2k_l + k_s + (1−ε)(k_l − k_s))
    /// ```
    ///
    /// # Why this one, and what it replaced
    ///
    /// It was `ε k_l + (1−ε) k_s` — the arithmetic mean, which is the **Voigt bound**. That is not a
    /// model of a packed bed; it is what you get if you do not think about it, and it is an *upper*
    /// bound: it is exact only when the two phases lie in parallel with the flux, which in a bed of
    /// spheres nothing does.
    ///
    /// The structural fact about a saturated bed is that **the liquid is the continuous phase** and the
    /// grains are dispersed in it. That is precisely what Maxwell–Eucken describes, and it is attained by
    /// a coated-sphere assemblage — a closed form, not a correlation. For coffee at 45% porosity the two
    /// differ by 11.0%: 0.38625 against 0.34811 W/m·K.
    ///
    /// Note that "liquid as host" is a statement about the geometry and not about which number is
    /// larger. For water in coffee it happens to coincide with the Hashin–Shtrikman *upper* bound,
    /// because the water conducts better than the grounds; for a metal powder it would be the lower one.
    /// `the_beds_conductivity.rs` asserts both the coincidence for these numbers and the bracketing in
    /// general, against `dualis_core::mixture::Mix`, which is checked against Maxwell–Garnett to
    /// `6.7e-16`. Two independently written forms of the same physics, in crates that do not depend on
    /// each other for it.
    ///
    /// # What the choice is worth, measured
    ///
    /// **Nothing, in every scenario this crate ships**, and that is why it went unexamined. Under flow
    /// the bed is isothermal, so conduction carries no heat and `λ` is multiplied by zero: swinging it
    /// over a factor of eight leaves the extraction yield identical to `1e-14`.
    ///
    /// It stops being free the moment there is a gradient. In a 20 °C basket the yield moves **4.9% per
    /// unit `ln λ`** on the grid `the_beds_conductivity.rs` measures it on, so the honest range the old
    /// rule left — Voigt to Reuss, a factor of 1.674 — was worth about **2.5% in extraction yield** on
    /// that grid, which is a taste-level difference in a cup. Maxwell–Eucken
    /// narrows the range that remains to 1.184 and costs 0.40% against the old value.
    pub fn bed_conductivity(liquid: f64, solid: f64, porosity: f64) -> f64 {
        let d = liquid - solid;
        let solids = 1.0 - porosity;
        let denominator = 2.0 * liquid + solid + solids * d;
        if denominator <= 0.0 {
            // A liquid of zero conductivity in a solid of zero conductivity. Not reachable from a
            // validated substance, and zero is the right answer rather than a `NaN`.
            return 0.0;
        }
        liquid * (2.0 * liquid + solid - 2.0 * solids * d) / denominator
    }

    /// The effective conductivity this cell is conducting with, in W/m·K.
    ///
    /// Exposed because it is a **modelling choice** rather than an implementation detail — see
    /// [`Puck::bed_conductivity`] for which choice and what the alternatives cost. A number that decides
    /// an answer and cannot be read from outside is a number nobody can check.
    ///
    /// A wall cell reports the wall's own conductivity, unmixed.
    pub fn conductivity_at(&self, i: usize, j: usize, k: usize) -> ThermalConductivity {
        let (nx, ny, _) = self.counts;
        ThermalConductivity::from_si(self.lambda[i + nx * (j + ny * k)])
    }

    /// How many cells the grid has, as `(nx, ny, nz)`. `y` is the flow axis.
    pub fn counts(&self) -> (usize, usize, usize) {
        self.counts
    }

    /// The cell spacing.
    pub fn spacing(&self) -> Length {
        Length::from_si(self.dx)
    }

    /// The bed's depth along the flow axis.
    pub fn depth(&self) -> Length {
        Length::from_si(self.counts.1 as f64 * self.dx)
    }

    /// The dry mass the basket holds.
    pub fn dose(&self) -> Mass {
        Mass::from_si(self.dose)
    }

    /// Whether a cell is packed bed rather than basket wall.
    pub fn is_packed(&self, i: usize, j: usize, k: usize) -> bool {
        self.index(i, j, k).map(|n| self.packed[n]).unwrap_or(false)
    }

    /// The flat index of a cell, or `None` out of range.
    pub fn index(&self, i: usize, j: usize, k: usize) -> Option<usize> {
        let (nx, ny, nz) = self.counts;
        (i < nx && j < ny && k < nz).then(|| i + nx * (j + ny * k))
    }

    fn clamped(&self, i: usize, j: usize, k: usize) -> usize {
        let (nx, ny, nz) = self.counts;
        i.min(nx - 1) + nx * (j.min(ny - 1) + ny * k.min(nz - 1))
    }

    /// Temperature at one cell.
    pub fn temperature_at(&self, i: usize, j: usize, k: usize) -> Temperature {
        Temperature::from_si(self.temperature[self.clamped(i, j, k)])
    }

    /// Pressure at one cell.
    pub fn pressure_at(&self, i: usize, j: usize, k: usize) -> Pressure {
        Pressure::from_si(self.pressure[self.clamped(i, j, k)])
    }

    /// How much of what could dissolve in this cell already has, from 0 to 1.
    pub fn extraction_at(&self, i: usize, j: usize, k: usize) -> f64 {
        let idx = self.clamped(i, j, k);
        if !self.packed[idx] {
            return 0.0;
        }
        let start = (1.0 - self.porosity[idx])
            * self.bed.solid_density.to_si()
            * self.dx.powi(3)
            * self.bed.soluble_fraction;
        if start <= 0.0 {
            0.0
        } else {
            ((start - self.solids[idx]) / start).clamp(0.0, 1.0)
        }
    }

    /// Solute concentration in this cell's pore liquid.
    pub fn concentration_at(&self, i: usize, j: usize, k: usize) -> Density {
        let idx = self.clamped(i, j, k);
        let pore = self.porosity[idx] * self.dx.powi(3);
        Density::from_si(if pore > 0.0 {
            self.dissolved[idx] / pore
        } else {
            0.0
        })
    }

    /// The **Darcy** velocity at a cell: volume flow per unit of total cross-section.
    ///
    /// This is what `−(k/μ)∇p` gives directly and it is not the speed anything moves at. A tracer
    /// travels at the *pore* velocity, which is this divided by the porosity — a factor of three
    /// for an espresso puck. Confusing the two is the classic error in porous transport, and it
    /// makes every arrival time three times too late.
    pub fn darcy_velocity_at(&self, i: usize, j: usize, k: usize) -> DVec3 {
        let (nx, ny, nz) = self.counts;
        let (i, j, k) = (i.min(nx - 1), j.min(ny - 1), k.min(nz - 1));
        let here = self.clamped(i, j, k);
        if !self.packed[here] {
            return DVec3::ZERO;
        }
        let grad = |lo: Option<usize>, hi: Option<usize>| -> f64 {
            match (lo, hi) {
                (Some(a), Some(b)) => (self.pressure[b] - self.pressure[a]) / (2.0 * self.dx),
                (None, Some(b)) => (self.pressure[b] - self.pressure[here]) / self.dx,
                (Some(a), None) => (self.pressure[here] - self.pressure[a]) / self.dx,
                (None, None) => 0.0,
            }
        };
        // A neighbour that is wall is not a neighbour: the gradient into it is not a driving
        // force, it is a no-flow face.
        let flowing = |n: Option<usize>| n.filter(|n| self.packed[*n]);
        let g = DVec3::new(
            grad(
                flowing(i.checked_sub(1).and_then(|a| self.index(a, j, k))),
                flowing(self.index(i + 1, j, k)),
            ),
            grad(
                flowing(j.checked_sub(1).and_then(|b| self.index(i, b, k))),
                flowing(self.index(i, j + 1, k)),
            ),
            grad(
                flowing(k.checked_sub(1).and_then(|c| self.index(i, j, c))),
                flowing(self.index(i, j, k + 1)),
            ),
        );
        -self.mobility(here) * g
    }

    /// The speed a tracer actually moves at: the Darcy velocity over the porosity.
    pub fn pore_velocity_at(&self, i: usize, j: usize, k: usize) -> DVec3 {
        let idx = self.clamped(i, j, k);
        let e = self.porosity[idx];
        if e <= 0.0 {
            DVec3::ZERO
        } else {
            self.darcy_velocity_at(i, j, k) / e
        }
    }

    /// The volumetric flow through the bed, from the last solve.
    pub fn flow_rate(&self) -> MassFlow {
        MassFlow::from_si(self.flow * self.liquid.density.to_si())
    }

    /// The liquid delivered so far.
    pub fn delivered(&self) -> Mass {
        Mass::from_si(self.delivered_volume * self.liquid.density.to_si())
    }

    /// The solute delivered so far — the dissolved mass in the cup.
    pub fn delivered_solute(&self) -> Mass {
        Mass::from_si(self.delivered_solute)
    }

    /// Extraction yield: solute in the cup as a fraction of the dry dose.
    ///
    /// The number a barista means by "20%". Its ceiling is [`Bed::soluble_fraction`], not one.
    pub fn yield_fraction(&self) -> f64 {
        if self.dose <= 0.0 {
            0.0
        } else {
            self.delivered_solute / self.dose
        }
    }

    /// Total dissolved solids: solute as a fraction of the beverage mass.
    ///
    /// The number a refractometer reads, and a different question from the yield: yield asks how
    /// much came out of the coffee, TDS asks how concentrated what came out is. A long shot can
    /// have a high yield and a low TDS.
    pub fn tds(&self) -> f64 {
        let beverage = self.delivered_volume * self.liquid.density.to_si() + self.delivered_solute;
        if beverage <= 0.0 {
            0.0
        } else {
            self.delivered_solute / beverage
        }
    }

    /// The temperature of the liquid leaving the basket, weighted by where it leaves.
    pub fn outlet_temperature(&self) -> Temperature {
        let (nx, ny, nz) = self.counts;
        let (mut num, mut den) = (0.0, 0.0);
        for k in 0..nz {
            for i in 0..nx {
                let idx = i + nx * (ny - 1 + ny * k);
                if !self.packed[idx] {
                    continue;
                }
                let q = self.outlet_conductance(idx) * self.pressure[idx];
                num += q * self.temperature[idx];
                den += q;
            }
        }
        Temperature::from_si(if den > 0.0 {
            num / den
        } else {
            self.inlet_temperature
        })
    }

    /// The mean temperature of the packed bed, ignoring the basket wall.
    pub fn bed_temperature(&self) -> Temperature {
        let (mut sum, mut n) = (0.0, 0usize);
        for (idx, packed) in self.packed.iter().enumerate() {
            if *packed {
                sum += self.temperature[idx];
                n += 1;
            }
        }
        Temperature::from_si(if n > 0 {
            sum / n as f64
        } else {
            self.inlet_temperature
        })
    }

    /// How unevenly the bed extracted, as the coefficient of variation of the per-cell extraction.
    ///
    /// # Read this with [`Puck::radial_contrast`], not instead of it
    ///
    /// A perfectly packed bed does **not** give zero here, and expecting it to is the trap. Fresh
    /// water enters at the inlet and is already loaded by the time it reaches the outlet, so an
    /// even puck extracts more at the top than at the bottom — measured, a spread of about 0.10 on
    /// a conventional shot, none of which is a fault.
    ///
    /// That axial gradient dominates this number, which makes it a **poor detector of
    /// channelling**: a wall channel that halves the yield moves it from 0.105 to 0.128. What it
    /// is good for is comparing two shots of the same geometry, where the gradient is common to
    /// both.
    pub fn unevenness(&self) -> f64 {
        let (nx, ny, nz) = self.counts;
        let (mut sum, mut sq, mut n) = (0.0, 0.0, 0usize);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if !self.is_packed(i, j, k) {
                        continue;
                    }
                    let e = self.extraction_at(i, j, k);
                    sum += e;
                    sq += e * e;
                    n += 1;
                }
            }
        }
        if n == 0 {
            return 0.0;
        }
        let mean = sum / n as f64;
        if mean <= 0.0 {
            return 0.0;
        }
        (sq / n as f64 - mean * mean).max(0.0).sqrt() / mean
    }

    /// How much more the outer ring of the bed extracted than its core.
    ///
    /// One for an evenly packed basket. Above one means the flow preferred the wall, which is what
    /// a channel *is* — and unlike [`Puck::unevenness`] this is blind to the axial gradient,
    /// because the ring and the core span the same depths.
    ///
    /// The ring is the packed cells with a wall neighbour in the radial plane; the core is the
    /// rest. That is the shape the defect actually takes when a puck shrinks away from the basket.
    pub fn radial_contrast(&self) -> f64 {
        let (nx, ny, nz) = self.counts;
        let (mut ring, mut ring_n, mut core, mut core_n) = (0.0, 0usize, 0.0, 0usize);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if !self.is_packed(i, j, k) {
                        continue;
                    }
                    let edge = i == 0
                        || k == 0
                        || i + 1 == nx
                        || k + 1 == nz
                        || !self.is_packed(i - 1, j, k)
                        || !self.is_packed(i + 1, j, k)
                        || !self.is_packed(i, j, k - 1)
                        || !self.is_packed(i, j, k + 1);
                    let e = self.extraction_at(i, j, k);
                    if edge {
                        ring += e;
                        ring_n += 1;
                    } else {
                        core += e;
                        core_n += 1;
                    }
                }
            }
        }
        if ring_n == 0 || core_n == 0 || core <= 0.0 {
            return 1.0;
        }
        (ring / ring_n as f64) / (core / core_n as f64)
    }

    /// How much the flow in disagrees with the flow out, relative to the flow.
    ///
    /// Zero for a converged solve. Measured at the two boundaries independently rather than taken
    /// from the solver's own residual, which is the same check
    /// [`Conductor::current_balance`](https://docs.rs/dualis-electrical) makes and for the same
    /// reason: an iterative solver's opinion of itself is not evidence.
    pub fn flow_balance(&self) -> f64 {
        let (a, b) = (self.boundary_flow(true), self.boundary_flow(false));
        let scale = a.abs().max(b.abs());
        if scale <= 0.0 {
            0.0
        } else {
            (a - b).abs() / scale
        }
    }

    /// The relative residual the last pressure solve reached.
    pub fn residual(&self) -> f64 {
        self.residual
    }

    /// Whether the last pressure solve met its tolerance.
    pub fn converged(&self) -> bool {
        self.converged
    }

    /// The effective dispersivity the advection scheme imposes, `dx/2`.
    ///
    /// Reported rather than hidden. First-order upwind carries a numerical diffusivity of
    /// `u·dx/2`, which is exactly the form of the mechanical dispersion a packed bed really has —
    /// `α_L·u`, with `α_L` of order the particle diameter. So the scheme *is* the dispersion model,
    /// and whether that model is right is the question of whether `dx/2` is near the grind.
    ///
    /// At 1 mm cells and a 250 µm grind it is twice too large, so a solute front is smeared about
    /// twice as much as it should be. Refining the grid reduces it; there is no term to switch off.
    pub fn numerical_dispersivity(&self) -> Length {
        Length::from_si(self.dx / 2.0)
    }

    /// One of the bed's fields, as something a sampler can read.
    pub fn field(&self, what: Observable) -> PuckField<'_> {
        PuckField { puck: self, what }
    }

    /// Solve for the pressure, to a relative residual of `tolerance`.
    pub fn solve(&mut self, tolerance: f64) -> bool {
        let budget = self
            .max_iterations
            .unwrap_or(ITERATION_BUDGET * self.pressure.len() + 32);
        self.solve_within(tolerance, budget)
    }

    /// Solve, spending at most `max_iterations`. Returns whether the tolerance was met.
    pub fn solve_within(&mut self, tolerance: f64, max_iterations: usize) -> bool {
        let n = self.pressure.len();
        self.refresh_faces();
        let faces = std::mem::take(&mut self.face_cache);
        let b = source(&faces, n, self.drive);
        let mut x = std::mem::take(&mut self.pressure);
        if x.len() != n {
            x = vec![0.0; n];
        }
        let mut r = b.clone();
        let ax = apply(&faces, &x);
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
        while rr.sqrt() / scale > tolerance && iterations < max_iterations {
            let ap = apply(&faces, &p);
            let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
            if pap <= 0.0 {
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

        self.pressure = x;
        self.face_cache = faces;
        self.residual = rr.sqrt() / scale;
        self.converged = self.residual <= tolerance;
        self.flow = self.boundary_flow(true);
        self.converged
    }

    // --- the discretisation ------------------------------------------------

    /// Mobility `k/μ` at a cell, m²·Pa⁻¹·s⁻¹. Zero for wall.
    fn mobility(&self, idx: usize) -> f64 {
        if !self.packed[idx] {
            return 0.0;
        }
        let k = self.grind[idx].permeability(self.porosity[idx]).to_si();
        let mu = self
            .liquid
            .viscosity(Temperature::from_si(self.temperature[idx]))
            .to_si();
        k / mu
    }

    /// Heat capacity of a cell, J·K⁻¹.
    fn heat_capacity(&self, idx: usize) -> f64 {
        self.capacity[idx]
    }

    /// Thermal conductivity of a cell, W·m⁻¹·K⁻¹.
    ///
    /// For the bed this is the volume-weighted mean of liquid and solid — the parallel model,
    /// which is an upper bound on the true effective conductivity of a packed bed. It is chosen
    /// because it errs towards the bed equilibrating with its wall faster rather than slower, and
    /// a thermal claim should be the conservative one.
    fn conductivity(&self, idx: usize) -> f64 {
        self.lambda[idx]
    }

    /// Hydraulic conductance from a boundary-adjacent cell to the boundary, m³·Pa⁻¹·s⁻¹.
    fn boundary_conductance(&self, idx: usize) -> f64 {
        self.mobility(idx) * self.dx * self.dx / (0.5 * self.dx)
    }

    fn outlet_conductance(&self, idx: usize) -> f64 {
        self.boundary_conductance(idx)
    }

    /// Every face carrying flow, as `(a, b, conductance)`.
    ///
    /// The conductance across a face between two cells is the **harmonic** mean of their
    /// mobilities, because two half-cells are in series. An arithmetic mean is the standard
    /// mistake and it is invisible in a uniform bed — which is why a wall cell, whose mobility is
    /// exactly zero, must give exactly zero here rather than half of its neighbour's.
    fn faces(&self) -> Vec<(Side, Side, f64)> {
        let (nx, ny, nz) = self.counts;
        let area = self.dx * self.dx;
        let mut out = Vec::new();
        for k in 0..nz {
            for i in 0..nx {
                // `j = 0` spelled out: the inlet is the low face of the flow axis, and reading
                // `nx * (ny * k)` back as that is harder than it needs to be.
                let inlet = i + nx * ny * k;
                out.push((
                    Side::Inlet,
                    Side::Cell(inlet),
                    self.boundary_conductance(inlet),
                ));
                let outlet = i + nx * (ny - 1 + ny * k);
                out.push((
                    Side::Cell(outlet),
                    Side::Outlet,
                    self.boundary_conductance(outlet),
                ));
            }
        }
        let interior = |a: usize, b: usize, out: &mut Vec<(Side, Side, f64)>| {
            let (ma, mb) = (self.mobility(a), self.mobility(b));
            let g = if ma <= 0.0 || mb <= 0.0 {
                0.0
            } else {
                area / (0.5 * self.dx / ma + 0.5 * self.dx / mb)
            };
            out.push((Side::Cell(a), Side::Cell(b), g));
        };
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx - 1 {
                    interior(i + nx * (j + ny * k), i + 1 + nx * (j + ny * k), &mut out);
                }
            }
        }
        for k in 0..nz {
            for j in 0..ny - 1 {
                for i in 0..nx {
                    interior(i + nx * (j + ny * k), i + nx * (j + 1 + ny * k), &mut out);
                }
            }
        }
        for k in 0..nz - 1 {
            for j in 0..ny {
                for i in 0..nx {
                    interior(i + nx * (j + ny * k), i + nx * (j + ny * (k + 1)), &mut out);
                }
            }
        }
        out
    }

    /// Rebuild [`Puck::face_cache`] from the current mobilities.
    fn refresh_faces(&mut self) {
        self.face_cache = self.faces();
    }

    fn pressure_of(&self, s: Side) -> f64 {
        match s {
            Side::Cell(i) => self.pressure[i],
            Side::Inlet => self.drive,
            Side::Outlet => 0.0,
        }
    }

    /// Total volumetric flow across the inlet (`true`) or the outlet.
    fn boundary_flow(&self, inlet: bool) -> f64 {
        let mut total = 0.0;
        for (a, b, g) in &self.face_cache {
            let (a, b, g) = (*a, *b, *g);
            match (a, b) {
                (Side::Inlet, _) if inlet => total += g * (self.drive - self.pressure_of(b)),
                (_, Side::Outlet) if !inlet => total += g * (self.pressure_of(a) - 0.0),
                _ => {}
            }
        }
        total
    }

    /// The largest step every cell stays positive at, from the last solve.
    ///
    /// The positive-coefficient bound for explicit upwind advection with explicit conduction: a
    /// cell may not lose more than it has in one step. That is `1/rate` where `rate` is the sum of
    /// everything leaving it per unit of what it holds — and it is a bound, not a guideline, so
    /// [`Domain::step`] refuses above it rather than producing a negative concentration.
    fn limiting_rate(&self) -> f64 {
        let (nx, ny, nz) = self.counts;
        let area = self.dx * self.dx;
        let mut heat_rate = vec![0.0f64; self.pressure.len()];
        let mut solute_rate = vec![0.0f64; self.pressure.len()];
        let rho_c_w = self.liquid.rho_c();

        for (a, b, g) in &self.face_cache {
            let (a, b, g) = (*a, *b, *g);
            let q = g * (self.pressure_of(a) - self.pressure_of(b));
            // Advection leaves the upwind side.
            match (a, b, q > 0.0) {
                (Side::Cell(i), _, true) | (_, Side::Cell(i), false) => {
                    heat_rate[i] += q.abs() * rho_c_w;
                    solute_rate[i] += q.abs();
                }
                _ => {}
            }
        }
        // Conduction across every face, both ways.
        let face_conduction = |a: usize, b: usize, rates: &mut Vec<f64>| {
            let (ka, kb) = (self.conductivity(a), self.conductivity(b));
            let g = if ka <= 0.0 || kb <= 0.0 {
                0.0
            } else {
                area / (0.5 * self.dx / ka + 0.5 * self.dx / kb)
            };
            rates[a] += g;
            rates[b] += g;
        };
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx - 1 {
                    face_conduction(
                        i + nx * (j + ny * k),
                        i + 1 + nx * (j + ny * k),
                        &mut heat_rate,
                    );
                }
            }
        }
        for k in 0..nz {
            for j in 0..ny - 1 {
                for i in 0..nx {
                    face_conduction(
                        i + nx * (j + ny * k),
                        i + nx * (j + 1 + ny * k),
                        &mut heat_rate,
                    );
                }
            }
        }
        for k in 0..nz - 1 {
            for j in 0..ny {
                for i in 0..nx {
                    face_conduction(
                        i + nx * (j + ny * k),
                        i + nx * (j + ny * (k + 1)),
                        &mut heat_rate,
                    );
                }
            }
        }

        let mut worst: f64 = 0.0;
        for idx in 0..self.pressure.len() {
            worst = worst.max(heat_rate[idx] / self.heat_capacity(idx));
            let pore = self.porosity[idx] * self.dx.powi(3);
            if pore > 0.0 {
                worst = worst.max(solute_rate[idx] / pore);
            }
        }
        worst
    }
}

/// `b` for the finite-volume system: what the driven boundary injects.
fn source(faces: &[(Side, Side, f64)], n: usize, drive: f64) -> Vec<f64> {
    let mut b = vec![0.0; n];
    for (a, c, g) in faces {
        match (a, c) {
            (Side::Inlet, Side::Cell(i)) => b[*i] += g * drive,
            (Side::Cell(i), Side::Inlet) => b[*i] += g * drive,
            _ => {}
        }
    }
    b
}

/// `A·x`: the sum over faces of `G·(x_a − x_b)`, with the boundaries on the diagonal.
fn apply(faces: &[(Side, Side, f64)], x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; x.len()];
    for (a, b, g) in faces {
        match (a, b) {
            (Side::Cell(i), Side::Cell(j)) => {
                let d = g * (x[*i] - x[*j]);
                y[*i] += d;
                y[*j] -= d;
            }
            (Side::Cell(i), _) | (_, Side::Cell(i)) => y[*i] += g * x[*i],
            _ => {}
        }
    }
    y
}

impl Domain for Puck {
    fn name(&self) -> &str {
        &self.name
    }

    /// [`Kind::Evolving`], even though the pressure is solved rather than marched.
    ///
    /// The flow has no state; the heat and the extraction do, and they are what sets the step.
    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    fn max_stable_dt(&self, _now: Time) -> Time {
        let rate = self.limiting_rate();
        Time::from_si(if rate > 0.0 {
            1.0 / rate
        } else {
            f64::INFINITY
        })
    }

    fn step(&mut self, _t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        // The viscosity moved with the temperature, so the flow field is not the one from last
        // step. Re-solve. Warm-started from the previous pressure, this is a few iterations.
        let budget = self
            .max_iterations
            .unwrap_or(ITERATION_BUDGET * self.pressure.len() + 32);
        if !self.solve_within(self.tolerance, budget) {
            return Err(Violation::at(
                self.name.clone(),
                "pressure residual",
                self.residual,
            ));
        }

        let h = dt.to_si();
        let limit = 1.0 / self.limiting_rate().max(f64::MIN_POSITIVE);
        if h > limit {
            return Err(Violation {
                quantity: "transport step".into(),
                site: self.name.clone(),
                before: limit,
                after: h,
                scale: limit,
                tolerance: 0.0,
            });
        }

        let n = self.pressure.len();
        let area = self.dx * self.dx;
        let rho_c_w = self.liquid.rho_c();
        let mut d_heat = vec![0.0f64; n];
        let mut d_solute = vec![0.0f64; n];
        let (mut out_volume, mut out_solute, mut out_enthalpy, mut in_enthalpy) =
            (0.0, 0.0, 0.0, 0.0);

        // --- flow, in flux form so that whatever leaves one cell arrives in another ---------
        for (a, b, g) in std::mem::take(&mut self.face_cache) {
            if g <= 0.0 {
                continue;
            }
            let q = g * (self.pressure_of(a) - self.pressure_of(b));
            let (temperature, concentration) = if q > 0.0 {
                match a {
                    Side::Cell(i) => (self.temperature[i], self.dissolved[i] / self.pore(i)),
                    Side::Inlet => (self.inlet_temperature, 0.0),
                    // Nothing flows in through the outlet at a positive q.
                    Side::Outlet => (self.inlet_temperature, 0.0),
                }
            } else {
                match b {
                    Side::Cell(i) => (self.temperature[i], self.dissolved[i] / self.pore(i)),
                    Side::Outlet => (self.inlet_temperature, 0.0),
                    Side::Inlet => (self.inlet_temperature, 0.0),
                }
            };
            let heat = q * rho_c_w * (temperature - T_REF);
            let solute = q * concentration;
            match (a, b) {
                (Side::Cell(i), Side::Cell(j)) => {
                    d_heat[i] -= heat * h;
                    d_heat[j] += heat * h;
                    d_solute[i] -= solute * h;
                    d_solute[j] += solute * h;
                }
                (Side::Inlet, Side::Cell(j)) => {
                    d_heat[j] += heat * h;
                    d_solute[j] += solute * h;
                    in_enthalpy += heat * h;
                }
                (Side::Cell(i), Side::Outlet) => {
                    d_heat[i] -= heat * h;
                    d_solute[i] -= solute * h;
                    out_enthalpy += heat * h;
                    out_solute += solute * h;
                    out_volume += q * h;
                }
                _ => {}
            }
        }

        // --- conduction, over the same faces, both directions -------------------------------
        let (nx, ny, nz) = self.counts;
        let conduct = |a: usize, b: usize, d_heat: &mut Vec<f64>, s: &Puck| {
            let (ka, kb) = (s.conductivity(a), s.conductivity(b));
            if ka <= 0.0 || kb <= 0.0 {
                return;
            }
            let g = area / (0.5 * s.dx / ka + 0.5 * s.dx / kb);
            let flux = g * (s.temperature[a] - s.temperature[b]) * h;
            d_heat[a] -= flux;
            d_heat[b] += flux;
        };
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx - 1 {
                    conduct(
                        i + nx * (j + ny * k),
                        i + 1 + nx * (j + ny * k),
                        &mut d_heat,
                        self,
                    );
                }
            }
        }
        for k in 0..nz {
            for j in 0..ny - 1 {
                for i in 0..nx {
                    conduct(
                        i + nx * (j + ny * k),
                        i + nx * (j + 1 + ny * k),
                        &mut d_heat,
                        self,
                    );
                }
            }
        }
        for k in 0..nz - 1 {
            for j in 0..ny {
                for i in 0..nx {
                    conduct(
                        i + nx * (j + ny * k),
                        i + nx * (j + ny * (k + 1)),
                        &mut d_heat,
                        self,
                    );
                }
            }
        }

        // --- apply -------------------------------------------------------------------------
        for idx in 0..n {
            self.temperature[idx] += d_heat[idx] / self.heat_capacity(idx);
            self.dissolved[idx] += d_solute[idx];
        }

        // --- dissolution, against the liquid already in the pore ----------------------------
        //
        // Within a step the cell is closed — the fluxes above have already been applied — so the
        // linear driving force has a closed solution over the step and there is no reason to take
        // a first-order one:
        //
        //   dm/dt = −K(m − m₀·C/C_sat),   C = (S − m)/V_pore,   S = m + dissolved  (fixed)
        //         = −K((1+β)m − βS),      β = m₀/(V_pore·C_sat)
        //   m(h)  = m_eq + (m − m_eq)·e^{−K(1+β)h},   m_eq = βS/(1+β)
        //
        // Exact, unconditionally positive, and it contributes no stability limit. It also means a
        // stalled cell *stays* at its equilibrium rather than oscillating around it, which an
        // explicit first-order step would do the moment `K(1+β)h` passed 2.
        let c_sat = self.bed.saturation_concentration.to_si();
        for idx in 0..n {
            if !self.packed[idx] || self.initial_solids[idx] <= 0.0 {
                continue;
            }
            let pore = self.porosity[idx] * self.dx.powi(3);
            if pore <= 0.0 {
                continue;
            }
            let rate = self.grind[idx]
                .extraction_rate(&self.bed, Temperature::from_si(self.temperature[idx]));
            let beta = self.initial_solids[idx] / (pore * c_sat);
            let total = self.solids[idx] + self.dissolved[idx];
            let equilibrium = beta * total / (1.0 + beta);
            let m =
                equilibrium + (self.solids[idx] - equilibrium) * (-rate * (1.0 + beta) * h).exp();
            // Clamped only against the arithmetic, not against the physics: `m` is already
            // between its start and its equilibrium by construction.
            let m = m.clamp(0.0, total);
            self.dissolved[idx] += self.solids[idx] - m;
            self.solids[idx] = m;
        }

        self.delivered_volume += out_volume;
        self.delivered_solute += out_solute;
        self.delivered_enthalpy += out_enthalpy;
        self.admitted_enthalpy += in_enthalpy;
        Ok(())
    }

    /// The bed's holdings, arranged so that both entries are **constant**.
    ///
    /// Mass is solute: what is still in the particles, plus what is dissolved in the pores, plus
    /// what has left in the cup. Energy is the bed's enthalpy plus what has left minus what was
    /// admitted. Nothing crosses the bus, so a change in either is a leak in the discretisation
    /// and not a coupling — which is what makes [`Domain::books_balance`] worth claiming here.
    fn ledger(&self) -> Ledger {
        let mut enthalpy = 0.0;
        for idx in 0..self.temperature.len() {
            enthalpy += self.heat_capacity(idx) * (self.temperature[idx] - T_REF);
        }
        let solute: f64 = self.solids.iter().sum::<f64>() + self.dissolved.iter().sum::<f64>();
        Ledger::new()
            .with(quantity::MASS, solute + self.delivered_solute)
            .with(
                quantity::ENERGY,
                enthalpy + self.delivered_enthalpy - self.admitted_enthalpy,
            )
    }

    /// **Yes.** Every transfer is in flux form: what one cell loses another gains, exactly, and
    /// what crosses a boundary is added to the delivered totals in the same statement that removes
    /// it from a cell. The solve's residual moves *where* the liquid goes, not *how much* of it
    /// there is.
    fn books_balance(&self) -> bool {
        true
    }

    fn readings(&self) -> Vec<Reading> {
        vec![
            Reading::new(&self.name, "flow", self.flow_rate().to_si() * 1000.0, "g/s"),
            Reading::new(
                &self.name,
                "delivered",
                self.delivered().to_si() * 1000.0,
                "g",
            ),
            Reading::new(&self.name, "yield", self.yield_fraction() * 100.0, "%"),
            Reading::new(&self.name, "TDS", self.tds() * 100.0, "%"),
            Reading::new(
                &self.name,
                "bed temperature",
                self.bed_temperature().to_si() - 273.15,
                "C",
            ),
            Reading::new(
                &self.name,
                "outlet temperature",
                self.outlet_temperature().to_si() - 273.15,
                "C",
            ),
            Reading::new(&self.name, "unevenness", self.unevenness(), ""),
            // The one that says *channelling* rather than merely *uneven*. One for a good puck;
            // `unevenness` beside it is 0.10 on a good puck and is mostly the axial gradient.
            Reading::new(&self.name, "ring over core", self.radial_contrast(), ""),
        ]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    /// **And mutably**, because a domain that can be read and not written is one a coupling can
    /// only fail at silently. `Simulation::domain_as_mut` returns `None` when this is not
    /// implemented, and a caller that wrote through it would do nothing and report nothing —
    /// which is how a whole coupled body came back reporting zero strain that read as *no stress*
    /// rather than as *not connected*.
    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }

    fn as_field(&self) -> Option<&dyn ScalarField> {
        Some(self)
    }

    fn checkpoint(&mut self) {
        self.saved = Some(Box::new(Saved {
            pressure: self.pressure.clone(),
            temperature: self.temperature.clone(),
            solids: self.solids.clone(),
            dissolved: self.dissolved.clone(),
            delivered_volume: self.delivered_volume,
            delivered_solute: self.delivered_solute,
            delivered_enthalpy: self.delivered_enthalpy,
            admitted_enthalpy: self.admitted_enthalpy,
        }));
    }

    fn restore(&mut self) {
        if let Some(s) = self.saved.take() {
            self.pressure = s.pressure;
            self.temperature = s.temperature;
            self.solids = s.solids;
            self.dissolved = s.dissolved;
            self.delivered_volume = s.delivered_volume;
            self.delivered_solute = s.delivered_solute;
            self.delivered_enthalpy = s.delivered_enthalpy;
            self.admitted_enthalpy = s.admitted_enthalpy;
            self.saved = None;
        }
    }
}

impl Puck {
    fn pore(&self, idx: usize) -> f64 {
        let v = self.porosity[idx] * self.dx.powi(3);
        if v > 0.0 {
            v
        } else {
            f64::INFINITY
        }
    }
}

/// The **extraction** field, which is the one [`Domain::as_field`] nominates.
///
/// Extraction rather than temperature, and the choice matters because a viewer that is handed one
/// field draws that one and nothing says what it missed. A bed under flow is *isothermal* unless
/// somebody deliberately cooled the basket, so a temperature panel is a flat rectangle on every
/// ordinary run — a picture that renders, looks fine, and carries no information. Extraction is
/// never flat: it has the axial gradient on a good puck and the radial one on a bad puck, which
/// are the two things anybody opens the picture to see.
///
/// The other four are [`Puck::field`] away, and `dualis_scene::sample_field` turns any of them
/// into a panel.
impl ScalarField for Puck {
    fn unit(&self) -> &'static str {
        ""
    }

    fn at(&self, p: LengthVec, t: Time) -> f64 {
        self.field(Observable::Extraction).at(p, t)
    }
}

/// One of a [`Puck`]'s fields, borrowed for sampling.
#[derive(Clone, Copy, Debug)]
pub struct PuckField<'a> {
    puck: &'a Puck,
    what: Observable,
}

impl ScalarField for PuckField<'_> {
    fn unit(&self) -> &'static str {
        match self.what {
            Observable::Temperature => "K",
            Observable::Pressure => "Pa",
            Observable::Speed => "m/s",
            Observable::Extraction => "",
            Observable::Concentration => "kg/m3",
        }
    }

    /// Nearest cell, not trilinear.
    ///
    /// Deliberate, and the reason is the basket wall. Interpolating across the boundary between a
    /// packed cell and a wall cell would produce a smooth ramp of "half-extracted steel" — values
    /// at points where the quantity does not exist. A viewer would render that as a soft edge and
    /// it would look better than the truth.
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let (nx, ny, nz) = self.puck.counts;
        let q = p.to_si() / self.puck.dx;
        if q.is_nan() {
            return 0.0;
        }
        let pick = |v: f64, n: usize| -> usize { (v.floor().max(0.0) as usize).min(n - 1) };
        let (i, j, k) = (pick(q.x, nx), pick(q.y, ny), pick(q.z, nz));
        match self.what {
            Observable::Temperature => self.puck.temperature_at(i, j, k).to_si(),
            Observable::Pressure => self.puck.pressure_at(i, j, k).to_si(),
            Observable::Speed => self.puck.pore_velocity_at(i, j, k).length(),
            Observable::Extraction => self.puck.extraction_at(i, j, k),
            Observable::Concentration => self.puck.concentration_at(i, j, k).to_si(),
        }
    }
}

/// Convenience: what a shot came out as, for a caller that wants the answer rather than the run.
#[derive(Clone, Copy, Debug)]
pub struct Shot {
    /// Beverage mass in the cup.
    pub beverage: Mass,
    /// Solute in the cup.
    pub solute: Mass,
    /// Solute over dry dose.
    pub yield_fraction: f64,
    /// Solute over beverage mass.
    pub tds: f64,
    /// How long it ran.
    pub duration: Time,
    /// Coefficient of variation of the per-cell extraction. Zero is even.
    pub unevenness: f64,
    /// Mean temperature of the liquid leaving the basket over the whole shot.
    pub outlet_temperature: Temperature,
}

impl Puck {
    /// What is in the cup so far.
    pub fn shot(&self, elapsed: Time) -> Shot {
        Shot {
            beverage: self.delivered(),
            solute: self.delivered_solute(),
            yield_fraction: self.yield_fraction(),
            tds: self.tds(),
            duration: elapsed,
            unevenness: self.unevenness(),
            outlet_temperature: Temperature::from_si(if self.delivered_volume > 0.0 {
                self.delivered_enthalpy / (self.delivered_volume * self.liquid.rho_c()) + T_REF
            } else {
                self.inlet_temperature
            }),
        }
    }

    /// The power the liquid is carrying out of the basket right now.
    pub fn thermal_power_out(&self) -> Power {
        Power::from_si(
            self.flow * self.liquid.rho_c() * (self.outlet_temperature().to_si() - T_REF),
        )
    }

    /// Peak pore speed anywhere in the bed — the number a channel shows up in.
    pub fn peak_pore_speed(&self) -> Velocity {
        let (nx, ny, nz) = self.counts;
        let mut worst: f64 = 0.0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    worst = worst.max(self.pore_velocity_at(i, j, k).length());
                }
            }
        }
        Velocity::from_si(worst)
    }

    /// Total extractable solute the basket started with.
    pub fn extractable(&self) -> Mass {
        Mass::from_si(self.extractable)
    }

    /// Energy delivered to the cup so far, measured above 273.15 K.
    pub fn delivered_enthalpy(&self) -> Energy {
        Energy::from_si(self.delivered_enthalpy)
    }
}
