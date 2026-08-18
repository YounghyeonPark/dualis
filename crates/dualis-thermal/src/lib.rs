//! dualis-thermal: heat, as a domain built on the `dualis-core` kernel.
//!
//! Four domains, spanning three dimensionalities and one graph. The first two sit either side of
//! the line that matters to a scheduler:
//!
//! - [`LumpedMass`] has one temperature and no internal structure. Its stability
//!   limit is its own time constant, which is seconds for a piece of glass in still
//!   air — so it takes one step per frame and costs nothing.
//! - [`Bar1D`] resolves a gradient on a grid, and pays the explicit diffusion limit
//!   `dt < dx²/2α` for it. On a millimetre grid that is about a second in N-BK7 and
//!   seven milliseconds in aluminium, and *that* two-orders-of-magnitude gap between
//!   two parts of the same instrument is what
//!   [`Schedule::Multirate`](dualis_core::Schedule::Multirate) exists for.
//!
//! The third answers a different question. Both of the above report *one* body:
//!
//! - [`ThermalNetwork`] is n lumped bodies joined by conductances — winding, stator,
//!   housing — and it carries the **drop across each joint**, which is the number a
//!   designer actually needs and the one a single lumped mass cannot give: it reports
//!   the temperature of the skin and the winding together. It also expresses a *contact*
//!   resistance between different materials, which [`Bar1D`]'s uniform grid cannot.
//!   A network of one node reduces to a [`LumpedMass`] bit for bit, so it inherits every
//!   check that domain already passes.
//!
//! And the fourth resolves what a bar cannot:
//!
//! - [`Solid3D`] is conduction in three dimensions on a cubic grid, which is what a hot spot
//!   needs: heat spreading *sideways* out of a spot is the whole job of a spreader plate and a
//!   fin, and a one-dimensional model has nowhere for it to go but along. It pays `dx²/6α`,
//!   a third of [`Bar1D`]'s limit, because the explicit limit tightens with every axis.
//!
//! # Where the heat comes from
//!
//! Neither domain generates any. They consume [`HEAT`] from the kernel's
//! [`Exchange`], and the thing that publishes it is optics:
//! `SurfaceOptics::absorptance` against a `SpectralPower` is a definite number of
//! watts, and those watts have to go somewhere. That is the whole coupling, and it
//! is auditable because it goes over the bus.
//!
//! # What is deliberately simple
//!
//! Both domains are explicit and both are linear in temperature except for the
//! radiative term. There is no implicit solver, no mesh, no natural convection
//! model — a convective loss is a coefficient the caller supplies, because
//! computing one honestly means solving a fluid problem and that is a different
//! crate. The point here is a domain that couples correctly and reports its own
//! stability limit, not a competitive thermal solver.

// Every public item carries a doc comment. Denied rather than warned: a public physics API
// whose `Length::mm` shows a blank summary in rustdoc is documented in the sense that a
// paragraph exists somewhere, and not in the sense a reader needs.
#![deny(missing_docs)]

pub mod network;
pub mod solid;

use dualis_core::conserved::quantity;
use dualis_core::{
    Domain, Exchange, Interface, Kind, Ledger, Reading, ScalarField, Substance, Violation,
};
use dualis_units::{
    Area, Energy, HeatCapacity, Length, LengthVec, Power, Temperature, Time, Volume,
    STEFAN_BOLTZMANN,
};
use glam::DVec3;
pub use network::{Node, SteadyState, ThermalNetwork};
pub use solid::{Face, GapPatch, Solid3D, STABLE_FOURIER_3D};

/// The bus channel heat arrives on, in joules.
///
/// Joules rather than watts, because a domain steps over an interval and what
/// crossed the interface is an amount, not a rate. The publisher multiplies by its
/// own `dt`, which is what makes the audit an equality rather than an approximation.
pub const HEAT: &str = quantity::ENERGY;

/// How a body loses heat to its surroundings.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Environment {
    /// The temperature the surroundings sit at, and what a body relaxes towards.
    pub ambient: Temperature,
    /// Convective coefficient `h`, W·m⁻²·K⁻¹. Still air is about 5 to 10; forced air
    /// 25 to 100; water hundreds. There is no right default, so there is no default.
    pub convection_w_per_m2_k: f64,
    /// Surface area available to lose heat through.
    pub area: Area,
}

impl Environment {
    /// A part sitting in still room air.
    pub fn still_air(ambient: Temperature, area: Area) -> Environment {
        Environment {
            ambient,
            convection_w_per_m2_k: 7.0,
            area,
        }
    }

    /// Power lost from a surface at `temperature`, convective plus radiative.
    ///
    /// The radiative term is `εσA(T⁴ - T_a⁴)`, and it is the only non-linearity in
    /// this crate. It is also not negligible: at room temperature a black surface
    /// radiates about 6 W·m⁻²·K⁻¹, which is the same order as still-air convection,
    /// so leaving it out halves the loss.
    pub fn loss_from(&self, temperature: Temperature, emissivity: f64) -> Power {
        let (t, ta) = (temperature.to_si(), self.ambient.to_si());
        let convective = self.convection_w_per_m2_k * self.area.to_si() * (t - ta);
        let radiative =
            emissivity * STEFAN_BOLTZMANN.to_si() * self.area.to_si() * (t.powi(4) - ta.powi(4));
        Power::from_si(convective + radiative)
    }
}

/// A body at one temperature.
///
/// The lumped approximation: valid when heat spreads through the body faster than it
/// escapes, which is the small-Biot-number condition `hL/k << 1`. For a 10 mm piece
/// of N-BK7 in still air that number is about 0.03, so it holds comfortably; for the
/// same piece in flowing water it does not, and [`Bar1D`] is the honest choice.
/// [`LumpedMass::biot_number`] says which situation you are in.
pub struct LumpedMass {
    name: String,
    substance: Substance,
    volume: Volume,
    /// Characteristic length for the Biot number — volume over surface area.
    thickness: Length,
    temperature: Temperature,
    environment: Environment,
    saved: Option<(Temperature, f64, f64)>,
    /// Joules taken from the bus over the run, for the books.
    absorbed: f64,
    /// Joules given up to the environment over the run.
    lost: f64,
}

impl LumpedMass {
    /// A body of one substance at one temperature, losing heat to its surroundings.
    ///
    /// `thickness` is the characteristic length for the Biot number, usually volume over
    /// surface area. It does not enter the dynamics — only
    /// [`LumpedMass::biot_number`], which is how you find out whether the lumped
    /// approximation was honest here.
    pub fn new(
        name: impl Into<String>,
        substance: Substance,
        volume: Volume,
        thickness: Length,
        initial: Temperature,
        environment: Environment,
    ) -> LumpedMass {
        LumpedMass {
            name: name.into(),
            substance,
            volume,
            thickness,
            temperature: initial,
            environment,
            saved: None,
            absorbed: 0.0,
            lost: 0.0,
        }
    }

    /// Absolute temperature of the whole body.
    pub fn temperature(&self) -> Temperature {
        self.temperature
    }

    /// Rise above ambient.
    pub fn rise(&self) -> Temperature {
        self.temperature - self.environment.ambient
    }

    /// `mc_p` for this body. Infinite if the substance has no specific heat recorded,
    /// which makes it refuse to warm rather than warm by a made-up amount.
    pub fn heat_capacity(&self) -> HeatCapacity {
        self.substance
            .heat_capacity(self.volume)
            .unwrap_or(HeatCapacity::from_si(f64::INFINITY))
    }

    /// Heat taken from the bus over the run.
    pub fn absorbed_energy(&self) -> Energy {
        Energy::from_si(self.absorbed)
    }

    /// Heat given up to the environment over the run.
    pub fn lost_energy(&self) -> Energy {
        Energy::from_si(self.lost)
    }

    /// `hL/k` — whether the lumped approximation is honest here. Under about 0.1 it
    /// is; well past that, the body has an internal gradient this domain cannot see.
    pub fn biot_number(&self) -> f64 {
        let Some(thermal) = self.substance.thermal else {
            return f64::INFINITY;
        };
        self.environment.convection_w_per_m2_k * self.thickness.to_si()
            / thermal.conductivity.to_si()
    }

    /// Time constant `C/(hA)` — how long it takes to get most of the way to
    /// equilibrium, and the step this domain must not much exceed.
    /// How long it takes to settle, linearised **at the temperature it is at now**.
    ///
    /// `C / (hA + 4εσA·T³)`. Both loss paths, because the crate's own
    /// [`Environment::loss_from`] says why: at room temperature a black surface radiates about
    /// 6 W·m⁻²·K⁻¹, the same order as still-air convection, so leaving it out roughly halves
    /// the conductance and doubles this number.
    ///
    /// It used to leave it out, and reported the same time constant for a polished surface and
    /// a blackbody one. Measured on a 1.12 kg box in still air under 21 W, against the time to
    /// reach 63% of its settled rise:
    ///
    /// ```text
    ///   emissivity    was    at rest    once hot    measured
    ///         0.05   53.0      50.8        48.3      49.2 min
    ///         0.09   53.0      49.2        45.4      46.6
    ///         0.50   53.0      37.1        30.1      32.2
    ///         0.90   53.0      29.9        23.8      25.6
    /// ```
    ///
    /// **Why the temperature it is at now, and not ambient.** The radiative conductance grows
    /// as `T³`, so a body running hot settles faster than the same body at rest. This is not a
    /// constant of the body; it is a property of its current state, and a function returning
    /// one number has to say which one. The two right-hand columns above are the same call on
    /// the same body cold and settled — and the measured figure lies *between* them, because a
    /// large-signal time constant is an average over the trajectory and the small-signal one
    /// bounds it at each end.
    ///
    /// [`LumpedMass::max_stable_dt`] re-reads this every step, so a scheduler tightens as the
    /// body warms: 179 s to 143 s over that last row, against 318 s before radiation was
    /// counted at all.
    ///
    /// Infinite when nothing carries heat away — no convection and no emissivity — which is a
    /// body that never settles rather than one that settles instantly.
    pub fn time_constant(&self) -> Time {
        let capacity = self.heat_capacity().to_si();
        let conductance = self.loss_conductance(self.temperature);
        if conductance <= 0.0 || !capacity.is_finite() {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(capacity / conductance)
    }

    /// `d(loss)/dT` at a temperature: convection plus the linearised radiative term.
    fn loss_conductance(&self, at: Temperature) -> f64 {
        linearised_loss_conductance(&self.environment, at, self.emissivity())
    }

    /// Steady-state rise for a constant absorbed power, **with radiation**.
    ///
    /// Solves `P = hA·ΔT + εσA((Tₐ+ΔT)⁴ − Tₐ⁴)` rather than `P/(hA)`. The answer a designer
    /// actually wants: not how it gets there, but where it ends up — and the radiative term is
    /// not a correction to it. On a 1.12 kg box in still air under 21 W, `P/(hA)` says 99.2 K
    /// whatever the surface is, against a measured 92.9 K at ε = 0.05 and **47.5 K at ε = 1.0**.
    /// It was over by 2.09× at the top of that range, and the error is always in the
    /// comfortable direction.
    ///
    /// One positive root, because the right-hand side is strictly increasing in `ΔT` above
    /// `−Tₐ`. Newton from the convective guess, which is an overestimate and therefore
    /// approaches from the side where the derivative is largest — three or four steps.
    ///
    /// Infinite only when nothing carries heat away at all.
    pub fn equilibrium_rise(&self, absorbed: Power) -> Temperature {
        let p = absorbed.to_si();
        let area = self.environment.area.to_si();
        let ha = self.environment.convection_w_per_m2_k * area;
        let er = self.emissivity() * STEFAN_BOLTZMANN.to_si() * area;
        let ta = self.environment.ambient.to_si();

        if !p.is_finite() || (ha <= 0.0 && er <= 0.0) {
            return Temperature::from_si(if p == 0.0 { 0.0 } else { f64::INFINITY });
        }
        if p == 0.0 {
            return Temperature::from_si(0.0);
        }
        if er <= 0.0 {
            return Temperature::from_si(p / ha);
        }

        // Start from whichever single-path answer exists; both are above the true root when
        // the other path is also carrying heat, so Newton descends onto it.
        let mut x = if ha > 0.0 {
            p / ha
        } else {
            (p / er + ta.powi(4)).max(0.0).powf(0.25) - ta
        };
        for _ in 0..64 {
            let t = (ta + x).max(0.0);
            let f = ha * x + er * (t.powi(4) - ta.powi(4)) - p;
            let df = ha + 4.0 * er * t * t * t;
            if df <= 0.0 || !f.is_finite() {
                break;
            }
            let step = f / df;
            x -= step;
            if step.abs() <= 1e-12 * (1.0 + x.abs()) {
                break;
            }
        }
        Temperature::from_si(x)
    }

    fn emissivity(&self) -> f64 {
        self.substance.thermal.map(|t| t.emissivity).unwrap_or(0.0)
    }
}

impl Domain for LumpedMass {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// A tenth of the time constant. Explicit Euler on `dT/dt = -(T-Ta)/τ` is stable
    /// up to `2τ` and *accurate* nowhere near it, so the limit reported is the
    /// accuracy one — a scheduler that honours it gets a curve rather than a
    /// staircase.
    fn max_stable_dt(&self, _now: Time) -> Time {
        self.time_constant() / 10.0
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let capacity = self.heat_capacity();
        if !capacity.to_si().is_finite() || capacity.to_si() <= 0.0 {
            return Err(Violation::at(
                &self.name,
                "substance has no heat capacity",
                capacity.to_si(),
            ));
        }

        // This step's share of the channel, not all of it. A lumped mass is subcycled under
        // `Schedule::Multirate`, and taking the whole outer step's joules on the first substep
        // deposits them all at its beginning — which stops the substep count from improving
        // anything. See `Exchange::take_share`.
        let gained = bus.take_share(HEAT, dt);
        self.absorbed += gained;

        let lost = self
            .environment
            .loss_from(self.temperature, self.emissivity());
        let lost_joules = lost.to_si() * dt.to_si();
        self.lost += lost_joules;

        let net = gained - lost_joules;
        self.temperature += Temperature::from_si(net / capacity.to_si());
        Ok(())
    }

    /// Energy this body accounts for: what is stored above ambient, plus what has
    /// already left to the environment.
    ///
    /// Not minus what it absorbed. `stored + lost` *is* what it absorbed, so
    /// subtracting that as well would cancel the entry against itself and leave the
    /// publisher's debt unmatched — the audit catches that immediately, which is how
    /// this convention got settled.
    fn ledger(&self) -> Ledger {
        let stored = self.heat_capacity().to_si() * self.rise().to_si();
        Ledger::new().with(quantity::ENERGY, stored + self.lost)
    }

    /// Everything the ledger reads, not only the temperature.
    ///
    /// `ledger()` is `stored + lost`, and `stored` follows the temperature while `lost` is a
    /// running total. Saving one and not the other means a sweep that gets rewound leaves its
    /// losses behind: under `Schedule::Iterative` the books then grow by one sweep of shed heat
    /// per iteration and the audit sees energy created out of nothing. Measured at 1567.6 J
    /// becoming 1600.9 J over forty advances of three sweeps each.
    ///
    /// It went unnoticed because nothing in this workspace has a residual, so `iterate` always
    /// converged on its first sweep and the restore branch never ran —
    /// `crates/dualis/tests/iterative_restore.rs` supplies the domain that makes it run.
    fn checkpoint(&mut self) {
        self.saved = Some((self.temperature, self.absorbed, self.lost));
    }

    fn restore(&mut self) {
        if let Some((t, absorbed, lost)) = self.saved {
            self.temperature = t;
            self.absorbed = absorbed;
            self.lost = lost;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    /// One temperature, which is the whole of what a lumped model claims to know.
    fn readings(&self) -> Vec<Reading> {
        vec![
            Reading::new(
                &self.name,
                "temperature",
                self.temperature.to_si() - 273.15,
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
}

/// One-dimensional explicit heat conduction on a uniform grid.
///
/// Exists to exercise the thing [`LumpedMass`] cannot: a real stability limit that a
/// scheduler has to subcycle around. The explicit update
/// `T'ᵢ = Tᵢ + α dt/dx² (Tᵢ₊₁ - 2Tᵢ + Tᵢ₋₁)` is stable only for
/// `α dt/dx² ≤ 1/2`, and exceeding it does not degrade the answer gracefully — it
/// oscillates and diverges within a few steps.
///
/// Ends are insulated, so the total heat is conserved exactly and the audit has
/// something sharp to check.
pub struct Bar1D {
    name: String,
    substance: Substance,
    /// Temperatures at the cell centres.
    cells: Vec<f64>,
    saved: Vec<f64>,
    dx: Length,
    /// Cross-sectional area, for turning joules into a temperature.
    area: Area,
    /// The boundary this bar offers to other domains, one face per cell, if it has one.
    boundary: Option<Interface>,
    absorbed: f64,
    /// The temperature the stored heat is measured from — see [`Bar1D::stored_heat`].
    reference: f64,
}

impl Bar1D {
    /// A bar of `cells` cells, each `dx` long, all starting at `initial`.
    pub fn new(
        name: impl Into<String>,
        substance: Substance,
        cells: usize,
        dx: Length,
        area: Area,
        initial: Temperature,
    ) -> Bar1D {
        let cells = cells.max(2);
        let temps = vec![initial.to_si(); cells];
        Bar1D {
            name: name.into(),
            substance,
            cells: temps.clone(),
            saved: temps,
            dx,
            area,
            boundary: None,
            absorbed: 0.0,
            reference: initial.to_si(),
        }
    }

    /// Expose the bar's long side as an [`Interface`], one face per cell.
    ///
    /// Without this the bar can only be heated lumpedly, and the heat lands in cell 0
    /// because that is where a surface absorbing light *would* put it if the bar knew where
    /// the surface was. It does not: `Exchange::publish` carries an amount and no place, so
    /// "the light hit the middle" is unsayable.
    ///
    /// With it, whoever illuminates the bar publishes a [`Flux`](dualis_core::Flux) over these faces and the
    /// heat appears where it landed. One face per cell deliberately: the two sides then
    /// share a discretisation, so nothing has to interpolate, and interpolation is where a
    /// coupling loses energy. A publisher on a different grid resamples explicitly with
    /// [`Flux::resample`](dualis_core::Flux::resample), which the bus insists on rather than doing quietly.
    ///
    /// `face_area` is the area of one cell's exposed side, which is not the bar's
    /// cross-section — a bar conducts along its length and is illuminated across it. It is
    /// only used to turn a lumped total into a distribution, so it does not enter the
    /// conduction at all.
    pub fn exposing(mut self, boundary: impl Into<String>, face_area: Area) -> Bar1D {
        self.boundary = Some(Interface::uniform(boundary, self.cells.len(), face_area));
        self
    }

    /// The boundary other domains publish onto, if [`Bar1D::exposing`] gave it one.
    pub fn boundary(&self) -> Option<&Interface> {
        self.boundary.as_ref()
    }

    /// Temperature of one cell, clamped to the ends of the bar.
    pub fn temperature_at(&self, index: usize) -> Temperature {
        Temperature::from_si(self.cells[index.min(self.cells.len() - 1)])
    }

    /// How many cells the bar is cut into.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Mean temperature along the bar.
    pub fn mean_temperature(&self) -> Temperature {
        Temperature::from_si(self.cells.iter().sum::<f64>() / self.cells.len() as f64)
    }

    /// Temperature difference between the ends — what a gradient looks like from
    /// outside, and what a lumped model reports as zero.
    pub fn end_to_end(&self) -> Temperature {
        Temperature::from_si(self.cells[self.cells.len() - 1] - self.cells[0])
    }

    fn cell_capacity(&self) -> f64 {
        let volume = Volume::from_si(self.area.to_si() * self.dx.to_si());
        self.substance
            .heat_capacity(volume)
            .map(|c| c.to_si())
            .unwrap_or(f64::INFINITY)
    }

    /// Heat held, measured from the temperature the bar started at.
    ///
    /// The reference point of an enthalpy is arbitrary, so it should be chosen for
    /// precision, and the natural-looking choice is the bad one. Against absolute zero a
    /// 20 mm aluminium bar of 1 cm² section holds 1.42 kJ, and a millijoule arriving is a
    /// change in the seventh significant figure — so differencing two such numbers leaves a
    /// rounding floor of a few times 10⁻¹² J whatever the transfer was, and the audit's
    /// *relative* check on a 1 mJ step is then asking for precision the arithmetic threw
    /// away. Refining the grid makes it worse rather than better: 1.6×10⁻¹² J at 41 cells
    /// against 7.3×10⁻¹² J at 161, because there are more absolute temperatures to add up.
    ///
    /// Measured from the initial temperature, the number being summed *is* the change, and
    /// the audit's precision tracks the heat that moved rather than the enthalpy it moved
    /// within.
    fn stored_heat(&self) -> f64 {
        self.cell_capacity() * self.cells.iter().map(|t| t - self.reference).sum::<f64>()
    }

    /// Heat taken from the bus over the run.
    pub fn absorbed_energy(&self) -> Energy {
        Energy::from_si(self.absorbed)
    }

    /// `α dt/dx²`, the number that must stay at or under 1/2.
    pub fn fourier_number(&self, dt: Time) -> f64 {
        let Some(alpha) = self.substance.diffusivity() else {
            return f64::INFINITY;
        };
        alpha.to_si() * dt.to_si() / (self.dx.to_si() * self.dx.to_si())
    }
}

impl Domain for Bar1D {
    fn books_balance(&self) -> bool {
        true
    }

    fn name(&self) -> &str {
        &self.name
    }

    /// `dx²/(2α)` — the explicit diffusion limit, exactly.
    fn max_stable_dt(&self, _now: Time) -> Time {
        let Some(alpha) = self.substance.diffusivity() else {
            return Time::from_si(f64::INFINITY);
        };
        Time::from_si(self.dx.to_si() * self.dx.to_si() / (2.0 * alpha.to_si()))
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let f = self.fourier_number(dt);
        if !f.is_finite() {
            return Err(Violation::at(&self.name, "substance has no diffusivity", f));
        }
        // Refuse rather than diverge. A scheduler honouring `max_stable_dt` never
        // sees this; one that ignores it gets told which limit it broke and by how
        // much, instead of a bar full of oscillating nonsense.
        if f > 0.5 + 1e-12 {
            return Err(Violation {
                quantity: "Fourier number".to_string(),
                site: format!("{} (explicit conduction)", self.name),
                before: 0.5,
                after: f,
                scale: 0.5,
                tolerance: 1e-12,
            });
        }

        let capacity = self.cell_capacity();

        // Lumped heat has no place, so it goes into the first cell — where a surface
        // absorbing light would put it, if the bar knew where the surface was. This step's
        // share of it: the bar subcycles hard under `Schedule::Multirate`, and taking the whole
        // interval at once would put every joule of it in the first substep.
        let gained = bus.take_share(HEAT, dt);
        self.absorbed += gained;
        self.cells[0] += gained / capacity;

        // Heat that does know where it landed. Taken before the conduction sweep so it
        // spreads on the same step it arrives, and taken out of the borrow of `boundary`
        // before the cells are written to.
        let arriving = match self.boundary.as_ref() {
            Some(boundary) => Some(bus.take_on(boundary, HEAT)?),
            None => None,
        };
        if let Some(flux) = arriving {
            for (cell, joules) in self.cells.iter_mut().zip(flux.per_face()) {
                *cell += joules / capacity;
            }
            self.absorbed += flux.total();
        }

        // Insulated ends: the boundary cell exchanges with its one neighbour only,
        // which is what makes the total conserved to the last bit.
        let previous = self.cells.clone();
        let last = previous.len() - 1;
        for i in 0..=last {
            let left = previous[i.saturating_sub(1)];
            let right = previous[(i + 1).min(last)];
            self.cells[i] = previous[i] + f * (left - 2.0 * previous[i] + right);
        }
        Ok(())
    }

    /// Heat gained since the start. The ends are insulated, so this is exactly what
    /// came in over the bus — see the note on [`LumpedMass::ledger`] for why the
    /// absorbed total is not subtracted here as well.
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.stored_heat())
    }

    fn checkpoint(&mut self) {
        self.saved = self.cells.clone();
    }

    fn restore(&mut self) {
        self.cells = self.saved.clone();
    }

    fn supports_restore(&self) -> bool {
        true
    }

    /// Mean and peak in celsius, and what it has absorbed.
    ///
    /// Both ends of the profile, because a bar's whole reason to exist rather than a lumped mass
    /// is that those two differ — reporting the mean alone would describe it as the thing it is
    /// not.
    fn readings(&self) -> Vec<Reading> {
        let peak = (0..self.cells.len())
            .map(|i| self.temperature_at(i).to_si())
            .fold(f64::MIN, f64::max);
        vec![
            Reading::new(
                &self.name,
                "mean",
                self.mean_temperature().to_si() - 273.15,
                "C",
            ),
            Reading::new(&self.name, "peak", peak - 273.15, "C"),
            Reading::new(&self.name, "absorbed", self.absorbed_energy().to_si(), "J"),
        ]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    /// The bar reads as a temperature field, so a renderer never has to know it is a bar.
    fn as_field(&self) -> Option<&dyn dualis_core::ScalarField> {
        Some(self)
    }
}

/// The bar as a temperature field, in kelvin.
///
/// The first implementation of [`ScalarField`] in the workspace, and it exists to answer a
/// question rather than to be used internally: the trait was written as the interface a
/// visualiser would read a simulation through, and until something implemented it, whether
/// it was the *right* interface was a guess. Two things came out of implementing it.
///
/// # The bar lies along x, and is uniform in y and z
///
/// Cell `i` is centred at `(i + ½)·dx`, so the bar occupies `0` to `n·dx`. Off the ends the
/// value is held constant, which is not a fudge: the ends are insulated, so the temperature
/// really does stop changing there. Off-axis it is uniform, which is not an approximation
/// being hidden either — a one-dimensional model *is* the claim that nothing varies across
/// the bar, and [`LumpedMass::biot_number`] is where you check whether that claim holds.
///
/// # `at` ignores the time it is given, and that is the interface's one rough edge
///
/// [`ScalarField::at`] takes a [`Time`], because a closed-form field like
/// [`Motion`](dualis_core::Motion) can answer for any instant. A marched domain cannot: it
/// holds *now* and nothing else. So the argument is ignored here, and a caller wanting a
/// different instant has to have recorded one.
///
/// The default [`ScalarField::rate`] would then read zero — it differences `at` across two
/// times — which would be wrong rather than merely unavailable, since the bar is visibly
/// heating. It is overridden below, and the fix is not a workaround: a diffusive field's
/// time derivative *is* `α∇²T`, so the governing equation supplies from the present state
/// exactly what the finite difference wanted history for.
///
/// # The derivatives use the domain's own stencil
///
/// [`ScalarField`] offers central differences over a step you choose, and its documentation
/// says a field that knows better should override them. This one does: `gradient` and
/// `laplacian` use the same mirrored three-point stencil that [`Domain::step`] integrates,
/// evaluated on the cells rather than by re-sampling the interpolated field. So the field
/// reports what the domain actually believes, and the `h` argument is ignored — asking for
/// a derivative on a scale finer than `dx` is asking for information the bar does not have.
///
/// They have to come from the cells rather than from `at`, and the reason is worth stating:
/// `at` interpolates linearly, so its exact second derivative is zero between nodes and
/// infinite at them. A Laplacian read off the interpolant would be useless. The cost is that
/// `gradient` is not quite the derivative of `at` — it is the derivative the *scheme* uses —
/// and that is an unavoidable property of sampling a discrete field, not a rough edge that
/// could be polished out.
impl ScalarField for Bar1D {
    /// **Kelvin**, because that is what the cells hold.
    ///
    /// Not celsius. `readings` reports celsius and a picture of a bar usually wants celsius, but
    /// both of those are *conversions a view chooses*; the field returns what it stores. Labelling
    /// this "C" would have put 293.15 under a degrees-celsius header, which is the failure a unit
    /// on a legend exists to prevent.
    fn unit(&self) -> &'static str {
        "K"
    }

    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let last = self.cells.len() - 1;
        // Position in cell-index space: cell centres land on the integers.
        let u = p.to_si().x / self.dx.to_si() - 0.5;
        // The NaN case is spelled out rather than folded into a negated comparison: it is a
        // real input a visualiser can hand over, and it must not reach the cast below.
        if u.is_nan() || u <= 0.0 {
            return self.cells[0];
        }
        if u >= last as f64 {
            return self.cells[last];
        }
        let i = u.floor() as usize;
        let f = u - i as f64;
        self.cells[i] * (1.0 - f) + self.cells[i + 1] * f
    }

    fn gradient(&self, p: LengthVec, _t: Time, _h: Length) -> DVec3 {
        let (left, _, right) = self.stencil_at(p);
        DVec3::new((right - left) / (2.0 * self.dx.to_si()), 0.0, 0.0)
    }

    fn laplacian(&self, p: LengthVec, _t: Time, _h: Length) -> f64 {
        let (left, centre, right) = self.stencil_at(p);
        let dx = self.dx.to_si();
        (left - 2.0 * centre + right) / (dx * dx)
    }

    /// `∂T/∂t = α∇²T` — the heat equation, evaluated rather than differenced.
    ///
    /// Conduction only. Heat arriving over the bus is a source term the field cannot see,
    /// so during illumination this reports how fast the bar is *spreading* what it has, not
    /// how fast it is warming. Away from the beam those are the same number.
    fn rate(&self, p: LengthVec, t: Time, _dt: Time) -> f64 {
        let Some(alpha) = self.substance.diffusivity() else {
            return 0.0;
        };
        alpha.to_si() * self.laplacian(p, t, self.dx)
    }
}

impl Bar1D {
    /// The three cell values the domain's own update uses at this point, with the ends
    /// mirrored exactly as [`Domain::step`] mirrors them.
    ///
    /// Past either end the stencil goes flat, so the derivatives agree with [`ScalarField::at`]
    /// holding its value out there. Inside, mirroring is what makes the ends insulated —
    /// but note what that does *not* mean. On a cell-centred grid the first sample sits half
    /// a cell inside the wall, so the gradient reported at `x = 0` is the mirrored estimate
    /// `(T₁ − T₀)/2dx` and not zero. Insulation shows up as no heat crossing the boundary,
    /// which the conservation audit checks; it does not show up as a zero slope at a
    /// position the grid cannot sample.
    fn stencil_at(&self, p: LengthVec) -> (f64, f64, f64) {
        let last = self.cells.len() - 1;
        let x = p.to_si().x / self.dx.to_si();
        // The bar occupies [0, L], so the walls are inside it and only points beyond them go
        // flat. NaN spelled out for the same reason as in `at`.
        if x.is_nan() || x < 0.0 {
            let v = self.cells[0];
            return (v, v, v);
        }
        if x >= self.cells.len() as f64 {
            let v = self.cells[last];
            return (v, v, v);
        }
        let i = (x as usize).min(last);
        (
            self.cells[i.saturating_sub(1)],
            self.cells[i],
            self.cells[(i + 1).min(last)],
        )
    }
}

/// `d(loss)/dT` for an environment at a temperature: convection plus the linearised radiative
/// term `4εσAT³`.
///
/// Shared by [`LumpedMass::time_constant`] and [`ThermalNetwork`] rather than written out twice,
/// so the two cannot drift apart. A network of one node has to reduce to a lumped mass *exactly*
/// — [`one_node_is_a_lumped_mass_bit_for_bit`] compares the two bit for bit over a whole
/// trajectory including the step limit — and a second copy of this expression is the obvious way
/// for that to quietly stop being true.
///
/// [`one_node_is_a_lumped_mass_bit_for_bit`]: https://github.com/YounghyeonPark/dualis-core/blob/main/crates/dualis-thermal/tests/network_closed_forms.rs
pub(crate) fn linearised_loss_conductance(
    environment: &Environment,
    at: Temperature,
    emissivity: f64,
) -> f64 {
    let area = environment.area.to_si();
    let t = at.to_si().max(0.0);
    environment.convection_w_per_m2_k * area
        + 4.0 * emissivity * STEFAN_BOLTZMANN.to_si() * area * t * t * t
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The box from the report that opened this: 1.122 kg, 96x60x60 mm, still air, 21 W.
    ///
    /// A real part from a downstream project rather than an invented one, and the emissivity is
    /// the parameter because that is exactly what the old functions were blind to.
    fn radiating_box(emissivity: f64) -> LumpedMass {
        use dualis_units::SpecificHeat;
        let (area, vol) = (0.030_24, 3.456e-4);
        let mut substance = Substance::aluminium_6061();
        substance.density = Density::kg_per_m3(1.122 / vol);
        if let Some(t) = substance.thermal.as_mut() {
            t.specific_heat = SpecificHeat::j_per_kg_k(600.0);
            t.emissivity = emissivity;
        }
        LumpedMass::new(
            "box",
            substance,
            Volume::from_si(vol),
            Length::from_si(vol / area),
            Temperature::celsius(25.0),
            Environment {
                ambient: Temperature::celsius(25.0),
                convection_w_per_m2_k: 7.0,
                area: Area::from_si(area),
            },
        )
    }

    /// Step to steady state under 21 W and report the settled rise.
    fn settle(body: &mut LumpedMass) -> f64 {
        let mut bus = Exchange::new();
        for k in 0..400_000 {
            bus.publish(HEAT, 21.0);
            body.step(Time::s(k as f64), Time::s(1.0), &mut bus)
                .unwrap();
        }
        body.rise().to_si()
    }

    use dualis_core::{Flux, Schedule, Simulation};
    use dualis_units::Density;

    fn lens_volume() -> Volume {
        // A 25 mm disc, 5 mm thick.
        Volume::from_si(std::f64::consts::PI * 0.0125f64.powi(2) * 0.005)
    }

    fn lens_area() -> Area {
        // Two faces plus the rim.
        let r = 0.0125f64;
        Area::from_si(2.0 * std::f64::consts::PI * r * r + std::f64::consts::TAU * r * 0.005)
    }

    fn lens(initial_c: f64) -> LumpedMass {
        LumpedMass::new(
            "lens",
            Substance::borosilicate_crown(),
            lens_volume(),
            Length::mm(5.0),
            Temperature::celsius(initial_c),
            Environment::still_air(Temperature::celsius(20.0), lens_area()),
        )
    }

    /// The lumped approximation is only honest when the Biot number is small, and the
    /// domain says which situation it is in rather than assuming.
    #[test]
    fn the_lumped_approximation_declares_when_it_applies() {
        let glass_in_air = lens(20.0);
        assert!(
            glass_in_air.biot_number() < 0.1,
            "still air over 5 mm of glass should be lumpable, Bi = {}",
            glass_in_air.biot_number()
        );

        // The same glass in flowing water is not lumpable: h is two orders larger.
        let mut wet = lens(20.0);
        wet.environment.convection_w_per_m2_k = 500.0;
        assert!(wet.biot_number() > 1.0, "Bi = {}", wet.biot_number());

        // A substance with no thermal data cannot answer at all, and says so.
        let unknown = LumpedMass::new(
            "unknown",
            Substance::bulk("x", Density::g_per_cm3(2.0)),
            lens_volume(),
            Length::mm(5.0),
            Temperature::celsius(20.0),
            Environment::still_air(Temperature::celsius(20.0), lens_area()),
        );
        assert!(!unknown.biot_number().is_finite());
    }

    /// Newton's law of cooling has a closed form: `ΔT(t) = ΔT₀ exp(-t/τ)`. Integrated
    /// with steps under the reported limit, the domain reproduces it — which is what
    /// says both the physics and the stability limit are right.
    #[test]
    fn cooling_follows_the_exponential_it_should() {
        // Radiation would add a second loss path and spoil the closed form, so this
        // one case turns it off to test the convective term alone.
        let mut body = lens(30.0);
        body.substance.thermal.as_mut().unwrap().emissivity = 0.0;
        let tau = body.time_constant();
        // 5.3 J/K of heat capacity against 9.6 mW/K of still-air conductance is
        // 549 s — nine minutes. A glass lens is a poor conductor with real heat
        // capacity, so it settles slowly, and that slowness is why its thermal domain
        // takes one step per video frame while a metal mount's does not.
        assert!(
            (tau.to_si() - 549.0).abs() < 5.0,
            "a lens in still air settles over about nine minutes, tau = {} s",
            tau.to_si()
        );

        let initial_rise = body.rise().to_si();

        // Run for one whole time constant, in `steps` equal pieces.
        let run = |steps: u32| {
            let mut body = lens(30.0);
            body.substance.thermal.as_mut().unwrap().emissivity = 0.0;
            let dt = tau / steps as f64;
            let mut bus = Exchange::new();
            for _ in 0..steps {
                body.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            body.rise().to_si()
        };

        // Against the *discrete* closed form first: explicit Euler on
        // `dT/dt = -(T-Ta)/tau` is exactly `(1 - h/tau)^n`, and the domain reproduces
        // it to machine precision. This is what says the update is the integrator it
        // claims to be, with no stray factor hiding in the heat capacity.
        for steps in [10u32, 100, 1000] {
            let discrete = initial_rise * (1.0 - 1.0 / steps as f64).powi(steps as i32);
            let got = run(steps);
            assert!(
                (got / discrete - 1.0).abs() < 1e-9,
                "{steps} steps: got {got:.6} K, discrete solution {discrete:.6} K"
            );
        }

        // Then against the *continuous* one, which is the physics. Euler undershoots
        // it, and the shortfall is first order in the step: 5.2% at tau/10, 0.51% at
        // tau/100, 0.05% at tau/1000 — a clean factor of ten each time, which is what
        // first order means and why `max_stable_dt` reports tau/10 rather than the
        // 2*tau that mere stability would allow.
        let exact = initial_rise * (-1.0f64).exp();
        let shortfall = |steps: u32| (exact - run(steps)) / exact;
        assert!((shortfall(10) - 0.0522).abs() < 1e-3, "{}", shortfall(10));
        assert!((shortfall(100) - 0.0051).abs() < 1e-3, "{}", shortfall(100));
        assert!(
            (shortfall(10) / shortfall(100) - 10.0).abs() < 1.0,
            "first order: ratio {}",
            shortfall(10) / shortfall(100)
        );
        assert!(run(10) < initial_rise, "it should have cooled");
    }

    /// Radiation is not a small correction at room temperature: for a black surface
    /// it is the same order as still-air convection, and leaving it out halves the
    /// loss. That is the mistake this test exists to make impossible.
    #[test]
    fn radiation_matters_as_much_as_convection() {
        let env = Environment::still_air(Temperature::celsius(20.0), lens_area());
        let hot = Temperature::celsius(30.0);
        let with_radiation = env.loss_from(hot, 0.90).to_si();
        let without = env.loss_from(hot, 0.0).to_si();
        let radiative = with_radiation - without;
        assert!(
            radiative / without > 0.6 && radiative / without < 1.0,
            "radiation is {:.2} of convection, not negligible",
            radiative / without
        );
        // At ambient nothing is lost either way.
        assert!(env.loss_from(Temperature::celsius(20.0), 0.9).to_si().abs() < 1e-12);
        // And a colder body gains, with the sign to prove it.
        assert!(env.loss_from(Temperature::celsius(10.0), 0.9).to_si() < 0.0);
    }

    /// **`equilibrium_rise` agrees with stepping the domain**, whatever the surface is.
    ///
    /// It used to be `P/(hA)`, which reports the same rise for a polished surface and a
    /// blackbody one, and was over by 2.09x at the top of that range — always in the
    /// comfortable direction. It solves the full balance now.
    ///
    /// Asserted against `step`, and that is the point: the crate had two public functions
    /// disagreeing with its own physics while `Environment::loss_from`'s documentation
    /// explained exactly why that would be wrong. Reported by a downstream project that built
    /// a conclusion on the old number and had to retract the mechanism.
    ///
    /// `quoted / settled` across the emissivity range was 1.07, 1.12, 1.37, 1.58, 1.99, 2.09.
    #[test]
    fn the_equilibrium_agrees_with_stepping_there_at_every_emissivity() {
        for e in [0.0, 0.05, 0.3, 0.9, 1.0] {
            let mut body = radiating_box(e);
            let quoted = body.equilibrium_rise(Power::w(21.0)).to_si();
            let settled = settle(&mut body);
            assert!(
                (quoted / settled - 1.0).abs() < 1e-6,
                "emissivity {e}: quoted {quoted:.4} K, settled {settled:.4} K"
            );
        }
        // And the surface matters, which is the whole complaint: a black box settles at about
        // half the rise of a polished one under the same load, where `P/(hA)` said 99.2 K for
        // both.
        let (mut black, mut shiny) = (radiating_box(1.0), radiating_box(0.05));
        let (b, s) = (settle(&mut black), settle(&mut shiny));
        assert!(b < 0.55 * s, "black {b:.1} K against polished {s:.1} K");
    }

    /// **The time constant tightens as the body warms**, and brackets the large-signal value.
    ///
    /// `4εσA·T³` grows with temperature, so this is a property of the state and not of the
    /// body. The measured figure — time to 63% of the settled rise — must lie between the value
    /// at rest and the value once hot, because a large-signal time constant is an average over
    /// the trajectory that the small-signal one bounds at each end.
    ///
    /// At ε = 0.9: 29.9 min at rest, 23.8 once settled, 25.6 large-signal. Before radiation was
    /// counted it was 53.0 whatever the surface was.
    #[test]
    fn the_time_constant_brackets_the_measured_one_and_tightens_when_hot() {
        let cold = radiating_box(0.9);
        let tau_cold = cold.time_constant().to_si();

        let mut body = radiating_box(0.9);
        let settled = settle(&mut body);
        let tau_hot = body.time_constant().to_si();
        assert!(
            tau_hot < tau_cold,
            "hot {tau_hot:.1} s against cold {tau_cold:.1} s"
        );

        let mut probe = radiating_box(0.9);
        let mut bus = Exchange::new();
        let mut t63 = f64::NAN;
        for k in 0..400_000 {
            bus.publish(HEAT, 21.0);
            probe
                .step(Time::s(k as f64), Time::s(1.0), &mut bus)
                .unwrap();
            let reached = probe.rise().to_si() >= settled * (1.0 - 1.0 / std::f64::consts::E);
            if t63.is_nan() && reached {
                t63 = k as f64;
            }
        }
        assert!(
            tau_hot < t63 && t63 < tau_cold,
            "the measured {:.1} min should lie between {:.1} and {:.1}",
            t63 / 60.0,
            tau_hot / 60.0,
            tau_cold / 60.0
        );
        // The scheduler follows it: a tenth of a time constant that is now shorter.
        assert!(body.max_stable_dt(Time::ZERO) < cold.max_stable_dt(Time::ZERO));
    }

    /// A body with no way to lose heat never settles, rather than settling instantly.
    #[test]
    fn a_body_that_cannot_lose_heat_has_no_equilibrium() {
        let sealed = LumpedMass::new(
            "sealed",
            Substance::aluminium_6061(),
            Volume::from_si(1e-4),
            Length::mm(10.0),
            Temperature::celsius(20.0),
            Environment {
                ambient: Temperature::celsius(20.0),
                convection_w_per_m2_k: 0.0,
                area: Area::from_si(0.0),
            },
        );
        assert!(sealed.equilibrium_rise(Power::w(1.0)).to_si().is_infinite());
        assert_eq!(sealed.equilibrium_rise(Power::w(0.0)).to_si(), 0.0);
        assert!(sealed.time_constant().to_si().is_infinite());
    }

    /// With no convection at all the root is the Stefan-Boltzmann one, exactly.
    ///
    /// The case the old `P/(hA)` divided by zero on and returned infinity for — a body in
    /// vacuum has an equilibrium, and it is a closed form.
    #[test]
    fn a_body_in_vacuum_settles_where_stefan_boltzmann_says() {
        let (area, vol) = (0.030_24, 3.456e-4);
        let emissivity = 0.8;
        let mut substance = Substance::aluminium_6061();
        if let Some(t) = substance.thermal.as_mut() {
            t.emissivity = emissivity;
        }
        let vacuum = LumpedMass::new(
            "vac",
            substance,
            Volume::from_si(vol),
            Length::from_si(vol / area),
            Temperature::celsius(25.0),
            Environment {
                ambient: Temperature::celsius(25.0),
                convection_w_per_m2_k: 0.0,
                area: Area::from_si(area),
            },
        );
        let rise = vacuum.equilibrium_rise(Power::w(21.0)).to_si();
        // Closed form, computed here: T = (P/(eps A sigma) + Ta^4)^(1/4).
        let ta = Temperature::celsius(25.0).to_si();
        let want =
            (21.0 / (emissivity * STEFAN_BOLTZMANN.to_si() * area) + ta.powi(4)).powf(0.25) - ta;
        assert!(
            (rise / want - 1.0).abs() < 1e-9,
            "vacuum: got {rise:.4} K, closed form {want:.4} K"
        );
    }

    /// The answer a designer wants: where does it end up.
    ///
    /// Asserted against stepping there rather than against a band, and the band is why. This
    /// test used to require 1 K to 20 K, which the convective-only formula satisfied and the
    /// real answer does not: 10 mW into this lens settles **0.60 K** up, because borosilicate
    /// radiates and the old formula pretended it did not. The band was wide enough to look
    /// generous and narrow enough to encode the bug.
    #[test]
    fn equilibrium_rise_is_the_number_that_matters() {
        let mut body = lens(20.0);
        let quoted = body.equilibrium_rise(Power::mw(10.0)).to_si();

        let mut bus = Exchange::new();
        for k in 0..200_000 {
            bus.publish(HEAT, 0.010);
            body.step(Time::s(k as f64), Time::s(1.0), &mut bus)
                .unwrap();
        }
        let settled = body.rise().to_si();
        assert!(
            (quoted / settled - 1.0).abs() < 1e-6,
            "quoted {quoted:.6} K, settled {settled:.6} K"
        );
        // Small enough to matter for a focus and not for the glass, which is the point of
        // the number — and that check lives on Substance.
        assert!(quoted > 0.1 && quoted < 2.0, "got {quoted} K");
        assert_eq!(
            Substance::borosilicate_crown().survives(Temperature::from_si(quoted)),
            Some(true)
        );
    }

    /// Explicit conduction refuses to run past its stability limit rather than
    /// diverging. This is the failure mode the whole `max_stable_dt` mechanism
    /// exists to prevent, and here it is caught even when a caller ignores it.
    #[test]
    fn explicit_conduction_refuses_an_unstable_step() {
        let mut bar = Bar1D::new(
            "bar",
            Substance::aluminium_6061(),
            20,
            Length::mm(1.0),
            Area::from_si(1e-4),
            Temperature::celsius(20.0),
        );
        let limit = bar.max_stable_dt(Time::ZERO);
        // Aluminium on a millimetre grid: about 7 ms.
        assert!(
            (limit.in_ms() - 7.2).abs() < 0.2,
            "limit {} ms",
            limit.in_ms()
        );
        assert!((bar.fourier_number(limit) - 0.5).abs() < 1e-12);

        let mut bus = Exchange::new();
        assert!(bar.step(Time::ZERO, limit, &mut bus).is_ok());
        let err = bar
            .step(Time::ZERO, limit * 1.5, &mut bus)
            .expect_err("past the limit must not be attempted");
        assert_eq!(err.quantity, "Fourier number");
        assert!(err.after > 0.5, "{err}");
    }

    /// Glass and aluminium differ by two orders of magnitude in how big a step they
    /// can take on the same grid — which is the concrete reason a multirate schedule
    /// is not a premature optimisation.
    #[test]
    fn two_materials_on_one_grid_need_different_steps() {
        let bar = |s: Substance| {
            Bar1D::new(
                "bar",
                s,
                10,
                Length::mm(1.0),
                Area::from_si(1e-4),
                Temperature::celsius(20.0),
            )
            .max_stable_dt(Time::ZERO)
            .to_si()
        };
        let glass = bar(Substance::borosilicate_crown());
        let metal = bar(Substance::aluminium_6061());
        assert!((glass - 0.967).abs() < 0.02, "glass {glass} s");
        assert!((metal - 0.0072).abs() < 0.001, "metal {metal} s");
        assert!(glass / metal > 100.0, "ratio {}", glass / metal);
    }

    /// Insulated conduction conserves heat exactly and flattens a gradient
    /// monotonically — the two things the heat equation is supposed to do.
    #[test]
    fn conduction_conserves_heat_and_flattens_a_gradient() {
        let mut bar = Bar1D::new(
            "bar",
            Substance::aluminium_6061(),
            21,
            Length::mm(1.0),
            Area::from_si(1e-4),
            Temperature::celsius(20.0),
        );
        // A hot spot in the middle.
        bar.cells[10] = Temperature::celsius(60.0).to_si();
        let total_before: f64 = bar.cells.iter().sum();
        let spread_before = bar.cells.iter().cloned().fold(0.0f64, f64::max)
            - bar.cells.iter().cloned().fold(f64::MAX, f64::min);

        let dt = bar.max_stable_dt(Time::ZERO);
        let mut bus = Exchange::new();
        for _ in 0..500 {
            bar.step(Time::ZERO, dt * 0.9, &mut bus).unwrap();
        }

        let total_after: f64 = bar.cells.iter().sum();
        assert!(
            (total_after / total_before - 1.0).abs() < 1e-12,
            "insulated ends must conserve heat exactly: {total_before} -> {total_after}"
        );
        let spread_after = bar.cells.iter().cloned().fold(0.0f64, f64::max)
            - bar.cells.iter().cloned().fold(f64::MAX, f64::min);
        assert!(
            spread_after < spread_before / 10.0,
            "the gradient should have flattened: {spread_before} -> {spread_after}"
        );
        // A lumped model would have reported zero gradient from the start; this one
        // still sees a little.
        assert!(bar.mean_temperature().in_celsius() > 21.0);
    }

    /// **What a place-aware flux buys, stated as the difference it makes.** The same
    /// joules, delivered to the same bar, once as a lumped total and once resolved over the
    /// boundary, land somewhere different — and the resolved one lands where the light
    /// actually was.
    ///
    /// This is the check that could not be written before: the lumped bus had no way to say
    /// "the middle", so both runs would have been the same run.
    #[test]
    fn heat_arrives_where_the_flux_says_it_did() {
        let build = || {
            Bar1D::new(
                "bar",
                Substance::aluminium_6061(),
                21,
                Length::mm(1.0),
                Area::from_si(1e-4),
                Temperature::celsius(20.0),
            )
        };
        let joules = 2.0;

        // Lumped: everything into cell 0, because there is nothing else it could mean.
        let mut lumped = build();
        let mut bus = Exchange::new();
        bus.publish(HEAT, joules);
        lumped
            .step(Time::ZERO, Time::from_si(1e-4), &mut bus)
            .unwrap();

        // Resolved: all of it onto face 10, the middle of the bar.
        let mut resolved = build().exposing("bar face", Area::from_si(1e-4));
        let boundary = resolved.boundary().expect("it was just given one").clone();
        assert_eq!(
            boundary.faces(),
            21,
            "one face per cell, so nothing interpolates"
        );
        let mut spot = vec![0.0; 21];
        spot[10] = joules;
        bus.publish_on(&boundary, HEAT, &Flux::from_faces(spot))
            .unwrap();
        resolved
            .step(Time::ZERO, Time::from_si(1e-4), &mut bus)
            .unwrap();

        // Same energy in, by the domain's own accounting.
        assert!(
            (resolved.absorbed_energy().to_si() - lumped.absorbed_energy().to_si()).abs() < 1e-15,
            "the two runs must differ in place, not in amount"
        );
        assert!(
            bus.unclaimed().next().is_none(),
            "and nothing was left on the bus"
        );

        // But a different bar. The lumped run is hot at the end it was told nothing about;
        // the resolved one is hot in the middle, where the beam was.
        let lumped_end = lumped.temperature_at(0).in_celsius();
        let lumped_middle = lumped.temperature_at(10).in_celsius();
        assert!(
            lumped_end > lumped_middle + 1.0,
            "lumped heat piled up at cell 0"
        );
        assert!(
            (lumped_middle - 20.0).abs() < 1e-9,
            "and never reached the middle"
        );

        let resolved_end = resolved.temperature_at(0).in_celsius();
        let resolved_middle = resolved.temperature_at(10).in_celsius();
        assert!(
            resolved_middle > resolved_end + 1.0,
            "resolved heat is in the middle"
        );
        assert!(
            (resolved_end - 20.0).abs() < 1e-9,
            "and the end is untouched"
        );
        // Symmetric about the spot, which a one-sided injection can never be.
        assert!(
            (resolved.temperature_at(9).in_celsius() - resolved.temperature_at(11).in_celsius())
                .abs()
                < 1e-12
        );
    }

    /// A distribution the bar cannot read is refused, and refusing it means the step fails
    /// rather than heating the wrong cells. The bar's grid is the discretisation; a
    /// publisher on a different one has to say so.
    #[test]
    fn a_flux_on_the_wrong_grid_stops_the_step() {
        let mut bar = Bar1D::new(
            "bar",
            Substance::aluminium_6061(),
            21,
            Length::mm(1.0),
            Area::from_si(1e-4),
            Temperature::celsius(20.0),
        )
        .exposing("bar face", Area::from_si(1e-4));
        let boundary = bar.boundary().unwrap().clone();

        // An illuminator on a 64-pixel grid, which is a perfectly reasonable thing to be.
        let camera = Interface::uniform("bar face", 64, Area::from_si(1e-4) * (21.0 / 64.0));
        let mut bus = Exchange::new();
        bus.publish_on(&camera, HEAT, &Flux::spread_over(2.0, &camera))
            .unwrap();

        let err = bar
            .step(Time::ZERO, Time::from_si(1e-4), &mut bus)
            .expect_err("a 64-face flux is not a 21-cell bar");
        assert!(err.quantity.contains("expected 21"), "{err}");
        assert!(
            (bar.temperature_at(10).in_celsius() - 20.0).abs() < 1e-9,
            "and nothing was heated"
        );

        // Resampling first is what works, and the bar then gets all of it.
        let crossed = bus
            .take_on(&camera, HEAT)
            .unwrap()
            .resample(&camera, &boundary)
            .unwrap();
        bus.publish_on(&boundary, HEAT, &crossed).unwrap();
        bar.step(Time::ZERO, Time::from_si(1e-4), &mut bus).unwrap();
        assert!((bar.absorbed_energy().to_si() - 2.0).abs() < 1e-12);
        assert!(bus.unclaimed().next().is_none());
    }

    /// The field reads the bar where the bar is, and holds its value past the insulated
    /// ends rather than running off the array.
    #[test]
    fn the_field_samples_the_bar_and_stops_at_its_ends() {
        let mut bar = Bar1D::new(
            "bar",
            Substance::aluminium_6061(),
            5,
            Length::mm(1.0),
            Area::from_si(1e-4),
            Temperature::celsius(20.0),
        );
        // A ramp, so every cell is distinguishable.
        for (i, cell) in bar.cells.iter_mut().enumerate() {
            *cell = 300.0 + i as f64;
        }
        let at = |mm: f64| bar.at(LengthVec::mm(mm, 0.0, 0.0), Time::ZERO);

        // Cell centres are at 0.5, 1.5, ... mm and read exactly.
        for i in 0..5 {
            assert!(
                (at(i as f64 + 0.5) - (300.0 + i as f64)).abs() < 1e-12,
                "cell {i}"
            );
        }
        // Halfway between two centres is halfway between their values.
        assert!((at(1.0) - 300.5).abs() < 1e-12);
        assert!((at(3.75) - 303.25).abs() < 1e-12);

        // Past either end the value is held. The ends are insulated, so the temperature
        // really does stop changing there — this is the physics, not a clamp for safety.
        assert!((at(-50.0) - 300.0).abs() < 1e-12);
        assert!((at(0.0) - 300.0).abs() < 1e-12);
        assert!((at(5.0) - 304.0).abs() < 1e-12);
        assert!((at(1e6) - 304.0).abs() < 1e-12);
        assert!(
            (at(f64::NAN) - 300.0).abs() < 1e-12,
            "a NaN must not index the array"
        );

        // Uniform across the bar, which is what a one-dimensional model asserts.
        assert_eq!(
            bar.at(LengthVec::mm(2.5, 0.0, 0.0), Time::ZERO),
            bar.at(LengthVec::mm(2.5, 40.0, -70.0), Time::ZERO)
        );
    }

    /// Gradient and Laplacian against closed forms. A linear ramp has a constant gradient
    /// and no curvature; a quadratic has a curvature the three-point stencil gets exactly,
    /// because a second difference of a quadratic is not an approximation.
    #[test]
    fn the_fields_derivatives_match_their_closed_forms() {
        let dx = 1e-3;
        let build = |f: &dyn Fn(f64) -> f64| {
            let mut bar = Bar1D::new(
                "bar",
                Substance::aluminium_6061(),
                21,
                Length::from_si(dx),
                Area::from_si(1e-4),
                Temperature::celsius(20.0),
            );
            for (i, cell) in bar.cells.iter_mut().enumerate() {
                *cell = f((i as f64 + 0.5) * dx);
            }
            bar
        };
        let probe = LengthVec::from_si(DVec3::new(10.5 * dx, 0.0, 0.0));
        let h = Length::from_si(dx);

        // T = 300 + 40x, so dT/dx = 40 K/m everywhere and the curvature is zero.
        let ramp = build(&|x| 300.0 + 40.0 * x);
        let g = ramp.gradient(probe, Time::ZERO, h);
        assert!((g.x - 40.0).abs() < 1e-9, "got {g}");
        assert!(
            g.y == 0.0 && g.z == 0.0,
            "a 1D bar has no transverse gradient"
        );
        assert!(
            ramp.laplacian(probe, Time::ZERO, h).abs() < 1e-6,
            "a ramp has no curvature"
        );

        // T = 300 + 5000x², so d²T/dx² = 10000 K/m² exactly, and dT/dx = 10000x.
        let curved = build(&|x| 300.0 + 5000.0 * x * x);
        let lap = curved.laplacian(probe, Time::ZERO, h);
        assert!((lap / 10_000.0 - 1.0).abs() < 1e-9, "got {lap}");
        let g = curved.gradient(probe, Time::ZERO, h);
        assert!((g.x / (10_000.0 * 10.5 * dx) - 1.0).abs() < 1e-9, "got {g}");

        // Past the ends the derivatives go to zero, because that is where `at` goes flat.
        // A field whose gradient disagreed with its own values would be worse than useless
        // to a visualiser, which reads both.
        for outside in [-5.0 * dx, 30.0 * dx] {
            let p = LengthVec::from_si(DVec3::new(outside, 0.0, 0.0));
            assert_eq!(
                curved.gradient(p, Time::ZERO, h),
                DVec3::ZERO,
                "at {outside}"
            );
            assert_eq!(curved.laplacian(p, Time::ZERO, h), 0.0, "at {outside}");
        }

        // But *at* the insulated wall the gradient is not zero, and pretending otherwise
        // would be a lie about what a cell-centred grid knows: its first sample is half a
        // cell inside, where the temperature really is still changing.
        let wall = LengthVec::ZERO;
        let expected = (curved.cells[1] - curved.cells[0]) / (2.0 * dx);
        assert!((curved.gradient(wall, Time::ZERO, h).x - expected).abs() < 1e-9);
        assert!(expected > 0.0, "the mirrored estimate is not zero");
    }

    /// **The check that makes the field worth having.** The rate it reports is exactly what
    /// the domain does on its next step — not approximately, to the last bit.
    ///
    /// That is not a coincidence and it is the reason `rate` is overridden. The explicit
    /// update *is* `T += α·dt·∇²T` on this stencil, so evaluating the governing equation and
    /// taking the step are the same arithmetic. A visualiser drawing `rate` is drawing what
    /// is about to happen, and the default finite difference in time — which would have
    /// needed history the domain does not keep — would have read zero.
    #[test]
    fn the_reported_rate_is_exactly_the_step_the_domain_takes() {
        let dx = 1e-3;
        let mut bar = Bar1D::new(
            "bar",
            Substance::aluminium_6061(),
            21,
            Length::from_si(dx),
            Area::from_si(1e-4),
            Temperature::celsius(20.0),
        );
        bar.cells[10] = Temperature::celsius(60.0).to_si();
        bar.cells[14] = Temperature::celsius(35.0).to_si();

        let dt = bar.max_stable_dt(Time::ZERO) * 0.5;
        let probes: Vec<LengthVec> = (0..21)
            .map(|i| LengthVec::from_si(DVec3::new((i as f64 + 0.5) * dx, 0.0, 0.0)))
            .collect();
        let predicted: Vec<f64> = probes
            .iter()
            .map(|p| bar.rate(*p, Time::ZERO, dt))
            .collect();
        let before: Vec<f64> = bar.cells.clone();

        bar.step(Time::ZERO, dt, &mut Exchange::new()).unwrap();

        for (i, p) in probes.iter().enumerate() {
            let observed = (bar.cells[i] - before[i]) / dt.to_si();
            let _ = p;
            if predicted[i].abs() < 1e-12 {
                assert!(
                    observed.abs() < 1e-9,
                    "cell {i}: {observed} against nothing"
                );
            } else {
                assert!(
                    (observed / predicted[i] - 1.0).abs() < 1e-12,
                    "cell {i}: predicted {} but the step did {observed}",
                    predicted[i]
                );
            }
        }
        // And the hot cell really was cooling while its neighbours warmed, so the test is
        // not passing on a bar where nothing happened.
        assert!(
            predicted[10] < -1.0,
            "the peak should be cooling: {}",
            predicted[10]
        );
        assert!(
            predicted[9] > 1.0,
            "and its neighbour warming: {}",
            predicted[9]
        );
    }

    /// A substance with no diffusivity has no conduction to report, rather than an
    /// infinity or a panic.
    #[test]
    fn a_field_with_no_diffusivity_reports_no_rate() {
        let bar = Bar1D::new(
            "bar",
            Substance::bulk("mystery", Density::from_si(1000.0)),
            5,
            Length::mm(1.0),
            Area::from_si(1e-4),
            Temperature::celsius(20.0),
        );
        assert_eq!(bar.rate(LengthVec::ZERO, Time::ZERO, Time::s(1.0)), 0.0);
    }

    /// The domain plugs into the kernel's scheduler and its books balance: heat taken
    /// from the bus is stored plus lost, to within the audit's tolerance, over a run
    /// long enough for both paths to matter.
    #[test]
    fn the_domain_balances_its_books_under_the_scheduler() {
        struct Heater {
            watts: f64,
            paid: f64,
        }
        impl Domain for Heater {
            fn name(&self) -> &str {
                "heater"
            }
            fn kind(&self) -> Kind {
                Kind::QuasiStatic
            }
            fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
                let joules = self.watts * dt.to_si();
                bus.publish(HEAT, joules);
                self.paid += joules;
                Ok(())
            }
            fn ledger(&self) -> Ledger {
                Ledger::new().with(quantity::ENERGY, -self.paid)
            }
            fn checkpoint(&mut self) {}
            fn restore(&mut self) {}
            fn supports_restore(&self) -> bool {
                true
            }
        }

        let mut sim = Simulation::new(Schedule::Multirate)
            .with(Heater {
                watts: 0.01,
                paid: 0.0,
            })
            .with(lens(20.0));

        for _ in 0..40 {
            sim.advance(Time::s(5.0)).expect("the books must balance");
        }
        // 40 windows of 5 s is 200 s, and the lens has warmed measurably.
        assert!((sim.time().to_si() - 200.0).abs() < 1e-9);
        // The audit already proved conservation; this is what it looks like.
        let total = sim.ledger().get(quantity::ENERGY).unwrap();
        assert!(total.abs() < 1e-9, "residual {total}");
    }

    /// A domain whose substance has no thermal data fails by name rather than
    /// silently doing nothing.
    #[test]
    fn a_substance_without_heat_capacity_is_refused() {
        let mut body = LumpedMass::new(
            "mystery",
            Substance::bulk("mystery", Density::g_per_cm3(3.0)),
            lens_volume(),
            Length::mm(5.0),
            Temperature::celsius(20.0),
            Environment::still_air(Temperature::celsius(20.0), lens_area()),
        );
        let mut bus = Exchange::new();
        let err = body.step(Time::ZERO, Time::s(1.0), &mut bus).unwrap_err();
        assert_eq!(err.site, "mystery");
        assert!(err.quantity.contains("heat capacity"), "{err}");
    }
}
