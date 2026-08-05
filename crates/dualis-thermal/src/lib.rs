//! dualis-thermal: heat, as a domain built on the `dualis-core` kernel.
//!
//! Two domains, chosen because they sit either side of the line that matters to a
//! scheduler:
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

use dualis_core::conserved::quantity;
use dualis_core::{Domain, Exchange, Kind, Ledger, Substance, Violation};
use dualis_units::{
    Area, Energy, HeatCapacity, Length, Power, Temperature, Time, Volume, STEFAN_BOLTZMANN,
};

/// The bus channel heat arrives on, in joules.
///
/// Joules rather than watts, because a domain steps over an interval and what
/// crossed the interface is an amount, not a rate. The publisher multiplies by its
/// own `dt`, which is what makes the audit an equality rather than an approximation.
pub const HEAT: &str = quantity::ENERGY;

/// How a body loses heat to its surroundings.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Environment {
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
    name: &'static str,
    substance: Substance,
    volume: Volume,
    /// Characteristic length for the Biot number — volume over surface area.
    thickness: Length,
    temperature: Temperature,
    environment: Environment,
    saved: Option<Temperature>,
    /// Joules taken from the bus over the run, for the books.
    absorbed: f64,
    /// Joules given up to the environment over the run.
    lost: f64,
}

impl LumpedMass {
    pub fn new(
        name: &'static str,
        substance: Substance,
        volume: Volume,
        thickness: Length,
        initial: Temperature,
        environment: Environment,
    ) -> LumpedMass {
        LumpedMass {
            name,
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

    pub fn temperature(&self) -> Temperature {
        self.temperature
    }

    /// Rise above ambient.
    pub fn rise(&self) -> Temperature {
        self.temperature - self.environment.ambient
    }

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
    pub fn time_constant(&self) -> Time {
        let capacity = self.heat_capacity().to_si();
        let conductance = self.environment.convection_w_per_m2_k * self.environment.area.to_si();
        if conductance <= 0.0 || !capacity.is_finite() {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(capacity / conductance)
    }

    /// Steady-state rise for a constant absorbed power, ignoring radiation.
    ///
    /// `ΔT = P/(hA)`. The answer a designer actually wants: not how it gets there,
    /// but where it ends up.
    pub fn equilibrium_rise(&self, absorbed: Power) -> Temperature {
        let conductance = self.environment.convection_w_per_m2_k * self.environment.area.to_si();
        if conductance <= 0.0 {
            return Temperature::from_si(f64::INFINITY);
        }
        Temperature::from_si(absorbed.to_si() / conductance)
    }

    fn emissivity(&self) -> f64 {
        self.substance.thermal.map(|t| t.emissivity).unwrap_or(0.0)
    }
}

impl Domain for LumpedMass {
    fn name(&self) -> &'static str {
        self.name
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
                self.name,
                "substance has no heat capacity",
                capacity.to_si(),
            ));
        }

        // Everything offered on the channel this step.
        let gained = bus.take(HEAT);
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

    fn checkpoint(&mut self) {
        self.saved = Some(self.temperature);
    }

    fn restore(&mut self) {
        if let Some(t) = self.saved {
            self.temperature = t;
        }
    }

    fn supports_restore(&self) -> bool {
        true
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
    name: &'static str,
    substance: Substance,
    /// Temperatures at the cell centres.
    cells: Vec<f64>,
    saved: Vec<f64>,
    dx: Length,
    /// Cross-sectional area, for turning joules into a temperature.
    area: Area,
    absorbed: f64,
    initial_heat: f64,
}

impl Bar1D {
    /// A bar of `cells` cells, each `dx` long, all starting at `initial`.
    pub fn new(
        name: &'static str,
        substance: Substance,
        cells: usize,
        dx: Length,
        area: Area,
        initial: Temperature,
    ) -> Bar1D {
        let cells = cells.max(2);
        let temps = vec![initial.to_si(); cells];
        let mut bar = Bar1D {
            name,
            substance,
            cells: temps.clone(),
            saved: temps,
            dx,
            area,
            absorbed: 0.0,
            initial_heat: 0.0,
        };
        bar.initial_heat = bar.stored_heat();
        bar
    }

    pub fn temperature_at(&self, index: usize) -> Temperature {
        Temperature::from_si(self.cells[index.min(self.cells.len() - 1)])
    }

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

    fn stored_heat(&self) -> f64 {
        // Against absolute zero, which is arbitrary but constant — only differences
        // enter the audit.
        self.cell_capacity() * self.cells.iter().sum::<f64>()
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
    fn name(&self) -> &'static str {
        self.name
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
            return Err(Violation::at(self.name, "substance has no diffusivity", f));
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

        // All the heat offered goes into the first cell, which is where a surface
        // absorbing light would put it.
        let gained = bus.take(HEAT);
        self.absorbed += gained;
        let capacity = self.cell_capacity();
        self.cells[0] += gained / capacity;

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
        Ledger::new().with(quantity::ENERGY, self.stored_heat() - self.initial_heat)
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use dualis_core::{Schedule, Simulation};
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

    /// The answer a designer wants: where does it end up. 10 mW into a lens in still
    /// air settles a few kelvin above ambient — enough to matter for a focus, not
    /// enough to break anything.
    #[test]
    fn equilibrium_rise_is_the_number_that_matters() {
        let body = lens(20.0);
        let rise = body.equilibrium_rise(Power::mw(10.0));
        assert!(
            rise.to_si() > 1.0 && rise.to_si() < 20.0,
            "10 mW should settle a few kelvin up, got {} K",
            rise.to_si()
        );
        // And the glass survives it comfortably — that check lives on Substance.
        assert_eq!(Substance::borosilicate_crown().survives(rise), Some(true));
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
            fn name(&self) -> &'static str {
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
