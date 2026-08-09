//! A Lennard-Jones fluid, as a domain the kernel can step.
//!
//! # What is being claimed, and how it is checked
//!
//! Every other domain in this workspace is checked against a value: a Fresnel coefficient, a
//! mode frequency, a temperature rise. This one mostly cannot be. A hundred atoms bouncing off
//! each other have no closed-form trajectory and are chaotic besides — two runs differing in
//! the last bit of one coordinate diverge completely inside a few hundred steps, and that is
//! the physics rather than a numerical failing.
//!
//! What *is* exact is what the trajectory averages to, and there is plenty of it:
//!
//! - **Equipartition.** `⟨KE⟩ = (3N − 3) k_BT / 2`, and the `−3` is not a rounding: momentum is
//!   conserved, so the centre of mass never moves and three degrees of freedom are frozen.
//!   Using `3N` overstates the temperature by `1/N`, which is a percent at a hundred atoms and
//!   is the classic way to be quietly wrong.
//! - **The ideal gas law.** At low density the virial term dies and `PV = Nk_BT` exactly.
//! - **The virial theorem.** At any density, `PV = Nk_BT + ⟨Σ f·r⟩/3`, which is the definition
//!   of the pressure rather than an approximation of it.
//! - **Energy conservation.** With no thermostat this is a symplectic integrator on a
//!   conservative system, so the total is bounded rather than drifting — the same statement
//!   [`velocity_verlet`](dualis_core::velocity_verlet) is checked for on a harmonic oscillator,
//!   holding here for a many-body potential.
//! - **Maxwell-Boltzmann.** Once thermostatted, the speeds settle into a distribution whose
//!   mean and variance are both fixed by one temperature.
//!
//! # Determinism, in a domain built out of random numbers
//!
//! The initial velocities and every Langevin kick come from
//! [`Rng::for_index`](dualis_core::Rng::for_index), keyed by step and particle. So the noise is
//! a pure function of `(seed, step, particle)` and a run is reproducible whatever order the
//! work is done in — which is the property the kernel's generator was built for and the one a
//! chaotic system makes impossible to recover any other way.

// Every public item carries a doc comment. See `dualis-units` for why this is denied.
#![deny(missing_docs)]

use dualis_core::conserved::quantity;
use dualis_core::{Domain, Exchange, Kind, Ledger, Rng, Violation};
use dualis_units::{Energy, Length, Mass, Pressure, Temperature, Time, Volume, BOLTZMANN};
use glam::DVec3;

use crate::box_::{CellList, PeriodicBox};
use crate::potential::LennardJones;

/// How the temperature is held, if it is.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Thermostat {
    /// None. The energy is conserved and the temperature is whatever it is — the NVE
    /// ensemble, and the only setting in which energy conservation means anything.
    None,
    /// Langevin: a drag proportional to velocity, and a random kick whose size is tied to the
    /// drag by the fluctuation-dissipation theorem.
    ///
    /// `damping` is `γ` in inverse seconds. Small is gentle and slow to equilibrate; large
    /// reaches the target quickly and distorts the dynamics on the way, since a thermostat
    /// strong enough to fix the temperature is also strong enough to change the diffusion
    /// coefficient. A tenth of the inverse collision time is the usual compromise.
    Langevin {
        /// The temperature the bath holds the fluid at, in the long run and on average.
        target: Temperature,
        /// `γ`, in inverse seconds. See the note on this variant for how to choose it.
        damping: f64,
    },
}

/// A box of atoms interacting through a pair potential.
pub struct Fluid {
    name: String,
    potential: LennardJones,
    bounds: PeriodicBox,
    mass: f64,
    positions: Vec<DVec3>,
    velocities: Vec<DVec3>,
    /// Accelerations from the last force evaluation, kept so a velocity-Verlet step costs one
    /// evaluation rather than two.
    accelerations: Vec<DVec3>,
    thermostat: Thermostat,
    seed: u64,
    /// Steps taken, which keys the noise so that it never repeats.
    step_count: u64,
    /// Potential energy at the current positions, from the same sweep that produced the
    /// accelerations — recomputing it separately would be a second `O(N)` pass and a second
    /// chance to disagree.
    potential_energy: f64,
    /// `Σ f·r` over pairs, likewise.
    virial: f64,
    /// Work the thermostat has done, in joules. Positive means it put energy in.
    thermostat_work: f64,
    saved: Option<(Vec<DVec3>, Vec<DVec3>, u64)>,
}

impl Fluid {
    /// Atoms on a face-centred cubic lattice, at rest.
    ///
    /// FCC rather than simple cubic because it is the packing argon actually freezes into, so
    /// a lattice at liquid density is close to a plausible configuration rather than a
    /// strained one. A simple cubic lattice at `ρ* = 0.8442` has neighbours well inside the
    /// repulsive wall and blows up on the first step.
    ///
    /// `cells` is the number of unit cells per side; each holds four atoms, so the count is
    /// `4·cells³` — 32, 108, 256, 500, 864 for two through six. Those are the numbers
    /// molecular-dynamics papers quote, and this is why.
    pub fn lattice(
        name: impl Into<String>,
        potential: LennardJones,
        mass: Mass,
        cells: usize,
        number_density: f64,
    ) -> Fluid {
        let cells = cells.max(1);
        let count = 4 * cells * cells * cells;
        let bounds = PeriodicBox::for_density(count, number_density);
        let a = bounds.length / cells as f64;
        // The four-atom basis of an FCC cell.
        let basis = [
            DVec3::ZERO,
            DVec3::new(0.5, 0.5, 0.0),
            DVec3::new(0.5, 0.0, 0.5),
            DVec3::new(0.0, 0.5, 0.5),
        ];
        let mut positions = Vec::with_capacity(count);
        for k in 0..cells {
            for j in 0..cells {
                for i in 0..cells {
                    let corner = DVec3::new(i as f64, j as f64, k as f64) * a;
                    for b in basis {
                        positions.push(corner + b * a);
                    }
                }
            }
        }
        let n = positions.len();
        let mut fluid = Fluid {
            name: name.into(),
            potential,
            bounds,
            mass: mass.to_si(),
            positions,
            velocities: vec![DVec3::ZERO; n],
            accelerations: vec![DVec3::ZERO; n],
            thermostat: Thermostat::None,
            seed: 0x_D05E_4D15,
            step_count: 0,
            potential_energy: 0.0,
            virial: 0.0,
            thermostat_work: 0.0,
            saved: None,
        };
        fluid.recompute_forces();
        fluid
    }

    /// Draw velocities from a Maxwell-Boltzmann distribution at `temperature`, then remove the
    /// centre-of-mass drift.
    ///
    /// Removing the drift matters twice. Physically, a box of gas sailing through its own
    /// periodic images is not the system anyone meant. Numerically, it is what makes the
    /// degrees of freedom `3N − 3` — and with the drift left in, the temperature computed from
    /// `3N − 3` would be systematically high by exactly the drift's share.
    ///
    /// Each component is an independent Gaussian of width `√(k_BT/m)`, which is what makes the
    /// *speed* Maxwell-Boltzmann. Drawn per particle from a keyed stream, so the configuration
    /// depends on the seed and not on the order.
    pub fn thermalised(mut self, temperature: Temperature, seed: u64) -> Fluid {
        self.seed = seed;
        let width = (BOLTZMANN.to_si() * temperature.to_si() / self.mass).sqrt();
        for (i, v) in self.velocities.iter_mut().enumerate() {
            let mut rng = Rng::for_index(seed, i as u64);
            *v = DVec3::new(
                rng.gaussian() * width,
                rng.gaussian() * width,
                rng.gaussian() * width,
            );
        }
        self.remove_drift();
        // The sampled temperature is not exactly the target — it is a draw from a
        // distribution — so rescale to land on it. Only the amplitude changes; the directions
        // are the ones the Gaussians chose.
        let sampled = self.temperature().to_si();
        if sampled > 0.0 {
            let factor = (temperature.to_si() / sampled).sqrt();
            for v in self.velocities.iter_mut() {
                *v *= factor;
            }
        }
        self
    }

    /// Attach a thermostat, or [`Thermostat::None`] for the NVE ensemble.
    ///
    /// Only without one does energy conservation mean anything, which is why the tests that
    /// check it never set one.
    pub fn with_thermostat(mut self, thermostat: Thermostat) -> Fluid {
        self.thermostat = thermostat;
        self
    }

    /// Subtract the mean velocity, so the centre of mass stands still.
    pub fn remove_drift(&mut self) {
        let n = self.velocities.len();
        if n == 0 {
            return;
        }
        let mean = self.velocities.iter().sum::<DVec3>() / n as f64;
        for v in self.velocities.iter_mut() {
            *v -= mean;
        }
    }

    /// How many atoms.
    pub fn count(&self) -> usize {
        self.positions.len()
    }

    /// The periodic box they live in.
    pub fn bounds(&self) -> PeriodicBox {
        self.bounds
    }

    /// Volume of the box. Fixed: there is no barostat, so this is the V of NVE and NVT.
    pub fn volume(&self) -> Volume {
        Volume::from_si(self.bounds.volume())
    }

    /// Particles per unit volume.
    pub fn number_density(&self) -> f64 {
        self.positions.len() as f64 / self.bounds.volume()
    }

    /// Position of one atom in metres, wrapped into the box. Clamped index.
    pub fn position(&self, i: usize) -> DVec3 {
        self.positions[i.min(self.positions.len() - 1)]
    }

    /// Velocity of one atom in m/s. Clamped index.
    pub fn velocity(&self, i: usize) -> DVec3 {
        self.velocities[i.min(self.velocities.len() - 1)]
    }

    /// Degrees of freedom: `3N − 3`.
    ///
    /// The three that are missing are the centre of mass, which the forces cannot move and
    /// [`Fluid::remove_drift`] set to zero. Counting them would make every reported temperature
    /// high by a factor of `3N/(3N−3)` — 1% at a hundred atoms, 0.1% at a thousand, and
    /// invisible either way unless you are checking against a closed form, which is exactly
    /// what this crate does.
    pub fn degrees_of_freedom(&self) -> f64 {
        (3 * self.positions.len()).saturating_sub(3) as f64
    }

    /// `½Σmv²`. With the drift removed this is entirely thermal, which is what makes
    /// [`Fluid::temperature`] a temperature.
    pub fn kinetic_energy(&self) -> Energy {
        let sum: f64 = self.velocities.iter().map(|v| v.length_squared()).sum();
        Energy::from_si(0.5 * self.mass * sum)
    }

    /// Summed pair energies, from the same sweep that produced the forces.
    ///
    /// Truncated **and shifted**, and the two need undoing in that order before this is
    /// comparable with a published number.
    ///
    /// [`LennardJones::energy_tail`] corrects the *truncated* potential, not the shifted one,
    /// so adding it to this sum directly is wrong — and wrong by more than the correction
    /// itself. The shift contributes `−n_pairs · u(rc)` to what is reported here, which at
    /// `rc = 2.5σ` and liquid density is about `+0.45 ε` per particle against a tail
    /// correction of `−0.45 ε`. Following the naive recipe leaves an answer that still drifts
    /// by ten times more with the cutoff than the uncorrected number should.
    ///
    /// Undo the shift first: add `n_pairs · LennardJones::shift()` back, then add the tail.
    /// Done that way the corrected energy per particle agrees across cutoffs from 2.5σ to
    /// 3.3σ to within 0.4%, where the raw values spread by 9%.
    pub fn potential_energy(&self) -> Energy {
        Energy::from_si(self.potential_energy)
    }

    /// Kinetic plus potential. Bounded rather than constant under Verlet, because a
    /// symplectic integrator conserves a shadow Hamiltonian and not the energy.
    pub fn total_energy(&self) -> Energy {
        self.kinetic_energy() + self.potential_energy()
    }

    /// Total momentum, which the pair forces cannot change.
    pub fn momentum(&self) -> DVec3 {
        self.velocities.iter().sum::<DVec3>() * self.mass
    }

    /// Temperature from equipartition: `T = 2·KE / (f·k_B)`.
    pub fn temperature(&self) -> Temperature {
        let f = self.degrees_of_freedom();
        if f <= 0.0 {
            return Temperature::from_si(0.0);
        }
        Temperature::from_si(2.0 * self.kinetic_energy().to_si() / (f * BOLTZMANN.to_si()))
    }

    /// Pressure by the virial theorem: `P = (N k_BT + Σf·r/3) / V`.
    ///
    /// The kinetic part alone is the ideal gas; the virial is everything the interactions add.
    /// It is negative where attraction dominates, which is why a real gas is easier to compress
    /// than an ideal one until the repulsive wall takes over.
    ///
    /// `N k_BT` and not `f k_BT/3`: the pressure counts all `N` particles, because the momentum
    /// constraint removes a degree of freedom from the energy and not from the volume
    /// derivative. Mixing the two conventions is a `1/N` error in the opposite direction from
    /// the temperature's.
    pub fn pressure(&self) -> Pressure {
        let n = self.positions.len() as f64;
        let ideal = n * BOLTZMANN.to_si() * self.temperature().to_si();
        Pressure::from_si((ideal + self.virial / 3.0) / self.bounds.volume())
    }

    /// `Σ f·r` over pairs, the interaction part of the pressure before it is divided by
    /// three and by the volume.
    pub fn virial(&self) -> f64 {
        self.virial
    }

    /// Energy the thermostat has added since the start; negative if it has taken some out.
    pub fn thermostat_work(&self) -> Energy {
        Energy::from_si(self.thermostat_work)
    }

    /// A step no larger than this keeps the fastest atom from crossing an appreciable fraction
    /// of the potential's width in one go.
    ///
    /// There is no CFL condition here — nothing propagates at a fixed speed — so the limit is
    /// a heuristic and is named as one: a hundredth of `σ√(m/ε)`, the time an atom with a
    /// typical thermal speed takes to cross the well. The usual reduced step of `0.005` is half
    /// of it, and the honest statement is that a molecular-dynamics step is chosen by watching
    /// the energy drift rather than by evaluating a formula.
    pub fn suggested_dt(&self) -> Time {
        let tau = self.potential.sigma * (self.mass / self.potential.epsilon).sqrt();
        Time::from_si(0.01 * tau)
    }

    /// Recompute accelerations, potential energy and virial in one sweep.
    fn recompute_forces(&mut self) {
        let list = CellList::build(self.bounds, self.potential.cutoff, &self.positions);
        for a in self.accelerations.iter_mut() {
            *a = DVec3::ZERO;
        }
        let mut energy = 0.0;
        let mut virial = 0.0;
        let (potential, mass) = (self.potential, self.mass);
        // Hoisted: `at_squared` would otherwise redo `shift()` and four products on every pair.
        let prepared = potential.prepared();
        let accelerations = &mut self.accelerations;
        list.for_each_pair(
            self.bounds,
            potential.cutoff,
            &self.positions,
            |i, j, d, r2| {
                if let Some(pair) = prepared.at_squared(r2) {
                    let force = d * pair.force_over_r;
                    // Equal and opposite, applied from one evaluation. This is what makes the
                    // momentum exact rather than nearly: the same bits are added to one
                    // particle and subtracted from the other.
                    //
                    // Divided **once** rather than once per particle. Two `force / mass` are six
                    // divisions on a `DVec3` where three will do, and computing the quotient a
                    // second time cannot give a different answer — so this is the same bits for
                    // half the divisions, which is the only kind of speed-up allowed to touch a
                    // pinned result.
                    let a = force / mass;
                    accelerations[i] += a;
                    accelerations[j] -= a;
                    energy += pair.energy;
                    virial += pair.force_over_r * r2;
                }
            },
        );
        self.potential_energy = energy;
        self.virial = virial;
    }
}

impl Domain for Fluid {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    fn max_stable_dt(&self, _now: Time) -> Time {
        self.suggested_dt()
    }

    /// One velocity-Verlet step, with the thermostat applied as a separate half-kick.
    ///
    /// Verlet by hand rather than through
    /// [`velocity_verlet`](dualis_core::velocity_verlet), for one reason worth stating: the
    /// kernel's version calls `acceleration` twice, and here the second evaluation of a step is
    /// the first evaluation of the next one. Keeping the accelerations halves the cost, which
    /// on the only expensive part of the calculation is not a micro-optimisation.
    fn step(&mut self, _t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        if !self.bounds.admits(self.potential.cutoff) {
            return Err(Violation {
                quantity: "minimum image".to_string(),
                site: format!("{} (periodic box)", self.name),
                before: self.bounds.length / 2.0,
                after: self.potential.cutoff,
                scale: self.bounds.length / 2.0,
                tolerance: 0.0,
            });
        }
        let h = dt.to_si();
        if h <= 0.0 {
            return Ok(());
        }

        // Half-kick, drift, recompute, half-kick.
        for (v, a) in self.velocities.iter_mut().zip(self.accelerations.iter()) {
            *v += *a * (h / 2.0);
        }
        for (p, v) in self.positions.iter_mut().zip(self.velocities.iter()) {
            *p = self.bounds.wrap(*p + *v * h);
        }
        self.recompute_forces();
        for (v, a) in self.velocities.iter_mut().zip(self.accelerations.iter()) {
            *v += *a * (h / 2.0);
        }

        self.apply_thermostat(h);
        self.step_count += 1;
        Ok(())
    }

    /// Energy held: kinetic plus potential, less whatever the thermostat put in.
    ///
    /// The subtraction is what makes a thermostatted run auditable at all. A Langevin bath is
    /// not part of the system, so without accounting for it the books would show energy
    /// appearing from nowhere — which is exactly what a thermostat does, and exactly what a
    /// conservation audit is supposed to refuse. Tracking the work turns it from a leak into a
    /// transfer.
    fn ledger(&self) -> Ledger {
        Ledger::new().with(
            quantity::ENERGY,
            self.total_energy().to_si() - self.thermostat_work,
        )
    }

    fn checkpoint(&mut self) {
        self.saved = Some((
            self.positions.clone(),
            self.velocities.clone(),
            self.step_count,
        ));
    }

    fn restore(&mut self) {
        if let Some((p, v, s)) = self.saved.clone() {
            self.positions = p;
            self.velocities = v;
            self.step_count = s;
            self.recompute_forces();
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }
}

impl Fluid {
    /// The Langevin half of a step: drag and a matched kick.
    ///
    /// ```text
    /// dv = −γ v dt + √(2 γ k_BT dt / m) · ξ
    /// ```
    ///
    /// The two terms are not independent. Drag alone cools to a standstill and noise alone
    /// heats without limit; the fluctuation-dissipation theorem fixes their ratio so that the
    /// steady state is a Boltzmann distribution at exactly `T`, and the `√γ` in the noise is
    /// that relation rather than a tuning constant. Getting it wrong does not fail loudly — it
    /// equilibrates to the wrong temperature, which is why the tests here check the temperature
    /// a long run *settles at* rather than the one it was asked for.
    ///
    /// Keyed by step and particle, so the noise is reproducible and independent between the
    /// two.
    fn apply_thermostat(&mut self, h: f64) {
        let Thermostat::Langevin { target, damping } = self.thermostat else {
            return;
        };
        if damping <= 0.0 {
            return;
        }
        let before = self.kinetic_energy().to_si();
        let kick = (2.0 * damping * BOLTZMANN.to_si() * target.to_si() * h / self.mass).sqrt();
        let decay = (-damping * h).exp();
        for (i, v) in self.velocities.iter_mut().enumerate() {
            // Exponential rather than `1 - γh`, so a large damping cannot overshoot into a
            // negative velocity — the same reason an implicit step is used for stiff decay
            // anywhere else.
            let mut rng = Rng::for_index(
                self.seed ^ 0x9E37_79B9_7F4A_7C15,
                self.step_count * self.positions.len() as u64 + i as u64,
            );
            let noise = DVec3::new(rng.gaussian(), rng.gaussian(), rng.gaussian());
            *v = *v * decay + noise * kick * ((1.0 - decay * decay) / (2.0 * damping * h)).sqrt();
        }
        // The bath does not push the box around either.
        self.remove_drift();
        self.thermostat_work += self.kinetic_energy().to_si() - before;
    }
}

/// Convenience: the reduced temperature `k_BT/ε` a fluid is at.
pub fn reduced_temperature(temperature: Temperature, potential: &LennardJones) -> f64 {
    BOLTZMANN.to_si() * temperature.to_si() / potential.epsilon
}

/// The reduced density `ρσ³`.
pub fn reduced_density(number_density: f64, potential: &LennardJones) -> f64 {
    number_density * potential.sigma.powi(3)
}

/// A temperature in reduced units, for setting up in the numbers papers quote.
pub fn temperature_from_reduced(reduced: f64, potential: &LennardJones) -> Temperature {
    Temperature::from_si(reduced * potential.epsilon / BOLTZMANN.to_si())
}

/// A mass, for a fluid in reduced units where `m = 1`.
pub fn unit_mass() -> Mass {
    Mass::from_si(1.0)
}

/// The distance `σ` as a length, for callers working in reduced units.
pub fn sigma_of(potential: &LennardJones) -> Length {
    Length::from_si(potential.sigma)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reduced units throughout: `ε = σ = m = 1`, so a temperature of 1.0 means `k_BT = ε` and
    /// the numbers are the ones papers quote.
    fn reduced(cells: usize, density: f64) -> Fluid {
        Fluid::lattice(
            "fluid",
            LennardJones::reduced(),
            unit_mass(),
            cells,
            density,
        )
    }

    fn kelvin(reduced: f64) -> Temperature {
        temperature_from_reduced(reduced, &LennardJones::reduced())
    }

    /// The lattice puts the right number of atoms in the right sized box, none of them on top
    /// of another.
    #[test]
    fn the_lattice_is_a_plausible_starting_configuration() {
        let fluid = reduced(3, 0.8442);
        assert_eq!(fluid.count(), 108, "four atoms to an fcc cell, 27 cells");
        assert!((fluid.number_density() - 0.8442).abs() < 1e-12);

        // Nearest neighbours on an fcc lattice sit at a/sqrt(2), where the cell edge follows
        // from four atoms a cell: a = (4/rho)^(1/3). Checked against that rather than against
        // a remembered number.
        let cell_edge = (4.0f64 / 0.8442).cbrt();
        let expected = cell_edge / 2f64.sqrt();
        let mut closest = f64::MAX;
        for i in 0..fluid.count() {
            for j in (i + 1)..fluid.count() {
                let d = fluid
                    .bounds()
                    .shortest(fluid.position(i), fluid.position(j))
                    .length();
                closest = closest.min(d);
            }
        }
        assert!(
            (closest / expected - 1.0).abs() < 1e-12,
            "nearest neighbour at {closest} sigma, closed form says {expected}"
        );
        // 1.188 sigma, which is just *outside* the well's minimum at 1.1225 — so every
        // neighbour starts in the attractive region and the lattice is under a slight tension
        // rather than being crushed. That is what lets the first step survive; a simple cubic
        // lattice at this density would put neighbours at 1.06 sigma, inside the wall.
        assert!(
            closest > LennardJones::reduced().minimum(),
            "the lattice starts strained"
        );
        assert!(
            closest < LennardJones::reduced().cutoff,
            "and neighbours are in range at all"
        );
    }

    /// **Equipartition, and the three degrees of freedom that are not there.**
    ///
    /// Momentum is conserved and the drift was removed, so the centre of mass is frozen and
    /// there are `3N − 3` degrees rather than `3N`. Counting all `3N` would report a
    /// temperature low by `1/N` — a percent here, and exactly the kind of error that survives
    /// forever because it looks like statistics.
    #[test]
    fn temperature_counts_the_degrees_of_freedom_that_move() {
        let fluid = reduced(3, 0.8442).thermalised(kelvin(1.2), 0xE011_5A17);
        assert_eq!(fluid.degrees_of_freedom(), 3.0 * 108.0 - 3.0);

        // Thermalising lands on the target exactly, because it rescales to it.
        let t = reduced_temperature(fluid.temperature(), &LennardJones::reduced());
        assert!((t - 1.2).abs() < 1e-12, "reduced temperature {t}");

        // Which is the same as saying the kinetic energy is (3N-3)kT/2.
        let expected = fluid.degrees_of_freedom() * BOLTZMANN.to_si() * kelvin(1.2).to_si() / 2.0;
        assert!((fluid.kinetic_energy().to_si() / expected - 1.0).abs() < 1e-12);

        // The drift really is gone, which is what licenses the -3.
        assert!(fluid.momentum().length() < 1e-15 * fluid.count() as f64);
    }

    /// **Momentum is exact, and that is what the half-neighbour loop bought.**
    ///
    /// Each pair is visited once and its force added to one particle and subtracted from the
    /// other, so the bits cancel rather than nearly cancelling. The same property `NBody` has
    /// and `TreeNBody` gives up, and here it is what licenses counting `3N − 3` degrees of
    /// freedom forever rather than only at the start.
    #[test]
    fn momentum_survives_a_long_run_to_the_last_bit() {
        let mut fluid = reduced(3, 0.7).thermalised(kelvin(1.5), 0x_1122_3344);
        let dt = fluid.max_stable_dt(Time::ZERO);
        let mut bus = Exchange::new();
        // The scale a drift is measured against: the total momentum a system with no
        // cancellation at all would carry, which is N times one particle's typical value. In
        // reduced units the mass is one, so that is N sqrt(2 KE / N).
        let scale = fluid.count() as f64
            * (2.0 * fluid.kinetic_energy().to_si() / fluid.count() as f64).sqrt();
        for _ in 0..500 {
            fluid.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        assert!(
            fluid.momentum().length() < 1e-12 * scale,
            "drift of {} after 500 steps, against a scale of {scale}",
            fluid.momentum().length()
        );
    }

    /// **Energy is bounded, not merely slow to drift.** Velocity Verlet is symplectic, so the
    /// total oscillates inside an envelope instead of walking away — the claim the kernel
    /// checks on a harmonic oscillator, here on a many-body potential with a truncated force.
    ///
    /// Measured as the spread over the run rather than as the endpoint difference, which is
    /// the distinction that matters: a dissipative integrator can land back near its start by
    /// luck, and a symplectic one cannot leave the band at all.
    #[test]
    fn energy_stays_inside_a_band_without_a_thermostat() {
        let mut fluid = reduced(3, 0.7).thermalised(kelvin(1.5), 0x_E0E0_E0E0);
        let dt = fluid.max_stable_dt(Time::ZERO);
        let mut bus = Exchange::new();
        let start = fluid.total_energy().to_si();
        let (mut lowest, mut highest) = (start, start);
        let mut late_lowest = f64::MAX;
        let mut late_highest = f64::MIN;
        for k in 0..4000 {
            fluid.step(Time::ZERO, dt, &mut bus).unwrap();
            let e = fluid.total_energy().to_si();
            lowest = lowest.min(e);
            highest = highest.max(e);
            if k >= 3000 {
                late_lowest = late_lowest.min(e);
                late_highest = late_highest.max(e);
            }
        }
        // Judged against the kinetic energy, which is the scale the fluctuation lives on. The
        // total is a difference of two larger numbers and passes near zero, so a relative
        // tolerance on it would mean nothing — the same trap the kernel's audit has a `scale`
        // field for.
        let scale = fluid.kinetic_energy().to_si();
        let band = (highest - lowest) / scale;
        assert!(band < 0.05, "the band was {band} of the kinetic energy");

        // Bounded rather than drifting: the last quarter sits inside the band the whole run
        // set, and its centre has not marched to one side.
        assert!(late_lowest >= lowest - 1e-15 && late_highest <= highest + 1e-15);
        let drift = (late_lowest + late_highest) / 2.0 - (lowest + highest) / 2.0;
        assert!(
            drift.abs() < 0.35 * (highest - lowest),
            "the band centre moved {drift} against a width of {}",
            highest - lowest
        );
    }

    /// **The ideal gas law**, which a real fluid obeys only when it is dilute enough to forget
    /// that it is real.
    ///
    /// The fluid has to be *equilibrated* first, and the reason is worth stating: a dilute
    /// lattice has no pairs inside the cutoff at all, so its virial is exactly zero and it
    /// reports `PV/Nk_BT = 1` to the last bit. That is not the ideal gas law being obeyed, it
    /// is a configuration with no physics in it, and a test passing on that would be measuring
    /// nothing. Melting the lattice puts neighbours at a spread of distances, which is what a
    /// gas is.
    ///
    /// The departure from one is then the second virial coefficient rather than an error:
    /// `PV/Nk_BT = 1 + B₂(T)ρ`. Negative at `T* = 2`, because attraction still dominates —
    /// which is why a real gas is easier to compress than an ideal one.
    ///
    /// # Why this averages over seeds
    ///
    /// Because one run does not settle it. The departure at `ρ* = 0.02` is 1.5%, and across
    /// four seeds it came out 1.68, 1.12, 1.81 and 1.41 percent — a spread of a third of the
    /// signal. The *ratio* between two densities is worse, landing anywhere from 1.35 to 2.92
    /// on a single seed while averaging to 2.06.
    ///
    /// A hundred atoms is a small sample and two thousand correlated snapshots are fewer
    /// independent ones than they look. The first version of this test used one seed, asserted
    /// the ratio was two, and passed — on the seed I happened to write down, and on two of the
    /// four above. That is the failure mode statistical tests have, and averaging is the fix
    /// rather than a wider tolerance.
    #[test]
    fn a_dilute_fluid_obeys_the_ideal_gas_law_and_departs_linearly() {
        let compressibility = |density: f64, seed: u64| {
            let target = kelvin(2.0);
            let mut fluid = reduced(3, density)
                .thermalised(target, seed)
                .with_thermostat(Thermostat::Langevin {
                    target,
                    damping: 1.0,
                });
            let dt = fluid.max_stable_dt(Time::ZERO);
            let mut bus = Exchange::new();
            let (mut sum, mut samples) = (0.0, 0.0);
            for k in 0..2500 {
                fluid.step(Time::ZERO, dt, &mut bus).unwrap();
                if k >= 800 {
                    let ideal =
                        fluid.count() as f64 * BOLTZMANN.to_si() * fluid.temperature().to_si();
                    sum += fluid.pressure().to_si() * fluid.volume().to_si() / ideal;
                    samples += 1.0;
                }
            }
            sum / samples
        };
        let seeds = [0x_1DEA_11A5u64, 0xAAAA, 0xBBBB, 0xCCCC];
        let mean = |density: f64| {
            seeds
                .iter()
                .map(|s| compressibility(density, *s))
                .sum::<f64>()
                / seeds.len() as f64
        };

        let thin = mean(0.02);
        let thick = mean(0.04);
        println!(
            "PV/NkT averaged over {} seeds: {thin:.5} at rho*=0.02, {thick:.5} at 0.04",
            seeds.len()
        );

        assert!(
            (thin - 1.0).abs() < 0.04,
            "a dilute gas should be nearly ideal, got PV/NkT = {thin}"
        );
        assert!(
            (thin - 1.0).abs() > 1e-6,
            "and not trivially so -- an exactly zero virial means nothing was tested"
        );
        assert!(thin < 1.0 && thick < 1.0, "B2 is negative at T* = 2");
        let ratio = (thick - 1.0) / (thin - 1.0);
        assert!(
            (ratio - 2.0).abs() < 0.35,
            "twice the density should be twice the departure, ratio {ratio}"
        );
    }

    /// **The virial**, checked against a pressure computed without the neighbour search.
    ///
    /// `Σ f·r` is accumulated during the force sweep from the cell list. Here it is recomputed
    /// by brute force over every pair, and the two have to agree. That crosses the neighbour
    /// search off the list of things the pressure could be silently wrong because of, which is
    /// the part of a molecular-dynamics code most likely to be subtly wrong.
    #[test]
    fn the_virial_agrees_with_a_brute_force_sum() {
        let fluid = reduced(3, 0.8442).thermalised(kelvin(1.0), 0x_1717_A115);
        let lj = LennardJones::reduced();
        let bounds = fluid.bounds();
        let mut brute = 0.0;
        for i in 0..fluid.count() {
            for j in (i + 1)..fluid.count() {
                let d = bounds.shortest(fluid.position(i), fluid.position(j));
                if let Some(pair) = lj.at_squared(d.length_squared()) {
                    brute += pair.force_over_r * d.length_squared();
                }
            }
        }
        assert!(brute.abs() > 1.0, "a dense fluid should have a real virial");
        assert!(
            (fluid.virial() / brute - 1.0).abs() < 1e-12,
            "cells gave {} and every pair gave {brute}",
            fluid.virial()
        );
    }

    /// **What a Langevin bath is for**: whatever the fluid starts at, it settles at the
    /// temperature the bath was given.
    ///
    /// Started deliberately cold, at a fifth of the target, so arriving there is the
    /// thermostat's doing rather than the initial condition's. Averaged over the second half
    /// of the run, since the first half is the approach.
    ///
    /// This is also the only test of the fluctuation-dissipation relation there can be. The
    /// drag and the noise are tied by `√(2γk_BT)`, and getting that wrong does not fail
    /// loudly — it equilibrates to the wrong temperature. So the check is on where it settles,
    /// never on where it was asked to go.
    #[test]
    fn a_langevin_bath_pulls_the_fluid_to_its_target() {
        let target = kelvin(1.4);
        let one_seed = |seed: u64| {
            let mut fluid = reduced(3, 0.6)
                .thermalised(kelvin(0.28), seed)
                .with_thermostat(Thermostat::Langevin {
                    target,
                    damping: 2.0,
                });
            let dt = fluid.max_stable_dt(Time::ZERO);
            let mut bus = Exchange::new();

            let mut sum = 0.0;
            let mut samples = 0.0;
            for k in 0..6000 {
                fluid.step(Time::ZERO, dt, &mut bus).unwrap();
                if k >= 3000 {
                    sum += fluid.temperature().to_si();
                    samples += 1.0;
                }
            }
            reduced_temperature(
                Temperature::from_si(sum / samples),
                &LennardJones::reduced(),
            )
        };

        // Four seeds, because one is not a measurement of a noisy quantity.
        //
        // The old comment reasoned that √(2/3N) = 7.9% a sample over three thousand samples
        // leaves 0.14% of error. That treats correlated samples as independent, and they are
        // not: the velocity correlation time is 1/gamma, so the run holds tens of independent
        // samples and not thousands. Measured rather than argued — across twelve seeds the
        // three-thousand-sample average ran from 1.376 to 1.420, a standard deviation of 0.96%
        // of the mean, which implies about fifty independent samples.
        //
        // So a single seed departs by up to 1.7%, and the 5% that stood here was five standard
        // deviations of slack. Averaging four seeds halves the spread to a standard error of
        // 0.48%. These four give 1.4207, 1.4170, 1.3891 and 1.4022 — a mean of 1.4072, 0.52%
        // high, which is the one sigma it should be. Two percent is three sigma of what is
        // left. More samples, not a wider tolerance: the tolerance got tighter, not looser.
        let seeds = [0x_BA77_0015, 0x_BA77_0016, 0x_BA77_0017, 0x_BA77_0018];
        let settled = seeds.iter().map(|s| one_seed(*s)).sum::<f64>() / seeds.len() as f64;
        assert!(
            (settled / 1.4 - 1.0).abs() < 0.02,
            "settled at a reduced temperature of {settled}, wanted 1.4"
        );
    }

    /// The bath pays for what it adds, so the books balance for a system that is explicitly
    /// not closed.
    ///
    /// It has to add energy here, because the fluid started cold. Without tracking that, the
    /// ledger would show joules appearing from nowhere — which is exactly what a conservation
    /// audit exists to refuse, and would refuse correctly.
    #[test]
    fn the_bath_pays_for_what_it_adds() {
        let mut fluid = reduced(3, 0.6)
            .thermalised(kelvin(0.3), 0x_BA77_0016)
            .with_thermostat(Thermostat::Langevin {
                target: kelvin(1.5),
                damping: 2.0,
            });
        let dt = fluid.max_stable_dt(Time::ZERO);
        let mut bus = Exchange::new();
        let opening = fluid.ledger().get(quantity::ENERGY).unwrap();
        let energy_before = fluid.total_energy().to_si();
        for _ in 0..2000 {
            fluid.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        assert!(
            fluid.thermostat_work().to_si() > 0.0,
            "warming a cold fluid means putting energy in"
        );
        // The two agree to the *integrator's* band and not exactly, and that is the honest
        // statement rather than a slack tolerance. Verlet conserves a shadow Hamiltonian
        // rather than the energy, so between thermostat kicks the total wanders inside a
        // envelope a few percent of the kinetic energy wide -- measured in
        // `energy_stays_inside_a_band_without_a_thermostat`. Asking for 1e-9 here would be
        // asking the dynamics to be something it is not.
        let gained = fluid.total_energy().to_si() - energy_before;
        let band = fluid.kinetic_energy().to_si();
        assert!(
            (gained - fluid.thermostat_work().to_si()).abs() < 0.05 * band,
            "gained {gained} but the bath logged {}, against a kinetic scale of {band}",
            fluid.thermostat_work().to_si()
        );
        // So the ledger has barely moved: what the fluid holds, less what the bath gave it.
        // Against the kinetic energy, because the quantity itself is a difference that passes
        // near zero and has no scale of its own.
        let closing = fluid.ledger().get(quantity::ENERGY).unwrap();
        assert!(
            (closing - opening).abs() < 0.05 * band,
            "the books moved from {opening} to {closing}, against a scale of {band}"
        );
    }

    /// **Determinism, in the one domain where nothing else could recover it.**
    ///
    /// The trajectory is chaotic: two runs differing in the last bit of one coordinate
    /// separate completely inside a few hundred steps, so "close enough" is not available as a
    /// notion. Either the arithmetic is identical or the answers are unrelated. Keying the
    /// noise on `(seed, step, particle)` makes it identical.
    #[test]
    fn a_seed_fixes_the_whole_trajectory() {
        let run = |seed: u64| {
            // Three cells rather than two: 32 atoms at this density give a box of 3.58
            // sigma, whose half is inside the 2.5 cutoff and which the domain refuses.
            let mut fluid = reduced(3, 0.7)
                .thermalised(kelvin(1.0), seed)
                .with_thermostat(Thermostat::Langevin {
                    target: kelvin(1.0),
                    damping: 1.0,
                });
            let dt = fluid.max_stable_dt(Time::ZERO);
            let mut bus = Exchange::new();
            for _ in 0..200 {
                fluid.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            (fluid.position(7), fluid.total_energy().to_si())
        };
        let (p1, e1) = run(0x_5EED_0001);
        let (p2, e2) = run(0x_5EED_0001);
        assert_eq!(p1, p2, "the same seed must give the same positions");
        assert_eq!(
            e1.to_bits(),
            e2.to_bits(),
            "and the same energy, bit for bit"
        );

        let (p3, _) = run(0x_5EED_0002);
        assert!(
            (p1 - p3).length() > 1e-6,
            "a different seed should give a different run"
        );
    }

    /// A box too small for the minimum image is refused rather than silently letting a
    /// particle interact with a second copy of its neighbour.
    #[test]
    fn a_box_smaller_than_twice_the_cutoff_is_refused() {
        // Four atoms at this density gives a box of 2.71 sigma against a cutoff of 2.5.
        let mut cramped = Fluid::lattice("cramped", LennardJones::reduced(), unit_mass(), 1, 0.2);
        assert!(!cramped.bounds().admits(2.5), "the test needs a small box");
        let err = cramped
            .step(Time::ZERO, Time::from_si(1e-3), &mut Exchange::new())
            .expect_err("a cutoff past L/2 must not be stepped");
        assert_eq!(err.quantity, "minimum image");

        // And one that is big enough runs.
        let mut roomy = reduced(4, 0.2);
        assert!(roomy.bounds().admits(2.5));
        let dt = roomy.max_stable_dt(Time::ZERO);
        assert!(roomy.step(Time::ZERO, dt, &mut Exchange::new()).is_ok());
    }
}
