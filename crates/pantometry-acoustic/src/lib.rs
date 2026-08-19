//! pantometry-acoustic: sound, as a domain on the `pantometry-core` kernel.
//!
//! The linear wave equation on a grid, in one dimension, two or three:
//!
//! - [`Tube`] — a duct, a rod, an organ pipe — with ends that can be open, closed, or matched to
//!   a different medium. Modes are a harmonic series, `n·c/2L`.
//! - [`Room`] — a floor plan, where the modes stop being a series and start being a lattice.
//! - [`Hall`] — with a ceiling. The vertical and oblique modes exist here and **do not exist at
//!   all** in a `Room`; a 2.4 m ceiling puts the first one at 71 Hz, well inside the range a
//!   room is judged on.
//!
//! The three are not a progression to be climbed. Each dimension costs the cells and `√d` in the
//! Courant limit on top, so a `Room` is the right model for a floor plan and a `Tube` for a duct.
//! `Hall` is for the question the other two cannot answer rather than the one they answer more
//! cheaply.
//!
//! ```text
//! ∂²p/∂t² = c² ∂²p/∂x²
//! ```
//!
//! # Why this domain and not a fluid one
//!
//! Sound is what a fluid does when nothing much is happening: small pressure
//! variations about a resting state, where the governing equations linearise and every
//! answer has a closed form to check against. Standing modes in a tube are `n c/2L`
//! exactly. A reflection off an impedance step is `(Z₂-Z₁)/(Z₂+Z₁)` exactly. A pulse
//! travels at `c` and arrives when it should.
//!
//! Full Navier-Stokes has none of that, and this workspace's README has said from the
//! start why: turbulence at any interesting Reynolds number is a question about
//! supercomputer budgets rather than about API design, and a solver that could not be
//! checked against anything would be decoration. So the fluid domain here is the linear
//! one, and it is complete rather than partial.
//!
//! # The fourth domain, and what it had to add to the kernel
//!
//! Nothing. That is the claim the crate split was for, and this is the fourth time it
//! has held: optics, heat, mechanics and now sound, none of them naming another and
//! none of them needing the kernel changed.
//!
//! What it does exercise for the first time is a **[`Kind::Evolving`] domain whose
//! stability limit is a genuine CFL condition** — `dt ≤ dx/c`, from a wave speed rather
//! than from a diffusion coefficient or a spring. Going past it does not degrade the
//! answer, it explodes within a handful of steps, so the limit is reported and the step
//! is refused rather than attempted.

// Every public item carries a doc comment. Denied rather than warned: a public physics API
// whose `Length::mm` shows a blank summary in rustdoc is documented in the sense that a
// paragraph exists somewhere, and not in the sense a reader needs.
#![deny(missing_docs)]
pub mod hall;
pub mod room;

pub use hall::Hall;
pub use room::Room;

use pantometry_core::conserved::quantity;
use pantometry_core::{Domain, Exchange, Kind, Ledger, Substance, Violation};
use pantometry_units::{Area, Density, Energy, Frequency, Length, Pressure, Time, Velocity};

/// Characteristic acoustic impedance, `ρc`, in Pa·s·m⁻¹.
///
/// What decides how much of a wave crosses a boundary and how much comes back. Two
/// media with the same impedance are acoustically the same medium however different
/// they are otherwise, which is why an ultrasound probe is coupled with gel and not with
/// air.
pub type Impedance = pantometry_units::Qty<-2, 1, -1, 0, 0, 0, 0>;

/// The impedance of a medium: `ρ c`.
pub fn impedance(density: Density, sound_speed: Velocity) -> Impedance {
    Impedance::from_si(density.to_si() * sound_speed.to_si())
}

/// Impedance of a substance, if its density and sound speed are both known.
pub fn impedance_of(substance: &Substance) -> Option<Impedance> {
    let acoustic = substance.acoustic?;
    Some(impedance(substance.density, acoustic.sound_speed))
}

/// Pressure reflection coefficient at a step from `from` to `into`:
/// `(Z₂ - Z₁)/(Z₂ + Z₁)`.
///
/// Positive into a harder medium and negative into a softer one, and the sign is the
/// physics: a wave hitting a closed end comes back the same way up, and one hitting an
/// open end comes back inverted. That inversion is why an open pipe's fundamental is an
/// octave below a closed one of the same length.
pub fn reflection_coefficient(from: Impedance, into: Impedance) -> f64 {
    let (z1, z2) = (from.to_si(), into.to_si());
    if z1 + z2 == 0.0 {
        return 0.0;
    }
    (z2 - z1) / (z2 + z1)
}

/// Fraction of *power* transmitted across an impedance step.
///
/// `1 - R²`, and it falls off brutally: air to water is a factor of 3600 in impedance
/// and lets through 0.1% of the power. That is why you cannot hear anything useful
/// through the surface of a pond, and why a stethoscope touches skin.
pub fn power_transmission(from: Impedance, into: Impedance) -> f64 {
    let r = reflection_coefficient(from, into);
    (1.0 - r * r).clamp(0.0, 1.0)
}

/// What the ends of a tube do to a wave.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum End {
    /// Rigid: pressure doubles, particle velocity is zero. A stopped pipe.
    Closed,
    /// Free: pressure is zero and the reflection is inverted. An open pipe's mouth.
    Open,
    /// Matched to an impedance, which absorbs a wave rather than returning it. `Z`
    /// equal to the tube's own is a perfect absorber — the numerical equivalent of an
    /// anechoic termination.
    Matched(Impedance),
}

/// A tube of gas, discretised along its length.
///
/// Pressure and particle velocity on a **staggered** grid: pressures at cell centres,
/// velocities at the faces between them, each leapfrogged over the other. That is the
/// standard acoustic finite-difference scheme, and it is worth the extra array for two
/// reasons that a pressure-only formulation cannot match.
///
/// The energy is honest. Acoustic energy is
/// `∫ [p²/(2ρc²) + ρu²/2] dV`, and both halves of it are stored rather than one being
/// reconstructed from a time difference. A pressure-only scheme conserves the energy of
/// the wave equation *in p*, which is the energy of the velocity potential and not the
/// energy of the sound — the two are different functionals and only one of them is in
/// joules.
///
/// And the boundaries are physical. An impedance is a statement relating pressure to
/// velocity, `u = p/Z`; with a velocity to apply it to, a matched end is that equation
/// and the power it absorbs is `p·u·A`, which is the definition of acoustic intensity
/// rather than an approximation of it.
pub struct Tube {
    name: String,
    /// Pressure at cell centres.
    pressure: Vec<f64>,
    /// Pressure one whole step earlier, so [`Tube::energy`] can report the quantity the
    /// scheme conserves exactly. See [`Room::energy`](crate::room::Room::energy) for why
    /// the obvious expression is not that quantity.
    pressure_prev: Vec<f64>,
    /// Particle velocity at the faces between cells, so one shorter.
    velocity: Vec<f64>,
    dx: f64,
    speed: f64,
    density: f64,
    area: f64,
    left: End,
    right: End,
    /// Whether `velocity` has been offset to the half step the scheme carries it at. See
    /// [`Room`](crate::room::Room)'s field of the same name; the defect and the fix are
    /// identical, because the two schemes are the same one in different dimensions.
    velocity_staggered: bool,
    /// The released state's energy minus the invariant, once staggered.
    energy_offset: f64,
    /// The physical energy at release, the datum the reported energy is measured against.
    energy_datum: f64,
    saved: Option<(Vec<f64>, Vec<f64>, bool)>,
    /// Energy handed to the bus by absorbing ends, in joules.
    radiated: f64,
}

impl Tube {
    /// A tube of the given length, at rest.
    pub fn new(
        name: impl Into<String>,
        length: Length,
        cells: usize,
        area: Area,
        density: Density,
        sound_speed: Velocity,
    ) -> Tube {
        let cells = cells.max(3);
        Tube {
            name: name.into(),
            pressure: vec![0.0; cells],
            pressure_prev: vec![0.0; cells],
            velocity_staggered: false,
            energy_offset: 0.0,
            energy_datum: 0.0,
            velocity: vec![0.0; cells - 1],
            dx: length.to_si() / (cells - 1) as f64,
            speed: sound_speed.to_si(),
            density: density.to_si(),
            area: area.to_si(),
            left: End::Closed,
            right: End::Closed,
            saved: None,
            radiated: 0.0,
        }
    }

    /// Air at 20 °C: 343 m/s and 1.204 kg/m³.
    pub fn of_air(name: impl Into<String>, length: Length, cells: usize, area: Area) -> Tube {
        Tube::new(
            name,
            length,
            cells,
            area,
            Density::kg_per_m3(1.204),
            Velocity::m_per_s(343.0),
        )
    }

    /// Set what each end does to a wave. Both are [`End::Closed`] until told otherwise.
    ///
    /// An absorbing end tightens [`Domain::max_stable_dt`], because it drains its half-width
    /// boundary cell on a timescale of its own.
    pub fn with_ends(mut self, left: End, right: End) -> Tube {
        self.left = left;
        self.right = right;
        self
    }

    /// Start with a pressure profile and no motion.
    ///
    /// At rest at `t = 0`, which is where an initial condition lives and *not* where a
    /// staggered scheme carries velocity. The first step makes up the missing half; see
    /// `velocity_staggered`. Storing the velocity explicitly rather than implying it from two
    /// equal time levels is what makes that correction possible at all.
    ///
    /// A pressure bump released from rest splits into two pulses going opposite ways, each
    /// half the height — worth knowing before reading a result off one of them.
    pub fn released_from(mut self, profile: impl Fn(Length) -> Pressure) -> Tube {
        for i in 0..self.pressure.len() {
            self.pressure[i] = profile(Length::from_si(i as f64 * self.dx)).to_si();
        }
        self.pressure_prev.copy_from_slice(&self.pressure);
        self.velocity.iter_mut().for_each(|u| *u = 0.0);
        self.velocity_staggered = false;
        self.energy_offset = 0.0;
        self.energy_datum = self.invariant_si();
        self
    }

    /// How many pressure samples along the tube.
    pub fn cells(&self) -> usize {
        self.pressure.len()
    }

    /// Length of the tube. Note `(n − 1)·dx`: the samples sit *on* the ends, so `n` of them
    /// span `n − 1` gaps.
    pub fn length(&self) -> Length {
        Length::from_si((self.pressure.len() - 1) as f64 * self.dx)
    }

    /// Speed of sound in the medium filling it.
    pub fn sound_speed(&self) -> Velocity {
        Velocity::m_per_s(self.speed)
    }

    /// The tube's own characteristic impedance.
    pub fn impedance(&self) -> Impedance {
        Impedance::from_si(self.density * self.speed)
    }

    /// Acoustic pressure at one sample, clamped to the ends. Signed: sound swings either
    /// side of the ambient pressure, which is not represented here at all.
    pub fn pressure_at(&self, cell: usize) -> Pressure {
        Pressure::from_si(self.pressure[cell.min(self.pressure.len() - 1)])
    }

    /// Pressure along the tube.
    pub fn pressure(&self) -> Vec<Pressure> {
        self.pressure
            .iter()
            .map(|p| Pressure::from_si(*p))
            .collect()
    }

    /// Particle velocity at the face between cell `i` and `i + 1`.
    pub fn velocity_at(&self, face: usize) -> Velocity {
        Velocity::from_si(self.velocity[face.min(self.velocity.len() - 1)])
    }

    /// Peak pressure anywhere.
    pub fn peak_pressure(&self) -> Pressure {
        Pressure::from_si(self.pressure.iter().fold(0.0f64, |a, p| a.max(p.abs())))
    }

    /// Acoustic energy: `∫ [p²/(2ρc²) + ρu²/2] dV`.
    ///
    /// Both halves are stored rather than one being inferred, so this is the energy of
    /// the sound and not of a related quantity. Conserved for rigid or pressure-release
    /// ends, and reduced by exactly what an absorbing end published to the bus.
    pub fn energy(&self) -> Energy {
        Energy::from_si(self.invariant_si() + self.energy_offset)
    }

    /// How much the discrete invariant differs from the released state's physical energy.
    ///
    /// Zero until the first step; `O(dt²)` and converging. See
    /// [`Room::startup_adjustment`](crate::room::Room::startup_adjustment).
    pub fn startup_adjustment(&self) -> Energy {
        Energy::from_si(self.energy_offset)
    }

    /// The scheme's invariant, before the release datum is applied.
    fn invariant_si(&self) -> f64 {
        let volume = self.area * self.dx;
        let rc2 = self.density * self.speed * self.speed;
        if rc2 <= 0.0 {
            return 0.0;
        }
        // End cells hold half a cell's worth, matching the weighting their update uses —
        // which is what keeps this conserved exactly rather than nearly. Every velocity face
        // is between two cells, so none of them is a special case.
        let last = self.pressure.len() - 1;
        let potential: f64 = self
            .pressure
            .iter()
            .zip(self.pressure_prev.iter())
            .enumerate()
            .map(|(i, (p, prev))| {
                let share = if i == 0 || i == last { 0.5 } else { 1.0 };
                share * p * prev / (2.0 * rc2) * volume
            })
            .sum();
        let kinetic: f64 = self
            .velocity
            .iter()
            .map(|u| self.density * u * u / 2.0 * volume)
            .sum();
        potential + kinetic
    }

    /// Joules handed to the bus by absorbing ends over the run.
    pub fn radiated_energy(&self) -> Energy {
        Energy::from_si(self.radiated)
    }

    /// The `n`th standing mode of this tube, exactly.
    ///
    /// A tube closed at both ends, or open at both, resonates at `n c/2L`. One closed
    /// and one open resonates at the odd multiples of `c/4L` — an octave lower for the
    /// same length, which is why a stopped organ pipe sounds an octave below an open one
    /// of the same height.
    pub fn mode_frequency(&self, n: u32) -> Frequency {
        let l = self.length().to_si();
        if l <= 0.0 || n == 0 {
            return Frequency::from_si(0.0);
        }
        let symmetric = matches!(
            (self.left, self.right),
            (End::Closed, End::Closed) | (End::Open, End::Open)
        );
        if symmetric {
            Frequency::from_si(n as f64 * self.speed / (2.0 * l))
        } else {
            // Odd harmonics of c/4L.
            Frequency::from_si((2 * n - 1) as f64 * self.speed / (4.0 * l))
        }
    }

    /// Courant number `c dt/dx`, which must not exceed 1.
    pub fn courant(&self, dt: Time) -> f64 {
        self.speed * dt.to_si() / self.dx
    }

    /// Velocity at a wall, given the end condition and the pressure just inside it.
    ///
    /// Rigid walls do not move. A pressure-release end has whatever velocity the wave
    /// arriving at it implies, which for a matched impedance is `p/Z` — and the general
    /// case interpolates between the two by the reflection coefficient, so `Z → ∞` is
    /// rigid and `Z → 0` is free.
    fn wall_velocity(&self, end: End, boundary_pressure: f64, outward: f64) -> f64 {
        match end {
            End::Closed => 0.0,
            // A free surface cannot hold pressure, so the wall moves as fast as the
            // incoming wave asks: u = p/Z with Z the tube's own.
            End::Open => outward * boundary_pressure / self.impedance().to_si(),
            End::Matched(z) => {
                let zz = z.to_si();
                if zz <= 0.0 {
                    return outward * boundary_pressure / self.impedance().to_si();
                }
                outward * boundary_pressure / zz
            }
        }
    }

    /// One staggered step: velocities from the pressure gradient, then pressures from
    /// the velocity divergence.
    fn advance(&mut self, dt: f64, bus: &mut Exchange) {
        let n = self.pressure.len();
        let rc2 = self.density * self.speed * self.speed;

        // Half a step the first time. The stored velocity starts where the initial condition
        // put it, at t = 0, and the scheme wants it at t = -dt/2; kicking it a whole step
        // instead is the leapfrog startup error, O(dt) and permanent.
        let vh = if self.velocity_staggered {
            dt
        } else {
            0.5 * dt
        };
        let starting = !self.velocity_staggered;
        self.velocity_staggered = true;

        // Faces first: rho du/dt = -dp/dx.
        for i in 0..self.velocity.len() {
            self.velocity[i] -=
                vh / (self.density * self.dx) * (self.pressure[i + 1] - self.pressure[i]);
        }

        // Then cells: dp/dt = -rho c^2 du/dx. The two walls supply the velocities just
        // outside the first and last faces. `outward` is +1 at the right-hand end and
        // -1 at the left, so that a positive pressure pushes the wall outwards in both
        // cases.
        let left_wall = self.wall_velocity(self.left, self.pressure[0], -1.0);
        let right_wall = self.wall_velocity(self.right, self.pressure[n - 1], 1.0);

        self.pressure_prev.copy_from_slice(&self.pressure);
        let mut absorbed = 0.0;
        for i in 0..n {
            let inflow = if i == 0 {
                left_wall
            } else {
                self.velocity[i - 1]
            };
            let outflow = if i == n - 1 {
                right_wall
            } else {
                self.velocity[i]
            };
            // The two end cells own half a cell each, so their divergence is divided by
            // half the spacing. See `Room`'s `wall_weight` for what leaving this out costs:
            // a second-order interior converging at first order.
            let share = if i == 0 || i == n - 1 { 2.0 } else { 1.0 };
            self.pressure[i] -= dt * rc2 / self.dx * share * (outflow - inflow);
        }

        // Power leaving through each wall is the acoustic intensity `p u A`, which is
        // the definition rather than an approximation of one. Rigid walls do not move,
        // so they radiate nothing and this is exactly zero for them.
        for (pressure, wall) in [
            (self.pressure[0], left_wall),
            (self.pressure[n - 1], right_wall),
        ] {
            // A wall moving outwards under positive pressure is doing work on whatever
            // is beyond it.
            absorbed += (pressure * wall).abs() * self.area * dt;
        }
        if absorbed > 0.0 {
            self.radiated += absorbed;
            bus.publish(quantity::ENERGY, absorbed);
        }

        // Startup, once. See `Room`'s step for the full account: converting the released
        // state's velocity from t = 0 to the half step the scheme carries it at moves the
        // invariant by O(dt²), which is discretisation and not a leak, so it is taken as a
        // datum rather than handed to the audit. Bounded, so a real first-step bug cannot
        // hide in it.
        if starting {
            let invariant = self.invariant_si();
            let offset = self.energy_datum - invariant;
            let scale = self.energy_datum.abs().max(invariant.abs());
            if scale > 0.0 && offset.abs() / scale <= 0.25 {
                self.energy_offset = offset;
            }
        }
    }
}

impl Domain for Tube {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// The CFL condition, `dx/c`, exactly.
    ///
    /// The first genuine CFL limit in the workspace: a wave must not cross more than one
    /// cell per step, because the three-point stencil only knows about its neighbours
    /// and information travelling further than that has nowhere to come from. Unlike the
    /// diffusion and contact limits elsewhere here, this one is not a factor of safety
    /// chosen for accuracy — at exactly `dx/c` the scheme is *non-dissipative*, and
    /// stepping shorter costs accuracy in the form of numerical dispersion rather than
    /// buying it.
    /// `dx/c` for a closed tube, and **half that** if either end absorbs.
    ///
    /// The extra factor is the impedance boundary's own limit rather than the wave's. An end
    /// with `u = p/Z` drains its cell exponentially, and because that cell owns only half a
    /// spacing the rate is `2ρc²/(Z·dx)` — so a step needs `dt ≤ Z·dx/(2ρc²)` for the decay
    /// to be monotone, which for a matched end is `dx/2c`.
    ///
    /// At the full `dx/c` the factor is exactly `−1`: the end inverts the wave instead of
    /// swallowing it. Not a divergence, so nothing here would have caught it — it looks like
    /// a perfectly stable run that reflects, and the only symptom is a duct that rings when
    /// it was meant to be anechoic.
    fn max_stable_dt(&self, _now: Time) -> Time {
        if self.speed <= 0.0 {
            return Time::from_si(f64::INFINITY);
        }
        let wave = self.dx / self.speed;
        let rc2 = self.density * self.speed * self.speed;
        let drain = |end: End| match end {
            End::Closed => f64::INFINITY,
            End::Open => self.dx * self.impedance().to_si() / (2.0 * rc2),
            End::Matched(z) => {
                let zz = if z.to_si() > 0.0 {
                    z.to_si()
                } else {
                    self.impedance().to_si()
                };
                self.dx * zz / (2.0 * rc2)
            }
        };
        Time::from_si(wave.min(drain(self.left)).min(drain(self.right)))
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let courant = self.courant(dt);
        if courant > 1.0 + 1e-12 {
            return Err(Violation {
                quantity: "Courant number".to_string(),
                site: format!("{} (explicit wave equation)", self.name),
                before: 1.0,
                after: courant,
                scale: 1.0,
                tolerance: 1e-12,
            });
        }
        if dt.to_si() <= 0.0 {
            return Ok(());
        }
        self.advance(dt.to_si(), bus);
        Ok(())
    }

    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.energy().to_si())
    }

    fn checkpoint(&mut self) {
        self.saved = Some((
            self.pressure.clone(),
            self.velocity.clone(),
            self.velocity_staggered,
        ));
    }

    fn restore(&mut self) {
        if let Some((pressure, velocity, staggered)) = self.saved.clone() {
            self.pressure_prev.copy_from_slice(&pressure);
            self.pressure = pressure;
            self.velocity = velocity;
            self.velocity_staggered = staggered;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn air_tube(cells: usize) -> Tube {
        Tube::of_air("tube", Length::m(1.0), cells, Area::from_si(1e-4))
    }

    /// **The end cells own half a cell, stated exactly.** One step from rest moves every
    /// pressure by `h²c²∇²p` on the mirrored three-point stencil, the closed ends included.
    ///
    /// Same check as `room::tests::one_step_from_rest_is_the_laplacian_the_field_reports`,
    /// and for the same reason: without the weighting an end cell moves by exactly half as
    /// much, which reads as a frequency a percent low rather than as anything obviously
    /// wrong. Here it is a mismatch at a named cell.
    #[test]
    fn one_step_from_rest_is_the_three_point_laplacian() {
        let n = 33;
        let mut tube = air_tube(n)
            .released_from(|x| Pressure::from_si((std::f64::consts::PI * x.to_si() / 1.0).cos()));
        let dx = tube.dx;
        let c = tube.speed;
        let h = tube.max_stable_dt(Time::ZERO).to_si() * 0.7;
        let before = tube.pressure.clone();
        tube.step(Time::ZERO, Time::from_si(h), &mut Exchange::new())
            .unwrap();

        let scale = before.iter().fold(0.0f64, |a, v| a.max(v.abs())) / (dx * dx);
        for i in 0..n {
            // Mirrored at the ends — the ghost outside a rigid wall is the *reflection* of
            // the cell inside it, `p₋₁ = p₁`, not a copy of the wall cell. Clamping instead
            // of reflecting halves the stencil there, which is exactly the bug being tested
            // for, so getting it wrong here would have hidden it.
            let left = if i == 0 { before[1] } else { before[i - 1] };
            let right = if i + 1 == n {
                before[n - 2]
            } else {
                before[i + 1]
            };
            let expected = (left - 2.0 * before[i] + right) / (dx * dx);
            // The ½ is Taylor's: from rest ṗ(0) = 0, so p(h) − p(0) = ½h²c²∇²p. See the
            // matching test in `room`, which asserted this without the half and passed —
            // against a scheme that kicked the velocity a whole step where it was owed half.
            let implied = (tube.pressure[i] - before[i]) / (0.5 * h * h * c * c);
            assert!(
                (implied - expected).abs() < scale * 1e-12,
                "cell {i}: the step implies {implied} but the stencil says {expected}"
            );
        }
    }

    /// An absorbing end restricts the step beyond what the wave does, and by a factor the
    /// impedance sets.
    ///
    /// A matched end drains its half-width cell at `2ρc²/(Z dx)`, so it needs `dt` under
    /// `Z dx / 2ρc²` — half the CFL limit when `Z` is the tube's own. At the full limit the
    /// factor is exactly `−1` and the end inverts the wave rather than absorbing it, which
    /// is stable, silent, and wrong.
    #[test]
    fn an_absorbing_end_tightens_the_step_and_a_closed_one_does_not() {
        let closed = air_tube(101);
        let wave = closed.max_stable_dt(Time::ZERO).to_si();
        assert!((wave - closed.dx / closed.speed).abs() < 1e-15);

        let matched = air_tube(101).with_ends(End::Closed, End::Matched(closed.impedance()));
        assert!(
            (matched.max_stable_dt(Time::ZERO).to_si() / (wave / 2.0) - 1.0).abs() < 1e-12,
            "a matched end should halve it, got {}",
            matched.max_stable_dt(Time::ZERO).to_si()
        );

        // Twice the impedance absorbs half as fast, so it costs half as much: the limit goes
        // back to the wave's own, and no further — the wave still has to be resolved.
        let stiff = air_tube(101).with_ends(
            End::Closed,
            End::Matched(Impedance::from_si(closed.impedance().to_si() * 4.0)),
        );
        assert!((stiff.max_stable_dt(Time::ZERO).to_si() - wave).abs() < 1e-15);

        // An open end is the same condition against the tube's own impedance.
        let open = air_tube(101).with_ends(End::Open, End::Closed);
        assert!((open.max_stable_dt(Time::ZERO).to_si() - wave / 2.0).abs() < 1e-15);
    }

    /// The impedances of the media everyone compares, against the published figures.
    #[test]
    fn impedances_match_the_published_figures() {
        let air = impedance(Density::kg_per_m3(1.204), Velocity::m_per_s(343.0));
        assert!(
            (air.to_si() - 413.0).abs() < 1.0,
            "air {} Pa s/m",
            air.to_si()
        );

        let water = impedance_of(&Substance::water()).expect("water has a sound speed");
        assert!(
            (water.to_si() / 1.48e6 - 1.0).abs() < 0.02,
            "water {} Pa s/m",
            water.to_si()
        );
        // Glass is stiffer still.
        let glass = impedance_of(&Substance::borosilicate_crown()).unwrap();
        assert!(glass > water && water > air);

        // A substance with no acoustic data says so rather than guessing.
        assert_eq!(
            impedance_of(&Substance::bulk("x", Density::kg_per_m3(1.0))),
            None
        );
    }

    /// The reason you cannot hear anything through the surface of a pond: air to water
    /// is a factor of 3600 in impedance and passes a thousandth of the power.
    #[test]
    fn an_impedance_mismatch_reflects_almost_everything() {
        let air = impedance(Density::kg_per_m3(1.204), Velocity::m_per_s(343.0));
        let water = impedance_of(&Substance::water()).unwrap();
        assert!(
            (water / air - 3585.0).abs() < 100.0,
            "ratio {}",
            water / air
        );

        let r = reflection_coefficient(air, water);
        assert!(r > 0.999, "almost all of it comes back, got {r}");
        let t = power_transmission(air, water);
        assert!(
            t > 0.0005 && t < 0.002,
            "about a thousandth gets through, got {t:e}"
        );

        // The sign carries the physics: into a harder medium the reflection keeps its
        // sign, into a softer one it inverts.
        assert!(reflection_coefficient(water, air) < -0.999);
        // And a matched pair reflects nothing at all.
        assert_eq!(reflection_coefficient(air, air), 0.0);
        assert!((power_transmission(air, air) - 1.0).abs() < 1e-15);
    }

    /// The CFL condition is a bound and not a suggestion. Past it the scheme does not
    /// degrade, it diverges — so the step is refused.
    #[test]
    fn the_courant_limit_is_enforced() {
        let mut tube = air_tube(101);
        // 1 m in 100 intervals is 10 mm a cell, and 343 m/s crosses that in 29.2 us.
        let limit = tube.max_stable_dt(Time::ZERO);
        assert!(
            (limit.in_us() - 29.15).abs() < 0.05,
            "limit {} us",
            limit.in_us()
        );
        assert!((tube.courant(limit) - 1.0).abs() < 1e-12);

        let mut bus = Exchange::new();
        assert!(tube.step(Time::ZERO, limit, &mut bus).is_ok());
        let err = tube
            .step(Time::ZERO, limit * 1.01, &mut bus)
            .expect_err("past the Courant limit must be refused");
        assert_eq!(err.quantity, "Courant number");
        assert!(err.after > 1.0);
    }

    /// A pulse travels at the speed of sound, and arrives when it should. The most
    /// direct statement the wave equation makes, and the easiest to get wrong by a
    /// factor.
    #[test]
    fn a_pulse_travels_at_the_speed_of_sound() {
        let cells = 401;
        let mut tube = Tube::of_air("tube", Length::m(4.0), cells, Area::from_si(1e-4))
            .with_ends(End::Closed, End::Closed)
            .released_from(|x| {
                // A narrow bump a quarter of the way along.
                let centre = 1.0;
                let width = 0.05;
                let u = (x.to_si() - centre) / width;
                Pressure::from_si(100.0 * (-u * u).exp())
            });
        let dt = tube.max_stable_dt(Time::ZERO);
        let dx = 4.0 / (cells - 1) as f64;

        // Released from rest, the bump splits into two halves going opposite ways. Follow
        // the right-going one for a known time and see where its crest is.
        let travel = Time::s(0.005); // 1.715 m at 343 m/s
        let steps = (travel.to_si() / dt.to_si()).round() as u32;
        let mut bus = Exchange::new();
        for _ in 0..steps {
            tube.step(Time::ZERO, dt, &mut bus).unwrap();
        }

        // Find the crest to the right of the start.
        let pressures = tube.pressure();
        let start_cell = (1.0 / dx).round() as usize;
        let (mut best, mut best_at) = (0.0f64, start_cell);
        for (i, p) in pressures.iter().enumerate().skip(start_cell + 5) {
            if p.to_si() > best {
                best = p.to_si();
                best_at = i;
            }
        }
        let distance = (best_at - start_cell) as f64 * dx;
        let expected = 343.0 * dt.to_si() * steps as f64;
        assert!(
            (distance / expected - 1.0).abs() < 0.02,
            "the crest moved {distance:.4} m where the speed of sound predicts \
             {expected:.4} m"
        );
        // And it is still recognisably a pulse rather than numerical mush.
        assert!(best > 30.0, "the crest should have survived, got {best}");
    }

    /// Standing modes at exactly `n c/2L` for a symmetric tube, and at the odd
    /// multiples of `c/4L` when one end is closed and the other open.
    ///
    /// The second is the more interesting claim: it is an octave below the first for the
    /// same length, which is why a stopped organ pipe is half the height of an open one
    /// at the same pitch.
    #[test]
    fn standing_modes_are_where_the_closed_form_puts_them() {
        let closed = air_tube(101).with_ends(End::Closed, End::Closed);
        // 343 / 2 = 171.5 Hz for a 1 m tube.
        assert!((closed.mode_frequency(1).to_si() - 171.5).abs() < 0.1);
        assert!((closed.mode_frequency(2).to_si() - 343.0).abs() < 0.1);
        assert!((closed.mode_frequency(3).to_si() - 514.5).abs() < 0.1);
        // Every harmonic is present.
        for n in 1..6u32 {
            let f = closed.mode_frequency(n).to_si();
            assert!((f / closed.mode_frequency(1).to_si() - n as f64).abs() < 1e-9);
        }

        // Open at both ends is the same series.
        let open = air_tube(101).with_ends(End::Open, End::Open);
        assert!((open.mode_frequency(1).to_si() - 171.5).abs() < 0.1);

        // One of each: the fundamental is an octave lower, and only odd harmonics.
        let stopped = air_tube(101).with_ends(End::Closed, End::Open);
        assert!(
            (stopped.mode_frequency(1).to_si() - 85.75).abs() < 0.1,
            "c/4L = 85.75 Hz, got {}",
            stopped.mode_frequency(1).to_si()
        );
        assert!(
            (closed.mode_frequency(1) / stopped.mode_frequency(1) - 2.0).abs() < 1e-9,
            "an octave, exactly"
        );
        // Odd multiples only: the second mode is three times the first, not twice.
        assert!((stopped.mode_frequency(2) / stopped.mode_frequency(1) - 3.0).abs() < 1e-9);
        assert!((stopped.mode_frequency(3) / stopped.mode_frequency(1) - 5.0).abs() < 1e-9);

        assert_eq!(closed.mode_frequency(0).to_si(), 0.0);
    }

    /// A closed tube's fundamental, measured rather than asserted: excite the mode
    /// shape, run for one predicted period, and it should have come back.
    ///
    /// This is what connects the closed form above to the code actually stepping.
    #[test]
    fn the_fundamental_oscillates_at_its_predicted_period() {
        let cells = 201;
        let length = 1.0;
        let mut tube = Tube::of_air("tube", Length::m(length), cells, Area::from_si(1e-4))
            .with_ends(End::Closed, End::Closed)
            .released_from(|x| {
                // The first mode of a tube with rigid ends: a half cosine in pressure.
                let k = std::f64::consts::PI / length;
                Pressure::from_si(100.0 * (k * x.to_si()).cos())
            });
        let period = 1.0 / tube.mode_frequency(1).to_si();
        assert!((period - 1.0 / 171.5).abs() < 1e-6);

        let dt = tube.max_stable_dt(Time::ZERO);
        let steps = (period / dt.to_si()).round() as u32;
        let start = tube.pressure_at(0).to_si();
        let mut bus = Exchange::new();

        // Half a period should have turned it upside down.
        for _ in 0..steps / 2 {
            tube.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        let half = tube.pressure_at(0).to_si();
        assert!(
            half < -0.9 * start,
            "half a period should invert the mode: {start:.2} to {half:.2}"
        );

        // And a whole period brings it back.
        for _ in 0..steps - steps / 2 {
            tube.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        let whole = tube.pressure_at(0).to_si();
        assert!(
            (whole / start - 1.0).abs() < 0.05,
            "a whole period should restore it: {start:.2} to {whole:.2}"
        );
    }

    /// Rigid ends neither absorb nor radiate, so the tube's energy stays put and there
    /// is nothing on the bus.
    #[test]
    fn a_closed_tube_keeps_its_energy_to_itself() {
        let mut tube = air_tube(201)
            .with_ends(End::Closed, End::Closed)
            .released_from(|x| {
                let u = (x.to_si() - 0.5) / 0.08;
                Pressure::from_si(200.0 * (-u * u).exp())
            });
        let start = tube.energy().to_si();
        assert!(start > 0.0, "a displaced tube holds energy");

        let dt = tube.max_stable_dt(Time::ZERO);
        let mut bus = Exchange::new();
        for _ in 0..2000 {
            tube.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        assert_eq!(tube.radiated_energy().to_si(), 0.0);
        assert!(bus.peek(quantity::ENERGY).abs() < 1e-30);
        // Exactly, not approximately: the leapfrog is lossless and `energy` reports the
        // time-centred quantity it conserves rather than one that wobbles.
        let now = tube.energy().to_si();
        assert!(
            (now / start - 1.0).abs() < 1e-9,
            "energy went from {start:e} to {now:e}"
        );
    }

    /// A matched end absorbs instead of reflecting, and the joules it takes go onto the
    /// bus — the same channel optics and mechanics publish heat on, so a thermal domain
    /// could pick up the sound a duct dumps into its walls.
    #[test]
    fn a_matched_end_absorbs_and_publishes_what_it_took() {
        let z = air_tube(3).impedance();
        let mut tube = air_tube(201)
            .with_ends(End::Closed, End::Matched(z))
            .released_from(|x| {
                let u = (x.to_si() - 0.3) / 0.05;
                Pressure::from_si(150.0 * (-u * u).exp())
            });
        let start = tube.energy().to_si();
        let dt = tube.max_stable_dt(Time::ZERO);
        let mut bus = Exchange::new();
        let mut collected = 0.0;
        for _ in 0..3000 {
            tube.step(Time::ZERO, dt, &mut bus).unwrap();
            collected += bus.take(quantity::ENERGY);
        }
        assert!(
            collected > 0.0,
            "a matched end should have absorbed something"
        );
        // The bus and the books agree trivially -- `step` does `self.radiated += absorbed`
        // and `bus.publish(..., absorbed)` from the same variable, so this is `a == a` and
        // catches nothing. Kept because it is cheap and would notice one of the two lines
        // being deleted, but stated as the tautology it is rather than dressed up with a
        // 1e-30 that looks like a precision claim.
        assert_eq!(collected, tube.radiated_energy().to_si());

        // **The check that has content**: what the end published against what the tube
        // actually lost. Nothing connected those two before, and multiplying the radiated
        // intensity by seven passed the entire workspace suite.
        let left = tube.energy().to_si();
        let lost = start - left;
        assert!(
            (collected / lost - 1.0).abs() < 2e-3,
            "published {collected:e} J against {lost:e} J that left the tube"
        );
        // The 2e-3 is the half-step stagger, not slack: `p·u·A` multiplies a pressure at a
        // whole step by a wall velocity half a step later, so the intensity carries an
        // O(omega dt) error that the time-centred energy functional does not.

        // And the tube drained rather than ringing, which is what a matched end is for.
        assert!(
            left < start * 0.01,
            "the tube should have drained: {start:e} to {left:e}"
        );
    }

    /// The fourth domain under the scheduler, subcycling to its own CFL limit.
    #[test]
    fn the_domain_runs_under_the_scheduler() {
        use pantometry_core::{Schedule, Simulation};

        let tube = air_tube(201)
            .with_ends(End::Closed, End::Closed)
            .released_from(|x| {
                let u = (x.to_si() - 0.5) / 0.08;
                Pressure::from_si(200.0 * (-u * u).exp())
            });
        let limit = tube.max_stable_dt(Time::ZERO);
        // 5 mm cells in air: 14.6 us.
        assert!((limit.in_us() - 14.58).abs() < 0.05, "{} us", limit.in_us());

        let mut sim = Simulation::new(Schedule::Multirate)
            // The scheme is lossless at any stable Courant number, and `energy` reports
            // what it conserves, so this can be as tight as the arithmetic allows.
            .conservation_tolerance(1e-9)
            .with(tube);

        // A millisecond of audio is 69 substeps.
        let report = sim.advance(Time::ms(1.0)).unwrap();
        assert_eq!(report.substeps[0].0, "tube");
        assert_eq!(report.substeps[0].1, 69);
        for _ in 0..20 {
            sim.advance(Time::ms(1.0)).expect("energy should hold");
        }
    }

    /// A degenerate tube is not a special case a caller has to avoid.
    #[test]
    fn degenerate_tubes_are_handled() {
        // Fewer than three cells is not a tube; the constructor rounds up.
        let tiny = air_tube(1);
        assert_eq!(tiny.cells(), 3);
        // A silent medium has no limit and does nothing.
        let mut frozen = Tube::new(
            "frozen",
            Length::m(1.0),
            11,
            Area::from_si(1e-4),
            Density::kg_per_m3(1.0),
            Velocity::m_per_s(0.0),
        );
        assert!(!frozen.max_stable_dt(Time::ZERO).to_si().is_finite());
        let mut bus = Exchange::new();
        assert!(frozen.step(Time::ZERO, Time::s(1.0), &mut bus).is_ok());
        assert_eq!(frozen.peak_pressure().to_si(), 0.0);
        // A zero step changes nothing.
        let mut quiet = air_tube(21).released_from(|_| Pressure::from_si(5.0));
        let before = quiet.pressure_at(3);
        quiet.step(Time::ZERO, Time::ZERO, &mut bus).unwrap();
        assert_eq!(quiet.pressure_at(3), before);
    }
}
