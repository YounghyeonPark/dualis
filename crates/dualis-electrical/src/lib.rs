//! dualis-electrical: resistive dissipation, as a domain built on the `dualis-core` kernel.
//!
//! The physics that *produces* the watts every other domain in this workspace has so far been
//! handed. `dualis-thermal` consumes heat, `dualis-optics` publishes it from absorbed light,
//! `dualis-mechanics` from a dashpot, `dualis-acoustic` from an absorbing duct end — and a motor
//! or a heater or a trace on a board gets hot for none of those reasons. It gets hot because
//! current went through resistance.
//!
//! Until this crate, the workspace's own examples papered over that with a source of a stated
//! number of watts. That is fine as a stand-in and it is not a model: nothing decides the number,
//! so nothing can be wrong about it.
//!
//! ```
//! use dualis_electrical::Winding;
//! use dualis_units::{Current, Length, Temperature};
//!
//! // A copper winding carrying 3 A. Its resistance is what decides the heat.
//! let coil = Winding::of_copper("coil", Length::m(24.0), 0.35e-6, Temperature::celsius(25.0))
//!     .driven_at(Current::a(3.0));
//! // rho*L/A at 20 C is 1.18217 ohm; five kelvin of copper adds 1.965%.
//! assert!((coil.resistance().to_si() - 1.205401).abs() < 1e-6);
//! assert!((coil.dissipation().to_si() - 10.848610).abs() < 1e-6);   // I^2 R
//! ```
//!
//! # The coupling this crate deliberately does not have
//!
//! Copper's resistivity rises about 0.393% per kelvin, so a winding that gets hot dissipates
//! more, which makes it hotter. That feedback is the whole reason thermal runaway is a thing a
//! designer worries about, and **it is not expressible here.**
//!
//! A domain would need to read another domain's *temperature* inside the step loop.
//! [`Exchange`] carries amounts — joules, coulombs — and not state, which
//! is exactly what makes the conservation audit an equality rather than an approximation. There
//! is no `peek_temperature`, and adding one would not be a small thing: a channel carrying state
//! is not conserved, cannot be audited, and is a short step from domains reading each other,
//! which is the property the crate split exists to hold.
//!
//! So this crate models the resistance at a temperature **you state**, and
//! [`Winding::at_temperature`] is how you state it. The number is right for that temperature and
//! the feedback is the caller's to close, between steps, with both temperatures in hand. What
//! that costs and whether the kernel should grow something is written up rather than decided
//! here — see the repository's `FRICTION.md`.
//!
//! # What the audit can and cannot see
//!
//! A [`Winding`] holds a finite [`reserve`](Winding::with_reserve) of joules, like every other
//! source in this workspace, because a source with an unlimited supply creates energy from
//! nothing every step.
//!
//! **The domain refuses that itself rather than leaving it to the audit**, and the reason is
//! worth stating: an infinite reserve does not fail the audit, it *disables* it. The ledger
//! reports `inf` before and `inf` after, `inf` compares equal to itself, and a winding pouring
//! joules into a plate runs green at any tolerance. This test was written expecting the audit to
//! catch it and it did not.
//!
//! What the audit cannot check is whether `I²R` is the right number of watts. Both sides of the
//! bus agree perfectly about whatever is published, so a resistance wrong by a factor of two
//! balances the books exactly. The tests are therefore against closed forms computed
//! independently: the resistivity of copper at a stated temperature, `P = I²R`, and the exact
//! equivalence of the constant-current and constant-voltage forms at the same operating point.

#![deny(missing_docs)]
#![forbid(unsafe_code)]

use dualis_core::conserved::quantity;
use dualis_core::{Domain, Exchange, Kind, Ledger, Reading, Violation};
use dualis_units::{Current, Length, Power, Resistance, Temperature, Time, Voltage};

/// The channel resistive loss is published on.
///
/// **`quantity::ENERGY`, not a string of this crate's own choosing.** The first version of this
/// crate declared `"heat"`, which reads correctly and is a different channel from the one
/// `dualis-thermal` takes from — so a winding published joules that nothing consumed, in a
/// coupling whose whole point was that the two crates already agree without naming each other.
///
/// The audit caught it on the first step, by name and with the amount: *heat, published but not
/// consumed, 0.045*. That is the library's own argument working on its author, and it is the
/// reason the channel is a shared constant rather than a literal in two places.
pub const HEAT: &str = quantity::ENERGY;

/// Resistivity of annealed copper at 20 °C, in ohm-metres.
///
/// The IACS reference value. Hard-drawn copper runs about 2% higher, and this is the number a
/// wire table is built from.
pub const COPPER_RESISTIVITY_20C: f64 = 1.724e-8;

pub mod conductor;

pub use conductor::Conductor;

/// Temperature coefficient of copper's resistivity, per kelvin, referenced to 20 °C.
///
/// `ρ(T) = ρ₂₀(1 + α(T − 20 °C))`. Linear, which is good to a per cent or so from about −50 °C
/// to 200 °C and is not good at cryogenic temperatures, where the residual resistivity of the
/// particular sample takes over and no coefficient describes it.
pub const COPPER_ALPHA: f64 = 0.00393;

/// How a winding is driven.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Drive {
    /// Constant current: dissipation is `I²R`, and rises with resistance.
    Current(f64),
    /// Constant voltage: dissipation is `V²/R`, and *falls* with resistance.
    Voltage(f64),
}

/// A length of conductor carrying current, dissipating `I²R` onto the heat channel.
///
/// The resistance is computed from the conductor's geometry and its temperature rather than
/// stated, so it is a model rather than a number: get the length or the cross-section wrong and
/// the dissipation is wrong in a way a closed form can catch.
pub struct Winding {
    name: String,
    /// Ohms at the reference temperature, before the coefficient is applied.
    resistance_20c: f64,
    alpha: f64,
    temperature: f64,
    drive: Drive,
    /// Joules left to spend. `f64::INFINITY` until [`Winding::with_reserve`] says otherwise,
    /// and `step` refuses while it is — an infinite ledger entry is not a large number, it is a
    /// number the audit cannot subtract.
    reserve: f64,
    dissipated: f64,
}

impl Winding {
    /// A winding of copper: `length` of wire with a given cross-section, at a temperature.
    ///
    /// `R = ρ(T)·L/A`. The cross-section is in square metres — 0.35 mm² of magnet wire is
    /// `0.35e-6`, which is roughly AWG 22.
    pub fn of_copper(
        name: impl Into<String>,
        length: Length,
        cross_section_m2: f64,
        at: Temperature,
    ) -> Winding {
        let r20 = COPPER_RESISTIVITY_20C * length.to_si() / cross_section_m2.max(f64::MIN_POSITIVE);
        Winding {
            name: name.into(),
            resistance_20c: r20,
            alpha: COPPER_ALPHA,
            temperature: at.to_si(),
            drive: Drive::Current(0.0),
            reserve: f64::INFINITY,
            dissipated: 0.0,
        }
    }

    /// A winding of stated resistance at 20 °C, for a conductor this crate has no geometry for.
    ///
    /// `alpha` is the temperature coefficient per kelvin; pass `0.0` for a resistor whose
    /// coefficient you do not want to model, which is the honest choice for a wirewound part
    /// whose datasheet quotes a tolerance band rather than a number.
    pub fn of_resistance(
        name: impl Into<String>,
        at_20c: Resistance,
        alpha: f64,
        at: Temperature,
    ) -> Winding {
        Winding {
            name: name.into(),
            resistance_20c: at_20c.to_si(),
            alpha,
            temperature: at.to_si(),
            drive: Drive::Current(0.0),
            reserve: f64::INFINITY,
            dissipated: 0.0,
        }
    }

    /// Drive it at a constant current. Dissipation is `I²R` and **rises** as it warms.
    pub fn driven_at(mut self, current: Current) -> Winding {
        self.drive = Drive::Current(current.to_si());
        self
    }

    /// Drive it from a constant voltage. Dissipation is `V²/R` and **falls** as it warms.
    ///
    /// The opposite sign of feedback from [`driven_at`](Winding::driven_at), and the difference
    /// decides whether a runaway is possible at all: a constant-current winding can run away, a
    /// constant-voltage one cannot.
    pub fn driven_from(mut self, voltage: Voltage) -> Winding {
        self.drive = Drive::Voltage(voltage.to_si());
        self
    }

    /// Joules it may dissipate before it goes quiet.
    ///
    /// **Not optional in practice.** Without it the reserve is infinite, and `step` refuses —
    /// see the comment there for why the refusal is the domain's job rather than the audit's.
    /// Saying where the energy comes from is the point.
    pub fn with_reserve(mut self, joules: f64) -> Winding {
        self.reserve = joules.max(0.0);
        self
    }

    /// Tell it what temperature it is now.
    ///
    /// **The manual half of a coupling this crate cannot close by itself.** A domain cannot read
    /// another's temperature inside the step loop — see the module documentation — so a caller
    /// wanting the electro-thermal feedback reads the thermal domain's temperature between
    /// steps and passes it here.
    pub fn at_temperature(&mut self, t: Temperature) {
        self.temperature = t.to_si();
    }

    /// Resistance at its current temperature: `R₂₀(1 + α(T − 20 °C))`.
    pub fn resistance(&self) -> Resistance {
        let dt = self.temperature - Temperature::celsius(20.0).to_si();
        Resistance::from_si(self.resistance_20c * (1.0 + self.alpha * dt))
    }

    /// Power it is dissipating right now.
    pub fn dissipation(&self) -> Power {
        self.dissipation_at(Temperature::from_si(self.temperature))
    }

    /// Resistance at an arbitrary temperature, without changing the winding.
    pub fn resistance_at(&self, at: Temperature) -> Resistance {
        let dt = at.to_si() - Temperature::celsius(20.0).to_si();
        Resistance::from_si(self.resistance_20c * (1.0 + self.alpha * dt))
    }

    /// What it would dissipate at a temperature, without changing the winding.
    ///
    /// **A pure function, and that is the point.** Everything else here is a `Domain` method,
    /// which means it can only be reached by stepping — and the electro-thermal feedback needs
    /// `P(T)` evaluated at a temperature the electrical domain has no way to learn, since
    /// [`Exchange`] carries amounts and not state.
    ///
    /// As a plain function it is composable by whoever *does* hold both sides. A caller with a
    /// thermal network and a winding can write `coil.dissipation_at(net.temperature(node))`
    /// between frames — which is what `dualis-world`'s scene 13 does — and any future in-loop
    /// coupling needs this same function rather than a different one. So it is correct under
    /// every answer to that design question, which is why it exists before the question is
    /// settled.
    pub fn dissipation_at(&self, at: Temperature) -> Power {
        let r = self.resistance_at(at).to_si();
        Power::from_si(match self.drive {
            Drive::Current(i) => i * i * r,
            // A short across an ideal source is not a physical answer, and returning an infinity
            // here would arrive at the audit as a NaN a step later, where its origin is lost.
            Drive::Voltage(v) if r > 0.0 => v * v / r,
            Drive::Voltage(_) => 0.0,
        })
    }

    /// The current at which this winding's feedback overcomes a heat path of conductance `g`.
    ///
    /// Thermal runaway is `dP/dT > dQ_out/dT`. For a constant-current winding `P = I²R₂₀(1+αΔT)`
    /// so `dP/dT = I²R₂₀α`, and against a path that sheds `g` watts per kelvin the threshold is
    /// exact:
    ///
    /// ```text
    ///     I_crit = √( g / (R₂₀ α) )
    /// ```
    ///
    /// `g` is the conductance of the **whole** path to ambient, which for anything with joints
    /// in it is not the surface's. A motor whose winding reaches air through 0.9 W/K and
    /// 2.4 W/K of joints and then 0.294 W/K of convection has a series `g` of 0.203 W/K, not
    /// 0.294 — and the threshold falls from 4.95 A to 4.11 A, a 17% margin a lumped model
    /// reports as present when it is not. That difference is the argument for
    /// [`ThermalNetwork`](https://docs.rs/dualis-thermal) over one body.
    ///
    /// **Do not assemble `g` by hand.** Those three numbers are a convection-only path, and the
    /// real one on that motor is 0.220 W/K because the housing also radiates at its operating
    /// temperature — so the hand-computed 4.11 A understates the true 4.28 A by 4%. Ask the
    /// network: `ThermalNetwork::path_conductance(node, at)` takes the slope of its own solved
    /// balance and therefore includes every path out, radiative terms and interior environments
    /// alike. The hand formula was found to be wrong by a sizing tool written against 0.6.0,
    /// which is what `FRICTION.md` 20 is about.
    ///
    /// Returns `None` for a voltage-driven winding, which cannot run away: `P = V²/R` *falls*
    /// as it warms, so the feedback has the opposite sign and there is no threshold to report.
    pub fn runaway_current(&self, path: dualis_units::Conductance) -> Option<Current> {
        match self.drive {
            Drive::Voltage(_) => None,
            Drive::Current(_) => {
                let denom = self.resistance_20c * self.alpha;
                if denom <= 0.0 || !path.to_si().is_finite() || path.to_si() <= 0.0 {
                    return None;
                }
                Some(Current::from_si((path.to_si() / denom).sqrt()))
            }
        }
    }

    /// The current through it, whichever way it is driven.
    pub fn current(&self) -> Current {
        let r = self.resistance().to_si();
        Current::from_si(match self.drive {
            Drive::Current(i) => i,
            Drive::Voltage(v) if r > 0.0 => v / r,
            Drive::Voltage(_) => 0.0,
        })
    }

    /// The voltage across it, whichever way it is driven.
    pub fn voltage(&self) -> Voltage {
        Voltage::from_si(match self.drive {
            Drive::Current(i) => i * self.resistance().to_si(),
            Drive::Voltage(v) => v,
        })
    }

    /// Joules it has put onto the bus over the run.
    pub fn dissipated_energy(&self) -> dualis_units::Energy {
        dualis_units::Energy::from_si(self.dissipated)
    }

    /// Joules it has left to spend.
    pub fn reserve(&self) -> dualis_units::Energy {
        dualis_units::Energy::from_si(self.reserve)
    }
}

impl Domain for Winding {
    fn name(&self) -> &str {
        &self.name
    }

    /// Nothing here is integrated: the dissipation is a closed-form function of the state, so
    /// there is no state to march and no stability limit to respect.
    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        // An infinite reserve does not merely go uncaught — it turns the audit *off*. The ledger
        // reports `inf` before and `inf` after, `inf` compares equal to itself, and a winding
        // pouring joules into a plate produces a run that is green at any tolerance. This
        // repository has shipped that failure once already, in `dualis-world`'s lamp, where a
        // skipped `with_reserve` left the reserve infinite and a scene audited clean at
        // tolerance zero with the lamp doing nothing.
        //
        // So the domain refuses rather than trusting the audit to notice, because the audit
        // structurally cannot.
        let watts = self.dissipation().to_si();
        if !watts.is_finite() || watts < 0.0 {
            return Err(Violation::at(
                &self.name,
                "dissipation is not a power",
                watts,
            ));
        }
        // Refused when it would actually publish, rather than on construction: a winding sitting
        // at zero current puts nothing on the bus and blinds nothing, so demanding a reserve of
        // it would be a rule about the API instead of about the physics. The drive is fixed at
        // construction — `driven_at` and `driven_from` consume `self` — so a winding that
        // publishes nothing on the first step publishes nothing ever, and this cannot be dodged
        // by starting at zero.
        if watts > 0.0 && !self.reserve.is_finite() {
            return Err(Violation::at(
                &self.name,
                "no reserve was set, so this winding supplies energy from nowhere and the audit \
                 cannot see it: an infinite ledger entry is not a large number, it is one that \
                 cannot be subtracted. Call with_reserve",
                watts,
            ));
        }
        let joules = (watts * dt.to_si()).min(self.reserve).max(0.0);
        self.reserve -= joules;
        self.dissipated += joules;
        bus.publish(HEAT, joules);
        Ok(())
    }

    /// **What is left to spend, not what has passed through.**
    ///
    /// The joules published are gone from here and are being reported by whoever took them.
    /// Adding them back would make the total grow by the heat that moved, which this workspace
    /// has written down more than once as the most common way to author a domain that audits
    /// green and is wrong.
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.reserve)
    }

    /// What it is dissipating, at what resistance, and what it has spent.
    ///
    /// The resistance is here because it is the *reason* the dissipation moves: at fixed current
    /// the power ratio is the resistance ratio, and a table with both columns shows that in a way
    /// a table with one cannot.
    fn readings(&self) -> Vec<Reading> {
        vec![
            Reading::new(&self.name, "dissipating", self.dissipation().to_si(), "W"),
            Reading::new(&self.name, "resistance", self.resistance().to_si(), "ohm"),
            Reading::new(&self.name, "spent", self.dissipated.max(0.0), "J"),
        ]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }
}
