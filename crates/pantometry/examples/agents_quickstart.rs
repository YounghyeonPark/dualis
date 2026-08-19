//! pantometry in one file, written for a reader with no context and a small budget for acquiring
//! any — an AI agent, most likely, but the shape suits a hurried human too.
//!
//! `AGENTS.md` quotes this file. It is run by CI, so the quickstart cannot drift away from
//! the library the way a hand-written snippet does.
//!
//! ```sh
//! cargo run --example agents_quickstart
//! ```

use pantometry::prelude::*;

// ---------------------------------------------------------------------------------------
// 1. Units are types. A whole class of mistake stops compiling.
// ---------------------------------------------------------------------------------------
//
// There is exactly one place a factor of a thousand may appear: a unit-bearing constructor.
// After that the dimensions travel with the value and the compiler checks the algebra.
//
//     let area: Area = Length::mm(10.0) * Length::mm(10.0);   // fine, Length x Length = Area
//     let wrong = Length::mm(10.0) + Time::s(1.0);            // does not compile
//
// `to_si()` is the only way back out to a bare f64, and it is deliberately noisy to write.

fn units_are_types() {
    let area: Area = Length::mm(10.0) * Length::mm(10.0);
    let absorbed: Power = Irradiance::mw_per_cm2(50.0) * area * 0.02;
    let capacity: HeatCapacity = Mass::g(2.0) * SpecificHeat::j_per_kg_k(858.0);
    let rise: Temperature = (absorbed * Time::s(1.0)) / capacity;

    println!("1. units");
    println!("   50 mW/cm^2 on a 1 cm^2 surface absorbing 2% warms 2 g of glass");
    println!("   by {:.3} mK in a second", rise.to_si() * 1e3);
}

// ---------------------------------------------------------------------------------------
// 2. A domain is anything that steps. Two methods are required; the rest have defaults.
// ---------------------------------------------------------------------------------------
//
// Domains never call each other. They meet on `Exchange`, a bus of named channels carrying
// SI amounts — joules, not watts, because a domain steps over an interval and what crossed
// is an amount. That is what makes the audit an equality rather than an approximation.

/// Publishes joules out of a finite tank.
struct Heater {
    /// Element power.
    watts: f64,
    /// Joules not yet spent. This is the ledger entry, and the reason the books close.
    reserve: f64,
}

impl Domain for Heater {
    fn name(&self) -> &'static str {
        "heater"
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let joules = (self.watts * dt.to_si()).min(self.reserve);
        self.reserve -= joules;
        bus.publish(HEAT, joules);
        Ok(())
    }

    /// **A ledger says what you are holding, not what has passed through you.**
    ///
    /// The joules published are gone from here and are being reported by whoever took
    /// them. Adding them back would make the total grow by the heat that was moved, which
    /// is the single most common way to write a domain that audits green and is wrong.
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.reserve)
    }
}

/// Takes whatever is on the bus and stores it.
struct Slab {
    stored: f64,
    /// Set for the second run below, to show what the kernel does about it.
    lossy: bool,
}

impl Domain for Slab {
    fn name(&self) -> &'static str {
        "slab"
    }

    fn step(&mut self, _t: Time, _dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let arrived = bus.take(HEAT);
        // The bug: every joule is taken off the bus, so nothing is left unclaimed, but only
        // nine tenths of them are booked. The transfer balances and the books do not.
        self.stored += if self.lossy { 0.9 * arrived } else { arrived };
        Ok(())
    }

    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.stored)
    }
}

// ---------------------------------------------------------------------------------------
// 3. The kernel audits, and a wrong model says so by name.
// ---------------------------------------------------------------------------------------
//
// This is the part worth knowing. Writing physics in numpy or a general-purpose engine, a
// wrong model runs happily and produces plausible output; there is no signal to act on.
// Here `advance` returns `Result<Report, Violation>`, and a `Violation` names the quantity,
// the site, what it was, what it became, and against what scale.

fn run(lossy: bool) -> Result<f64, Violation> {
    let mut sim = Simulation::new(Schedule::Staggered)
        // Declaration order is execution order under Staggered, so the publisher goes first
        // and the consumer sees this step's joules rather than the last step's.
        .with(Heater {
            watts: 100.0,
            reserve: 500.0,
        })
        .with(Slab { stored: 0.0, lossy });

    for _ in 0..100 {
        sim.advance(Time::ms(50.0))?;
    }
    // `Simulation::ledger` merges every domain's books; `get` returns None for a quantity
    // nobody reported, which is a different thing from zero and is not conflated with it.
    Ok(sim
        .ledger()
        .get(quantity::ENERGY)
        .expect("both domains report energy"))
}

/// A winding inside a case: the thing that fails is not the thing you can measure.
///
/// `AGENTS.md` quotes this function. The reason it is here rather than written by hand in the
/// document is that CI runs this file, so the snippet an agent copies is a snippet that
/// compiled this morning.
fn junction_to_case() {
    let mut motor = ThermalNetwork::new("motor");
    let winding = motor.node(
        "winding",
        Substance::copper(),
        Volume::from_si(18e-6),
        Length::mm(2.0),
        Temperature::celsius(25.0),
    );
    let case = motor.node_losing_to(
        "case",
        Substance::aluminium_6061(),
        Volume::from_si(220e-6),
        Length::mm(4.0),
        Temperature::celsius(25.0),
        Environment::still_air(Temperature::celsius(25.0), Area::from_si(0.042)),
    );
    motor
        .link(winding, case, Conductance::w_per_k(0.9))
        .expect("two distinct nodes and a positive conductance");
    // Where heat arriving on the bus lands. A network that is never told leaves it unclaimed,
    // which the audit refuses — joules that arrived nowhere are joules that went missing.
    motor.absorbing(winding).expect("winding is a node of this");

    let mut sim = Simulation::new(Schedule::Staggered)
        .with(Heater {
            watts: 6.0,
            reserve: 6_000.0,
        })
        .with(motor);
    for _ in 0..900 {
        sim.advance(Time::s(1.0)).expect("the books close");
    }

    let motor = sim
        .domain_as::<ThermalNetwork>("motor")
        .expect("it is still there");
    let (hot, cold) = (
        motor.node_named("winding").unwrap(),
        motor.node_named("case").unwrap(),
    );
    let drop = motor.temperature(hot).to_si() - motor.temperature(cold).to_si();
    println!(
        "   winding {:.1} C, case {:.1} C — a drop of {drop:.1} K across the joint",
        motor.temperature(hot).to_si() - 273.15,
        motor.temperature(cold).to_si() - 273.15,
    );
    println!("   a LumpedMass would have reported the case's number for both");

    // 6 W across 0.9 W/K is 6.67 K at steady state. Measured 6.2 K at fifteen minutes: the
    // joint itself equilibrates in about 70 s, and the shortfall is the winding's own 62 J/K
    // still filling while the whole assembly warms on its 2000 s constant. Bounded below the
    // steady-state value because exceeding it would mean more power crossing the joint than
    // ever entered.
    assert!(drop > 5.5 && drop < 6.67, "drop {drop:.3} K");

    // And pinned to one decimal, because `AGENTS.md` quotes this number. Nothing here consults
    // a clock or a random source, so the run is reproducible to the bit and a tight pin costs
    // nothing — while the band above would let the documented figure go stale in silence. This
    // repository has shipped stale numbers before; the fix is to make the prose's claim be an
    // assertion somewhere.
    assert!(
        (drop * 10.0).round() == 62.0,
        "AGENTS.md says 6.2 K; this run gives {drop:.3} K. Update both or neither."
    );
}

fn main() {
    units_are_types();

    println!("\n2. a correct pair of domains");
    let total = run(false).expect("the books close");
    // 500 J in the tank at the start, and 500 J between the tank and the slab at the end.
    assert!((total - 500.0).abs() < 1e-9, "total {total}");
    println!("   500.0 J started in the heater's tank");
    println!("   {total:.1} J is still accounted for after five seconds");

    println!("\n3. the same pair, with the slab booking 90% of what it takes");
    match run(true) {
        Ok(total) => panic!("the audit should have refused this: {total}"),
        Err(v) => {
            println!("   advance() returned Err, and this is the whole message:");
            println!("\n     {v}\n");
            // The point for an agent: this is machine-readable, not just printable.
            assert_eq!(v.quantity, quantity::ENERGY);
            println!("   v.quantity = {:?}", v.quantity);
            println!("   v.site     = {:?}", v.site);
            println!("   v.before   = {:.4} J", v.before);
            println!("   v.after    = {:.4} J", v.after);
        }
    }

    println!("\n4. several bodies and the drop between them");
    junction_to_case();

    // ---------------------------------------------------------------------------------
    // What the audit does NOT catch, which matters as much as what it does.
    // ---------------------------------------------------------------------------------
    //
    // It catches quantities appearing or vanishing, amounts left unclaimed on the bus, and
    // fluxes that disagree face by face across a shared boundary. It does not catch a model
    // that is internally consistent and physically wrong: publish the power instead of the
    // energy, forgetting the factor of dt, and the publisher and consumer agree perfectly
    // about a number that is off by 1/dt.
    //
    // That class is caught by checking against something the code did not compute — a closed
    // form, an exact limit, or a convergence rate. Every test in this workspace does one of
    // those, and the five examples do it in public. If you are adding physics here, that is
    // the standard; see CONTRIBUTING.md.
    println!("\nRead next: AGENTS.md for the API surface, examples/ for physics that is");
    println!("checked in public, CONTRIBUTING.md if you are adding to it.");
}
