//! dualis in one file, written for a reader with no context and a small budget for acquiring
//! any — an AI agent, most likely, but the shape suits a hurried human too.
//!
//! `AGENTS.md` quotes this file. It is run by CI, so the quickstart cannot drift away from
//! the library the way a hand-written snippet does.
//!
//! ```sh
//! cargo run --example agents_quickstart
//! ```

use dualis::prelude::*;

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
    /// Joules not yet spent. This is the ledger entry, and the reason the books close.
    reserve: f64,
}

impl Domain for Heater {
    fn name(&self) -> &'static str {
        "heater"
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let joules = (100.0 * dt.to_si()).min(self.reserve); // a 100 W element
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
        .with(Heater { reserve: 500.0 })
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
