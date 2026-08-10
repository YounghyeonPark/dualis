//! A small domain's leak, hidden inside a large domain's total.
//!
//! The whole-simulation audit sums every ledger before comparing, so the scale it measures against
//! is the *sum*. A domain holding a microjoule beside one holding a kilojoule can lose everything
//! it has without moving that sum by anything a tolerance could catch — and no tolerance fixes it,
//! because the problem is the scale rather than the number.
//!
//! `Domain::books_balance` is the opt-in that says "check me on my own", and this file is the
//! demonstration that it makes a difference.

use dualis_core::conserved::quantity;
use dualis_core::units::Time;
use dualis_core::{Domain, Exchange, Kind, Ledger, Schedule, Simulation, Violation};

/// A domain holding energy, optionally leaking a fraction of it each step, optionally claiming
/// its books are exact.
struct Holder {
    name: &'static str,
    energy: f64,
    leak: f64,
    claims_exact: bool,
}

impl Domain for Holder {
    fn name(&self) -> &str {
        self.name
    }
    fn kind(&self) -> Kind {
        Kind::Evolving
    }
    fn step(&mut self, _t: Time, _dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        self.energy *= 1.0 - self.leak;
        Ok(())
    }
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.energy)
    }
    fn books_balance(&self) -> bool {
        self.claims_exact
    }
}

fn holder(name: &'static str, energy: f64, leak: f64, claims_exact: bool) -> Holder {
    Holder {
        name,
        energy,
        leak,
        claims_exact,
    }
}

/// **The failure: a small domain loses a fifth of itself and the sum does not notice.**
///
/// A kilojoule beside a microjoule. The small one leaks 20% a step — catastrophic, obvious, the
/// kind of bug a conservation audit exists for — and the *total* moves by 2×10⁻¹⁰, which is inside
/// even a `1e-9` tolerance.
///
/// This is not a tolerance that was set too loosely. Tightening it to `1e-12` would refuse the run
/// for floating-point noise long before it could see a leak of this shape.
#[test]
fn a_small_domains_leak_hides_inside_a_large_domains_total() {
    let mut sim = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .with(holder("big", 1000.0, 0.0, false))
        .with(holder("small", 1e-6, 0.2, false));

    sim.advance(Time::from_si(1.0))
        .expect("the sum barely moves, so the audit sees nothing");

    let after = sim
        .domain("small")
        .expect("still there")
        .ledger()
        .get(quantity::ENERGY)
        .unwrap_or(0.0);
    assert!(
        (after / 1e-6 - 0.8).abs() < 1e-12,
        "the small domain really did lose a fifth: {after:.6e}"
    );

    // And the reason the audit could not see it, stated as arithmetic: the leak is a fifth of a
    // microjoule against a total of a kilojoule.
    let relative = 0.2 * 1e-6 / 1000.0;
    assert!(
        relative < 1e-9,
        "{relative:.3e} is inside any tolerance a real simulation would use"
    );
}

/// **The same leak, caught, because the domain claims its books are exact.**
///
/// Nothing else changes: the same two domains, the same sizes, the same tolerance. The small one
/// says `books_balance`, so it is checked against **its own** holdings rather than against the sum,
/// and a fifth is a fifth at any scale.
#[test]
fn claiming_exact_books_catches_it() {
    let mut sim = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .with(holder("big", 1000.0, 0.0, false))
        .with(holder("small", 1e-6, 0.2, true));

    let caught = sim
        .advance(Time::from_si(1.0))
        .expect_err("a domain checked on its own scale cannot lose a fifth quietly");
    assert_eq!(caught.quantity, "energy");
    assert!(
        caught.site.starts_with("small"),
        "the violation should name the domain that leaked, said {:?}",
        caught.site
    );
    // And the numbers are the domain's own, not the simulation's.
    assert!(
        (caught.scale - 1e-6).abs() < 1e-18,
        "the scale should be the small domain's, was {:.3e}",
        caught.scale
    );
}

/// **A domain that legitimately exchanges is not accused of leaking.**
///
/// The check is not "your ledger did not change" — a consumer's ledger is *supposed* to change.
/// It is "your ledger changed by exactly what you took minus what you published", which is a much
/// sharper statement and the reason this is worth doing at all.
#[test]
fn a_domain_that_takes_from_the_bus_is_measured_against_what_it_took() {
    struct Source {
        left: f64,
    }
    impl Domain for Source {
        fn name(&self) -> &str {
            "source"
        }
        fn kind(&self) -> Kind {
            Kind::QuasiStatic
        }
        fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
            let joules = (4.0 * dt.to_si()).min(self.left);
            self.left -= joules;
            bus.publish(quantity::ENERGY, joules);
            Ok(())
        }
        fn ledger(&self) -> Ledger {
            Ledger::new().with(quantity::ENERGY, self.left)
        }
        fn books_balance(&self) -> bool {
            true
        }
    }

    struct Consumer {
        held: f64,
        /// A fraction of what it takes that never reaches its books.
        skim: f64,
    }
    impl Domain for Consumer {
        fn name(&self) -> &str {
            "consumer"
        }
        fn kind(&self) -> Kind {
            Kind::Evolving
        }
        fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
            self.held += bus.take_share(quantity::ENERGY, dt) * (1.0 - self.skim);
            Ok(())
        }
        fn ledger(&self) -> Ledger {
            Ledger::new().with(quantity::ENERGY, self.held)
        }
        fn books_balance(&self) -> bool {
            true
        }
    }

    // Honest: everything taken is kept, and both domains' books balance against the bus.
    let mut sim = Simulation::new(Schedule::Staggered)
        .with(Source { left: 100.0 })
        .with(Consumer {
            held: 0.0,
            skim: 0.0,
        });
    sim.advance(Time::from_si(1.0))
        .expect("a source paying and a consumer keeping is exactly balanced books");
    let held = sim
        .domain("consumer")
        .unwrap()
        .ledger()
        .get(quantity::ENERGY)
        .unwrap_or(0.0);
    assert!((held - 4.0).abs() < 1e-12, "4 J moved, got {held}");

    // Dishonest by one part in ten thousand: it takes and does not record all of it. The
    // whole-simulation audit *would* also catch this one, because both domains are the same size
    // — the point here is that the per-domain check names the culprit rather than the total.
    let mut sim = Simulation::new(Schedule::Staggered)
        .with(Source { left: 100.0 })
        .with(Consumer {
            held: 0.0,
            skim: 1e-4,
        });
    let caught = sim
        .advance(Time::from_si(1.0))
        .expect_err("a consumer that drops part of what it took has not balanced its books");
    assert!(
        caught.site.starts_with("consumer"),
        "it should name the consumer, said {:?}",
        caught.site
    );
}

/// **A domain that does not claim exact books is left alone, and the per-quantity tolerance
/// applies to the ones that do.**
///
/// Two things at once, because they are the same mechanism. The opt-in has to cost nothing to
/// decline, or every domain modelling a boundary the bus does not carry — a lumped mass losing
/// heat to still air — would be accused of leaking. And the number a claiming domain is held to
/// is the *quantity's* tolerance, so the two halves of this change compose rather than fight.
///
/// The pairing with a large domain is what keeps the whole-simulation audit out of it: at these
/// sizes the sum does not move, so anything that fires here fired because of the claim.
#[test]
fn declining_the_claim_costs_nothing_and_claiming_uses_the_quantity_tolerance() {
    let paired = |claims: bool| {
        Simulation::new(Schedule::Staggered)
            .conservation_tolerance(1e-9)
            .with(holder("big", 1000.0, 0.0, false))
            .with(holder("radiating", 1e-6, 0.05, claims))
    };

    paired(false)
        .advance(Time::from_si(1.0))
        .expect("a domain that made no claim is not held to one");
    paired(true)
        .advance(Time::from_si(1.0))
        .expect_err("and one that made it is");

    // The per-quantity override reaches the per-domain check too: energy at 10% allows a 5% leak
    // that energy at 1e-9 refuses. Same domain, same claim, same run — only the number differs.
    Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .conservation_tolerance_for(quantity::ENERGY, 0.1)
        .with(holder("big", 1000.0, 0.0, false))
        .with(holder("radiating", 1e-6, 0.05, true))
        .advance(Time::from_si(1.0))
        .expect("energy at 10% allows a 5% leak, however small the domain");

    // And a tolerance on a *different* quantity does not silence it, which is the mistake a
    // single shared number would make.
    Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .conservation_tolerance_for(quantity::MOMENTUM, 0.1)
        .with(holder("big", 1000.0, 0.0, false))
        .with(holder("radiating", 1e-6, 0.05, true))
        .advance(Time::from_si(1.0))
        .expect_err("loosening momentum must not loosen energy");
}
