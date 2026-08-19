//! What the bus's second-consumer check must and must not refuse.
//!
//! `two_consumers.rs` pins the two endpoints — a robbery is refused, a subcycling consumer is
//! not. This file pins the arithmetic between them, because a review that mutated the check
//! found the middle unguarded: deleting its per-sweep datum, restricting it to one channel, or
//! putting an unearned floor on it all left the suite green.
//!
//! Every scene here is small and synthetic on purpose. The claims are about the *bus*, and a
//! real domain would bring its own physics into the answer.

use pantometry::core::conserved::quantity;
use pantometry::core::{Exchange, Kind, Ledger, Violation};
use pantometry::prelude::*;

/// A domain that publishes a fixed amount each step, takes a fixed number of times, or both.
///
/// Deliberately not a physics: what is under test is the bus's bookkeeping, and a `LumpedMass`
/// would answer with a temperature as well.
struct Part {
    name: String,
    publishes: Vec<(&'static str, f64)>,
    takes: Vec<&'static str>,
    held: f64,
}

impl Part {
    fn new(name: &str) -> Part {
        Part {
            name: name.to_string(),
            publishes: Vec::new(),
            takes: Vec::new(),
            held: 0.0,
        }
    }
    fn publishing(mut self, channel: &'static str, amount: f64) -> Part {
        self.publishes.push((channel, amount));
        self
    }
    fn taking(mut self, channel: &'static str) -> Part {
        self.takes.push(channel);
        self
    }
}

impl Domain for Part {
    fn name(&self) -> &str {
        &self.name
    }
    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }
    fn step(&mut self, _t: Time, _dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        for (channel, amount) in &self.publishes {
            bus.publish(channel, *amount);
        }
        for channel in &self.takes {
            self.held += bus.take(channel);
        }
        Ok(())
    }
    fn ledger(&self) -> Ledger {
        // Everything published here comes out of an unmodelled reserve, so no ledger: these
        // scenes are about the take check and not about conservation.
        Ledger::new()
    }
}

fn run(parts: Vec<Part>) -> Result<(), Violation> {
    let mut sim = Simulation::new(Schedule::Staggered);
    for p in parts {
        sim = sim.with(p);
    }
    sim.advance(Time::s(1.0)).map(|_| ())
}

/// **Two takers of a channel nobody fed are not a defect.** The check counted takes once, so an
/// empty channel taken from twice was refused exactly as a full one — and two blocks side by
/// side with no source is the first scene anybody assembling parts writes.
#[test]
fn two_takers_of_an_empty_channel_are_fine() {
    run(vec![
        Part::new("a").taking(quantity::ENERGY),
        Part::new("b").taking(quantity::ENERGY),
    ])
    .expect("nothing was published, so nothing could have been mis-split");
}

/// **A second taker that finds an emptied channel is refused, and named.** The case the check
/// exists for, restated here at the bus level.
#[test]
fn a_second_taker_of_a_fed_channel_is_refused() {
    let why = run(vec![
        Part::new("source").publishing(quantity::ENERGY, 5.0),
        Part::new("first").taking(quantity::ENERGY),
        Part::new("second").taking(quantity::ENERGY),
    ])
    .expect_err("the second taker found nothing");
    assert!(why.site.starts_with("second"), "{}", why.site);
    assert_eq!(why.before, 5.0, "the message should carry what moved");
    assert_eq!(why.after, 0.0, "and what this domain got");
}

/// **A producer between two consumers is allowed, and this is the narrower promise stated.**
/// Both consumers received a real amount and nothing went missing; which arrangement was meant
/// cannot be read from a bus that carries amounts and an order.
#[test]
fn a_producer_between_two_consumers_is_allowed() {
    run(vec![
        Part::new("source").publishing(quantity::ENERGY, 5.0),
        Part::new("first").taking(quantity::ENERGY),
        Part::new("again").publishing(quantity::ENERGY, 3.0),
        Part::new("second").taking(quantity::ENERGY),
    ])
    .expect("both consumers received a real amount");
}

/// **A domain that publishes and takes its own offer back received nothing, and cannot use that
/// to mask a robbery.**
///
/// The check reads the **net** — `taken − published` — which is the quantity `attribute` uses a
/// few lines away in the same function. Counting the gross let this scene pass: `patch` takes
/// its own token back, its gross is nonzero, and the lamp's five joules quietly went to `near`
/// alone.
#[test]
fn a_self_publisher_cannot_mask_a_robbery() {
    let why = run(vec![
        Part::new("lamp").publishing(quantity::ENERGY, 5.0),
        Part::new("near").taking(quantity::ENERGY),
        Part::new("patch")
            .publishing(quantity::ENERGY, 0.001)
            .taking(quantity::ENERGY),
    ])
    .expect_err("patch received its own offer and nothing else");
    assert!(why.site.starts_with("patch"), "{}", why.site);
}

/// **And the mirror: a self-publisher's own receipts are not somebody else's loss.** Counting
/// the gross made `bystander` the accused for a channel `busy` had only cycled through itself.
#[test]
fn a_self_publisher_does_not_accuse_the_next_taker() {
    run(vec![
        Part::new("busy")
            .publishing(quantity::ENERGY, 3.0)
            .taking(quantity::ENERGY),
        Part::new("bystander").taking(quantity::ENERGY),
    ])
    .expect("busy moved nothing on net, so bystander lost nothing");
}

/// **A negative offer still moves something.** `publish` places no sign restriction, and a
/// domain modelling a net exchange flips sign with its state. Comparing against zero rather
/// than against a magnitude let a chiller's five joules leave the second taker empty in silence.
#[test]
fn a_negative_transfer_is_still_a_transfer() {
    let why = run(vec![
        Part::new("chiller").publishing(quantity::ENERGY, -5.0),
        Part::new("near").taking(quantity::ENERGY),
        Part::new("far").taking(quantity::ENERGY),
    ])
    .expect_err("far found the channel empty");
    assert!(why.site.starts_with("far"), "{}", why.site);
    assert_eq!(why.before, 5.0, "the magnitude of what moved");
}

/// **Two earlier takers whose receipts cancel still moved something**, and must not disarm the
/// check for everybody after them. Ten joules crossed; a net of zero says otherwise.
#[test]
fn cancelling_receipts_do_not_disarm_the_check() {
    let why = run(vec![
        Part::new("up").publishing(quantity::ENERGY, 5.0),
        Part::new("a").taking(quantity::ENERGY),
        Part::new("down").publishing(quantity::ENERGY, -5.0),
        Part::new("b").taking(quantity::ENERGY),
        Part::new("c").taking(quantity::ENERGY),
    ])
    .expect_err("c found the channel empty after ten joules had crossed");
    assert!(why.site.starts_with("c"), "{}", why.site);
}

/// **Every channel is checked, not the first one in name order.**
///
/// A mutation that inspected only one channel survived the whole suite, because no scene had a
/// robbery on any channel but the alphabetically first. Here `energy` is clean and `mass` is
/// not, so a check that stopped early would pass this and should not.
#[test]
fn a_robbery_on_the_second_channel_is_seen() {
    let why = run(vec![
        Part::new("source").publishing(quantity::MASS, 2.0),
        Part::new("first")
            .taking(quantity::ENERGY)
            .taking(quantity::MASS),
        Part::new("second")
            .taking(quantity::ENERGY)
            .taking(quantity::MASS),
    ])
    .expect_err("the mass channel was emptied before the second taker");
    assert_eq!(why.quantity, quantity::MASS);
    assert!(why.site.starts_with("second"), "{}", why.site);
}

/// **The check's sensitivity is the sweep's own transfers, not the run's history.**
///
/// Reading the bus's running totals made the smallest visible mis-split `2⁻⁵²` of everything the
/// channel had ever carried: after a gigajoule, a microjoule robbery was invisible, and the same
/// scene passed or failed according to how long it had been running. Here a large history is
/// built first and a small robbery is committed after it.
#[test]
fn a_long_history_does_not_blind_the_check() {
    let mut sim = Simulation::new(Schedule::Staggered)
        .with(Part::new("bulk").publishing(quantity::ENERGY, 1e9))
        .with(Part::new("sink").taking(quantity::ENERGY));
    for _ in 0..4 {
        sim.advance(Time::s(1.0)).expect("one taker, no defect");
    }

    // Now the same bus, with a robbery four orders below the rounding of that history.
    let mut sim = Simulation::new(Schedule::Staggered)
        .with(Part::new("bulk").publishing(quantity::ENERGY, 1e9))
        .with(Part::new("sink").taking(quantity::ENERGY))
        .with(Part::new("tiny").publishing(quantity::ENERGY, 1e-8))
        .with(Part::new("first").taking(quantity::ENERGY))
        .with(Part::new("second").taking(quantity::ENERGY));
    for _ in 0..4 {
        let _ = sim.advance(Time::s(1.0));
    }
    let why = sim
        .advance(Time::s(1.0))
        .expect_err("a microjoule robbery must survive a gigajoule of history");
    assert!(why.site.starts_with("second"), "{}", why.site);
}

/// **A taker that only starts taking partway through a run is not accused of a past sweep's
/// traffic.**
///
/// The scene that flips when the check's per-sweep datum is removed: a source with a reserve
/// that gives out, a consumer that always takes, and a second that joins later — a part that
/// only draws while it is in contact. Reading the run's totals instead of the sweep's makes the
/// later one look like a thief for what the first received on an earlier step.
#[test]
fn a_taker_that_joins_later_is_not_accused_of_an_earlier_sweep() {
    /// Publishes on the first step only, like a reserve that gives out.
    struct Once {
        spent: bool,
    }
    impl Domain for Once {
        fn name(&self) -> &str {
            "once"
        }
        fn kind(&self) -> Kind {
            Kind::QuasiStatic
        }
        fn step(&mut self, _t: Time, _dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
            if !self.spent {
                bus.publish(quantity::ENERGY, 5.0);
                self.spent = true;
            }
            Ok(())
        }
    }

    let mut sim = Simulation::new(Schedule::Staggered)
        .with(Once { spent: false })
        .with(Part::new("always").taking(quantity::ENERGY))
        .with(Part::new("later").taking(quantity::ENERGY));

    // The first advance is the one with something to share, and it is refused — correctly.
    sim.advance(Time::s(1.0))
        .expect_err("both takers want the one offer");

    // With nothing left to share, the same two coexist for as long as the run lasts.
    let mut sim = Simulation::new(Schedule::Staggered)
        .with(Once { spent: true })
        .with(Part::new("always").taking(quantity::ENERGY))
        .with(Part::new("later").taking(quantity::ENERGY));
    for step in 0..3 {
        sim.advance(Time::s(1.0))
            .unwrap_or_else(|v| panic!("step {step} refused with nothing on the channel: {v}"));
    }
}
