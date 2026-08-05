//! Running several domains at once.
//!
//! A domain is a piece of physics that can be stepped: heat in a block of glass,
//! a rigid body under contact, light through a train of surfaces. Each one knows
//! its own equations and nothing about the others. This module is how they share
//! a clock and a budget without knowing about each other.
//!
//! # The timescale problem, which is the real one
//!
//! Domains do not agree on how big a step is. An explicit FDTD electromagnetic
//! solver on a nanometre grid is stable to about 10⁻¹⁷ s; heat conduction to about
//! 10⁻⁹ s; rigid contact to 10⁻⁴ s; and a thermal drift that defocuses an
//! instrument plays out over seconds. Stepping all of them at the smallest limit
//! integrates the slow ones ten billion times for nothing.
//!
//! Two mechanisms deal with that, and they are the reason this module is not just
//! a `for` loop over domains:
//!
//! - **[`Kind::QuasiStatic`]** — a domain with no state to roll forward, which is
//!   re-solved on demand instead of stepped. Light crosses an instrument in
//!   nanoseconds; against a thermal timescale that is zero, so optics is not
//!   integrated at all. This is the largest single saving available, and it is
//!   what the closed-form [`Motion`](crate::motion::Motion) and the instantaneous
//!   `SurfaceOptics` were already doing before there was a scheduler to notice.
//! - **[`Schedule::Multirate`]** — each evolving domain takes as many equal
//!   substeps of the shared window as its own stability limit requires, so the
//!   slow domain is not dragged down to the fast one's step.
//!
//! # Coupling, and why it goes through a bus
//!
//! Domains never touch each other. They publish to and consume from an
//! [`Exchange`], which is a set of named channels carrying SI amounts. That is not
//! only a borrow-checker convenience: it is what makes the transfer *auditable*.
//! Each domain conserves energy internally, but the interface between two
//! discretisations of the same surface — ray hits on one side, mesh nodes on the
//! other — is exactly where interpolation quietly loses or invents some. The bus
//! compares what was published against what was consumed and refuses to let the
//! difference pass silently.
//!
//! # What the schedules cost
//!
//! [`Schedule::OneWay`] is unconditionally stable and embarrassingly parallel,
//! because nothing feeds back. [`Schedule::Staggered`] costs one exchange per
//! step and is stable only while the coupling is weak — and *not* fixable by
//! shrinking `dt`, since some strongly coupled systems (the standard example is
//! fluid-structure interaction at comparable densities, the added-mass effect)
//! become more unstable as the step shrinks. That is what
//! [`Schedule::Iterative`] is for, and why it is worth its cost.

use std::collections::BTreeMap;

use dualis_units::Time;

use crate::conserved::{audit, Ledger, Violation};
use crate::integrator::substeps_for;

/// Whether a domain has state to roll forward.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// Has state, and a stability limit on how far it can be stepped at once.
    Evolving,
    /// Has no state: solved from its inputs whenever asked, in zero time. Optics,
    /// a static load, an equilibrium reaction. Never subcycled — a solve is a
    /// solve.
    QuasiStatic,
}

/// One piece of physics.
///
/// The only required methods are the name and the step; the rest have defaults
/// that describe a well-behaved evolving domain with no stability limit and no
/// books to keep.
pub trait Domain {
    fn name(&self) -> &'static str;

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// The largest step this domain can take from `now` and stay stable — a CFL
    /// condition, a diffusion limit, a contact penetration budget.
    ///
    /// Infinite means "no limit", which is the honest answer for a quasi-static
    /// domain and for a linear one being solved implicitly.
    fn max_stable_dt(&self, now: Time) -> Time {
        let _ = now;
        Time::from_si(f64::INFINITY)
    }

    /// Advance by `dt` from `t`, reading inputs from `bus` and publishing outputs
    /// to it. A quasi-static domain ignores `dt`.
    ///
    /// Must be a pure function of its state and its inputs: no wall clock, no
    /// unordered reduction, no shared generator. [`Rng::for_index`](crate::Rng::for_index)
    /// is how a domain gets randomness without giving that up.
    fn step(&mut self, t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation>;

    /// How far this domain still is from agreeing with its neighbours, for
    /// [`Schedule::Iterative`]. Zero means converged.
    fn residual(&self) -> f64 {
        0.0
    }

    /// What this domain is holding, for the conservation audit.
    fn ledger(&self) -> Ledger {
        Ledger::new()
    }

    /// Save state so an iterative sweep can be re-run from the same starting
    /// point. A domain that does not implement this cannot take part in
    /// [`Schedule::Iterative`], and [`Simulation::advance`] says so rather than
    /// silently iterating from the wrong state.
    fn checkpoint(&mut self) {}

    /// Restore the last [`Domain::checkpoint`].
    fn restore(&mut self) {}

    fn supports_restore(&self) -> bool {
        false
    }
}

/// The channel between domains: named quantities, in SI base units.
///
/// A domain publishes what it produced and consumes what it needs. Nothing else
/// crosses between domains, which means every transfer is in one place and can be
/// checked in one place.
#[derive(Clone, Debug, Default)]
pub struct Exchange {
    published: BTreeMap<&'static str, f64>,
    consumed: BTreeMap<&'static str, f64>,
}

impl Exchange {
    pub fn new() -> Exchange {
        Exchange::default()
    }

    /// Offer an amount on a channel. Repeated publishes accumulate, so several
    /// surfaces can each contribute to one heat load.
    pub fn publish(&mut self, channel: &'static str, si_amount: f64) {
        *self.published.entry(channel).or_insert(0.0) += si_amount;
    }

    /// Take everything on a channel, recording that it was taken. The channel is
    /// left empty: an amount consumed twice would be an amount doubled.
    pub fn take(&mut self, channel: &'static str) -> f64 {
        let amount = self.published.insert(channel, 0.0).unwrap_or(0.0);
        *self.consumed.entry(channel).or_insert(0.0) += amount;
        amount
    }

    /// Look without taking.
    pub fn peek(&self, channel: &'static str) -> f64 {
        self.published.get(channel).copied().unwrap_or(0.0)
    }

    /// Channels that were published to but never taken from, with what is left on
    /// them. Energy sitting here at the end of a step is energy that left one
    /// domain and arrived nowhere.
    pub fn unclaimed(&self) -> impl Iterator<Item = (&'static str, f64)> + '_ {
        self.published
            .iter()
            .filter(|(_, v)| v.abs() > 0.0)
            .map(|(k, v)| (*k, *v))
    }

    /// Fail if anything published was not consumed.
    ///
    /// This is the check that catches a coupling whose two sides disagree — a
    /// surface that absorbed 3.7 mW handing it to a mesh that received 3.4 mW
    /// because the interpolation between their discretisations lost the rest.
    pub fn audit_transfers(&self, site: &str, abs_tol: f64) -> Result<(), Violation> {
        for (channel, left) in self.published.iter() {
            if left.abs() > abs_tol {
                return Err(Violation {
                    quantity: (*channel).to_string(),
                    site: format!("{site} (published but not consumed)"),
                    before: *left,
                    after: 0.0,
                    // An absolute check: the amount left on the channel *is* the
                    // scale, because all of it went missing.
                    scale: left.abs(),
                    tolerance: abs_tol,
                });
            }
        }
        Ok(())
    }

    /// Total taken from a channel over the run, for reporting.
    pub fn total_consumed(&self, channel: &str) -> f64 {
        self.consumed.get(channel).copied().unwrap_or(0.0)
    }

    /// Empty the offers, keeping the running consumption totals.
    pub fn clear_offers(&mut self) {
        self.published.clear();
    }
}

/// How the domains are interleaved.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Schedule {
    /// One pass in declared order, no feedback expected. Unconditionally stable;
    /// the only schedule whose domains could safely run concurrently.
    OneWay,
    /// One pass in declared order, with each domain seeing the previous ones'
    /// output from this step and the later ones' from the last. Cheap, and stable
    /// only while the coupling is weak.
    Staggered,
    /// Repeat the pass until every domain's residual is under `tol`, or fail.
    ///
    /// The cost is `max_iter` passes; the benefit is stability where a staggered
    /// scheme diverges no matter how small the step. Failing to converge is
    /// reported as a [`Violation`] rather than accepted, because an unconverged
    /// coupling that is allowed through is the most expensive kind of wrong
    /// answer: it looks like physics.
    Iterative { max_iter: u32, tol: f64 },
    /// As [`Schedule::Staggered`], but each evolving domain takes as many equal
    /// substeps as its own stability limit needs.
    Multirate,
}

/// What one [`Simulation::advance`] actually did.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Report {
    /// Substeps taken, per domain, in declared order.
    pub substeps: Vec<(&'static str, u32)>,
    /// Coupling iterations used. One for every schedule but `Iterative`.
    pub iterations: u32,
    /// Largest residual left at the end.
    pub residual: f64,
}

/// A set of domains sharing a clock.
pub struct Simulation {
    domains: Vec<Box<dyn Domain>>,
    schedule: Schedule,
    bus: Exchange,
    t: Time,
    transfer_tol: f64,
    conservation_tol: f64,
}

impl Simulation {
    /// Domains are stepped in the order they are added. That order is part of the
    /// physics under a staggered schedule — put the quasi-static producers before
    /// the evolving consumers — and it is fixed rather than discovered, so two
    /// runs take the same path.
    pub fn new(schedule: Schedule) -> Simulation {
        Simulation {
            domains: Vec::new(),
            schedule,
            bus: Exchange::new(),
            t: Time::ZERO,
            transfer_tol: 1e-12,
            conservation_tol: 1e-9,
        }
    }

    pub fn with(mut self, domain: impl Domain + 'static) -> Simulation {
        self.domains.push(Box::new(domain));
        self
    }

    /// Absolute tolerance on the bus audit, in SI units of whatever is on the
    /// channel. Default 1e-12.
    pub fn transfer_tolerance(mut self, tol: f64) -> Simulation {
        self.transfer_tol = tol;
        self
    }

    /// Relative tolerance on the whole-simulation conservation audit across a
    /// step. Default 1e-9.
    pub fn conservation_tolerance(mut self, tol: f64) -> Simulation {
        self.conservation_tol = tol;
        self
    }

    pub fn time(&self) -> Time {
        self.t
    }

    pub fn bus(&self) -> &Exchange {
        &self.bus
    }

    pub fn domain(&self, name: &str) -> Option<&dyn Domain> {
        self.domains
            .iter()
            .find(|d| d.name() == name)
            .map(|d| d.as_ref())
    }

    /// Every domain's books, summed.
    pub fn ledger(&self) -> Ledger {
        self.domains
            .iter()
            .fold(Ledger::new(), |total, d| total.merged(&d.ledger()))
    }

    /// Advance every domain by `dt`.
    ///
    /// Fails without advancing the clock if a domain fails, if the bus does not
    /// balance, if an iterative coupling does not converge, or if the totalled
    /// ledgers moved by more than the conservation tolerance.
    pub fn advance(&mut self, dt: Time) -> Result<Report, Violation> {
        let before = self.ledger();
        let report = match self.schedule {
            Schedule::OneWay | Schedule::Staggered => self.sweep(dt, false)?,
            Schedule::Multirate => self.sweep(dt, true)?,
            Schedule::Iterative { max_iter, tol } => self.iterate(dt, max_iter, tol)?,
        };

        self.bus.audit_transfers("bus", self.transfer_tol)?;
        let after = self.ledger();
        if !before.is_empty() || !after.is_empty() {
            audit("simulation", &before, &after, self.conservation_tol)?;
        }
        self.t += dt;
        Ok(report)
    }

    /// One pass over the domains in declared order.
    fn sweep(&mut self, dt: Time, multirate: bool) -> Result<Report, Violation> {
        let now = self.t;
        let mut substeps = Vec::with_capacity(self.domains.len());
        for domain in self.domains.iter_mut() {
            // A quasi-static domain has no state to march, so subdividing its
            // step would just solve the same problem several times.
            let n = if multirate && domain.kind() == Kind::Evolving {
                substeps_for(dt, domain.max_stable_dt(now))
            } else {
                1
            };
            let h = dt / n as f64;
            let mut t = now;
            for _ in 0..n {
                domain.step(t, h, &mut self.bus)?;
                t += h;
            }
            substeps.push((domain.name(), n));
        }
        let residual = self
            .domains
            .iter()
            .map(|d| d.residual())
            .fold(0.0f64, f64::max);
        Ok(Report {
            substeps,
            iterations: 1,
            residual,
        })
    }

    /// Repeat the pass from the same starting state until the residuals settle.
    fn iterate(&mut self, dt: Time, max_iter: u32, tol: f64) -> Result<Report, Violation> {
        if let Some(bad) = self.domains.iter().find(|d| !d.supports_restore()) {
            return Err(Violation::at(
                bad.name(),
                "iterative coupling needs a restorable domain",
                0.0,
            ));
        }
        for domain in self.domains.iter_mut() {
            domain.checkpoint();
        }

        let mut last = Report::default();
        for iteration in 1..=max_iter {
            if iteration > 1 {
                for domain in self.domains.iter_mut() {
                    domain.restore();
                }
                self.bus.clear_offers();
            }
            let mut report = self.sweep(dt, true)?;
            report.iterations = iteration;
            last = report;
            if last.residual <= tol {
                return Ok(last);
            }
        }

        // Not converged. Reporting this rather than proceeding is the whole point:
        // an unconverged coupling produces plausible numbers, which is worse than
        // producing none.
        Err(Violation {
            quantity: "coupling residual".to_string(),
            site: format!("simulation (after {max_iter} iterations)"),
            before: 0.0,
            after: last.residual,
            scale: last.residual.abs(),
            tolerance: tol,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conserved::quantity;

    /// A quasi-static source: converts an input into watts on the bus without any
    /// state of its own. This is the shape optics has — solved, never stepped.
    struct Lamp {
        watts: f64,
        delivered: f64,
    }

    impl Domain for Lamp {
        fn name(&self) -> &'static str {
            "lamp"
        }
        fn kind(&self) -> Kind {
            Kind::QuasiStatic
        }
        fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
            let joules = self.watts * dt.to_si();
            bus.publish(quantity::ENERGY, joules);
            self.delivered += joules;
            Ok(())
        }
        fn ledger(&self) -> Ledger {
            // Energy that has left the lamp is still in the system's books until
            // something else takes it, so the lamp reports what it has paid out.
            Ledger::new().with(quantity::ENERGY, -self.delivered)
        }
        fn checkpoint(&mut self) {}
        fn restore(&mut self) {}
        fn supports_restore(&self) -> bool {
            true
        }
    }

    /// An evolving sink with a stability limit: a lumped thermal mass that must
    /// not be stepped past a fraction of its time constant.
    struct Block {
        joules: f64,
        limit: Time,
        saved: f64,
    }

    impl Domain for Block {
        fn name(&self) -> &'static str {
            "block"
        }
        fn max_stable_dt(&self, _now: Time) -> Time {
            self.limit
        }
        fn step(&mut self, _t: Time, _dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
            self.joules += bus.take(quantity::ENERGY);
            Ok(())
        }
        fn ledger(&self) -> Ledger {
            Ledger::new().with(quantity::ENERGY, self.joules)
        }
        fn checkpoint(&mut self) {
            self.saved = self.joules;
        }
        fn restore(&mut self) {
            self.joules = self.saved;
        }
        fn supports_restore(&self) -> bool {
            true
        }
    }

    fn lamp_and_block(schedule: Schedule, limit: Time) -> Simulation {
        Simulation::new(schedule)
            .with(Lamp {
                watts: 0.01,
                delivered: 0.0,
            })
            .with(Block {
                joules: 0.0,
                limit,
                saved: 0.0,
            })
    }

    /// The chain works end to end: a quasi-static producer hands energy across
    /// the bus to an evolving consumer, the books balance, and the clock moves.
    #[test]
    fn energy_crosses_the_bus_and_the_books_balance() {
        let mut sim = lamp_and_block(Schedule::Staggered, Time::s(1.0));
        let report = sim.advance(Time::s(2.0)).expect("a balanced step");
        assert_eq!(report.iterations, 1);
        assert!((sim.time().to_si() - 2.0).abs() < 1e-15);
        // 10 mW for 2 s is 20 mJ, and all of it arrived.
        assert!((sim.bus().total_consumed(quantity::ENERGY) - 0.02).abs() < 1e-15);
        // The system as a whole is where it started: the lamp is down what the
        // block is up.
        assert_eq!(sim.ledger().get(quantity::ENERGY), Some(0.0));
    }

    /// Energy published and not consumed is caught. This is the interpolation bug
    /// at a coupling interface, in its simplest possible form: a producer with no
    /// consumer.
    #[test]
    fn energy_that_arrives_nowhere_is_a_violation() {
        let mut sim = Simulation::new(Schedule::Staggered).with(Lamp {
            watts: 0.01,
            delivered: 0.0,
        });
        let err = sim.advance(Time::s(1.0)).expect_err("nothing consumed it");
        assert_eq!(err.quantity, "energy");
        assert!(err.site.contains("not consumed"), "{err}");
        // And the clock did not move, so the failure is not half-applied.
        assert_eq!(sim.time(), Time::ZERO);
    }

    /// Multirate: the domain with the tight limit subcycles, and the quasi-static
    /// one does not, because there is nothing to subdivide.
    #[test]
    fn only_evolving_domains_subcycle() {
        let mut sim = lamp_and_block(Schedule::Multirate, Time::s(0.3));
        let report = sim.advance(Time::s(1.0)).unwrap();
        assert_eq!(
            report.substeps,
            vec![("lamp", 1), ("block", 4)],
            "the block needs ceil(1.0/0.3) = 4 substeps; the lamp needs none"
        );
        // Subcycling must not change the total that crossed.
        assert!((sim.bus().total_consumed(quantity::ENERGY) - 0.01).abs() < 1e-15);
    }

    /// A domain with no stability limit is not subcycled at all, however long the
    /// step.
    #[test]
    fn an_unlimited_domain_takes_one_step() {
        let mut sim = lamp_and_block(Schedule::Multirate, Time::from_si(f64::INFINITY));
        let report = sim.advance(Time::s(1e6)).unwrap();
        assert_eq!(report.substeps, vec![("lamp", 1), ("block", 1)]);
    }

    /// Iterative coupling converges and reports how many passes it took.
    struct Settling {
        residual: f64,
        saved: f64,
    }

    impl Domain for Settling {
        fn name(&self) -> &'static str {
            "settling"
        }
        fn step(&mut self, _t: Time, _dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
            // Each pass halves the disagreement with the neighbour.
            self.residual /= 2.0;
            Ok(())
        }
        fn residual(&self) -> f64 {
            self.residual
        }
        fn checkpoint(&mut self) {
            self.saved = self.residual;
        }
        fn restore(&mut self) {
            // The restore puts the state back but keeps the improved coupling
            // guess, which is what makes the iteration converge rather than loop.
            let improved = self.residual;
            self.residual = self.saved.min(improved);
        }
        fn supports_restore(&self) -> bool {
            true
        }
    }

    #[test]
    fn an_iterative_coupling_converges_and_says_how_long_it_took() {
        let mut sim = Simulation::new(Schedule::Iterative {
            max_iter: 20,
            tol: 1e-3,
        })
        .with(Settling {
            residual: 1.0,
            saved: 0.0,
        });
        let report = sim.advance(Time::s(1.0)).unwrap();
        // 1.0 halved ten times is 9.8e-4, the first value under 1e-3.
        assert_eq!(report.iterations, 10);
        assert!(report.residual <= 1e-3);
    }

    /// Not converging is a failure, not a result. An unconverged coupling gives
    /// numbers that look like physics, which is the worst thing it could do.
    #[test]
    fn failing_to_converge_is_reported_not_accepted() {
        let mut sim = Simulation::new(Schedule::Iterative {
            max_iter: 3,
            tol: 1e-9,
        })
        .with(Settling {
            residual: 1.0,
            saved: 0.0,
        });
        let err = sim
            .advance(Time::s(1.0))
            .expect_err("three halvings is not 1e-9");
        assert_eq!(err.quantity, "coupling residual");
        assert!(err.site.contains("after 3 iterations"), "{err}");
        assert_eq!(sim.time(), Time::ZERO);
    }

    /// A domain that cannot put itself back cannot be iterated, and is told so by
    /// name rather than being iterated from the wrong state.
    #[test]
    fn iteration_refuses_a_domain_that_cannot_rewind() {
        struct NoRewind;
        impl Domain for NoRewind {
            fn name(&self) -> &'static str {
                "no-rewind"
            }
            fn step(&mut self, _t: Time, _dt: Time, _b: &mut Exchange) -> Result<(), Violation> {
                Ok(())
            }
        }
        let mut sim = Simulation::new(Schedule::Iterative {
            max_iter: 5,
            tol: 1e-6,
        })
        .with(NoRewind);
        let err = sim.advance(Time::s(1.0)).unwrap_err();
        assert_eq!(err.site, "no-rewind");
        assert!(err.quantity.contains("restorable"), "{err}");
    }

    /// The whole scheduler is deterministic: same domains, same schedule, same
    /// numbers, down to the substep counts.
    #[test]
    fn advancing_is_reproducible() {
        let run = || {
            let mut sim = lamp_and_block(Schedule::Multirate, Time::s(0.07));
            let mut reports = Vec::new();
            for _ in 0..5 {
                reports.push(sim.advance(Time::s(0.25)).unwrap());
            }
            (reports, sim.bus().total_consumed(quantity::ENERGY))
        };
        let (a, ea) = run();
        let (b, eb) = run();
        assert_eq!(a, b);
        assert_eq!(ea.to_bits(), eb.to_bits(), "not bit-identical");
        assert_eq!(a[0].substeps, vec![("lamp", 1), ("block", 4)]);
    }

    /// Taking from a channel empties it, so an amount cannot be consumed twice.
    #[test]
    fn a_channel_cannot_be_drained_twice() {
        let mut bus = Exchange::new();
        bus.publish(quantity::ENERGY, 5.0);
        bus.publish(quantity::ENERGY, 3.0);
        assert_eq!(bus.peek(quantity::ENERGY), 8.0);
        assert_eq!(bus.take(quantity::ENERGY), 8.0);
        assert_eq!(bus.take(quantity::ENERGY), 0.0);
        assert_eq!(bus.total_consumed(quantity::ENERGY), 8.0);
        assert!(bus.unclaimed().next().is_none());
    }
}
