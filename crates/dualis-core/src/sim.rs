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

use std::any::Any;
use std::collections::BTreeMap;

use dualis_units::Time;

use crate::bodies::Bodies;
use crate::conserved::{audit_with, Ledger, Tolerances, Violation};
use crate::field::ScalarField;
use crate::integrator::substeps_for;
use crate::scene::{mismatch, Flux, Interface};

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
    /// What this domain is called. Used to look it up and to name it in a violation.
    ///
    /// Borrowed rather than `&'static str`, so a name can come from a scene file. That was
    /// the first thing the workspace's own application could not do: every constructor
    /// wanted a compile-time name and the name it had was a `String` read off disk, so it
    /// leaked one per domain to get past the signature.
    fn name(&self) -> &str;

    /// Whether it has state to roll forward. Defaults to [`Kind::Evolving`].
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

    /// Whether [`Domain::checkpoint`] and [`Domain::restore`] actually do something.
    ///
    /// [`Schedule::Iterative`] refuses to run a domain that says no, rather than iterating
    /// from the wrong state and reporting a residual that means nothing.
    fn supports_restore(&self) -> bool {
        false
    }

    /// This domain as [`Any`], so a caller can get the concrete type back out of a
    /// [`Simulation`] — see [`Simulation::domain_as`].
    ///
    /// Opt-in, and returning `None` by default, because it cannot be automatic. Deriving it
    /// from the trait would need `Domain: Any` plus upcasting `dyn Domain` to `dyn Any`,
    /// which is a newer Rust than this crate promises. A domain that wants to be inspected
    /// writes `fn as_any(&self) -> Option<&dyn Any> { Some(self) }` and is done.
    ///
    /// The coupling never needs this: domains meet through [`Exchange`] and nothing else,
    /// which is the property the whole design rests on. What needs it is everything *around*
    /// the simulation — a test asserting a temperature profile, a visualiser drawing one —
    /// and that is a reader, not a participant.
    fn as_any(&self) -> Option<&dyn Any> {
        None
    }

    /// The same, mutably, so a caller can *write* to a domain between steps.
    ///
    /// **This does not weaken "domains never read each other."** That rule is about what happens
    /// inside [`Domain::step`], where the only channel is [`Exchange`]. This is the owner of the
    /// simulation, outside the step loop, holding `&mut Simulation` already — it could drop the
    /// domain and rebuild it, so denying it a write was never protecting anything.
    ///
    /// What needs it is a feedback loop the bus cannot carry. A copper winding's resistance rises
    /// with its temperature, and that temperature lives in a thermal domain: neither can see the
    /// other's state, and neither should. A caller between frames can see both, and until this
    /// existed it could read one and not write the other, which made the loop unclosable from
    /// anywhere at all.
    ///
    /// Opt-in and `None` by default, like [`Domain::as_any`] — and that default is a hazard this
    /// workspace has been bitten by twice, in `FRICTION.md` findings 7 and 12: a domain that
    /// forgets it is not broken, it is silently absent from whatever asks. If you implement
    /// `as_any`, implement this beside it.
    fn as_any_mut(&mut self) -> Option<&mut dyn Any> {
        None
    }

    /// This domain as a [`ScalarField`], if it has one to show.
    ///
    /// Opt-in and `None` by default, in the same style as [`Domain::as_any`] and for a
    /// sharper reason than that one. `ScalarField` was written as the interface a visualiser
    /// would read a simulation through, and then a visualiser found it unreachable: it holds
    /// `&dyn Domain`, and there was no way to ask that for a field. So it downcast to
    /// concrete types instead and knew every domain by name — precisely what the interface
    /// existed to avoid.
    ///
    /// A domain with a field writes `fn as_field(&self) -> Option<&dyn ScalarField>
    /// { Some(self) }`. See [`Simulation::field`].
    fn as_field(&self) -> Option<&dyn ScalarField> {
        None
    }

    /// The named scalars this domain reports, for a table, a chart or a caption.
    ///
    /// **The number a domain has when it has no picture.** A source has a remaining tank, a
    /// winding has a dissipation, a thermal network has a temperature per node — and for several
    /// of those the scalar *is* the result. `as_field` covers the domains that are continua and
    /// there was no counterpart for the rest, so a caller wanting them had to know every domain
    /// by name and downcast to each.
    ///
    /// That is what makes this a trait method rather than a function somewhere above: a layer
    /// that collects readings by matching on domain types has to be edited every time a physics
    /// is added, which is the one thing this workspace's structure exists to avoid.
    ///
    /// Return what the domain is *for* rather than a uniform summary. A mean over a pressure
    /// field is zero by symmetry and would be a column of noise; the peak is the number a reader
    /// wants. Nobody but the domain knows which.
    ///
    /// Empty by default, and opt-in like [`as_any`](Domain::as_any) and
    /// [`as_field`](Domain::as_field) — with the hazard `as_any` has already taught once: four
    /// mechanics domains never opted into it, and an orbit scene ran, conserved, and drew nothing
    /// at all. A domain that forgets this one is silently absent from every table, not broken.
    fn readings(&self) -> Vec<Reading> {
        Vec::new()
    }

    /// This domain as a countable set of bodies, if that is what it is.
    ///
    /// The counterpart to [`as_field`](Domain::as_field), and between them they cover both kinds
    /// of thing a domain can be. A caller wanting to draw, measure or export no longer has to
    /// name `NBody`, `ContactSystem` or `Fluid` — which it did for months, recorded as
    /// `FRICTION.md` finding 11, until splitting the layers made it unpayable.
    ///
    /// Opt-in and `None` by default, with the hazard that default has now taught three times: a
    /// domain that forgets is silently absent rather than broken.
    fn as_bodies(&self) -> Option<&dyn Bodies> {
        None
    }
}

/// Delegation, so a domain chosen at run time can be added like any other.
///
/// Without this a caller holding `Box<dyn Domain>` — which is what building from data
/// produces — could not hand it to [`Simulation::with`], even though the simulation stores
/// exactly that internally. Prefer [`Simulation::with_boxed`], which avoids boxing the box;
/// this impl is here so that generic code over `impl Domain` works on a boxed one too.
impl Domain for Box<dyn Domain> {
    fn name(&self) -> &str {
        (**self).name()
    }
    fn kind(&self) -> Kind {
        (**self).kind()
    }
    fn max_stable_dt(&self, now: Time) -> Time {
        (**self).max_stable_dt(now)
    }
    fn step(&mut self, t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        (**self).step(t, dt, bus)
    }
    fn residual(&self) -> f64 {
        (**self).residual()
    }
    fn ledger(&self) -> Ledger {
        (**self).ledger()
    }
    fn checkpoint(&mut self) {
        (**self).checkpoint()
    }
    fn restore(&mut self) {
        (**self).restore()
    }
    fn supports_restore(&self) -> bool {
        (**self).supports_restore()
    }
    fn as_any_mut(&mut self) -> Option<&mut dyn Any> {
        (**self).as_any_mut()
    }
    fn readings(&self) -> Vec<Reading> {
        (**self).readings()
    }
    fn as_bodies(&self) -> Option<&dyn Bodies> {
        (**self).as_bodies()
    }
    fn as_any(&self) -> Option<&dyn Any> {
        (**self).as_any()
    }
    fn as_field(&self) -> Option<&dyn ScalarField> {
        (**self).as_field()
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
    /// Channels that carry a place as well as an amount, keyed by
    /// `(interface name, channel)` so the audit reports them in a fixed order.
    spatial: BTreeMap<(String, &'static str), Flux>,
    spatial_consumed: BTreeMap<(String, &'static str), f64>,
    /// The outer step the current sweep is covering, in seconds. Zero when nobody has said —
    /// a bare `Exchange` in a test — and [`Exchange::take_share`] falls back to taking
    /// everything, which is the honest answer when the interval is unknown.
    interval: f64,
    /// How much of `interval` is still unclaimed, per channel. See `take_share`.
    unclaimed_time: BTreeMap<&'static str, f64>,
    /// How many separate `take` calls each channel saw this step.
    ///
    /// Counted because the conservation audit structurally cannot see the failure it detects.
    /// [`Exchange::take`] empties a channel, so a *second* consumer of the same channel gets
    /// zero — and the books balance perfectly, because everything published was taken. Two
    /// plates under one lamp warm at the rate of one plate, and the audit reports it clean.
    ///
    /// Every scene and every integration test in this workspace had at most one consumer per
    /// channel, which is why this went unnoticed until a world with six domains was attempted.
    takers: BTreeMap<&'static str, u32>,
}

impl Exchange {
    /// An empty bus.
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
        *self.takers.entry(channel).or_insert(0) += 1;
        amount
    }

    /// Look without taking.
    pub fn peek(&self, channel: &'static str) -> f64 {
        self.published.get(channel).copied().unwrap_or(0.0)
    }

    /// Take the share of a channel that belongs to a substep of length `dt`.
    ///
    /// For a domain that subcycles. [`Exchange::take`] empties the channel, which is right for
    /// a domain stepping once per interval and wrong for one stepping many times: a publisher
    /// offers a whole outer step's worth at once, so the first substep would take all of it and
    /// the rest would find the channel dark. Every joule of the interval then lands at its
    /// beginning, and **refining the substep stops improving the answer** — see
    /// [`Schedule::Multirate`], where the measured error is 26% at a 300 s outer step whatever
    /// the substep count.
    ///
    /// The share is taken against the time *remaining*, not against the whole interval. That is
    /// what makes it exact: after handing out `A·dt/T` and reducing both, `A/T` is unchanged, so
    /// the last substep — which asks for at least what is left — receives the remainder and the
    /// channel ends empty to the last bit. Apportioning against the whole interval instead
    /// leaves `O(n·ε·A)` stranded, and [`Exchange::audit_transfers`] uses an absolute tolerance
    /// that would eventually refuse it.
    ///
    /// Falls back to [`Exchange::take`] when the interval is unknown, so a domain written
    /// against this works unchanged under a bare `Exchange` and under
    /// [`Schedule::Staggered`], where it steps once and the share is the whole.
    pub fn take_share(&mut self, channel: &'static str, dt: Time) -> f64 {
        let h = dt.to_si();
        if self.interval <= 0.0 || !h.is_finite() || h <= 0.0 {
            return self.take(channel);
        }
        let left = *self.unclaimed_time.entry(channel).or_insert(self.interval);
        // The last substep asks for everything that is left, and gets it. Compared with a
        // slack of `1e-12` of the interval rather than exactly, because `n` substeps of `dt/n`
        // do not sum to `dt` in binary: three of a third leave a residue one ulp wide, and an
        // exact comparison misses the final share and strands it on the channel.
        if h >= left || left - h <= self.interval * 1e-12 {
            self.unclaimed_time.insert(channel, 0.0);
            return self.take(channel);
        }
        let amount = self.published.get(channel).copied().unwrap_or(0.0);
        let share = amount * h / left;
        self.unclaimed_time.insert(channel, left - h);
        *self.published.entry(channel).or_insert(0.0) -= share;
        *self.consumed.entry(channel).or_insert(0.0) += share;
        share
    }

    /// Tell the bus what interval the current sweep covers, so [`Exchange::take_share`] can
    /// apportion. Called by [`Simulation::advance`]; a standalone `Exchange` need not.
    pub fn covering(&mut self, dt: Time) {
        self.interval = dt.to_si().max(0.0);
        self.unclaimed_time.clear();
        self.takers.clear();
    }

    /// Offer an amount that knows where on a boundary it landed.
    ///
    /// The spatial counterpart of [`publish`](Exchange::publish), and the reason
    /// [`scene`](crate::scene) exists: a coating absorbs where the beam is, and a lumped
    /// number cannot say that. Repeated publishes accumulate face by face, so two
    /// mechanisms heating the same surface add up in place.
    ///
    /// Refuses a [`Flux`] whose face count does not match the interface. Silently padding
    /// or truncating would put energy on the wrong part of the boundary, which is worse
    /// than losing it — losing it the audit would catch.
    pub fn publish_on(
        &mut self,
        interface: &Interface,
        channel: &'static str,
        flux: &Flux,
    ) -> Result<(), Violation> {
        if flux.faces() != interface.faces() {
            return Err(mismatch(
                &format!("publish on {}/{channel}", interface.name()),
                interface.faces(),
                flux.faces(),
            ));
        }
        let key = (interface.name().to_string(), channel);
        match self.spatial.get_mut(&key) {
            Some(existing) => existing.add(flux),
            None => {
                self.spatial.insert(key, flux.clone());
                Ok(())
            }
        }
    }

    /// Take everything offered on an interface's channel, leaving it empty.
    ///
    /// Returns zeros rather than an error when nothing was published, because a consumer
    /// stepping a boundary that happens to be dark this step is not a fault. A face-count
    /// disagreement *is*, and is reported: the two sides do not share a discretisation, and
    /// the fix is [`Flux::resample`] at whichever side owns the decision.
    pub fn take_on(
        &mut self,
        interface: &Interface,
        channel: &'static str,
    ) -> Result<Flux, Violation> {
        let key = (interface.name().to_string(), channel);
        // Removed rather than zeroed. A drained channel is empty, and an empty channel
        // should not go on pinning a face count for the rest of the step — the next
        // publisher on that boundary is entitled to its own discretisation.
        let Some(offered) = self.spatial.remove(&key) else {
            return Ok(Flux::zeros(interface.faces()));
        };
        if offered.faces() != interface.faces() {
            // Put it back: a consumer that could not read it has not consumed it, and the
            // audit should still see the energy sitting there unclaimed.
            let found = offered.faces();
            self.spatial.insert(key, offered);
            return Err(mismatch(
                &format!("take from {}/{channel}", interface.name()),
                interface.faces(),
                found,
            ));
        }
        *self.spatial_consumed.entry(key).or_insert(0.0) += offered.total();
        Ok(offered)
    }

    /// Look at a spatial channel without taking it.
    pub fn peek_on(&self, interface: &Interface, channel: &'static str) -> Option<&Flux> {
        self.spatial.get(&(interface.name().to_string(), channel))
    }

    /// Channels that were published to but never taken from, with what is left on
    /// them. Energy sitting here at the end of a step is energy that left one
    /// domain and arrived nowhere.
    ///
    /// Spatial channels appear as `"interface/channel"`, with the total left on them.
    pub fn unclaimed(&self) -> impl Iterator<Item = (String, f64)> + '_ {
        self.published
            .iter()
            .filter(|(_, v)| v.abs() > 0.0)
            .map(|(k, v)| ((*k).to_string(), *v))
            .chain(
                self.spatial
                    .iter()
                    .filter(|(_, f)| f.total().abs() > 0.0)
                    .map(|((i, c), f)| (format!("{i}/{c}"), f.total())),
            )
    }

    /// Fail if anything published was not consumed.
    ///
    /// This is the check that catches a coupling whose two sides disagree — a
    /// surface that absorbed 3.7 mW handing it to a mesh that received 3.4 mW
    /// because the interpolation between their discretisations lost the rest.
    ///
    /// The original design said that, and then could not check it: with one number per
    /// channel there was no discretisation to disagree about. Spatial channels close that
    /// gap, and they are audited **face by face** rather than on their total — a
    /// redistribution that moves heat from one side of a mirror to the other keeps the sum
    /// exactly right, so a total-only check would pass the one bug the spatial coupling
    /// exists to prevent. The failure names the face.
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
        for ((interface, channel), flux) in self.spatial.iter() {
            for (face, left) in flux.per_face().iter().enumerate() {
                if left.abs() > abs_tol {
                    return Err(Violation {
                        quantity: format!("{interface}/{channel} face {face}"),
                        site: format!("{site} (published but not consumed)"),
                        before: *left,
                        after: 0.0,
                        scale: left.abs(),
                        tolerance: abs_tol,
                    });
                }
            }
        }
        Ok(())
    }

    /// Total taken from a channel over the run, for reporting.
    pub fn total_consumed(&self, channel: &str) -> f64 {
        self.consumed.get(channel).copied().unwrap_or(0.0)
    }

    /// Total taken from a spatial channel over the run, summed over its faces.
    pub fn total_consumed_on(&self, interface: &Interface, channel: &'static str) -> f64 {
        self.spatial_consumed
            .get(&(interface.name().to_string(), channel))
            .copied()
            .unwrap_or(0.0)
    }

    /// Empty the offers, keeping the running consumption totals.
    pub fn clear_offers(&mut self) {
        self.published.clear();
        self.spatial.clear();
        self.unclaimed_time.clear();
        self.takers.clear();
    }

    /// How many times each channel has been taken from this sweep.
    ///
    /// Raw counts, because the bus cannot interpret them: a domain subcycling ten times takes
    /// ten times, and ten domains taking once each also takes ten times. Only
    /// [`Simulation`] knows whose turn it was, and it compares this between turns — see
    /// `Simulation::sweep`, where the check that a channel had at most one *consumer* lives.
    pub fn takes_per_channel(&self) -> impl Iterator<Item = (&'static str, u32)> + '_ {
        self.takers.iter().map(|(c, n)| (*c, *n))
    }
}

/// One named scalar from one domain at one instant.
///
/// Deliberately flat and owned: it crosses a layer boundary, gets written to a CSV column and a
/// chart legend, and neither of those wants a borrow into a running simulation.
#[derive(Clone, Debug, PartialEq)]
pub struct Reading {
    /// Which domain it came from. Filled in by the domain, because only it knows its own name.
    pub domain: String,
    /// What it is — `"mean"`, `"peak"`, `"reserve"`, a node's name.
    pub label: String,
    /// The value, in SI, with one exception this workspace has already made everywhere else:
    /// temperatures are celsius, because that is the unit a column of them is read in.
    pub value: f64,
    /// The unit, for a header row or an axis. `&'static str` because a unit is a compile-time
    /// fact about the quantity, not data — unlike a domain's name, which comes from a file.
    pub unit: &'static str,
}

impl Reading {
    /// A reading, named.
    pub fn new(
        domain: impl Into<String>,
        label: impl Into<String>,
        value: f64,
        unit: &'static str,
    ) -> Reading {
        Reading {
            domain: domain.into(),
            label: label.into(),
            value,
            unit,
        }
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
    Iterative {
        /// Give up after this many sweeps. Reaching it is a [`Violation`], not a result.
        max_iter: u32,
        /// The residual every domain must fall under for the step to be accepted.
        tol: f64,
    },
    /// As [`Schedule::Staggered`], but each evolving domain takes as many equal
    /// substeps as its own stability limit needs.
    ///
    /// # It does not refine a coupled quantity, and the audit cannot tell you
    ///
    /// Read this before choosing it for accuracy, because that is the obvious reason to and it
    /// is the wrong one.
    ///
    /// One domain is stepped to completion before the next. A quasi-static publisher is never
    /// subcycled, so it puts a whole outer step's worth on the bus once; a subcycling consumer
    /// then calls [`Exchange::take`] on its **first** substep and takes all of it. So every
    /// joule of the interval is deposited at its beginning and decays for the rest of it, and
    /// refining the substep does not move the answer toward the truth. Taking the limit of
    /// `u ← u·gⁿ + (P·dt/C)·g^(n−1)` with `g = 1 − h/τ` as `n → ∞` gives
    /// `u·e^(−dt/τ) + (P·dt/C)·e^(−dt/τ)`, which is not the solution: the error is first order
    /// in the **outer** step and independent of the substep entirely.
    ///
    /// Measured on a lumped plate under a steady lamp, against the closed form: 26.2% low at a
    /// 300 s outer step, 13.8% at 150 s, 7.1% at 75 s — *whatever* the substep count. At the
    /// same outer step it is not reliably better than [`Schedule::Staggered`] and at a coarse
    /// one it is worse, with the errors on opposite sides.
    ///
    /// **Every one of those runs passes the conservation audit at around 1e-12.** The total
    /// that crossed is exactly right; only its distribution in time is wrong, and a [`Ledger`]
    /// has no representation for *when*. This is the time-domain twin of the reason
    /// [`Exchange::audit_transfers`] had to become a per-face check in space — a quantity moved
    /// to the wrong part of an interval keeps its total, and conservation is blind to it.
    ///
    /// So: choose this for **stability**, which is what it delivers — a domain whose limit is a
    /// hundredth of the frame no longer forces the frame to shrink. Choose the outer step for
    /// **accuracy**, because that is what sets it. `crates/dualis/tests/multirate_timing.rs`
    /// pins the consequence.
    Multirate,
}

/// What one [`Simulation::advance`] actually did.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Report {
    /// Substeps taken, per domain, in declared order.
    ///
    /// Owned names, because [`Domain::name`] is borrowed from the domain and this report
    /// outlives the borrow — the same consequence of names being data rather than
    /// constants that shows up everywhere else in this module.
    pub substeps: Vec<(String, u32)>,
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
    conservation_tol: Tolerances,
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
            conservation_tol: Tolerances::default(),
        }
    }

    /// Add a domain whose type was chosen at run time.
    ///
    /// What [`Simulation::with`] cannot do: building a domain from a scene file produces a
    /// `Box<dyn Domain>`, and `with` wants a concrete type. The simulation has always stored
    /// boxes internally, so this is the shorter path and not a wider one.
    pub fn with_boxed(mut self, domain: Box<dyn Domain>) -> Simulation {
        self.domains.push(domain);
        self
    }

    /// Add a domain. Order matters for [`Schedule::Staggered`] and its relatives: a domain
    /// sees the output of those declared before it from this step, and of those after it from
    /// the last one.
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
    /// step, for every quantity that has no override. Default 1e-9.
    pub fn conservation_tolerance(mut self, tol: f64) -> Simulation {
        let overrides: Vec<(&'static str, f64)> = self.conservation_tol.overrides().collect();
        self.conservation_tol = overrides
            .into_iter()
            .fold(Tolerances::uniform(tol), |t, (q, v)| t.with(q, v));
        self
    }

    /// Relative tolerance for **one** quantity, overriding the default.
    ///
    /// The reason this exists: a Barnes-Hut N-body gives up exact momentum by construction, and
    /// energy in a rigid room is exact to `1e-15`. Under one number either the momentum check
    /// refuses a correct run or the energy check stops being able to see anything. A quantity's
    /// achievable accuracy is a property of the scheme carrying it.
    ///
    /// ```
    /// # use dualis_core::{Schedule, Simulation};
    /// # use dualis_core::conserved::quantity;
    /// let sim = Simulation::new(Schedule::Staggered)
    ///     .conservation_tolerance(1e-12)
    ///     .conservation_tolerance_for(quantity::MOMENTUM, 1e-6);
    /// assert_eq!(sim.tolerances().for_quantity(quantity::ENERGY), 1e-12);
    /// assert_eq!(sim.tolerances().for_quantity(quantity::MOMENTUM), 1e-6);
    /// ```
    pub fn conservation_tolerance_for(mut self, quantity: &'static str, tol: f64) -> Simulation {
        self.conservation_tol = std::mem::take(&mut self.conservation_tol).with(quantity, tol);
        self
    }

    /// What this simulation checks each quantity against.
    pub fn tolerances(&self) -> &Tolerances {
        &self.conservation_tol
    }

    /// How far the simulation has been advanced.
    pub fn time(&self) -> Time {
        self.t
    }

    /// The coupling bus, for reading what crossed between domains.
    pub fn bus(&self) -> &Exchange {
        &self.bus
    }

    /// Every domain, in the order they were added.
    ///
    /// `domain` answers by name, which is right for a caller that knows what it is looking for
    /// and useless for one that must visit them all. A layer capturing a run has to enumerate,
    /// and without this it had to be handed the list by whoever built the simulation — which
    /// means the layer above knows the composition rather than asking.
    ///
    /// Order is declaration order, which is also execution order under the staggered schedules,
    /// so a caller iterating this sees domains in the order they act.
    pub fn domains(&self) -> impl Iterator<Item = &dyn Domain> + '_ {
        self.domains.iter().map(|d| &**d as &dyn Domain)
    }

    /// A domain by name, through the trait. For the concrete type, see
    /// [`Simulation::domain_as`].
    pub fn domain(&self, name: &str) -> Option<&dyn Domain> {
        self.domains
            .iter()
            .find(|d| d.name() == name)
            .map(|d| d.as_ref())
    }

    /// A domain's [`ScalarField`], if it has one and opted in.
    ///
    /// The domain-agnostic counterpart of [`Simulation::domain_as`]: a renderer can sample
    /// every field in a simulation without knowing what any of them are. That was the whole
    /// point of `ScalarField` and it was not reachable until [`Domain::as_field`] existed.
    pub fn field(&self, name: &str) -> Option<&dyn ScalarField> {
        self.domain(name)?.as_field()
    }

    /// A domain by name and concrete type, for a caller that needs more than the
    /// [`Domain`] trait exposes — a temperature profile, a body's position.
    ///
    /// Returns `None` if the name is not here, if the type is wrong, or if that domain did
    /// not implement [`Domain::as_any`]. Prefer [`Simulation::field`] when what is wanted is
    /// a field to sample: that one does not need the concrete type at all.
    pub fn domain_as<T: Any>(&self, name: &str) -> Option<&T> {
        self.domain(name)?.as_any()?.downcast_ref::<T>()
    }

    /// The same, mutably, for a caller closing a feedback loop between steps.
    ///
    /// `None` if there is no such domain, if it is not a `T`, or if it does not implement
    /// [`Domain::as_any_mut`] — three different reasons that look alike from here, which is why
    /// that method's documentation asks for it to be implemented beside `as_any`.
    pub fn domain_as_mut<T: Any>(&mut self, name: &str) -> Option<&mut T> {
        self.domains
            .iter_mut()
            .find(|d| d.name() == name)?
            .as_any_mut()?
            .downcast_mut::<T>()
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
        // What a substep's share is measured against. Set here rather than in `sweep`, because
        // `iterate` sweeps repeatedly over the same interval.
        self.bus.covering(dt);
        let report = match self.schedule {
            Schedule::OneWay | Schedule::Staggered => self.sweep(dt, false)?,
            Schedule::Multirate => self.sweep(dt, true)?,
            Schedule::Iterative { max_iter, tol } => self.iterate(dt, max_iter, tol)?,
        };

        self.bus.audit_transfers("bus", self.transfer_tol)?;
        let after = self.ledger();
        if !before.is_empty() || !after.is_empty() {
            audit_with("simulation", &before, &after, &self.conservation_tol)?;
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
            // Which channels had already been drawn on before this domain's turn.
            let before: Vec<(&'static str, u32)> = self.bus.takes_per_channel().collect();
            for _ in 0..n {
                domain.step(t, h, &mut self.bus)?;
                t += h;
            }

            // A channel this domain took from that an *earlier* domain had already emptied.
            //
            // `Exchange::take` empties a channel, so the second consumer gets zero — and every
            // total agrees, because everything published was consumed. Two plates under one lamp
            // warm at the rate of one plate and the books balance to the bit. The conservation
            // audit structurally cannot see it.
            //
            // Counted per *turn* rather than per call, because a subcycling domain takes once
            // per substep and that is one consumer collecting its own interval in pieces.
            //
            // Refused rather than apportioned: splitting needs a rule the kernel has no way to
            // choose — equally, by heat capacity, by area? — and any rule it picked would be
            // silently wrong for someone, which is the failure being fixed rather than a fresh
            // one. A caller who knows the answer can publish on channels of their own.
            for (channel, now_taken) in self.bus.takes_per_channel() {
                let was = before
                    .iter()
                    .find(|(c, _)| *c == channel)
                    .map_or(0, |(_, n)| *n);
                if was > 0 && now_taken > was {
                    return Err(Violation {
                        quantity: channel.to_string(),
                        site: format!(
                            "{} (a second domain took from a channel already emptied)",
                            domain.name()
                        ),
                        before: was as f64,
                        after: now_taken as f64,
                        scale: now_taken as f64,
                        tolerance: 0.0,
                    });
                }
            }
            substeps.push((domain.name().to_string(), n));
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
    use dualis_units::Area;

    /// A quasi-static source: converts an input into watts on the bus without any
    /// state of its own. This is the shape optics has — solved, never stepped.
    struct Lamp {
        watts: f64,
        delivered: f64,
    }

    impl Domain for Lamp {
        fn name(&self) -> &str {
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
        fn name(&self) -> &str {
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
            vec![("lamp".to_string(), 1), ("block".to_string(), 4)],
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
        assert_eq!(
            report.substeps,
            vec![("lamp".to_string(), 1), ("block".to_string(), 1)]
        );
    }

    /// Iterative coupling converges and reports how many passes it took.
    struct Settling {
        residual: f64,
        saved: f64,
    }

    impl Domain for Settling {
        fn name(&self) -> &str {
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
            fn name(&self) -> &str {
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
        assert_eq!(
            a[0].substeps,
            vec![("lamp".to_string(), 1), ("block".to_string(), 4)]
        );
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

    /// A spatial channel behaves like a lumped one — accumulate, drain once — but face by
    /// face, so two mechanisms heating the same mirror add up *where* each of them did.
    #[test]
    fn a_spatial_channel_accumulates_and_drains_in_place() {
        let mirror = Interface::uniform("mirror", 4, Area::from_si(1e-4));
        let mut bus = Exchange::new();

        // Absorption in the coating, on the two faces the beam covers.
        bus.publish_on(
            &mirror,
            quantity::ENERGY,
            &Flux::from_faces(vec![0.0, 2.0, 3.0, 0.0]),
        )
        .unwrap();
        // And a mount conducting into one edge, which is a different mechanism on the same
        // boundary. It must land on face 0, not be averaged in.
        bus.publish_on(
            &mirror,
            quantity::ENERGY,
            &Flux::from_faces(vec![1.0, 0.0, 0.0, 0.0]),
        )
        .unwrap();

        assert_eq!(
            bus.peek_on(&mirror, quantity::ENERGY).unwrap().per_face(),
            &[1.0, 2.0, 3.0, 0.0]
        );

        let taken = bus.take_on(&mirror, quantity::ENERGY).unwrap();
        assert_eq!(taken.per_face(), &[1.0, 2.0, 3.0, 0.0]);
        assert!((bus.total_consumed_on(&mirror, quantity::ENERGY) - 6.0).abs() < 1e-15);
        // Emptied, so it cannot be consumed twice.
        assert_eq!(bus.take_on(&mirror, quantity::ENERGY).unwrap().total(), 0.0);
        assert!(bus.unclaimed().next().is_none());

        // A channel nobody published to reads as zeros over the right boundary, not an
        // error: a mirror that happens to be dark this step is not a fault.
        let dark = bus.take_on(&mirror, "photons").unwrap();
        assert_eq!(dark.faces(), 4);
        assert_eq!(dark.total(), 0.0);
    }

    /// **The bug the spatial audit exists to catch.** A consumer that keeps the total but
    /// moves it to the wrong part of the boundary is invisible to a total-only check, and
    /// is exactly the failure a shared discretisation is supposed to prevent.
    #[test]
    fn the_audit_names_the_face_that_was_left_holding_something() {
        let mirror = Interface::uniform("mirror", 8, Area::from_si(1e-4));
        let mut bus = Exchange::new();

        // Ten joules on face 6.
        let mut absorbed = vec![0.0; 8];
        absorbed[6] = 10.0;
        bus.publish_on(&mirror, quantity::ENERGY, &Flux::from_faces(absorbed))
            .unwrap();

        // A consumer takes it and puts back the same total in the wrong place. The sum is
        // exactly right, and the sum is not what is being checked.
        let taken = bus.take_on(&mirror, quantity::ENERGY).unwrap();
        let mut misplaced = vec![0.0; 8];
        misplaced[1] = -taken.total();
        misplaced[2] = taken.total();
        bus.publish_on(&mirror, quantity::ENERGY, &Flux::from_faces(misplaced))
            .unwrap();

        assert!(
            bus.peek_on(&mirror, quantity::ENERGY)
                .unwrap()
                .total()
                .abs()
                < 1e-12,
            "the total balances, which is the whole point of the example"
        );
        let err = bus
            .audit_transfers("mirror coupling", 1e-9)
            .expect_err("a redistribution that keeps the total must still be caught");
        assert!(err.quantity.contains("face 1"), "{err}");
        assert!(err.quantity.contains("mirror/energy"), "{err}");
    }

    /// Two sides that do not share a discretisation are refused rather than resampled
    /// behind the caller's back, on both the publishing and the consuming side.
    #[test]
    fn a_discretisation_disagreement_is_refused_at_the_bus() {
        let coarse = Interface::uniform("mirror", 4, Area::from_si(1e-4));
        let fine = Interface::uniform("mirror", 16, Area::from_si(0.25e-4));
        let mut bus = Exchange::new();

        // Publishing 16 faces onto a 4-face boundary.
        let err = bus
            .publish_on(&coarse, quantity::ENERGY, &Flux::zeros(16))
            .expect_err("16 faces is not 4 faces");
        assert!(err.quantity.contains("expected 4"), "{err}");
        assert!(err.site.contains("mirror/energy"), "{err}");

        // And a consumer whose own boundary is finer than what was published. Note both
        // interfaces are named "mirror": the channel matches, the discretisation does not,
        // and it is the face count that decides.
        bus.publish_on(&coarse, quantity::ENERGY, &Flux::from_faces(vec![1.0; 4]))
            .unwrap();
        let err = bus
            .take_on(&fine, quantity::ENERGY)
            .expect_err("a 16-cell mesh must not read a 4-face flux");
        assert!(err.quantity.contains("expected 16"), "{err}");
        assert!(err.quantity.contains("found 4"), "{err}");

        // A refused take consumed nothing, so the energy is still there to be found.
        assert!((bus.peek_on(&coarse, quantity::ENERGY).unwrap().total() - 4.0).abs() < 1e-15);
        assert_eq!(bus.total_consumed_on(&coarse, quantity::ENERGY), 0.0);
        assert!(bus.audit_transfers("mirror", 1e-9).is_err());

        // Saying it explicitly is what works, and it conserves.
        let crossed = bus
            .take_on(&coarse, quantity::ENERGY)
            .unwrap()
            .resample(&coarse, &fine)
            .unwrap();
        assert_eq!(crossed.faces(), 16);
        assert!((crossed.total() - 4.0).abs() < 1e-12);
    }
}
