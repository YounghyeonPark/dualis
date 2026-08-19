//! What `Schedule::Multirate` does with a coupled quantity, and what the audit cannot see.
//!
//! Found from outside, by an application that compared a lumped plate under a lamp against the
//! closed form of its own scheme. The total energy is right to 1e-12 every step; the *time* at
//! which it arrives is not, and `Ledger` has no representation for when.

use pantometry::prelude::*;

/// A steady source with a finite tank, quasi-static: it publishes once per outer step.
struct Lamp {
    watts: f64,
    reserve: f64,
    saved: Option<f64>,
}

impl Domain for Lamp {
    fn name(&self) -> &str {
        "lamp"
    }
    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }
    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let j = (self.watts * dt.to_si()).min(self.reserve);
        self.reserve -= j;
        bus.publish(HEAT, j);
        Ok(())
    }
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.reserve)
    }
    fn checkpoint(&mut self) {
        self.saved = Some(self.reserve);
    }
    fn restore(&mut self) {
        if let Some(r) = self.saved {
            self.reserve = r;
        }
    }
    fn supports_restore(&self) -> bool {
        true
    }
}

/// Aluminium with the radiation turned off, so the loss is exactly `hA·(T − T_a)` and the
/// closed form below is the closed form of the problem rather than of a linearisation.
///
/// This is the difference between a reference and a plausible number. With
/// `aluminium_6061`'s real emissivity the `T⁴` term makes the loss larger than `hA·ΔT`, the
/// equilibrium lower, and the discrete solution sits *below* a linear closed form — which is
/// the wrong side for explicit Euler, since `(1 − h/τ)^n < e^(−t/τ)` means it must overshoot.
/// The first version of this test compared against the linear form anyway, and read the
/// mismatch as a convergence failure. Refining the step then made the disagreement *worse*,
/// which is the signal that the reference was wrong and not the scheme.
///
/// Buildable from the prelude only because `ThermalProps` and the two unit types it needs were
/// added to it — which was itself a finding, from a consumer who wanted exactly this.
fn grey_aluminium() -> Substance {
    let mut s = Substance::aluminium_6061();
    if let Some(t) = s.thermal.as_mut() {
        t.emissivity = 0.0;
    }
    s
}

fn plate() -> LumpedMass {
    LumpedMass::new(
        "plate",
        grey_aluminium(),
        Volume::from_si(60e-3 * 60e-3 * 3e-3),
        Length::mm(1.5),
        Temperature::celsius(20.0),
        Environment::still_air(
            Temperature::celsius(20.0),
            Area::from_si(2.0 * 60e-3 * 60e-3),
        ),
    )
}

/// Run to 600 s in `outer` steps and report the rise.
fn rise_after(schedule: Schedule, outer: Time, steps: usize) -> f64 {
    let mut sim = Simulation::new(schedule)
        .conservation_tolerance(1e-6)
        .with(Lamp {
            watts: 2.0,
            reserve: f64::INFINITY,
            saved: None,
        })
        .with(plate());
    for _ in 0..steps {
        sim.advance(outer).expect("the books close");
    }
    sim.domain_as::<LumpedMass>("plate")
        .expect("the plate is there")
        .rise()
        .to_si()
}

/// **Subcycling improves a coupled quantity, which it did not until `take_share` existed.**
///
/// `Simulation::sweep` steps one domain to completion before the next, so a quasi-static
/// publisher offers a whole outer step's joules at once. A subcycling consumer calling
/// [`Exchange::take`] would empty the channel on its *first* substep, depositing every joule of
/// the interval at its beginning — and then refining the substep stops improving anything. The
/// limit of `u ← u·gⁿ + (P·dt/C)·g^(n−1)` as `n → ∞` is `u·e^(−dt/τ) + (P·dt/C)·e^(−dt/τ)`,
/// which is not the solution: first order in the **outer** step, independent of the substep.
///
/// Measured, before the fix and after, on this plate at a 300 s outer step:
///
/// ```text
///                 staggered   multirate   analytic   multi err
///   take          303.670     300.033     301.920      1.89
///   take_share      -- same --  301.616    301.920      0.304
/// ```
///
/// The schedule chosen *for* accuracy used to be the worse of the two. It is now fourteen times
/// better than the alternative at the same outer step.
///
/// **What the audit says throughout: nothing.** Both versions pass at around 1e-12, because the
/// total that crossed was always exactly right and only its distribution in time was wrong — and
/// a `Ledger` has no representation for *when*. This is the time-domain twin of the reason
/// `Exchange::audit_transfers` had to become a per-face check in space, and it is why this test
/// compares against a closed form rather than against conservation.
#[test]
fn multirate_now_beats_staggered_because_a_substep_takes_only_its_share() {
    // The closed form of the continuous problem, for reference.
    let p = plate();
    let c = p.heat_capacity().to_si();
    let tau = p.time_constant().to_si();
    let ha = c / tau;
    let analytic = 2.0 / ha * (1.0 - (-600.0 / tau).exp());

    // At a coarse outer step the plate subcycles under multirate and does not under staggered.
    let coarse = Time::from_si(300.0);
    let stag = (rise_after(Schedule::Staggered, coarse, 2) - analytic).abs();
    let multi = (rise_after(Schedule::Multirate, coarse, 2) - analytic).abs();
    assert!(
        multi * 5.0 < stag,
        "subcycling should now be worth something: staggered off by {stag:.4} K, multirate by {multi:.4} K"
    );

    // What the fix does *not* change: the substep, not the outer step, now sets the answer.
    // `substeps_for` gives `ceil(dt / max_stable_dt)`, so a 300 s and a 150 s outer step both
    // land on a 50 s substep — and produce bit-identical results, which they should, because a
    // substep receiving its own share does not care how the outer steps were grouped. That is
    // the property the front-loading version could not have: there, the grouping was the whole
    // story.
    let at_300 = rise_after(Schedule::Multirate, Time::from_si(300.0), 2);
    let at_150 = rise_after(Schedule::Multirate, Time::from_si(150.0), 4);
    assert_eq!(
        at_300, at_150,
        "the same substep should give the same answer however the outer steps are grouped"
    );

    // And a genuinely finer substep is genuinely better, which is what the schedule is for.
    let fine = (rise_after(Schedule::Multirate, Time::from_si(75.0), 8) - analytic).abs();
    assert!(
        fine < multi,
        "a 37.5 s substep should beat a 50 s one: {fine:.4} K against {multi:.4} K"
    );
}

/// A share taken over many substeps empties the channel exactly, leaving nothing stranded.
///
/// The reason `take_share` apportions against the time *remaining* rather than against the whole
/// interval. With `A·dt/T` and both reduced, `A/T` is unchanged, so the last substep receives the
/// remainder and the channel ends empty to the last bit. Against the whole interval instead,
/// `n` shares of `dt/n` leave `O(n·ε·A)` behind — and `audit_transfers` uses an *absolute*
/// tolerance, so at a large enough amount it would refuse a run that was arithmetically fine.
#[test]
fn a_channel_apportioned_over_substeps_ends_exactly_empty() {
    for n in [3usize, 7, 64, 1000] {
        let mut bus = Exchange::new();
        let dt = Time::from_si(1.0);
        bus.covering(dt);
        // A large amount, so an absolute residue would show.
        bus.publish(HEAT, 1.234_567_890_123e9);
        let h = Time::from_si(1.0 / n as f64);
        let mut got = 0.0;
        for _ in 0..n {
            got += bus.take_share(HEAT, h);
        }
        assert_eq!(
            bus.peek(HEAT),
            0.0,
            "{n} substeps left {} on the channel",
            bus.peek(HEAT)
        );
        // The *channel* empties exactly — that is the sharp claim, and it is what stops a
        // residue being stranded where an absolute `audit_transfers` would refuse it. The
        // running total cannot be exact: it is a sum of `n` terms each carrying a division, so
        // it drifts by a few times `n·ε`. Measured 1.9e-11 relative at n = 1000, against
        // 1000·ε = 2.2e-13 for the additions alone; the rest is the divisions.
        assert!(
            (got / 1.234_567_890_123e9 - 1.0).abs() < 1e-9,
            "{n} substeps collected {got}"
        );
        assert!(
            bus.unclaimed().next().is_none(),
            "{n}: something is unclaimed"
        );
    }
}

/// Without an interval the share is the whole, so a domain written against `take_share` works
/// unchanged on a bare `Exchange` and under a schedule that does not subcycle.
#[test]
fn an_unknown_interval_hands_over_everything() {
    let mut bus = Exchange::new();
    bus.publish(HEAT, 42.0);
    assert_eq!(bus.take_share(HEAT, Time::from_si(0.1)), 42.0);
    assert_eq!(bus.peek(HEAT), 0.0);
}

/// **A hypothesis from a design review, tested.**  saves only the
/// temperature, and  reports . Under  the
/// simulation restores before each sweep after the first, so a second sweep would add its
/// losses on top of the first while the temperature is rewound.
#[test]
fn a_lumped_mass_survives_an_iterative_sweep() {
    let mut sim = Simulation::new(Schedule::Iterative {
        max_iter: 3,
        tol: 0.0,
    })
    .conservation_tolerance(1e-9)
    .with(Lamp {
        watts: 2.0,
        reserve: 1e9,
        saved: None,
    })
    .with(plate());
    for _ in 0..4 {
        sim.advance(Time::from_si(10.0))
            .expect("an iterative sweep must not create energy");
    }
}
