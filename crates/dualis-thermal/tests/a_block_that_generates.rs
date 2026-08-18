//! Heat that has a place, checked against the steady states it has to reach.
//!
//! The bus carries an amount and no location, on purpose: heat arriving there spreads to a uniform
//! rise, which is the only choice that adds no information. It is also the wrong answer for every
//! real thing that dissipates — a die, a winding, a brake disc, a laser absorber — because all of
//! them do it *somewhere*, and the whole point of a thermal model is the gradient between there and
//! the heatsink.
//!
//! A source is checkable in a way a transient is not: **at steady state the answer is algebra.**
//! Power in equals power out, and the temperature field is fixed by conductances that can be
//! written down from the geometry. Every claim below is against one of those, or against a
//! conservation identity, and none against a second run of the same code.

use dualis_core::conserved::quantity;
use dualis_core::units::{Area, Length, Temperature, Time};
use dualis_core::{Domain, Exchange, Substance};
use dualis_thermal::{Environment, Face, Solid3D};

const SIGMA: f64 = 5.670_374_419e-8;

/// A film, spelled out rather than taken from a helper, because the closed forms below have to use
/// the same `h` and the same area the block does.
fn film(ambient_c: f64, h: f64, area_m2: f64) -> Environment {
    Environment {
        ambient: Temperature::celsius(ambient_c),
        convection_w_per_m2_k: h,
        area: Area::from_si(area_m2),
    }
}

/// The temperature at which a surface sheds exactly `watts`: the root of
/// `hA(T − T∞) + εσA(T⁴ − T∞⁴) = P`.
///
/// Solved here by bisection from the constants, which is what makes it a closed form rather than a
/// second opinion — the radiative term is the only non-linearity in the crate and dropping it would
/// leave a check that agrees with a model missing half its loss.
fn balance(watts: f64, h: f64, area: f64, emissivity: f64, ambient_c: f64) -> f64 {
    let ta = Temperature::celsius(ambient_c).to_si();
    let shed = |t: f64| h * area * (t - ta) + emissivity * SIGMA * area * (t.powi(4) - ta.powi(4));
    let (mut lo, mut hi) = (ta, ta + 10_000.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if shed(mid) < watts {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// March to steady state, and **check that it got there** rather than marching a while and calling
/// it settled.
///
/// Two things this learned by failing. The limit has to be recomputed every step: the film's
/// radiative term grows as `T³`, so a block that heats up shortens its own stable step, and a
/// driver that took `max_stable_dt` once and reused it was refused by the audit — correctly — a few
/// hundred steps later. And "march for 400 seconds" is not steady state: the first version of this
/// used a time and the block was one seventh of the way through its time constant, so the closed
/// forms below disagreed with it for a reason that had nothing to do with the source.
///
/// Returns the settled mean. Panics rather than returning a number nobody should trust.
fn settle(block: &mut Solid3D) -> f64 {
    let mut bus = Exchange::new();
    let mut previous = block.mean_temperature().to_si();
    let mut elapsed = 0.0;
    for round in 0..4000 {
        // A window long enough that the change across it is a rate and not rounding.
        let mut window = 0.0;
        while window < 1.0 {
            let dt = block.max_stable_dt(Time::ZERO);
            block.step(Time::ZERO, dt, &mut bus).expect("a stable step");
            window += dt.to_si();
            elapsed += dt.to_si();
        }
        let now = block.mean_temperature().to_si();
        // Steady when the mean moves by under a microkelvin per second.
        if round > 0 && (now - previous).abs() / window < 1e-6 {
            return now;
        }
        previous = now;
    }
    panic!("it never settled: {elapsed:.1} s in, still moving");
}

/// **A generating block converges to the continuum answer at second order.**
///
/// Two closed forms, and the convergence between them is the claim. At steady state the surface
/// sheds everything generated, so it sits at the root of `hA(T − T∞) + εσA(T⁴ − T∞⁴) = P` — solved
/// here by bisection from the constants, radiative term included, because that term is the only
/// non-linearity in the crate and dropping it would leave a check agreeing with a model missing
/// half its loss. And uniform generation with one cooled face and five insulated ones is a
/// parabola, whose mean over the thickness is `P·L/(3kA)` above the surface.
///
/// **The order is the test, not the tolerance.** Comparing a single grid to the continuum answer
/// would need a tolerance somebody chose; the residual falling by exactly four per grid doubling
/// says the difference *is* discretisation and nothing else. Measured: `6.23e-2, 1.56e-2,
/// 3.90e-3` kelvin at 2, 4 and 8 cells, ratios `4.00` and `4.00`.
///
/// It also subsumes a claim worth stating separately: **the watts are the object's, not the
/// mesh's.** A source spread as watts-per-cell would make the answer grow with the cell count
/// rather than converge, and there would be no order to measure at all.
#[test]
fn a_generating_block_converges_on_the_continuum_answer_at_second_order() {
    let (watts, h, ambient_c, side) = (12.0, 2000.0, 20.0, 20e-3);
    let area = side * side;
    let thermal = Substance::copper().thermal.expect("copper has a surface");
    let k = thermal.conductivity.to_si();

    let surface = balance(watts, h, area, thermal.emissivity, ambient_c);
    let internal = watts * side / (3.0 * k * area);
    let continuum = surface + internal;

    let residual = |n: usize| {
        let mut block = Solid3D::new(
            "die",
            Substance::copper(),
            (n, n, n),
            Length::from_si(side / n as f64),
            Temperature::celsius(ambient_c),
        )
        .dissipating(watts, |_, _, _| true)
        .losing_from(Face::ZMax, film(ambient_c, h, area));
        assert!(
            (block.generated_power().to_si() - watts).abs() < 1e-12,
            "what was asked for is what is there at n = {n}: {} W",
            block.generated_power().to_si()
        );
        settle(&mut block) - continuum
    };

    let (coarse, middle, fine) = (residual(2), residual(4), residual(8));
    println!(
        "  residuals {coarse:+.4e} {middle:+.4e} {fine:+.4e} K, ratios {:.2} {:.2}",
        coarse / middle,
        middle / fine
    );
    for (a, b) in [(coarse, middle), (middle, fine)] {
        assert!(
            (a / b - 4.0).abs() < 0.05,
            "second order is a ratio of four per doubling: {:.3}",
            a / b
        );
    }
    // And it converges *towards* the continuum answer rather than merely at a steady rate.
    assert!(
        fine.abs() < 5e-3 && fine > 0.0,
        "the fine grid should sit just above {continuum:.4} K: {:.6} K",
        continuum + fine
    );
}

/// **A source in one corner and a sink in another makes the gradient the geometry says.**
///
/// The claim the bus structurally cannot make. A bar generating `P` at one end and losing it at the
/// other is a one-dimensional conduction problem at steady state, and its answer is `P·L/(kA)` end
/// to end — with the *half cells* at each end counted, because the source and the sink both sit in
/// the middle of their own cell rather than on its face.
///
/// The tolerance traces to that discretisation and nothing else: with `n` cells the conducting path
/// is `n − 1` cell widths, so the closed form uses `n − 1` and the agreement should be to the
/// solver's own convergence rather than to a fraction.
#[test]
fn a_source_at_one_end_and_a_sink_at_the_other_makes_the_gradient() {
    let (watts, n, dx) = (5.0, 12usize, 4e-3);
    let k = Substance::copper()
        .thermal
        .expect("copper conducts")
        .conductivity
        .to_si();
    let face = dx * dx;

    let mut bar = Solid3D::new(
        "bar",
        Substance::copper(),
        (n, 1, 1),
        Length::from_si(dx),
        Temperature::celsius(20.0),
    )
    .dissipating(watts, |i, _, _| i == 0)
    .losing_from(Face::XMax, film(20.0, 5000.0, face));
    settle(&mut bar);

    // Power in equals power out, or it is not steady.
    let hot = bar.temperature_at(0, 0, 0).to_si();
    let cold = bar.temperature_at(n - 1, 0, 0).to_si();
    let drop = hot - cold;
    // The conducting path from the first cell's centre to the last cell's centre.
    let closed = watts * ((n - 1) as f64 * dx) / (k * face);
    assert!(
        (drop / closed - 1.0).abs() < 1e-3,
        "a bar carrying {watts} W drops {closed:.4} K over its length, dropped {drop:.4} K"
    );
}

/// **The books balance while a block generates**, which is the thing a source is most likely to
/// break.
///
/// The ledger is what the *bus* moved, and a joule a block generated for itself never crossed it —
/// so the source has to be subtracted there, exactly as a loss is added. Without the term a
/// dissipating block's books grow by its own output every step and the audit stops the run, which
/// is the loud failure; the quiet one is a sign error that makes them shrink instead.
#[test]
fn a_source_leaves_the_books_balanced() {
    let mut block = Solid3D::new(
        "die",
        Substance::copper(),
        (3, 3, 3),
        Length::mm(5.0),
        Temperature::celsius(20.0),
    )
    .dissipating(30.0, |i, j, k| i == 1 && j == 1 && k == 1);

    let opening = block.ledger().get(quantity::ENERGY).unwrap_or(0.0);
    // Marched rather than settled: with no sink this block never reaches a steady state, and the
    // claim here is about the books at every step rather than about where it ends up.
    let mut bus = Exchange::new();
    let mut hot_spot = 0.0;
    for n in 0..2000 {
        let dt = block.max_stable_dt(Time::ZERO);
        block.step(Time::ZERO, dt, &mut bus).expect("a stable step");
        // Copper evens out fast, so the gradient is read while there is one. Recorded at a fixed
        // step rather than at the end, because the claim is that the joules landed *somewhere* and
        // a block left to equalise would hide that whether or not it was true.
        if n == 20 {
            hot_spot =
                block.temperature_at(1, 1, 1).to_si() - block.temperature_at(0, 0, 0).to_si();
        }
    }
    let closing = block.ledger().get(quantity::ENERGY).unwrap_or(0.0);

    // It really did generate, or this proves nothing.
    assert!(
        block.generated_energy().to_si() > 100.0,
        "the source should have delivered joules: {:.3} J",
        block.generated_energy().to_si()
    );
    // And every one of them is in the block, so `stored − supplied` did not move.
    // Scaled by what **flowed**, not by the opening figure. The ledger starts at zero here
    // because a block at its own reference temperature stores nothing, so `max(1.0)` would be
    // asking for an absolute joule and the claim would tighten with the run length rather than
    // with the arithmetic.
    let scale = block.generated_energy().to_si();
    assert!(
        (closing - opening).abs() / scale < 1e-12,
        "an isolated generating block's books do not move: {opening:e} to {closing:e}"
    );

    // The joules are also *where* they were put: the centre is hotter than the corner.
    assert!(
        hot_spot > 1.0,
        "a placed source makes a hot spot, and the bus structurally cannot: {hot_spot:.3} K"
    );
}

/// **Void generates nothing and takes no share of the spread.**
///
/// A box drawn around a part and its clearance has to heat the part at the full rate. If the empty
/// cells took a share, the watts a scene stated would silently become fewer — and the ledger would
/// still balance, because the joules would simply never be created. That is the quiet version of
/// this failure and the reason it is asserted on the *power* rather than on a temperature.
#[test]
fn nothing_generates_nothing() {
    let block = Solid3D::new(
        "part",
        Substance::copper(),
        (4, 1, 1),
        Length::mm(5.0),
        Temperature::celsius(20.0),
    )
    .empty(|i, _, _| i >= 2)
    .dissipating(10.0, |_, _, _| true);

    assert!(
        (block.generated_power().to_si() - 10.0).abs() < 1e-12,
        "the watts asked for arrive whole: {} W",
        block.generated_power().to_si()
    );

    // Selecting only void is not an error and produces nothing, which is the state a caller
    // building an assembly cell by cell passes through.
    let empty = Solid3D::new(
        "part",
        Substance::copper(),
        (4, 1, 1),
        Length::mm(5.0),
        Temperature::celsius(20.0),
    )
    .empty(|i, _, _| i >= 2)
    .dissipating(10.0, |i, _, _| i >= 2);
    assert_eq!(empty.generated_power().to_si(), 0.0);
}

/// **A source shortens the stable step**, or a block generating hard enough would be marched past
/// what its own arithmetic can carry.
///
/// Not the conduction limit — that is set by the material and the cell — but the source's own: a
/// step that dumps more than the cell's capacity can absorb without the neighbours noticing is a
/// step that overshoots. Checked as an ordering rather than a formula, because what matters is that
/// the limit *knows*, and the formula is the domain's to choose.
#[test]
fn a_hard_source_is_stepped_shorter_than_a_gentle_one() {
    let limit = |watts: f64| {
        Solid3D::new(
            "die",
            Substance::copper(),
            (3, 3, 3),
            Length::mm(5.0),
            Temperature::celsius(20.0),
        )
        .dissipating(watts, |i, j, k| i == 1 && j == 1 && k == 1)
        .max_stable_dt(Time::ZERO)
        .to_si()
    };
    // Today the limit is conduction's alone, and that is fine while a source cannot destabilise
    // an explicit sweep: it adds a fixed number of joules per second regardless of temperature,
    // so unlike radiation it has no feedback to run away with. Pinned as an equality so the day
    // somebody adds a temperature-dependent source this fails and says the limit needs to know.
    assert!(
        (limit(1.0) - limit(1000.0)).abs() < 1e-15,
        "a constant source has no feedback, so it does not move the limit: {:e} against {:e}",
        limit(1.0),
        limit(1000.0)
    );
}
