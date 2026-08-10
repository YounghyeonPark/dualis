//! What current can this joint carry, and how many units pass?
//!
//! ```text
//! cargo run --release --example busbar_rating            # the study
//! cargo run --release --example busbar_rating out.svg    # and the plots
//! ```
//!
//! The other examples ask whether the physics is right. This one asks the question an engineer is
//! actually paid to answer, and it runs the whole way from geometry to a production yield:
//!
//! ```text
//!   1  geometry  → resistance      a bolted joint is a constriction, and no formula gives its R
//!   2  network   → thermal path    W/K out of the joint, from the solved balance
//!   3  fixed point                 R rises with T, T rises with I²R — solve them together
//!   4  bisection → the rating      the current that lands the joint on its limit
//!   5  runaway                     the current at which the feedback stops converging at all
//!   6  Monte Carlo → yield         with the tolerances a factory actually holds
//! ```
//!
//! # The one thing that makes this trustworthy
//!
//! Every step above has a closed form to check it against, and the last two share one:
//!
//! ```text
//!   ΔT(I) = I²R₂₀ / (g − I²R₂₀α)          the fixed point, solved exactly
//!   I_runaway = √( g / (R₂₀α) )           where that denominator reaches zero
//! ```
//!
//! So the iteration is checked against the algebra, the algebra's pole is checked against
//! `Winding::runaway_current`, and the resistance that feeds both is checked against `ρL/A` for a
//! plain bar and against the classical constriction resistance `ρ/2a` for the joint. A design
//! study whose numbers cannot be checked is a design study that reports whatever it computed.
//!
//! # Where a reduction is used, and why that is the professional choice
//!
//! The joint's resistance is solved as a **field**, because a constriction has no `ρL/A`. The
//! thermal path is a **four-node network**, because the sweep below evaluates it about 900 times
//! and a 3D conduction solve of a body whose Biot number is 0.003 would be spending an hour to
//! confirm what a lumped model gets right in microseconds.
//!
//! Solve the hard geometry once; sweep the reduction. Using the expensive model everywhere is not
//! rigour, it is a failure to know which question is being asked.

use dualis::prelude::*;

mod common;
use common::svg::{rgb, ticks, Plot};
use common::{check, check_between, heading};

/// Copper at 20 °C.
const RHO_20: f64 = 1.724e-8;
/// Per kelvin. Copper's resistivity coefficient.
const ALPHA: f64 = 0.00393;

/// The joint: 25 mm along the current, on a 20 x 10 mm busbar section, at 1 mm cells.
const NX: usize = 25;
const NY: usize = 20;
const NZ: usize = 10;
const DX: f64 = 1e-3;
/// Radius of the contact patch through the interface, in millimetres. A bolted lap joint touches
/// over a fraction of its overlap, and that fraction is the joint.
const CONTACT_R: f64 = 4.0;

/// What the joint is allowed to reach. A tinned joint under a plated bolt is usually held to this.
const LIMIT_C: f64 = 105.0;
/// Design ambient.
const AMBIENT_C: f64 = 40.0;

fn main() {
    // ================================================================ 1. geometry -> resistance
    heading("1. The joint's resistance, solved from its geometry");

    let plain = conductor(false);
    let bulk = plain.resistance().to_si();
    let length = NX as f64 * DX;
    let area = (NY * NZ) as f64 * DX * DX;
    check(
        "a plain bar is rho L / A",
        bulk,
        RHO_20 * length / area,
        1e-12,
        "ohm",
    );

    // **Verify the solver against a known limit before trusting it on the real part.** Maxwell's
    // constriction resistance for a circular aperture joining two half-spaces is `rho/2a`, and it
    // is a *limit*: it assumes the conductor is unbounded either side. A finite section confines
    // the current before the aperture does, so the aperture adds less — and the deficit has to
    // vanish as the section grows. That is the check, and it is a convergence rather than a value.
    // Two discretisation limits are in play and they push opposite ways, so this is a band and a
    // *direction* rather than an equality:
    //
    //   section >> a    a finite section confines the current before the aperture does, so the
    //                   aperture adds less than the half-space formula — low, and rising
    //   a >> dx         a staircase aperture of a few cells adds numerical resistance — high
    //
    // At 3 mm on a 1 mm grid they partly cancel. Claiming a tight equality here would be claiming
    // the cancellation is exact, which it is not; what is checked is that the approach is monotone
    // toward the limit and lands near it. A solver with a wrong face conductance or a missing
    // transverse coupling misses that by far more than the band.
    const PROBE_R: f64 = 3.0;
    let maxwell = RHO_20 / (2.0 * PROBE_R * DX);
    println!(
        "  {:<34} {:>10.3} uohm   for a {:.0} mm patch between half-spaces",
        "Maxwell's rho / 2a",
        maxwell * 1e6,
        2.0 * PROBE_R
    );
    let mut approach = Vec::new();
    for section in [12usize, 16, 24, 32] {
        let widened = constriction(section, PROBE_R);
        approach.push((section, widened / maxwell));
        println!(
            "    {:>2} x {:>2} mm section ({:>4.1}a)   {:>7.4} uohm   {:>6.3}x the limit",
            section,
            section,
            section as f64 / PROBE_R,
            widened * 1e6,
            widened / maxwell
        );
    }
    for pair in approach.windows(2) {
        assert!(
            pair[1].1 > pair[0].1,
            "widening the section must raise the constriction toward the limit:              {:?} then {:?}",
            pair[0],
            pair[1]
        );
    }
    let (_, closest) = *approach.last().unwrap();
    check_between("and the widest lands on it", closest, 0.85, 1.15, "x");

    // Now the part that is actually being designed.
    let joint = conductor(true);
    let r20 = joint.resistance().to_si();
    assert!(joint.converged(), "residual {:.2e}", joint.residual());
    println!(
        "  {:<34} {:>10.3} uohm   {NX} x {NY} x {NZ} mm bar, bulk only",
        "the busbar alone",
        bulk * 1e6
    );
    println!(
        "  {:<34} {:>10.3} uohm   with the bolted interface",
        "the joint as built",
        r20 * 1e6
    );
    println!(
        "  {:<34} {:>10.1} %      of the joint's resistance is the contact",
        "so the contact is",
        100.0 * (r20 - bulk) / r20
    );
    check_between(
        "the contact is a real part of the budget",
        100.0 * (r20 - bulk) / r20,
        10.0,
        80.0,
        "%",
    );

    // ================================================================ 2. the thermal path
    heading("2. The path out, from the network's own solved balance");
    let (net, hot) = thermal_path();
    // At a representative loss, because a radiative term makes the path mildly non-linear and
    // `path_conductance` takes the slope where it is asked. FRICTION 20: the hand-assembled
    // convection-only figure was 4% wrong on a real motor.
    let g = net
        .path_conductance(hot, Power::w(0.5))
        .expect("the joint has a path to ambient")
        .to_si();
    println!("  {:<34} {:>10.4} W/K  at 0.5 W", "joint to ambient", g);
    println!(
        "  {:<34} {:>10.4}      well under 0.1, so a lumped node is honest",
        "Biot number",
        net.biot_number(hot).unwrap_or(f64::NAN)
    );
    check_between(
        "the lumped model is valid here",
        net.biot_number(hot).unwrap_or(1.0),
        0.0,
        0.1,
        "Bi",
    );

    // ================================================================ 3. the fixed point
    heading("3. The electro-thermal fixed point, and the algebra it must match");
    let rise_at = |amps: f64, r20: f64, g: f64| -> (f64, u32) {
        // R rises with T; T rises with I²R. Neither is known without the other, so iterate.
        // `FRICTION.md` 21: this loop is eight lines every consumer of this library writes.
        let (mut rise, mut iterations) = (0.0f64, 0);
        for _ in 0..200 {
            let r = r20 * (1.0 + ALPHA * (AMBIENT_C + rise - 20.0));
            let next = amps * amps * r / g;
            iterations += 1;
            if (next - rise).abs() < 1e-12 * next.max(1.0) {
                rise = next;
                break;
            }
            rise = next;
        }
        (rise, iterations)
    };
    // The same fixed point, solved exactly. R(T) is linear in T, so the loop above is a linear
    // recurrence and its limit is available in closed form.
    let rise_exact = |amps: f64, r20: f64, g: f64| {
        let r_ambient = r20 * (1.0 + ALPHA * (AMBIENT_C - 20.0));
        let p = amps * amps;
        p * r_ambient / (g - p * r20 * ALPHA)
    };

    for amps in [100.0, 200.0, 300.0] {
        let (rise, iters) = rise_at(amps, r20, g);
        println!(
            "  {:>5.0} A   rise {:>7.2} K   joint at {:>6.1} C   ({iters} iterations)",
            amps,
            rise,
            AMBIENT_C + rise
        );
        check(
            &format!("{amps:.0} A against the closed form"),
            rise,
            rise_exact(amps, r20, g),
            1e-9,
            "K",
        );
    }

    // ================================================================ 4. the rating
    heading("4. The rating: bisection onto the limit");
    let target = LIMIT_C - AMBIENT_C;
    let rating = bisect(|a| rise_at(a, r20, g).0 - target, 1.0, 2000.0);
    let density = rating / (area * 1e6);
    println!("  {:<34} {:>10.1} A", "continuous rating", rating);
    println!(
        "  {:<34} {:>10.2} A/mm2   ({LIMIT_C:.0} C limit, {AMBIENT_C:.0} C ambient)",
        "current density", density
    );
    check(
        "the rating lands on the limit",
        AMBIENT_C + rise_at(rating, r20, g).0,
        LIMIT_C,
        1e-9,
        "C",
    );
    check_between(
        "and it is a plausible busbar density",
        density,
        1.0,
        8.0,
        "A/mm2",
    );

    // ================================================================ 5. runaway
    heading("5. How much margin there is before the feedback wins");
    // A winding is a resistor with a temperature coefficient, which is exactly what this joint is.
    let coil = Winding::of_resistance(
        "joint",
        Resistance::from_si(r20),
        ALPHA,
        Temperature::celsius(AMBIENT_C),
    )
    .driven_at(Current::a(rating));
    let runaway = coil
        .runaway_current(Conductance::from_si(g))
        .expect("a current-driven joint has a threshold");
    // Where the closed form's denominator reaches zero. Computed here, not taken from the domain.
    let exact = (g / (r20 * ALPHA)).sqrt();
    println!("  {:<34} {:>10.1} A", "thermal runaway at", runaway.to_si());
    println!(
        "  {:<34} {:>10.2}x the rating",
        "margin",
        runaway.to_si() / rating
    );
    check(
        "runaway = sqrt(g / R20 alpha)",
        runaway.to_si(),
        exact,
        1e-12,
        "A",
    );
    check_between(
        "the design has real margin",
        runaway.to_si() / rating,
        1.3,
        5.0,
        "x",
    );

    // ================================================================ 6. yield
    heading("6. Yield, over the tolerances a factory actually holds");
    // A bolted joint's contact resistance is the loose number in any busbar: torque, plating
    // thickness and surface finish all move it, and +-15% is a tight assembly rather than a
    // typical one. Ambient is a duty-cycle question rather than a manufacturing one.
    const SPREAD: f64 = 0.15;
    const AMBIENT_LO: f64 = 35.0;
    const AMBIENT_HI: f64 = 50.0;
    let units = 20_000;
    let study = Ensemble::new(0x8115_ba12, units);

    // One unit, at a given current. Two independent draws from a generator keyed by the unit's
    // index — `Rng::for_index`, which `Ensemble` supplies — so this study is bit-identical
    // however many threads run it and whatever order they finish in.
    let unit_at = |amps: f64, mut rng: Rng| {
        let r_unit = r20 * (1.0 + rng.range(-SPREAD, SPREAD));
        let ambient = rng.range(AMBIENT_LO, AMBIENT_HI);
        let mut rise = 0.0f64;
        for _ in 0..200 {
            let r = r_unit * (1.0 + ALPHA * (ambient + rise - 20.0));
            let next = amps * amps * r / g;
            if (next - rise).abs() < 1e-12 * next.max(1.0) {
                rise = next;
                break;
            }
            rise = next;
        }
        ambient + rise
    };

    let temps = study.run(|_, rng| unit_at(rating, rng));
    let passed = temps.iter().filter(|t| **t <= LIMIT_C).count();
    let yield_pct = 100.0 * passed as f64 / units as f64;
    let hottest = temps.iter().copied().fold(f64::MIN, f64::max);
    let estimate = study
        .estimate(|_, rng| unit_at(rating, rng))
        .expect("20 000 samples is not empty");

    println!(
        "  {:<34} {:>10} units, R +-{:.0}%, ambient {AMBIENT_LO:.0}-{AMBIENT_HI:.0} C",
        "sampled",
        units,
        SPREAD * 100.0
    );
    println!(
        "  {:<34} {:>10.2} C +- {:.3} (standard error)",
        "mean joint temperature", estimate.mean, estimate.standard_error
    );
    println!("  {:<34} {:>10.2} C", "hottest unit", hottest);
    println!(
        "  {:<34} {:>10.2} %",
        "pass at the rated current", yield_pct
    );

    // **The rating was set at nominal, so about half the units miss it.** That is the finding, not
    // a defect: a rating computed from nominal values is a coin toss in production, and this is
    // what the Monte Carlo is for.
    check_between("yield at the nominal rating", yield_pct, 20.0, 60.0, "%");

    // So derate until the yield is what a specification would ask for.
    // Bisected on the *same* 20 000 units at every trial current, not on a fresh draw each time.
    // A resampled objective is noisy, and bisection on a noisy objective converges to wherever
    // the noise happened to change sign — which looks like an answer and is a coin toss.
    let yield_at = |amps: f64| {
        study
            .run(|_, rng| unit_at(amps, rng))
            .into_iter()
            .filter(|t| *t <= LIMIT_C)
            .count() as f64
            / units as f64
    };
    // Negated, because yield **falls** with current while the temperature rise climbs — see
    // `bisect`, which now refuses the wrong direction instead of returning a bound.
    let derated = bisect(|a| 0.999 - yield_at(a), 1.0, rating);
    println!(
        "  {:<34} {:>10.1} A   ({:.0}% of the nominal rating)",
        "for 99.9% yield, derate to",
        derated,
        100.0 * derated / rating
    );
    check_between(
        "the derating is real but not ruinous",
        derated / rating,
        0.75,
        0.98,
        "x",
    );

    if let Some(path) = common::output_path() {
        common::write(
            &path,
            &draw(r20, g, rating, runaway.to_si(), derated, &temps, &joint),
        );
    }
}

/// The constriction resistance alone, for a square section of `n` mm — the joint's total minus
/// the same block without an interface. Used to show the solver converging on Maxwell's limit.
fn constriction(n: usize, radius: f64) -> f64 {
    let build = |joint: bool| {
        let mut c = Conductor::new(
            "probe",
            (NX, n, n),
            Length::from_si(DX),
            Resistivity::ohm_m(RHO_20),
            Voltage::mv(1.0),
        );
        if joint {
            let centre = (n as f64 - 1.0) / 2.0;
            c.set_region(
                |i, j, k| {
                    let r = ((j as f64 - centre).powi(2) + (k as f64 - centre).powi(2)).sqrt();
                    i == NX / 2 && r > radius
                },
                Resistivity::ohm_m(RHO_20 * 1e12),
            );
            c.solve(1e-13);
        }
        c.resistance().to_si()
    };
    build(true) - build(false)
}

/// A copper block carrying current along x, optionally with a bolted interface through it.
fn conductor(with_joint: bool) -> Conductor {
    let mut c = Conductor::new(
        "joint",
        (NX, NY, NZ),
        Length::from_si(DX),
        Resistivity::ohm_m(RHO_20),
        Voltage::mv(1.0),
    );
    if with_joint {
        // An insulating sheet across the middle with a circular contact patch through it, which
        // is what a bolted joint is: two bars touching over a fraction of their overlap.
        let (cy, cz) = ((NY as f64 - 1.0) / 2.0, (NZ as f64 - 1.0) / 2.0);
        c.set_region(
            |i, j, k| {
                let r = ((j as f64 - cy).powi(2) + (k as f64 - cz).powi(2)).sqrt();
                i == NX / 2 && r > CONTACT_R
            },
            Resistivity::ohm_m(RHO_20 * 1e12),
        );
        c.solve(1e-13);
    }
    c
}

/// The joint as a lumped node losing to still air through a mounting bar.
fn thermal_path() -> (ThermalNetwork, Node) {
    let mut net = ThermalNetwork::new("busbar");
    let volume = Volume::from_si(NX as f64 * NY as f64 * NZ as f64 * DX * DX * DX);
    let skin =
        Area::from_si(2.0 * ((NX * NY) as f64 + (NY * NZ) as f64 + (NX * NZ) as f64) * DX * DX);
    let hot = net.node_losing_to(
        "joint",
        Substance::copper(),
        volume,
        Length::from_si(NZ as f64 * DX / 2.0),
        Temperature::celsius(AMBIENT_C),
        Environment::still_air(Temperature::celsius(AMBIENT_C), skin),
    );
    net.absorbing(hot).expect("the joint takes the I2R");
    (net, hot)
}

/// Bisect an **increasing** function for its root. Deterministic and bounded.
///
/// The bracket is asserted rather than assumed. Handed a decreasing function this used to walk
/// the wrong way and return `lo` — a number in range, plausible, and wrong, which is the failure
/// mode a silent helper has. The yield curve falls with current where the temperature curve
/// rises, and that is exactly how it happened.
fn bisect(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    assert!(
        f(lo) < 0.0 && f(hi) > 0.0,
        "bisect needs f(lo) < 0 < f(hi); got f({lo}) = {} and f({hi}) = {}.          A decreasing function has to be negated.",
        f(lo),
        f(hi)
    );
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// The design curve, the current density in the joint, and the yield histogram.
fn draw(
    r20: f64,
    g: f64,
    rating: f64,
    runaway: f64,
    derated: f64,
    temps: &[f64],
    joint: &Conductor,
) -> String {
    let (w, h) = (980.0, 690.0);
    let i_max = runaway * 0.97;
    let curve: Vec<(f64, f64)> = (1..=200)
        .map(|k| {
            let a = i_max * k as f64 / 200.0;
            let r_ambient = r20 * (1.0 + ALPHA * (AMBIENT_C - 20.0));
            let p = a * a;
            (a, AMBIENT_C + p * r_ambient / (g - p * r20 * ALPHA))
        })
        .collect();

    let mut top = Plot::new(w, h, (0.0, i_max), (AMBIENT_C, LIMIT_C * 2.2)).viewport(
        78.0,
        58.0,
        w - 140.0,
        280.0,
    );
    top.title("Busbar joint: temperature against current, and where the rating falls");
    top.polyline([(0.0, LIMIT_C), (i_max, LIMIT_C)], &rgb(200, 60, 40), 1.5);
    top.label(
        i_max * 0.12,
        LIMIT_C + 8.0,
        &format!("{LIMIT_C:.0} C limit"),
        11.0,
        &rgb(200, 60, 40),
        "start",
    );
    top.polyline(curve.iter().copied(), &rgb(255, 82, 33), 2.5);
    for (x, colour, text) in [
        (
            derated,
            rgb(60, 130, 90),
            format!("{derated:.0} A  99.9% yield"),
        ),
        (rating, rgb(40, 60, 110), format!("{rating:.0} A  nominal")),
        (
            runaway,
            rgb(120, 120, 120),
            format!("{runaway:.0} A  runaway"),
        ),
    ] {
        top.polyline([(x, AMBIENT_C), (x, LIMIT_C * 2.2)], &colour, 1.0);
        top.label(x, LIMIT_C * 2.0, &text, 10.0, &colour, "middle");
    }
    top.axes(
        &ticks(0.0, i_max, 6),
        &ticks(AMBIENT_C, LIMIT_C * 2.2, 5),
        |v| format!("{v:.0} A"),
        |v| format!("{v:.0} C"),
    );
    top.caption("the curve is the closed-form fixed point; it has a pole at the runaway current, which is why the margin is finite");

    // The current density through the contact patch, on the mid-plane.
    let (nx, ny, nz) = joint.counts();
    let mid = nz / 2;
    let dens: Vec<f64> = (0..ny)
        .flat_map(|j| (0..nx).map(move |i| (i, j)).collect::<Vec<_>>())
        .map(|(i, j)| joint.current_density_magnitude(i, j, mid).to_si())
        .collect();
    let peak = dens.iter().fold(0.0f64, |m, v| m.max(*v)).max(1e-30);
    let mut left =
        Plot::new(w, h, (0.0, nx as f64), (0.0, ny as f64)).viewport(78.0, 404.0, 380.0, 190.0);
    left.raster(nx, ny, (0.0, nx as f64), (0.0, ny as f64), |i, j| {
        common::svg::heat(dens[j * nx + i] / peak)
    });
    left.caption("|J| through the contact patch: the crowding is the constriction resistance");

    // Yield.
    let lo = temps.iter().copied().fold(f64::MAX, f64::min);
    let hi = temps.iter().copied().fold(f64::MIN, f64::max);
    let bins = 44;
    let mut counts = vec![0.0f64; bins];
    for t in temps {
        let k = (((t - lo) / (hi - lo).max(1e-12)) * (bins - 1) as f64).round() as usize;
        counts[k.min(bins - 1)] += 1.0;
    }
    let tallest = counts.iter().fold(0.0f64, |m, v| m.max(*v));
    let mut right =
        Plot::new(w, h, (lo, hi), (0.0, tallest * 1.1)).viewport(540.0, 404.0, w - 600.0, 190.0);
    for (k, n) in counts.iter().enumerate() {
        let a = lo + (hi - lo) * k as f64 / bins as f64;
        let b = lo + (hi - lo) * (k + 1) as f64 / bins as f64;
        let colour = if b <= LIMIT_C {
            rgb(60, 130, 90)
        } else {
            rgb(200, 60, 40)
        };
        right.cell((a, b), (0.0, *n), &colour);
    }
    right.polyline(
        [(LIMIT_C, 0.0), (LIMIT_C, tallest * 1.1)],
        &rgb(40, 40, 40),
        1.5,
    );
    right.axes(
        &ticks(lo, hi, 5),
        &ticks(0.0, tallest * 1.1, 4),
        |v| format!("{v:.0} C"),
        |v| format!("{v:.0}"),
    );
    right.caption("20 000 units at the nominal rating: red is over the limit");

    let mut foot = Plot::new(w, h, (0.0, 1.0), (0.0, 1.0)).viewport(78.0, 668.0, w - 140.0, 1.0);
    foot.footnote(
        "resistance solved as a field, thermal path as a four-node network, yield by a deterministic ensemble — same answer on any thread count",
    );
    common::svg::document(
        w,
        h,
        [
            top.into_body(),
            left.into_body(),
            right.into_body(),
            foot.into_body(),
        ],
    )
}
