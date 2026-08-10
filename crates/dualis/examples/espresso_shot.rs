//! An espresso machine, from the pump to the cup, with the puck's inside visible.
//!
//! ```text
//! cargo run --release --example espresso_shot              # the study
//! cargo run --release --example espresso_shot shot.html    # and the cross-sections, in 3D
//! ```
//!
//! The question a barista is actually asking is not "what happens" but "why did *that* happen and
//! which knob fixes it". So this runs the whole way from the three settings to the thing in the
//! cup, and then opens the basket up:
//!
//! ```text
//!   1  the machine        pump, boiler, and a portafilter with a temperature of its own
//!   2  one shot           17.6 g in, 38 g out, and the flow curve that got there
//!   3  the cross-section  temperature, speed and extraction on a vertical cut
//!   4  the three knobs    grind, temperature, pressure — each swept on its own
//!   5  the fault          a gap at the basket wall, and what it looks like from outside
//!   6  the cold group     what a portafilter left on the counter costs
//! ```
//!
//! # What makes this checkable rather than merely plausible
//!
//! Every number below is either a closed form or a measurement against one:
//!
//! ```text
//!   Q = kA dp / (mu L)             Darcy — reproduced to 1e-15 by the solve
//!   k ∝ d^2,  K ∝ 1/d^2           Kozeny-Carman and sphere diffusion, independently derived
//!   t_transit = eps L / u          the pore velocity, not the Darcy velocity
//!   solute in = solute out         the ledger, to 1e-14 over 17 000 steps
//!   dT <= C_wall dT_wall / (m c)   the cold group cannot take more than it can hold
//! ```
//!
//! # The one thing here that is fitted
//!
//! Two numbers, and they are named: [`Grind::FINES_RATIO`] and `Bed::solute_diffusivity`. Both
//! were set from a single conventional shot. Everything else — how the flow splits, how the front
//! moves, how a channel starves the bed behind it — is a consequence, and step 4 sweeps ranges
//! neither of them was fitted over.

mod common;

use common::{check, check_between, heading};
use dualis::prelude::*;
use dualis_view::report;

/// A 58 mm basket in a 66 mm box on a 2 mm grid: a 20 mm bed with a 4 mm metal jacket.
///
/// The jacket is the point of the extra four millimetres. A basket that fills its grid has metal
/// only in the corners — the right heat capacity in the wrong shape — and a cut through the axis
/// crosses none of it, so the temperature panel comes out flat and says nothing.
const NX: usize = 33;
const NY: usize = 10;
const NZ: usize = 33;
const DX: f64 = 2e-3;
const RADIUS: f64 = 29e-3;
/// How many steps between re-reads of the stability limit. A countdown rather than a modulo,
/// because `usize::is_multiple_of` is Rust 1.87 and this workspace builds on 1.78.
const RECHECK: usize = 256;

/// What a lockstep run of several baskets produced.
struct Together {
    /// When each basket reached its weight, in the order they were given.
    finished: Vec<f64>,
    /// A cross-section of every basket in every field, on a shared clock.
    frames: Vec<Frame>,
    /// `(t, g/s, grams)` for the first basket, sampled twice a second.
    curve: Vec<(f64, f64, f64)>,
}

/// What the pump holds across the puck, once it is up.
const PUMP_BAR: f64 = 9.0;
/// What the group delivers.
const BREW_C: f64 = 93.0;
/// Inter-particle porosity of a tamped bed. With a 600 kg/m³ particle this is a 330 kg/m³ puck,
/// which is what 17.6 g in 20 mm of a 58 mm basket weighs.
const POROSITY: f64 = 0.45;
/// The target in the cup.
const TARGET_G: f64 = 38.0;

fn machine(grind: Grind, bar: f64, brew_c: f64, porosity: f64) -> Puck {
    Puck::new(
        "basket",
        Basket {
            counts: (NX, NY, NZ),
            cell: Length::from_si(DX),
            radius: Length::from_si(RADIUS),
            grind,
            porosity,
            pressure: Pressure::from_si(bar * 1e5),
            temperature: Temperature::celsius(brew_c),
            ..Basket::espresso()
        },
    )
}

/// The vertical plane through the axis: `nz = 1` at the mid-depth.
///
/// **A cross-section, not a volume**, and that is what makes it affordable to carry one on every
/// frame. A 29x10x29 volume is 8410 numbers; nine of them on fifty frames is 45 MB of a document
/// nobody can open, to animate three baskets. The plane is 290, and a plane through the axis is
/// what the question was about.
fn cut() -> Extent {
    let mid = (NZ as f64 * 0.5) * DX;
    Extent::new(
        LengthVec::from_si(glam::DVec3::new(0.0, 0.0, mid)),
        LengthVec::from_si(glam::DVec3::new(NX as f64 * DX, NY as f64 * DX, mid)),
        NX,
        NY,
        1,
    )
}

/// Pull three baskets to the same weight at once, capturing a cross-section as they go.
///
/// Together rather than one after another, on a **shared clock and a shared step**, so a frame
/// compares three baskets at the same instant. A basket that has reached its weight stops — which
/// is not a modelling convenience but the thing being shown: the channelled one finishes first,
/// and after that its picture is the picture of a shot that is over.
///
/// The step is the smallest any of them will take. Three different steps would put the three on
/// three different clocks and there would be nothing to compare.
fn pull_together(
    baskets: &mut [(&'static str, Puck)],
    target_g: f64,
    cutoff_s: f64,
    every_s: f64,
) -> Together {
    let mut bus = Exchange::new();
    let mut t = 0.0;
    let mut done = vec![f64::NAN; baskets.len()];
    let mut frames = Vec::new();
    let mut curve = Vec::new();
    let mut next_frame = 0.0;
    let mut next_sample = 0.0;
    let mut dt = baskets
        .iter()
        .map(|(_, p)| p.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5)
        .fold(f64::INFINITY, f64::min);
    let mut until_recheck = RECHECK;

    loop {
        if t >= next_frame {
            let mut panels = Vec::new();
            for (label, p) in baskets.iter() {
                for (what, name) in [
                    (Observable::Extraction, "extraction"),
                    (Observable::Concentration, "dissolved"),
                    (Observable::Temperature, "temperature"),
                ] {
                    panels.push(sample_field(
                        format!("{label} - {name}"),
                        &p.field(what),
                        cut(),
                        Pose::IDENTITY,
                        Time::from_si(t),
                    ));
                }
            }
            let readings = baskets
                .iter()
                .flat_map(|(label, p)| {
                    vec![
                        Reading::new(*label, "in the cup", p.delivered().to_si() * 1000.0, "g"),
                        Reading::new(*label, "yield", p.yield_fraction() * 100.0, "%"),
                    ]
                })
                .collect();
            frames.push(Frame {
                time_s: t,
                panels,
                readings,
            });
            next_frame += every_s;
        }
        if t >= next_sample {
            let (_, first) = &baskets[0];
            curve.push((
                t,
                first.flow_rate().to_si() * 1000.0,
                first.delivered().to_si() * 1000.0,
            ));
            next_sample += 0.5;
        }

        let mut any = false;
        for (n, (_, p)) in baskets.iter_mut().enumerate() {
            if done[n].is_nan() {
                if p.delivered().to_si() * 1000.0 >= target_g {
                    done[n] = t;
                    continue;
                }
                p.step(Time::from_si(t), Time::from_si(dt), &mut bus)
                    .expect("stable");
                any = true;
            }
        }
        if !any || t >= cutoff_s {
            for d in done.iter_mut() {
                if d.is_nan() {
                    *d = t;
                }
            }
            break;
        }
        t += dt;
        until_recheck -= 1;
        if until_recheck == 0 {
            dt = baskets
                .iter()
                .map(|(_, p)| p.max_stable_dt(Time::from_si(t)).to_si() * 0.5)
                .fold(f64::INFINITY, f64::min);
            until_recheck = RECHECK;
        }
    }
    Together {
        finished: done,
        frames,
        curve,
    }
}

/// Pull one to a weight, sampling the flow as it goes. Returns the elapsed time and the curve.
fn pull(p: &mut Puck, target_g: f64, cutoff_s: f64) -> (f64, Vec<(f64, f64, f64)>) {
    let mut bus = Exchange::new();
    let mut t = 0.0;
    let mut dt = Time::from_si(p.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let mut curve = Vec::new();
    let mut next = 0.0;
    let mut until_recheck = RECHECK;
    while p.delivered().to_si() * 1000.0 < target_g && t < cutoff_s {
        p.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
        until_recheck -= 1;
        if until_recheck == 0 {
            // The step limit moves as the bed heats and the viscosity with it.
            dt = Time::from_si(p.max_stable_dt(Time::from_si(t)).to_si() * 0.5);
            until_recheck = RECHECK;
        }
        if t >= next {
            curve.push((
                t,
                p.flow_rate().to_si() * 1000.0,
                p.delivered().to_si() * 1000.0,
            ));
            next += 0.5;
        }
    }
    (t, curve)
}

fn main() {
    // ======================================================================== 1. the machine
    heading("1. The machine, and what each part of it is");

    let basket = machine(Grind::espresso(), PUMP_BAR, BREW_C, POROSITY);
    let dose = basket.dose().to_si() * 1000.0;
    println!(
        "  {:<38} {:>10.2} g      from (1-eps) x 600 kg/m3 x the bed volume",
        "dose, which nothing sets directly", dose
    );
    check_between("dose for a 58 mm basket 20 mm deep", dose, 16.0, 19.0, "g");
    println!(
        "  {:<38} {:>10.2} g      at 30% soluble, the ceiling on any yield",
        "extractable",
        basket.extractable().to_si() * 1000.0
    );

    // Darcy in closed form, against the solve. This is the anchor: if it holds, the geometry, the
    // harmonic face means and the half-cell boundaries are all right, and everything downstream
    // rides on a flow field that is correct rather than merely smooth.
    let mu = Liquid::water()
        .viscosity(Temperature::celsius(BREW_C))
        .to_si();
    let k = Grind::espresso().permeability(POROSITY).to_si();
    let open = open_area(&basket);
    let closed = k * open * PUMP_BAR * 1e5 / (mu * NY as f64 * DX);
    println!(
        "  {:<38} {:>10.1} mm2   the {:.0} mm basket, as this grid resolves it against pi r^2 = {:.0}",
        "open cross-section",
        open * 1e6,
        RADIUS * 2000.0,
        std::f64::consts::PI * RADIUS * RADIUS * 1e6
    );
    // In mL/s, because `m3/s` at this scale prints as four zeros and a check nobody can read is
    // a check nobody will notice failing.
    check(
        "Darcy's Q = kA dp / (mu L)",
        basket.flow_rate().to_si() / Liquid::water().density.to_si() * 1e6,
        closed * 1e6,
        1e-9,
        "mL/s",
    );
    println!(
        "  {:<38} {:>10.3e} m2   about {:.0}x below what the sieve diameter predicts",
        "bed permeability",
        k,
        1.0 / (Grind::FINES_RATIO * Grind::FINES_RATIO)
    );
    println!(
        "  {:<38} {:>10.1} s      how long a water molecule is in the bed",
        "pore residence time",
        (NY as f64 * DX)
            / (basket.flow_rate().to_si() / Liquid::water().density.to_si() / open / POROSITY)
    );

    // ============================================================================ 2. one shot
    heading("2. One shot, pulled to a weight the way anybody pulls one");

    // The fault and the cold group are set up now and run **beside** the good one, on the same
    // clock, so that the cross-sections compare three baskets at the same instant rather than
    // three separate runs stitched together afterwards.
    let mut faulted = machine(Grind::espresso(), PUMP_BAR, BREW_C, POROSITY);
    let ring = wall_ring(&faulted);
    faulted.repack(0.60, |i, j, k| ring[i + NX * (j + NY * k)]);
    let mut cold = machine(Grind::espresso(), PUMP_BAR, BREW_C, POROSITY);
    cold.set_wall_temperature(Temperature::celsius(45.0));
    cold.set_inlet_temperature(Temperature::celsius(BREW_C));
    let wall_capacity = wall_heat_capacity(&cold);

    let opening = basket.ledger();
    let mut trio = vec![
        ("even", basket),
        ("wall gap", faulted),
        ("cold basket", cold),
    ];
    let together = pull_together(&mut trio, TARGET_G, 90.0, 2.0);
    let (frames, curve) = (together.frames, together.curve);
    let (elapsed, bad_t, cold_t) = (
        together.finished[0],
        together.finished[1],
        together.finished[2],
    );
    let cold = trio.pop().expect("three").1;
    let faulted = trio.pop().expect("three").1;
    let basket = trio.pop().expect("three").1;
    let shot = basket.shot(Time::from_si(elapsed));
    println!(
        "  {:.1} g in, {:.1} g out, {:.1} s     1:{:.1}",
        dose,
        shot.beverage.to_si() * 1000.0,
        elapsed,
        shot.beverage.to_si() * 1000.0 / dose
    );
    check_between(
        "extraction yield",
        shot.yield_fraction * 100.0,
        17.0,
        23.0,
        "%",
    );
    check_between(
        "TDS, what a refractometer reads",
        shot.tds * 100.0,
        6.0,
        11.0,
        "%",
    );
    check_between(
        "mean temperature in the cup",
        shot.outlet_temperature.to_si() - 273.15,
        88.0,
        94.0,
        "C",
    );

    // The flow is flat, and that is a real prediction rather than an assumption: nothing here
    // holds it flat. Permeability is fixed and the viscosity barely moves at constant temperature,
    // so Darcy gives a constant. A real shot slows because the puck compacts and the fines
    // migrate — neither of which is modelled, and this is where that shows.
    let (first, last) = (curve[1].1, curve[curve.len() - 1].1);
    println!(
        "  {:<38} {:>10.3} g/s   and {:.3} g/s at the end — flat, because nothing here clogs",
        "flow at one second", first, last
    );
    check_between("the flow held steady", last / first, 0.97, 1.03, "x");

    // The books, **measured on the run that just happened** rather than asserted from the type.
    // `books_balance` is a claim a domain makes about itself, and a claim is worth exactly as much
    // as the measurement beside it.
    let closing = basket.ledger();
    for q in [quantity::MASS, quantity::ENERGY] {
        let (a, b) = (opening.get(q).unwrap_or(0.0), closing.get(q).unwrap_or(0.0));
        let scale = opening.scale_of(q).unwrap_or(0.0).max(a.abs()).max(1e-300);
        let drift = (b - a).abs() / scale;
        println!(
            "  {:<38} {drift:>10.2e}      over the whole shot",
            format!("{q} ledger drift")
        );
        assert!(
            drift < 1e-11,
            "{q} drifted by {drift:.3e}, which is a leak rather than rounding"
        );
    }
    println!(
        "  {:<38} {:>10.2e}      flow in against flow out at the two boundaries",
        "flow balance",
        basket.flow_balance()
    );
    assert!(
        basket.flow_balance() < 1e-9,
        "an unconverged solve puts liquid where it did not come from: {:.3e}",
        basket.flow_balance()
    );

    // ===================================================================== 3. the cross-section
    heading("3. Inside the basket, on a vertical cut through the axis");

    // `y` is the flow axis, so a slice at fixed `k` is a vertical cut. Depth-by-depth, because
    // that is the profile the shot is built out of.
    println!("  depth      T (C)    speed (mm/s)   extraction   TDS in the pore (%)");
    for j in 0..NY {
        let mid = NX / 2;
        println!(
            "  {:>4.0} mm   {:>7.2}   {:>11.3}   {:>10.3}   {:>13.2}",
            (j as f64 + 0.5) * DX * 1000.0,
            basket.temperature_at(mid, j, NZ / 2).to_si() - 273.15,
            basket.pore_velocity_at(mid, j, NZ / 2).length() * 1000.0,
            basket.extraction_at(mid, j, NZ / 2),
            basket.concentration_at(mid, j, NZ / 2).to_si() / 10.0,
        );
    }
    let top = basket.extraction_at(NX / 2, 0, NZ / 2);
    let bottom = basket.extraction_at(NX / 2, NY - 1, NZ / 2);
    println!(
        "  {:<38} {:>10.3}      the top gives up {:.0}% more than the bottom",
        "top over bottom",
        top / bottom,
        (top / bottom - 1.0) * 100.0
    );
    assert!(
        top > bottom,
        "water enters clean and leaves loaded, so the inlet end must extract more: {top:.3} \
         against {bottom:.3}"
    );
    println!(
        "  {:<38} {:>10.4}      and none of it is a fault — see radial contrast below",
        "spread across the whole bed",
        basket.unevenness()
    );
    check_between(
        "radial contrast on a good puck",
        basket.radial_contrast(),
        0.97,
        1.03,
        "x",
    );

    // ====================================================================== 4. the three knobs
    heading("4. The three knobs, each swept on its own");

    println!("  grind        time      g/s     yield     TDS");
    let mut grinds = Vec::new();
    for microns in [175.0, 210.0, 250.0, 300.0, 350.0] {
        let mut p = machine(
            Grind::sieved(Length::from_si(microns * 1e-6)),
            PUMP_BAR,
            BREW_C,
            POROSITY,
        );
        let (t, _) = pull(&mut p, TARGET_G, 400.0);
        let s = p.shot(Time::from_si(t));
        println!(
            "  {microns:>4.0} um   {t:>7.1} s   {:>5.2}   {:>6.2}%   {:>5.2}%",
            TARGET_G / t,
            s.yield_fraction * 100.0,
            s.tds * 100.0
        );
        grinds.push((microns, t, s.yield_fraction * 100.0));
    }
    // The time to a fixed weight goes as `1/k`, and `k` goes as `d²`. So halving the grind should
    // quadruple the time — checked against the two ends of the sweep, which span a factor of two.
    let (fine, coarse) = (grinds[0], grinds[4]);
    check(
        "time to 38 g goes as 1/d^2",
        fine.1 / coarse.1,
        (coarse.0 / fine.0).powi(2),
        0.05,
        "x",
    );
    assert!(
        fine.2 > coarse.2,
        "and the finer bed, given that time, extracts more: {:.2}% against {:.2}%",
        fine.2,
        coarse.2
    );

    println!("\n  temperature   time      yield     TDS");
    let mut temps = Vec::new();
    for c in [85.0, 89.0, 93.0, 96.0] {
        let mut p = machine(Grind::espresso(), PUMP_BAR, c, POROSITY);
        let (t, _) = pull(&mut p, TARGET_G, 120.0);
        let s = p.shot(Time::from_si(t));
        println!(
            "  {c:>7.0} C   {t:>7.1} s   {:>6.2}%   {:>5.2}%",
            s.yield_fraction * 100.0,
            s.tds * 100.0
        );
        temps.push((c, t, s.yield_fraction * 100.0));
    }
    // Temperature moves two things at once and they push opposite ways on the *time*: hotter water
    // is thinner, so it runs faster and has less contact — and hotter water extracts faster. The
    // yield rises anyway, which says the Arrhenius term wins. That is worth stating because the
    // naive expectation is that a faster shot extracts less.
    assert!(
        temps[3].1 < temps[0].1,
        "hotter water is thinner and runs faster: {:.1} s at 96 C against {:.1} s at 85 C",
        temps[3].1,
        temps[0].1
    );
    assert!(
        temps[3].2 > temps[0].2,
        "and extracts more anyway, because Arrhenius beats the shorter contact: {:.2}% against \
         {:.2}%",
        temps[3].2,
        temps[0].2
    );
    println!(
        "  {:<38} {:>10.2}%     over 11 K, which is why brew temperature is a fine adjustment",
        "yield across the whole range",
        temps[3].2 - temps[0].2
    );

    println!("\n  pressure      time      yield     TDS");
    let mut pressures = Vec::new();
    for bar in [6.0, 9.0, 12.0] {
        let mut p = machine(Grind::espresso(), bar, BREW_C, POROSITY);
        let (t, _) = pull(&mut p, TARGET_G, 120.0);
        let s = p.shot(Time::from_si(t));
        println!(
            "  {bar:>7.0} bar {t:>7.1} s   {:>6.2}%   {:>5.2}%",
            s.yield_fraction * 100.0,
            s.tds * 100.0
        );
        pressures.push((bar, t, s.yield_fraction * 100.0));
    }
    // Pressure is the one knob that is purely a time knob: it moves the flow linearly and touches
    // nothing else. So the time to a fixed weight is exactly inversely proportional to it, and the
    // yield follows only through the time. That is why "more pressure" is not a way to extract
    // more — it is a way to extract less, faster.
    check(
        "time to 38 g goes as 1/dp",
        pressures[0].1 / pressures[2].1,
        12.0 / 6.0,
        0.03,
        "x",
    );
    assert!(
        pressures[0].2 > pressures[2].2,
        "so raising the pressure lowers the yield: {:.2}% at 6 bar against {:.2}% at 12 bar",
        pressures[0].2,
        pressures[2].2
    );

    // ============================================================================= 5. the fault
    heading("5. A gap at the basket wall, which is what channelling actually is");

    let bad = faulted.shot(Time::from_si(bad_t));
    println!(
        "  {:<38} {:>10.1} s      against {:.1} s for the same weight",
        "time to 38 g", bad_t, elapsed
    );
    println!(
        "  {:<38} {:>10.2}%     against {:.2}%",
        "yield",
        bad.yield_fraction * 100.0,
        shot.yield_fraction * 100.0
    );
    println!(
        "  {:<38} {:>10.3}      one for an even puck; this is the diagnosis",
        "ring over core extraction",
        faulted.radial_contrast()
    );
    assert!(
        bad.yield_fraction < shot.yield_fraction,
        "a channel gets less out of the basket: {:.2}% against {:.2}%",
        bad.yield_fraction * 100.0,
        shot.yield_fraction * 100.0
    );
    check_between(
        "radial contrast with a wall gap",
        faulted.radial_contrast(),
        1.10,
        1.60,
        "x",
    );
    // **The number a barista could actually have seen is the time.** The shot ran short, and that
    // is the only symptom visible without a refractometer — which is why "it gushed" and "it was
    // sour" are the same diagnosis.
    println!(
        "  {:<38} {:>10.0}%     the whole of what is visible from outside the basket",
        "how much sooner it finished",
        (1.0 - bad_t / elapsed) * 100.0
    );

    // ======================================================================== 6. the cold group
    heading("6. A portafilter that was left on the counter");

    let cold_shot = cold.shot(Time::from_si(cold_t));
    let drop = BREW_C - (cold_shot.outlet_temperature.to_si() - 273.15);
    let bound = wall_capacity * (BREW_C - 45.0)
        / (cold_shot.beverage.to_si() * Liquid::water().specific_heat.to_si());
    // **What the wall is here, and what it is not.** The cells outside the inscribed cylinder are
    // the corners of the grid box — a shell a few millimetres thick, 34 J/K. A real portafilter is
    // a few hundred grams of brass at ten times that, and it sits below the basket rather than
    // around it. So this understates the effect by about an order of magnitude, and the bound
    // below is the honest statement: whatever the wall's capacity, the cup cannot lose more than
    // it can absorb.
    println!(
        "  {:<38} {:>10.1} J/K    the basket shell on this grid, not a whole portafilter",
        "what the wall holds", wall_capacity
    );
    println!(
        "  {:<38} {:>10.2} C      and the bound it cannot exceed is {:.2} C",
        "how much cooler the cup came out", drop, bound
    );
    check_between("the cold basket cost the cup", drop, 0.2, bound, "C");
    println!(
        "  {:<38} {:>10.2}%     against {:.2}% into a hot one",
        "yield",
        cold_shot.yield_fraction * 100.0,
        shot.yield_fraction * 100.0
    );
    assert!(
        cold_shot.yield_fraction < shot.yield_fraction,
        "a cold basket under-extracts: {:.3}% against {:.3}%",
        cold_shot.yield_fraction * 100.0,
        shot.yield_fraction * 100.0
    );

    // ================================================================================== output
    if let Some(path) = common::output_path() {
        let mut frames = frames;
        settle_framing(&mut frames);
        println!(
            "\n  {} frames, each a vertical cut through three baskets in nine fields",
            frames.len()
        );
        common::write(
            &path,
            &report::html("An espresso shot, from the pump to the cup", &frames),
        );
    } else {
        println!(
            "\n  Pass a filename ending in .html for the cross-sections:\n    cargo run --release \
             --example espresso_shot shot.html"
        );
    }
}

/// The open cross-section the inscribed cylinder actually has on this grid.
fn open_area(p: &Puck) -> f64 {
    let mut n = 0;
    for k in 0..NZ {
        for i in 0..NX {
            if p.is_packed(i, 0, k) {
                n += 1;
            }
        }
    }
    n as f64 * DX * DX
}

/// The packed cells with a wall neighbour in the radial plane.
fn wall_ring(p: &Puck) -> Vec<bool> {
    let mut ring = vec![false; NX * NY * NZ];
    for k in 0..NZ {
        for j in 0..NY {
            for i in 0..NX {
                if !p.is_packed(i, j, k) {
                    continue;
                }
                ring[i + NX * (j + NY * k)] = i == 0
                    || k == 0
                    || i + 1 == NX
                    || k + 1 == NZ
                    || !p.is_packed(i - 1, j, k)
                    || !p.is_packed(i + 1, j, k)
                    || !p.is_packed(i, j, k - 1)
                    || !p.is_packed(i, j, k + 1);
            }
        }
    }
    ring
}

fn wall_heat_capacity(p: &Puck) -> f64 {
    let al = Substance::aluminium_6061();
    let c = al.thermal.as_ref().expect("aluminium has thermal props");
    let mut cells = 0;
    for k in 0..NZ {
        for j in 0..NY {
            for i in 0..NX {
                if !p.is_packed(i, j, k) {
                    cells += 1;
                }
            }
        }
    }
    cells as f64 * DX.powi(3) * al.density.to_si() * c.specific_heat.to_si()
}
