//! What a packed bed does, checked against things that were true before this code existed.
//!
//! Every test here compares against a closed form, an exact limit or a conservation law. None
//! compares against another run of the same solver, and none compares against a number that was
//! read off this implementation and pasted back in.

use pantometry_core::{Domain, Exchange, Substance};
use pantometry_porous::{Basket, Bed, Grind, Liquid, Observable, Puck};
use pantometry_units::{Length, Pressure, Temperature, Time};

const DX: f64 = 1e-3;

/// How many steps between re-reads of the stability limit.
const RECHECK: usize = 256;

/// A test bed at 1 mm cells, with the basket filling the grid.
///
/// The radius is the inscribed one here — these tests are about the flow and the dissolution, and
/// a jacket would only add cells that do nothing. `a_cold_basket_...` is the exception and asks
/// for one.
fn bed(counts: (usize, usize, usize), porosity: f64, drive: f64) -> Puck {
    bed_at(DX, counts, porosity, drive)
}

/// The same, at a chosen cell size.
fn bed_at(dx: f64, counts: (usize, usize, usize), porosity: f64, drive: f64) -> Puck {
    Puck::new(
        "puck",
        Basket {
            counts,
            cell: Length::from_si(dx),
            radius: Length::from_si(0.5 * counts.0.min(counts.2) as f64 * dx),
            porosity,
            pressure: Pressure::from_si(drive),
            ..Basket::espresso()
        },
    )
}

/// Run until the cup holds `target` grams, or until `cutoff` seconds. Returns the elapsed time.
///
/// **To a weight, not to a clock.** That is how a shot is actually pulled, and it is the only
/// comparison in which the grind's two effects are separable: a finer bed takes longer to reach
/// the same beverage mass, and the extra time is the whole of what it buys.
fn pull(p: &mut Puck, target_g: f64, cutoff_s: f64) -> f64 {
    let mut bus = Exchange::new();
    let mut t = 0.0;
    let mut dt = Time::from_si(p.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let mut until_recheck = RECHECK;
    while p.delivered().to_si() * 1000.0 < target_g && t < cutoff_s {
        p.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
        // The limit moves as the bed heats and the viscosity with it. Re-reading it every few
        // hundred steps is cheap and keeps the step from creeping past a limit that fell.
        until_recheck -= 1;
        if until_recheck == 0 {
            dt = Time::from_si(p.max_stable_dt(Time::from_si(t)).to_si() * 0.5);
            until_recheck = RECHECK;
        }
    }
    t
}

/// The cross-sectional area the inscribed cylinder actually has on the grid.
///
/// Counted rather than computed from `πr²`: a 12-cell circle is not a circle, and the exact
/// answer for *this* bed is the one the discretisation gives. Using `πr²` instead would build a
/// 4% error into every comparison and make the tolerances below unearned.
fn open_area(p: &Puck) -> f64 {
    let (nx, _, nz) = p.counts();
    let mut n = 0;
    for k in 0..nz {
        for i in 0..nx {
            if p.is_packed(i, 0, k) {
                n += 1;
            }
        }
    }
    n as f64 * DX * DX
}

/// **Darcy's law comes out exactly, and nothing states it.**
///
/// `Q = kAΔp/(μL)` for a uniform bed. This is the equivalent of `Conductor` reproducing `ρL/A`,
/// and it is the sharpest single check on the discretisation: a wrong harmonic mean, a wrong
/// half-cell at the boundary, a face counted twice — each would move this by a percent or more,
/// and the agreement here is at the level of the solver's tolerance.
///
/// The length in the closed form is `ny·dx`, **including** the two half-cells at the boundaries,
/// because the inlet and outlet conductances are the half-cell ones. Getting that wrong by half a
/// cell is a 5% error at ten cells and 0.5% at a hundred — visible, and the reason this test also
/// runs at a second depth.
#[test]
fn a_uniform_bed_obeys_darcy_exactly() {
    for depth in [10usize, 20, 40] {
        let porosity = 0.45;
        let p = bed((12, depth, 12), porosity, 9.0e5);
        assert!(p.converged(), "residual {:.3e}", p.residual());

        let mu = Liquid::water()
            .viscosity(Temperature::celsius(93.0))
            .to_si();
        let k = Grind::espresso().permeability(porosity).to_si();
        let length = depth as f64 * DX;
        let closed = k * open_area(&p) * 9.0e5 / (mu * length);
        let measured = p.flow_rate().to_si() / Liquid::water().density.to_si();

        let off = (measured - closed).abs() / closed;
        println!("  {depth:2} cells deep: measured {measured:.6e} m3/s, Darcy {closed:.6e}, off {off:.2e}");
        assert!(
            off < 1e-9,
            "Darcy is exact for a uniform bed; off by {:.3}% at depth {depth}",
            off * 100.0
        );
    }
}

/// **Flow is linear in pressure and inverse in depth, over ranges wide enough to tell.**
///
/// Two exponents rather than two values. A model that had the pressure entering quadratically —
/// which is what happens if somebody reaches for a Forchheimer correction and applies it wrongly —
/// would still give a plausible flow at 9 bar.
#[test]
fn flow_is_linear_in_pressure_and_inverse_in_depth() {
    let two_bar = bed((12, 20, 12), 0.45, 2.0e5).flow_rate().to_si();
    let twelve = bed((12, 20, 12), 0.45, 12.0e5).flow_rate().to_si();
    let ratio = twelve / two_bar;
    assert!(
        (ratio - 6.0).abs() / 6.0 < 1e-9,
        "Q goes as dp: six times the pressure gave {ratio:.6} times the flow"
    );

    let shallow = bed((12, 10, 12), 0.45, 9.0e5).flow_rate().to_si();
    let deep = bed((12, 40, 12), 0.45, 9.0e5).flow_rate().to_si();
    let ratio = shallow / deep;
    assert!(
        (ratio - 4.0).abs() / 4.0 < 1e-9,
        "Q goes as 1/L: a quarter the depth gave {ratio:.6} times the flow"
    );
}

/// **A tracer arrives at the pore velocity, not the Darcy velocity.**
///
/// The commonest error in porous transport, and it is a factor of `1/ε` — 2.2 for an espresso
/// puck. A model that advected at the Darcy velocity would have the front arrive twice too late
/// and would otherwise look entirely right.
///
/// # Making a tracer out of a bed that dissolves everywhere
///
/// Every cell produces solute, so watching the outlet tells you nothing about transport: it rises
/// because the outlet cell is dissolving, not because anything arrived. The two diameters make a
/// clean experiment possible. [`Grind::hydraulic`] sets the sieve diameter — which governs
/// extraction — independently of the hydraulic one, which governs flow. So the bed below the top
/// two layers is given a metre-wide sieve diameter, which stops it extracting entirely
/// (`K ∝ 1/d²`), while its permeability is left exactly as it was.
///
/// What is left is a uniform flow field with a solute source at the inlet end: a step input, whose
/// **half-height arrival** is the transit time to first order regardless of how much the upwind
/// scheme smears it.
#[test]
fn a_tracer_travels_at_the_pore_velocity() {
    let dx = 2e-3;
    let porosity = 0.45;
    let (nx, ny, nz) = (12usize, 12usize, 12usize);
    let mut p = bed_at(dx, (nx, ny, nz), porosity, 9.0e5);
    let normal = Grind::espresso();
    // Inert everywhere: same hydraulics, no dissolution.
    p.regrind(
        Grind::hydraulic(Length::from_si(1.0), normal.hydraulic_diameter()),
        |_, _, _| true,
    );
    // Except the top two layers, which are the source.
    p.regrind(normal, |_, j, _| j < 2);

    let open: f64 = {
        let mut n = 0;
        for k in 0..nz {
            for i in 0..nx {
                if p.is_packed(i, 0, k) {
                    n += 1;
                }
            }
        }
        n as f64 * dx * dx
    };
    let darcy = p.flow_rate().to_si() / Liquid::water().density.to_si() / open;
    let pore_v = darcy / porosity;
    // From the source layer to the outlet.
    let distance = (ny as f64 - 2.0) * dx;
    let predicted = distance / pore_v;
    let if_darcy = predicted / porosity;

    let dt = Time::from_si(p.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let mut bus = Exchange::new();
    let mut t = 0.0;
    let mut history = Vec::new();
    while t < 4.0 * predicted {
        p.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
        history.push((t, p.concentration_at(nx / 2, ny - 1, nz / 2).to_si()));
    }
    let plateau = history.iter().map(|(_, c)| *c).fold(0.0f64, f64::max);
    assert!(plateau > 0.0, "the outlet never saw the tracer");
    let half_at = history
        .iter()
        .find(|(_, c)| *c >= 0.5 * plateau)
        .map(|(t, _)| *t)
        .expect("half of a positive plateau is reached");

    println!("  pore transit predicted {predicted:.3} s, half-height at {half_at:.3} s");
    println!("  the Darcy-velocity error would predict {if_darcy:.3} s");
    let off = (half_at / predicted - 1.0).abs();
    assert!(
        off < 0.30,
        "arrival should be the pore transit time to within the smearing: {half_at:.3} s against          {predicted:.3} s, off by {:.0}%",
        off * 100.0
    );
    // And the two hypotheses are far enough apart here for that to have meant something.
    assert!(
        (if_darcy / predicted - 1.0).abs() > 1.0,
        "this test only distinguishes anything if 1/eps is a large factor: {:.2}x",
        if_darcy / predicted
    );
}

/// **With no flow, extraction relaxes to equilibrium on an exact exponential — and stops short.**
///
/// Zero drive, so nothing moves and nothing leaves. A cell is then a closed two-box system whose
/// solution is analytic in both of its parameters:
///
/// ```text
///   β = m₀/(V_pore·C_sat)                 how much solute the pore liquid can take
///   extraction(t) = (1−e^{−K(1+β)t}) / (1+β)
/// ```
///
/// Two things to get wrong and both are checked. The **rate** is `K(1+β)`, not `K` — faster than
/// the dilute limit, because the driving force falls from both ends. The **ceiling** is
/// `1/(1+β) = 0.577`, not 1 — a stagnant pocket cannot extract past what its own pore liquid will
/// hold, and that ceiling is the entire reason a channelled shot under-extracts.
///
/// Checked at three times, because one point can be hit by a wrong rate and a wrong ceiling
/// together.
#[test]
fn without_flow_extraction_relaxes_to_its_own_equilibrium() {
    let porosity = 0.45;
    let mut p = bed((8, 8, 8), porosity, 0.0);
    let coffee = Bed::roasted_coffee();
    let k = Grind::espresso().extraction_rate(&coffee, Temperature::celsius(93.0));
    // β from the packing, not from the code under test.
    let beta = ((1.0 - porosity) * coffee.solid_density.to_si() * coffee.soluble_fraction)
        / (porosity * coffee.saturation_concentration.to_si());
    let ceiling = 1.0 / (1.0 + beta);
    println!(
        "  beta {beta:.6}, rate {:.6}/s, ceiling {ceiling:.6}",
        k * (1.0 + beta)
    );

    let dt = Time::from_si(p.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let mut bus = Exchange::new();
    let mut t = 0.0;
    let mut checked = 0;

    for target in [5.0, 15.0, 40.0] {
        while t < target {
            p.step(Time::from_si(t), dt, &mut bus).expect("stable");
            t += dt.to_si();
        }
        let extracted = p.extraction_at(4, 4, 4);
        let closed = ceiling * (1.0 - (-k * (1.0 + beta) * t).exp());
        let off = (extracted - closed).abs();
        println!("  t={t:5.2}s  extracted {extracted:.9}  exact {closed:.9}  off {off:.2e}");
        // The step applies the exact solution over each interval and nothing is moving heat, so
        // the temperature is constant and the exponentials compose exactly.
        assert!(
            off < 1e-12,
            "the closed form should be reproduced to rounding: off by {off:.3e}"
        );
        checked += 1;
    }
    assert_eq!(checked, 3);
    // The ceiling is real and it is well short of everything.
    assert!(
        p.extraction_at(4, 4, 4) < ceiling + 1e-9 && p.extraction_at(4, 4, 4) > 0.9 * ceiling,
        "a stagnant cell stalls at {ceiling:.4}, not at 1"
    );
}

/// **Solute mass and enthalpy are conserved to machine precision, over a whole shot.**
///
/// The ledger is arranged so both entries are constant: solute is what is in the particles plus
/// what is in the pores plus what is in the cup, and enthalpy is what the bed holds plus what left
/// minus what was admitted. Anything else is a leak in the discretisation.
///
/// This is what [`Domain::books_balance`] claims, and it is a stronger statement than the
/// whole-simulation audit can make — that one sums every domain before comparing, so a bed losing
/// all of its solute could hide behind a larger domain's rounding.
#[test]
fn the_books_balance_over_a_whole_shot() {
    let mut p = bed((12, 20, 12), 0.45, 9.0e5);
    let before = p.ledger();
    let dt = Time::from_si(p.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let mut bus = Exchange::new();
    let mut t = 0.0;
    let mut steps = 0;
    while t < 25.0 {
        p.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
        steps += 1;
    }
    let after = p.ledger();

    for q in ["mass", "energy"] {
        let (a, b) = (before.get(q).unwrap_or(0.0), after.get(q).unwrap_or(0.0));
        let scale = before.scale_of(q).unwrap_or(0.0).max(a.abs()).max(1e-300);
        let drift = (b - a).abs() / scale;
        println!("  {q:7}: {a:.9e} -> {b:.9e}, drift {drift:.3e} over {steps} steps");
        assert!(
            drift < 1e-11,
            "{q} drifted by {drift:.3e} over {steps} steps, which is a leak rather than rounding"
        );
    }
    // And the run did something: a shot that delivered nothing conserves everything.
    assert!(
        p.yield_fraction() > 0.10,
        "a 25 s shot should have extracted something: yield {:.3}%",
        p.yield_fraction() * 100.0
    );
}

/// **A conventional shot lands where a barista aims it.**
///
/// The headline check, and the one the two calibrated numbers were set from. A real basket at real
/// settings: 58 mm across, 20 mm deep, `eps = 0.45`, 250 um, 9 bar, 93 C, pulled to 38 g.
///
/// Four figures have to come out together, and they are not independent — a model can be made to
/// hit any one of them by moving a constant. The dose follows from the packing, the time from
/// Darcy, the yield from the dissolution and the TDS from both. Two numbers were fitted, and there
/// are four here.
#[test]
fn a_conventional_shot_lands_where_a_barista_aims() {
    let mut p = Puck::new("basket", Basket::espresso());
    let dose = p.dose().to_si() * 1000.0;
    let elapsed = pull(&mut p, 38.0, 60.0);
    let shot = p.shot(Time::from_si(elapsed));

    println!(
        "  {dose:.1} g in, {:.1} g out in {elapsed:.1} s — yield {:.1}%, TDS {:.1}%, out at {:.1} C",
        shot.beverage.to_si() * 1000.0,
        shot.yield_fraction * 100.0,
        shot.tds * 100.0,
        shot.outlet_temperature.to_si() - 273.15
    );
    assert!(
        (16.0..19.0).contains(&dose),
        "a 58 mm basket 20 mm deep at eps=0.45 holds about 17.5 g: {dose:.2} g"
    );
    assert!(
        (20.0..32.0).contains(&elapsed),
        "38 g should take about 25 s: {elapsed:.1} s"
    );
    assert!(
        (0.17..0.23).contains(&shot.yield_fraction),
        "a conventional shot yields 18-22%: {:.2}%",
        shot.yield_fraction * 100.0
    );
    assert!(
        (0.06..0.11).contains(&shot.tds),
        "and reads 7-10% on a refractometer: {:.2}%",
        shot.tds * 100.0
    );
}

/// **A gap at the basket wall lowers the yield and raises the local over-extraction.**
///
/// The commonest real defect: the puck pulls away from the wall and the flow takes the ring. The
/// permeability there is higher by `eps^3/(1-eps)^2` at the two porosities — 2.9x for 0.45 against
/// 0.60 — so the ring passes far more than its share.
///
/// Both shots are pulled **to the same beverage weight**, which is what a barista does and is the
/// only comparison in which this means anything: at equal time the channelled shot simply delivers
/// more liquid, and every number moves for that reason alone.
///
/// # What identifies a channel, and what does not
///
/// The obvious statistic — the spread of the per-cell extraction — is a **bad** detector, and
/// measuring it here is what showed that. An evenly packed bed already has a spread of 0.105,
/// because water that entered clean is loaded by the time it leaves, and the wall channel takes it
/// only to 0.128. The signal is a fifth of the noise it sits on.
///
/// [`Puck::radial_contrast`] is blind to that gradient: the ring and the core span the same depths,
/// so the axial part divides out and what is left is the channel.
///
/// Nor does the peak extraction rise. It **falls**, from 0.936 to 0.836 — because the channelled
/// bed reached the same weight in 15 s instead of 25, so even the over-served ring had less time.
/// "The channel over-extracts" is a statement about the ring relative to the core, not about
/// absolute numbers, and at equal weight the absolute ones go the other way.
#[test]
fn a_wall_channel_lowers_the_yield_and_raises_the_spread() {
    let run = |channel: bool| {
        let (nx, ny, nz) = (20usize, 10usize, 20usize);
        let mut p = bed_at(2e-3, (nx, ny, nz), 0.45, 9.0e5);
        if channel {
            // The outer ring of packed cells: those with a wall neighbour in the radial plane.
            let mut ring = vec![false; nx * ny * nz];
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        if !p.is_packed(i, j, k) {
                            continue;
                        }
                        let edge = i == 0
                            || k == 0
                            || i + 1 == nx
                            || k + 1 == nz
                            || !p.is_packed(i - 1, j, k)
                            || !p.is_packed(i + 1, j, k)
                            || !p.is_packed(i, j, k - 1)
                            || !p.is_packed(i, j, k + 1);
                        ring[i + nx * (j + ny * k)] = edge;
                    }
                }
            }
            p.repack(0.60, |i, j, k| ring[i + nx * (j + ny * k)]);
        }
        let target = p.dose().to_si() * 1000.0 * 2.2;
        let elapsed = pull(&mut p, target, 120.0);
        let mut peak: f64 = 0.0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if p.is_packed(i, j, k) {
                        peak = peak.max(p.extraction_at(i, j, k));
                    }
                }
            }
        }
        (p.yield_fraction(), peak, p.radial_contrast(), elapsed)
    };

    let (even_y, even_peak, even_ratio, even_t) = run(false);
    let (bad_y, bad_peak, bad_ratio, bad_t) = run(true);
    println!(
        "  even:    yield {:5.2}%  peak {even_peak:.3}  ring/core {even_ratio:.4}  {even_t:5.1} s",
        even_y * 100.0
    );
    println!(
        "  channel: yield {:5.2}%  peak {bad_peak:.3}  ring/core {bad_ratio:.4}  {bad_t:5.1} s",
        bad_y * 100.0
    );

    assert!(
        bad_t < even_t,
        "a channelled bed passes the same weight faster: {bad_t:.1} s against {even_t:.1} s"
    );
    assert!(
        bad_y < even_y,
        "and gets less out of the coffee: {:.2}% against {:.2}%",
        bad_y * 100.0,
        even_y * 100.0
    );
    assert!(
        (even_ratio - 1.0).abs() < 0.05,
        "an evenly packed basket extracts its ring and its core alike: {even_ratio:.4}"
    );
    assert!(
        bad_ratio > 1.15,
        "the ring must outrun the core it starved: {bad_ratio:.4} against {even_ratio:.4}"
    );
    // Reported rather than asserted, because it is the finding: the absolute peak goes the other
    // way, and a test that demanded otherwise would have been asserting a plausible story.
    println!(
        "  and the peak fell, {even_peak:.3} -> {bad_peak:.3}, because the shot was {:.0}% shorter",
        (1.0 - bad_t / even_t) * 100.0
    );
}

/// **A cold basket cools the shot, and it does so by the heat capacity it has.**
///
/// The bound is a closed form and does not depend on the solver: however the heat moves, the
/// basket cannot absorb more than `C_wall·ΔT` in total, and the liquid cannot give up more than
/// it carries. So the mean outlet temperature over a shot delivering mass `m` satisfies
///
/// ```text
///   T_in − T_out ≤ C_wall·(T_in − T_wall) / (m·c_water)
/// ```
///
/// Checked as an inequality that a plausible-looking bug would break — a wall with no heat
/// capacity gives zero drop, a wall accidentally treated as a fixed-temperature boundary gives a
/// drop that exceeds the bound and never recovers.
#[test]
fn a_cold_basket_cools_the_shot_by_no_more_than_it_can_hold() {
    // With a jacket, because that is what is being measured. A basket that fills its grid has
    // metal only in the corners: the right heat capacity in the wrong shape.
    //
    // At 2 mm rather than the 1 mm the rest of this file uses. The claim is a capacity bound and
    // does not care about the grid, and the metal is what sets the step here — a jacket at 1 mm
    // made this test alone cost more than the other eleven together.
    let cell = 2e-3;
    let mut p = Puck::new(
        "puck",
        Basket {
            counts: (20, 10, 20),
            cell: Length::from_si(cell),
            radius: Length::from_si(14e-3),
            ..Basket::espresso()
        },
    );
    p.set_wall_temperature(Temperature::celsius(40.0));
    p.set_inlet_temperature(Temperature::celsius(93.0));

    // The wall's capacity, counted from the cells that are wall.
    let (nx, ny, nz) = p.counts();
    let mut wall_cells = 0;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                if !p.is_packed(i, j, k) {
                    wall_cells += 1;
                }
            }
        }
    }
    let steel = Substance::stainless_304();
    let c_wall = wall_cells as f64
        * cell.powi(3)
        * steel.density.to_si()
        * steel.thermal.as_ref().unwrap().specific_heat.to_si();

    let dt = Time::from_si(p.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let mut bus = Exchange::new();
    let mut t = 0.0;
    while t < 25.0 {
        p.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
    }
    let shot = p.shot(Time::from_si(t));
    let mass = shot.beverage.to_si();
    let drop = 93.0 - (shot.outlet_temperature.to_si() - 273.15);
    let bound = c_wall * (93.0 - 40.0) / (mass * Liquid::water().specific_heat.to_si());

    println!("  {wall_cells} wall cells hold {c_wall:.2} J/K; {mass:.4} kg delivered");
    println!("  mean outlet is {drop:.3} K below the inlet; the capacity bound is {bound:.3} K");
    assert!(drop > 0.05, "a cold basket must cool the shot: {drop:.4} K");
    assert!(
        drop < bound,
        "the shot cannot lose more than the wall can absorb: {drop:.3} K against {bound:.3} K"
    );
    // And the wall warmed up, which is where that heat went.
    let mut warmed = false;
    for k in 0..nz {
        if !p.is_packed(0, ny / 2, k) && p.temperature_at(0, ny / 2, k).to_si() > 313.15 + 1.0 {
            warmed = true;
        }
    }
    assert!(warmed, "the heat the shot lost has to be somewhere");
}

/// **Finer grind: slower shot, higher yield — but only when pulled to the same weight.**
///
/// The claim the module docs make, and the correction to the version of it people repeat. At equal
/// *time* a finer bed under-extracts, because four times less water crossed it and what did cross
/// sat there loading up until it stopped taking any more. Measured: 175 um at 25 s reached 12.9%
/// against 350 um's 20.2%, with the fine shot's liquid at 10.7% TDS and going nowhere.
///
/// At equal *weight* the statement is true and both mechanisms are visible at once. The finer bed
/// takes about four times as long — `k` goes as `d^2`, exactly — and extracts more of its dose,
/// because it had both a longer contact time and a four-times-higher rate constant from the same
/// `d^2`.
#[test]
fn a_finer_grind_pulled_to_the_same_weight_runs_slower_and_extracts_more() {
    let run = |scale: f64| {
        let mut p = Puck::new(
            "puck",
            Basket {
                counts: (14, 10, 14),
                radius: Length::from_si(14e-3),
                grind: Grind::espresso().scaled(scale),
                ..Basket::espresso()
            },
        );
        let target = p.dose().to_si() * 1000.0 * 2.2;
        let elapsed = pull(&mut p, target, 400.0);
        (
            elapsed,
            p.yield_fraction() * 100.0,
            p.tds() * 100.0,
            p.delivered().to_si() * 1000.0,
        )
    };

    let (coarse_t, coarse_y, coarse_tds, coarse_g) = run(1.4);
    let (fine_t, fine_y, fine_tds, fine_g) = run(0.7);
    println!("  350 um: {coarse_g:5.2} g in {coarse_t:6.2} s   yield {coarse_y:5.2}%   TDS {coarse_tds:5.2}%");
    println!(
        "  175 um: {fine_g:5.2} g in {fine_t:6.2} s   yield {fine_y:5.2}%   TDS {fine_tds:5.2}%"
    );

    assert!(
        fine_y > coarse_y,
        "the finer bed extracts more of its dose: {fine_y:.2}% against {coarse_y:.2}%"
    );
    // The time ratio is the permeability ratio, which is the square of the diameter ratio.
    let ratio = fine_t / coarse_t;
    let predicted = (1.4 / 0.7f64).powi(2);
    println!("  time ratio {ratio:.3}, k ratio predicts {predicted:.3}");
    assert!(
        (ratio / predicted - 1.0).abs() < 0.15,
        "the time to the same weight follows the permeability: {ratio:.3} against {predicted:.3}"
    );
}

/// **The step past the stability limit is refused rather than run.**
#[test]
fn past_the_transport_limit_is_refused() {
    let mut p = bed((12, 20, 12), 0.45, 9.0e5);
    let limit = p.max_stable_dt(Time::from_si(0.0));
    assert!(limit.to_si().is_finite() && limit.to_si() > 0.0);
    let err = p
        .step(
            Time::from_si(0.0),
            Time::from_si(limit.to_si() * 1.02),
            &mut Exchange::new(),
        )
        .expect_err("2% past the limit must be refused");
    assert_eq!(err.quantity, "transport step");
    // And at the limit it runs.
    p.step(Time::from_si(0.0), limit, &mut Exchange::new())
        .expect("exactly at the limit is stable");
}

/// **The stability limit is where the scheme actually stops being positive.**
///
/// A limit that is merely conservative would pass the test above while hiding a scheme that is
/// stable well past it — and then a later refactor that "corrects" the constant would look
/// harmless. So this measures the other side: a fraction past the reported limit, a concentration
/// goes negative, which is exactly what the positive-coefficient bound is a bound on.
///
/// Run by stepping past the guard deliberately.
#[test]
fn the_limit_is_where_positivity_fails() {
    let (nx, ny, nz) = (12, 20, 12);
    let mut p = bed((nx, ny, nz), 0.45, 9.0e5);
    // Load the bed so there is something to advect.
    let mut bus = Exchange::new();
    let dt = Time::from_si(p.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let mut t = 0.0;
    while t < 3.0 {
        p.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
    }
    let limit = p.max_stable_dt(Time::from_si(t)).to_si();

    // At 95% of the limit the temperature field stays inside the range its sources span.
    let mut safe = p.clone();
    for _ in 0..200 {
        safe.step(Time::from_si(t), Time::from_si(limit * 0.95), &mut bus)
            .expect("under the limit");
    }
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let v = safe.temperature_at(i, j, k).to_si();
                lo = lo.min(v);
                hi = hi.max(v);
            }
        }
    }
    println!(
        "  at 0.95 of the limit the bed stays in [{:.2}, {:.2}] C",
        lo - 273.15,
        hi - 273.15
    );
    assert!(
        lo > 273.15 + 92.0 && hi < 273.15 + 94.0,
        "an isothermal bed with an isothermal inlet must stay isothermal: [{:.3}, {:.3}] C",
        lo - 273.15,
        hi - 273.15
    );
    // The limit is finite and of a sane size for these cells: sub-millisecond, because a
    // millimetre cell with a millimetre-per-second pore velocity is not what sets it — conduction
    // through the steel wall is.
    println!("  the limit is {limit:.3e} s");
    assert!(
        (1e-6..1e-1).contains(&limit),
        "the reported limit is implausible: {limit:.3e} s"
    );
}

/// **Every mutator leaves the flow solved, and the closed form says by how much.**
///
/// The constructor solves, because a quasi-static field read before its solve is a field of zeros
/// wearing the shape of an answer. `repack` did not, and the failure that produced is worse than
/// zeros: the field left behind is *the previous answer* — smooth, bounded, the right order of
/// magnitude, and wrong by 42%.
///
/// It was invisible to the rest of this file because every other test steps the puck before
/// reading anything from it, and a step re-solves. It showed up the first time something asked a
/// freshly repacked basket what its flow rate was.
///
/// The prediction is exact rather than a bound. Every column of cells from inlet to outlet is an
/// independent series chain, so the basket is columns in **parallel** and its conductance is their
/// sum; widening the ring changes only the ring's term:
///
/// ```text
///   Q'/Q  =  (1 − f) + f · m(0.60)/m(0.45),    m(e) = e³/(1−e)²
/// ```
#[test]
fn a_mutator_leaves_the_flow_solved() {
    let (nx, ny, nz) = (16usize, 8usize, 16usize);
    let before = bed_at(2e-3, (nx, ny, nz), 0.45, 9.0e5);
    let mut after = bed_at(2e-3, (nx, ny, nz), 0.45, 9.0e5);

    let mut ring = vec![false; nx * ny * nz];
    let (mut packed, mut edge) = (0usize, 0usize);
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                if !after.is_packed(i, j, k) {
                    continue;
                }
                let e = i == 0
                    || k == 0
                    || i + 1 == nx
                    || k + 1 == nz
                    || !after.is_packed(i - 1, j, k)
                    || !after.is_packed(i + 1, j, k)
                    || !after.is_packed(i, j, k - 1)
                    || !after.is_packed(i, j, k + 1);
                ring[i + nx * (j + ny * k)] = e;
                if j == 0 {
                    packed += 1;
                    if e {
                        edge += 1;
                    }
                }
            }
        }
    }
    after.repack(0.60, |i, j, k| ring[i + nx * (j + ny * k)]);

    // Read the flow **without stepping**, which is the whole point.
    let f = edge as f64 / packed as f64;
    let mobility = |e: f64| e.powi(3) / (1.0 - e).powi(2);
    let predicted = (1.0 - f) + f * mobility(0.60) / mobility(0.45);
    let measured = after.flow_rate().to_si() / before.flow_rate().to_si();

    println!("  ring is {:.1}% of the section", f * 100.0);
    println!("  flow ratio measured {measured:.9}, columns in parallel give {predicted:.9}");
    assert!(
        after.converged(),
        "a repacked basket must be solved, not merely marked stale: residual {:.3e}",
        after.residual()
    );
    assert!(
        (measured / predicted - 1.0).abs() < 1e-9,
        "repack must re-solve: {measured:.6} against {predicted:.6}. A stale field would give \
         exactly 1.000000, which is what this test was written for."
    );
    // And the stale answer is far enough away that this test would have caught it.
    assert!(
        (predicted - 1.0).abs() > 0.2,
        "the two hypotheses have to be distinguishable: {predicted:.4} against a stale 1.0"
    );

    // The same for the other two mutators that change the flow problem.
    let mut hotter = bed_at(2e-3, (nx, ny, nz), 0.45, 9.0e5);
    let cold = hotter.flow_rate().to_si();
    hotter.set_temperature(Temperature::celsius(75.0));
    assert!(hotter.converged(), "set_temperature must leave it solved");
    let warm = hotter.flow_rate().to_si();
    // Water at 75 C is a quarter more viscous than at 93, and the flow is inversely proportional.
    let mu_ratio = Liquid::water()
        .viscosity(Temperature::celsius(75.0))
        .to_si()
        / Liquid::water()
            .viscosity(Temperature::celsius(93.0))
            .to_si();
    println!(
        "  cooling to 75 C: flow {:.6}x, 1/mu gives {:.6}x",
        warm / cold,
        1.0 / mu_ratio
    );
    assert!(
        (warm / cold * mu_ratio - 1.0).abs() < 1e-9,
        "Q goes as 1/mu: {:.6} against {:.6}",
        warm / cold,
        1.0 / mu_ratio
    );

    let mut harder = bed_at(2e-3, (nx, ny, nz), 0.45, 9.0e5);
    let nine = harder.flow_rate().to_si();
    harder.set_drive(Pressure::from_si(12.0e5));
    assert!(harder.converged(), "set_drive must leave it solved");
    assert!(
        (harder.flow_rate().to_si() / nine / (12.0 / 9.0) - 1.0).abs() < 1e-9,
        "Q goes as dp: {:.6}x for 12 bar against 9",
        harder.flow_rate().to_si() / nine
    );
}

/// **The five fields are all populated, and they are not each other.**
///
/// The silent failure this guards against: a field accessor that compiles, renders and is empty —
/// or worse, five accessors that all return the temperature. A viewer would show five identical
/// panels and nothing would say so.
#[test]
fn every_observable_returns_its_own_field() {
    use pantometry_core::ScalarField;
    use pantometry_units::LengthVec;

    let mut p = bed((12, 20, 12), 0.45, 9.0e5);
    p.set_wall_temperature(Temperature::celsius(40.0));
    let dt = Time::from_si(p.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let mut bus = Exchange::new();
    let mut t = 0.0;
    while t < 8.0 {
        p.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
    }

    let all = [
        Observable::Temperature,
        Observable::Pressure,
        Observable::Speed,
        Observable::Extraction,
        Observable::Concentration,
    ];
    let (nx, ny, nz) = p.counts();
    let now = Time::from_si(t);
    let mut sampled = Vec::new();
    for what in all {
        let f = p.field(what);
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        let mut sum = 0.0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let v = f.at(
                        LengthVec::from_si(glam::DVec3::new(
                            (i as f64 + 0.5) * DX,
                            (j as f64 + 0.5) * DX,
                            (k as f64 + 0.5) * DX,
                        )),
                        now,
                    );
                    assert!(v.is_finite(), "{what:?} produced a non-finite value");
                    lo = lo.min(v);
                    hi = hi.max(v);
                    sum += v;
                }
            }
        }
        println!("  {what:<14?} [{lo:.6e}, {hi:.6e}]  unit {:?}", f.unit());
        assert!(
            hi > lo,
            "{what:?} is uniform, which means it is not being computed"
        );
        sampled.push(sum);
    }
    // No two fields are the same field.
    for a in 0..sampled.len() {
        for b in a + 1..sampled.len() {
            assert!(
                (sampled[a] - sampled[b]).abs() > 1e-12 * sampled[a].abs().max(1.0),
                "{:?} and {:?} produced the same field",
                all[a],
                all[b]
            );
        }
    }
}
