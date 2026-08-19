//! `ThermalNetwork::steady_state` against closed forms, and against the marching it replaces.
//!
//! The solver and the step loop are two independent routes to the same balance, so agreement
//! between them is worth something — but only in one direction. Marching is the *reference*
//! here: it is the thing already checked against series-resistance formulae, exponential decay
//! and a measured convergence rate. So every test below either compares against a formula
//! computed in this file, or against a march that has been run long enough to have arrived.

use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::{Environment, ThermalNetwork, HEAT};
use pantometry_units::{Area, Conductance, Length, Power, Temperature, Time, Volume};

/// Aluminium with radiation switched off, so the balance is linear and the closed forms below
/// are the closed forms of the problem rather than of a linearisation.
fn grey() -> Substance {
    let mut s = Substance::aluminium_6061();
    if let Some(t) = s.thermal.as_mut() {
        t.emissivity = 0.0;
    }
    s
}

fn ambient() -> Temperature {
    Temperature::celsius(20.0)
}

/// The three-node ladder from the closed-form suite: winding, stator, housing.
fn ladder(substance: fn() -> Substance) -> (ThermalNetwork, [pantometry_thermal::Node; 3], f64) {
    let area = Area::from_si(0.03);
    let mut net = ThermalNetwork::new("ladder");
    let n1 = net.node(
        "winding",
        substance(),
        Volume::from_si(1e-5),
        Length::mm(3.0),
        ambient(),
    );
    let n2 = net.node(
        "stator",
        substance(),
        Volume::from_si(5e-5),
        Length::mm(3.0),
        ambient(),
    );
    let n3 = net.node_losing_to(
        "housing",
        substance(),
        Volume::from_si(2e-4),
        Length::mm(3.0),
        ambient(),
        Environment::still_air(ambient(), area),
    );
    net.link(n1, n2, Conductance::w_per_k(1.5)).unwrap();
    net.link(n2, n3, Conductance::w_per_k(0.8)).unwrap();
    net.absorbing(n1).unwrap();
    // still_air is h = 7 W/m^2/K; with emissivity zero that is the whole conductance to ambient.
    (net, [n1, n2, n3], 7.0 * area.to_si())
}

/// **Resistances in series, to machine precision rather than to a march's tolerance.**
///
/// `T₁ − Tₐ = P(1/K₁ + 1/K₂ + 1/G₃)`, and each joint carries the whole power because nothing
/// leaves in between. This is a direct solve of a linear system, so "close" is not the standard:
/// the residual is driven to zero and the answer should agree with the formula at 1e-12, not at
/// the 1e-6 the marched version of this test settles for after two million steps.
#[test]
fn a_series_ladder_solves_to_the_formula_at_machine_precision() {
    let (net, [n1, n2, n3], g3) = ladder(grey);
    let power = 4.0;
    let settled = net
        .steady_state(Power::w(power))
        .expect("it has a solution");

    let rise = |n| settled.temperature(n).to_si() - ambient().to_si();
    let want_1 = power * (1.0 / 1.5 + 1.0 / 0.8 + 1.0 / g3);
    assert!(
        (rise(n1) / want_1 - 1.0).abs() < 1e-12,
        "winding: {:.9} K against {want_1:.9} K",
        rise(n1)
    );
    assert!(
        ((rise(n1) - rise(n2)) / (power / 1.5) - 1.0).abs() < 1e-12,
        "first joint: {:.9} K",
        rise(n1) - rise(n2)
    );
    assert!(
        ((rise(n2) - rise(n3)) / (power / 0.8) - 1.0).abs() < 1e-12,
        "second joint: {:.9} K",
        rise(n2) - rise(n3)
    );
    assert!(
        (rise(n3) / (power / g3) - 1.0).abs() < 1e-12,
        "housing: {:.9} K",
        rise(n3)
    );
}

/// **The solver and the step loop agree, and marching is the one being trusted.**
///
/// Two independent routes to the same balance: one solves it, one converges to it. The step
/// loop is the reference — it is what the closed-form suite already checks — so this asks
/// whether the solver lands where marching arrives, run out to thirty time constants where the
/// remaining transient is `e⁻³⁰ ≈ 1e-13`.
///
/// Run with **radiation on**, which is the case that makes this test worth having: the balance
/// is `T⁴` and a single linear solve would answer a different question. If Newton were dropped
/// for one solve from ambient the two would disagree in the third figure.
#[test]
fn the_solver_lands_where_marching_arrives_with_radiation_on() {
    let power = 4.0;
    let (net, nodes, _) = ladder(Substance::aluminium_6061);
    let solved: Vec<f64> = {
        let s = net
            .steady_state(Power::w(power))
            .expect("it has a solution");
        nodes.iter().map(|n| s.temperature(*n).to_si()).collect()
    };

    let (mut net, nodes, _) = ladder(Substance::aluminium_6061);
    let (dt, mut bus) = (0.05, Exchange::new());
    for k in 0..2_000_000 {
        bus.publish(HEAT, power * dt);
        net.step(Time::s(k as f64 * dt), Time::s(dt), &mut bus)
            .unwrap();
    }
    let marched: Vec<f64> = nodes.iter().map(|n| net.temperature(*n).to_si()).collect();

    for ((s, m), n) in solved.iter().zip(&marched).zip(nodes.iter()) {
        assert!(
            (s - m).abs() < 1e-6,
            "{}: solved {s:.9} K against marched {m:.9} K",
            net.label(*n)
        );
    }

    // And radiation is actually doing something here, or this compared the linear case twice:
    // with emissivity 0.09 the housing sheds partly by radiation, so it settles cooler than the
    // convection-only answer the grey ladder gives.
    let (grey_net, grey_nodes, _) = ladder(grey);
    let grey_housing = grey_net
        .steady_state(Power::w(power))
        .unwrap()
        .temperature(grey_nodes[2])
        .to_si();
    assert!(
        grey_housing - marched[2] > 1.0,
        "radiation should cool the housing meaningfully: grey {grey_housing:.3} K against \
         radiating {:.3} K",
        marched[2]
    );
}

/// **One node reduces to `LumpedMass::equilibrium_rise`.**
///
/// The reduction that ties this to physics already checked. `equilibrium_rise` Newton-solves
/// `P = hA·ΔT + εσA((Tₐ+ΔT)⁴ − Tₐ⁴)` on a single body, and a network of one node is that body,
/// so the two must agree — with radiation on, which is the only case where the claim has
/// content.
#[test]
fn one_node_agrees_with_the_lumped_equilibrium_rise() {
    let volume = Volume::from_si(60e-3 * 60e-3 * 3e-3);
    let area = Area::from_si(2.0 * 60e-3 * 60e-3);
    let power = 12.0;

    let lump = pantometry_thermal::LumpedMass::new(
        "plate",
        Substance::aluminium_6061(),
        volume,
        Length::mm(1.5),
        ambient(),
        Environment::still_air(ambient(), area),
    );
    let want = lump.equilibrium_rise(Power::w(power)).to_si();

    let mut net = ThermalNetwork::new("plate");
    let only = net.node_losing_to(
        "plate",
        Substance::aluminium_6061(),
        volume,
        Length::mm(1.5),
        ambient(),
        Environment::still_air(ambient(), area),
    );
    net.absorbing(only).unwrap();
    let got = net
        .steady_state(Power::w(power))
        .expect("one node with an environment has a solution")
        .temperature(only)
        .to_si()
        - ambient().to_si();

    assert!(
        (got / want - 1.0).abs() < 1e-9,
        "network {got:.9} K against lumped {want:.9} K"
    );
    // Radiation matters at this power, so the agreement is not two linear solves agreeing: the
    // convection-only answer is far away.
    let linear = power / (7.0 * area.to_si());
    assert!(
        linear - want > 5.0,
        "the radiative term should be worth several kelvin here: linear {linear:.3} K against \
         {want:.3} K"
    );
}

/// **A network with nowhere for heat to go is refused, not answered.**
///
/// It has no steady state — it warms without limit — and marching it is perfectly well defined,
/// so the temptation is to return the linear solve's answer, which would be a plausible finite
/// number for a temperature that does not exist. The singular matrix is the physics saying so.
#[test]
fn a_network_with_no_environment_is_refused_rather_than_answered() {
    let mut net = ThermalNetwork::new("sealed");
    let a = net.node(
        "a",
        grey(),
        Volume::from_si(1e-5),
        Length::mm(2.0),
        ambient(),
    );
    let b = net.node(
        "b",
        grey(),
        Volume::from_si(4e-5),
        Length::mm(2.0),
        ambient(),
    );
    net.link(a, b, Conductance::w_per_k(1.0)).unwrap();
    net.absorbing(a).unwrap();

    let refused = net.steady_state(Power::w(3.0)).unwrap_err();
    assert!(
        refused.quantity.contains("nowhere to go"),
        "{}",
        refused.quantity
    );

    // With no power there *is* a solution — everything at its starting temperature — and it is
    // returned rather than refused, because the refusal is about the physics and not about the
    // shape of the network.
    let idle = net
        .steady_state(Power::w(0.0))
        .expect("no power, no problem");
    assert_eq!(idle.temperature(a).to_si(), idle.temperature(b).to_si());

    // A node with no path to *any* environment is the same failure with a subtler cause: the
    // network as a whole has one, this component does not.
    let mut split = ThermalNetwork::new("split");
    let joined = split.node_losing_to(
        "joined",
        grey(),
        Volume::from_si(1e-5),
        Length::mm(2.0),
        ambient(),
        Environment::still_air(ambient(), Area::from_si(0.01)),
    );
    let orphan = split.node(
        "orphan",
        grey(),
        Volume::from_si(1e-5),
        Length::mm(2.0),
        ambient(),
    );
    split.absorbing(joined).unwrap();
    let err = split.steady_state(Power::w(1.0)).unwrap_err();
    assert!(err.quantity.contains("singular"), "{}", err.quantity);
    let _ = orphan;
}

/// **The solve does not disturb the network.**
///
/// It is a question, not a step. Asking it must leave every temperature and both running totals
/// exactly where they were, or a caller who asks mid-run has silently changed their simulation.
#[test]
fn asking_where_it_settles_does_not_move_anything() {
    let (mut net, nodes, _) = ladder(grey);
    let (dt, mut bus) = (0.05, Exchange::new());
    for k in 0..1_000 {
        bus.publish(HEAT, 4.0 * dt);
        net.step(Time::s(k as f64 * dt), Time::s(dt), &mut bus)
            .unwrap();
    }

    let before: Vec<u64> = nodes
        .iter()
        .map(|n| net.temperature(*n).to_si().to_bits())
        .collect();
    let (absorbed, lost) = (net.absorbed_energy().to_si(), net.lost_energy().to_si());

    let settled = net.steady_state(Power::w(4.0)).unwrap();
    assert_eq!(settled.nodes(), 3);

    let after: Vec<u64> = nodes
        .iter()
        .map(|n| net.temperature(*n).to_si().to_bits())
        .collect();
    assert_eq!(before, after, "steady_state moved the network");
    assert_eq!(net.absorbed_energy().to_si(), absorbed);
    assert_eq!(net.lost_energy().to_si(), lost);

    // And it answered about the settled state rather than about the current one, which after a
    // thousand steps is nowhere near it.
    let now = net.temperature(nodes[0]).to_si();
    let end = settled.temperature(nodes[0]).to_si();
    assert!(
        end - now > 5.0,
        "fifty seconds in, the winding is at {now:.3} K and settles at {end:.3} K"
    );
}

/// **A radiation-dominated balance, against a root found by bisection rather than by Newton.**
///
/// The test that was missing. `NEWTON_STEPS` was first set to eight — twice the worst case the
/// tests above exercised — and that refused **a kilowatt**, because none of them loads a node
/// hard enough for the `T⁴` term to dominate. At ambient the radiative slope `4εσAT³` is tiny
/// against what the balance needs, so the first solve overshoots enormously and Newton then
/// walks down at the `3/4` ratio a quartic gives: twelve iterations here, sixty-six at a
/// terawatt.
///
/// Checked against **bisection**, which shares no arithmetic with the solver and cannot fail in
/// the same direction: the balance is monotone in `T`, so a sign change brackets the root and
/// eighty halvings pin it to the last bit.
#[test]
fn a_kilowatt_on_one_node_solves_where_bisection_says() {
    let area = Area::from_si(0.01);
    let volume = Volume::from_si(1e-5);
    let power = 1_000.0;

    let mut net = ThermalNetwork::new("hot");
    let n = net.node_losing_to(
        "n",
        Substance::aluminium_6061(),
        volume,
        Length::mm(1.0),
        ambient(),
        Environment::still_air(ambient(), area),
    );
    net.absorbing(n).unwrap();
    let got = net
        .steady_state(Power::w(power))
        .expect("a kilowatt is not an unreasonable question")
        .temperature(n)
        .to_si();

    // The balance, written out here from the constants rather than read off the library.
    let ta = ambient().to_si();
    let emissivity = Substance::aluminium_6061().thermal.unwrap().emissivity;
    let sigma = 5.670_374_419e-8;
    let residual = |t: f64| {
        power
            - 7.0 * area.to_si() * (t - ta)
            - emissivity * sigma * area.to_si() * (t.powi(4) - ta.powi(4))
    };
    let (mut lo, mut hi) = (ta, 1.0e6);
    assert!(
        residual(lo) > 0.0 && residual(hi) < 0.0,
        "the root is bracketed"
    );
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let want = 0.5 * (lo + hi);

    assert!(
        (got / want - 1.0).abs() < 1e-9,
        "solved {got:.6} K against bisection's {want:.6} K"
    );
    // And radiation really is carrying this, or the many-iteration path was never taken: the
    // convection-only answer is an order of magnitude away.
    let convection_only = ta + power / (7.0 * area.to_si());
    assert!(
        convection_only / want > 5.0,
        "this should be radiation-dominated: convection alone gives {convection_only:.0} K \
         against {want:.0} K"
    );
}

/// **The path conductance is the series conductance, exactly, when nothing radiates.**
///
/// `1/g = 1/K₁ + 1/K₂ + 1/G₃` — resistances in series, the same formula the ladder settles to,
/// so this is checked against arithmetic in this file rather than against another solve. Exact
/// because the balance is linear: the slope `ΔP/ΔT` is the same wherever it is taken, which the
/// test confirms by taking it at two operating points a factor of a hundred apart.
///
/// It exists because a caller was assembling this by hand out of numbers the network already
/// held, to hand to `Winding::runaway_current`. `FRICTION.md` 20.
#[test]
fn the_path_conductance_is_resistances_in_series() {
    let (net, [n1, n2, n3], g3) = ladder(grey);
    let series = 1.0 / (1.0 / 1.5 + 1.0 / 0.8 + 1.0 / g3);

    let g = net.path_conductance(n1, Power::w(4.0)).unwrap().to_si();
    assert!(
        (g / series - 1.0).abs() < 1e-9,
        "{g:.9} W/K against {series:.9} W/K"
    );

    // Linear, so the operating point does not matter.
    let far = net.path_conductance(n1, Power::w(400.0)).unwrap().to_si();
    assert!((far / g - 1.0).abs() < 1e-9, "{far:.9} against {g:.9}");

    // A node further down the path sees less resistance, and exactly how much less.
    let g2 = net.path_conductance(n2, Power::w(4.0)).unwrap().to_si();
    let g3n = net.path_conductance(n3, Power::w(4.0)).unwrap().to_si();
    assert!((g2 - 1.0 / (1.0 / 0.8 + 1.0 / g3)).abs() < 1e-9, "{g2:.9}");
    assert!(
        (g3n - g3).abs() < 1e-9,
        "{g3n:.9} against the environment's {g3:.9}"
    );
    assert!(g < g2 && g2 < g3n, "further from the air is a worse path");

    // With radiation it is the local slope and must be *better* than the grey path, because
    // radiating is a second way out. Sign, not size: the size is the whole point of it being
    // local rather than a constant.
    let (radiating, nodes, _) = ladder(Substance::aluminium_6061);
    let hot = radiating
        .path_conductance(nodes[0], Power::w(40.0))
        .unwrap()
        .to_si();
    assert!(hot > g, "radiating {hot:.6} should beat grey {g:.6}");
}
