//! `ThermalNetwork` against closed forms, and against the failures conservation cannot see.
//!
//! A link contributes `+q` to one node and `−q` to another in the same sum, so the ledger is
//! blind to links by construction: a sign error, a transposed index or a dropped link passes the
//! conservation audit at machine precision. Every check here is therefore per node or against a
//! formula computed in this file — never on the total.

use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::{Environment, Node, ThermalNetwork, HEAT};
use pantometry_units::{Area, Conductance, Length, Temperature, Time, Volume};

/// Aluminium with the radiation switched off, so environment loss is exactly `hA·ΔT` and the
/// closed forms below are the closed forms of the problem rather than of a linearisation.
fn grey() -> Substance {
    let mut s = Substance::aluminium_6061();
    if let Some(t) = s.thermal.as_mut() {
        t.emissivity = 0.0;
    }
    s
}

fn run(net: &mut ThermalNetwork, dt: f64, steps: usize) {
    let mut bus = Exchange::new();
    for k in 0..steps {
        net.step(Time::s(k as f64 * dt), Time::s(dt), &mut bus)
            .unwrap();
    }
}

/// **A network of one node is a `LumpedMass`, bit for bit, over a whole trajectory.**
///
/// The reduction that makes the new domain trustworthy: `LumpedMass` is checked against
/// exponential decay, radiative equilibrium and a measured time constant, and this inherits all
/// of it at zero cost — but only if the agreement is *exact*. "Close" would mean the two have
/// separate arithmetic that happens to be near, which is the state the shared
/// `linearised_loss_conductance` exists to prevent.
///
/// Run with radiation **on**, so the nonlinear `esA(T^4 - Ta^4)` is on trial too, and with heat
/// arriving on the bus, so the source path is compared as well as the loss path.
#[test]
fn one_node_is_a_lumped_mass_bit_for_bit() {
    let substance = Substance::aluminium_6061();
    let volume = Volume::from_si(60e-3 * 60e-3 * 3e-3);
    let thickness = Length::mm(1.5);
    let start = Temperature::celsius(85.0);
    let env = || {
        Environment::still_air(
            Temperature::celsius(20.0),
            Area::from_si(2.0 * 60e-3 * 60e-3),
        )
    };

    let mut lump = pantometry_thermal::LumpedMass::new(
        "plate",
        substance.clone(),
        volume,
        thickness,
        start,
        env(),
    );
    let mut net = ThermalNetwork::new("plate");
    let only = net.node_losing_to("plate", substance, volume, thickness, start, env());
    net.absorbing(only).unwrap();

    assert_eq!(
        net.max_stable_dt(Time::s(0.0)).to_si().to_bits(),
        lump.max_stable_dt(Time::s(0.0)).to_si().to_bits(),
        "the step limits must agree exactly, or the two are computing different conductances"
    );

    let dt = Time::s(0.25);
    for k in 0..4_000 {
        // Time-varying, so the comparison is not one repeated arithmetic coincidence. Heat in for
        // the first half, then it coasts down.
        let watts = if k < 2_000 { 6.0 } else { 0.0 };

        let mut bus_a = Exchange::new();
        bus_a.publish(HEAT, watts * dt.to_si());
        lump.step(Time::s(k as f64 * 0.25), dt, &mut bus_a).unwrap();

        let mut bus_b = Exchange::new();
        bus_b.publish(HEAT, watts * dt.to_si());
        net.step(Time::s(k as f64 * 0.25), dt, &mut bus_b).unwrap();

        assert_eq!(
            net.temperature(only).to_si().to_bits(),
            lump.temperature().to_si().to_bits(),
            "diverged at step {k}: {} against {}",
            net.temperature(only).to_si(),
            lump.temperature().to_si()
        );
    }

    // It has to have gone somewhere, or this compared a constant with a constant.
    let rise = net.temperature(only).to_si() - start.to_si();
    assert!(
        rise.abs() > 5.0,
        "the trajectory should have moved, got {rise:.3} K"
    );
    assert_eq!(
        net.lost_energy().to_si().to_bits(),
        lump.lost_energy().to_si().to_bits()
    );
    assert_eq!(
        net.absorbed_energy().to_si().to_bits(),
        lump.absorbed_energy().to_si().to_bits()
    );
}

/// **Two nodes and nothing else: one exponential, and its rate is a pure function of the link.**
///
/// With no environments and no source, `C₁T₁' = K(T₂−T₁)` and `C₂T₂' = K(T₁−T₂)`, so the
/// difference obeys `D' = −K·D·(C₁+C₂)/(C₁C₂)` exactly:
///
/// ```text
///   τ = C₁C₂ / ( K (C₁ + C₂) )
/// ```
///
/// The cleanest statement of what a conductance *is*, and it isolates the link completely — no
/// environment, no radiation, no nonlinearity. Checked first against the **discrete** form
/// `(1 − h/τ)ⁿ`, which pins that no stray factor hides in the capacities, and then against the
/// continuous one by its *rate*: the Euler shortfall must fall by ten when the step does.
/// Measured 10.011; the band is wide enough that a scheme gone half-order still fails it.
#[test]
fn two_nodes_relax_toward_each_other_at_the_rate_the_conductance_sets() {
    let volume = Volume::from_si(1e-4);
    let build = || {
        let mut net = ThermalNetwork::new("pair");
        let a = net.node(
            "a",
            grey(),
            volume,
            Length::mm(5.0),
            Temperature::celsius(80.0),
        );
        let b = net.node(
            "b",
            grey(),
            Volume::from_si(3e-4),
            Length::mm(5.0),
            Temperature::celsius(20.0),
        );
        net.link(a, b, Conductance::w_per_k(0.8)).unwrap();
        (net, a, b)
    };

    let c1 = grey().heat_capacity(volume).unwrap().to_si();
    let c2 = grey().heat_capacity(Volume::from_si(3e-4)).unwrap().to_si();
    let tau = c1 * c2 / (0.8 * (c1 + c2));

    // The discrete form, to nine figures: this is the scheme's own solution, so anything but
    // agreement means a factor is wrong somewhere in the capacities or the flux.
    let (mut net, a, b) = build();
    let h = 0.5;
    let steps = 200;
    run(&mut net, h, steps);
    let got = net.temperature(a).to_si() - net.temperature(b).to_si();
    let discrete = 60.0 * (1.0 - h / tau).powi(steps as i32);
    assert!(
        (got / discrete - 1.0).abs() < 1e-9,
        "discrete: got {got:.9}, expected {discrete:.9}"
    );

    // And the rate against the continuum. Explicit Euler is first order, so refining the step
    // by ten must cut the shortfall by about ten — the claim with content, since a wrong `tau`
    // would show as a constant offset that refinement does not touch.
    let err = |h: f64| {
        let (mut net, a, b) = build();
        let steps = (100.0 / h) as usize;
        run(&mut net, h, steps);
        let got = net.temperature(a).to_si() - net.temperature(b).to_si();
        (got - 60.0 * (-100.0 / tau).exp()).abs()
    };
    let (coarse, fine) = (err(0.5), err(0.05));
    let ratio = coarse / fine;
    assert!(
        (7.0..13.0).contains(&ratio),
        "first order: {coarse:.6e} at h=0.5 against {fine:.6e} at h=0.05, ratio {ratio:.2}"
    );
}

/// **A series ladder settles where resistances in series say.**
///
/// Chain `1 –K– 2 –K– 3`, environment on the last node only, power into the first:
///
/// ```text
///   T₁ − Tₐ = P · ( 1/K₁ + 1/K₂ + 1/G₃ ),   Tᵢ − Tᵢ₊₁ = P / Kᵢ
/// ```
///
/// Textbook, independent of this code, and it is exactly the number the report that asked for
/// this domain could not get: the drop **across the joint** between what fails and what you can
/// measure. A single lumped mass reports the skin temperature and calls it the winding.
#[test]
fn a_series_ladder_settles_where_resistances_in_series_say() {
    let area = Area::from_si(0.03);
    let mut net = ThermalNetwork::new("ladder");
    let n1 = net.node(
        "winding",
        grey(),
        Volume::from_si(1e-5),
        Length::mm(3.0),
        Temperature::celsius(20.0),
    );
    let n2 = net.node(
        "stator",
        grey(),
        Volume::from_si(5e-5),
        Length::mm(3.0),
        Temperature::celsius(20.0),
    );
    let n3 = net.node_losing_to(
        "housing",
        grey(),
        Volume::from_si(2e-4),
        Length::mm(3.0),
        Temperature::celsius(20.0),
        Environment::still_air(Temperature::celsius(20.0), area),
    );
    net.link(n1, n2, Conductance::w_per_k(1.5)).unwrap();
    net.link(n2, n3, Conductance::w_per_k(0.8)).unwrap();
    net.absorbing(n1).unwrap();

    let power = 4.0;
    let dt = 0.05;
    let mut bus = Exchange::new();
    for k in 0..2_000_000 {
        bus.publish(HEAT, power * dt);
        net.step(Time::s(k as f64 * dt), Time::s(dt), &mut bus)
            .unwrap();
    }

    // Computed here, from the conductances this test set.
    let g3 = 7.0 * area.to_si(); // still_air is h = 7 W/m^2/K, and emissivity is zero
    let want_1 = power * (1.0 / 1.5 + 1.0 / 0.8 + 1.0 / g3);
    let rise_1 = net.temperature(n1).to_si() - Temperature::celsius(20.0).to_si();
    assert!(
        (rise_1 / want_1 - 1.0).abs() < 1e-6,
        "winding: {rise_1:.6} K against {want_1:.6} K"
    );

    // And the drops across each joint, which is the thing a single lumped mass cannot report.
    let d12 = net.temperature(n1).to_si() - net.temperature(n2).to_si();
    let d23 = net.temperature(n2).to_si() - net.temperature(n3).to_si();
    assert!(
        (d12 / (power / 1.5) - 1.0).abs() < 1e-6,
        "across the first joint: {d12:.4} K"
    );
    assert!(
        (d23 / (power / 0.8) - 1.0).abs() < 1e-6,
        "across the second: {d23:.4} K"
    );

    // The steady flux through every joint is the whole input, since nothing leaves in between.
    assert!((net.heat_flow(n1, n2).to_si() / power - 1.0).abs() < 1e-6);
}

/// **Reciprocity: inject at one node and read at another, then swap. The rises are equal.**
///
/// The steady-state conductance matrix is symmetric, so its inverse is. This is the single most
/// valuable test here, because it catches exactly the class conservation cannot: a link applied
/// asymmetrically — flux added to one node and not subtracted from the other, or a transposed
/// index — keeps every total right and breaks this identity.
#[test]
fn injecting_at_a_and_reading_at_b_equals_the_other_way_round() {
    let settle = |into: usize| -> f64 {
        let mut net = ThermalNetwork::new("web");
        let a = net.node(
            "a",
            grey(),
            Volume::from_si(1e-5),
            Length::mm(2.0),
            Temperature::celsius(20.0),
        );
        let b = net.node(
            "b",
            grey(),
            Volume::from_si(7e-5),
            Length::mm(2.0),
            Temperature::celsius(20.0),
        );
        let c = net.node_losing_to(
            "c",
            grey(),
            Volume::from_si(3e-5),
            Length::mm(2.0),
            Temperature::celsius(20.0),
            Environment::still_air(Temperature::celsius(20.0), Area::from_si(0.02)),
        );
        // A triangle, so the test also covers a cycle — where a careless per-link loop over a
        // shared buffer double-counts.
        net.link(a, b, Conductance::w_per_k(0.9)).unwrap();
        net.link(b, c, Conductance::w_per_k(0.4)).unwrap();
        net.link(a, c, Conductance::w_per_k(0.25)).unwrap();

        let nodes = [a, b, c];
        net.absorbing(nodes[into]).unwrap();
        let read = nodes[1 - into.min(1)];
        let (dt, mut bus) = (0.02, Exchange::new());
        for k in 0..1_500_000 {
            bus.publish(HEAT, 3.0 * dt);
            net.step(Time::s(k as f64 * dt), Time::s(dt), &mut bus)
                .unwrap();
        }
        net.temperature(read).to_si() - Temperature::celsius(20.0).to_si()
    };

    let (a_to_b, b_to_a) = (settle(0), settle(1));
    assert!(
        a_to_b > 1.0,
        "the test needs a rise to compare, got {a_to_b:.4} K"
    );
    assert!(
        (a_to_b / b_to_a - 1.0).abs() < 1e-6,
        "reciprocity: {a_to_b:.6} K against {b_to_a:.6} K"
    );
}

/// **A network already at one temperature does not move, bit for bit.**
///
/// Every link flux and every loss is identically zero, so the state must be unchanged after any
/// number of steps — not close, *equal*. Catches a sign error or a transposed index that a
/// total-only check cannot see, because those keep the sum right while moving individual nodes.
#[test]
fn a_uniform_network_is_a_fixed_point_exactly() {
    let ambient = Temperature::celsius(25.0);
    let mut net = ThermalNetwork::new("still");
    let a = net.node("a", grey(), Volume::from_si(1e-5), Length::mm(2.0), ambient);
    let b = net.node_losing_to(
        "b",
        grey(),
        Volume::from_si(4e-5),
        Length::mm(2.0),
        ambient,
        Environment::still_air(ambient, Area::from_si(0.01)),
    );
    let c = net.node("c", grey(), Volume::from_si(2e-5), Length::mm(2.0), ambient);
    net.link(a, b, Conductance::w_per_k(1.2)).unwrap();
    net.link(b, c, Conductance::w_per_k(0.7)).unwrap();
    net.link(c, a, Conductance::w_per_k(0.3)).unwrap();

    let before: Vec<u64> = [a, b, c]
        .iter()
        .map(|n| net.temperature(*n).to_si().to_bits())
        .collect();
    run(&mut net, 0.1, 5_000);
    let after: Vec<u64> = [a, b, c]
        .iter()
        .map(|n| net.temperature(*n).to_si().to_bits())
        .collect();
    assert_eq!(before, after, "a uniform network moved");
    assert_eq!(net.lost_energy().to_si(), 0.0);
}

/// **A closed network conserves, and the ledger says so against a scale that means something.**
///
/// No environments and nothing on the bus, so `Σ Cᵢ(Tᵢ − Tᵢ(0))` must stay at zero — measured
/// against the largest single node's holding rather than against the net, which is what
/// `Ledger`'s scale exists for. Includes a cycle, where a per-link loop over a shared buffer
/// would double-count.
#[test]
fn a_closed_network_conserves_against_a_scale_that_is_not_zero() {
    let mut net = ThermalNetwork::new("closed");
    let a = net.node(
        "a",
        grey(),
        Volume::from_si(1e-5),
        Length::mm(2.0),
        Temperature::celsius(90.0),
    );
    let b = net.node(
        "b",
        grey(),
        Volume::from_si(4e-5),
        Length::mm(2.0),
        Temperature::celsius(20.0),
    );
    let c = net.node(
        "c",
        grey(),
        Volume::from_si(2e-5),
        Length::mm(2.0),
        Temperature::celsius(50.0),
    );
    net.link(a, b, Conductance::w_per_k(1.2)).unwrap();
    net.link(b, c, Conductance::w_per_k(0.7)).unwrap();
    net.link(c, a, Conductance::w_per_k(0.3)).unwrap();

    let holding = |net: &ThermalNetwork, n: Node| {
        let v = if n == a {
            1e-5
        } else if n == b {
            4e-5
        } else {
            2e-5
        };
        grey().heat_capacity(Volume::from_si(v)).unwrap().to_si() * net.rise(n).to_si()
    };

    let steps = 200_000;
    run(&mut net, 0.05, steps);
    let net_total = holding(&net, a) + holding(&net, b) + holding(&net, c);
    let scale = [a, b, c]
        .iter()
        .map(|n| holding(&net, *n).abs())
        .fold(0.0f64, f64::max);
    assert!(
        scale > 100.0,
        "the test needs heat to have moved, scale {scale:.3} J"
    );
    // The tolerance is the accumulation bound, not a number chosen after looking: each step
    // rounds the temperatures at `eps` relative, and the worst case is that every one of the
    // `steps` roundings lands the same way. `2e5 * 2.2e-16 = 4.4e-11`. The measured drift is
    // `1.0e-12`, about 2% of the bound and roughly ten times the `sqrt(N)*eps` a random walk
    // would give — the right order for float accumulation and far below anything a dropped or
    // one-sided link would produce. It scales with the run, so lengthening the test does not
    // silently make it stricter than the arithmetic allows.
    let bound = steps as f64 * f64::EPSILON;
    assert!(
        net_total.abs() / scale < bound,
        "closed network drifted by {net_total:.6e} J against a scale of {scale:.3} J, \
         relative {:.3e} against an accumulation bound of {bound:.3e}",
        net_total.abs() / scale
    );

    // And it really did level out, so the conservation above is not the conservation of nothing.
    let spread = net.temperature(a).to_si() - net.temperature(b).to_si();
    assert!(
        spread.abs() < 1e-6,
        "the network should have levelled, spread {spread:.6} K"
    );
}

/// The mistakes a caller makes are refused, and named.
#[test]
fn the_mistakes_are_refused_by_name() {
    let mut net = ThermalNetwork::new("net");
    let a = net.node(
        "a",
        grey(),
        Volume::from_si(1e-5),
        Length::mm(2.0),
        Temperature::celsius(20.0),
    );
    let b = net.node(
        "b",
        grey(),
        Volume::from_si(1e-5),
        Length::mm(2.0),
        Temperature::celsius(20.0),
    );

    let self_link = net.link(a, a, Conductance::w_per_k(1.0)).unwrap_err();
    assert!(
        self_link.quantity.contains("itself"),
        "{}",
        self_link.quantity
    );

    let negative = net.link(a, b, Conductance::w_per_k(-1.0)).unwrap_err();
    assert!(
        negative.quantity.contains("not negative"),
        "{}",
        negative.quantity
    );

    // A handle from a different network is refused rather than addressing whatever sits at that
    // index — the one hole a plain index could not close.
    let mut other = ThermalNetwork::new("other");
    let elsewhere = other.node(
        "a",
        grey(),
        Volume::from_si(1e-5),
        Length::mm(2.0),
        Temperature::celsius(20.0),
    );
    let foreign = net
        .link(a, elsewhere, Conductance::w_per_k(1.0))
        .unwrap_err();
    assert!(
        foreign.quantity.contains("different network"),
        "{}",
        foreign.quantity
    );

    // Parallel conductances add, following the bus's convention for repeated offers.
    net.link(a, b, Conductance::w_per_k(0.5)).unwrap();
    net.link(a, b, Conductance::w_per_k(0.5)).unwrap();
    let (t_a, t_b) = (net.temperature(a).to_si(), net.temperature(b).to_si());
    assert_eq!(t_a, t_b); // still uniform, so the flow is zero and only the count is on trial
    assert_eq!(net.nodes(), 2);
    assert_eq!(net.node_named("b"), Some(b));
    assert_eq!(net.node_named("nope"), None);
}
