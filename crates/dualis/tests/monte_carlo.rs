//! A Monte Carlo study run through the facade, against a closed form and against itself.
//!
//! `Ensemble` claims two things: that many independent samples can be spread over cores, and
//! that doing so cannot change the answer. The unit tests check the second on a bare generator.
//! This checks both on the thing the workspace is actually for — physics — because a sample
//! closure that reaches for shared state is the way the claim gets broken, and only a realistic
//! closure would do that.

use dualis::prelude::*;

/// **Photon shot noise, against the closed form, and identical on any number of threads.**
///
/// A detector collecting `n` photoelectrons on average has variance `n` — that is what makes
/// Poisson statistics the floor a camera cannot get under, and the signal-to-noise ratio `√n`.
/// Both are known exactly, so the Monte Carlo is checked against arithmetic done here rather
/// than against a second run of itself.
#[test]
fn shot_noise_is_root_n_and_does_not_care_how_many_cores_found_it() {
    let mean_signal = 10_000.0;
    let draw = move |_: u64, mut rng: Rng| rng.poisson(mean_signal) as f64;

    let single = Ensemble::new(31, 50_000)
        .estimate(draw)
        .expect("many samples");
    let many = Ensemble::new(31, 50_000)
        .with_threads(12)
        .estimate(draw)
        .expect("many samples");

    // Bit-identical, not merely statistically indistinguishable. A shared generator would give
    // two answers that both look plausible and differ, which is the failure that never gets
    // caught because it is indistinguishable from noise by inspection.
    assert_eq!(single.mean.to_bits(), many.mean.to_bits());
    assert_eq!(
        single.standard_error.to_bits(),
        many.standard_error.to_bits()
    );

    // The mean is the rate, within the error it reports on itself.
    assert!(
        many.within(4.0, mean_signal),
        "{mean_signal} is {:.2} standard errors from {:.2}",
        (many.mean - mean_signal).abs() / many.standard_error,
        many.mean
    );

    // And the *spread* is √n, which is the physics rather than the estimator. Computed here.
    let want = mean_signal.sqrt();
    let got = many.standard_deviation();
    assert!(
        (got / want - 1.0).abs() < 0.02,
        "shot noise {got:.2} against √n = {want:.2}"
    );

    // The signal-to-noise ratio a datasheet quotes, from the same two numbers.
    assert!(
        ((many.mean / got) / want - 1.0).abs() < 0.03,
        "SNR {:.1} against √n = {want:.1}",
        many.mean / got
    );
}

/// **A parameter sweep: many simulations, not many draws.**
///
/// The other Monte Carlo shape, and the more expensive one — each sample here *runs a
/// simulation to a settled state* rather than drawing a number. A winding whose length varies
/// by a manufacturing tolerance, asked how hot it gets.
///
/// The claim is the same and matters more: twelve threads each stepping their own `Simulation`
/// must give the same answers as one thread doing them in turn. `Domain` is not `Send`, so each
/// sample builds its own and nothing crosses a thread boundary — which the type system enforces
/// here rather than leaving to discipline.
#[test]
fn a_sweep_of_whole_simulations_is_reproducible_across_threads() {
    let sweep = |_: u64, mut rng: Rng| {
        // ±3% on the wire length, which is a real tolerance on a wound part.
        let length = 62.0 * (1.0 + 0.03 * (2.0 * rng.unit() - 1.0));
        let ambient = Temperature::celsius(25.0);

        let mut net = ThermalNetwork::new("motor");
        let w = net.node(
            "winding",
            Substance::copper(),
            Volume::cm3(18.0),
            Length::mm(2.0),
            ambient,
        );
        let h = net.node_losing_to(
            "housing",
            Substance::aluminium_6061(),
            Volume::cm3(220.0),
            Length::mm(4.0),
            ambient,
            Environment::still_air(ambient, Area::cm2(420.0)),
        );
        net.link(w, h, Conductance::w_per_k(0.9)).unwrap();
        net.absorbing(w).unwrap();

        let coil = Winding::of_copper("coil", Length::m(length), 0.35e-6, ambient)
            .driven_at(Current::a(2.0));
        // The electro-thermal fixed point, closed by the caller as it must be.
        let mut watts = coil.dissipation().to_si();
        for _ in 0..30 {
            let settled = net.steady_state(Power::w(watts)).unwrap();
            watts = coil.dissipation_at(settled.temperature(w)).to_si();
        }
        net.steady_state(Power::w(watts))
            .unwrap()
            .temperature(w)
            .to_si()
            - 273.15
    };

    let one = Ensemble::new(1234, 400).run(sweep);
    let twelve = Ensemble::new(1234, 400).with_threads(12).run(sweep);
    for (i, (a, b)) in one.iter().zip(&twelve).enumerate() {
        assert_eq!(a.to_bits(), b.to_bits(), "sample {i}: {a} against {b}");
    }

    // The sweep says something: a ±3% tolerance on wire length moves the winding temperature by
    // a few kelvin, which is the number a designer wants and cannot get from one run.
    let hottest = one.iter().cloned().fold(f64::MIN, f64::max);
    let coolest = one.iter().cloned().fold(f64::MAX, f64::min);
    assert!(
        hottest - coolest > 1.0,
        "the tolerance should matter; spread {:.3} K",
        hottest - coolest
    );
    assert!(
        (60.0..110.0).contains(&coolest),
        "the coolest sample settles at {coolest:.1} C"
    );
}
