//! Every domain reports its own scalars, so no layer has to know them by name.
//!
//! `Domain::readings` exists because the alternative was a `match` over domain types in whatever
//! layer wanted a table — which has to be edited every time a physics is added, and that is the
//! one thing this workspace's structure exists to avoid.

use pantometry::prelude::*;
use pantometry_electrical::Winding;

/// **A domain names its own numbers, and a boxed one still does.**
///
/// `Simulation` stores `Vec<Box<dyn Domain>>`, so the forwarding impl is the path every layer
/// above will actually take. A `Box` that forgot to forward would make every reading vanish for
/// exactly the callers who need them, and nothing would look broken.
#[test]
fn a_boxed_domain_still_reports_its_readings() {
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

    // Directly, and through the box the simulation will hold it in.
    let direct = net.readings();
    let boxed: Box<dyn Domain> = Box::new(net);
    assert_eq!(direct, boxed.readings(), "the box swallowed the readings");

    assert_eq!(direct.len(), 2);
    assert_eq!(direct[0].label, "winding");
    assert_eq!(direct[1].label, "housing");
    assert!(direct.iter().all(|r| r.domain == "motor" && r.unit == "C"));
}

/// **Each domain reports what it is *for*, not a uniform summary.**
///
/// The reason this is a trait method and not a generic reduction over whatever a domain exposes.
/// A room's mean pressure is zero by symmetry; a bar's mean hides the gradient that is its whole
/// point; a winding's resistance is the *reason* its dissipation moves. Only the domain knows.
#[test]
fn a_domain_reports_the_number_it_exists_to_produce() {
    // A bar reports both ends, because the difference is why it is not a lumped mass.
    let bar = Bar1D::new(
        "bar",
        Substance::aluminium_6061(),
        41,
        Length::mm(0.5),
        Area::mm2(100.0),
        Temperature::celsius(20.0),
    );
    let got = bar.readings();
    let labels: Vec<&str> = got.iter().map(|r| r.label.as_str()).collect();
    assert_eq!(labels, ["mean", "peak", "absorbed"]);

    // A room reports its peak, not its mean — released in a mode, the mean is zero and the peak
    // is the amplitude.
    let room = Room::of_air("room", Length::m(4.4), Length::m(3.1), 41).released_in_mode(
        1,
        1,
        Pressure::from_si(1.0),
    );
    let r = room.readings();
    assert_eq!(r.len(), 1);
    assert_eq!(r[0].label, "peak");
    assert!(
        (r[0].value - 1.0).abs() < 1e-9,
        "a mode released at 1 Pa should peak at 1 Pa, got {}",
        r[0].value
    );

    // A winding reports the resistance beside the power, because at fixed current the power
    // ratio *is* the resistance ratio and one column cannot show that.
    let coil = Winding::of_copper("coil", Length::m(24.0), 0.35e-6, Temperature::celsius(20.0))
        .driven_at(Current::a(2.0));
    let got = coil.readings();
    let labels: Vec<&str> = got.iter().map(|r| r.label.as_str()).collect();
    assert_eq!(labels, ["dissipating", "resistance", "spent"]);
}

/// A domain that says nothing says nothing, rather than something wrong.
#[test]
fn the_default_is_empty_and_that_is_visible() {
    struct Quiet;
    impl Domain for Quiet {
        fn name(&self) -> &str {
            "quiet"
        }
        fn step(&mut self, _t: Time, _dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
            Ok(())
        }
    }
    assert!(Quiet.readings().is_empty());
    // Which is the same hazard `as_any` and `as_field` have taught twice: a domain that forgets
    // this is silently absent from every table rather than broken. Absent columns are at least
    // visible in a header; that is the whole mitigation and it is worth saying out loud.
}
