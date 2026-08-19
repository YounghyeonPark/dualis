//! The second coupling, and the one that shows the first was not a special case.
//!
//! `pantometry-optics` publishes absorbed light as heat and `pantometry-thermal` consumes it.
//! That could have been an arrangement between two crates written together. So here
//! the publisher is `pantometry-mechanics` instead: a bouncing ball's dashpot dissipates
//! mechanical energy, and the same [`LumpedMass`] picks the joules up off the same
//! channel, with nothing changed on the thermal side and nothing added to the kernel.
//!
//! Three crates now meet on one string and a conservation audit. None of them names
//! another.

use glam::DVec3;
use pantometry::prelude::*;

fn steel_ball() -> Substance {
    // A small aluminium ball, because the kernel already knows its properties.
    Substance::aluminium_6061()
}

fn ball_volume() -> Volume {
    // 20 mm diameter.
    Volume::from_si(4.0 / 3.0 * std::f64::consts::PI * 0.01f64.powi(3))
}

fn ball_area() -> Area {
    Area::from_si(4.0 * std::f64::consts::PI * 0.01f64.powi(2))
}

fn ball_mass() -> Mass {
    steel_ball().mass_of(ball_volume())
}

/// A ball dropped onto a floor, with a dashpot stiff enough to lose most of the drop
/// in a handful of bounces.
fn bouncing() -> ContactSystem {
    ContactSystem::new(
        "ball",
        &[Body::new(
            ball_mass(),
            LengthVec::m(0.0, 0.0, 0.5),
            VelocityVec::ZERO,
        )],
        AccelerationVec::from_si(-DVec3::Z * G0.to_si()),
        Ground::floor(),
        Stiffness::from_si(2e5),
        // Critical damping for this ball on this spring is 2*sqrt(k m) = 95 N s/m, so
        // 20 is a damping ratio of 0.21 and a restitution of about 0.51 — a bounce that
        // visibly loses energy without stopping dead.
        Damping::from_si(20.0),
    )
}

/// The ball warms up, and it warms up by exactly the mechanical energy it lost.
fn warming_ball() -> LumpedMass {
    LumpedMass::new(
        "ball-temperature",
        steel_ball(),
        ball_volume(),
        Length::mm(10.0),
        Temperature::celsius(20.0),
        Environment::still_air(Temperature::celsius(20.0), ball_area()),
    )
}

/// Mechanics publishes, thermal consumes, and the kernel audits the crossing — with
/// neither domain crate mentioning the other.
#[test]
fn a_bouncing_ball_warms_itself() {
    let contact = bouncing();
    let start_energy = contact.mechanical_energy();
    let restitution = contact.restitution();
    assert!(
        restitution > 0.1 && restitution < 0.8,
        "the test needs a bounce that visibly loses energy, got e = {restitution:.3}"
    );

    let mut sim = Simulation::new(Schedule::Multirate)
        // First order at twenty steps per contact oscillation — see the note in
        // pantometry-mechanics on why this number is the integrator's and not the audit's.
        .conservation_tolerance(2e-2)
        .with(contact)
        .with(warming_ball());

    for _ in 0..200 {
        sim.advance(Time::ms(10.0))
            .expect("mechanical energy plus heat must be conserved across the bus");
    }
    assert!((sim.time().to_si() - 2.0).abs() < 1e-9);

    // Heat crossed, all of it was claimed, and none is sitting on the bus.
    let crossed = sim.bus().total_consumed(quantity::ENERGY);
    assert!(
        crossed > 0.0,
        "the dashpot should have dissipated something"
    );
    assert!(sim.bus().unclaimed().next().is_none());

    // And the joules that crossed are the mechanical energy that went missing.
    //
    // This used to read `start_energy - 0.0f64.max(0.0)`, which is start_energy: the final
    // energy was never computed, and what remained was the one-sided claim that the dashpot
    // cannot dissipate more than the ball began with. True, and nearly free — it would have
    // held for any dissipation between zero and everything. The comment promised an equality
    // and the assertion delivered a bound.
    //
    // `ContactSystem` now answers `as_any`, so the energy actually left is reachable and the
    // equality can be the thing that is checked.
    //
    // The tolerance, from both ends. The dashpot's work is accumulated as -F·dx with the force
    // taken at the old velocity and the displacement at the new one, which is first order, and
    // the semi-implicit step has a first-order energy error of its own — `ContactSystem`'s doc
    // measures that at about 2% over a second of bouncing, which is where this simulation's
    // `conservation_tolerance(2e-2)` comes from. Measured here: 0.054617 J crossed against
    // 0.055455 J lost, 1.51% apart. From the other end, publishing 5% more heat than the
    // dashpot dissipated lands at 3.41% and must fail. So 2.5e-2 — the legitimate error uses
    // 60% of it and the sabotage overruns it by a third. Both numbers are in the comment
    // because a tolerance with only one of them written down is a tolerance nobody can revise.
    let end_energy = sim
        .domain_as::<ContactSystem>("ball")
        .expect("the contact system is still in the simulation")
        .mechanical_energy();
    let lost = start_energy.to_si() - end_energy.to_si();
    assert!(
        lost > 0.0,
        "the ball must have lost mechanical energy, got {lost:.4} J"
    );
    assert!(
        (crossed / lost - 1.0).abs() < 2.5e-2,
        "the joules that crossed the bus are the joules that went missing: \
         {crossed:.4} J crossed against {lost:.4} J lost"
    );
}

/// The temperature rise, followed through to a number — and the number is tiny.
///
/// Worth asserting precisely because it is small and surprising. Mechanical energy
/// feels large and heat capacity is larger: an 11 g aluminium ball dropped half a
/// metre carries 55 mJ, its heat capacity is 10 J/K, and turning the whole fall into
/// heat warms it five thousandths of a kelvin. A model that reported a warm ball would
/// be wrong in a way nobody would question.
///
/// And the rise does not depend on the mass at all. `mgh/(m c_p)` cancels to
/// `gh/c_p`, so a bearing ball and a wrecking ball dropped the same distance heat by
/// the same amount — which is the sort of thing dimensions make visible and a pile of
/// f64 does not.
#[test]
fn the_rise_from_a_fall_is_tiny_and_independent_of_mass() {
    let mass = ball_mass();
    assert!(
        (mass.to_si() * 1e3 - 11.31).abs() < 0.05,
        "a 20 mm aluminium ball is 11.3 g, got {} g",
        mass.to_si() * 1e3
    );

    let drop = Length::m(0.5);
    let potential: Energy = Energy::from_si(mass.to_si() * G0.to_si() * drop.to_si());
    assert!(
        (potential.to_si() * 1e3 - 55.4).abs() < 0.5,
        "half a metre is 55 mJ, got {} mJ",
        potential.to_si() * 1e3
    );

    let capacity: HeatCapacity = steel_ball()
        .heat_capacity(ball_volume())
        .expect("aluminium has a specific heat");
    assert!(
        (capacity.to_si() - 10.13).abs() < 0.05,
        "the ball holds 10.1 J/K, got {:?}",
        capacity.to_si()
    );

    let rise: Temperature = potential / capacity;
    assert!(
        (rise.to_si() * 1e3 - 5.47).abs() < 0.05,
        "the whole drop is 5.5 mK, got {} mK",
        rise.to_si() * 1e3
    );
    assert!(
        rise.to_si() < 0.01,
        "a bouncing ball does not get warm, and the model should say so"
    );

    // g h / c_p, with the mass gone.
    let specific_heat = steel_ball()
        .thermal
        .expect("aluminium is thermal")
        .specific_heat;
    let closed_form = G0.to_si() * drop.to_si() / specific_heat.to_si();
    assert!(
        (rise.to_si() / closed_form - 1.0).abs() < 1e-12,
        "the rise should be g h / c_p exactly: {} vs {closed_form}",
        rise.to_si()
    );
}

/// An undamped contact publishes nothing, so a simulation with no consumer is still
/// balanced — the bus is only unbalanced when something was actually dropped.
#[test]
fn a_lossless_contact_publishes_nothing_to_lose() {
    let lossless = ContactSystem::new(
        "ball",
        &[Body::new(
            ball_mass(),
            LengthVec::m(0.0, 0.0, 0.5),
            VelocityVec::ZERO,
        )],
        AccelerationVec::from_si(-DVec3::Z * G0.to_si()),
        Ground::floor(),
        Stiffness::from_si(2e5),
        Damping::from_si(0.0),
    );
    assert!((lossless.restitution() - 1.0).abs() < 1e-15);

    let mut sim = Simulation::new(Schedule::Multirate)
        // 2e-2 for the reason given on ContactSystem::max_stable_dt: a penalty contact
        // is a non-smooth potential, so the symplectic energy bound does not apply and
        // each bounce shifts the total by O(h).
        .conservation_tolerance(2e-2)
        .with(lossless);
    for _ in 0..100 {
        sim.advance(Time::ms(10.0))
            .expect("nothing is dissipated, so nothing is unclaimed");
    }
    assert_eq!(sim.bus().total_consumed(quantity::ENERGY), 0.0);
}

/// Gravity as a domain alongside the rest, and the thing it brings that nothing else
/// does: a second conserved quantity for the kernel's audit to police.
#[test]
fn the_kernel_polices_momentum_as_readily_as_energy() {
    let mut sim = Simulation::new(Schedule::Multirate)
        .conservation_tolerance(1e-11)
        .with(NBody::new(
            "binary",
            &[
                Body::new(
                    Mass::kg(5e12),
                    LengthVec::m(-500.0, 0.0, 0.0),
                    VelocityVec::m_per_s(0.0, 0.4, 0.0),
                ),
                Body::new(
                    Mass::kg(5e12),
                    LengthVec::m(500.0, 0.0, 0.0),
                    VelocityVec::m_per_s(0.0, -0.4, 0.0),
                ),
            ],
        ));

    // A symmetric binary: the momentum is exactly zero and must stay exactly zero.
    for _ in 0..50 {
        sim.advance(Time::s(10.0))
            .expect("pairwise forces cancel, so momentum cannot move");
    }
    let ledger = sim.ledger();
    for axis in [
        pantometry::mechanics::conserved::MOMENTUM_X,
        pantometry::mechanics::conserved::MOMENTUM_Y,
        pantometry::mechanics::conserved::MOMENTUM_Z,
    ] {
        let p = ledger.get(axis).expect("the domain reports every axis");
        assert!(p.abs() < 1e-6, "{axis} drifted to {p:e}");
    }
}
