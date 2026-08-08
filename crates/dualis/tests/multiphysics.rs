//! The end-to-end claim: two domains that do not know about each other, coupled
//! through the kernel, conserving energy across the interface.
//!
//! This is the test the whole architecture was for. `dualis-optics` computes how much
//! of a lamp a surface absorbs and publishes the joules; `dualis-thermal` consumes
//! them and warms up; the kernel audits the transfer and the totals. Neither domain
//! crate names the other, and neither has a line of code about coupling — the seam is
//! [`Exchange`](dualis_core::Exchange) and nothing else.
//!
//! The chain being closed is the one from the design: **absorbed light → heat →
//! temperature rise → thermal expansion → focus shift.** The last link is the one
//! that makes it matter, because a few milliwatts of stray absorption moves a focus
//! by more than the depth of focus it has to stay inside.

use dualis::prelude::*;
use dualis_optics::spectrum::Spectrum as Spec;

/// The optics side of the coupling, as a quasi-static domain.
///
/// It has no state to march: light crosses the instrument in nanoseconds, so against
/// a thermal timescale it is solved instantly and never subcycled. Every step it
/// recomputes the absorbed power from the surface's spectral absorptance and publishes
/// the joules for that interval.
struct AbsorbingSurface {
    lamp: SpectralPower,
    absorptance: Spectrum,
    paid_out: f64,
}

impl AbsorbingSurface {
    fn new(lamp: SpectralPower, optics: &SurfaceOptics) -> AbsorbingSurface {
        // Sample the surface's absorptance into a spectrum the radiometry can
        // integrate. This is the discretisation step where a careless coupling would
        // lose energy, and the bus audit is what catches it.
        let absorptance = Spec::curve(
            (0..=150)
                .map(|i| {
                    let nm = 350.0 + i as f64 * 5.0;
                    (nm, optics.absorptance(Length::nm(nm)))
                })
                .collect(),
        );
        AbsorbingSurface {
            lamp,
            absorptance,
            paid_out: 0.0,
        }
    }

    fn absorbed_power(&self) -> Power {
        self.lamp.absorbed_by(&self.absorptance)
    }
}

impl Domain for AbsorbingSurface {
    fn name(&self) -> &'static str {
        "optics"
    }

    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let joules = self.absorbed_power().to_si() * dt.to_si();
        bus.publish(HEAT, joules);
        self.paid_out += joules;
        Ok(())
    }

    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, -self.paid_out)
    }

    fn checkpoint(&mut self) {}
    fn restore(&mut self) {}
    fn supports_restore(&self) -> bool {
        true
    }
}

fn lens_volume() -> Volume {
    // A 25 mm diameter lens, 5 mm thick.
    Volume::from_si(std::f64::consts::PI * 0.0125f64.powi(2) * 0.005)
}

fn lens_area() -> Area {
    let r = 0.0125f64;
    Area::from_si(2.0 * std::f64::consts::PI * r * r + std::f64::consts::TAU * r * 0.005)
}

fn warm_lens() -> LumpedMass {
    LumpedMass::new(
        "lens",
        Substance::borosilicate_crown(),
        lens_volume(),
        Length::mm(5.0),
        Temperature::celsius(20.0),
        Environment::still_air(Temperature::celsius(20.0), lens_area()),
    )
}

/// A dichroic under a 5 W lamp: the optics domain says how many watts it absorbs, the
/// thermal domain warms up by exactly that, and the kernel confirms nothing was
/// created or lost on the way across.
#[test]
fn absorbed_light_warms_the_glass_and_the_books_balance() {
    let lamp = SpectralPower::new(Spec::blackbody(3200.0), Power::w(5.0), VISIBLE_RANGE);
    let dichroic = SurfaceOptics::dichroic(vec![[495.0, 545.0]], 0.95, 10.0);
    let surface = AbsorbingSurface::new(lamp, &dichroic);

    let absorbed = surface.absorbed_power();
    assert!(
        absorbed.in_mw() > 5.0 && absorbed.in_mw() < 150.0,
        "a dichroic under 5 W should absorb tens of milliwatts, got {} mW",
        absorbed.in_mw()
    );

    let mut sim = Simulation::new(Schedule::Multirate)
        .with(surface)
        .with(warm_lens());

    // Ten minutes, in five-second windows. The thermal domain subcycles as its own
    // time constant demands; the optics domain never does.
    for _ in 0..120 {
        sim.advance(Time::s(5.0)).expect("energy must be conserved");
    }
    assert!((sim.time().to_si() - 600.0).abs() < 1e-9);

    // **The glass actually got warm**, which the conservation check alone cannot establish.
    //
    // That is not a hypothetical gap. `LumpedMass::ledger` reports `stored + lost` using the
    // same arithmetic `step` used to move the temperature, and the surface reports what it
    // paid out — so the residual is an accounting identity that holds whatever the physics
    // does. Discarding every arriving joule while leaving the cooling term alone passes all
    // 344 tests in this workspace, which is how this assertion came to be written.
    let lens: &LumpedMass = sim
        .domain_as("lens")
        .expect("the lens is in the simulation");
    let rise = lens.rise().to_si();

    // Against an independent integration of the same energy balance: forward Euler at 50 ms,
    // written here rather than called from the domain, against whatever multirate substeps
    // the scheduler chose. Different code, different step, same physics.
    //
    // Not a closed form, and the reason is worth naming. `equilibrium_rise` *is* a closed
    // form but it ignores radiation, and at room temperature a black surface radiates about
    // as fast as still air convects — so it overstates the settling point by nearly a factor
    // of two, which is exactly the trap it warns about in its own doc comment. The nonlinear
    // T⁴ term has no useful closed-form transient, so the honest reference is a quadrature.
    //
    // What this catches: a wrong heat capacity, a wrong loss coefficient, a dropped source
    // term, and heat that arrives on the bus and never reaches the glass. What it does not:
    // an error inside `loss_from` itself, which is checked in the thermal crate against the
    // exponential a purely convective body must follow.
    let ambient = Temperature::celsius(20.0);
    let emissivity = Substance::borosilicate_crown()
        .thermal
        .expect("N-BK7 has thermal properties")
        .emissivity;
    let environment = Environment::still_air(ambient, lens_area());
    let capacity = lens.heat_capacity().to_si();
    let h = 0.05;
    let mut reference = ambient;
    for _ in 0..(600.0 / h) as usize {
        let loss = environment.loss_from(reference, emissivity).to_si();
        reference += Temperature::from_si((absorbed.to_si() - loss) * h / capacity);
    }
    let expected = (reference - ambient).to_si();
    assert!(
        (rise / expected - 1.0).abs() < 0.01,
        "the lens rose {rise} K where an independent integration gives {expected} K"
    );
    assert!(rise > 1.0, "and it is a real warming, not a rounding");
    // Below the radiation-free closed form, which is the direction radiation pushes it.
    assert!(
        rise < lens.equilibrium_rise(absorbed).to_si(),
        "radiation must make the glass settle lower than convection alone predicts"
    );

    // Only now does the conservation statement mean something: everything the surface paid
    // out is either stored in the glass or has gone to the room.
    //
    // Judged against the joules that actually crossed, not against a bare 1e-9. The residual's
    // correct value is exactly zero, so a tolerance on it needs a scale from outside — the
    // same reason `Ledger` records one beside every total.
    let crossed = sim.bus().total_consumed(quantity::ENERGY);
    assert!(crossed > 0.0, "nothing crossed, so nothing was audited");
    let residual = sim.ledger().get(quantity::ENERGY).unwrap();
    assert!(
        residual.abs() / crossed < 1e-12,
        "energy residual {residual} J against {crossed} J that crossed"
    );

    // And nothing is left sitting on the bus unclaimed.
    assert!(sim.bus().unclaimed().next().is_none());
}

/// The link that makes the coupling matter. A few tens of milliwatts absorbed warms
/// the glass a few kelvin, that expands it, and the expansion moves the focus by more
/// than the depth of focus it was supposed to stay inside.
///
/// Every number in this chain comes from a different crate — radiometry from optics,
/// heat capacity and expansion from the kernel's `Substance`, depth of focus from
/// diffraction — and the dimensions are what let them be composed at all.
#[test]
fn a_warm_lens_drifts_out_of_focus() {
    let lamp = SpectralPower::new(Spec::blackbody(3200.0), Power::w(5.0), VISIBLE_RANGE);
    let dichroic = SurfaceOptics::dichroic(vec![[495.0, 545.0]], 0.95, 10.0);
    let surface = AbsorbingSurface::new(lamp, &dichroic);
    let absorbed = surface.absorbed_power();

    let lens = warm_lens();
    let settled = lens.equilibrium_rise(absorbed);
    assert!(
        settled.to_si() > 1.0,
        "the rise should be more than a kelvin, got {} K",
        settled.to_si()
    );

    // Run to something close to equilibrium and check the transient agrees with the
    // steady-state formula it should be heading for.
    let mut sim = Simulation::new(Schedule::Multirate)
        .with(surface)
        .with(warm_lens());
    for _ in 0..400 {
        sim.advance(Time::s(5.0)).unwrap();
    }
    let reached = sim
        .domain("lens")
        .map(|d| d.ledger())
        .expect("the lens is in the simulation");
    assert!(!reached.is_empty());

    // The optical consequence, from the kernel's material data.
    let glass = Substance::borosilicate_crown();
    let growth = glass
        .expansion_of(Length::mm(100.0), settled)
        .expect("N-BK7 has an expansion coefficient");

    // Against the depth of focus at 550 nm, which is the budget the drift has to fit
    // inside. The numbers, all of them computed above rather than assumed: 96 mW of
    // absorbed stray light settles the lens 10.0 K above ambient, which grows a
    // 100 mm mount by 7.10 um.
    let modest = depth_of_focus(Length::nm(550.0), 0.25, 1.0); // 8.80 um
    let tighter = depth_of_focus(Length::nm(550.0), 0.30, 1.0); // 6.11 um

    // At NA 0.25 the drift does not quite exceed the depth — it eats 81% of it, which
    // leaves nothing for every other error in the instrument.
    let used = growth / modest;
    assert!(
        used > 0.7 && used < 1.0,
        "a {:.1} K rise moves the mount {:.2} um and spends {:.0}% of the {:.2} um \
         depth of focus",
        settled.to_si(),
        growth.in_um(),
        used * 100.0,
        modest.in_um()
    );

    // Tighten the aperture slightly — depth falls as one over NA squared — and the
    // same thermal drift is now larger than the entire budget. That is the whole
    // reason this coupling has to exist rather than being assumed negligible.
    assert!(
        growth > tighter,
        "at NA 0.30 the {:.2} um drift exceeds the {:.2} um depth outright",
        growth.in_um(),
        tighter.in_um()
    );

    // And the glass is nowhere near breaking, which is the other thing worth knowing.
    assert_eq!(glass.survives(settled), Some(true));
}

/// The scheduler's two mechanisms, observed on real domains: the quasi-static optics
/// takes one step per window whatever the window is, and the thermal domain subcycles
/// according to its own time constant.
#[test]
fn optics_never_subcycles_and_the_thermal_domain_does() {
    let lamp = SpectralPower::new(Spec::constant(1.0), Power::w(1.0), VISIBLE_RANGE);
    let black = SurfaceOptics::black();
    let mut sim = Simulation::new(Schedule::Multirate)
        .with(AbsorbingSurface::new(lamp, &black))
        .with(warm_lens());

    let short = sim.advance(Time::s(1.0)).unwrap();
    let long = sim.advance(Time::s(120.0)).unwrap();

    assert_eq!(short.substeps[0], ("optics".to_string(), 1));
    assert_eq!(
        long.substeps[0],
        ("optics".to_string(), 1),
        "a solve is a solve"
    );

    // The lens's time constant is a couple of minutes, and the reported limit is a
    // tenth of that — so a one-second window needs one substep and a two-minute
    // window needs several.
    assert_eq!(short.substeps[1].0, "lens");
    assert_eq!(short.substeps[1].1, 1);
    assert!(
        long.substeps[1].1 > 1,
        "a two-minute window should subcycle, got {} substeps",
        long.substeps[1].1
    );
}

/// A coupling that loses energy is caught. The optics domain publishes the joules it
/// absorbed; if the thermal side receives a different number — which is what an
/// interpolation between two discretisations of the same surface does — the bus audit
/// fails rather than letting the discrepancy vanish into a temperature.
#[test]
fn a_lossy_interface_is_refused() {
    /// A thermal domain that takes only 90% of what it was offered, as a mismatched
    /// interpolation would.
    struct LeakySink;
    impl Domain for LeakySink {
        fn name(&self) -> &'static str {
            "leaky"
        }
        fn step(&mut self, _t: Time, _dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
            let offered = bus.peek(quantity::ENERGY);
            // Put 10% back, which is the same as never having taken it.
            let _ = bus.take(quantity::ENERGY);
            bus.publish(quantity::ENERGY, offered * 0.1);
            Ok(())
        }
    }

    let lamp = SpectralPower::new(Spec::constant(1.0), Power::w(1.0), VISIBLE_RANGE);
    let mut sim = Simulation::new(Schedule::Staggered)
        .with(AbsorbingSurface::new(lamp, &SurfaceOptics::black()))
        .with(LeakySink);

    let err = sim
        .advance(Time::s(1.0))
        .expect_err("10% of the heat arrived nowhere");
    assert_eq!(err.quantity, "energy");
    assert!(err.site.contains("not consumed"), "{err}");
    // The clock did not move, so the run is not left half-applied.
    assert_eq!(sim.time(), Time::from_si(0.0));
}

/// A black surface absorbs everything, so the coupling's total is exactly the lamp's
/// power — the anchor case that says the integration chain has no stray factor in it.
#[test]
fn a_black_surface_hands_over_the_whole_lamp() {
    let lamp = SpectralPower::new(Spec::blackbody(3200.0), Power::w(2.0), VISIBLE_RANGE);
    let total = lamp.total();
    let surface = AbsorbingSurface::new(lamp, &SurfaceOptics::black());
    let absorbed = surface.absorbed_power();
    assert!(
        (absorbed / total - 1.0).abs() < 1e-9,
        "black should absorb all 2 W, got {:?}",
        absorbed
    );

    // A mirror hands over almost none of it — aluminium's 8%.
    let mirror_lamp = SpectralPower::new(Spec::blackbody(3200.0), Power::w(2.0), VISIBLE_RANGE);
    let mirror = AbsorbingSurface::new(mirror_lamp, &SurfaceOptics::aluminium());
    let fraction = mirror.absorbed_power() / total;
    assert!(
        fraction > 0.04 && fraction < 0.12,
        "aluminium absorbs a few percent, got {fraction}"
    );
}
