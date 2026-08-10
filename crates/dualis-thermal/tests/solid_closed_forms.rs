//! `Solid3D` against closed forms, not against a second copy of the same sweep.
//!
//! Three dimensions is where a conduction scheme goes wrong quietly. A swapped axis, a stencil
//! missing one of its six arms, or a spacing squared in the wrong place all produce a block that
//! still cools, still conserves, and still looks like heat spreading — so the checks here are
//! chosen to be sensitive to exactly those, and every one of them is against something the sweep
//! did not compute.

use dualis_core::conserved::quantity;
use dualis_core::units::{Energy, Length, Temperature, Time};
use dualis_core::{Domain, Exchange, Schedule, Simulation, Substance};
use dualis_thermal::{Solid3D, STABLE_FOURIER_3D};

fn block(counts: (usize, usize, usize), dx_mm: f64) -> Solid3D {
    Solid3D::new(
        "block",
        Substance::aluminium_6061(),
        counts,
        Length::mm(dx_mm),
        Temperature::celsius(20.0),
    )
}

/// Step a block on its own, with nothing on the bus.
fn run(block: &mut Solid3D, dt: Time, steps: usize) {
    let mut bus = Exchange::new();
    let mut t = Time::from_si(0.0);
    for _ in 0..steps {
        block.step(t, dt, &mut bus).expect("stable and isolated");
        t += dt;
    }
}

/// **A separable cosine mode decays by exactly its own eigenvalue, every step.**
///
/// This is the load-bearing test, and it is an *equality* rather than a tolerance. The mode
/// `cos(aπ(i+½)/nx)·cos(bπ(j+½)/ny)·cos(cπ(k+½)/nz)` is an exact eigenvector of the seven-point
/// stencil with mirrored faces, so after `n` steps its amplitude is `A₀·(1 + F·λdx²)ⁿ` with
/// nothing dropped — no discretisation error, no truncation, no leading order.
///
/// Three modes, chosen so they cannot pass for each other. `(1,0,0)` and `(0,1,0)` differ only by
/// which axis varies and would be identical on a cubic block, so the block is **7 × 5 × 3** —
/// deliberately unequal, so a stencil that read y where it meant x predicts the wrong number.
/// `(2,1,1)` varies along all three at once, which is the only case where a missing arm shows.
#[test]
fn a_separable_mode_decays_at_its_exact_discrete_rate() {
    for mode in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 1, 1)] {
        let mut b = block((7, 5, 3), 2.0);
        b.release_mode(mode, Temperature::celsius(20.0), 4.0);

        let dt = Time::from_si(b.max_stable_dt(Time::from_si(0.0)).to_si() * 0.8);
        let per_step = b.mode_amplification(mode, dt);
        assert!(
            per_step.abs() <= 1.0,
            "mode {mode:?} amplifies by {per_step}, which is the instability itself"
        );

        let steps = 25;
        run(&mut b, dt, steps);

        let want = 4.0 * per_step.powi(steps as i32);
        let got = b.mode_amplitude(mode);
        // Machine precision, scaled by the amplitude the mode started at. Nothing here is
        // approximated, so anything above rounding is a wrong operator.
        assert!(
            (got - want).abs() < 1e-12 * 4.0,
            "mode {mode:?}: amplitude {got:.15e} against the exact {want:.15e}"
        );
        // And it really decayed, or the check above is satisfied by a block that did nothing.
        assert!(got.abs() < 0.9 * 4.0, "mode {mode:?} barely moved: {got}");
    }
}

/// **The three axes are interchangeable, and the block is not cubic.**
///
/// The same physical block described three ways: `7×5×3` released in `(1,0,0)`, `5×3×7` in
/// `(0,0,1)`, `3×7×5` in `(0,1,0)`. Each is the same mode along the axis of length 7, so all
/// three must decay identically.
///
/// This is the check a swapped index survives every *other* test in this file: a stencil that
/// consistently read y for x still conserves, still decays, and still matches its own
/// eigenvalue, because the eigenvalue would be computed with the same mistake.
#[test]
fn the_axes_are_the_same_physics_in_a_different_order() {
    let amplitude_after = |counts, mode| {
        let mut b = block(counts, 2.0);
        b.release_mode(mode, Temperature::celsius(20.0), 4.0);
        // Near the stability limit and for long enough that the mode is mostly gone. At the
        // 1e-4 s this first said, forty steps decayed it by 1.3% — three orientations agreeing
        // to 1e-12 about a block that had barely moved, which is agreement about nothing. The
        // guard below is what caught that, and it is why the guard is there.
        let dt = Time::from_si(5e-3);
        run(&mut b, dt, 200);
        b.mode_amplitude(mode)
    };

    let along_x = amplitude_after((7, 5, 3), (1, 0, 0));
    let along_z = amplitude_after((5, 3, 7), (0, 0, 1));
    let along_y = amplitude_after((3, 7, 5), (0, 1, 0));

    assert!(
        (along_x - along_z).abs() < 1e-12 && (along_x - along_y).abs() < 1e-12,
        "the long axis should not matter: x {along_x:.12e}, y {along_y:.12e}, z {along_z:.12e}"
    );
    assert!(
        along_x.abs() < 0.2 * 4.0,
        "the mode should be mostly gone, or three orientations agree about nothing: {along_x}"
    );
}

/// **The discrete decay rate approaches the continuum one at second order.**
///
/// The exact test above pins the scheme to itself. This one pins the scheme to the *physics*: the
/// slowest mode of an insulated block of size `L` decays at `α·π²/L²`, and the discrete
/// eigenvalue must converge to that as the grid refines.
///
/// The **rate** and not the value. A first-order scheme also converges, and its error also gets
/// small — the thing that separates them is that halving `dx` quarters a second-order error and
/// only halves a first-order one. That distinction is what found the acoustic boundary defect in
/// this workspace, where an error of 1.4% looked exactly like ordinary discretisation.
#[test]
fn the_rate_converges_to_the_continuum_at_second_order() {
    let length = 20e-3;
    let alpha = Substance::aluminium_6061()
        .diffusivity()
        .expect("aluminium conducts")
        .to_si();
    let exact = alpha * std::f64::consts::PI.powi(2) / (length * length);

    let error_at = |n: usize| {
        let dx = length / n as f64;
        let b = block((n, 1, 1), dx * 1e3);
        let dt = Time::from_si(1e-6);
        // The decay rate the discrete operator actually has: -ln(amplification)/dt.
        let rate = -b.mode_amplification((1, 0, 0), dt).ln() / dt.to_si();
        (rate / exact - 1.0).abs()
    };

    let (coarse, fine) = (error_at(8), error_at(16));
    let ratio = coarse / fine;
    assert!(
        (ratio - 4.0).abs() < 0.15,
        "second order means the error quarters on refinement: {coarse:.3e} -> {fine:.3e}, \
         ratio {ratio:.4} against 4"
    );
    // And it is converging *to* the right number, not merely converging.
    assert!(fine < 0.01, "16 cells should be within 1%, was {fine:.4e}");
}

/// **An insulated block conserves every joule, judged against a scale.**
///
/// The total is not compared against itself — a block that started uniform and stayed uniform
/// would pass that. It is compared against the largest cell-to-cell difference the run ever held,
/// so the tolerance is measured against the size of the thing that had to cancel.
#[test]
fn an_insulated_block_conserves_exactly() {
    let mut b = block((9, 7, 5), 1.5);
    b.release_mode((2, 1, 1), Temperature::celsius(20.0), 30.0);

    let capacity = Substance::aluminium_6061()
        .heat_capacity(b.volume())
        .expect("aluminium has a specific heat")
        .to_si();
    let cells = b.counts().0 * b.counts().1 * b.counts().2;
    let mean_at = |b: &Solid3D| b.mean_temperature().to_si();

    let before = mean_at(&b);
    let spread = b.peak_temperature().to_si() - b.coldest_temperature().to_si();
    assert!(
        spread > 50.0,
        "the mode should be a real gradient: {spread} K"
    );

    run(&mut b, Time::from_si(1e-4), 600);

    let after = mean_at(&b);
    // The scale is the energy that was moving about, not the energy that is there: the block
    // holds about 1 kJ above absolute zero and the mode is worth a few joules of that, so a
    // relative check against the total would hide a leak of the entire gradient.
    let moved = capacity * spread;
    let lost = capacity * (after - before).abs();
    assert!(
        lost < 1e-9 * moved,
        "insulated: {lost:.3e} J lost against {moved:.3e} J of gradient over {cells} cells"
    );
    // The gradient really did relax, or nothing was asked of the sweep.
    let left = b.peak_temperature().to_si() - b.coldest_temperature().to_si();
    assert!(
        left < 0.2 * spread,
        "the gradient should be mostly gone, or nothing was asked of the sweep:          {spread:.3} K -> {left:.3} K"
    );
}

/// **The stability limit is one sixth, it is a third of the bar's, and it is refused.**
///
/// The number that makes three dimensions expensive, stated three ways: the constant, the
/// reported `max_stable_dt`, and what happens to a caller who ignores it.
#[test]
fn the_third_dimension_costs_a_factor_of_three_and_the_limit_is_enforced() {
    let b = block((5, 5, 5), 2.0);
    let dx = 2e-3;
    let alpha = Substance::aluminium_6061().diffusivity().unwrap().to_si();

    let limit = b.max_stable_dt(Time::from_si(0.0)).to_si();
    let exact = dx * dx / (6.0 * alpha);
    assert!(
        (limit / exact - 1.0).abs() < 1e-12,
        "dx²/6α: {limit:.9e} against {exact:.9e}"
    );
    assert!((b.fourier_number(Time::from_si(limit)) - STABLE_FOURIER_3D).abs() < 1e-12);

    // A third of the 1D limit of dx²/2α, which is the whole cost of the extra two axes.
    let bar_limit = dx * dx / (2.0 * alpha);
    assert!(
        (bar_limit / limit - 3.0).abs() < 1e-9,
        "three dimensions should cost exactly 3× in steps, got {:.6}",
        bar_limit / limit
    );

    // Refuse rather than diverge, and name the number that was broken.
    let mut b = block((5, 5, 5), 2.0);
    let mut bus = Exchange::new();
    let over = Time::from_si(limit * 1.05);
    let err = b
        .step(Time::from_si(0.0), over, &mut bus)
        .expect_err("5% past the limit must be refused");
    assert_eq!(err.quantity, "Fourier number");
    assert!(
        (err.after / STABLE_FOURIER_3D - 1.05).abs() < 1e-9,
        "the violation should say by how much: {}",
        err.after
    );

    // And just inside it is accepted, so the check is a limit and not a blanket refusal.
    let mut b = block((5, 5, 5), 2.0);
    b.step(
        Time::from_si(0.0),
        Time::from_si(limit),
        &mut Exchange::new(),
    )
    .expect("exactly at the limit is stable");
}

/// **Heat with no place goes nowhere in particular.**
///
/// `Bar1D` puts placeless heat in its first cell, which is defensible for a bar: it has an end,
/// and a surface absorbing light would plausibly be there. A block has six faces and no
/// distinguished cell, so choosing one would invent a location the bus never carried — and a hot
/// spot that came out of a tie-break is worse than none, because it looks like physics.
///
/// So a uniform block fed from the plain channel must stay uniform, exactly, and rise by the
/// amount the whole block's heat capacity says.
#[test]
fn placeless_heat_leaves_the_block_uniform() {
    let mut b = block((4, 3, 2), 3.0);
    let mut bus = Exchange::new();
    let joules = 12.0;
    bus.publish(quantity::ENERGY, joules);
    b.step(Time::from_si(0.0), Time::from_si(1e-5), &mut bus)
        .expect("stable");

    let (nx, ny, nz) = b.counts();
    let first = b.temperature_at(0, 0, 0).to_si();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let t = b.temperature_at(i, j, k).to_si();
                assert!(
                    (t - first).abs() < 1e-12,
                    "cell ({i},{j},{k}) is at {t} against {first} — placeless heat found a place"
                );
            }
        }
    }

    let capacity = Substance::aluminium_6061()
        .heat_capacity(b.volume())
        .unwrap()
        .to_si();
    let rise = first - Temperature::celsius(20.0).to_si();
    assert!(
        (rise - joules / capacity).abs() < 1e-12,
        "{joules} J into {capacity:.4} J/K should be {:.9} K, was {rise:.9}",
        joules / capacity
    );
    assert!((b.absorbed_energy().to_si() - joules).abs() < 1e-12);
}

/// **A hot spot spreads sideways, which is the entire reason this domain exists.**
///
/// One cell heated in the middle of a plate. A one-dimensional model of the same plate has
/// nowhere for the heat to go but along, so the check is that the *transverse* neighbours warm
/// too, and by the same amount as the in-line ones — six equal arms is what an isotropic stencil
/// on a cubic grid means, and an unequal one is a spacing used in the wrong place.
#[test]
fn a_hot_spot_spreads_the_same_way_along_every_axis() {
    let mut b = block((9, 9, 9), 1.0);
    b.deposit(4, 4, 4, Energy::from_si(2.0));
    let hot = b.temperature_at(4, 4, 4).to_si();
    assert!(hot > 20.0 + 273.15, "the deposit should have warmed it");

    run(&mut b, Time::from_si(2e-6), 5);

    let arms = [
        b.temperature_at(3, 4, 4).to_si(),
        b.temperature_at(5, 4, 4).to_si(),
        b.temperature_at(4, 3, 4).to_si(),
        b.temperature_at(4, 5, 4).to_si(),
        b.temperature_at(4, 4, 3).to_si(),
        b.temperature_at(4, 4, 5).to_si(),
    ];
    let ambient = Temperature::celsius(20.0).to_si();
    for (n, arm) in arms.iter().enumerate() {
        assert!(*arm > ambient + 1e-6, "arm {n} never warmed: {arm}");
        assert!(
            (arm - arms[0]).abs() < 1e-12 * (arms[0] - ambient),
            "arm {n} at {arm} against arm 0 at {}: the stencil is not isotropic",
            arms[0]
        );
    }
    // And the corner, which is only reachable by two steps, is warmer than nothing and colder
    // than a face — a plain check that heat is diffusing rather than being copied about.
    let corner = b.temperature_at(3, 3, 4).to_si();
    assert!(
        corner > ambient && corner < arms[0],
        "diagonal {corner} should sit between ambient {ambient} and the face {}",
        arms[0]
    );
}

/// **A block one cell thick in two directions is a bar, exactly.**
///
/// The reduction that makes the three-dimensional operator checkable against a one-dimensional
/// closed form. An axis of one cell is a mirror against itself, contributes nothing to the
/// eigenvalue, and must leave the remaining axis behaving as the 1D scheme does — including its
/// looser stability limit being *unavailable*, since this domain still reports `dx²/6α`.
#[test]
fn one_cell_thick_reduces_to_a_bar() {
    let n = 12;
    let mut b = block((n, 1, 1), 2.0);
    b.release_mode((1, 0, 0), Temperature::celsius(20.0), 5.0);

    let dt = Time::from_si(2e-5);
    let f = b.fourier_number(dt);
    // The 1D eigenvalue, written out here rather than taken from the domain.
    let lambda = -4.0 * (std::f64::consts::PI / (2.0 * n as f64)).sin().powi(2);
    let want_per_step = 1.0 + f * lambda;
    assert!(
        (b.mode_amplification((1, 0, 0), dt) - want_per_step).abs() < 1e-15,
        "the flat axes should contribute nothing to the eigenvalue"
    );

    run(&mut b, dt, 30);
    let want = 5.0 * want_per_step.powi(30);
    assert!(
        (b.mode_amplitude((1, 0, 0)) - want).abs() < 1e-12 * 5.0,
        "1D limit: {:.15e} against {want:.15e}",
        b.mode_amplitude((1, 0, 0))
    );
}

/// **The field reads back the cells, and does not extrapolate past a face.**
///
/// A trilinear sampler is four lerps and an index, and it is very easy to write one that is
/// subtly wrong at a face or that transposes two axes on the way in. Cell centres must come back
/// exactly; a point outside must clamp rather than continue the gradient, because outside an
/// insulated face the temperature is not defined and extrapolating would draw a block hotter
/// than any cell in it.
#[test]
fn the_field_reads_back_the_cells_and_clamps_at_the_faces() {
    use dualis_core::ScalarField;
    let mut b = block((5, 4, 3), 2.0);
    b.release_mode((1, 1, 1), Temperature::celsius(20.0), 10.0);
    let t = Time::from_si(0.0);

    let (nx, ny, nz) = b.counts();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let at = b.at(b.centre_of(i, j, k), t);
                let cell = b.temperature_at(i, j, k).to_si();
                assert!(
                    (at - cell).abs() < 1e-9,
                    "centre ({i},{j},{k}): field {at} against cell {cell}"
                );
            }
        }
    }

    // Outside, in every direction, is the nearest cell and never past it.
    let hottest = b.peak_temperature().to_si();
    let coldest = b.coldest_temperature().to_si();
    for p in [
        dualis_core::units::LengthVec::m(-1.0, 0.0, 0.0),
        dualis_core::units::LengthVec::m(0.0, -1.0, 0.0),
        dualis_core::units::LengthVec::m(0.0, 0.0, -1.0),
        dualis_core::units::LengthVec::m(1.0, 1.0, 1.0),
    ] {
        let v = b.at(p, t);
        assert!(
            v <= hottest + 1e-9 && v >= coldest - 1e-9,
            "sampling outside gave {v}, past the block's own range [{coldest}, {hottest}]"
        );
    }
    assert_eq!(b.unit(), "K", "the cells hold kelvin");
}

/// **The block audits clean inside a `Simulation`, subcycling to its own limit.**
///
/// The domain in the place it will actually be used: a source publishing joules, a block taking
/// them, and the kernel's audit live for the whole run. A frame interval a hundred times the
/// block's stability limit, so `Multirate` has to subcycle it — which is the schedule that had a
/// first-order defect in this workspace and now has to hand a subcycling consumer its share.
#[test]
fn it_runs_and_audits_inside_a_simulation() {
    // A source with **books**. The first version of this published joules and reported an empty
    // ledger, so the audit saw the total go 0 → 4 J and refused the run — correctly. A domain
    // that pays out has to say what it has left, or every joule it hands over looks created.
    struct Source {
        watts: f64,
        reserve: f64,
    }
    impl Domain for Source {
        fn name(&self) -> &str {
            "lamp"
        }
        fn kind(&self) -> dualis_core::Kind {
            dualis_core::Kind::QuasiStatic
        }
        fn step(
            &mut self,
            _t: Time,
            dt: Time,
            bus: &mut Exchange,
        ) -> Result<(), dualis_core::Violation> {
            let joules = (self.watts * dt.to_si()).min(self.reserve);
            self.reserve -= joules;
            bus.publish(quantity::ENERGY, joules);
            Ok(())
        }
        fn ledger(&self) -> dualis_core::Ledger {
            dualis_core::Ledger::new().with(quantity::ENERGY, self.reserve)
        }
    }

    let seconds = 0.5;
    let watts = 8.0;
    let mut sim = Simulation::new(Schedule::Multirate)
        .with(Source {
            watts,
            // Comfortably more than the run spends, so the source never runs dry and the
            // absorbed total below is the full `watts × seconds`.
            reserve: watts * seconds * 10.0,
        })
        .with(block((6, 6, 6), 2.0));

    sim.advance(Time::from_si(seconds))
        .expect("an insulated block taking every joule cannot leak");

    let b = sim
        .domain_as::<Solid3D>("block")
        .expect("the block is still there");
    let paid = watts * seconds;
    assert!(
        (b.absorbed_energy().to_si() - paid).abs() < 1e-9 * paid,
        "{paid} J published, {} J absorbed",
        b.absorbed_energy().to_si()
    );

    // Uniform, because the source publishes on the plain channel and that carries no place.
    let spread = b.peak_temperature().to_si() - b.coldest_temperature().to_si();
    assert!(
        spread < 1e-9,
        "placeless heat made a gradient of {spread} K"
    );

    // And the rise is what the block's own heat capacity says.
    let capacity = Substance::aluminium_6061()
        .heat_capacity(b.volume())
        .unwrap()
        .to_si();
    let rise = b.mean_temperature().to_si() - Temperature::celsius(20.0).to_si();
    assert!(
        (rise - paid / capacity).abs() < 1e-9 * (paid / capacity),
        "{rise:.6} K against {:.6} K",
        paid / capacity
    );
}

/// **At the reported limit the worst-resolved mode does not grow, and past it, it does.**
///
/// The limit stated as physics rather than as a constant. `STABLE_FOURIER_3D` being `1/6` is a
/// claim about what the *sharpest* thing on the grid does — the mode that alternates fastest,
/// whose eigenvalue approaches `-12/dx²` when all three axes contribute. A block released in a
/// smooth mode is stable at three times the limit and says nothing about it.
///
/// Written because a mutation exposed the hole: replacing `1/6` with `1/2` left nine of the ten
/// tests here green, since none of the others ever excited a mode sharp enough to care.
#[test]
fn the_limit_is_where_the_sharpest_mode_stops_growing() {
    let n = 5;
    // The highest mode this grid can represent — very nearly one alternation per cell.
    let sharpest = (n - 1, n - 1, n - 1);

    let mut b = block((n, n, n), 2.0);
    b.release_mode(sharpest, Temperature::celsius(20.0), 1.0);
    let limit = b.max_stable_dt(Time::from_si(0.0));

    // Measured by stepping, not predicted. Sixty steps is far more than enough for a growing
    // mode to be obvious: at three times this dt the amplitude would be past 10⁶.
    run(&mut b, limit, 60);
    let after = b.mode_amplitude(sharpest);
    assert!(
        after.abs() <= 1.0 + 1e-9,
        "sixty steps at the reported limit turned amplitude 1 into {after:.6e}"
    );

    // **Marginal, not conservative**, which is what makes it *the* limit rather than a safe
    // number below it. The sharpest mode should sit just inside the unit circle and flip sign
    // every step — a scheme that damped it comfortably would be leaving stability unused.
    let per_step = block((n, n, n), 2.0).mode_amplification(sharpest, limit);
    assert!(
        (-1.0..-0.5).contains(&per_step),
        "the limit should be marginal and oscillatory: amplification {per_step:.6}"
    );

    // And three times that step — the 1D limit applied to a 3D block — runs away. Predicted
    // rather than stepped, because `step` refuses it, which is the guard doing its job. The
    // point of the pair is that the refusal is placed where divergence actually begins.
    let over = Time::from_si(limit.to_si() * 3.0);
    let runaway = block((n, n, n), 2.0).mode_amplification(sharpest, over);
    assert!(
        runaway.abs() > 1.0 && runaway.abs().powi(60) > 1e6,
        "three times the limit should diverge: amplification {runaway:.6}"
    );
    let mut b = block((n, n, n), 2.0);
    b.step(Time::from_si(0.0), over, &mut Exchange::new())
        .expect_err("and it is refused rather than run");
}
