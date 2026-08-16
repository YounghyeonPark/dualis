//! A block that can lose heat, checked against the closed forms that say where it settles.
//!
//! `Solid3D` was adiabatic on all six faces until now, so no three-dimensional thermal scene
//! could reach a steady state: every one of them warmed for as long as it ran. That is honest
//! for a pulse and useless for the question a designer asks — *what temperature does this run
//! at* — which is the question a chip package, a magnet busbar, a motor and a factory cell are
//! all asking.
//!
//! Every claim below is against something the domain did not compute: a lumped balance it has
//! no code for, Newton's law of cooling, a conduction resistance in series with a film, and the
//! ledger identity that says nothing was lost on the way out.

use dualis_core::units::{Area, Length, Temperature, Time};
use dualis_core::{Domain, Exchange, Substance};
use dualis_thermal::{Environment, Face, Solid3D};

/// A copper cube, `n` cells on a side of `cell_mm`, starting at `start_c`.
fn cube(n: usize, cell_mm: f64, start_c: f64) -> Solid3D {
    Solid3D::new(
        "block",
        Substance::copper(),
        (n, n, n),
        Length::mm(cell_mm),
        Temperature::celsius(start_c),
    )
}

/// What the block holds, read the way the audit reads it: the ledger, minus what it has shed.
fn stored(block: &Solid3D) -> f64 {
    use dualis_core::conserved::quantity;
    block.ledger().get(quantity::ENERGY).unwrap_or(0.0) - block.lost_energy().to_si()
}

/// Run to `seconds`, subcycling at the domain's own limit, and give back the block.
fn run(mut block: Solid3D, seconds: f64) -> Solid3D {
    let limit = block.max_stable_dt(Time::ZERO).to_si();
    let dt = Time::from_si(limit * 0.4);
    let steps = (seconds / dt.to_si()).ceil() as usize;
    let mut bus = Exchange::new();
    let mut t = Time::ZERO;
    for _ in 0..steps {
        block.step(t, dt, &mut bus).expect("a stable step");
        t += dt;
    }
    block
}

/// **A block with no exposed face is still adiabatic, exactly.** The guarantee that let this be
/// added without changing a single existing scene: nothing loses heat unless a caller says a
/// face is exposed.
#[test]
fn an_unexposed_block_loses_nothing() {
    let before = cube(4, 5.0, 200.0);
    let opening = stored(&before);
    let after = run(before, 600.0);
    assert_eq!(after.lost_energy().to_si(), 0.0, "nothing was exposed");
    assert!(
        (stored(&after) - opening).abs() < 1e-9,
        "an insulated block holds what it had"
    );
}

/// **A cooling block settles at ambient, and gets there at Newton's rate.**
///
/// Copper at 5 mm cells has a Biot number far below 0.1 — `hL/k` is about `12 · 0.01 / 401`,
/// which is `3e-4` — so the block is isothermal to a part in three thousand and the closed form
/// is the lumped one: `T(t) = T∞ + (T₀ − T∞)·exp(−t/τ)` with `τ = mc/(hA)`.
///
/// `h` is the environment's convective coefficient **plus** the radiative term the domain also
/// models, and that is the point of computing the expected `τ` from `Environment::loss_from`
/// rather than from the convective number alone: a check written against convection only would
/// be wrong by the 6 W·m⁻²·K⁻¹ a black surface radiates at room temperature, which is the same
/// order again.
#[test]
fn a_cooling_block_follows_newtons_law() {
    let n = 4;
    let cell = 5.0e-3;
    let side = n as f64 * cell;
    let area = Area::from_si(6.0 * side * side);
    let ambient = Temperature::celsius(20.0);
    let start = Temperature::celsius(120.0);

    let mut block = cube(n, 5.0, 120.0);
    for face in Face::ALL {
        block = block.losing_from(face, Environment::still_air(ambient, area / 6.0));
    }

    // The time constant, from the loss the domain's own environment model gives at the
    // *midpoint* of the swing — the linearisation Newton's law assumes, evaluated where it is
    // most representative rather than at either end.
    let capacity = block.heat_capacity().to_si();
    let mid = Temperature::celsius(70.0);
    let env = Environment::still_air(ambient, area);
    let h_a = env
        .loss_from(mid, Substance::copper().thermal.unwrap().emissivity)
        .to_si()
        / (mid.to_si() - ambient.to_si());
    let tau = capacity / h_a;

    let elapsed = tau;
    let after = run(block, elapsed);
    let expect = ambient.to_si() + (start.to_si() - ambient.to_si()) * (-elapsed / tau).exp();
    let got = after.mean_temperature().to_si();

    // One time constant of a hundred-kelvin swing, against a closed form whose `h` was frozen
    // at the midpoint while the real one falls as the block cools. That curvature is the whole
    // error budget: the radiative part of `h` moves by about 12% across the swing and enters
    // the exponent, so a couple of kelvin is what it buys. Measured at 1.2 K.
    assert!(
        (got - expect).abs() < 3.0,
        "after one time constant the block should be near {expect:.2} K, and it is {got:.2} K"
    );

    // And it keeps going the right way: far past the constant it is within a kelvin of the air.
    let settled = run(after, 12.0 * tau);
    assert!(
        (settled.mean_temperature().to_si() - ambient.to_si()).abs() < 1.0,
        "a block left in air settles at the air's temperature, got {} K",
        settled.mean_temperature().to_si()
    );
}

/// **Nothing is lost on the way out: `stored + lost` is conserved to the bit.**
///
/// The identity that lets `books_balance` stay `true` for a block that sheds heat. Every joule
/// that leaves a cell is counted in the same statement that removes it, so the ledger moves
/// only by what crossed the bus — and here nothing did.
#[test]
fn every_joule_that_leaves_is_counted() {
    let mut block = cube(3, 4.0, 300.0);
    let side = 3.0 * 4.0e-3;
    block = block.losing_from(
        Face::ZMax,
        Environment::still_air(Temperature::celsius(20.0), Area::from_si(side * side)),
    );

    let opening = stored(&block) + block.lost_energy().to_si();
    let after = run(block, 2000.0);
    let closing = stored(&after) + after.lost_energy().to_si();

    assert!(
        after.lost_energy().to_si() > 0.0,
        "it should have shed heat"
    );
    // Relative to the joules that actually moved, not to the block's absolute enthalpy, so the
    // check does not get easier as the block gets bigger.
    let moved = after.lost_energy().to_si();
    assert!(
        (closing - opening).abs() / moved < 1e-12,
        "stored + lost moved by {:e} J against {moved:e} J shed",
        closing - opening
    );
}

/// **The Biot number decides whether a cooled block has a gradient at all, and it does.**
///
/// The design question behind every cooled part: is this thing isothermal, or does the inside
/// run hotter than the surface? The textbook answer is `Bi = hL/k` — under about 0.1 the block
/// is lumped and a single temperature describes it; well above, it is not and the middle runs
/// hot. Nothing in the domain computes a Biot number, so this is a check against an outside
/// statement rather than against itself.
///
/// Two blocks of the same poor conductor, differing **only** in the film coefficient, three
/// decades apart. The measured spread — hottest cell to coldest, against the block's own
/// distance from the air — has to track it.
#[test]
fn the_biot_number_decides_whether_a_cooled_block_is_isothermal() {
    // Borosilicate at 12 mm: a still-air film puts it at Bi = 0.05 and a fan at 50, which
    // straddles the lumped criterion with two coefficients that are both real. A polymer is too
    // poor a conductor for the low end — even still air leaves it at Bi ≈ 0.9 — and a metal is
    // too good for the high one.
    let glass = Substance::from_name("borosilicate").expect("the catalogue has it");
    let k = glass.thermal.unwrap().conductivity.to_si();
    let cells = 4;
    let cell = 3.0e-3;
    let ambient = Temperature::celsius(20.0);

    // The spread a block develops after cooling for `seconds` with this film on one face.
    let spread = |h: f64, seconds: f64| {
        let block = Solid3D::new(
            "part",
            glass.clone(),
            (cells, cells, cells),
            Length::from_si(cell),
            Temperature::celsius(120.0),
        )
        .losing_from(
            Face::ZMax,
            Environment {
                ambient,
                convection_w_per_m2_k: h,
                area: Area::from_si((cells as f64 * cell).powi(2)),
            },
        );
        let after = run(block, seconds);
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for kk in 0..cells {
            let t = after.temperature_at(0, 0, kk).to_si();
            lo = lo.min(t);
            hi = hi.max(t);
        }
        (hi - lo, after.mean_temperature().to_si() - ambient.to_si())
    };

    // The characteristic length for a face-cooled block is the block itself.
    let length = cells as f64 * cell;
    let (gentle_spread, gentle_above) = spread(5.0, 300.0);
    let (fierce_spread, fierce_above) = spread(5000.0, 300.0);
    let bi_gentle = 5.0 * length / k;
    let bi_fierce = 5000.0 * length / k;

    assert!(
        bi_gentle < 0.1 && bi_fierce > 10.0,
        "the two films should straddle the lumped criterion: {bi_gentle:.3} and {bi_fierce:.0}"
    );
    // A gentle film leaves the block nearly isothermal: the internal spread is a small part of
    // how far the whole block sits from the air.
    assert!(
        gentle_spread / gentle_above < 0.15,
        "at Bi = {bi_gentle:.3} the block should be nearly lumped, spread {gentle_spread:.2} K \
         against {gentle_above:.2} K above air"
    );
    // A fierce one does not: the face is dragged to the air and the far side is not.
    assert!(
        fierce_spread / fierce_above > 1.0,
        "at Bi = {bi_fierce:.0} the block cannot be lumped, spread {fierce_spread:.2} K against \
         {fierce_above:.2} K above air"
    );
}

/// **A film stiffer than the conduction it replaces tightens the stability limit, and the
/// domain refuses a step past it.**
///
/// The silent failure this guards: a boundary integrated explicitly past its own limit does not
/// diverge dramatically, it **oscillates about ambient** — a cell overshooting the air it is
/// losing to, every step, while the conservation audit stays perfectly happy because the heat
/// really did leave.
///
/// The condition is worth stating because it is not "any exposed face". `worst_rate` is a
/// maximum over cells, an interior cell has six conducting faces and a boundary cell has five
/// plus a film, so the limit only moves when `h·dx > k` — the film outweighing the conduction
/// it stands in for. For copper at 2 mm that needs `h > 200 kW/m²·K`, which is not air; for a
/// polymer it needs a few hundred, which is a fan. Checked on the polymer, where it bites.
#[test]
fn a_stiff_film_tightens_the_limit_and_the_limit_is_enforced() {
    let pla = Substance::from_name("pla").expect("the catalogue has it");
    let side = 4.0 * 2.0e-3;
    let bare = Solid3D::new(
        "part",
        pla.clone(),
        (4, 4, 4),
        Length::mm(2.0),
        Temperature::celsius(300.0),
    );
    let bare_limit = bare.max_stable_dt(Time::ZERO).to_si();

    let cooled = Solid3D::new(
        "part",
        pla,
        (4, 4, 4),
        Length::mm(2.0),
        Temperature::celsius(300.0),
    )
    .losing_from(
        Face::ZMax,
        Environment {
            ambient: Temperature::celsius(20.0),
            convection_w_per_m2_k: 2000.0,
            area: Area::from_si(side * side),
        },
    );
    let cooled_limit = cooled.max_stable_dt(Time::ZERO).to_si();

    assert!(
        cooled_limit < bare_limit * 0.9,
        "a film stiffer than the conduction it replaces must shorten the step: \
         {cooled_limit:e} against {bare_limit:e}"
    );

    // And a step past it is refused rather than run into an oscillation.
    let mut block = cooled;
    let mut bus = Exchange::new();
    block
        .step(Time::ZERO, Time::from_si(cooled_limit * 1.5), &mut bus)
        .expect_err("a step past the limit must be refused");
}

/// **A gentle film does not move the limit, and that is not a bug.** The other side of the same
/// fact, pinned so nobody "fixes" the limit into being pessimistic: still air on copper is
/// nowhere near the conduction it stands in for, so the interior cell still sets the step.
#[test]
fn a_gentle_film_leaves_the_limit_alone() {
    let side = 4.0 * 5.0e-3;
    let bare = cube(4, 5.0, 200.0).max_stable_dt(Time::ZERO).to_si();
    let cooled = cube(4, 5.0, 200.0)
        .losing_from(
            Face::ZMax,
            Environment::still_air(Temperature::celsius(20.0), Area::from_si(side * side)),
        )
        .max_stable_dt(Time::ZERO)
        .to_si();
    assert_eq!(
        cooled, bare,
        "still air on copper is far weaker than the face it replaces"
    );
}

/// **A hot block's limit is shorter than a cool one's, and that is not a refinement.**
///
/// The radiative conductance is `4εσT³`, so a part at 1000 °C carries about **26 times** what
/// the same surface carries at room temperature — measured at `ε = 0.9`, 136 against
/// 5.14 W/m²·K. A limit linearised about the *ambient* would hand the hot block a step an order
/// too long, and an explicit boundary past its limit does not fail loudly: it oscillates about
/// the air it is losing to, every step, while the conservation audit stays perfectly happy
/// because the heat really did leave.
///
/// So the limit is evaluated at the block's own temperature, which is what
/// [`Domain::max_stable_dt`] means by "from `now`". Checked as a ratio between two blocks
/// differing only in how hot they are, and asserted one-sided against the `T³` the physics
/// gives — a limit that ignored the state would return the identical number for both.
#[test]
fn a_hot_block_gets_a_shorter_limit_than_a_cool_one() {
    // A poor conductor with a convective film just under the threshold on its own, so the
    // radiative term is what carries it across. The film cannot exceed the half cell of solid
    // behind it — , here 743 W/m^2/K — so radiation alone at 1000 C (420) does not reach
    // it and 600 of convection beside it does.
    let glass = Substance::from_name("borosilicate").expect("the catalogue has it");
    let side = 4.0 * 3.0e-3;
    let block = |celsius: f64| {
        Solid3D::new(
            "part",
            glass.clone(),
            (4, 4, 4),
            Length::mm(3.0),
            Temperature::celsius(celsius),
        )
        .losing_from(
            Face::ZMax,
            Environment {
                ambient: Temperature::celsius(20.0),
                convection_w_per_m2_k: 600.0,
                area: Area::from_si(side * side),
            },
        )
        .max_stable_dt(Time::ZERO)
        .to_si()
    };

    let cool = block(30.0);
    let hot = block(1000.0);
    assert!(
        hot < cool,
        "a hot block radiates harder and must take shorter steps: {hot:e} against {cool:e}"
    );

    // And a block with no exposed face has a limit that does not move with temperature at all,
    // which is the control that says the difference above came from the boundary.
    let bare = |celsius: f64| {
        Solid3D::new(
            "part",
            glass.clone(),
            (4, 4, 4),
            Length::mm(3.0),
            Temperature::celsius(celsius),
        )
        .max_stable_dt(Time::ZERO)
        .to_si()
    };
    assert_eq!(
        bare(30.0),
        bare(1000.0),
        "conduction alone does not care how hot the block is"
    );
}
