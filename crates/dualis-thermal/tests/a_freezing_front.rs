//! A freezing front, against Neumann's exact solution of Stefan's problem.
//!
//! The second of `ARCHITECTURE.md`'s depth entries — every domain was single-phase — and like the
//! first it is the second half of a crate that exists rather than a new one. `Substance` carries
//! [`FusionProps`](dualis_core::substance::FusionProps) now and `Solid3D` accounts for it.
//!
//! # Why this problem and not another
//!
//! Because it has an **exact solution**, which almost nothing with a moving boundary does. A
//! semi-infinite liquid at its melting point, its surface dropped to `T_s` below it, freezes to a
//! depth
//!
//! ```text
//!   X(t) = 2λ√(αt)      where  λ e^{λ²} erf(λ) = St/√π,   St = c(T_m − T_s)/L
//! ```
//!
//! That is a position, not a rate and not a limit. So the test is the front's *place* at several
//! times against a transcendental equation solved here, and the temperature everywhere behind it
//! against `erf(x/2√(αt))/erf(λ)`.
//!
//! Measured, for ice under a 20 K undercooling:
//!
//! ```text
//!   dx = 1.0 mm    t = 100 s    5.311 mm against 5.288 mm     0.431%
//!                  t = 400 s   10.570 mm against 10.567 mm    0.032%
//!                  t = 900 s   15.843 mm against 15.848 mm    0.033%
//!   dx = 0.5 mm    t = 100 s    5.285 mm against 5.283 mm     0.032%
//!                  t = 400 s   10.566 mm against 10.565 mm    0.009%
//!                  t = 900 s   15.847 mm against 15.848 mm    0.012%
//! ```
//!
//! The error is **not** monotone in time — it falls by a factor of thirteen between 5 mm and 10 mm of
//! front and then stops — so from the coarse run alone it looks like a floor that refining could not
//! move. It is not one: halving `dx` improves all three, by 13.3×, 3.5× and 2.8×.
//!
//! And nine times the time gives 2.9984× the depth, against `√9 = 3`. That is the `√t` law with no
//! closed form involved at all, which is the claim a scheme with the right coefficient and the wrong
//! power would fail.
//!
//! # The scheme, and the one it is not
//!
//! Enthalpy, bookkept as a temperature and a melted fraction. The state of a cell is one monotone
//! number — `T − T_m` below, `φ·ℓ` inside, `ℓ + T − T_m` above — so the sweep adds energy to it and
//! inverts, and the front is wherever the fraction happens to be between nought and one. Nothing
//! tracks it, nothing iterates, and energy is conserved because energy *is* the state.
//!
//! The method it is not is **apparent heat capacity**, which smears `L` over a narrow temperature
//! interval and multiplies `c_p` up inside it. That method has a failure this one cannot: a step
//! large enough to cross the interval skips the latent heat entirely, and the front then runs fast
//! by up to the factor `L/(cΔT_interval)`. Here the inverse map says where an overshoot's energy
//! goes instead of a branch guessing, which `a_step_that_crosses_the_whole_latent_heat` measures at
//! machine precision.

use dualis_core::conserved::quantity;
use dualis_core::{
    units::{Energy, Length, Temperature, Time, Volume},
    Domain, Exchange, Substance,
};
use dualis_thermal::Solid3D;

/// Ice, from `Substance::ice`. Spelled out here so the closed form is not built from the same call
/// the domain makes.
const K: f64 = 2.22;
const RHO: f64 = 917.0;
const CP: f64 = 2050.0;
const LATENT: f64 = 333_550.0;
/// `k/(ρc)`, m²/s.
const ALPHA: f64 = K / (RHO * CP);

/// `erf`, by its Maclaurin series.
///
/// Written out rather than taken from a crate, and the series rather than a rational fit, because
/// every argument this test needs is under `0.4`: the similarity variable runs from nought at the
/// cold face to `λ ≈ 0.24` at the front, and beyond the front the answer is the melting point by
/// definition. For `|x| ≤ 1` the terms fall like `x^{2n}/n!` and forty of them are past machine
/// precision. It would be wrong for large `x` and is never asked.
fn erf(x: f64) -> f64 {
    assert!(x.abs() <= 1.0, "this series is only used on |x| <= 1: {x}");
    let mut term = x;
    let mut sum = x;
    for n in 1..40 {
        term *= -x * x / n as f64;
        sum += term / (2 * n + 1) as f64;
    }
    sum * 2.0 / std::f64::consts::PI.sqrt()
}

/// Solve `λ e^{λ²} erf(λ) = St/√π` for `λ`, by bisection.
///
/// The left side is zero at zero and increases without bound, so bisection cannot miss it and needs
/// no derivative. Sixty halvings of `[0, 1]` is `1e-18`, which is past what `f64` distinguishes.
fn lambda(stefan: f64) -> f64 {
    let target = stefan / std::f64::consts::PI.sqrt();
    let f = |l: f64| l * (l * l).exp() * erf(l) - target;
    let (mut lo, mut hi) = (0.0, 1.0);
    assert!(f(hi) > 0.0, "the root is not bracketed for St = {stefan}");
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Forty millimetres of ice at one millimetre, and the front reaches sixteen — so the far face is
/// never involved and the column is semi-infinite in the only sense that matters.
const DEPTH: f64 = 40e-3;
const DX: f64 = 1e-3;
/// Depths behind the front, in metres, that the profile is checked at. Cell centres at both
/// resolutions, and all well behind the front at 600 s so none is inside the mush.
const PROBES: [f64; 4] = [2e-3, 4e-3, 6e-3, 8e-3];
/// The column's cell count at [`DX`], for the tests that walk cells by index.
const CELLS_1MM: usize = 40;

/// A column of liquid at exactly its melting point, with the cold face about to be applied.
///
/// `set_melted_fraction` is what makes "liquid at exactly zero" sayable. Starting a hair above the
/// melting point instead would put `ε·c_p` of sensible heat into an initial condition that Neumann's
/// solution does not have, and the front would arrive early by an amount nothing would explain.
fn column(dx: f64) -> Solid3D {
    let cells = (DEPTH / dx).round() as usize;
    let mut c = Solid3D::new(
        "ice",
        Substance::ice(),
        (1, 1, cells),
        Length::from_si(dx),
        Temperature::celsius(0.0),
    );
    for k in 0..cells {
        c.set_melted_fraction(0, 0, k, 1.0);
    }
    c
}

/// How far the solid reaches from the cold face, in metres.
///
/// From the melted *volume*, so the partial cell the front is inside counts for its fraction and the
/// answer is not quantised to `dx`. The half cell is the Dirichlet boundary: clamping cell zero puts
/// the surface at its **centre**, so half of that cell lies before `x = 0` and does not count as
/// frozen depth.
fn frozen_depth(c: &Solid3D, dx: f64) -> f64 {
    let area = dx * dx;
    let solid = DEPTH - c.melted_volume().to_si() / area;
    solid - 0.5 * dx
}

/// March a column with its cold face clamped, stopping at each of `until`, and report the front.
fn march(dx: f64, undercooling: f64, until: &[f64]) -> Vec<(f64, f64)> {
    let cold = Temperature::celsius(-undercooling);
    let mut c = column(dx);
    let dt = Time::from_si(c.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
    let mut out = Vec::new();
    let mut t = 0.0;
    for want_at in until {
        while t < *want_at {
            c.set_temperature(0, 0, 0, cold);
            c.step(Time::from_si(t), dt, &mut Exchange::new())
                .expect("stable");
            t += dt.to_si();
        }
        c.set_temperature(0, 0, 0, cold);
        out.push((t, frozen_depth(&c, dx)));
    }
    out
}

/// **The front is where Neumann's solution puts it, at three times and two resolutions.**
///
/// The load-bearing test, and it checks a *position* against a transcendental closed form rather
/// than a rate against a fit. Three times because `√t` and `t` differ by everything over a factor of
/// nine and by little over a factor of 1.1 — one time would pass for a scheme with the wrong
/// exponent.
///
/// Two cell sizes because that is what says where the remaining error is, and the answer was not
/// what a first draft assumed. The error is **not** monotone in time: it falls from 0.43% to 0.03%
/// between 5 mm and 10 mm of front and then stops, so something other than resolving the front
/// limits it there. Halving `dx` is how to find out which.
#[test]
fn the_front_advances_as_neumanns_solution_says() {
    let undercooling = 20.0;
    let stefan = CP * undercooling / LATENT;
    let l = lambda(stefan);
    let coefficient = 2.0 * l * ALPHA.sqrt();
    println!("  St = {stefan:.5}, lambda = {l:.6}, so X(t) = {coefficient:.4e} sqrt(t) m");

    let times = [100.0, 400.0, 900.0];
    let mut by_dx = Vec::new();
    for dx in [DX, DX / 2.0] {
        let mut errors = Vec::new();
        for (t, got) in march(dx, undercooling, &times) {
            let closed = coefficient * t.sqrt();
            let off = (got / closed - 1.0).abs();
            println!(
                "  dx = {:.1} mm, t = {t:6.1} s   {:.3} mm against {:.3} mm    {:.3}%",
                dx * 1e3,
                got * 1e3,
                closed * 1e3,
                off * 100.0
            );
            // The bound is `(dx/X)²`, which is second order in how well the front is resolved: 3.6%,
            // 0.9% and 0.4% at a millimetre. Per time rather than one blanket number, because a
            // blanket number sized for the shallowest front would say nothing about the deepest.
            let allowed = (dx / closed).powi(2);
            assert!(
                off < allowed,
                "t = {t}: {:.3}% is past the {:.3}% that resolving the front to {:.1} mm costs",
                off * 100.0,
                allowed * 100.0,
                dx * 1e3
            );
            errors.push(off);
        }
        by_dx.push(errors);
    }

    // Halving the cells improves **every** time, which is what says the remaining error is spatial —
    // a first draft asserted the opposite from the coarse run alone, where the error stops falling
    // with *time* and looked like a floor. It is not a floor: 13.3x, 3.5x and 2.8x.
    println!(
        "  halving dx changes the errors by {:.2}x, {:.2}x, {:.2}x",
        by_dx[0][0] / by_dx[1][0],
        by_dx[0][1] / by_dx[1][1],
        by_dx[0][2] / by_dx[1][2]
    );
    for (n, (coarse, fine)) in by_dx[0].iter().zip(&by_dx[1]).enumerate() {
        assert!(
            fine < coarse,
            "time {n} must improve on refinement: {:.4}% -> {:.4}%",
            coarse * 100.0,
            fine * 100.0
        );
    }

    // The exponent, from the two ends of the finer run: X ∝ √t means nine times the time is three
    // times the depth. This is the claim a scheme with the right coefficient and the wrong power
    // would fail, and it needs no closed form at all.
    let fine = march(DX / 2.0, undercooling, &times);
    let ratio = fine[2].1 / fine[0].1;
    println!(
        "  and nine times the time is {ratio:.4}x the depth, against sqrt(9) = 3 — off {:.2e}",
        (ratio / 3.0 - 1.0).abs()
    );
    // The bound is the two front errors this ratio is built from: each is under 0.05% at the finer
    // resolution, so their quotient is under 0.1%. Measured 5.5e-4, which is inside it.
    assert!(
        (ratio / 3.0 - 1.0).abs() < 1e-3,
        "the front goes as sqrt(t): {ratio:.4} against 3"
    );
}

/// **Behind the front the temperature is the error function, and ahead of it there is no gradient.**
///
/// The other half of Neumann's solution, and the half that would catch a scheme moving the front
/// correctly for the wrong reason. `T(x,t) = T_s + (T_m − T_s)·erf(η)/erf(λ)` with `η = x/2√(αt)`.
///
/// Ahead of the front the claim is stronger than "small": the liquid is at the melting point
/// **exactly**, to the last bit, because every one of its faces has zero temperature difference
/// across it and zero times a conductance is zero. That is not an approximation to the one-phase
/// problem; it *is* the one-phase problem, and it is the reason ice is the right catalogue entry to
/// check this with.
///
/// Two resolutions, and the tolerance is the *rate* rather than a value. A first draft asserted
/// 0.1 K, measured 0.158, and 0.1 was a number chosen before anything was run.
#[test]
fn behind_the_front_the_profile_is_an_error_function() {
    let undercooling = 20.0;
    let l = lambda(CP * undercooling / LATENT);
    let melting = Temperature::celsius(0.0).to_si();
    let cold = Temperature::celsius(-undercooling);

    // Averaged over five instants, not read off one. The front advances in a staircase across a
    // fixed grid, so the field at a single moment depends on where the front happens to sit between
    // two cell centres — and that phase is not a smooth function of `dx`. Measured from one snapshot
    // the "order" came out 1.21 and then 3.31, which is the crossing phase and not convergence.
    let sample_times = [300.0, 450.0, 600.0, 750.0, 900.0];
    let mut means = Vec::new();
    for dx in [DX, DX / 2.0, DX / 4.0] {
        let cells = (DEPTH / dx).round() as usize;
        let mut c = column(dx);
        let dt = Time::from_si(c.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
        let mut t = 0.0;
        let mut errors = Vec::new();
        let mut ahead = 0;
        for want_at in sample_times {
            while t < want_at {
                c.set_temperature(0, 0, 0, cold);
                c.step(Time::from_si(t), dt, &mut Exchange::new())
                    .expect("stable");
                t += dt.to_si();
            }
            c.set_temperature(0, 0, 0, cold);
            let front = frozen_depth(&c, dx);
            let scale = 2.0 * (ALPHA * t).sqrt();
            // Probed at **fixed depths**, not at every cell behind the front. "Every cell" measures a
            // different set of points at each resolution — the nearest one to the front moves inward
            // as `dx` shrinks, into exactly the region where a cell-centred profile is least like a
            // continuum one, and that alone read as order 0.66.
            let mut worst: f64 = 0.0;
            for x in PROBES {
                let k = (x / dx).round() as usize;
                assert!(
                    (k as f64 * dx - x).abs() < 1e-12,
                    "{x} m has to be a cell centre at dx = {dx}"
                );
                let got = c.temperature_at(0, 0, k).to_si();
                let want = melting - undercooling + undercooling * erf(x / scale) / erf(l);
                worst = worst.max((got - want).abs());
            }
            errors.push(worst);

            // Ahead of the front, exact rather than close. A single bit here would mean heat had
            // crossed a face with no temperature difference across it.
            ahead = 0;
            for k in 1..cells {
                if k as f64 * dx > front + dx {
                    let got = c.temperature_at(0, 0, k).to_si();
                    assert!(
                        got == melting,
                        "ahead of the front nothing has happened at all: cell {k} is {got:.17e}, \
                         not {melting:.17e}"
                    );
                    assert!(
                        c.melted_fraction_at(0, 0, k) == 1.0,
                        "and it is still entirely liquid: cell {k} is {}",
                        c.melted_fraction_at(0, 0, k)
                    );
                    ahead += 1;
                }
            }
        }
        let mean = errors.iter().sum::<f64>() / errors.len() as f64;
        println!(
            "  dx = {:.2} mm: mean worst error {mean:.4} K over five instants ({:.4} to {:.4}), and \
             {ahead} cells still untouched to the last bit",
            dx * 1e3,
            errors.iter().copied().fold(f64::MAX, f64::min),
            errors.iter().copied().fold(0.0, f64::max)
        );
        assert!(ahead > 10, "there should be a liquid region: {ahead} cells");
        means.push(mean);
    }

    // **The claim is the fourfold refinement, not the consecutive ratios, and that is a finding.**
    //
    // Between neighbouring resolutions the order comes out 0.78 and then 1.84, and averaging five
    // instants did not settle it. The reason is that a fixed grid makes the front advance in a
    // **staircase**: the field at any moment depends on where the front sits between two cell
    // centres, and that phase is not a smooth function of `dx`. An order read from two adjacent
    // resolutions here is reading that phase.
    //
    // Over the full fourfold refinement it is a factor of 6.2 — between first and second order, and
    // first order is what a cell-centred enthalpy method should give: the cell the front is inside
    // holds the melting point across its whole width where the continuum has a gradient, so the
    // boundary condition is misplaced by `O(dx)` and the temperature field inherits it.
    //
    // The front *position* does better than the field it comes out of, and the reason is that
    // `melted_volume` is an integral of a conserved quantity while a temperature is a point sample.
    // A first draft asserted second order here on no evidence, and then first order on five instants
    // that were not enough to earn it.
    let orders: Vec<f64> = means.windows(2).map(|w| (w[0] / w[1]).log2()).collect();
    let overall = (means[0] / means[2]).log2() / 2.0;
    println!(
        "  {:.4} K, {:.4} K, {:.4} K — adjacent orders {:.2} and {:.2}, overall {overall:.2}",
        means[0], means[1], means[2], orders[0], orders[1]
    );
    for pair in means.windows(2) {
        assert!(
            pair[1] < pair[0],
            "every refinement must improve it: {means:?}"
        );
    }
    assert!(
        overall > 1.0 && overall < 2.0,
        "four times the cells is between first and second order: {means:?} gives {overall:.3}"
    );
    assert!(
        means[0] < 0.01 * undercooling,
        "and the coarse run is already under a percent of the undercooling: {:.4} K of \
         {undercooling}",
        means[0]
    );
}

/// **The plateau is exactly `Λ/P` long, and it is the whole of the latent heat.**
///
/// The simplest closed form phase change has, and the one that would catch a latent heat that was
/// merely *large* rather than right: a cell at its melting point fed constant power holds there for
/// `ρLV/P` seconds and not a moment more.
///
/// For a cubic millimetre of ice that is 306 mJ, against 1.88 mJ per kelvin of warming — so the
/// plateau is worth **163 K**, and a scheme that dropped the latent heat would not be slightly wrong.
#[test]
fn a_cell_at_its_melting_point_holds_there_for_exactly_the_latent_heat() {
    let mut c = Solid3D::new(
        "cube",
        Substance::ice(),
        (1, 1, 1),
        Length::from_si(DX),
        Temperature::celsius(0.0),
    );
    let volume = Volume::from_si(DX * DX * DX);
    let latent = Substance::ice()
        .latent_energy(volume)
        .expect("ice melts")
        .to_si();
    let capacity = Substance::ice()
        .heat_capacity(volume)
        .expect("ice has a specific heat")
        .to_si();
    println!(
        "  {:.1} mJ of latent heat against {:.3} mJ/K of capacity — a plateau worth {:.1} K",
        latent * 1e3,
        capacity * 1e3,
        latent / capacity
    );
    assert!(
        (latent / capacity - LATENT / CP).abs() < 1e-9,
        "the plateau is L/c_p: {:.4} against {:.4}",
        latent / capacity,
        LATENT / CP
    );

    // A milliwatt, delivered in a thousand equal parcels. The plateau should be 306 s.
    let power = 1e-3;
    let dt = 0.1;
    let want = latent / power;
    let mut t = 0.0;
    let mut left_at = None;
    let melting = Temperature::celsius(0.0).to_si();
    while t < 2.0 * want {
        c.deposit(0, 0, 0, Energy::from_si(power * dt));
        t += dt;
        if left_at.is_none() && c.temperature_at(0, 0, 0).to_si() > melting {
            left_at = Some(t);
        }
    }
    let left = left_at.expect("it must warm up eventually");
    println!("  it left the melting point at {left:.1} s, against L/P = {want:.1} s");
    // Within one parcel, which is the resolution the delivery has and nothing to do with the scheme.
    assert!(
        (left - want).abs() <= dt * 1.001,
        "the plateau is the latent heat over the power: {left:.3} s against {want:.3} s"
    );
    // And it did hold: no temperature rise at all until then.
    assert!(
        c.temperature_at(0, 0, 0).to_si() > melting,
        "past the plateau it warms"
    );
}

/// **A single step that crosses the whole latent heat loses none of it.**
///
/// The failure the apparent-heat-capacity method has and this one does not. That method spreads `L`
/// over a temperature interval and raises `c_p` inside it, so a step big enough to jump the interval
/// pays only the sensible part and the rest is gone — silently, and in the direction that makes a
/// front run fast.
///
/// Here the state is a monotone function of energy, so an overshoot is not a case to handle: the
/// inverse map puts the remainder on the far side because that is where that much energy is. Asserted
/// at machine precision with **ten times** the latent heat in one delivery, which is a hundred and
/// sixty times the interval any apparent capacity would use.
#[test]
fn a_step_that_crosses_the_whole_latent_heat_deposits_the_remainder() {
    let volume = Volume::from_si(DX * DX * DX);
    let latent = Substance::ice()
        .latent_energy(volume)
        .expect("ice melts")
        .to_si();
    let capacity = Substance::ice()
        .heat_capacity(volume)
        .expect("ice has c_p")
        .to_si();
    let melting = Temperature::celsius(0.0).to_si();

    for multiple in [0.5, 1.0, 1.0000001, 3.0, 10.0] {
        let mut c = Solid3D::new(
            "cube",
            Substance::ice(),
            (1, 1, 1),
            Length::from_si(DX),
            Temperature::celsius(0.0),
        );
        let joules = multiple * latent;
        c.deposit(0, 0, 0, Energy::from_si(joules));

        // Where that energy must be: all of it melting until the cell is liquid, the rest warming it.
        let (want_phi, want_t) = if multiple >= 1.0 {
            (1.0, melting + (joules - latent) / capacity)
        } else {
            (multiple, melting)
        };
        let phi = c.melted_fraction_at(0, 0, 0);
        let temp = c.temperature_at(0, 0, 0).to_si();
        println!(
            "  {multiple:>9} x L: melted {phi:.9}, at {:.6} C — want {want_phi:.9} and {:.6} C",
            temp - melting,
            want_t - melting
        );
        assert!(
            (phi - want_phi).abs() < 1e-15,
            "{multiple} x L melts {want_phi}, got {phi}"
        );
        assert!(
            (temp - want_t).abs() < 1e-12 * (1.0 + (want_t - melting).abs()),
            "{multiple} x L leaves it at {want_t:.9}, got {temp:.9}"
        );
        // And the ledger says exactly what went in, which is the statement that nothing was skipped.
        let held = c.ledger().get(quantity::ENERGY).expect("energy is booked");
        assert!(
            (held / joules - 1.0).abs() < 1e-14,
            "{multiple} x L: {held:.9e} J on the books against {joules:.9e} delivered"
        );
    }
}

/// **A freezing column's books balance, and the latent heat is most of what is in them.**
///
/// A conservation check that could not exist before: the block loses heat through a clamped face, so
/// `absorbed_energy` does not see it, and the audit is against what the *state* says it holds. If the
/// ledger read temperature alone it would miss the 306 mJ per cubic millimetre the front has released
/// and call the run a leak of exactly that.
///
/// The split is the point. At 600 s the front is 13 mm down and the latent heat released is several
/// times the sensible heat, so a ledger that dropped it would be wrong by more than a factor of two
/// rather than by a correction.
///
/// # The comparison is a *change*, and getting that wrong cost a run
///
/// `set_melted_fraction` declaring the column liquid is an **opening balance**, exactly as
/// `set_temperature` is: it moves what the block holds and not what it has absorbed, so the ledger
/// opens at `Σ C·ℓ` — 12.235 J for forty cells of ice. A first draft compared the closing balance
/// against the state measured from the melting point and was out by precisely that, forty cells'
/// worth, which is what said the two were measuring from different places rather than disagreeing.
#[test]
fn the_ledger_counts_the_latent_heat_and_not_only_the_cooling() {
    let undercooling = 20.0;
    let cold = Temperature::celsius(-undercooling);
    let mut c = column(DX);
    let opening = c.ledger().get(quantity::ENERGY).expect("energy is booked");
    let dt = Time::from_si(c.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
    let mut t = 0.0;
    while t < 600.0 {
        c.set_temperature(0, 0, 0, cold);
        c.step(Time::from_si(t), dt, &mut Exchange::new())
            .expect("stable");
        t += dt.to_si();
    }
    c.set_temperature(0, 0, 0, cold);

    // The **change** in the books, against the change in state. Both sides start from the column as
    // it was declared: liquid, everywhere, at exactly the melting point.
    let held = c.ledger().get(quantity::ENERGY).expect("energy is booked") - opening;
    // What the state says, summed cell by cell from the two things a cell can hold.
    let cell_capacity = RHO * CP * DX * DX * DX;
    let cell_latent = RHO * LATENT * DX * DX * DX;
    let melting = Temperature::celsius(0.0).to_si();
    let mut sensible = 0.0;
    let mut fusion = 0.0;
    for k in 0..CELLS_1MM {
        sensible += cell_capacity * (c.temperature_at(0, 0, k).to_si() - melting);
        fusion += cell_latent * (c.melted_fraction_at(0, 0, k) - 1.0);
    }
    println!(
        "  {:.4} J on the books: {:.4} J of cooling and {:.4} J of freezing",
        held, sensible, fusion
    );
    assert!(
        (held - (sensible + fusion)).abs() < 1e-12 * fusion.abs(),
        "the ledger is the state: {held:.9} against {:.9}",
        sensible + fusion
    );
    println!(
        "  the latent heat is {:.1}x the sensible, so a ledger that dropped it would be wrong by          more than a factor of two",
        fusion.abs() / sensible.abs()
    );
    assert!(
        fusion.abs() > 2.0 * sensible.abs(),
        "the latent heat should dominate, or this is not a phase-change problem: {:.4} against \
         {:.4} J",
        fusion.abs(),
        sensible.abs()
    );
}

/// **A checkpoint carries the phase, because a temperature is not a state.**
///
/// The silent failure this guards, and it is a live path: `Schedule::Iterative` restores every
/// iteration and the conservation audit restores on a violation. A checkpoint of temperatures alone
/// would bring a half-frozen column back as a fully liquid one at the same temperatures, losing the
/// entire latent heat released so far with nothing to say it had gone — and the ledger would then
/// balance, because the ledger reads the same state that was corrupted.
#[test]
fn a_checkpoint_carries_the_phase_and_not_only_the_temperature() {
    let cold = Temperature::celsius(-20.0);
    let mut c = column(DX);
    let dt = Time::from_si(c.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
    let mut t = 0.0;
    while t < 120.0 {
        c.set_temperature(0, 0, 0, cold);
        c.step(Time::from_si(t), dt, &mut Exchange::new())
            .expect("stable");
        t += dt.to_si();
    }
    c.checkpoint();
    let (depth, held) = (
        frozen_depth(&c, DX),
        c.ledger().get(quantity::ENERGY).expect("booked"),
    );
    assert!(
        depth > 3e-3,
        "there should be a front to lose: {depth:.4} m"
    );

    // Carry on, then go back.
    for n in 0..200 {
        c.set_temperature(0, 0, 0, cold);
        c.step(
            Time::from_si(t + n as f64 * dt.to_si()),
            dt,
            &mut Exchange::new(),
        )
        .expect("stable");
    }
    assert!(
        frozen_depth(&c, DX) > depth,
        "it should have advanced before being restored"
    );
    c.restore();

    let back = frozen_depth(&c, DX);
    println!(
        "  the front was at {:.3} mm, advanced, and came back to {:.3} mm",
        depth * 1e3,
        back * 1e3
    );
    assert!(
        back == depth,
        "a restore returns the phase exactly: {back:.17e} against {depth:.17e}"
    );
    assert!(
        c.ledger().get(quantity::ENERGY).expect("booked") == held,
        "and therefore the books"
    );
}

/// **A substance that does not melt gives the block it gave before, bit for bit.**
///
/// The reduction, and it has to be exact for the same reason the multi-material one did: the phase
/// machinery is on the same sweep every existing closed form is checked through, so if it moved a
/// non-melting answer then every one of those was quietly reinterpreted.
///
/// Aluminium *does* melt, at 855 K, and this is the honest statement of what absence means:
/// `Substance::aluminium_6061` carries no `FusionProps`, so the domain models it as never changing
/// phase — which is correct for a heat sink and would be wrong for a casting.
#[test]
fn a_block_of_something_that_does_not_melt_is_unchanged() {
    let build = || {
        let mut b = Solid3D::new(
            "b",
            Substance::aluminium_6061(),
            (4, 3, 5),
            Length::mm(2.0),
            Temperature::celsius(20.0),
        );
        b.deposit(2, 1, 2, Energy::from_si(3.0));
        b
    };
    let mut plain = build();
    assert_eq!(plain.melted_fraction_at(2, 1, 2), 0.0, "nothing melts");
    assert_eq!(
        plain.melted_volume().to_si(),
        0.0,
        "and nothing is liquid to begin with"
    );

    // `set_melted_fraction` on a substance with no fusion is ignored rather than half-applied: there
    // is no state for it to mean, and inventing one would make a block that melts by being asked to.
    plain.set_melted_fraction(2, 1, 2, 1.0);
    assert_eq!(plain.melted_fraction_at(2, 1, 2), 0.0);

    let mut reference = build();
    let dt = Time::from_si(reference.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    for n in 0..60 {
        let t = Time::from_si(n as f64 * dt.to_si());
        plain.step(t, dt, &mut Exchange::new()).expect("stable");
        reference.step(t, dt, &mut Exchange::new()).expect("stable");
    }
    let (a, b) = (
        plain.peak_temperature().to_si(),
        reference.peak_temperature().to_si(),
    );
    assert!(a - b == 0.0, "{a:.17e} against {b:.17e}");
    assert_eq!(plain.melted_volume().to_si(), 0.0, "and still none");
}
