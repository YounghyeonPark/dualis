//! A superheated liquid slows the front, against the two-phase Neumann solution.
//!
//! `a_freezing_front.rs` is the **one-phase** problem: the liquid sits at the melting point, so no
//! heat flows through it and its properties cannot matter. That is Stefan's original problem and it is
//! exact.
//!
//! It stops being the right model the moment the liquid is warmer than the melting point, because then
//! the liquid conducts heat *toward* the interface and the front has to remove that as well as the
//! latent heat. Water conducts a quarter of what ice does and holds twice as much, so the correction is
//! not small: **20 K of superheat slows the front 16%**, from 15.85 mm to 13.33 mm at 900 s.
//!
//! # The closed form, and the check on the algebra
//!
//! The energy balance at `x = X(t) = 2λ√(α_s t)` is `k_s ∂T_s/∂x − k_l ∂T_l/∂x = ρL dX/dt`, with an
//! `erf` profile in the solid and an `erfc` one in the liquid. Working it through:
//!
//! ```text
//!   exp(−λ²)/erf(λ) − (k_l/k_s)·ν·((T_i−T_m)/(T_m−T_s))·exp(−ν²λ²)/erfc(νλ) = λ√π/St
//!
//!   ν = √(α_s/α_l),   St = c_s(T_m − T_s)/L
//! ```
//!
//! Setting `T_i = T_m` kills the middle term and leaves `λ e^{λ²} erf(λ) = St/√π`, which is the
//! one-phase condition the other file uses. **That reduction is the check on the algebra**, and it is
//! asserted here rather than trusted: the two-phase `λ` converges on the one-phase one linearly in the
//! superheat, at `9.45e-3` per kelvin.
//!
//! # What it costs, and why two phases are not the default
//!
//! The cell holding the interface is part solid and part liquid, so its conductivity is a mixture and
//! the heat reaching the interface sees less conductance than the exact solution gives it. That is a
//! **first-order** error at the front, and it is why `Substance::ice` carries `liquid: None`: switching
//! it on took the *one-phase* answer from 0.43% out to 6.9% at forty cells, for a problem whose physics
//! had not changed. A default that costs a caller a factor of sixteen is the wrong default.
//!
//! So the bound here is `dx/X` and not `(dx/X)²`, and it is stated as the mush's price rather than
//! discovered as a disappointment.

use pantometry_core::substance::ThermalProps;
use pantometry_core::{
    units::{Length, Temperature, Time},
    Domain, Exchange, Substance,
};
use pantometry_thermal::Solid3D;

const K_S: f64 = 2.22;
const RHO: f64 = 917.0;
const C_S: f64 = 2050.0;
const LATENT: f64 = 333_550.0;
const ALPHA_S: f64 = K_S / (RHO * C_S);
/// Water, at the same density as the ice — a density change at the front is a separate term and this
/// is the standard treatment without it.
const K_L: f64 = 0.598;
const C_L: f64 = 4182.0;
const ALPHA_L: f64 = K_L / (RHO * C_L);

const DEPTH: f64 = 40e-3;
const DROP: f64 = 20.0;

/// Ice with water named as its liquid phase, which is what makes a block two-phase.
fn two_phase_ice() -> Substance {
    let base = Substance::ice();
    let fusion = base.fusion.expect("ice melts");
    base.with_fusion(fusion.with_liquid(ThermalProps {
        conductivity: pantometry_core::units::ThermalConductivity::w_per_m_k(K_L),
        specific_heat: pantometry_core::units::SpecificHeat::j_per_kg_k(C_L),
        expansion: pantometry_core::units::ThermalExpansion::ppm_per_k(207.0),
        emissivity: 0.96,
    }))
}

fn erf(x: f64) -> f64 {
    // The Maclaurin series, valid for the arguments this file uses. `erfc(νλ)` reaches νλ ≈ 0.67 at
    // the largest ν here, so nothing is asked beyond 1.
    assert!(x.abs() <= 3.0, "the series loses precision past 3: {x}");
    let mut term = x;
    let mut sum = x;
    for n in 1..60 {
        term *= -x * x / n as f64;
        sum += term / (2 * n + 1) as f64;
    }
    sum * 2.0 / std::f64::consts::PI.sqrt()
}

fn erfc(x: f64) -> f64 {
    1.0 - erf(x)
}

/// The largest step an all-**solid** block of this material allows, which is the tightest state a
/// freezing run reaches.
fn solid_limit() -> f64 {
    let c = Solid3D::new(
        "probe",
        Substance::ice(),
        (1, 1, 8),
        Length::from_si(1e-3),
        Temperature::celsius(-1.0),
    );
    c.max_stable_dt(Time::from_si(0.0)).to_si()
}

/// `λ` for a solid growing into a liquid `superheat` kelvin above the melting point.
fn two_phase_lambda(superheat: f64) -> f64 {
    let nu = (ALPHA_S / ALPHA_L).sqrt();
    let stefan = C_S * DROP / LATENT;
    let rhs = std::f64::consts::PI.sqrt() / stefan;
    let f = |lam: f64| {
        let first = (-lam * lam).exp() / erf(lam);
        let second =
            (K_L / K_S) * nu * (superheat / DROP) * (-nu * nu * lam * lam).exp() / erfc(nu * lam);
        first - second - rhs * lam
    };
    // `hi` is 0.45 rather than 1: lambda is under 0.3 for every case here, and a larger bracket
    // would ask `erf` for `nu*hi` = 2.75, where the series is still convergent but no longer sharp.
    let (mut lo, mut hi) = (1e-9, 0.45);
    assert!(f(hi) < 0.0, "the root is bracketed");
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// The one-phase `λ`, from `λ e^{λ²} erf(λ) = St/√π`.
fn one_phase_lambda() -> f64 {
    let target = C_S * DROP / LATENT / std::f64::consts::PI.sqrt();
    let g = |lam: f64| lam * (lam * lam).exp() * erf(lam) - target;
    let (mut lo, mut hi) = (0.0, 0.45);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if g(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// **The two-phase condition reduces to the one-phase one, linearly in the superheat.**
///
/// The check on the algebra, and it is available before any solver runs. If the transcendental were
/// wrong, everything measured against it would be wrong in the same direction and nothing downstream
/// would say so.
///
/// Linear convergence rather than mere agreement: the extra term is proportional to the superheat, so
/// the gap should be too, and it is — `9.45e-12`, `9.45e-9`, `9.45e-6` at `1e-9`, `1e-6` and `1e-3`
/// kelvin. Three decades of exactly proportional error is a much stronger statement than one equality.
#[test]
fn the_two_phase_condition_reduces_to_the_one_phase_one() {
    let one = one_phase_lambda();
    println!("  one-phase lambda {one:.9}");
    let mut previous: Option<(f64, f64)> = None;
    for superheat in [1e-9, 1e-6, 1e-3] {
        let two = two_phase_lambda(superheat);
        let gap = (two / one - 1.0).abs();
        println!("  superheat {superheat:<8} lambda {two:.9}  gap {gap:.3e}");
        if let Some((prev_sh, prev_gap)) = previous {
            let ratio = (gap / prev_gap) / (superheat / prev_sh);
            assert!(
                (ratio - 1.0).abs() < 0.05,
                "the gap is proportional to the superheat: {ratio:.4} per decade"
            );
        }
        previous = Some((superheat, gap));
    }
    // And the coefficient, so a change of sign or of scale has something to fail against.
    let per_kelvin = (two_phase_lambda(1e-3) / one - 1.0).abs() / 1e-3;
    println!("  which is {per_kelvin:.4e} per kelvin of superheat");
    assert!(
        (per_kelvin / 9.45e-3 - 1.0).abs() < 0.02,
        "9.45e-3 per kelvin: {per_kelvin:.4e}"
    );
}

/// **A warmer liquid freezes more slowly, by the amount the closed form says.**
///
/// Five superheats, and the front position at 900 s across them. The one-phase model gives 15.85 mm
/// whatever the liquid is doing; the two-phase closed form takes that to 13.33 mm at 20 K, and the
/// point of the test is that the solver follows the second curve and not the first.
#[test]
fn a_superheated_liquid_slows_the_front_by_what_the_closed_form_says() {
    let one = one_phase_lambda();
    let cells = (DEPTH / 1e-3).round() as usize;
    let mut worst_two: f64 = 0.0;
    let mut worst_one: f64 = 0.0;
    for superheat in [0.0f64, 2.0, 5.0, 10.0, 20.0] {
        let lam = two_phase_lambda(superheat.max(1e-12));
        let want = 2.0 * lam * (ALPHA_S * 900.0).sqrt();
        let one_phase_says = 2.0 * one * (ALPHA_S * 900.0).sqrt();

        let mut c = Solid3D::new(
            "column",
            two_phase_ice(),
            (1, 1, cells),
            Length::from_si(1e-3),
            Temperature::celsius(superheat),
        );
        for k in 0..cells {
            c.set_melted_fraction(0, 0, k, 1.0);
        }
        let cold = Temperature::celsius(-DROP);
        // **From the solid's limit, not the block's current one.** Freezing *tightens* the limit:
        // water's is 7.57x ice's, so a step sized on an all-liquid block is refused a few hundred
        // steps later once there is ice in it. The guard catches that rather than diverging, which is
        // how this was found — but a caller either sizes `dt` on the tightest phase, as here, or uses
        // `Schedule::Multirate`, which re-asks every step.
        let dt = Time::from_si(solid_limit() * 0.9);
        let mut t = 0.0;
        while t < 900.0 {
            c.set_temperature(0, 0, 0, cold);
            c.step(Time::from_si(t), dt, &mut Exchange::new())
                .expect("stable");
            t += dt.to_si();
        }
        c.set_temperature(0, 0, 0, cold);
        let solid = DEPTH - c.melted_volume().to_si() / 1e-6;
        let got = solid - 0.5e-3;

        let off_two = (got / want - 1.0).abs();
        let off_one = (got / one_phase_says - 1.0).abs();
        println!(
            "  superheat {superheat:5.1} K   front {:6.2} mm   two-phase says {:6.2} ({:5.2}% out)   \
             one-phase says {:6.2} ({:5.2}% out)",
            got * 1e3,
            want * 1e3,
            off_two * 100.0,
            one_phase_says * 1e3,
            off_one * 100.0
        );
        worst_two = worst_two.max(off_two);
        if superheat >= 20.0 {
            worst_one = off_one;
        }
    }
    // `dx/X` and not its square: the interface cell's conductivity is a mixture, which is a
    // first-order error at the front and is what two phases cost.
    let allowed = 1e-3 / (2.0 * one * (ALPHA_S * 900.0).sqrt());
    println!(
        "  worst against the two-phase form {:.2}%, against a bound of {:.2}% = dx/X",
        worst_two * 100.0,
        allowed * 100.0
    );
    assert!(
        worst_two < allowed,
        "the solver follows the two-phase curve: worst {:.3}%",
        worst_two * 100.0
    );
    // And the one-phase model is wrong by far more than that at 20 K, which is why two phases exist.
    println!(
        "  and at 20 K the one-phase model is {:.1}% out — {:.0}x the bound",
        worst_one * 100.0,
        worst_one / allowed
    );
    assert!(
        worst_one > 3.0 * allowed,
        "the one-phase model should be clearly wrong here, not marginally: {:.2}%",
        worst_one * 100.0
    );
}

/// **The two models disagree by one cell's worth, and refining closes it.**
///
/// In the *continuum* one-phase problem a liquid at the melting point cannot influence anything: every
/// face in it has zero temperature difference, and zero times a conductance is zero. I asserted that of
/// the solver and it is **false**, by 2.1% at twenty cells.
///
/// The reason is the discretisation and not the physics. A fixed grid puts a cell that is *wholly*
/// liquid where the continuum has an interface partway through it, and that cell conducts at water's
/// 0.598 rather than ice's 2.22 — so the conductance from the front to the cold face is throttled while
/// the front is crossing it. One cell, mis-modelled, every step.
///
/// So the honest claim is that the gap is `O(dx)` and closes on refinement, which is what a first-order
/// interface error means. Asserting equality would have been asserting the continuum's property of a
/// scheme that does not have it.
#[test]
fn the_two_models_differ_by_one_cell_and_refining_closes_it() {
    let run = |substance: Substance, cells: usize, dx: f64| {
        let mut c = Solid3D::new(
            "column",
            substance,
            (1, 1, cells),
            Length::from_si(dx),
            Temperature::celsius(0.0),
        );
        for k in 0..cells {
            c.set_melted_fraction(0, 0, k, 1.0);
        }
        let cold = Temperature::celsius(-DROP);
        // Scaled with the grid so both resolutions are the same fraction of their own limit, and from
        // the solid phase for the reason above.
        let dt = Time::from_si(0.4 * dx * dx / (2.0 * ALPHA_S) * 0.9);
        let mut t = 0.0;
        while t < 120.0 {
            c.set_temperature(0, 0, 0, cold);
            c.step(Time::from_si(t), dt, &mut Exchange::new())
                .expect("stable");
            t += dt.to_si();
        }
        c.melted_volume().to_si() / (dx * dx)
    };
    let mut gaps = Vec::new();
    for (cells, dx) in [(20, 1e-3), (40, 5e-4), (80, 2.5e-4)] {
        let one = run(Substance::ice(), cells, dx);
        let two = run(two_phase_ice(), cells, dx);
        let gap = (two / one - 1.0).abs();
        println!(
            "  dx = {:.2} mm: one-phase {:.4} mm melted, two-phase {:.4} — gap {:.3}%",
            dx * 1e3,
            one * 1e3,
            two * 1e3,
            gap * 100.0
        );
        gaps.push(gap);
    }
    let overall = gaps[0] / gaps[2];
    println!("  a fourfold refinement closes it {overall:.2}x, against 4 for first order");
    for pair in gaps.windows(2) {
        assert!(
            pair[1] < pair[0],
            "every refinement must close it: {gaps:?}"
        );
    }
    assert!(
        overall > 2.5,
        "first order or better over four times the cells: {gaps:?} closes {overall:.2}x"
    );
}
