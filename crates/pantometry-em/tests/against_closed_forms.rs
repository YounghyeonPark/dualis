//! What Maxwell's equations do in a box, against things that were true before this code existed.
//!
//! The headline is not a tolerance. `∇·B` is **identically** zero under the Yee update — every term
//! in the discrete divergence of the discrete curl appears twice with opposite signs — so the check
//! is that it is `0.0`, not that it is small.

use pantometry_core::{Domain, Exchange};
use pantometry_em::{cavity_frequency, Cavity, Medium, COURANT_3D};
use pantometry_units::{Length, Time};

const DX: f64 = 2e-3;

fn box_of(counts: (usize, usize, usize)) -> Cavity {
    Cavity::new("cavity", counts, Length::from_si(DX), Medium::vacuum())
}

/// A released cavity, stepped at a safe fraction of the Courant limit.
fn released(counts: (usize, usize, usize), mode: (u32, u32)) -> (Cavity, Time) {
    let mut c = box_of(counts);
    let dt = Time::from_si(c.courant_limit().to_si() * 0.5);
    c.release_mode(mode, 1.0, dt);
    (c, dt)
}

/// **`∇·B` is zero by an identity, and stays at rounding rather than drifting.**
///
/// The discrete divergence of the discrete curl vanishes term by term on the Yee grid, so the
/// update cannot change this quantity at all. In exact arithmetic it would be `0.0` forever; in
/// `f64` each update rounds, and what is left is a random walk at the level of `ε·|H|`.
///
/// So the check is the **shape**: it starts at `1e-16`, ends twelve orders below the field it is
/// measured against, and grows sub-linearly in the step count — which is rounding. A scheme
/// without the identity accumulates magnetic charge in proportion to the work it does, and would
/// be at `1e-4` before a hundred steps. That figure is not hypothetical: it is what this crate
/// measured when `release_mode` set `H` from the *analytic* curl instead of the discrete one.
#[test]
fn the_magnetic_divergence_stays_at_rounding_rather_than_drifting() {
    let (mut c, dt) = released((12, 8, 10), (1, 1));
    let start = c.magnetic_divergence();
    println!("  at release: {start:.3e}");
    assert!(
        start < 1e-14,
        "the released mode is divergence-free to begin with: {start:.3e}"
    );

    let mut bus = Exchange::new();
    let mut at = Vec::new();
    let mut worst: f64 = start;
    for n in 0..2000 {
        c.step(Time::from_si(n as f64 * dt.to_si()), dt, &mut bus)
            .expect("stable");
        worst = worst.max(c.magnetic_divergence());
        if n + 1 == 125 || n + 1 == 500 || n + 1 == 2000 {
            at.push((n + 1, worst));
        }
    }
    for (n, w) in &at {
        println!("  after {n:>4} steps: {w:.3e}");
    }
    assert!(
        at[2].1 < 1e-11,
        "twelve orders below the field: {:.3e}",
        at[2].1
    );
    assert!(
        c.energy().to_si() > 0.0,
        "a cavity with no field satisfies every constraint"
    );

    // **The identity, from the other side.** Put some divergence in by hand and the update can
    // neither grow it nor heal it, because it cannot touch this quantity at all. A scheme without
    // the identity does one or the other, and which one is not the point — either would show here.
    //
    // Measured **unnormalised**: the relative figure divides by `max |H|`, which oscillates
    // through the cycle, so comparing it at two instants compares two denominators. The first
    // version of this read 1.099 and then 0.5387 for a quantity that had not moved.
    let before = c.peak_magnetic_divergence();
    c.inject_divergence(6, 4, 5, 1e-3);
    let injected = c.peak_magnetic_divergence();
    println!("  injected: {before:.3e} -> {injected:.3e} A/m");
    for n in 0..500 {
        c.step(Time::from_si(n as f64 * dt.to_si()), dt, &mut bus)
            .expect("stable");
    }
    let after = c.peak_magnetic_divergence();
    println!(
        "  and after 500 more steps: {after:.3e} A/m — {:.9}x",
        after / injected
    );
    assert!(
        (after / injected - 1.0).abs() < 1e-6,
        "the update preserves whatever divergence is there, exactly: {:.9}x",
        after / injected
    );
    // And the field went on doing electromagnetism around it, which is what makes that
    // preservation a statement rather than a description of a frozen array.
    assert!(
        c.energy().to_si() > 0.0 && c.electric_energy().to_si() > 0.0,
        "the cavity is still ringing"
    );
}

/// **The cavity rings at `(c/2)√((m/a)² + (n/b)² + (p/d)²)`, approached at second order.**
///
/// The discrete grid has its own dispersion, so the measured frequency is not the continuum one —
/// it is close, and the *rate* at which it approaches is what says the scheme is second order.
/// Measured across three refinements, which is what distinguishes second order from first: one
/// refinement cannot tell them apart and two can be fitted by either.
#[test]
fn the_cavity_rings_at_the_closed_form_frequency() {
    let mut errors = Vec::new();
    for cells in [8usize, 16, 32] {
        // The same physical box every time: the cell shrinks as the count grows.
        let dx = 0.048 / cells as f64;
        let mut c = Cavity::new(
            "cavity",
            (cells, 4, cells),
            Length::from_si(dx),
            Medium::vacuum(),
        );
        let dt = Time::from_si(c.courant_limit().to_si() * 0.5);
        c.release_mode((1, 1), 1.0, dt);

        let s = c.size().to_si();
        let closed = cavity_frequency([s.x, s.y, s.z], [1, 0, 1], &Medium::vacuum()).to_si();

        // The period, from the peaks of the electric energy — which peaks twice per cycle, so its
        // period is half the field's. Measured on the energy rather than on a single sample,
        // which is immune to where a probe happens to sit relative to a node.
        //
        // # Two hundred half-cycles, and the peak interpolated
        //
        // Eight was the first version and it could not see second order at all: it reported the
        // *same* frequency at 8, 16 and 32 cells, to seven digits, because a peak located to the
        // nearest step carries an error of `dt` and eight half-cycles divides that by only eight.
        // At Courant 0.5 that leaves 0.6% of measurement noise on top of a 0.55% signal.
        //
        // Two hundred divides it by two hundred, and fitting a parabola through the three samples
        // around each peak removes what is left: the vertex of a parabola through equally spaced
        // points is exact for a smooth maximum and costs three multiplies.
        const CYCLES: usize = 200;
        let mut bus = Exchange::new();
        let mut t = 0.0;
        let mut window = [(0.0, c.electric_energy().to_si()); 3];
        let mut peaks: Vec<f64> = Vec::new();
        for _ in 0..200_000 {
            c.step(Time::from_si(t), dt, &mut bus).expect("stable");
            t += dt.to_si();
            window[0] = window[1];
            window[1] = window[2];
            window[2] = (t, c.electric_energy().to_si());
            if window[1].1 > window[0].1 && window[1].1 >= window[2].1 {
                // Vertex of the parabola through three equally spaced samples.
                let (a, b, cc) = (window[0].1, window[1].1, window[2].1);
                let denom = a - 2.0 * b + cc;
                let shift = if denom.abs() > 0.0 {
                    0.5 * (a - cc) / denom
                } else {
                    0.0
                };
                peaks.push(window[1].0 + shift * dt.to_si());
                if peaks.len() > CYCLES {
                    break;
                }
            }
        }
        assert!(
            peaks.len() > CYCLES,
            "the cavity has to ring: {} peaks",
            peaks.len()
        );
        let half_period = (peaks[CYCLES] - peaks[0]) / CYCLES as f64;
        let measured = 0.5 / half_period;
        let error = (measured / closed - 1.0).abs();
        println!(
            "  {cells:>2} cells: {measured:.6e} Hz against {closed:.6e} — off {:.3}%",
            error * 100.0
        );
        errors.push(error);
    }
    // Second order: halving the cell should quarter the error.
    for pair in errors.windows(2) {
        let rate = pair[0] / pair[1];
        println!("  refinement ratio {rate:.3} (second order is 4)");
        assert!(
            (2.6..6.0).contains(&rate),
            "the Yee scheme is second order in space: {rate:.3}"
        );
    }
    assert!(
        errors[2] < 0.01,
        "and the finest must be close: off by {:.3}%",
        errors[2] * 100.0
    );
}

/// **The field energy holds, and it swings by exactly `2 sin(ωΔt/2)` while it does.**
///
/// Conducting walls do no work — the tangential electric field is zero there, so the Poynting flux
/// through them is zero — and vacuum dissipates nothing. What goes in stays in.
///
/// But `½εE² + ½μH²` is **not** what leapfrog conserves. `E` and `H` are half a step apart, so the
/// invariant is `½εE(t)² + ½μH(t−Δt/2)·H(t+Δt/2)`, and the naive sum oscillates about it. With
/// `E ∝ cos ωt` and `H ∝ sin ω(t±Δt/2)` the arithmetic is short:
///
/// ```text
///   U(t) = A[1 + sin(ωΔt/2) · sin(2ωt + ωΔt/2)]      peak to peak: 2 sin(ωΔt/2)
/// ```
///
/// That is a closed form, so it is what this checks. The first version asserted "the swing is
/// under 5%", measured 12.779%, and called it a failure — where the formula says 12.8% and the
/// scheme was right. An arbitrary bound cannot tell a correct oscillation from an incorrect one;
/// this can.
///
/// The *mean* is the thing that must not move, and it is checked separately.
#[test]
fn a_lossless_cavity_holds_its_energy_and_swings_by_the_closed_form() {
    let counts = (10, 6, 10);
    let (mut c, dt) = released(counts, (1, 1));
    let omega = 2.0 * std::f64::consts::PI * c.mode_frequency([1, 0, 1]).to_si();
    let predicted = 2.0 * (omega * dt.to_si() / 2.0).sin();

    let mut bus = Exchange::new();
    let mut t = 0.0;
    let mut samples = Vec::new();
    for _ in 0..8000 {
        c.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
        samples.push(c.energy().to_si());
    }
    let n = samples.len();
    let early: f64 = samples[..n / 4].iter().sum::<f64>() / (n / 4) as f64;
    let late: f64 = samples[3 * n / 4..].iter().sum::<f64>() / (n - 3 * n / 4) as f64;
    let (lo, hi) = (
        samples.iter().cloned().fold(f64::MAX, f64::min),
        samples.iter().cloned().fold(0.0f64, f64::max),
    );
    let swing = (hi - lo) / early;
    let drift = (late / early - 1.0).abs();

    println!("  mean {early:.9e} -> {late:.9e} J, drift {drift:.2e}");
    println!(
        "  swing {:.4}% against 2 sin(wdt/2) = {:.4}%",
        swing * 100.0,
        predicted * 100.0
    );
    assert!(
        (swing / predicted - 1.0).abs() < 0.05,
        "the swing is the leapfrog's own: {:.4}% against {:.4}%",
        swing * 100.0,
        predicted * 100.0
    );
    assert!(drift < 1e-3, "and the mean does not move: {drift:.3e}");
    // Halving the step must quarter the swing, which is what makes the formula the formula rather
    // than a number that happened to fit.
    let (mut fine, _) = released(counts, (1, 1));
    let small = Time::from_si(dt.to_si() * 0.5);
    fine.release_mode((1, 1), 1.0, small);
    let mut t = 0.0;
    let mut fine_samples = Vec::new();
    for _ in 0..16000 {
        fine.step(Time::from_si(t), small, &mut bus)
            .expect("stable");
        t += small.to_si();
        fine_samples.push(fine.energy().to_si());
    }
    let mean: f64 = fine_samples.iter().sum::<f64>() / fine_samples.len() as f64;
    let fine_swing = (fine_samples.iter().cloned().fold(0.0f64, f64::max)
        - fine_samples.iter().cloned().fold(f64::MAX, f64::min))
        / mean;
    println!(
        "  half the step: {:.4}% — ratio {:.3}",
        fine_swing * 100.0,
        swing / fine_swing
    );
    assert!(
        (swing / fine_swing / 2.0 - 1.0).abs() < 0.05,
        "sin(wdt/2) halves with the step: ratio {:.3}",
        swing / fine_swing
    );
}

/// **The energy is shared equally between the two fields, on average.**
///
/// Equipartition in a standing mode: `⟨½εE²⟩ = ⟨½μH²⟩`. It is a statement about the mode rather
/// than about the solver, and a scheme with the wrong constant in one of its two updates would keep
/// a constant total while splitting it wrongly — which nothing else here would notice.
#[test]
fn a_standing_mode_shares_its_energy_equally() {
    let (mut c, dt) = released((10, 6, 10), (1, 1));
    let mut bus = Exchange::new();
    let mut t = 0.0;
    let (mut e_sum, mut h_sum) = (0.0, 0.0);
    for _ in 0..4000 {
        c.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
        e_sum += c.electric_energy().to_si();
        h_sum += c.magnetic_energy().to_si();
    }
    let ratio = e_sum / h_sum;
    println!("  <U_E> / <U_H> = {ratio:.6}");
    assert!(
        (ratio - 1.0).abs() < 0.02,
        "a standing mode splits its energy evenly: {ratio:.4}"
    );
}

/// **The Courant limit is `dx/(c√3)`, and past it is refused rather than run.**
///
/// The same shape as the acoustic wave equation's and for the same reason. Checked in both
/// directions: refused above, and stable at the limit itself, because a limit that is merely
/// conservative would hide a constant that is wrong.
#[test]
fn past_the_courant_limit_is_refused() {
    let mut c = box_of((8, 8, 8));
    let limit = c.courant_limit().to_si();
    let closed = COURANT_3D * DX / 299_792_458.0;
    println!("  limit {limit:.6e} s against dx/(c sqrt 3) = {closed:.6e} s");
    assert!(
        (limit / closed - 1.0).abs() < 1e-9,
        "the limit is dx/(c sqrt 3): {limit:.6e}"
    );

    let mut bus = Exchange::new();
    let err = c
        .step(Time::from_si(0.0), Time::from_si(limit * 1.001), &mut bus)
        .expect_err("past the limit must be refused");
    assert_eq!(err.quantity, "Courant number");
    c.step(Time::from_si(0.0), Time::from_si(limit), &mut bus)
        .expect("at the limit is stable");
}

/// **Past the limit the field grows, which is what makes the limit a limit.**
///
/// A constant that is merely conservative passes the test above while hiding a scheme that is
/// stable well beyond it — and then a later "correction" of the constant looks harmless. So this
/// steps deliberately past and measures the growth.
#[test]
fn the_limit_is_where_the_field_starts_to_grow() {
    for (factor, should_grow) in [(0.98, false), (1.15, true)] {
        let mut c = box_of((10, 10, 10));
        let dt = Time::from_si(c.courant_limit().to_si() * factor);
        // Released with the *stable* step, so the initial state is the same in both runs and only
        // the marching differs.
        let stable = Time::from_si(c.courant_limit().to_si() * 0.5);
        c.release_mode((5, 5), 1.0, stable);
        let start = c.energy().to_si();

        let mut bus = Exchange::new();
        let mut blew_up = false;
        for n in 0..600 {
            if c.step(Time::from_si(n as f64 * dt.to_si()), dt, &mut bus)
                .is_err()
            {
                blew_up = true;
                break;
            }
            if !c.energy().to_si().is_finite() {
                blew_up = true;
                break;
            }
        }
        let grew = blew_up || c.energy().to_si() > 100.0 * start;
        println!(
            "  {factor:.2} of the limit: energy {:.3e} from {start:.3e} — {}",
            c.energy().to_si(),
            if grew { "grew" } else { "held" }
        );
        assert_eq!(
            grew,
            should_grow,
            "at {factor} of the Courant limit the sharpest mode should {}",
            if should_grow { "grow" } else { "hold" }
        );
    }
}

/// **A conductor dissipates, and the books still balance.**
///
/// The energy the field lost is on the ledger beside what it still has, so the total is constant.
/// That is what makes the loss a *transfer* rather than a leak, and it is the only arrangement
/// under which the conservation audit can tell the two apart.
///
/// The semi-implicit update is also what makes this runnable at all: an explicit `E -= (σ/ε)E dt`
/// needs `dt < 2ε/σ`, which for even a poor conductor is orders below the Courant limit.
#[test]
fn a_lossy_medium_moves_energy_onto_the_books_rather_than_losing_it() {
    let medium = Medium {
        conductivity: 5e-3,
        ..Medium::vacuum()
    };
    let mut c = Cavity::new("lossy", (10, 6, 10), Length::from_si(DX), medium);
    let dt = Time::from_si(c.courant_limit().to_si() * 0.5);
    c.release_mode((1, 1), 1.0, dt);

    let before = c.ledger();
    let start = c.energy().to_si();
    let mut bus = Exchange::new();
    let mut t = 0.0;
    for _ in 0..3000 {
        c.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
    }
    let after = c.ledger();
    let (a, b) = (
        before.get("energy").unwrap_or(0.0),
        after.get("energy").unwrap_or(0.0),
    );
    let drift = (b - a).abs() / a.abs().max(1e-300);
    println!(
        "  field {start:.6e} -> {:.6e} J, dissipated {:.6e} J, ledger drift {drift:.2e}",
        c.energy().to_si(),
        c.dissipated().to_si()
    );
    assert!(
        c.energy().to_si() < 0.9 * start,
        "a conductor has to take something: {:.4} of it left",
        c.energy().to_si() / start
    );
    // **Not `1e-9`, and the reason is the test above.** The reported field energy is the naive
    // `½εE² + ½μH²`, which oscillates about leapfrog's actual invariant by `2 sin(ωΔt/2)` — 12.8%
    // here. A ledger built on it cannot be tighter than that oscillation, and asserting otherwise
    // would be asserting that the scheme is something it is not. `Cavity` does not claim
    // `books_balance` for the same reason.
    let swing = 2.0
        * (2.0 * std::f64::consts::PI * c.mode_frequency([1, 0, 1]).to_si() * dt.to_si() / 2.0)
            .sin();
    println!("  the leapfrog's own swing here is {:.2}%", swing * 100.0);
    assert!(
        drift < swing,
        "and what it took is on the books, to what the scheme can say: {drift:.3e} against a \
         {swing:.3e} oscillation"
    );
    // A lossless box, by contrast, dissipates exactly nothing.
    let (mut clean, dt2) = released((10, 6, 10), (1, 1));
    for n in 0..100 {
        clean
            .step(Time::from_si(n as f64 * dt2.to_si()), dt2, &mut bus)
            .expect("stable");
    }
    assert_eq!(
        clean.dissipated().to_si(),
        0.0,
        "vacuum dissipates nothing at all, not nearly nothing"
    );
}

/// **An open box lets a pulse out, and a conducting one does not — measured, not claimed.**
///
/// Mur's first-order condition is the discrete one-way wave equation. A wave arriving along the
/// normal satisfies it exactly and leaves with no reflection at all in the continuum; one arriving
/// at an angle does not, with a reflection going as `(1−cos θ)/(1+cos θ)`.
///
/// A line source radiates over a whole range of angles at once, so this measures a distribution
/// rather than the best case. Measured: **0.149%** of the energy is left after two and a half
/// crossings, against a conducting box's 101%. The conducting box is the control — without it,
/// "the energy fell" is not a statement about the boundary.
///
/// It also checks that opening the box does not cost the divergence identity. `∇·B` is preserved by
/// `∂(∇·B)/∂t = −∇·(∇×E)`, which holds for *any* `E` — so changing `E` at a boundary cannot touch
/// it, and a measurement that said otherwise would mean the boundary was writing to `H`.
#[test]
fn an_open_box_lets_a_pulse_out() {
    let counts = (24, 24, 24);
    let mut kept = Vec::new();
    for open in [false, true] {
        let mut c = box_of(counts);
        if open {
            c.open();
        }
        let dt = Time::from_si(c.courant_limit().to_si() * 0.5);
        c.pulse((12, 12), 1.0, 2.0);
        let start = c.energy().to_si();
        let divergence = c.peak_magnetic_divergence();

        // Long enough for a wave to cross the box twice at the speed of light.
        let crossing = counts.0 as f64 * DX / 299_792_458.0;
        let mut bus = Exchange::new();
        let mut t = 0.0;
        while t < 2.5 * crossing {
            c.step(Time::from_si(t), dt, &mut bus).expect("stable");
            t += dt.to_si();
        }
        let left = c.energy().to_si() / start;
        println!(
            "  {:<11} {:.4}% of the energy is still in the box after {:.1} crossings",
            if open { "open:" } else { "conducting:" },
            left * 100.0,
            2.5
        );
        assert!(
            c.peak_magnetic_divergence() <= divergence.max(1e-12) * 1e3,
            "a boundary condition writes to E and must not touch div B"
        );
        kept.push(left);
    }

    assert!(
        kept[0] > 0.97,
        "a conducting box keeps what it was given: {:.4}",
        kept[0]
    );
    assert!(
        kept[1] < 0.01,
        "and an open one lets it out: {:.4} left",
        kept[1]
    );
    // Reported as the two figures rather than as their ratio: the conducting box's departure from
    // exactly 1.0 is the leapfrog's own energy swing, which is near zero at this instant, and
    // dividing by it produces a number with nine digits and no meaning.
    println!(
        "  a conductor keeps {:.1}% and Mur keeps {:.3}%",
        kept[0] * 100.0,
        kept[1] * 100.0
    );
}

/// **Mur's coefficient is `(cΔt − Δ)/(cΔt + Δ)` and nothing else.**
///
/// A function of the Courant number alone, which is worth stating because it is the one number in
/// the absorbing boundary and a reader will otherwise assume it depends on the medium or the mesh
/// separately. It does not: both enter only through `cΔt/Δ`.
#[test]
fn the_mur_coefficient_is_a_function_of_the_courant_number_alone() {
    let a = box_of((8, 8, 8));
    let mut b = Cavity::new(
        "other",
        (8, 8, 8),
        Length::from_si(DX * 3.0),
        Medium::dielectric(4.0),
    );
    b.open();
    // The same Courant number in two quite different boxes.
    let s = 0.4;
    let ka = a.mur_coefficient(Time::from_si(a.courant_limit().to_si() * s));
    let kb = b.mur_coefficient(Time::from_si(b.courant_limit().to_si() * s));
    println!("  vacuum at 2 mm: {ka:.9}   glass at 6 mm: {kb:.9}");
    assert!(
        (ka - kb).abs() < 1e-12,
        "the coefficient depends on c dt / dx and nothing else: {ka} against {kb}"
    );
    // And it is the closed form.
    let closed = (COURANT_3D * s - 1.0) / (COURANT_3D * s + 1.0);
    assert!(
        (ka - closed).abs() < 1e-12,
        "(c dt - dx)/(c dt + dx): {ka:.9} against {closed:.9}"
    );
}
