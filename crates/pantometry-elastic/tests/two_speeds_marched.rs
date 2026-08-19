//! `Waves` marched, against the two speeds `Elastic` states in closed form.
//!
//! The speeds themselves are algebra and are checked in `crates/pantometry/tests/two_wave_speeds.rs`.
//! This is the solver: does a body given inertia actually carry a compression at `√((λ+2μ)/ρ)` and a
//! shear at `√(μ/ρ)`.
//!
//! # How a frequency is got out of a march
//!
//! A clamped-clamped span released in its `n`-th half-wave oscillates, and the closed form is
//! `f = n·c/(2L)`. Measured from the mode's projected amplitude — sign changes located by linear
//! interpolation and averaged over many of them — so it is the coefficient of a shape and not a
//! reading off one node.
//!
//! **The leapfrog's own dispersion is removed rather than tolerated.** Central differences on a single
//! eigenmode oscillate at `Ω = 2·arcsin(ω dt/2)` per step, not at `ω dt`. So the measured `Ω` is
//! inverted through `ω = 2·sin(Ω/2)/dt` before anything is compared, leaving the *spatial* error
//! alone. Without that step the two errors add and the sum gets called an accuracy.
//!
//! # Measured
//!
//! ```text
//!   elements   c_p error   c_s error
//!         16      0.161%     0.161%
//!         32      0.040%     0.040%
//!         64      0.010%     0.010%
//!
//!   ratio against sqrt(2(1-v)/(1-2v)):   9.4e-9 at nu = 0.2
//!                                        8.2e-8 at nu = 0.33
//!                                        4.7e-8 at nu = 0.45
//! ```
//!
//! Second order in both — a factor of four per halving, to two figures — and the **ratio holds to
//! about `1e-8`**, four orders tighter than either speed alone. That is because the
//! two modes have the same shape and therefore the same discretisation error, which cancels in a
//! quotient. That is the sharpest statement here: `E` and `ρ` cancel out of the ratio algebraically
//! and the mesh error cancels out of it numerically, so what is left is the operator.

use pantometry_core::{
    units::{Length, Time},
    Domain, Exchange,
};
use pantometry_elastic::{Axis, Elastic, Waves};

/// Elements along the span; one across the other two axes, which makes it the one-dimensional
/// problem the closed forms are about and costs nothing to resolve.
fn column(elements: usize, hold: [Axis; 2]) -> Waves {
    let mut w = Waves::new(
        "column",
        (1, 1, elements),
        Length::mm(1.0),
        Elastic::aluminium_6061(),
    );
    for axis in hold {
        w.hold(axis);
    }
    w.clamp_ends(Axis::Z);
    w
}

/// March and return the frequency of the released mode, in hertz, with the leapfrog's dispersion
/// taken out.
fn measured_frequency(w: &mut Waves, mode: usize, along: Axis, dt: Time, steps: usize) -> f64 {
    let h = dt.to_si();
    let mut crossings = Vec::new();
    let mut previous = w.mode_amplitude(mode, Axis::Z, along);
    assert!(
        previous.abs() > 0.0,
        "the mode has to be present before it can be measured"
    );
    for n in 0..steps {
        w.step(Time::from_si(n as f64 * h), dt, &mut Exchange::new())
            .expect("stable");
        let now = w.mode_amplitude(mode, Axis::Z, along);
        if (previous <= 0.0) != (now <= 0.0) {
            // Linear interpolation between the two samples, so the period is not quantised to a step.
            let frac = previous / (previous - now);
            crossings.push((n as f64 + frac) * h);
        }
        previous = now;
    }
    assert!(
        crossings.len() >= 8,
        "not enough oscillation to measure: {} crossings",
        crossings.len()
    );
    // Half a period between consecutive zero crossings, fitted over all of them.
    let first = crossings[0];
    let last = crossings[crossings.len() - 1];
    let halves = (crossings.len() - 1) as f64;
    let period_steps = 2.0 * (last - first) / halves;
    // Invert the leapfrog relation: the recurrence turns at Ω = 2 arcsin(ω dt / 2) per step.
    let omega_discrete = 2.0 * std::f64::consts::PI / period_steps;
    let omega = 2.0 * (omega_discrete * h / 2.0).sin() / h;
    omega / (2.0 * std::f64::consts::PI)
}

/// **A constrained compression travels at `√((λ+2μ)/ρ)`, converging at second order.**
///
/// Holding `x` and `y` at zero everywhere leaves `ρü_z = (λ+2μ)u_z''`, which is the closed form's own
/// problem: the sides cannot bulge because they are held. That is `λ+2μ` and not `E`, and for
/// aluminium the two differ by 22% — so a scheme that used `E` fails this by far more than any mesh
/// error.
///
/// Three resolutions, because the claim is a **rate**. A single mesh would pass for a scheme that was
/// accidentally right at that spacing, and this workspace has shipped a first-order scheme claiming
/// second once.
#[test]
fn a_constrained_compression_travels_at_the_p_wave_speed() {
    let material = Elastic::aluminium_6061();
    let closed_speed = material.p_wave_speed();
    println!("  c_p = {:.1} m/s", closed_speed.to_si());
    let mut errors = Vec::new();
    for elements in [16, 32, 64] {
        let mut w = column(elements, [Axis::X, Axis::Y]);
        w.release_mode(1, Axis::Z, Axis::Z, Length::from_si(1e-9));
        let dt = Time::from_si(w.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
        let want = w.mode_frequency(1, Axis::Z, closed_speed).to_si();
        // Enough steps for a dozen or so half periods.
        let steps = (12.0 / (want * dt.to_si())) as usize;
        let got = measured_frequency(&mut w, 1, Axis::Z, dt, steps);
        let off = (got / want - 1.0).abs();
        println!(
            "  {elements:3} elements: {got:.4e} Hz against n c/2L {want:.4e} — off {:.3}%",
            off * 100.0
        );
        errors.push(off);
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            order > 1.8,
            "the frequency converges at second order: {errors:?} gives order {order:.3}"
        );
    }
    assert!(
        errors[2] < 1e-3,
        "and the finest mesh is inside a tenth of a percent: {:.4}%",
        errors[2] * 100.0
    );
}

/// **A shear wave travels at `√(μ/ρ)`, and it is slower.**
///
/// Holding `y` and `z` leaves `ρü_x = μ u_x''` — a transverse displacement varying along the span,
/// which is the one wave a fluid cannot carry. Same mesh, same time scheme, same element: only the
/// component that is free has changed, and the speed that comes out is a different one.
#[test]
fn a_shear_wave_travels_at_the_s_wave_speed() {
    let material = Elastic::aluminium_6061();
    let closed_speed = material.s_wave_speed();
    println!(
        "  c_s = {:.1} m/s, against c_p = {:.1}",
        closed_speed.to_si(),
        material.p_wave_speed().to_si()
    );
    let mut errors = Vec::new();
    for elements in [16, 32, 64] {
        let mut w = column(elements, [Axis::Y, Axis::Z]);
        w.release_mode(1, Axis::Z, Axis::X, Length::from_si(1e-9));
        let dt = Time::from_si(w.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
        let want = w.mode_frequency(1, Axis::Z, closed_speed).to_si();
        let steps = (12.0 / (want * dt.to_si())) as usize;
        let got = measured_frequency(&mut w, 1, Axis::X, dt, steps);
        let off = (got / want - 1.0).abs();
        println!(
            "  {elements:3} elements: {got:.4e} Hz against n c/2L {want:.4e} — off {:.3}%",
            off * 100.0
        );
        errors.push(off);
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            order > 1.8,
            "the frequency converges at second order: {errors:?} gives order {order:.3}"
        );
    }
}

/// **The marched ratio of the two speeds is `√(2(1−ν)/(1−2ν))`, and the mesh error cancels out of it.**
///
/// The sharpest thing in this file. `E` and `ρ` cancel out of the ratio algebraically; the *spatial
/// discretisation error* cancels out of it numerically, because the two modes have the same shape on
/// the same mesh and therefore the same error. So this is far tighter than either speed on its own —
/// and it is the check a scheme with the wrong stiffness and a compensating wrong mass cannot pass.
///
/// Three values of Poisson's ratio, because the identity is a *function* of `ν` and one value would
/// pass for a scheme that had it as a constant.
#[test]
fn the_marched_ratio_is_poissons_ratio_and_the_mesh_error_cancels() {
    for nu in [0.2, 0.33, 0.45] {
        let material = Elastic::new(
            pantometry_core::units::Pressure::from_si(68.9e9),
            nu,
            pantometry_core::units::Density::from_si(2700.0),
        )
        .expect("representable");
        let mut ratios = Vec::new();
        for elements in [16, 32] {
            let mut got = [0.0; 2];
            for (which, (hold, along)) in
                [([Axis::X, Axis::Y], Axis::Z), ([Axis::Y, Axis::Z], Axis::X)]
                    .into_iter()
                    .enumerate()
            {
                let mut w = Waves::new("c", (1, 1, elements), Length::mm(1.0), material);
                for axis in hold {
                    w.hold(axis);
                }
                w.clamp_ends(Axis::Z);
                w.release_mode(1, Axis::Z, along, Length::from_si(1e-9));
                let dt = Time::from_si(w.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
                let speed = if which == 0 {
                    material.p_wave_speed()
                } else {
                    material.s_wave_speed()
                };
                let want = w.mode_frequency(1, Axis::Z, speed).to_si();
                let steps = (12.0 / (want * dt.to_si())) as usize;
                got[which] = measured_frequency(&mut w, 1, along, dt, steps);
            }
            ratios.push(got[0] / got[1]);
        }
        let closed = material.speed_ratio();
        let worst = ratios
            .iter()
            .map(|r| (r / closed - 1.0).abs())
            .fold(0.0f64, f64::max);
        println!(
            "  nu = {nu}: marched ratio {:.6} and {:.6} against sqrt(2(1-v)/(1-2v)) {closed:.6} \
             — worst {worst:.2e}",
            ratios[0], ratios[1]
        );
        assert!(
            worst < 1e-5,
            "nu = {nu}: the ratio is the identity, mesh error and all: {ratios:?} against {closed}"
        );
    }
}

/// **A leapfrog's energy swings by `2 sin(ωΔt/2)` and does not drift, and that is not a defect.**
///
/// The kinetic term is built from a backward difference, so it is the velocity half a step ago while
/// the strain term is the displacement now. The two are offset by half a step, which makes the total
/// oscillate at twice the mode frequency by exactly that factor. `Room` records the same thing.
///
/// What matters is that it does not **drift**: the swing at the end is the swing at the beginning. A
/// scheme leaking energy would show a falling envelope, and one gaining it a rising one, and either is
/// invisible in a single period.
#[test]
fn the_energy_swings_by_the_leapfrog_factor_and_does_not_drift() {
    let material = Elastic::aluminium_6061();
    let mut w = column(32, [Axis::X, Axis::Y]);
    w.release_mode(1, Axis::Z, Axis::Z, Length::from_si(1e-9));
    let dt = Time::from_si(w.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let want = w
        .mode_frequency(1, Axis::Z, material.p_wave_speed())
        .to_si();
    let omega = 2.0 * std::f64::consts::PI * want;
    let predicted = 2.0 * (omega * dt.to_si() / 2.0).sin();

    let mut early = (f64::MAX, 0.0f64);
    let mut late = (f64::MAX, 0.0f64);
    let steps = (40.0 / (want * dt.to_si())) as usize;
    for n in 0..steps {
        w.step(
            Time::from_si(n as f64 * dt.to_si()),
            dt,
            &mut Exchange::new(),
        )
        .expect("stable");
        let e = w.total_energy(dt).to_si();
        let bucket = if n < steps / 8 {
            Some(&mut early)
        } else if n > steps - steps / 8 {
            Some(&mut late)
        } else {
            None
        };
        if let Some(b) = bucket {
            b.0 = b.0.min(e);
            b.1 = b.1.max(e);
        }
    }
    let swing = |b: (f64, f64)| (b.1 - b.0) / (0.5 * (b.1 + b.0));
    println!(
        "  swing {:.4} early and {:.4} late, against 2 sin(omega dt/2) = {predicted:.4}",
        swing(early),
        swing(late)
    );
    println!(
        "  and the mean energy went {:.6e} -> {:.6e} J",
        0.5 * (early.0 + early.1),
        0.5 * (late.0 + late.1)
    );
    assert!(
        (swing(late) / predicted - 1.0).abs() < 0.05,
        "the swing is the leapfrog factor: {:.5} against {predicted:.5}",
        swing(late)
    );
    // No drift: the mean over the last eighth is the mean over the first, to a part in a thousand.
    let drift = (0.5 * (late.0 + late.1)) / (0.5 * (early.0 + early.1)) - 1.0;
    assert!(
        drift.abs() < 1e-3,
        "a symplectic scheme does not drift: {:.3e} over {steps} steps",
        drift
    );
}

/// **Past the limit is refused, and the limit is computed from the operator rather than assumed.**
///
/// `2/√λ_max(M⁻¹K)` by Gershgorin, so it is a bound that holds for this element and this mass
/// lumping rather than a Courant number borrowed from a stencil this is not. Reported against the
/// naive `dx/(c_p√3)` so the gap is on the record.
#[test]
fn the_stability_limit_is_the_operators_and_it_is_enforced() {
    let material = Elastic::aluminium_6061();
    let w = Waves::new("cube", (4, 4, 4), Length::mm(1.0), material);
    let limit = w.max_stable_dt(Time::from_si(0.0)).to_si();
    let courant = 1e-3 / (material.p_wave_speed().to_si() * 3.0f64.sqrt());
    println!(
        "  limit {limit:.4e} s against dx/(c_p sqrt 3) = {courant:.4e} — {:.3}x",
        limit / courant
    );
    assert!(limit.is_finite() && limit > 0.0);

    let mut over = Waves::new("cube", (4, 4, 4), Length::mm(1.0), material);
    over.release_mode(1, Axis::Z, Axis::Z, Length::from_si(1e-9));
    let err = over
        .step(
            Time::from_si(0.0),
            Time::from_si(limit * 1.05),
            &mut Exchange::new(),
        )
        .expect_err("5% past the limit must be refused");
    assert_eq!(err.quantity, "leapfrog stability");

    // And just inside it runs, so the check is a limit and not a blanket refusal.
    let mut at = Waves::new("cube", (4, 4, 4), Length::mm(1.0), material);
    at.release_mode(1, Axis::Z, Axis::Z, Length::from_si(1e-9));
    at.step(
        Time::from_si(0.0),
        Time::from_si(limit),
        &mut Exchange::new(),
    )
    .expect("exactly at the limit is accepted");
}

/// **A body nobody has touched sits still, and a checkpoint carries both halves of the state.**
///
/// The first half is the trivial solution of `ρü = ∇·σ` and a scheme that drifts off it has a source
/// nobody asked for. The second is the leapfrog's own requirement: the state is `u` **and** `u`
/// one step ago, because a single displacement does not say which way it is going. A checkpoint of
/// one of them would restore a body with the wrong velocity and nothing would say so.
#[test]
fn an_undisturbed_body_sits_still_and_a_checkpoint_carries_the_velocity() {
    let mut w = Waves::new("still", (3, 3, 3), Length::mm(2.0), Elastic::steel());
    let dt = Time::from_si(w.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
    for n in 0..50 {
        w.step(
            Time::from_si(n as f64 * dt.to_si()),
            dt,
            &mut Exchange::new(),
        )
        .expect("stable");
    }
    assert_eq!(w.strain_energy().to_si(), 0.0, "nothing should have moved");
    assert_eq!(w.displacement_at(1, 1, 1), [0.0, 0.0, 0.0]);

    let mut w = column(16, [Axis::X, Axis::Y]);
    w.release_mode(1, Axis::Z, Axis::Z, Length::from_si(1e-9));
    let dt = Time::from_si(w.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    // Get it moving, so that `u` and `prev` differ and a half-saved state would be detectable.
    for n in 0..40 {
        w.step(
            Time::from_si(n as f64 * dt.to_si()),
            dt,
            &mut Exchange::new(),
        )
        .expect("stable");
    }
    w.checkpoint();
    let mark = w.total_energy(dt).to_si();
    let moving = w.mode_amplitude(1, Axis::Z, Axis::Z);
    for n in 0..40 {
        w.step(
            Time::from_si(n as f64 * dt.to_si()),
            dt,
            &mut Exchange::new(),
        )
        .expect("stable");
    }
    assert!(
        (w.mode_amplitude(1, Axis::Z, Axis::Z) - moving).abs() > 1e-14,
        "it should have moved on before being restored"
    );
    w.restore();
    assert_eq!(
        w.mode_amplitude(1, Axis::Z, Axis::Z),
        moving,
        "a restore returns the displacement exactly"
    );
    assert_eq!(
        w.total_energy(dt).to_si(),
        mark,
        "and the velocity with it, which is what the second array is for"
    );
}

/// **The field the analysis layer draws is the displacement magnitude, and it is not empty.**
///
/// A domain that offers `as_field` and returns nothing useful is the failure this workspace hunts by
/// name: the panel renders, the scale reads "in 1 shades", and nothing anywhere says the picture is of
/// no data. So the field is checked at the nodes it is built from, at the clamped ends where it must
/// be zero, and outside the body where it must clamp rather than extrapolate.
///
/// `|u|` costs the sign, which is stated in the method's own docs and is visible here: a standing
/// half-wave has one antinode and the magnitude has one maximum, but a *full* wave has two antinodes
/// of opposite sign and the magnitude shows both as bright. That is the trade for one scalar.
#[test]
fn the_drawn_field_is_the_displacement_magnitude_and_it_is_populated() {
    use pantometry_core::units::LengthVec;

    let elements = 16;
    let mut w = column(elements, [Axis::X, Axis::Y]);
    w.release_mode(1, Axis::Z, Axis::Z, Length::from_si(1e-9));
    let field = w
        .as_field()
        .expect("a body carrying waves has a field to draw");
    assert_eq!(field.unit(), "m", "a displacement is in metres");

    let dx = 1e-3;
    let at = |k: usize| {
        field.at(
            LengthVec::from_si(glam::DVec3::new(0.0, 0.0, k as f64 * dx)),
            Time::from_si(0.0),
        )
    };
    // Zero at both clamped ends, largest in the middle, and every node in between nonzero.
    assert_eq!(at(0), 0.0, "the clamped end has not moved");
    assert_eq!(at(elements), 0.0, "nor the other one");
    let middle = at(elements / 2);
    assert!(
        (middle - 1e-9).abs() < 1e-21,
        "the antinode is the amplitude: {middle:.6e} against 1e-9"
    );
    let populated = (1..elements).filter(|k| at(*k) > 0.0).count();
    assert_eq!(
        populated,
        elements - 1,
        "every interior node should be displaced, not just the one that was looked at"
    );

    // Halfway between two nodes it interpolates rather than snapping, which is what makes a picture
    // smooth and is the part a nearest-neighbour sampler would get wrong without ever looking empty.
    let between = field.at(
        LengthVec::from_si(glam::DVec3::new(0.0, 0.0, 3.5 * dx)),
        Time::from_si(0.0),
    );
    let (lo, hi) = (at(3), at(4));
    assert!(
        (between - 0.5 * (lo + hi)).abs() < 1e-24,
        "trilinear between nodes: {between:.9e} against {:.9e}",
        0.5 * (lo + hi)
    );

    // Outside the body it clamps. Extrapolating would draw material moving where there is none.
    let past = field.at(
        LengthVec::from_si(glam::DVec3::new(0.0, 0.0, 5.0 * elements as f64 * dx)),
        Time::from_si(0.0),
    );
    assert_eq!(
        past,
        at(elements),
        "past the face it clamps rather than continuing"
    );
    let before = field.at(
        LengthVec::from_si(glam::DVec3::new(0.0, 0.0, -3.0 * dx)),
        Time::from_si(0.0),
    );
    assert_eq!(before, at(0), "and before it, too");
}
