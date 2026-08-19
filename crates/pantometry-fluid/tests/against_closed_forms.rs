//! What incompressible flow does, against the three exact solutions that exist.
//!
//! Each is chosen to be blind to a different mistake, and the module docs say which. Poiseuille and
//! Couette are steady and unidirectional, so the advection term is **identically zero** in both and
//! neither can check it. Taylor–Green can: it solves the complete nonlinear equations. Beside them
//! sit two machine-precision statements a decay rate is too coarse to see — a uniform flow must
//! stay exactly uniform, and momentum in a periodic box must not move at all.

use glam::DVec3;
use pantometry_core::{Domain, Exchange};
use pantometry_fluid::{
    poiseuille_mean_speed, taylor_green_rate, Channel, Fluid, Walls, CELL_REYNOLDS_LIMIT,
};
use pantometry_units::{Density, Diffusivity, Length, Time};

/// A fluid whose viscosity is convenient rather than real, so a small box runs at a Reynolds
/// number the scheme is honest at.
fn syrup(nu: f64) -> Fluid {
    Fluid::new(Density::from_si(1000.0), Diffusivity::from_si(nu))
}

fn run(c: &mut Channel, seconds: f64, safety: f64) -> f64 {
    let mut bus = Exchange::new();
    let mut t = 0.0;
    while t < seconds {
        let dt = (c.max_stable_dt(Time::from_si(t)).to_si() * safety).min(seconds - t);
        if dt <= 0.0 {
            break;
        }
        c.step(Time::from_si(t), Time::from_si(dt), &mut bus)
            .expect("stable");
        t += dt;
    }
    t
}

/// **Poiseuille flow comes out exactly — against the *discrete* parabola, which is not the
/// continuum one.**
///
/// A body force between two no-slip walls settles to `u(y) = (g/2ν)y(h−y)` in the continuum, whose
/// mean is `gh²/12ν`. The discrete answer is a different parabola, and the difference is a closed
/// form rather than an error term:
///
/// ```text
///   u_j = (g/2ν) [ (h² + Δ²)/4 − (y_j − h/2)² ]        mean = (gh²/12ν)(1 + 2/n²)
/// ```
///
/// # Where the `2/n²` comes from
///
/// The interior equation is satisfied by **any** quadratic — a second difference of one is exact —
/// so the discretisation contributes nothing there. The wall does. A no-slip condition imposed by
/// reflecting the first cell, `u₋₁ = −u₀`, makes the *linear interpolation* between the two
/// vanish at the wall, and a parabola is not its own linear interpolation: the ghost is off by the
/// curvature, `(g/4ν)Δ²`.
///
/// So the test is two statements. The discrete form holds to **machine precision at any mesh**,
/// which says the scheme is right; and the gap to the continuum falls as `1/n²`, which says the
/// boundary treatment is second order. Comparing to the continuum alone and calling the difference
/// a tolerance would have hidden both.
///
/// The first version of this test did exactly that, and measured 19/18 at six cells — which is
/// `1 + 2/36` to the digit, and was the scheme telling it so.
#[test]
fn poiseuille_flow_is_exact_against_the_discrete_parabola() {
    let nu = 1e-3;
    let g = 0.02;
    let mut gaps = Vec::new();
    for ny in [6usize, 10, 16] {
        let h = 0.02;
        let dx = h / ny as f64;
        let mut c = Channel::new(
            "channel",
            (4, ny, 4),
            Length::from_si(dx),
            syrup(nu),
            Walls::Sliding {
                low: 0.0,
                high: 0.0,
            },
        );
        c.drive(DVec3::new(g, 0.0, 0.0));
        // Long enough to be steady: the viscous time across the gap is `h²/ν`.
        run(&mut c, 8.0 * h * h / nu, 0.4);

        let continuum = poiseuille_mean_speed(g, h, nu);
        let discrete = continuum * (1.0 + 2.0 / (ny * ny) as f64);
        let mean = c.mean_speed().to_si();
        let mut worst: f64 = 0.0;
        for j in 0..ny {
            let y = (j as f64 + 0.5) * dx;
            let exact = g / (2.0 * nu) * ((h * h + dx * dx) / 4.0 - (y - h / 2.0).powi(2));
            worst = worst.max((c.layer_speed(j).to_si() - exact).abs() / exact);
        }
        let gap = mean / continuum - 1.0;
        println!(
            "  {ny:>2} cells: mean {mean:.9e} against the discrete {discrete:.9e} — off {:.2e}; \
             worst layer {worst:.2e}; above the continuum by {:.4}%",
            (mean / discrete - 1.0).abs(),
            gap * 100.0
        );
        assert!(
            (mean / discrete - 1.0).abs() < 1e-9,
            "the discrete parabola is exact: {mean:.9e} against {discrete:.9e}"
        );
        assert!(
            worst < 1e-9,
            "and so is every layer of it: worst off by {worst:.3e}"
        );
        gaps.push(gap);
    }
    // And the gap to the continuum closes at second order in the mesh.
    for pair in gaps.windows(2) {
        let rate = pair[0] / pair[1];
        println!("  the gap to the continuum shrank {rate:.3}x");
    }
    let predicted = (16.0f64 / 6.0).powi(2);
    let measured = gaps[0] / gaps[2];
    println!("  6 cells to 16: {measured:.3}x against (16/6)^2 = {predicted:.3}");
    assert!(
        (measured / predicted - 1.0).abs() < 1e-6,
        "the wall's error is exactly 2/n^2: {measured:.4} against {predicted:.4}"
    );
}

/// **Couette flow is a straight line, and it is exact for the same reason.**
///
/// One wall slides, no force, and the steady profile is `u = U·y/h`. Linear, so the second
/// difference is exactly zero and the answer is exact. It is also completely blind to advection —
/// there is nothing to advect — which is what makes it a check on the wall condition and on
/// nothing else.
#[test]
fn couette_flow_is_a_straight_line() {
    let (nu, u_wall, ny) = (1e-3, 0.01, 12usize);
    let h = 0.02;
    let dx = h / ny as f64;
    let mut c = Channel::new(
        "couette",
        (4, ny, 4),
        Length::from_si(dx),
        syrup(nu),
        Walls::Sliding {
            low: 0.0,
            high: u_wall,
        },
    );
    run(&mut c, 10.0 * h * h / nu, 0.4);

    let mut worst: f64 = 0.0;
    for j in 0..ny {
        let y = (j as f64 + 0.5) * dx;
        let closed = u_wall * y / h;
        worst = worst.max((c.layer_speed(j).to_si() - closed).abs() / u_wall);
    }
    println!("  worst layer off by {worst:.3e} of the wall speed");
    assert!(worst < 1e-9, "a linear profile is exact: {worst:.3e}");
    // And the mean is half the wall speed, which is the integral of the line.
    assert!(
        (c.mean_speed().to_si() / (0.5 * u_wall) - 1.0).abs() < 1e-9,
        "the mean of a straight line is its midpoint: {:.6e}",
        c.mean_speed().to_si()
    );
}

/// **Taylor–Green decays at `2νk²`, which is the only check here that sees the advection.**
///
/// An exact solution of the **complete** nonlinear equations: the advection term is not absent, it
/// is balanced by the pressure gradient. So this is the one test in the file that would notice a
/// wrong advection scheme, and the rate is a closed form rather than a fitted number.
///
/// Measured across two meshes, because the decay of a *discrete* mode is not exactly `2νk²` — the
/// discrete Laplacian's eigenvalue is, and it approaches the continuum one at second order.
///
/// # Two errors that partly cancel, and both are accounted for
///
/// Against the discrete Laplacian's own rate `ν·4(1−cos kΔ)/Δ²` the measurement comes out **high**
/// by 1.02% and 0.26% — which is forward Euler's `σΔt/2` to the digit, since `Δt` goes as `Δ²`
/// here. Against the continuum it comes out **low** by 0.27% and 0.07%, because the spatial error
/// is the larger of the two and has the other sign.
///
/// Knowing both is what makes the continuum comparison meaningful. A scheme with a genuine defect
/// would not land between two errors that each have a formula.
#[test]
fn a_taylor_green_vortex_decays_at_the_closed_form_rate() {
    let nu = 2e-3;
    let mut errors = Vec::new();
    for n in [16usize, 32] {
        let l = 0.04;
        let dx = l / n as f64;
        let k = 2.0 * std::f64::consts::PI / l;
        let mut c = Channel::new(
            "vortex",
            (n, n, 4),
            Length::from_si(dx),
            syrup(nu),
            Walls::None,
        );
        let u0 = 0.02;
        c.set_velocity(|p| {
            DVec3::new(
                u0 * (k * p.x).cos() * (k * p.y).sin(),
                -u0 * (k * p.x).sin() * (k * p.y).cos(),
                0.0,
            )
        });

        let rate = taylor_green_rate(k, nu);
        let start = c.kinetic_energy().to_si();
        // One e-folding of the energy, which decays at twice the velocity's rate.
        let target = 1.0 / (2.0 * rate);
        let t = run(&mut c, target, 0.4);
        let end = c.kinetic_energy().to_si();
        let measured = -(end / start).ln() / (2.0 * t);
        let error = (measured / rate - 1.0).abs();
        println!(
            "  {n:>2}^2: decay {measured:.6e} /s against 2nu k^2 {rate:.6e} — off {:.3}%",
            error * 100.0
        );
        assert!(
            end < 0.5 * start,
            "the vortex has to actually decay: {:.4} of it left",
            end / start
        );
        errors.push(error);
    }
    let rate = errors[0] / errors[1];
    println!("  refinement ratio {rate:.2} (second order is 4)");
    assert!(
        (2.6..6.0).contains(&rate),
        "the scheme is second order: {rate:.3}"
    );
    assert!(
        errors[1] < 2e-3,
        "and the finer mesh must be close: off by {:.4}%",
        errors[1] * 100.0
    );
}

/// **A uniform flow stays exactly uniform**, which catches almost every interpolation mistake.
///
/// Free-stream preservation. Advection of a constant is zero, diffusion of a constant is zero, and
/// the divergence of a constant is zero, so nothing should happen — and a scheme that interpolated
/// a velocity to the wrong place, or differenced across a periodic seam incorrectly, produces a
/// ripple immediately. It costs nothing to check and it is the sharpest test in the file.
#[test]
fn a_uniform_flow_is_left_alone() {
    let mut c = Channel::new(
        "stream",
        (8, 6, 5),
        Length::from_si(1e-3),
        syrup(1e-4),
        Walls::None,
    );
    let flow = DVec3::new(0.03, -0.02, 0.01);
    c.set_velocity(|_| flow);
    let before = c.momentum_x();

    let mut worst: f64 = 0.0;
    let mut bus = Exchange::new();
    let dt = Time::from_si(c.max_stable_dt(Time::from_si(0.0)).to_si() * 0.4);
    for n in 0..500 {
        c.step(Time::from_si(n as f64 * dt.to_si()), dt, &mut bus)
            .expect("stable");
        let (nx, ny, nz) = c.counts();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let v = c.velocity_at(i, j, k);
                    worst = worst.max((v - flow).length() / flow.length());
                }
            }
        }
    }
    println!("  after 500 steps the worst cell is off by {worst:.3e} of the stream");
    assert!(
        worst < 1e-12,
        "a uniform flow is a solution and must be left exactly alone: {worst:.3e}"
    );
    assert!(
        (c.momentum_x() / before - 1.0).abs() < 1e-12,
        "and its momentum with it"
    );
}

/// **Momentum in a periodic box is conserved exactly.**
///
/// Flux form: every face's contribution to the advection appears twice with opposite signs, and
/// the pressure gradient of a periodic field sums to zero. So the total is unchanged to machine
/// precision — not to the truncation error, which is what the convective form would give.
///
/// Checked on a flow that is genuinely doing something, or conservation is free.
#[test]
fn momentum_in_a_periodic_box_does_not_move() {
    let mut c = Channel::new(
        "box",
        (16, 16, 4),
        Length::from_si(2.5e-3),
        syrup(2e-3),
        Walls::None,
    );
    let k = 2.0 * std::f64::consts::PI / 0.04;
    // A vortex plus a mean flow: the mean is what is conserved, the vortex is what makes the
    // advection non-trivial while it is conserved.
    c.set_velocity(|p| {
        DVec3::new(
            0.01 + 0.02 * (k * p.x).cos() * (k * p.y).sin(),
            -0.02 * (k * p.x).sin() * (k * p.y).cos(),
            0.0,
        )
    });
    let before = c.momentum_x();
    let energy_before = c.kinetic_energy().to_si();

    let mut bus = Exchange::new();
    let dt = Time::from_si(c.max_stable_dt(Time::from_si(0.0)).to_si() * 0.4);
    for n in 0..400 {
        c.step(Time::from_si(n as f64 * dt.to_si()), dt, &mut bus)
            .expect("stable");
    }
    let drift = (c.momentum_x() / before - 1.0).abs();
    println!(
        "  momentum {before:.9e} -> {:.9e}, drift {drift:.3e}",
        c.momentum_x()
    );
    println!(
        "  while the energy fell to {:.4} of its start, so the flow was not idle",
        c.kinetic_energy().to_si() / energy_before
    );
    assert!(
        drift < 1e-11,
        "flux form conserves momentum exactly: {drift:.3e}"
    );
    assert!(
        c.kinetic_energy().to_si() < 0.95 * energy_before,
        "the vortex has to be doing something"
    );
}

/// **The divergence after a projection is the pressure solve's residual, and it is reported.**
///
/// Weaker than electromagnetism's `∇·B`, which is an identity of the update: here it holds to
/// whatever the conjugate-gradient solve was asked for, and knowing which of the two you have is
/// the difference between a guarantee and a hope.
#[test]
fn the_projection_leaves_the_flow_divergence_free_to_the_solve() {
    let mut c = Channel::new(
        "box",
        (12, 10, 6),
        Length::from_si(2e-3),
        syrup(1e-3),
        Walls::Sliding {
            low: 0.0,
            high: 0.02,
        },
    );
    c.drive(DVec3::new(0.01, 0.0, 0.0));
    let mut bus = Exchange::new();
    let dt = Time::from_si(c.max_stable_dt(Time::from_si(0.0)).to_si() * 0.4);
    let mut worst: f64 = 0.0;
    for n in 0..300 {
        c.step(Time::from_si(n as f64 * dt.to_si()), dt, &mut bus)
            .expect("stable");
        worst = worst.max(c.divergence());
    }
    println!(
        "  worst divergence {worst:.3e} m/s against a peak speed of {:.3e} m/s",
        c.peak_speed()
    );
    println!("  the last solve's residual was {:.3e}", c.residual());
    assert!(c.converged(), "the pressure solve has to converge");
    assert!(
        worst < 1e-9 * c.peak_speed().max(1e-12),
        "the projection is only as good as its solve: {worst:.3e}"
    );
}

/// **The cell Reynolds number is a property of the mesh, and past two the step is refused.**
///
/// The one limit here that shortening the step cannot fix. Checked in both directions: a mesh
/// under the limit runs, and the same flow on a coarser mesh is refused — with the *same* time
/// step, so that nothing about the refusal is about time.
#[test]
fn too_coarse_a_mesh_is_refused_however_short_the_step() {
    let flow = 0.05;
    let nu = 1e-4;
    let mut bus = Exchange::new();
    for (dx, should_run) in [(1e-3, true), (1e-2, false)] {
        let mut c = Channel::new(
            "box",
            (8, 8, 4),
            Length::from_si(dx),
            syrup(nu),
            Walls::None,
        );
        c.set_velocity(|_| DVec3::new(flow, 0.0, 0.0));
        let re = c.cell_reynolds();
        // A step far below every other limit, so the only thing that can refuse is the mesh.
        let dt = Time::from_si(c.max_stable_dt(Time::from_si(0.0)).to_si() * 1e-4);
        let ok = c.step(Time::from_si(0.0), dt, &mut bus).is_ok();
        println!(
            "  dx {dx:.0e}: cell Reynolds {re:.2}, step {}",
            if ok { "taken" } else { "refused" }
        );
        assert_eq!(
            ok, should_run,
            "at Re_cell = {re:.2} against a limit of {CELL_REYNOLDS_LIMIT}"
        );
        if !ok {
            let err = c
                .step(Time::from_si(0.0), dt, &mut bus)
                .expect_err("refused");
            assert_eq!(err.quantity, "cell Reynolds number");
        }
    }
}

/// **Past the viscous limit the flow is refused rather than run.**
#[test]
fn past_the_viscous_limit_is_refused() {
    let mut c = Channel::new(
        "box",
        (6, 6, 6),
        Length::from_si(1e-3),
        syrup(1e-4),
        Walls::None,
    );
    let closed = 1e-6 / (6.0 * 1e-4);
    println!(
        "  limit {:.6e} s against dx^2/6nu = {closed:.6e} s",
        c.viscous_limit().to_si()
    );
    assert!(
        (c.viscous_limit().to_si() / closed - 1.0).abs() < 1e-12,
        "the viscous limit is the Fourier one"
    );
    let mut bus = Exchange::new();
    let err = c
        .step(Time::from_si(0.0), Time::from_si(closed * 1.01), &mut bus)
        .expect_err("past the limit must be refused");
    assert_eq!(err.quantity, "flow step");
    c.step(Time::from_si(0.0), Time::from_si(closed), &mut bus)
        .expect("at the limit is stable");
}
