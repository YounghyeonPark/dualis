//! The GPU against the domain it accelerates.
//!
//! Not "does it run" — every wrong stencil runs. The question is how far a single-precision port
//! lands from the `f64` domain that is the reference, and whether that distance is the shape of
//! rounding or the shape of a bug.
//!
//! # These tests skip when there is no GPU, and say so
//!
//! A CI runner usually has no adapter. Skipping is the honest outcome — the alternative is a
//! software rasteriser, which would be checking a different implementation than anyone runs. What
//! is **not** acceptable is skipping quietly, so every skip prints why.

use dualis_core::units::{Energy, Length, Temperature, Time};
use dualis_core::{Domain, Exchange, Substance};
use dualis_gpu::GpuSolid;
use dualis_thermal::Solid3D;

const N: usize = 16;
const DX: f64 = 1e-3;

/// Build the pair. `None` when this machine has no GPU, with a reason on stdout.
fn pair() -> Option<(Solid3D, GpuSolid)> {
    let cpu = Solid3D::new(
        "cpu",
        Substance::aluminium_6061(),
        (N, N, N),
        Length::from_si(DX),
        Temperature::celsius(20.0),
    );
    match GpuSolid::new(
        "gpu",
        Substance::aluminium_6061(),
        (N, N, N),
        Length::from_si(DX),
        Temperature::celsius(20.0),
    ) {
        Ok(gpu) => Some((cpu, gpu)),
        Err(why) => {
            println!("skipped: {why}. Nothing here can run without one, and a software adapter");
            println!("         would be checking a different implementation than anyone uses.");
            None
        }
    }
}

/// Deposit the same joules in the same cell on both.
fn seed(cpu: &mut Solid3D, gpu: &mut GpuSolid) {
    cpu.deposit(N / 2, N / 2, N / 2, Energy::from_si(2.0));
    gpu.deposit(N / 2, N / 2, N / 2, Energy::from_si(2.0));
}

fn run(cpu: &mut Solid3D, gpu: &mut GpuSolid, steps: usize, dt: Time) {
    let mut bus = Exchange::new();
    let mut t = 0.0;
    for _ in 0..steps {
        cpu.step(Time::from_si(t), dt, &mut bus).expect("stable");
        gpu.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
    }
}

/// The largest relative difference between the two grids, measured against the range the field
/// actually spans rather than against each cell — a cell at ambient is 293 K and dividing by it
/// would report every real difference as tiny.
fn divergence(cpu: &Solid3D, gpu: &mut GpuSolid) -> f64 {
    let ambient = Temperature::celsius(20.0).to_si();
    let mut worst = 0.0f64;
    let mut scale: f64 = 1e-30;
    let cells = gpu.cells();
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                let a = cpu.temperature_at(i, j, k).to_si();
                scale = scale.max((a - ambient).abs());
            }
        }
    }
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                let a = cpu.temperature_at(i, j, k).to_si();
                let b = cells[i + N * (j + N * k)];
                worst = worst.max((a - b).abs() / scale);
            }
        }
    }
    worst
}

/// **The GPU reproduces the reference to single precision, and the figure is reported.**
///
/// The claim is not that they agree. WGSL has no `f64`, so they cannot: the port is a
/// lower-precision arithmetic and the only useful question is how much lower.
///
/// The tolerance is earned rather than tried. `f32` has about 7 decimal digits, so one update —
/// a sum of seven terms and two multiplies — loses of order `1e-7` relative. Over `k` steps the
/// error walks as roughly `√k`, so 200 steps is about `1.4e-6`. Anything at that order is
/// rounding; anything above it by orders is a different stencil.
#[test]
fn the_gpu_reproduces_the_cpu_to_single_precision() {
    let Some((mut cpu, mut gpu)) = pair() else {
        return;
    };
    seed(&mut cpu, &mut gpu);
    let dt = Time::from_si(cpu.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);

    // Sixty, not two hundred. At half the stability limit the diffusion length after `k` steps is
    // `dx·√(k/3)`, so 200 steps reaches 8 mm on a 16 mm block — the spot has hit the walls and
    // levelled, and comparing two nearly-uniform grids is agreement for free. Sixty puts it at
    // 4.5 mm, a quarter of the block, where there is still structure to disagree about.
    let steps = 60;
    run(&mut cpu, &mut gpu, steps, dt);
    let worst = divergence(&cpu, &mut gpu);
    let expected = 1e-7 * (steps as f64).sqrt();

    println!("  after {steps} steps: worst relative difference {worst:.3e}");
    println!("  single precision over that many steps predicts about {expected:.3e}");
    assert!(
        worst < 20.0 * expected,
        "the two have diverged by more than rounding: {worst:.3e} against {expected:.3e}"
    );
    // And they are genuinely both computing something, or agreeing is free. Measured against the
    // deposit rather than against a round number: 2 J into one cell of this block is an 827 K
    // rise, and after 200 steps it has spread but is nowhere near level.
    let spread = cpu.peak_temperature().to_si() - cpu.coldest_temperature().to_si();
    let uniform = 2.0
        / Substance::aluminium_6061()
            .heat_capacity(cpu.volume())
            .unwrap()
            .to_si();
    println!(
        "  the spot is still {:.0}x the levelled rise",
        spread / uniform
    );
    assert!(
        spread > 5.0 * uniform,
        "the spot must still be a spot: {spread:.4} K against {uniform:.4} K levelled"
    );
}

/// **The divergence is the shape of rounding, not of a bug.**
///
/// A wrong stencil — a missing arm, a swapped axis, a mirror that is a zero — does not drift; it
/// is wrong immediately and stays wrong. Rounding grows slowly. So the check is the *shape*: after
/// ten steps the two must be far closer than after two hundred, and both far below what a real
/// disagreement would give.
#[test]
fn the_difference_grows_like_rounding_rather_than_appearing_at_once() {
    let Some((mut cpu, mut gpu)) = pair() else {
        return;
    };
    seed(&mut cpu, &mut gpu);
    let dt = Time::from_si(cpu.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);

    run(&mut cpu, &mut gpu, 10, dt);
    let early = divergence(&cpu, &mut gpu);
    run(&mut cpu, &mut gpu, 50, dt);
    let late = divergence(&cpu, &mut gpu);

    println!("  10 steps {early:.3e}   60 steps {late:.3e}");
    assert!(
        late > early,
        "rounding accumulates; a wrong stencil would be wrong at step one"
    );
    assert!(
        early < 1e-6,
        "ten steps of f32 should be near-exact, was {early:.3e}"
    );
}

/// **Energy is conserved on the GPU too, to what single precision can hold.**
///
/// The faces are insulated and nothing is on the bus, so the mean is fixed. In `f64` the CPU holds
/// that to about `1e-12`; `f32` cannot, and the gap is the honest cost of the port.
///
/// This is why `GpuSolid` declines `books_balance` and why a scene using it needs
/// `conservation_tolerance_for(ENERGY, ..)` loosened: `Simulation`'s default `1e-9` is below what
/// the arithmetic can deliver, and a run would be refused for being single precision rather than
/// for being wrong.
#[test]
fn the_gpu_conserves_to_what_f32_can_hold() {
    let Some((mut cpu, mut gpu)) = pair() else {
        return;
    };
    seed(&mut cpu, &mut gpu);
    let dt = Time::from_si(cpu.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);

    let before_cpu = cpu.mean_temperature().to_si();
    let before_gpu = gpu.mean_temperature().to_si();
    run(&mut cpu, &mut gpu, 60, dt);
    let cpu_drift = (cpu.mean_temperature().to_si() - before_cpu).abs() / before_cpu;
    let gpu_drift = (gpu.mean_temperature().to_si() - before_gpu).abs() / before_gpu;

    println!("  mean drift: cpu {cpu_drift:.3e}   gpu {gpu_drift:.3e}");
    assert!(
        cpu_drift < 1e-12,
        "the f64 reference is exact: {cpu_drift:.3e}"
    );
    assert!(
        gpu_drift < 1e-6,
        "even f32 should conserve to a millionth: {gpu_drift:.3e}"
    );
    // The point of the pair: the GPU is measurably looser, and by how much is the finding.
    println!(
        "  so the accelerator is about {:.0}x looser than the domain it accelerates",
        (gpu_drift / cpu_drift.max(1e-16)).max(1.0)
    );
}

/// **The stability limit is refused on the GPU exactly as on the CPU.**
#[test]
fn past_the_limit_is_refused() {
    let Some((cpu, mut gpu)) = pair() else {
        return;
    };
    let limit = cpu.max_stable_dt(Time::from_si(0.0));
    assert!(
        (gpu.max_stable_dt(Time::from_si(0.0)).to_si() / limit.to_si() - 1.0).abs() < 1e-12,
        "the two report the same limit; the scheme is the same and only the precision differs"
    );
    let err = gpu
        .step(
            Time::from_si(0.0),
            Time::from_si(limit.to_si() * 1.05),
            &mut Exchange::new(),
        )
        .expect_err("5% past the limit must be refused");
    assert_eq!(err.quantity, "Fourier number");
}

/// **It is actually faster, and by how much is measured rather than assumed.**
///
/// Reported rather than asserted. A timing threshold on somebody else's machine fails for reasons
/// that have nothing to do with the code, and a GPU that is slower than a CPU on a small grid is a
/// true fact about small grids rather than a defect.
#[test]
fn how_much_faster() {
    let Some((mut cpu, mut gpu)) = pair() else {
        return;
    };
    seed(&mut cpu, &mut gpu);
    let dt = Time::from_si(cpu.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let steps = 400;
    let mut bus = Exchange::new();

    let start = std::time::Instant::now();
    for _ in 0..steps {
        cpu.step(Time::from_si(0.0), dt, &mut bus).expect("stable");
    }
    let cpu_time = start.elapsed().as_secs_f64();

    let start = std::time::Instant::now();
    for _ in 0..steps {
        gpu.step(Time::from_si(0.0), dt, &mut bus).expect("stable");
    }
    // Force the queue to drain, or this times the submission and not the work.
    let _ = gpu.mean_temperature();
    let gpu_time = start.elapsed().as_secs_f64();

    println!(
        "  {}^3 cells, {steps} steps: cpu {:.3} s, gpu {:.3} s — {:.2}x",
        N,
        cpu_time,
        gpu_time,
        cpu_time / gpu_time.max(1e-9)
    );
    println!("  a readback is one transfer, so a run that audits every step pays for it");
}
