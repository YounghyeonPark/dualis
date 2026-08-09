//! Where the time actually goes, before anyone optimises anything.
//!
//! This workspace has never measured itself. That is defensible — it makes no speed claims —
//! but it means any conversation about SIMD or a GPU starts from a guess about which loop
//! matters, and a guess is how effort goes into the wrong one.
//!
//! Dependency-free on purpose. A benchmark harness is a public commitment and a licence to
//! audit, for a number that is only ever read by a person deciding what to work on.
//!
//! ```sh
//! cargo run --release --example where_the_time_goes
//! ```
//!
//! **It uses a clock, which nothing else in this workspace may do.** That is why it is an
//! example rather than a test: the determinism rule is about results, and a timing that varies
//! run to run is the one number here that is allowed to. Nothing it prints is asserted.

use dualis::prelude::*;
use std::time::Instant;

/// Run `f` until at least `seconds` have passed, and report the per-step cost.
fn time_it(label: &str, work: &str, count: usize, mut f: impl FnMut()) -> f64 {
    // One untimed pass, so a first-touch page fault or a lazily built cell list is not counted
    // as the steady-state cost of a step.
    f();
    let start = Instant::now();
    let mut steps = 0u64;
    while start.elapsed().as_secs_f64() < 0.6 {
        f();
        steps += 1;
    }
    let per_step = start.elapsed().as_secs_f64() / steps as f64;
    let per_unit = per_step / count as f64;
    println!(
        "  {label:<24} {:>9.3} ms/step   {:>8.1} ns per {work}   ({count} of them)",
        per_step * 1e3,
        per_unit * 1e9
    );
    per_step
}

fn main() {
    println!("where the time goes — one step of each, release build\n");
    let dt = Time::ms(1.0);

    // ---- Molecular: the densest inner loop in the workspace.
    println!("dualis-molecular");
    let mut totals = Vec::new();
    for cells in [4usize, 6, 8] {
        let atoms = 4 * cells * cells * cells;
        let mut fluid = Fluid::lattice(
            "argon",
            LennardJones::reduced(),
            dualis_molecular::unit_mass(),
            cells,
            0.85,
        )
        .thermalised(Temperature::from_si(1.0), 7);
        let mut bus = Exchange::new();
        let t = time_it(&format!("Fluid, {atoms} atoms"), "atom", atoms, || {
            fluid
                .step(Time::from_si(1e-3), Time::from_si(1e-3), &mut bus)
                .unwrap();
            bus.take(quantity::ENERGY);
        });
        totals.push(("fluid", atoms, t));
    }

    // ---- Mechanics: O(N²) against Barnes-Hut. Both sizes on purpose: the tree *loses* at
    // 1024 and wins at 4096, so a single size would say the wrong thing about it. Measured
    // crossover is near N = 2000, and at 16384 the tree is 3.6x ahead.
    println!("\ndualis-mechanics");
    for n in [1024usize, 4096] {
        let bodies: Vec<Body> = (0..n)
            .map(|i| {
                let mut r = Rng::for_index(11, i as u64);
                Body::new(
                    Mass::kg(1e20),
                    LengthVec::m(
                        (r.unit() - 0.5) * 1e10,
                        ({
                            let mut q = Rng::for_index(12, i as u64);
                            q.unit()
                        } - 0.5)
                            * 1e10,
                        ({
                            let mut q = Rng::for_index(13, i as u64);
                            q.unit()
                        } - 0.5)
                            * 1e10,
                    ),
                    VelocityVec::m_per_s(0.0, 0.0, 0.0),
                )
            })
            .collect();
        let mut direct = NBody::new("direct", &bodies);
        let mut bus = Exchange::new();
        time_it(&format!("NBody, {n} bodies"), "body", n, || {
            direct.step(Time::s(0.0), Time::s(1.0), &mut bus).unwrap();
        });
        let mut tree = TreeNBody::new("tree", &bodies);
        time_it(&format!("TreeNBody, {n} bodies"), "body", n, || {
            tree.step(Time::s(0.0), Time::s(1.0), &mut bus).unwrap();
        });
    }

    // ---- Acoustics: a 2D stencil, which is the shape SIMD is usually best at.
    println!("\ndualis-acoustic");
    for cells in [201usize, 401] {
        let mut room = Room::of_air("room", Length::m(4.4), Length::m(3.1), cells);
        let (nx, ny) = room.cells();
        let n = nx * ny;
        let mut bus = Exchange::new();
        let h = room.max_stable_dt(Time::s(0.0));
        time_it(&format!("Room, {cells} across"), "cell", n, || {
            room.step(Time::s(0.0), h, &mut bus).unwrap();
        });
    }

    // ---- Thermal: the 1D stencil and the dense solve, for scale.
    println!("\ndualis-thermal");
    let mut bar = Bar1D::new(
        "bar",
        Substance::aluminium_6061(),
        4001,
        Length::mm(0.005),
        Area::mm2(100.0),
        Temperature::celsius(20.0),
    );
    let mut bus = Exchange::new();
    let h = bar.max_stable_dt(Time::s(0.0));
    time_it("Bar1D, 4001 cells", "cell", 4001, || {
        bar.step(Time::s(0.0), h, &mut bus).unwrap();
    });
    let _ = dt;

    println!(
        "\nread it as: which loop would have to get faster for a run to get faster.\n\
         nothing here is asserted, and nothing here is a promise."
    );
}
