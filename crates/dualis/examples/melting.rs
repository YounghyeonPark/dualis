//! An argon crystal melting, watched through its own structure.
//!
//! ```text
//! cargo run --release --example melting            # numbers, checked
//! cargo run --release --example melting out.svg    # and a picture
//! ```
//!
//! Nothing here decides that the solid has melted. The radial distribution function is
//! accumulated at three temperatures and the shape says which phase it is — sharp spikes at the
//! lattice shells for a crystal, two or three broad humps for a liquid, one bump at contact for
//! a gas. That is also what a neutron diffraction experiment measures, so this is one of the
//! few places a simulation and a real instrument produce the same curve.
//!
//! The crystal panel is the one that can be checked exactly. A face-centred cubic lattice has
//! neighbour shells at `1 : √2 : √3 : 2` times the nearest-neighbour distance, holding 12, 6, 24
//! and 12 atoms — combinatorics rather than measurement — and the histogram reproduces all four
//! before anything is asked of it at a temperature where the answer is not known.
//!
//! Release mode on purpose. This is the first example that is genuinely expensive: five hundred
//! atoms, three temperatures, and tens of thousands of force evaluations each.

use dualis::prelude::*;
use dualis_molecular::{fcc_shells, RadialDistribution};

mod common;
use common::svg::{document, rgb, ticks, Plot};
use common::{check, check_between, heading};

/// Reduced units throughout: `ε = σ = m = 1`, which is how the literature quotes state points.
const DENSITY: f64 = 0.8442;
/// Five cells a side, so 500 atoms — enough for a box of 8.4 σ and a `g(r)` reaching 4.2.
const CELLS: usize = 5;

fn kelvin(reduced: f64) -> Temperature {
    dualis_molecular::temperature_from_reduced(reduced, &LennardJones::reduced())
}

fn fluid_at(reduced_t: f64, density: f64, seed: u64) -> Fluid {
    let target = kelvin(reduced_t);
    Fluid::lattice(
        "argon",
        LennardJones::reduced(),
        dualis_molecular::unit_mass(),
        CELLS,
        density,
    )
    .thermalised(target, seed)
    .with_thermostat(Thermostat::Langevin {
        target,
        damping: 1.0,
    })
}

/// One state point: equilibrate, then accumulate structure while sampling the thermodynamics.
struct Sample {
    label: &'static str,
    reduced_t: f64,
    density: f64,
    curve: Vec<(f64, f64)>,
    peak: (f64, f64),
    coordination: f64,
    settled_t: f64,
    pressure: f64,
}

fn measure(label: &'static str, reduced_t: f64, density: f64, seed: u64) -> Sample {
    let mut fluid = fluid_at(reduced_t, density, seed);
    let dt = fluid.max_stable_dt(Time::ZERO);
    let mut bus = Exchange::new();
    for _ in 0..3000 {
        fluid
            .step(Time::ZERO, dt, &mut bus)
            .expect("inside the box");
    }
    let mut rdf = RadialDistribution::new(&fluid, 4.5, 180);
    let (mut t_sum, mut p_sum, mut samples) = (0.0, 0.0, 0.0);
    for k in 0..6000 {
        fluid.step(Time::ZERO, dt, &mut bus).unwrap();
        if k % 10 == 0 {
            rdf.accumulate(&fluid);
        }
        t_sum += fluid.temperature().to_si();
        p_sum += fluid.pressure().to_si();
        samples += 1.0;
    }
    let lj = LennardJones::reduced();
    Sample {
        label,
        reduced_t,
        density,
        curve: rdf.curve(),
        peak: rdf.first_peak(),
        // Out to the first minimum of a dense fluid, which is where "first shell" ends.
        coordination: rdf.coordination(1.55),
        settled_t: dualis_molecular::reduced_temperature(
            Temperature::from_si(t_sum / samples),
            &lj,
        ),
        pressure: p_sum / samples / lj.epsilon * lj.sigma.powi(3),
    }
}

fn main() {
    heading("The crystal, whose answer is combinatorics");
    // A cold lattice, checked against the fcc shells before anything harder is attempted. Four
    // atoms to a cell fixes the edge at (4/rho)^(1/3), and everything else follows from it.
    let cold = Fluid::lattice(
        "crystal",
        LennardJones::reduced(),
        dualis_molecular::unit_mass(),
        CELLS,
        DENSITY,
    );
    let edge = (4.0 / DENSITY).cbrt();
    let shells = fcc_shells(edge);
    let mut exact = RadialDistribution::new(&cold, shells[3].0 * 1.05, 700);
    exact.accumulate(&cold);
    let mut running = 0.0;
    for (radius, count) in shells {
        running += count as f64;
        let midpoint = radius * 1.02;
        check(
            &format!("neighbours within {:.3} sigma", midpoint),
            exact.coordination(midpoint),
            running,
            1e-12,
            "",
        );
    }
    check(
        "shell ratio, 2nd over 1st",
        shells[1].0 / shells[0].0,
        2f64.sqrt(),
        1e-12,
        "",
    );
    check(
        "shell ratio, 4th over 1st",
        shells[3].0 / shells[0].0,
        2.0,
        1e-12,
        "",
    );

    heading("Three state points, each equilibrated and then sampled");
    let samples = [
        measure("crystal, T* = 0.15", 0.15, DENSITY, 0x_C01D_0001),
        measure("liquid,  T* = 0.85", 0.85, DENSITY, 0x_119A_1D02),
        measure("gas,     T* = 3.00", 3.00, 0.05, 0x_9A50_0003),
    ];
    for s in &samples {
        println!(
            "  {:<20} rho* = {:<7.4} T* settled {:.3}   P* {:>8.3}   g peak {:.2} at {:.3} sigma   first shell {:.1}",
            s.label, s.density, s.settled_t, s.pressure, s.peak.1, s.peak.0, s.coordination
        );
    }

    heading("What the structure says, without anyone deciding");
    let (crystal, liquid, gas) = (&samples[0], &samples[1], &samples[2]);

    // Every phase peaks where the potential does, because that is where an atom wants to sit.
    let well = LennardJones::reduced().minimum();
    for s in &samples[..2] {
        check_between(
            &format!("{}: first peak, in sigma", s.label.trim_end()),
            s.peak.0,
            well - 0.2,
            well + 0.25,
            "",
        );
    }

    // The height of that peak is the phase. A lattice concentrates its neighbours into one
    // shell; a liquid smears them over a range; a gas has none to concentrate.
    check_between("crystal: peak height", crystal.peak.1, 4.0, 7.0, "x");
    check_between("liquid:  peak height", liquid.peak.1, 2.3, 3.5, "x");
    check_between("gas:     peak height", gas.peak.1, 1.1, 1.9, "x");
    assert!(
        crystal.peak.1 > liquid.peak.1 && liquid.peak.1 > gas.peak.1,
        "sharpness must fall as the order does"
    );

    // But the *number* of first neighbours barely moves between solid and liquid: melting costs
    // the order, not the packing. That is the whole reason a liquid is nearly as dense as its
    // solid while a gas is a thousand times thinner.
    check_between(
        "crystal: first-shell count",
        crystal.coordination,
        11.0,
        13.5,
        "",
    );
    check_between(
        "liquid:  first-shell count",
        liquid.coordination,
        11.5,
        14.0,
        "",
    );
    assert!(
        (crystal.coordination - liquid.coordination).abs() < 2.0,
        "melting should keep the neighbours: {} against {}",
        crystal.coordination,
        liquid.coordination
    );
    check_between("gas:     first-shell count", gas.coordination, 0.3, 1.2, "");

    // Long-range order is what melting destroys, and the *third* shell is where to look. Not
    // the fourth, which was the first thing tried and reads 1.045 for the crystal against 0.999
    // for the liquid -- no discrimination at all. The reason is worth knowing: an fcc lattice's
    // shells hold 12, 6, 24 and 12 atoms, and thermal broadening washes out a thin shell at a
    // large radius long before a fat one at a smaller radius. The third is the most populous of
    // the four, so it is the last to go.
    let g_at = |s: &Sample, r: f64| {
        s.curve
            .iter()
            .filter(|(x, _)| (*x - r).abs() < 0.1)
            .map(|(_, g)| *g)
            .fold(0.0f64, f64::max)
    };
    let (third, fourth) = (shells[2].0, shells[3].0);
    check_between(
        &format!("crystal: g at the 3rd shell, {third:.2} sigma"),
        g_at(crystal, third),
        2.0,
        4.0,
        "x",
    );
    check_between(
        "liquid:  g at the same radius",
        g_at(liquid, third),
        1.0,
        1.6,
        "x",
    );
    assert!(
        g_at(crystal, third) > 1.8 * g_at(liquid, third),
        "the third shell should separate them clearly"
    );
    // And the fourth genuinely does not, which is worth stating rather than quietly not testing.
    println!(
        "  the 4th shell, by contrast: crystal {:.3}, liquid {:.3} -- indistinguishable",
        g_at(crystal, fourth),
        g_at(liquid, fourth)
    );

    heading("And the thermodynamics agrees it is the same fluid");
    // The thermostat landed each run on its target, which is the only reason the phases can be
    // compared at all.
    for s in &samples {
        check(
            &format!("{}: settled temperature", s.label.trim_end()),
            s.settled_t,
            s.reduced_t,
            0.06,
            "T*",
        );
    }
    // The dilute gas is nearly ideal: P* = rho* T* to a few percent.
    let ideal = gas.density * gas.settled_t;
    check("gas: P* against rho* T*", gas.pressure, ideal, 0.06, "");
    // The cold crystal is under tension -- its neighbours sit past the well's minimum, so they
    // pull inwards and the pressure is negative. A gas can never do that.
    //
    // Asserted as a band that excludes zero, rather than as `pressure.max(0.0)` against a
    // scale, which is how this was written first. That form says "the pressure is not
    // positive", and a virial that came out identically zero -- from a cutoff that caught no
    // pairs, say, which is a mistake this workspace has actually made -- satisfies it
    // perfectly. The tension has to be *there*, not merely not-absent.
    //
    // The bounds are the physics and not the measurement: at this density every particle has
    // twelve neighbours a little past r_min, so the tension is some units of epsilon per
    // sigma^3, and a lattice that had melted or exploded would leave the band on one side or
    // the other. Measured: -4.82.
    check_between(
        "crystal: pressure is negative",
        crystal.pressure,
        -8.0,
        -1.0,
        "P*",
    );
    println!(
        "  crystal P* = {:.3}, so the lattice is pulling itself together",
        crystal.pressure
    );

    let Some(path) = common::output_path() else {
        println!(
            "\npass a path to write an SVG, e.g. `cargo run --release --example melting out.svg`"
        );
        return;
    };
    common::write(&path, &draw(&samples, &shells));
}

/// Three stacked panels, one per phase, on a shared radius axis.
fn draw(samples: &[Sample], shells: &[(f64, usize); 4]) -> String {
    let (w, h) = (880.0, 560.0);
    let top = samples
        .iter()
        .flat_map(|s| s.curve.iter().map(|(_, g)| *g))
        .fold(0.0f64, f64::max)
        .min(9.0)
        * 1.08;
    let far = samples
        .iter()
        .flat_map(|s| s.curve.iter().map(|(r, _)| *r))
        .fold(0.0f64, f64::max);

    let colours = [rgb(48, 108, 186), rgb(158, 40, 32), rgb(120, 120, 120)];
    let mut parts = Vec::new();
    for (k, s) in samples.iter().enumerate() {
        let y = 56.0 + k as f64 * 158.0;
        let mut panel = Plot::new(w, h, (0.0, far), (0.0, top)).viewport(72.0, y, 760.0, 118.0);
        panel.axes(
            &ticks(0.0, far, 9),
            &[0.0, 1.0, 3.0, 6.0, 9.0],
            |v| format!("{v:.0}"),
            |v| format!("{v:.0}"),
        );
        // One is the ideal gas, and where a curve sits relative to it is the whole reading.
        panel.polyline([(0.0, 1.0), (far, 1.0)], "#00000040", 1.0);
        // The lattice shells, so the crystal's spikes can be seen landing on them.
        for (radius, _) in shells.iter() {
            if *radius < far {
                panel.polyline([(*radius, 0.0), (*radius, top)], "#00000018", 1.0);
            }
        }
        panel.polyline(s.curve.iter().copied(), &colours[k], 1.9);
        panel.label(
            76.0,
            y - 8.0,
            &format!(
                "{}   rho* {:.4}   peak {:.1} at {:.2} sigma   first shell {:.1} neighbours",
                s.label.trim_end(),
                s.density,
                s.peak.1,
                s.peak.0,
                s.coordination
            ),
            12.0,
            "#3a3a3a",
            "start",
        );
        parts.push(panel.into_body());
    }

    let mut frame = Plot::new(w, h, (0.0, 1.0), (0.0, 1.0));
    frame.label(
        72.0,
        26.0,
        "A Lennard-Jones fluid through melting, read off its own radial distribution",
        15.0,
        "#1b1b1b",
        "start",
    );
    frame.label(
        808.0,
        26.0,
        "500 atoms   grey lines: the fcc neighbour shells",
        12.0,
        "#6a6a6a",
        "end",
    );
    frame.label(
        452.0,
        h - 14.0,
        "separation in sigma  --  the flat line at one is an ideal gas, and distance from it is structure",
        11.5,
        "#6a6a6a",
        "middle",
    );
    parts.push(frame.into_body());

    document(w, h, parts)
}
