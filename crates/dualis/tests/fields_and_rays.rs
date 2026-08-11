//! Where the wave picture and the ray picture have to agree.
//!
//! `dualis-optics` computes a surface's reflectance from Fresnel's equations — algebra on two
//! refractive indices. `dualis-em` marches Maxwell on a grid and knows nothing about either. The
//! two crates share no code at all; they do not even depend on each other.
//!
//! So a reflectance that comes out the same from both is a real check, and it is the only kind of
//! check that says the two are limits of one physics rather than two unrelated models that happen
//! to be in the same workspace.
//!
//! # What makes a plane wave statable in a box
//!
//! A wave along `z` polarised along `y` has `Ey` **tangential** to the `x` faces, so a conducting
//! wall there would force it to zero and turn the plane wave into a waveguide mode — with a cutoff,
//! a different phase velocity, and a reflectance that is a function of both. Conductors on the `y`
//! faces and *magnetic* conductors on the `x` ones leave the wave exactly uniform, and the `z`
//! faces absorb so nothing comes back to be counted twice.

use dualis::em::{Boundary, Cavity, Medium, Wall};
use dualis::optics::fresnel_reflectance;
use dualis::prelude::*;

/// Cells per wavelength, where a test does not vary it.
///
/// Forty rather than twenty, and the difference is not cosmetic: at twenty the reflectance of a
/// glass surface comes out **10.5%** low, at forty 2.6% and at eighty 0.63%. That is the
/// second-order convergence `the_grid_and_the_algebra_agree` measures, and it is why a single
/// resolution with a loose tolerance would have been a worse test than a tight one at three.
const PER_WAVELENGTH: usize = 40;

/// A one-dimensional plane-wave testbed in a three-dimensional grid.
///
/// **Two cells across and one deep**, which is the smallest grid that still has an interior. Every
/// transverse cell holds the same field — that is what the boundary conditions are for — so a
/// 4 x 4 cross-section is sixteen copies of one answer and costs sixteen times as much. It was
/// 4 x 4 first, and this file took **486 s** in a debug build; the four jobs of CI's OS matrix run
/// debug builds.
///
/// Two rather than one because `advance_electric` updates the interior `x` faces, `1..nx`, and at
/// `nx = 1` that range is empty: the magnetic walls would mirror a field nothing ever wrote.
fn testbed(cells_along_z: usize, dx: f64) -> Cavity {
    let mut c = Cavity::new(
        "slab",
        (2, 1, cells_along_z),
        Length::from_si(dx),
        Medium::vacuum(),
    );
    c.set_boundary(Wall::XLow, Boundary::Magnetic);
    c.set_boundary(Wall::XHigh, Boundary::Magnetic);
    c.set_boundary(Wall::ZLow, Boundary::Open);
    c.set_boundary(Wall::ZHigh, Boundary::Open);
    c
}

/// A Gaussian-modulated pulse: narrow enough in time to separate the incident and reflected
/// passes at the monitor, wide enough in frequency to be centred on a wavelength the grid carries.
fn pulse(centre: f64, width: f64, wavelength: f64) -> impl Fn(f64, f64, f64) -> f64 {
    let k = 2.0 * std::f64::consts::PI / wavelength;
    move |_x: f64, _y: f64, z: f64| {
        let u = (z - centre) / width;
        (-u * u).exp() * (k * (z - centre)).cos()
    }
}

/// Run and record `Ey` at a monitor plane, returning the time series.
fn watch(c: &mut Cavity, monitor: usize, steps: usize, dt: Time) -> Vec<f64> {
    let mut bus = Exchange::new();
    let mut out = Vec::with_capacity(steps);
    for n in 0..steps {
        c.step(Time::from_si(n as f64 * dt.to_si()), dt, &mut bus)
            .expect("stable");
        out.push(c.electric_at(1, 0, monitor).y);
    }
    out
}

/// The energy in a window of the series, which is what a reflectance is a ratio of.
fn energy(series: &[f64], from: usize, to: usize) -> f64 {
    series[from.min(series.len())..to.min(series.len())]
        .iter()
        .map(|v| v * v)
        .sum()
}

/// Where a pulse's energy sits in time, `Σ t·E² / Σ E²`.
///
/// The centroid rather than the peak, because the peak of an oscillating pulse is a peak of its
/// **carrier** and moves in steps of half a period — 0.83 fs here, against delays of a few
/// femtoseconds. The first version of `a_slab_delays_a_pulse` used the peak and read 9% high for
/// that reason and no other.
fn centroid(series: &[f64], from: usize, to: usize, dt: f64) -> f64 {
    let slice = &series[from.min(series.len())..to.min(series.len())];
    let total: f64 = slice.iter().map(|v| v * v).sum();
    if total <= 0.0 {
        return 0.0;
    }
    let weighted: f64 = slice
        .iter()
        .enumerate()
        .map(|(i, v)| (from + i) as f64 * v * v)
        .sum();
    weighted / total * dt
}

/// Where everything sits, in wavelengths, so the pulse and the monitor cannot overlap.
///
/// # The geometry is the measurement
///
/// The incident and reflected passes are separated at the monitor by `2(z_i − z_m)/c`, and the
/// pulse has to be shorter than that *and* clear of the monitor when it starts. The first version
/// put a pulse of half-width 1.6λ at 4λ and the monitor at 5λ — inside the pulse — so the
/// "incident" window opened with the pulse already sitting on it and the reflectance came out 34.7%
/// where Fresnel says 14.8%.
struct Layout {
    /// Where the pulse starts.
    source: f64,
    /// Its Gaussian half-width.
    width: f64,
    /// Where `Ey` is watched.
    monitor: f64,
    /// Where the dielectric begins.
    interface: f64,
    /// How long the box is.
    length: f64,
}

impl Layout {
    fn cells(&self, what: f64, dx: f64, wavelength: f64) -> usize {
        (what * wavelength / dx).round() as usize
    }
}

/// Set a slab up, run it, and return `(reflectance, the series' two windows)`.
fn reflectance(
    layout: &Layout,
    wavelength: f64,
    per_wavelength: usize,
    build: impl Fn(&mut Cavity, usize),
) -> f64 {
    let dx = wavelength / per_wavelength as f64;
    let nz = layout.cells(layout.length, dx, wavelength);
    let interface = layout.cells(layout.interface, dx, wavelength);
    let monitor = layout.cells(layout.monitor, dx, wavelength);

    let mut c = testbed(nz, dx);
    build(&mut c, interface);
    let dt = Time::from_si(c.courant_limit().to_si() * 0.5);
    c.launch_along_z(
        dt,
        pulse(
            layout.source * wavelength,
            layout.width * wavelength,
            wavelength,
        ),
    );

    let step_distance = 299_792_458.0 * dt.to_si();
    let there = ((layout.monitor - layout.source) * wavelength / step_distance) as usize;
    let back = ((2.0 * layout.interface - layout.monitor - layout.source) * wavelength
        / step_distance) as usize;
    let series = watch(&mut c, monitor, back + (back - there), dt);
    let split = (there + back) / 2;
    energy(&series, split, series.len()) / energy(&series, 0, split)
}

/// **FDTD reproduces Fresnel's reflectance, and neither crate knows the other exists.**
///
/// A pulse crosses a monitor plane going forward, meets a dielectric half-space, and crosses the
/// monitor again coming back. The two passes are separated in time by `2(z_interface − z_monitor)/c`,
/// so the reflectance is the ratio of the energy in the two windows — no fitting, no phase
/// unwrapping, just two integrals.
///
/// # The agreement is a convergence, and saying which matters
///
/// A grid does not reproduce Fresnel exactly at any resolution — the discrete dispersion and the
/// averaged permittivity at the interface both contribute, and both are second order. Measured at
/// `n = 2`, where Fresnel says 11.111%:
///
/// ```text
///   20 cells per wavelength    −10.543%
///   40                          −2.556%      ratio 4.13
///   80                          −0.634%      ratio 4.03
/// ```
///
/// So the claim is the *rate*, not a tolerance somebody chose. A single resolution with a loose
/// bound would pass for a scheme converging at first order, and first order is what an interface
/// placed half a cell wrong would give.
///
/// Then three indices at the finest, because one is a number and three is a curve: 1.5 gives 4%,
/// 2.0 gives 11.1% and 3.5 gives 30.9%.
#[test]
fn the_grid_and_the_algebra_agree_on_a_surfaces_reflectance() {
    let wavelength = 500e-9;
    // As tight as the separation allows: the two passes at the monitor are `2(14−8) = 12λ` apart
    // and the pulse is about 4λ long, so they do not touch. Every wavelength of box is cells and
    // steps at three resolutions, and the widest of them costs sixteen times the narrowest.
    let layout = Layout {
        source: 4.0,
        width: 1.0,
        monitor: 8.0,
        interface: 14.0,
        length: 20.0,
    };

    let closed = fresnel_reflectance(1.0, 2.0, 1.0);
    let mut errors = Vec::new();
    for per in [20usize, 40, 80] {
        let measured = reflectance(&layout, wavelength, per, |c, interface| {
            c.fill(Medium::dielectric(4.0), |_, _, k| k >= interface);
        });
        let error = (measured / closed - 1.0).abs();
        println!(
            "  {per:>2} cells per wavelength: {:.4}% against Fresnel's {:.4}% — off {:.3}%",
            measured * 100.0,
            closed * 100.0,
            error * 100.0
        );
        errors.push(error);
    }
    for pair in errors.windows(2) {
        let rate = pair[0] / pair[1];
        println!("  refinement ratio {rate:.2} (second order is 4)");
        assert!(
            (2.6..6.0).contains(&rate),
            "the agreement converges at second order: {rate:.3}"
        );
    }

    for n2 in [1.5, 2.0, 3.5] {
        let measured = reflectance(&layout, wavelength, 80, |c, interface| {
            c.fill(Medium::dielectric(n2 * n2), |_, _, k| k >= interface);
        });
        let closed = fresnel_reflectance(1.0, n2, 1.0);
        println!(
            "  n = {n2}: FDTD {:.4}% against Fresnel {:.4}% — off {:.3}%",
            measured * 100.0,
            closed * 100.0,
            (measured / closed - 1.0).abs() * 100.0
        );
        assert!(
            (measured / closed - 1.0).abs() < 0.02,
            "and lands on it at eighty cells: {measured:.5} against {closed:.5}"
        );
    }
}

/// **A quarter-wave layer of `√n` makes the reflection vanish**, which is interference and not an
/// interface.
///
/// The anti-reflection coating: a layer of index `√(n₁n₃)` and optical thickness `λ/4` reflects
/// **nothing** at that wavelength, because the two reflections it creates are equal in size and
/// opposite in phase. Fresnel alone cannot say this — it is a statement about two surfaces at a
/// particular separation, and it is exactly the sort of thing a field solver is for.
///
/// So this is the other half of the test above: there, the grid agreed with the algebra; here, it
/// says something the algebra by itself does not.
#[test]
fn a_quarter_wave_coating_cancels_its_own_reflection() {
    let wavelength = 500e-9;
    // A wider pulse than the test above, because a coating is a quarter wave at **one**
    // wavelength: a short pulse is a broad spectrum and most of it would not be cancelled.
    let layout = Layout {
        source: 10.0,
        width: 3.0,
        monitor: 22.0,
        interface: 34.0,
        length: 50.0,
    };
    let n_glass = 2.25f64;
    let n_coat = n_glass.sqrt();
    let dx = wavelength / PER_WAVELENGTH as f64;
    let thickness_cells = (wavelength / (4.0 * n_coat) / dx).round() as usize;

    let mut measured = Vec::new();
    for coated in [false, true] {
        let r = reflectance(&layout, wavelength, PER_WAVELENGTH, |c, interface| {
            c.fill(Medium::dielectric(n_glass * n_glass), |_, _, k| {
                k >= interface
            });
            if coated {
                c.fill(Medium::dielectric(n_coat * n_coat), |_, _, k| {
                    k >= interface - thickness_cells && k < interface
                });
            }
        });
        measured.push(r);
        println!(
            "  {:<8} reflectance {:.4}%",
            if coated { "coated:" } else { "bare:" },
            r * 100.0
        );
    }

    let bare_closed = fresnel_reflectance(1.0, n_glass, 1.0);
    println!(
        "  bare against Fresnel {:.4}%, and the coating removed {:.1}x of it",
        bare_closed * 100.0,
        measured[0] / measured[1].max(1e-12)
    );
    assert!(
        (measured[0] / bare_closed - 1.0).abs() < 0.05,
        "the bare surface is the Fresnel one: {:.5} against {bare_closed:.5}",
        measured[0]
    );
    assert!(
        measured[1] < measured[0] / 10.0,
        "a quarter wave of the geometric mean must very nearly cancel it: {:.5} against {:.5}",
        measured[1],
        measured[0]
    );
}

/// **Light slows by exactly the index**, which is the one thing the two pictures share by
/// definition rather than by derivation.
///
/// `n = √εᵣ` in the medium constructor, and `c/n` out of the grid. Measured from the delay a slab
/// adds to a pulse's arrival, so it is a transit time rather than a restatement of the constructor.
#[test]
fn a_slab_delays_a_pulse_by_exactly_its_index() {
    let wavelength = 500e-9;
    let dx = wavelength / PER_WAVELENGTH as f64;
    let nz = 40 * PER_WAVELENGTH;
    let n = 2.0f64;
    let slab = (10 * PER_WAVELENGTH, 20 * PER_WAVELENGTH);
    let monitor = 32 * PER_WAVELENGTH;

    let mut arrivals = Vec::new();
    for with_slab in [false, true] {
        let mut c = testbed(nz, dx);
        if with_slab {
            c.fill(Medium::dielectric(n * n), |_, _, k| {
                k >= slab.0 && k < slab.1
            });
        }
        let dt = Time::from_si(c.courant_limit().to_si() * 0.5);
        c.launch_along_z(dt, pulse(4.0 * wavelength, 1.2 * wavelength, wavelength));
        // Long enough for the slowest case to arrive and pass, and no longer: after that the
        // internal reflections trickle out and would drag the centroid.
        let step_distance = 299_792_458.0 * dt.to_si();
        let steps = (1.4 * (monitor as f64 * dx) / step_distance) as usize;
        let series = watch(&mut c, monitor, steps, dt);
        arrivals.push(centroid(&series, 0, series.len(), dt.to_si()));
    }

    let extra = arrivals[1] - arrivals[0];
    let thickness = (slab.1 - slab.0) as f64 * dx;
    let closed = thickness * (n - 1.0) / 299_792_458.0;
    println!(
        "  a {:.2} um slab of n = {n} delayed the pulse {:.4} fs against (n-1)d/c = {:.4} fs",
        thickness * 1e6,
        extra * 1e15,
        closed * 1e15
    );
    assert!(
        (extra / closed - 1.0).abs() < 0.06,
        "the delay is (n-1)d/c: {:.4e} against {closed:.4e}",
        extra
    );
}
