//! A laser on a mirror, and the hot spot a lumped model cannot see.
//!
//! ```text
//! cargo run --example beam_hot_spot            # numbers, checked
//! cargo run --example beam_hot_spot out.svg    # and a picture
//! ```
//!
//! Two domains that have never heard of each other. `dualis-optics` integrates a 100 W green
//! laser against a mirror's spectral absorptance and gets a real number of watts;
//! `dualis-thermal` conducts them through aluminium. They meet over a boundary cut into
//! faces both address, so the heat arrives *where the beam is* rather than as a total —
//! and the kernel audits the crossing face by face.
//!
//! The point of the picture is the gap between two curves. The mean rise is what a lumped
//! model predicts, and it is exactly right: `Q/C`, with insulated ends. The peak is four
//! times it, sitting under a 1 mm beam on a 20 mm bar, and a coating fails at its peak.

use dualis::prelude::*;
use dualis_optics::spectrum::Spectrum as Spec;

mod common;
use common::svg::{heat, rgb, ticks, Plot};
use common::{check, check_between, check_zero, heading};

/// 20 mm of aluminium, side-on to the beam, cut into 81 cells.
const SPAN: f64 = 20e-3;
const CELLS: usize = 81;
const CROSS_SECTION: f64 = 1e-4;
/// 1/e² irradiance radius. A twentieth of the bar, so the beam is genuinely a spot.
const WAIST: f64 = 1e-3;

/// The optics side: absorbed watts, published with the profile they landed with.
struct Beam {
    lamp: SpectralPower,
    absorptance: Spectrum,
    boundary: Interface,
    paid_out: f64,
}

impl Beam {
    fn new(boundary: &Interface) -> Beam {
        let mirror = SurfaceOptics::mirror(0.99);
        Beam {
            lamp: SpectralPower::new(
                Spec::gaussian(Length::nm(532.0), Length::nm(1.0), 1.0),
                Power::w(100.0),
                (Length::nm(400.0), Length::nm(700.0)),
            ),
            absorptance: Spec::curve(
                (0..=60)
                    .map(|i| {
                        let nm = 400.0 + i as f64 * 5.0;
                        (nm, mirror.absorptance(Length::nm(nm)))
                    })
                    .collect(),
            ),
            boundary: boundary.clone(),
            paid_out: 0.0,
        }
    }

    fn absorbed_power(&self) -> Power {
        self.lamp.absorbed_by(&self.absorptance)
    }
}

impl Domain for Beam {
    fn name(&self) -> &'static str {
        "beam"
    }
    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }
    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let joules = self.absorbed_power().to_si() * dt.to_si();
        // exp(-2r²/w²), written the way the physics is written. The kernel normalises it so
        // the faces sum to `joules` exactly, and leaves the shape alone.
        let flux = Flux::profiled(joules, &self.boundary, |u| {
            (-2.0 * ((u - 0.5) * SPAN / WAIST).powi(2)).exp()
        });
        bus.publish_on(&self.boundary, HEAT, &flux)?;
        self.paid_out += joules;
        Ok(())
    }
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, -self.paid_out)
    }
    fn checkpoint(&mut self) {}
    fn restore(&mut self) {}
    fn supports_restore(&self) -> bool {
        true
    }
}

fn bar() -> Bar1D {
    let dx = Length::from_si(SPAN / CELLS as f64);
    Bar1D::new(
        "bar",
        Substance::aluminium_6061(),
        CELLS,
        dx,
        Area::from_si(CROSS_SECTION),
        Temperature::celsius(20.0),
    )
    .exposing("bar face", Area::from_si(10e-3 * dx.to_si()))
}

/// One snapshot: the time it was taken, the sampled profile in (mm, °C), and the mean and
/// peak rise above ambient.
type Frame = (f64, Vec<(f64, f64)>, f64, f64);

/// Read the bar through [`ScalarField`] rather than through its own accessors.
///
/// Deliberately: the trait is the interface a visualiser sees a simulation through, and
/// nothing had implemented it until now. Sampling at 400 points across a bar of 81 cells
/// exercises the interpolation, and the fact that this function never mentions `Bar1D`
/// is the property that matters — the same code will draw any scalar field.
fn sample(field: &dyn ScalarField, samples: usize) -> Vec<(f64, f64)> {
    (0..samples)
        .map(|i| {
            let x = SPAN * i as f64 / (samples - 1) as f64;
            let t = field.at(
                LengthVec::from_si(glam::DVec3::new(x, 0.0, 0.0)),
                Time::ZERO,
            );
            (x * 1e3, t - 273.15)
        })
        .collect()
}

fn main() {
    let plate = bar();
    let boundary = plate.boundary().expect("the bar exposes a face").clone();
    let beam = Beam::new(&boundary);
    let absorbed = beam.absorbed_power().to_si();

    let mut sim = Simulation::new(Schedule::Multirate).with(beam).with(plate);

    heading("A 100 W laser on a 99% mirror, 20 mm of aluminium behind it");
    // 1% of 100 W, from a real spectral integral rather than from the number 0.01.
    check("absorbed power", absorbed, 1.0, 0.06, "W");

    // Snapshots on a roughly logarithmic ladder: the spot is sharpest early and the profile
    // stops changing shape once conduction has crossed the bar.
    let stops = [0.01f64, 0.05, 0.2, 1.0, 5.0];
    let window = Time::ms(0.5);
    let mut frames: Vec<Frame> = Vec::new();
    let mut elapsed = 0.0;
    for stop in stops {
        while elapsed < stop - 1e-12 {
            sim.advance(window).expect("energy must be conserved");
            elapsed += window.to_si();
        }
        let plate: &Bar1D = sim.domain_as("bar").expect("the bar is in the simulation");
        let mean = plate.mean_temperature().in_celsius() - 20.0;
        let peak = plate.temperature_at(CELLS / 2).in_celsius() - 20.0;
        frames.push((stop, sample(plate, 400), mean, peak));
    }

    let capacity = Substance::aluminium_6061()
        .heat_capacity(Volume::from_si(CROSS_SECTION * SPAN))
        .unwrap()
        .to_si();

    heading("The lumped answer is exact, and it is not the answer you need");
    let (_, _, mean_early, peak_early) = frames[1];
    check(
        "mean rise at 50 ms  (closed form Q/C)",
        mean_early,
        absorbed * 0.05 / capacity,
        1e-9,
        "K",
    );
    check_between("peak rise at 50 ms", peak_early, 0.03, 0.06, "K");
    check_between(
        "peak / mean at 50 ms  (concentrated)",
        peak_early / mean_early,
        3.0,
        7.0,
        "x",
    );
    let (_, _, mean_late, peak_late) = frames[4];
    // Not 1.0, and never will be: the excess settles to a fixed number of kelvin while the
    // mean climbs forever, so it is the *ratio* that decays. tests/beam_heats_where_it_lands
    // pins the settled excess against a quadrature solution; here it only has to be small.
    check_between("peak / mean at 5 s", peak_late / mean_late, 1.0, 1.2, "x");

    heading("What the field reports, read through ScalarField");
    let plate: &Bar1D = sim.domain_as("bar").unwrap();
    let at = |mm: f64| LengthVec::from_si(glam::DVec3::new(mm * 1e-3, 0.0, 0.0));
    let h = Length::from_si(SPAN / CELLS as f64);
    let mm_of = |i: usize| (i as f64 + 0.5) * SPAN * 1e3 / CELLS as f64;

    // Symmetry, which costs nothing to check and would catch an off-by-one in the stencil.
    check(
        "dT/dx at the beam centre  (symmetry)",
        plate.gradient(at(10.0), Time::ZERO, h).x,
        0.0,
        1e-9,
        "K/m",
    );
    // Off the end the field is flat, because the ends are insulated and the temperature
    // stops changing there. Not the same claim as "zero slope at the wall", which would be
    // false — the first cell centre is half a cell inside, where it is still changing.
    check(
        "dT/dx beyond the end of the bar",
        plate.gradient(at(-5.0), Time::ZERO, h).x,
        0.0,
        1e-9,
        "K/m",
    );
    // The steepest slope is somewhere between the beam and the end, not at either — the
    // shape of the answer, and a stencil bug would move it. Measured as a distance from the
    // centre because there are two of them, one on each side.
    let steepest = (0..CELLS)
        .max_by(|a, b| {
            let g = |i: &usize| plate.gradient(at(mm_of(*i)), Time::ZERO, h).x.abs();
            g(a).total_cmp(&g(b))
        })
        .unwrap();
    check_between(
        "steepest gradient, distance from the beam",
        (mm_of(steepest) - 10.0).abs(),
        0.5,
        5.0,
        "mm",
    );

    // What the beam puts into its own cell, against what conduction takes straight back out.
    // Most of it leaves sideways: that is why the peak settles instead of running away.
    let cell_capacity = capacity / CELLS as f64;
    let centre_watts = Flux::profiled(absorbed, &boundary, |u| {
        (-2.0 * ((u - 0.5) * SPAN / WAIST).powi(2)).exp()
    })
    .at(CELLS / 2);
    let source = centre_watts / cell_capacity;
    let conduction = plate.rate(at(10.0), Time::ZERO, window);
    println!("  absorbed into the centre cell         {source:>12.4} K/s");
    println!("  conducted back out of it              {conduction:>12.4} K/s");
    check_between(
        "  fraction leaving sideways",
        -conduction / source,
        0.7,
        0.95,
        "",
    );

    // **The exact one.** Conduction moves heat about; it does not make any. Summed over an
    // insulated bar the reported rates cancel to the last bit, because the mirrored stencil
    // telescopes — every cell's gain is its neighbour's loss, and the end cells have no
    // neighbour to lose to. A boundary handled wrongly would show up here and nowhere else.
    let net: f64 = (0..CELLS)
        .map(|i| plate.rate(at(mm_of(i)), Time::ZERO, window))
        .sum();
    let largest = (0..CELLS)
        .map(|i| plate.rate(at(mm_of(i)), Time::ZERO, window).abs())
        .fold(0.0f64, f64::max);
    println!("  net conduction over the whole bar     {net:>12.3e} K/s   (largest single cell {largest:.3})");
    assert!(
        net.abs() < largest * 1e-12,
        "conduction must not create heat: net {net} against a scale of {largest}"
    );

    heading("Conservation, audited face by face");
    let residual = sim.ledger().get(quantity::ENERGY).unwrap();
    // Judged against the joules that actually crossed, which is the only scale that means
    // anything here: the books balance, so the total is zero by design.
    check_zero("energy residual", residual, absorbed * 5.0, 1e-9, "J");
    assert!(
        sim.bus().unclaimed().next().is_none(),
        "nothing may be left on the bus"
    );
    println!("  nothing left unclaimed on any face of the boundary");

    let Some(path) = common::output_path() else {
        println!("\npass a path to write an SVG, e.g. `cargo run --example beam_hot_spot out.svg`");
        return;
    };
    common::write(&path, &draw(&frames, absorbed));
}

/// Two panels: the profile at five instants, and the same data as a heat strip.
fn draw(frames: &[Frame], absorbed: f64) -> String {
    let hottest = frames
        .iter()
        .flat_map(|(_, pts, _, _)| pts.iter().map(|(_, t)| *t))
        .fold(f64::MIN, f64::max);
    let coldest = 20.0;
    let top = hottest + (hottest - coldest) * 0.12;

    let mut p = Plot::new(760.0, 460.0, (0.0, SPAN * 1e3), (coldest, top));
    p.axes(
        &ticks(0.0, SPAN * 1e3, 8),
        &ticks(coldest, top, 6),
        |v| format!("{v:.0}"),
        |v| format!("{v:.2}"),
    );
    p.title("Where a 1 mm beam puts its heat in a 20 mm aluminium bar");
    p.caption(&format!("{:.0} mW absorbed", absorbed * 1e3));

    // The heat strip along the bottom, from the last frame — the horizontal extent of the
    // colour is what "the beam is here" looks like without reading the axis.
    let (_, last, _, _) = frames.last().unwrap();
    let (lo, hi) = last.iter().fold((f64::MAX, f64::MIN), |(a, b), (_, t)| {
        (a.min(*t), b.max(*t))
    });
    let strip_top = coldest + (top - coldest) * 0.06;
    p.raster(
        last.len() - 1,
        1,
        (0.0, SPAN * 1e3),
        (coldest, strip_top),
        |i, _| heat(((last[i].1 - lo) / (hi - lo).max(1e-12)).clamp(0.0, 1.0)),
    );

    let colours = [
        rgb(160, 196, 232),
        rgb(96, 152, 214),
        rgb(48, 108, 186),
        rgb(198, 96, 48),
        rgb(158, 40, 32),
    ];
    for (i, (t, pts, mean, _)) in frames.iter().enumerate() {
        let colour = &colours[i.min(colours.len() - 1)];
        p.polyline(pts.iter().copied(), colour, 1.9);
        // The lumped model's whole answer, as a flat line, so the gap is the picture.
        p.polyline([(0.0, 20.0 + mean), (SPAN * 1e3, 20.0 + mean)], colour, 0.6);
        let label = if *t < 1.0 {
            format!("{:.0} ms", t * 1e3)
        } else {
            format!("{t:.0} s")
        };
        let peak = pts[pts.len() / 2].1;
        p.text(SPAN * 1e3 * 0.52, peak, &label, 12.0, colour, "start");
    }
    p.label(
        64.0,
        446.0,
        "thick: temperature along the bar   thin: the lumped mean, which every curve has exactly",
        11.5,
        "#6a6a6a",
        "start",
    );
    p.finish()
}
