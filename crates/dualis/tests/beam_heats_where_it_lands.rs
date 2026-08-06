//! The claim a lumped coupling could not make: heat arrives *where* the light did.
//!
//! `multiphysics.rs` proved energy crosses between two domains that do not know about each
//! other. It could not prove anything about *place*, because the bus carried one number per
//! channel per step — so "a 1 mm beam on a 20 mm mirror" and "1 W spread over the whole
//! mirror" were the same message, and the second one is wrong in the way that matters. A
//! coating fails at its hot spot, not at its average.
//!
//! With `Flux` and `Interface` the message has a place in it, and this file is what that
//! buys:
//!
//! - The **mean** temperature is still exactly the lumped answer, `Q/C`. The spatial model
//!   contains the lumped one; it does not replace it.
//! - The **peak** is several times that mean, and a designer sizing a coating from `Q/C`
//!   would have under-predicted the thing that breaks.
//! - The answer does not depend on the grid, which is the property that makes a spatial
//!   coupling a model rather than an artefact of its mesh.
//! - The settled shape matches a solution computed by quadrature from the steady-state
//!   energy balance, which never went through the coupling at all.
//! - And the coupling *interval* turns out to be the largest error source here. It is also
//!   the one no conservation check can see: a 10 ms window delivers exactly the right joules
//!   and reads the peak 12% low.
//!
//! # One expectation this file had to correct
//!
//! An insulated bar under a continuous beam does not flatten out. It settles into a fixed
//! shape and then rides upward on a mean that climbs forever — so the hot spot is permanent,
//! and the lumped model's error in kelvin never shrinks. What shrinks is the *fraction*,
//! because the mean grows without bound. The intuition that says "given long enough it
//! evens out" is about a bar heated *once*, and it does not survive leaving the source on.

use dualis::prelude::*;
use dualis_optics::spectrum::Spectrum as Spec;

/// The absorbed power of a beam, published with the profile it landed with.
///
/// Everything except the last line is what `AbsorbingSurface` in `multiphysics.rs` already
/// did: a real spectral absorptance integrated against a real lamp, giving watts. The
/// difference is that those watts are handed over as a distribution over the boundary the
/// consumer owns, so the consumer does not have to guess.
struct ProfiledBeam {
    lamp: SpectralPower,
    absorptance: Spectrum,
    /// The consumer's discretisation, taken as given. This domain publishes onto the grid
    /// the boundary has rather than onto one of its own, which is why nothing interpolates.
    boundary: Interface,
    /// `w`, the 1/e² irradiance radius.
    waist: Length,
    /// How long the illuminated boundary is, so a waist can be compared to it.
    span: Length,
    paid_out: f64,
}

impl ProfiledBeam {
    fn new(
        lamp: SpectralPower,
        optics: &SurfaceOptics,
        boundary: &Interface,
        waist: Length,
        span: Length,
    ) -> ProfiledBeam {
        let absorptance = Spec::curve(
            (0..=60)
                .map(|i| {
                    let nm = 400.0 + i as f64 * 5.0;
                    (nm, optics.absorptance(Length::nm(nm)))
                })
                .collect(),
        );
        ProfiledBeam {
            lamp,
            absorptance,
            boundary: boundary.clone(),
            waist,
            span,
            paid_out: 0.0,
        }
    }

    fn absorbed_power(&self) -> Power {
        self.lamp.absorbed_by(&self.absorptance)
    }
}

impl Domain for ProfiledBeam {
    fn name(&self) -> &'static str {
        "beam"
    }

    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let joules = self.absorbed_power().to_si() * dt.to_si();
        // A centred Gaussian, written as the physics writes it: `exp(-2r²/w²)` in the
        // boundary's own coordinate. The kernel normalises it to `joules` exactly.
        let (span, waist) = (self.span.to_si(), self.waist.to_si());
        let flux = Flux::profiled(joules, &self.boundary, |u| {
            (-2.0 * ((u - 0.5) * span / waist).powi(2)).exp()
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

/// 20 mm of aluminium, side-on to the beam.
const SPAN: f64 = 20e-3;
/// A square 10 mm bar, so a cell's cross-section is 1 cm².
const CROSS_SECTION: f64 = 1e-4;

fn bar(cells: usize) -> Bar1D {
    let dx = Length::from_si(SPAN / cells as f64);
    Bar1D::new(
        "bar",
        Substance::aluminium_6061(),
        cells,
        dx,
        Area::from_si(CROSS_SECTION),
        Temperature::celsius(20.0),
    )
    // The illuminated side of one cell: 10 mm across by `dx` along.
    .exposing("bar face", Area::from_si(10e-3 * dx.to_si()))
}

/// A 100 W green laser on a mirror that absorbs about a percent of it.
fn beam_for(boundary: &Interface, waist_mm: f64) -> ProfiledBeam {
    let laser = SpectralPower::new(
        Spec::gaussian(Length::nm(532.0), Length::nm(1.0), 1.0),
        Power::w(100.0),
        (Length::nm(400.0), Length::nm(700.0)),
    );
    let mirror = SurfaceOptics::mirror(0.99);
    ProfiledBeam::new(
        laser,
        &mirror,
        boundary,
        Length::mm(waist_mm),
        Length::from_si(SPAN),
    )
}

/// Total heat capacity of the bar, which is the denominator of the lumped answer and does
/// not depend on how many cells it was cut into.
fn bar_capacity() -> f64 {
    let al = Substance::aluminium_6061();
    let volume = Volume::from_si(CROSS_SECTION * SPAN);
    al.heat_capacity(volume).unwrap().to_si()
}

/// **The test this whole layer exists for.** A 1 mm beam on a 20 mm bar: the mean rise is
/// exactly the lumped answer, and the peak is many times it.
#[test]
fn a_focused_beam_leaves_a_hot_spot_a_lumped_model_cannot_see() {
    let cells = 81;
    let plate = bar(cells);
    let boundary = plate.boundary().unwrap().clone();
    let beam = beam_for(&boundary, 1.0);

    let absorbed = beam.absorbed_power();
    assert!(
        absorbed.in_mw() > 500.0 && absorbed.in_mw() < 1500.0,
        "a 1% mirror under 100 W absorbs about a watt, got {} mW",
        absorbed.in_mw()
    );

    let mut sim = Simulation::new(Schedule::Multirate).with(beam).with(plate);

    // 50 ms. Long enough for the deposited heat to start spreading — aluminium's diffusion
    // length √(αt) is 1.9 mm here, twice the waist — and short enough that it has not
    // reached the ends of a 20 mm bar. Which is the regime a real mirror lives in under a
    // real beam, and the regime a lumped model has nothing to say about.
    for _ in 0..50 {
        sim.advance(Time::ms(1.0))
            .expect("energy must be conserved");
    }

    // The audit inside `advance` already proved conservation face by face. State it.
    let residual = sim.ledger().get(quantity::ENERGY).unwrap();
    assert!(residual.abs() < 1e-9, "energy residual {residual}");
    assert!(sim.bus().unclaimed().next().is_none());

    let plate: &Bar1D = sim.domain_as("bar").unwrap();
    let deposited = absorbed.to_si() * 0.050;
    assert!(
        (plate.absorbed_energy().to_si() / deposited - 1.0).abs() < 1e-12,
        "{} J arrived of {deposited} J",
        plate.absorbed_energy().to_si()
    );

    // The closed form, and it is exact: insulated ends, so the mean rise is Q/C whatever
    // the heat did on the way. This is the lumped model's entire answer.
    let lumped_rise = deposited / bar_capacity();
    let mean_rise = plate.mean_temperature().in_celsius() - 20.0;
    assert!(
        (mean_rise / lumped_rise - 1.0).abs() < 1e-9,
        "the mean must be the lumped answer exactly: {mean_rise} against {lumped_rise}"
    );

    // And the part the lumped model has no access to. The peak is in the middle, where the
    // beam was, and it is several times the mean.
    let middle = cells / 2;
    let peak_rise = plate.temperature_at(middle).in_celsius() - 20.0;
    let hottest = (0..cells)
        .max_by(|a, b| {
            plate
                .temperature_at(*a)
                .to_si()
                .total_cmp(&plate.temperature_at(*b).to_si())
        })
        .unwrap();
    assert_eq!(hottest, middle, "the hot spot is where the beam is");
    assert!(
        peak_rise / mean_rise > 3.0,
        "a 1 mm beam on a 20 mm bar must concentrate: peak {peak_rise} K against mean {mean_rise} K"
    );

    // The ends are barely touched — a lumped model would have warmed them by the full mean.
    let end_rise = plate.temperature_at(0).in_celsius() - 20.0;
    assert!(
        end_rise < mean_rise * 0.2,
        "the far end should still be cold: {end_rise} K against a mean of {mean_rise} K"
    );

    // Symmetric about the beam, since the illumination and the conduction both are.
    for offset in 1..middle {
        let left = plate.temperature_at(middle - offset).in_celsius() - 20.0;
        let right = plate.temperature_at(middle + offset).in_celsius() - 20.0;
        assert!(
            (left - right).abs() < 1e-9 * peak_rise.max(1e-30),
            "asymmetric at offset {offset}: {left} against {right}"
        );
    }
}

/// **The property that makes it a model and not a mesh artefact.** Two bars cut into
/// different numbers of cells, under the same physical beam, reach the same peak.
///
/// This is the check the spatial coupling has to pass before any of its numbers mean
/// anything, and it is only possible because the publisher takes the consumer's
/// discretisation as given rather than imposing one.
#[test]
fn the_answer_does_not_depend_on_the_grid() {
    let run = |cells: usize| -> (f64, f64) {
        let plate = bar(cells);
        let boundary = plate.boundary().unwrap().clone();
        let mut sim = Simulation::new(Schedule::Multirate)
            .with(beam_for(&boundary, 2.0))
            .with(plate);
        for _ in 0..50 {
            sim.advance(Time::ms(1.0))
                .expect("energy must be conserved");
        }
        let plate: &Bar1D = sim.domain_as("bar").unwrap();
        (
            plate.temperature_at(cells / 2).in_celsius() - 20.0,
            plate.mean_temperature().in_celsius() - 20.0,
        )
    };

    // 41, 81 and 161 cells: a 4× range in grid spacing and a 16× range in stable step.
    let (coarse_peak, coarse_mean) = run(41);
    let (medium_peak, medium_mean) = run(81);
    let (fine_peak, fine_mean) = run(161);

    // The mean is exact on every grid, because it is Q/C and C does not care about cells.
    for (label, mean) in [("41", coarse_mean), ("81", medium_mean), ("161", fine_mean)] {
        assert!(
            (mean / coarse_mean - 1.0).abs() < 1e-9,
            "{label} cells gave a different mean: {mean} against {coarse_mean}"
        );
    }

    // The peak converges. It is a sampled value of a continuous field, so it is not
    // identical across grids — a coarse cell averages over more of the profile and reads
    // low. Refining must reduce the disagreement, and it does: doubling the grid roughly
    // halves it, which is the first-order convergence a midpoint-sampled profile has.
    let coarse_gap = (coarse_peak / fine_peak - 1.0).abs();
    let medium_gap = (medium_peak / fine_peak - 1.0).abs();
    assert!(
        coarse_gap < 0.10,
        "41 cells should already be within 10%: {coarse_gap}"
    );
    assert!(
        medium_gap < coarse_gap * 0.75,
        "refining must converge: {coarse_gap} at 41 cells, {medium_gap} at 81"
    );
}

/// The quasi-steady temperature profile an insulated bar reaches under continuous heating,
/// solved by quadrature rather than by marching — the independent answer this file's
/// simulations are checked against.
///
/// # Where it comes from
///
/// Once the transient has died, every part of the bar rises at the same rate `ṫ = Q̇/C`,
/// because the ends are insulated and there is nowhere else for the heat to go. Balance the
/// energy of the segment from 0 to `x`: the source deposits `S(x)`, conduction carries
/// `−kA·T'(x)` out through `x`, and what is left raises `ρc_pAx` of material at `ṫ`. So
///
/// ```text
/// T'(x) = (ṫ/α) · [x − L·F(x)],     F(x) = S(x)/Q̇
/// ```
///
/// with `F` the fraction of the source falling left of `x`. Integrating once gives the
/// shape and integrating again gives its mean, both by trapezoid here.
///
/// For a point source at the centre this has the closed form `peak − mean = ṫL²/12α`, which
/// is asserted below as the limit of this quadrature — so the quadrature is itself checked
/// against an exact result before anything else is checked against it.
///
/// Note what makes this independent rather than circular: it solves the *steady* problem as
/// an integral where the simulation marches the *transient* problem as a difference
/// equation. Neither shares a discretisation or a line of code with the other.
fn quasi_steady_peak_above_mean(
    rise_rate: f64,
    alpha: f64,
    length: f64,
    profile: impl Fn(f64) -> f64,
) -> f64 {
    // Odd, so one sample sits exactly at the centre where the peak is.
    const N: usize = 40_001;
    let h = 1.0 / (N - 1) as f64;

    let source: Vec<f64> = (0..N).map(|i| profile(i as f64 * h)).collect();
    let mut fraction = vec![0.0; N];
    for i in 1..N {
        fraction[i] = fraction[i - 1] + 0.5 * (source[i - 1] + source[i]) * h;
    }
    let total = fraction[N - 1];

    // In u = x/L the slope is (ṫL²/α)·[u − F(u)].
    let k = rise_rate * length * length / alpha;
    let slope: Vec<f64> = (0..N)
        .map(|i| k * (i as f64 * h - fraction[i] / total))
        .collect();
    let mut shape = vec![0.0; N];
    for i in 1..N {
        shape[i] = shape[i - 1] + 0.5 * (slope[i - 1] + slope[i]) * h;
    }

    let peak = shape[(N - 1) / 2];
    let mean = h * (shape.iter().sum::<f64>() - 0.5 * (shape[0] + shape[N - 1]));
    peak - mean
}

/// The quadrature above, against the one case with an exact answer: a point source at the
/// centre gives `peak − mean = ṫL²/12α`.
///
/// First, because everything else in this file leans on it.
#[test]
fn the_reference_solution_reproduces_its_own_closed_form() {
    let (rate, alpha, length) = (0.2, 6.9e-5, 20e-3);
    let exact = rate * length * length / (12.0 * alpha);

    // How far a beam of a given width falls short of the point-source peak.
    let deficit = |width: f64| {
        let value = quasi_steady_peak_above_mean(rate, alpha, length, |u| {
            (-2.0 * ((u - 0.5) / width).powi(2)).exp()
        });
        1.0 - value / exact
    };

    // A beam a thousandth of the bar wide is nearly a point source, but not quite: rounding
    // off the kink in the parabola over a width `w` costs an amount linear in `w`, and here
    // that is a quarter of a percent.
    let narrow = deficit(1e-3);
    assert!(
        narrow > 0.0 && narrow < 4e-3,
        "a near-point source should be a fraction of a percent short, got {narrow}"
    );

    // Asserting the *scaling* rather than picking a tolerance, which is what makes this a
    // check on the quadrature instead of a check on my patience. Doubling the beam doubles
    // the shortfall; ten times the beam is ten times the shortfall.
    assert!(
        (deficit(2e-3) / narrow - 2.0).abs() < 0.05,
        "the shortfall must be linear in the width: {} against {narrow}",
        deficit(2e-3)
    );
    assert!(
        (deficit(1e-2) / narrow - 10.0).abs() < 0.5,
        "and still linear a decade out: {} against {narrow}",
        deficit(1e-2)
    );

    // Monotonic all the way up, and a beam a *tenth* of the bar wide still peaks at three
    // quarters of the point-source value. The concentration is not fragile, which is why it
    // matters in practice rather than only in the limit.
    let mut previous = 0.0;
    for width in [0.01, 0.02, 0.05, 0.1] {
        let d = deficit(width);
        assert!(
            d > previous,
            "width {width} did not reduce the peak further"
        );
        previous = d;
    }
    assert!(previous < 0.3, "a tenth of the bar cost {previous}");
}

/// **Where the spatial coupling's numbers come from, checked against that reference.**
///
/// Under continuous heating an insulated bar does *not* flatten out. It settles into a
/// quasi-steady shape — hottest under the beam, coldest at the ends — and that shape then
/// rides upward on a mean that climbs forever. So the hot spot never goes away, and the
/// lumped model stays wrong by a fixed number of kelvin for as long as the beam is on.
///
/// What does vanish is the *fractional* error, because the mean grows without bound. That is
/// the honest statement of when a lumped model is good enough, and it is a statement about
/// the ratio rather than about the difference. Getting that backwards is easy: a bar that is
/// heated *once* really does flatten, and the intuition does not survive the source being
/// left on.
#[test]
fn the_profile_settles_to_a_shape_that_matches_the_reference_solution() {
    let cells = 81;
    let waist = 1e-3;
    let plate = bar(cells);
    let boundary = plate.boundary().unwrap().clone();
    let mut sim = Simulation::new(Schedule::Multirate)
        .with(beam_for(&boundary, waist * 1e3))
        .with(plate);

    let sample = |sim: &Simulation| -> (f64, f64) {
        let plate: &Bar1D = sim.domain_as("bar").unwrap();
        let mean = plate.mean_temperature().in_celsius() - 20.0;
        (
            mean,
            plate.temperature_at(cells / 2).in_celsius() - 20.0 - mean,
        )
    };

    // A coupling window comfortably under the bar's own 0.44 ms step, so the transfer
    // interval is not the limiting error — the test below is about what happens when it is.
    let window = Time::ms(0.5);
    // Eight seconds, against a bar diffusion time L²/α of 5.8 s. The excess was measured
    // constant to twelve figures from six seconds on, so this is well past settled.
    let steps = 16_000;
    for _ in 0..steps {
        sim.advance(window).expect("energy must be conserved");
    }
    let (first_mean, first_excess) = sample(&sim);

    let alpha = Substance::aluminium_6061().diffusivity().unwrap().to_si();
    let rate = first_mean / sim.time().to_si();
    let reference = quasi_steady_peak_above_mean(rate, alpha, SPAN, |u| {
        (-2.0 * ((u - 0.5) * SPAN / waist).powi(2)).exp()
    });
    assert!(
        (first_excess / reference - 1.0).abs() < 0.02,
        "the settled shape must match the reference: {first_excess} K against {reference} K"
    );

    // Settled, and it stays settled: another eight seconds doubles the mean and leaves the
    // excess exactly where it was. Two claims in one — the shape is stationary, and the hot spot is
    // permanent.
    for _ in 0..steps {
        sim.advance(window).unwrap();
    }
    let (second_mean, second_excess) = sample(&sim);
    assert!(
        (second_mean / first_mean - 2.0).abs() < 1e-6,
        "the mean must keep climbing linearly: {first_mean} then {second_mean}"
    );
    assert!(
        (second_excess / first_excess - 1.0).abs() < 1e-6,
        "the shape must be stationary: {first_excess} K then {second_excess} K"
    );

    // So the lumped model's absolute error never shrinks and its fractional error halves as
    // the mean doubles. Which is the whole answer to "when can I use a lumped model": when
    // you care about the average, never when you care about the hottest point.
    let first_fraction = first_excess / first_mean;
    let second_fraction = second_excess / second_mean;
    assert!(
        (second_fraction / first_fraction - 0.5).abs() < 1e-5,
        "the fraction must halve: {first_fraction} then {second_fraction}"
    );

    let residual = sim.ledger().get(quantity::ENERGY).unwrap();
    assert!(residual.abs() < 1e-9, "energy residual {residual}");
    assert!(sim.bus().unclaimed().next().is_none());
}

/// **How often the domains talk is itself an error source, and here it is the largest one.**
///
/// The bar's own stability limit is 0.44 ms on this grid, so a 10 ms coupling window means it
/// takes the whole window's heat on its first substep and then conducts for twenty-two
/// substeps with nothing arriving. The energy is exactly right — the audit says so — but the
/// profile spends most of each window relaxing, and the peak reads 12% low.
///
/// Worth a test of its own because it is invisible to every check a conservative coupling
/// normally gets. Nothing is lost, nothing is created, both domains are inside their
/// stability limits, and the answer is still wrong. The only thing that catches it is
/// comparing against a solution that never went through the coupling at all.
#[test]
fn the_coupling_interval_is_an_error_source_and_it_converges() {
    let cells = 81;
    let waist = 1e-3;
    let alpha = Substance::aluminium_6061().diffusivity().unwrap().to_si();

    let settled = |window_ms: f64| -> (f64, f64) {
        let plate = bar(cells);
        let boundary = plate.boundary().unwrap().clone();
        let mut sim = Simulation::new(Schedule::Multirate)
            .with(beam_for(&boundary, waist * 1e3))
            .with(plate);
        // Eight seconds in every case, so only the window differs.
        for _ in 0..(8.0 / (window_ms * 1e-3)).round() as usize {
            sim.advance(Time::ms(window_ms))
                .expect("energy is conserved whatever the window");
        }
        let plate: &Bar1D = sim.domain_as("bar").unwrap();
        let mean = plate.mean_temperature().in_celsius() - 20.0;
        let excess = plate.temperature_at(cells / 2).in_celsius() - 20.0 - mean;
        (mean, excess)
    };

    let (mean, coarse) = settled(10.0);
    let (mean_1ms, medium) = settled(1.0);
    let (mean_02ms, fine) = settled(0.2);

    // Every window deposits the same total, to the last bit the audit can see. Which is the
    // point of the test: the conserved quantity is right in all three and the distribution
    // is not, so conservation was never going to catch this.
    for (label, m) in [("1 ms", mean_1ms), ("0.2 ms", mean_02ms)] {
        assert!(
            (m / mean - 1.0).abs() < 1e-9,
            "{label} must deposit the same total: {m} against {mean}"
        );
    }

    let reference = quasi_steady_peak_above_mean(mean / 8.0, alpha, SPAN, |u| {
        (-2.0 * ((u - 0.5) * SPAN / waist).powi(2)).exp()
    });
    let error = |v: f64| (v / reference - 1.0).abs();

    assert!(
        error(coarse) > 0.08 && error(coarse) < 0.16,
        "a 10 ms window should read about 12% low, got {:.4}",
        error(coarse)
    );
    assert!(
        error(medium) < error(coarse) / 3.0,
        "1 ms must be much better than 10 ms: {:.4} against {:.4}",
        error(medium),
        error(coarse)
    );
    assert!(
        error(fine) < 0.005,
        "a 0.2 ms window should be inside half a percent, got {:.4}",
        error(fine)
    );
}
