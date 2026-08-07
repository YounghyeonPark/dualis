//! The radial distribution function: how much more likely a neighbour is at `r` than chance.
//!
//! ```text
//! g(r) = ⟨pairs in the shell at r⟩ / ⟨pairs an ideal gas would have there⟩
//! ```
//!
//! One is the ideal gas. Above one means the structure prefers that separation, below one means
//! it avoids it, and zero inside about `0.8σ` means atoms cannot overlap.
//!
//! One only in the thermodynamic limit, strictly. The denominator counts `N²/2` pairs where a
//! finite sample has `N(N−1)/2`, so a structureless system plateaus at `(N−1)/N` — 0.93% low at
//! 108 atoms and 3.1% low at 32. That is the standard convention and almost every code uses it,
//! but it is worth knowing before reading a small box's plateau as a discrepancy.
//! [`RadialDistribution::coordination`] is free of the offset: it divides by `N` and doubles the
//! pair count, so it returns neighbours per atom exactly. It is the single most
//! useful thing to look at in a molecular simulation, because it says what *phase* you have
//! without anyone having to decide: a crystal gives sharp spikes at the lattice shells, a liquid
//! gives two or three broad humps that die out, and a gas gives a bump at contact and then one
//! flat line forever.
//!
//! It is also measurable. `g(r)` is what a neutron or X-ray diffraction experiment returns, by
//! Fourier transform of the structure factor, so this is one of the few quantities where a
//! simulation and a real instrument produce the same curve.
//!
//! # Checked on a lattice first
//!
//! A face-centred cubic lattice has an exactly known answer, which is what makes it the right
//! thing to validate against before pointing this at a fluid. Its neighbour shells sit at
//!
//! ```text
//! a/√2,  a,  a√(3/2),  a√2,  ...
//! ```
//!
//! — ratios of `1 : √2 : √3 : 2` — holding `12, 6, 24, 12` neighbours. Those are combinatorial
//! facts about the lattice rather than measurements, so an implementation either reproduces them
//! or is wrong.
//!
//! # The range is half the box, not the potential's cutoff
//!
//! Two different limits that are easy to confuse. The potential stops at `rc` because further
//! pairs contribute nothing worth computing; `g(r)` stops at `L/2` because further than that the
//! minimum image is ambiguous and the sample is no longer a sphere. So the interesting structure
//! past `rc` — the second and third shells of a liquid — is visible here even though no force
//! acts across it.

use crate::box_::CellList;
use crate::fluid::Fluid;

/// Pair separations accumulated into shells, and the ideal-gas count to divide by.
///
/// Accumulated over snapshots rather than computed from one, because `g(r)` is an average and a
/// single configuration of a few hundred atoms is far too noisy to read a second peak off.
pub struct RadialDistribution {
    counts: Vec<f64>,
    bin_width: f64,
    range: f64,
    /// Snapshots added, so the counts can be turned into an average.
    snapshots: f64,
    particles: usize,
    number_density: f64,
}

impl RadialDistribution {
    /// Shells out to `range`, or to half the box if that is nearer.
    ///
    /// Clamped rather than refused: asking for the structure of a liquid past `L/2` is a
    /// reasonable thing to want and an impossible thing to have, and truncating is what every
    /// molecular-dynamics code does. [`RadialDistribution::range`] reports what was actually
    /// used, so a caller plotting the result knows where it stops.
    pub fn new(fluid: &Fluid, range: f64, bins: usize) -> RadialDistribution {
        let range = range
            .min(fluid.bounds().length / 2.0)
            .max(f64::MIN_POSITIVE);
        let bins = bins.max(1);
        RadialDistribution {
            counts: vec![0.0; bins],
            bin_width: range / bins as f64,
            range,
            snapshots: 0.0,
            particles: fluid.count(),
            number_density: fluid.number_density(),
        }
    }

    /// How far out the histogram reaches, after clamping to half the box.
    pub fn range(&self) -> f64 {
        self.range
    }

    /// How many shells the range is cut into.
    pub fn bins(&self) -> usize {
        self.counts.len()
    }

    /// Centre of a shell, which is where its `g` should be plotted.
    pub fn radius(&self, bin: usize) -> f64 {
        (bin as f64 + 0.5) * self.bin_width
    }

    /// Add one configuration.
    pub fn accumulate(&mut self, fluid: &Fluid) {
        let bounds = fluid.bounds();
        let positions: Vec<glam::DVec3> = (0..fluid.count()).map(|i| fluid.position(i)).collect();
        let list = CellList::build(bounds, self.range, &positions);
        let bin_width = self.bin_width;
        let counts = &mut self.counts;
        let bins = counts.len();
        list.for_each_pair(bounds, self.range, &positions, |_, _, _, r2| {
            let bin = (r2.sqrt() / bin_width) as usize;
            if bin < bins {
                counts[bin] += 1.0;
            }
        });
        self.snapshots += 1.0;
    }

    /// `g(r)` for one shell.
    ///
    /// The denominator is the number of pairs an ideal gas of the same density would put in the
    /// same shell: `(N/2)·ρ·V_shell`, with the shell's exact volume rather than `4πr²Δr`. The
    /// `N/2` is there because each pair was counted once, and using `N` is the standard way to
    /// end up with a `g(r)` that plateaus at two.
    pub fn g(&self, bin: usize) -> f64 {
        if bin >= self.counts.len() || self.snapshots <= 0.0 {
            return 0.0;
        }
        let inner = bin as f64 * self.bin_width;
        let outer = inner + self.bin_width;
        let shell = 4.0 / 3.0 * std::f64::consts::PI * (outer.powi(3) - inner.powi(3));
        let ideal = 0.5 * self.particles as f64 * self.number_density * shell;
        if ideal <= 0.0 {
            return 0.0;
        }
        self.counts[bin] / self.snapshots / ideal
    }

    /// `(r, g(r))` for every shell, ready to plot.
    pub fn curve(&self) -> Vec<(f64, f64)> {
        (0..self.bins())
            .map(|b| (self.radius(b), self.g(b)))
            .collect()
    }

    /// Average number of neighbours within `radius` of a particle.
    ///
    /// ```text
    /// n(r) = ∫₀ʳ 4πs²ρ g(s) ds
    /// ```
    ///
    /// The integral of the structure rather than the structure itself, and the more robust thing
    /// to compare against a lattice: a shell's *position* moves with the bin width but the count
    /// inside it does not.
    pub fn coordination(&self, radius: f64) -> f64 {
        let mut total = 0.0;
        for bin in 0..self.bins() {
            let inner = bin as f64 * self.bin_width;
            if inner >= radius {
                break;
            }
            // Two pairs per neighbour relationship, counted once, so twice the pairs over N.
            total += 2.0 * self.counts[bin] / self.snapshots / self.particles as f64;
        }
        total
    }

    /// Where the tallest peak is, which for a liquid is the first neighbour shell.
    pub fn first_peak(&self) -> (f64, f64) {
        (0..self.bins())
            .map(|b| (self.radius(b), self.g(b)))
            .fold((0.0, 0.0), |best, c| if c.1 > best.1 { c } else { best })
    }
}

/// The neighbour shells of a face-centred cubic lattice with cell edge `a`.
///
/// Radius and count, in order. Exact combinatorics rather than a measurement, which is what
/// makes them worth testing a histogram against.
pub fn fcc_shells(cell_edge: f64) -> [(f64, usize); 4] {
    let a = cell_edge;
    [
        (a / 2f64.sqrt(), 12),
        (a, 6),
        (a * 1.5f64.sqrt(), 24),
        (a * 2f64.sqrt(), 12),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid::{temperature_from_reduced, unit_mass, Thermostat};
    use crate::potential::LennardJones;
    use dualis_core::{Domain, Exchange};
    use dualis_units::{Temperature, Time};

    fn kelvin(reduced: f64) -> Temperature {
        temperature_from_reduced(reduced, &LennardJones::reduced())
    }

    /// **The exact case.** An fcc lattice's shells are at `1 : √2 : √3 : 2` times the nearest
    /// neighbour, holding 12, 6, 24 and 12 atoms, and the histogram has to find all four.
    ///
    /// Nothing statistical about it — one snapshot of a perfect lattice, and the answer is
    /// combinatorics. If this fails, no amount of averaging over a fluid would have revealed it.
    #[test]
    fn a_lattice_shows_its_exact_neighbour_shells() {
        let density = 1.2;
        let fluid = Fluid::lattice("ice", LennardJones::reduced(), unit_mass(), 4, density);
        // Four atoms to a cell fixes the edge: a = (4/rho)^(1/3).
        let a = (4.0 / density).cbrt();
        let shells = fcc_shells(a);

        let mut rdf = RadialDistribution::new(&fluid, shells[3].0 * 1.05, 600);
        rdf.accumulate(&fluid);

        // Each shell's atoms are all at exactly one radius, so its whole count lands in one bin
        // and the coordination number steps up by exactly that count. Probed *between* shells
        // rather than a hair either side of one: a histogram cannot resolve finer than a bin,
        // and asking it to is a test of the bin width instead of the physics.
        let between = |a: f64, b: f64| (a + b) / 2.0;
        let probes = [
            (shells[0].0 * 0.9, 0.0),
            (between(shells[0].0, shells[1].0), 12.0),
            (between(shells[1].0, shells[2].0), 18.0),
            (between(shells[2].0, shells[3].0), 42.0),
            (shells[3].0 * 1.02, 54.0),
        ];
        for (radius, expected) in probes {
            let found = rdf.coordination(radius);
            assert!(
                (found - expected).abs() < 1e-9,
                "within {radius}: {found} neighbours, expected {expected}"
            );
        }

        // The tallest peak is the twelve nearest neighbours.
        let (peak_r, peak_g) = rdf.first_peak();
        assert!(
            (peak_r / shells[0].0 - 1.0).abs() < 0.01,
            "peak at {peak_r}, first shell at {}",
            shells[0].0
        );
        assert!(
            peak_g > 20.0,
            "a lattice spike should be tall, got {peak_g}"
        );

        // The shell ratios, stated as the ratios rather than as four radii.
        assert!((shells[1].0 / shells[0].0 - 2f64.sqrt()).abs() < 1e-12);
        assert!((shells[2].0 / shells[0].0 - 3f64.sqrt()).abs() < 1e-12);
        assert!((shells[3].0 / shells[0].0 - 2.0).abs() < 1e-12);
    }

    /// **A dilute gas is structureless**: `g(r) = 1` beyond the core, because a neighbour at any
    /// distance is exactly as likely as chance.
    ///
    /// That is the normalisation being right and nothing else, which is why it is worth its own
    /// test — a `g` that plateaued at two would look plausible on a liquid and is the standard
    /// consequence of forgetting the factor of two in the pair count.
    #[test]
    fn a_dilute_gas_has_no_structure_to_speak_of() {
        let target = kelvin(3.0);
        let mut fluid = Fluid::lattice("gas", LennardJones::reduced(), unit_mass(), 3, 0.05)
            .thermalised(target, 0x_9A50_0001)
            .with_thermostat(Thermostat::Langevin {
                target,
                damping: 1.0,
            });
        let dt = fluid.max_stable_dt(Time::ZERO);
        let mut bus = Exchange::new();
        for _ in 0..1000 {
            fluid.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        let mut rdf = RadialDistribution::new(&fluid, 6.0, 60);
        for k in 0..2000 {
            fluid.step(Time::ZERO, dt, &mut bus).unwrap();
            if k % 5 == 0 {
                rdf.accumulate(&fluid);
            }
        }

        // Past two sigma there is no force and, at this density, no correlation either.
        let far: Vec<f64> = rdf
            .curve()
            .into_iter()
            .filter(|(r, _)| *r > 2.5)
            .map(|(_, g)| g)
            .collect();
        assert!(far.len() > 10, "the test needs bins out there");
        let mean = far.iter().sum::<f64>() / far.len() as f64;
        assert!(
            (mean - 1.0).abs() < 0.05,
            "a structureless gas should sit at one, got {mean}"
        );

        // And nothing gets inside the repulsive wall.
        for (r, g) in rdf.curve() {
            if r < 0.7 {
                assert_eq!(g, 0.0, "an atom at {r} sigma should be impossible");
            }
        }
    }

    /// A liquid is between the two: a tall first peak, a second hump, and then flat.
    ///
    /// Not an exact claim, which is why it is asserted qualitatively — but the *shape* is what
    /// distinguishes a liquid from a gas and from a crystal, and the first peak of a
    /// Lennard-Jones liquid sitting near the potential minimum is a fact worth pinning.
    #[test]
    fn a_liquid_peaks_near_the_potential_minimum() {
        let target = kelvin(0.85);
        let mut fluid = Fluid::lattice("liquid", LennardJones::reduced(), unit_mass(), 4, 0.8442)
            .thermalised(target, 0x_119A_1D00)
            .with_thermostat(Thermostat::Langevin {
                target,
                damping: 1.0,
            });
        let dt = fluid.max_stable_dt(Time::ZERO);
        let mut bus = Exchange::new();
        for _ in 0..2000 {
            fluid.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        let mut rdf = RadialDistribution::new(&fluid, 5.0, 100);
        for k in 0..3000 {
            fluid.step(Time::ZERO, dt, &mut bus).unwrap();
            if k % 5 == 0 {
                rdf.accumulate(&fluid);
            }
        }

        let (peak_r, peak_g) = rdf.first_peak();
        // Near the well's minimum at 1.1225 sigma, which is where an atom most wants to sit.
        assert!(
            (peak_r - LennardJones::reduced().minimum()).abs() < 0.15,
            "first peak at {peak_r}, the well bottoms at {}",
            LennardJones::reduced().minimum()
        );
        // A liquid's first peak is a few times the ideal density -- not a lattice spike, and
        // not the flat line a gas gives.
        assert!(
            (2.0..6.0).contains(&peak_g),
            "a liquid first peak should be a few times one, got {peak_g}"
        );
        // Structure dies out. Measured over the outer third of whatever range the box allowed
        // rather than past a fixed radius: at this density 256 atoms give a box of 6.7 sigma, so
        // `g` stops at 3.36 and a filter for `r > 4` selects nothing and averages to NaN. Which
        // it did, on the first attempt.
        let curve = rdf.curve();
        let tail: Vec<f64> = curve
            .iter()
            .filter(|(r, _)| *r > rdf.range() * 2.0 / 3.0)
            .map(|(_, g)| *g)
            .collect();
        assert!(
            tail.len() > 5,
            "the tail needs bins in it, got {}",
            tail.len()
        );
        let mean = tail.iter().sum::<f64>() / tail.len() as f64;
        assert!((mean - 1.0).abs() < 0.1, "the tail sat at {mean}");

        // And the first shell holds about twelve neighbours, which is close packing -- a liquid
        // keeps its neighbours and loses only the order in which it has them.
        let first_shell = rdf.coordination(1.55);
        assert!(
            (8.0..14.0).contains(&first_shell),
            "first-shell coordination was {first_shell}"
        );
    }

    /// The range is clamped to half the box rather than being taken on trust.
    #[test]
    fn the_range_stops_at_half_the_box() {
        let fluid = Fluid::lattice("small", LennardJones::reduced(), unit_mass(), 3, 0.8442);
        let half = fluid.bounds().length / 2.0;
        let rdf = RadialDistribution::new(&fluid, 1000.0, 50);
        assert!(
            (rdf.range() - half).abs() < 1e-12,
            "range was {}",
            rdf.range()
        );
        // A request inside the box is honoured as asked.
        let modest = RadialDistribution::new(&fluid, 2.0, 50);
        assert!((modest.range() - 2.0).abs() < 1e-12);
    }
}
