//! A periodic box, and finding neighbours inside one in linear time.
//!
//! # Why periodic
//!
//! A thousand atoms in a vacuum is mostly surface: the cube root of 1000 is 10, so about half
//! of them are in the outer shell, and whatever they do dominates the average. Periodic
//! boundaries make the box tile space instead, so every atom has a full neighbourhood and none
//! of them is at an edge. It is not a trick — it is the statement that the sample stands for a
//! bulk, and it is what makes a thousand atoms say anything about a litre of argon.
//!
//! What it costs is honest to name. Nothing in the box can be larger than the box, so a
//! correlation longer than `L/2` is cut off, and a fluctuation whose wavelength exceeds `L` is
//! not represented at all. Near a critical point, where correlations diverge, that is fatal
//! and no amount of averaging fixes it.
//!
//! # The minimum image convention
//!
//! Each pair interacts through the nearest of the infinitely many copies. That is only
//! unambiguous while the cutoff is under `L/2`, because past that a particle would see two
//! images of the same neighbour — and worse, itself. [`PeriodicBox::admits`] refuses that
//! rather than silently double-counting.
//!
//! # Cell lists
//!
//! Sorting particles into cells at least a cutoff wide means a neighbour can only be in the
//! same cell or one of the 26 around it, so the work is proportional to the number of
//! particles rather than to its square. At 500 particles that is already a factor of ten; at
//! 10 000 it is the difference between a second and a minute.

use glam::DVec3;

/// A cubic box that tiles space.
///
/// Cubic rather than general: an orthorhombic or triclinic cell is a small generalisation of
/// the arithmetic and a large one of the testing, and nothing here needs a shape yet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeriodicBox {
    /// Edge length in metres. The box runs from 0 to this in each direction.
    pub length: f64,
}

impl PeriodicBox {
    /// A cube of the given edge length, in metres.
    pub fn cubic(length: f64) -> PeriodicBox {
        PeriodicBox {
            length: length.max(f64::MIN_POSITIVE),
        }
    }

    /// A box holding `count` particles at a given number density.
    pub fn for_density(count: usize, number_density: f64) -> PeriodicBox {
        let volume = count as f64 / number_density.max(f64::MIN_POSITIVE);
        PeriodicBox::cubic(volume.cbrt())
    }

    /// Volume in cubic metres.
    pub fn volume(&self) -> f64 {
        self.length.powi(3)
    }

    /// Whether the minimum image convention is unambiguous at this cutoff.
    ///
    /// It needs `rc < L/2`. At exactly `L/2` a particle sits equidistant from two images and
    /// the answer depends on a rounding mode, which is not a basis for physics.
    pub fn admits(&self, cutoff: f64) -> bool {
        cutoff < self.length / 2.0
    }

    /// Fold a position back into the box.
    ///
    /// `rem_euclid` rather than a subtract-in-a-loop: a particle that has been kicked several
    /// box lengths in one step is a sign of trouble elsewhere, but wrapping it should still
    /// terminate.
    pub fn wrap(&self, p: DVec3) -> DVec3 {
        DVec3::new(
            p.x.rem_euclid(self.length),
            p.y.rem_euclid(self.length),
            p.z.rem_euclid(self.length),
        )
    }

    /// The shortest vector from `b` to `a` over all images.
    ///
    /// Each component is brought into `[−L/2, L/2)`, which is the minimum image. The rounding
    /// is done with `round`, so a component at exactly `L/2` goes to `−L/2` — consistent, and
    /// the reason [`PeriodicBox::admits`] keeps the cutoff strictly inside.
    pub fn shortest(&self, a: DVec3, b: DVec3) -> DVec3 {
        let d = a - b;
        d - self.length
            * DVec3::new(
                (d.x / self.length).round(),
                (d.y / self.length).round(),
                (d.z / self.length).round(),
            )
    }
}

/// Particles sorted into cells, so that a neighbour search visits a constant number of them.
///
/// Rebuilt every step rather than maintained incrementally. Incremental updates are faster and
/// are how a production code does it, but they carry state that can disagree with the
/// positions, and a neighbour list that is quietly stale produces physics that looks right.
/// Rebuilding is `O(N)` and the pair loop is `O(N)` too, so the constant is what is being
/// traded, not the order.
pub struct CellList {
    /// Cells per side.
    divisions: usize,
    cell_size: f64,
    /// Particle indices, grouped by cell and contiguous.
    entries: Vec<u32>,
    /// Where each cell's run starts in `entries`, with a final sentinel.
    starts: Vec<u32>,
}

impl CellList {
    /// Sort `positions` into cells at least `cutoff` wide.
    ///
    /// Falls back to a single cell — and therefore to every pair — when the box cannot hold
    /// three divisions. Below three, a cell's 26 neighbours start repeating and the
    /// bookkeeping to avoid double-counting costs more than the loop it saves.
    pub fn build(bounds: PeriodicBox, cutoff: f64, positions: &[DVec3]) -> CellList {
        let divisions = if cutoff > 0.0 {
            ((bounds.length / cutoff).floor() as usize).max(1)
        } else {
            1
        };
        let divisions = if divisions < 3 { 1 } else { divisions };
        let cell_size = bounds.length / divisions as f64;
        let cells = divisions * divisions * divisions;

        // Counting sort: one pass to count, one to place. Stable in index order, so the
        // result is a function of the input and not of any hashing.
        let mut counts = vec![0u32; cells + 1];
        let index_of = |p: DVec3| -> usize {
            let coord = |v: f64| {
                let k = (v.rem_euclid(bounds.length) / cell_size) as usize;
                k.min(divisions - 1)
            };
            (coord(p.z) * divisions + coord(p.y)) * divisions + coord(p.x)
        };
        for p in positions {
            counts[index_of(*p) + 1] += 1;
        }
        for k in 0..cells {
            counts[k + 1] += counts[k];
        }
        let starts = counts.clone();
        let mut cursor = counts;
        let mut entries = vec![0u32; positions.len()];
        for (i, p) in positions.iter().enumerate() {
            let cell = index_of(*p);
            entries[cursor[cell] as usize] = i as u32;
            cursor[cell] += 1;
        }

        CellList {
            divisions,
            cell_size,
            entries,
            starts,
        }
    }

    /// Cells per side. One means the list degenerated to every pair, which happens in a box
    /// too small to divide — see [`CellList::build`].
    pub fn divisions(&self) -> usize {
        self.divisions
    }

    /// Edge length of one cell, which is never less than the cutoff it was built for.
    pub fn cell_size(&self) -> f64 {
        self.cell_size
    }

    /// Call `visit` once for every pair within the cutoff, each pair exactly once.
    ///
    /// The half-neighbour trick: of the 27 cells around one, only 13 and a half are visited —
    /// the cell itself with `i < j`, and thirteen of its neighbours chosen so that the
    /// opposite thirteen are reached from the other side. That halves the work and, more
    /// usefully, means a pair is seen once so Newton's third law can be applied by adding a
    /// force to one particle and subtracting it from the other. Momentum is then conserved bit
    /// for bit rather than nearly.
    ///
    /// `visit` receives the two indices and the minimum-image separation from `j` to `i`.
    pub fn for_each_pair<F>(
        &self,
        bounds: PeriodicBox,
        cutoff: f64,
        positions: &[DVec3],
        mut visit: F,
    ) where
        F: FnMut(usize, usize, DVec3, f64),
    {
        let rc2 = cutoff * cutoff;
        let length = bounds.length;
        let consider = |i: usize, j: usize, visit: &mut F| {
            // The minimum image, with the wrap skipped when it provably does not apply.
            //
            // `PeriodicBox::shortest` costs three divisions and three `round()` calls, and it
            // runs on every *candidate* — roughly four hundred per atom, of which about fifty
            // are inside the cutoff. The cost is in the ones that fail, which measurement put at
            // ~10 ns each and made the dominant term in `Fluid::step`.
            //
            // When every component satisfies `2|dₖ| < L` the nearest image is the direct
            // difference. `2 * x` is exact in binary floating point, so that comparison is the
            // exact `|dₖ| < L/2` and not an approximation of it.
            //
            // No pair is *missed*: a component at or beyond `L/2` takes the slow path, so the
            // wrap still happens wherever it might matter.
            //
            // And where the fast path is taken it returns the same bits, with one degenerate
            // exception worth naming rather than glossing. `round(dₖ/L)` is zero whenever
            // `|dₖ| < L/2`, *unless* the quotient rounds to exactly `0.5` — which needs `|dₖ|`
            // within an ulp of half a box. There the two answers are `d` and `d − L`, the two
            // periodic images of a pair exactly half a box apart, and they are equidistant: the
            // minimum image is genuinely ambiguous and either is as correct as the other.
            //
            // It is also deterministic. The same input takes the same branch on every platform,
            // which is the property this workspace actually promises — not that `shortest` and
            // this agree on a tie, but that a run gives the same answer everywhere.
            let direct = positions[i] - positions[j];
            let d = if 2.0 * direct.x.abs() < length
                && 2.0 * direct.y.abs() < length
                && 2.0 * direct.z.abs() < length
            {
                direct
            } else {
                bounds.shortest(positions[i], positions[j])
            };
            let r2 = d.length_squared();
            if r2 < rc2 {
                visit(i, j, d, r2);
            }
        };

        if self.divisions == 1 {
            for i in 0..positions.len() {
                for j in (i + 1)..positions.len() {
                    consider(i, j, &mut visit);
                }
            }
            return;
        }

        let n = self.divisions as isize;
        let cell_at = |x: isize, y: isize, z: isize| -> usize {
            let w = |v: isize| v.rem_euclid(n) as usize;
            (w(z) * self.divisions + w(y)) * self.divisions + w(x)
        };
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let here = cell_at(x, y, z);
                    let mine =
                        &self.entries[self.starts[here] as usize..self.starts[here + 1] as usize];
                    // Inside the cell, ordered so each pair appears once.
                    for (a, &i) in mine.iter().enumerate() {
                        for &j in &mine[a + 1..] {
                            consider(i as usize, j as usize, &mut visit);
                        }
                    }
                    // And the thirteen forward neighbours.
                    for (dx, dy, dz) in HALF_NEIGHBOURS {
                        let there = cell_at(x + dx, y + dy, z + dz);
                        let theirs = &self.entries
                            [self.starts[there] as usize..self.starts[there + 1] as usize];
                        for &i in mine {
                            for &j in theirs {
                                consider(i as usize, j as usize, &mut visit);
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Thirteen of the twenty-six surrounding cells, chosen so that every unordered pair of
/// adjacent cells is generated exactly once when the whole grid is swept.
///
/// The rule is lexicographic: keep an offset if it is "greater than zero" reading `(dz, dy,
/// dx)` in that order. Its negation is then excluded, and the cell on the other side reaches
/// this one instead.
const HALF_NEIGHBOURS: [(isize, isize, isize); 13] = [
    (1, 0, 0),
    (-1, 1, 0),
    (0, 1, 0),
    (1, 1, 0),
    (-1, -1, 1),
    (0, -1, 1),
    (1, -1, 1),
    (-1, 0, 1),
    (0, 0, 1),
    (1, 0, 1),
    (-1, 1, 1),
    (0, 1, 1),
    (1, 1, 1),
];

#[cfg(test)]
mod tests {
    use super::*;
    use dualis_core::Rng;
    use std::collections::HashSet;

    #[test]
    fn wrapping_and_the_minimum_image_agree_about_where_things_are() {
        let b = PeriodicBox::cubic(10.0);
        assert_eq!(
            b.wrap(DVec3::new(11.0, -1.0, 25.0)),
            DVec3::new(1.0, 9.0, 5.0)
        );
        assert_eq!(
            b.wrap(DVec3::new(0.0, 10.0, 9.999)),
            DVec3::new(0.0, 0.0, 9.999)
        );

        // Two particles a whisker either side of a face are neighbours, not nine apart.
        let d = b.shortest(DVec3::new(9.5, 0.0, 0.0), DVec3::new(0.5, 0.0, 0.0));
        assert_eq!(d, DVec3::new(-1.0, 0.0, 0.0));
        assert!((d.length() - 1.0).abs() < 1e-15);

        // Antisymmetric, which is what makes Newton's third law expressible.
        let (p, q) = (DVec3::new(1.0, 2.0, 3.0), DVec3::new(8.0, 9.0, 0.5));
        assert_eq!(b.shortest(p, q), -b.shortest(q, p));
        // And never longer than the half-diagonal.
        assert!(b.shortest(p, q).length() <= b.length * 3f64.sqrt() / 2.0 + 1e-12);
    }

    /// The convention has a domain, and the box says where it ends rather than leaving a
    /// caller to double-count.
    #[test]
    fn a_cutoff_past_half_the_box_is_refused() {
        let b = PeriodicBox::cubic(10.0);
        assert!(b.admits(4.999));
        assert!(
            !b.admits(5.0),
            "at exactly L/2 the nearest image is a coin toss"
        );
        assert!(!b.admits(7.0));

        // A density constructor that agrees with itself.
        let sized = PeriodicBox::for_density(864, 0.8442);
        assert!((864.0 / sized.volume() - 0.8442).abs() < 1e-12);
    }

    /// **The property the whole cell list exists to have.** It finds exactly the pairs a
    /// brute-force sweep finds, each of them once.
    ///
    /// Checked against the `O(N²)` loop rather than against itself, at a density high enough
    /// that cells are genuinely full, and with a box small enough that images matter.
    #[test]
    fn cells_find_the_same_pairs_as_every_pair() {
        let bounds = PeriodicBox::cubic(12.0);
        let cutoff = 2.5;
        let mut rng = Rng::new(0x0CE1_1157);
        let positions: Vec<DVec3> = (0..400)
            .map(|_| {
                DVec3::new(
                    rng.range(0.0, 12.0),
                    rng.range(0.0, 12.0),
                    rng.range(0.0, 12.0),
                )
            })
            .collect();

        let mut by_cells: Vec<(usize, usize)> = Vec::new();
        let list = CellList::build(bounds, cutoff, &positions);
        assert!(
            list.divisions() >= 4,
            "the test needs real cells, got {}",
            list.divisions()
        );
        assert!(
            list.cell_size() >= cutoff,
            "a cell must be at least a cutoff wide"
        );
        list.for_each_pair(bounds, cutoff, &positions, |i, j, _, _| {
            by_cells.push((i.min(j), i.max(j)));
        });

        let mut brute: Vec<(usize, usize)> = Vec::new();
        for i in 0..positions.len() {
            for j in (i + 1)..positions.len() {
                if bounds.shortest(positions[i], positions[j]).length_squared() < cutoff * cutoff {
                    brute.push((i, j));
                }
            }
        }

        assert!(
            !brute.is_empty(),
            "nothing was in range, so nothing was tested"
        );
        let seen: HashSet<_> = by_cells.iter().copied().collect();
        assert_eq!(seen.len(), by_cells.len(), "a pair was visited twice");
        assert_eq!(seen, brute.iter().copied().collect::<HashSet<_>>());
    }

    /// A box too small to divide falls back to every pair rather than to a wrong answer.
    #[test]
    fn a_box_under_three_cells_wide_falls_back_to_every_pair() {
        let bounds = PeriodicBox::cubic(5.0);
        let positions: Vec<DVec3> = (0..40)
            .map(|k| {
                let f = k as f64 / 40.0;
                DVec3::new(f * 5.0, (f * 3.7) % 5.0, (f * 1.3) % 5.0)
            })
            .collect();
        // 5 / 2.4 floors to 2, which is under three, so one cell.
        let list = CellList::build(bounds, 2.4, &positions);
        assert_eq!(list.divisions(), 1);

        let mut count = 0;
        list.for_each_pair(bounds, 2.4, &positions, |_, _, _, _| count += 1);
        let brute = (0..40)
            .flat_map(|i| ((i + 1)..40).map(move |j| (i, j)))
            .filter(|(i, j)| {
                bounds
                    .shortest(positions[*i], positions[*j])
                    .length_squared()
                    < 2.4 * 2.4
            })
            .count();
        assert_eq!(count, brute);
    }

    /// The separation handed to `visit` is the minimum image one, pointing from `j` to `i`.
    /// A sign error here would invert every force in the simulation.
    #[test]
    fn the_separation_points_from_the_second_to_the_first() {
        let bounds = PeriodicBox::cubic(9.0);
        let positions = vec![DVec3::new(1.0, 1.0, 1.0), DVec3::new(2.5, 1.0, 1.0)];
        let list = CellList::build(bounds, 2.0, &positions);
        let mut seen = 0;
        list.for_each_pair(bounds, 2.0, &positions, |i, j, d, r2| {
            seen += 1;
            let expected = bounds.shortest(positions[i], positions[j]);
            assert_eq!(d, expected);
            assert!((r2 - 1.5 * 1.5).abs() < 1e-12);
        });
        assert_eq!(seen, 1);
    }
}
