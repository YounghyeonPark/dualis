//! Sound in two dimensions: a room rather than a tube.
//!
//! [`Tube`](crate::Tube) is one-dimensional, which covers a duct and an organ pipe and
//! nothing with a shape. A room has modes in every direction at once, and they interact
//! in a way a one-dimensional model cannot show: the frequencies are not a harmonic
//! series, they crowd together as they rise, and the low ones are far enough apart to be
//! heard individually as the boom a small room has at one particular note.
//!
//! # The closed form, which is exact and not obvious
//!
//! A rectangular room with rigid walls resonates at
//!
//! ```text
//! f(nx, ny) = (c/2) √((nx/Lx)² + (ny/Ly)²)
//! ```
//!
//! for every pair of non-negative integers. Two things follow that a tube does not show.
//! The series is not harmonic — `f(1,1)` is not a whole multiple of `f(1,0)` — so a room
//! does not ring on a note, it rings on a chord that is not one. And the modes get denser
//! as the frequency rises, going as `f²` in two dimensions rather than staying evenly
//! spaced, which is why a room's colouration is audible at the bottom of the spectrum and
//! not at the top.
//!
//! # What this shares with the tube, and what it does not
//!
//! The same staggered scheme, extended: pressures at cell centres, velocities on the
//! faces between them in each direction. The CFL condition tightens, and by a factor that
//! is worth knowing — `dt ≤ dx/(c√2)` in two dimensions rather than `dx/c`, because a
//! wave travelling diagonally covers `√2` cells while crossing one. In three dimensions
//! it would be `√3`. Every explicit wave solver pays that, and it is the reason
//! three-dimensional acoustics is expensive rather than merely large.

use dualis_core::conserved::quantity;
use dualis_core::{Domain, Exchange, Kind, Ledger, Violation};
use dualis_units::{Area, Density, Energy, Frequency, Length, Pressure, Time, Velocity};

/// A rectangular room, discretised on a uniform grid with rigid walls.
pub struct Room {
    name: &'static str,
    /// Pressure at cell centres, row-major.
    pressure: Vec<f64>,
    /// Pressure one whole step earlier, kept only so that [`Room::energy`] can be the
    /// quantity the scheme actually conserves rather than one that wobbles.
    pressure_prev: Vec<f64>,
    /// Velocity on the faces between horizontally adjacent cells: `(nx - 1)` per row.
    vx: Vec<f64>,
    /// Velocity on the faces between vertically adjacent cells: `(ny - 1)` rows.
    vy: Vec<f64>,
    nx: usize,
    ny: usize,
    dx: f64,
    speed: f64,
    density: f64,
    /// Depth, so a two-dimensional room still has a volume and an energy in joules.
    depth: f64,
    saved: Option<(Vec<f64>, Vec<f64>, Vec<f64>)>,
}

impl Room {
    /// A room of the given size, at rest, with rigid walls.
    ///
    /// The grid spacing is taken from the width; the height is quantised to the nearest
    /// whole number of the same cells, so the cells stay square and the CFL limit stays
    /// isotropic.
    pub fn new(
        name: &'static str,
        width: Length,
        height: Length,
        cells_across: usize,
        depth: Length,
        density: Density,
        sound_speed: Velocity,
    ) -> Room {
        let nx = cells_across.max(3);
        let dx = width.to_si() / (nx - 1) as f64;
        let ny = ((height.to_si() / dx).round() as usize + 1).max(3);
        Room {
            name,
            pressure: vec![0.0; nx * ny],
            pressure_prev: vec![0.0; nx * ny],
            vx: vec![0.0; (nx - 1) * ny],
            vy: vec![0.0; nx * (ny - 1)],
            nx,
            ny,
            dx,
            speed: sound_speed.to_si(),
            density: density.to_si(),
            depth: depth.to_si(),
            saved: None,
        }
    }

    /// A room of air at 20 °C, one metre deep.
    pub fn of_air(name: &'static str, width: Length, height: Length, cells_across: usize) -> Room {
        Room::new(
            name,
            width,
            height,
            cells_across,
            Length::m(1.0),
            Density::kg_per_m3(1.204),
            Velocity::m_per_s(343.0),
        )
    }

    /// Start with a pressure field and no motion.
    pub fn released_from(mut self, profile: impl Fn(Length, Length) -> Pressure) -> Room {
        for j in 0..self.ny {
            for i in 0..self.nx {
                let x = Length::from_si(i as f64 * self.dx);
                let y = Length::from_si(j as f64 * self.dx);
                self.pressure[j * self.nx + i] = profile(x, y).to_si();
            }
        }
        self.pressure_prev = self.pressure.clone();
        self.vx.iter_mut().for_each(|v| *v = 0.0);
        self.vy.iter_mut().for_each(|v| *v = 0.0);
        self
    }

    /// Excite one rigid-wall mode exactly: `cos(nx π x/Lx) cos(ny π y/Ly)`.
    ///
    /// The shape whose oscillation frequency the closed form predicts, so that the
    /// prediction can be checked against the integration rather than against itself.
    pub fn released_in_mode(self, nx: u32, ny: u32, amplitude: Pressure) -> Room {
        let (lx, ly) = (self.width().to_si(), self.height().to_si());
        let a = amplitude.to_si();
        self.released_from(|x, y| {
            let px = if lx > 0.0 {
                (nx as f64 * std::f64::consts::PI * x.to_si() / lx).cos()
            } else {
                1.0
            };
            let py = if ly > 0.0 {
                (ny as f64 * std::f64::consts::PI * y.to_si() / ly).cos()
            } else {
                1.0
            };
            Pressure::from_si(a * px * py)
        })
    }

    pub fn width(&self) -> Length {
        Length::from_si((self.nx - 1) as f64 * self.dx)
    }

    pub fn height(&self) -> Length {
        Length::from_si((self.ny - 1) as f64 * self.dx)
    }

    pub fn cells(&self) -> (usize, usize) {
        (self.nx, self.ny)
    }

    pub fn pressure_at(&self, i: usize, j: usize) -> Pressure {
        let i = i.min(self.nx - 1);
        let j = j.min(self.ny - 1);
        Pressure::from_si(self.pressure[j * self.nx + i])
    }

    pub fn peak_pressure(&self) -> Pressure {
        Pressure::from_si(self.pressure.iter().fold(0.0f64, |a, p| a.max(p.abs())))
    }

    /// Acoustic energy, in the form the scheme actually conserves.
    ///
    /// The obvious expression — `p²` at one time level plus `u²` at another — is not it.
    /// On a staggered grid the pressures sit at whole steps and the velocities half a
    /// step later, so reading both at once carries an `O(ω dt)` wobble: measured at 4% for
    /// a pulse at five substeps a millisecond, which does not drift but is more than
    /// enough to make a conservation audit meaningless.
    ///
    /// The time-centred form uses the *product* of the pressure either side of the step,
    /// `pⁿ pⁿ⁺¹`, against the velocity at the half step between them. That combination is
    /// exactly conserved by the lossless leapfrog rather than approximately, which is what
    /// lets the audit be tight. It costs one extra array.
    ///
    /// A consequence worth knowing: a single cell's contribution can be negative, where
    /// the pressure changed sign across the step. The total cannot.
    pub fn energy(&self) -> Energy {
        let volume = self.dx * self.dx * self.depth;
        let rc2 = self.density * self.speed * self.speed;
        if rc2 <= 0.0 {
            return Energy::from_si(0.0);
        }
        let potential: f64 = self
            .pressure
            .iter()
            .zip(self.pressure_prev.iter())
            .map(|(p, prev)| p * prev / (2.0 * rc2) * volume)
            .sum();
        let kinetic: f64 = self
            .vx
            .iter()
            .chain(self.vy.iter())
            .map(|u| self.density * u * u / 2.0 * volume)
            .sum();
        Energy::from_si(potential + kinetic)
    }

    /// Cross-sectional area of one cell, for reporting.
    pub fn cell_area(&self) -> Area {
        Area::from_si(self.dx * self.depth)
    }

    /// The `(nx, ny)` rigid-wall mode frequency: `(c/2)√((nx/Lx)² + (ny/Ly)²)`.
    ///
    /// Exact, and not a harmonic series: `f(1,1)` of a square room is `√2` times
    /// `f(1,0)`, which is a tritone and not a note anyone would call in tune.
    pub fn mode_frequency(&self, nx: u32, ny: u32) -> Frequency {
        let (lx, ly) = (self.width().to_si(), self.height().to_si());
        if lx <= 0.0 || ly <= 0.0 {
            return Frequency::from_si(0.0);
        }
        let a = nx as f64 / lx;
        let b = ny as f64 / ly;
        Frequency::from_si(self.speed / 2.0 * (a * a + b * b).sqrt())
    }

    /// Every mode below a frequency, in ascending order.
    ///
    /// The count grows as `f²` in two dimensions — the modes fill an area of the `(nx,
    /// ny)` lattice — which is why a room's individual resonances are audible at the
    /// bottom of the range and merge into a statistical hiss further up.
    pub fn modes_below(&self, limit: Frequency) -> Vec<(u32, u32, Frequency)> {
        let mut found = Vec::new();
        // Bound the search from the largest index either axis could contribute.
        let max_index = |l: f64| {
            if l <= 0.0 {
                0
            } else {
                (2.0 * limit.to_si() * l / self.speed).floor() as u32
            }
        };
        let (mx, my) = (
            max_index(self.width().to_si()),
            max_index(self.height().to_si()),
        );
        for nx in 0..=mx {
            for ny in 0..=my {
                if nx == 0 && ny == 0 {
                    continue;
                }
                let f = self.mode_frequency(nx, ny);
                if f <= limit {
                    found.push((nx, ny, f));
                }
            }
        }
        found.sort_by(|a, b| a.2.to_si().total_cmp(&b.2.to_si()));
        found
    }

    /// Courant number, `c dt √2/dx` in two dimensions.
    pub fn courant(&self, dt: Time) -> f64 {
        self.speed * dt.to_si() * 2f64.sqrt() / self.dx
    }
}

impl Domain for Room {
    fn name(&self) -> &'static str {
        self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// `dx/(c√2)`, tighter than the tube's `dx/c` by exactly the diagonal.
    ///
    /// A wave going corner to corner covers `√2` cells while crossing one, and the
    /// three-point stencil in each direction cannot know about it. In three dimensions the
    /// factor is `√3`, so the same room at the same resolution costs 22% more steps again
    /// on top of the extra dimension — which is the real reason three-dimensional
    /// acoustics is dear.
    fn max_stable_dt(&self, _now: Time) -> Time {
        if self.speed <= 0.0 {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(self.dx / (self.speed * 2f64.sqrt()))
    }

    fn step(&mut self, _t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let courant = self.courant(dt);
        if courant > 1.0 + 1e-12 {
            return Err(Violation {
                quantity: "Courant number".to_string(),
                site: format!("{} (two-dimensional wave equation)", self.name),
                before: 1.0,
                after: courant,
                scale: 1.0,
                tolerance: 1e-12,
            });
        }
        let h = dt.to_si();
        if h <= 0.0 {
            return Ok(());
        }
        let rc2 = self.density * self.speed * self.speed;
        let (nx, ny) = (self.nx, self.ny);

        // Faces: rho du/dt = -grad p.
        for j in 0..ny {
            for i in 0..nx - 1 {
                let face = j * (nx - 1) + i;
                self.vx[face] -= h / (self.density * self.dx)
                    * (self.pressure[j * nx + i + 1] - self.pressure[j * nx + i]);
            }
        }
        for j in 0..ny - 1 {
            for i in 0..nx {
                let face = j * nx + i;
                self.vy[face] -= h / (self.density * self.dx)
                    * (self.pressure[(j + 1) * nx + i] - self.pressure[j * nx + i]);
            }
        }

        // Cells: dp/dt = -rho c^2 div u. Rigid walls mean no velocity through them, so
        // the flux outside the domain is zero — which is the boundary condition and needs
        // no ghost cells.
        self.pressure_prev.copy_from_slice(&self.pressure);
        for j in 0..ny {
            for i in 0..nx {
                let left = if i == 0 {
                    0.0
                } else {
                    self.vx[j * (nx - 1) + i - 1]
                };
                let right = if i == nx - 1 {
                    0.0
                } else {
                    self.vx[j * (nx - 1) + i]
                };
                let below = if j == 0 {
                    0.0
                } else {
                    self.vy[(j - 1) * nx + i]
                };
                let above = if j == ny - 1 {
                    0.0
                } else {
                    self.vy[j * nx + i]
                };
                self.pressure[j * nx + i] -= h * rc2 / self.dx * ((right - left) + (above - below));
            }
        }
        Ok(())
    }

    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.energy().to_si())
    }

    fn checkpoint(&mut self) {
        self.saved = Some((self.pressure.clone(), self.vx.clone(), self.vy.clone()));
    }

    fn restore(&mut self) {
        if let Some((p, vx, vy)) = self.saved.clone() {
            self.pressure_prev.copy_from_slice(&p);
            self.pressure = p;
            self.vx = vx;
            self.vy = vy;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square(cells: usize) -> Room {
        Room::of_air("room", Length::m(4.0), Length::m(4.0), cells)
    }

    /// The closed form, and the two things about it a tube cannot show.
    #[test]
    fn room_modes_are_not_a_harmonic_series() {
        let room = Room::of_air("room", Length::m(5.0), Length::m(3.0), 51);
        assert!((room.width().to_si() - 5.0).abs() < 1e-12);
        // The height is quantised to whole cells, so it lands near 3 m rather than on it.
        assert!(
            (room.height().to_si() - 3.0).abs() < 0.11,
            "{}",
            room.height().to_si()
        );

        // (c/2) sqrt((1/5)^2) = 34.3 Hz along the length.
        assert!((room.mode_frequency(1, 0).to_si() - 34.3).abs() < 0.1);
        // And (c/2)(1/3) = 57.2 Hz across the width.
        assert!(
            (room.mode_frequency(0, 1).to_si() - 57.2).abs() < 1.5,
            "{}",
            room.mode_frequency(0, 1).to_si()
        );

        // The oblique mode is the quadrature sum, not the arithmetic one — which is what
        // makes the series inharmonic.
        let f10 = room.mode_frequency(1, 0).to_si();
        let f01 = room.mode_frequency(0, 1).to_si();
        let f11 = room.mode_frequency(1, 1).to_si();
        assert!(
            (f11 - (f10 * f10 + f01 * f01).sqrt()).abs() < 1e-9,
            "the oblique mode is the quadrature sum"
        );
        assert!(f11 < f10 + f01, "and so it is below the arithmetic sum");

        // In a square room the diagonal mode is exactly root two times the axial one: a
        // tritone, which is why a square room sounds wrong rather than merely loud.
        let sq = square(41);
        assert!((sq.mode_frequency(1, 1) / sq.mode_frequency(1, 0) - 2f64.sqrt()).abs() < 1e-12);
        assert_eq!(sq.mode_frequency(0, 0).to_si(), 0.0);
    }

    /// The modes crowd together as they rise, going as `f²` in two dimensions. That is
    /// why a small room's colouration is heard at the bottom of the range and turns into
    /// an even hiss further up.
    #[test]
    fn the_modes_get_denser_as_the_frequency_rises() {
        let room = Room::of_air("room", Length::m(5.0), Length::m(4.0), 51);
        let count = |f: f64| room.modes_below(Frequency::hz(f)).len();

        // Doubling the frequency should roughly quadruple the count.
        let (low, high) = (count(200.0), count(400.0));
        let ratio = high as f64 / low as f64;
        assert!(
            (ratio - 4.0).abs() < 1.0,
            "the count should go as f squared: {low} below 200 Hz and {high} below 400, \
             a ratio of {ratio:.2}"
        );

        // They come out in ascending order and the first is the longest dimension's.
        let modes = room.modes_below(Frequency::hz(120.0));
        assert!(modes.len() > 4);
        for pair in modes.windows(2) {
            assert!(pair[0].2 <= pair[1].2, "modes should be sorted");
        }
        assert_eq!((modes[0].0, modes[0].1), (1, 0));
        // Nothing below the fundamental, and the zero mode is not a mode.
        assert!(room.modes_below(Frequency::hz(30.0)).is_empty());
    }

    /// The CFL condition tightens by exactly the diagonal, and going past it is refused.
    #[test]
    fn two_dimensions_tighten_the_courant_limit_by_root_two() {
        let room = square(41);
        // 4 m over 40 intervals is 100 mm; 343 m/s crosses that in 292 us, and the
        // two-dimensional limit is that over root two.
        let limit = room.max_stable_dt(Time::ZERO);
        assert!(
            (limit.in_us() - 206.2).abs() < 0.5,
            "limit {} us",
            limit.in_us()
        );
        let tube_limit = 0.1 / 343.0;
        assert!(
            (tube_limit / limit.to_si() - 2f64.sqrt()).abs() < 1e-9,
            "exactly root two tighter than one dimension"
        );
        assert!((room.courant(limit) - 1.0).abs() < 1e-12);

        let mut room = square(41);
        let mut bus = Exchange::new();
        assert!(room.step(Time::ZERO, limit, &mut bus).is_ok());
        let err = room
            .step(Time::ZERO, limit * 1.01, &mut bus)
            .expect_err("past the limit must be refused");
        assert_eq!(err.quantity, "Courant number");
    }

    /// A mode excited exactly should oscillate at exactly its predicted frequency: half a
    /// period inverts it, a whole one restores it.
    ///
    /// This is what ties the closed form to the code actually stepping, and it is the test
    /// that would fail if the stencil, the staggering or the wall condition were wrong.
    #[test]
    fn an_excited_mode_oscillates_at_its_predicted_frequency() {
        for (nx, ny) in [(1u32, 0u32), (1, 1), (2, 1)] {
            let mut room = square(81).released_in_mode(nx, ny, Pressure::from_si(100.0));
            let f = room.mode_frequency(nx, ny).to_si();
            let period = 1.0 / f;
            // Well inside the limit, so the numerical dispersion stays small.
            let dt = room.max_stable_dt(Time::ZERO) * 0.5;
            let steps = (period / dt.to_si()).round() as u32;
            let start = room.pressure_at(0, 0).to_si();
            assert!(
                start.abs() > 50.0,
                "the corner is an antinode of every mode"
            );

            let mut bus = Exchange::new();
            for _ in 0..steps / 2 {
                room.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            let half = room.pressure_at(0, 0).to_si();
            assert!(
                half < -0.85 * start,
                "mode ({nx},{ny}): half a period should invert it, {start:.1} to {half:.1}"
            );

            for _ in 0..steps - steps / 2 {
                room.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            let whole = room.pressure_at(0, 0).to_si();
            assert!(
                (whole / start - 1.0).abs() < 0.1,
                "mode ({nx},{ny}): a whole period should restore it, {start:.1} to \
                 {whole:.1}"
            );
        }
    }

    /// Rigid walls keep the energy in, so it holds rather than draining.
    #[test]
    fn a_closed_room_keeps_its_energy() {
        let mut room = square(61).released_from(|x, y| {
            let (dx, dy) = (x.to_si() - 2.0, y.to_si() - 2.0);
            let r2 = dx * dx + dy * dy;
            Pressure::from_si(200.0 * (-r2 / 0.09).exp())
        });
        let start = room.energy().to_si();
        assert!(start > 0.0);

        let dt = room.max_stable_dt(Time::ZERO) * 0.5;
        let mut bus = Exchange::new();
        for _ in 0..3000 {
            room.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        let now = room.energy().to_si();
        assert!(
            (now / start - 1.0).abs() < 0.02,
            "energy went from {start:e} to {now:e}"
        );
        // Rigid walls publish nothing.
        assert!(bus.peek(quantity::ENERGY).abs() < 1e-30);
    }

    /// A pulse in the middle spreads outwards in a circle rather than along the axes,
    /// which is the thing a pair of one-dimensional tubes could not produce.
    #[test]
    fn a_pulse_spreads_isotropically() {
        let mut room = square(81).released_from(|x, y| {
            let (dx, dy) = (x.to_si() - 2.0, y.to_si() - 2.0);
            let r2 = dx * dx + dy * dy;
            Pressure::from_si(200.0 * (-r2 / 0.02).exp())
        });
        let dt = room.max_stable_dt(Time::ZERO) * 0.5;
        let mut bus = Exchange::new();
        // Long enough for the front to travel 0.8 m, well short of the walls.
        let steps = (0.8 / 343.0 / dt.to_si()).round() as u32;
        for _ in 0..steps {
            room.step(Time::ZERO, dt, &mut bus).unwrap();
        }

        // Sample the same distance from the centre along an axis and along the diagonal.
        let (cx, cy) = (40usize, 40usize);
        let offset = 16usize; // 0.8 m at 50 mm a cell
        let diag = (offset as f64 / 2f64.sqrt()).round() as usize;
        let along_axis = room.pressure_at(cx + offset, cy).to_si().abs();
        let along_diagonal = room.pressure_at(cx + diag, cy + diag).to_si().abs();
        assert!(
            along_axis > 1.0,
            "the front should have arrived: {along_axis}"
        );
        assert!(
            (along_axis - along_diagonal).abs() / along_axis.max(along_diagonal) < 0.25,
            "a circular front should reach both at once: {along_axis:.3} along the axis \
             and {along_diagonal:.3} along the diagonal"
        );
    }

    /// Under the scheduler, subcycling to the two-dimensional limit.
    #[test]
    fn the_domain_runs_under_the_scheduler() {
        use dualis_core::{Schedule, Simulation};

        let room = square(41).released_from(|x, y| {
            let (dx, dy) = (x.to_si() - 2.0, y.to_si() - 2.0);
            Pressure::from_si(200.0 * (-(dx * dx + dy * dy) / 0.5).exp())
        });
        // Tight, because `Room::energy` reports the time-centred quantity the leapfrog
        // conserves exactly rather than the one that wobbles by a few percent. The first
        // version of that function needed 5e-2 here and was still failing.
        let mut sim = Simulation::new(Schedule::Multirate)
            .conservation_tolerance(1e-9)
            .with(room);
        let report = sim.advance(Time::ms(1.0)).unwrap();
        assert_eq!(report.substeps[0].0, "room");
        // 1 ms at a 206 us limit is 5 substeps.
        assert_eq!(report.substeps[0].1, 5);
        for _ in 0..20 {
            sim.advance(Time::ms(1.0)).expect("energy should hold");
        }
    }

    /// Degenerate rooms are handled.
    #[test]
    fn degenerate_rooms_are_handled() {
        let tiny = Room::of_air("tiny", Length::m(0.1), Length::m(0.1), 1);
        let (nx, ny) = tiny.cells();
        assert!(nx >= 3 && ny >= 3);

        let mut silent = Room::new(
            "silent",
            Length::m(1.0),
            Length::m(1.0),
            11,
            Length::m(1.0),
            Density::kg_per_m3(1.0),
            Velocity::m_per_s(0.0),
        );
        assert!(!silent.max_stable_dt(Time::ZERO).to_si().is_finite());
        let mut bus = Exchange::new();
        assert!(silent.step(Time::ZERO, Time::s(1.0), &mut bus).is_ok());
        assert_eq!(silent.peak_pressure().to_si(), 0.0);

        // A zero step changes nothing.
        let mut room = square(21).released_from(|_, _| Pressure::from_si(7.0));
        let before = room.pressure_at(5, 5);
        room.step(Time::ZERO, Time::ZERO, &mut bus).unwrap();
        assert_eq!(room.pressure_at(5, 5), before);
    }
}
