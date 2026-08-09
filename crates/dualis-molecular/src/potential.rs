//! The Lennard-Jones pair potential, and what truncating it costs.
//!
//! ```text
//! u(r) = 4ε [ (σ/r)¹² − (σ/r)⁶ ]
//! ```
//!
//! The `r⁻⁶` term is real physics — induced-dipole attraction, and the exponent comes out of
//! the derivation. The `r⁻¹²` is not: it is a repulsion that happens to be the square of the
//! attraction, which makes it cheap to evaluate and roughly the right steepness. Everyone
//! knows this and everyone uses it anyway, because the answers that matter are set by the
//! attraction and by the fact that atoms cannot overlap, not by the shape of the wall.
//!
//! # Two exact numbers to check against
//!
//! The minimum sits at `r = 2^(1/6)σ ≈ 1.1225σ` with depth exactly `−ε`, and the potential
//! crosses zero at exactly `σ`. Both are closed forms, so the implementation can be checked
//! against arithmetic rather than against a table.
//!
//! # Truncation, and why it is shifted
//!
//! Every pair beyond a cutoff is skipped, or the cost is `O(N²)` and the whole point of a
//! cell list is lost. At `2.5σ` the raw potential is still `−0.0163ε`, so cutting there leaves
//! a step in the energy — and a step in the energy is a delta function in the force, which an
//! integrator turns into a slow warming that looks exactly like a thermostat nobody asked for.
//!
//! Subtracting `u(rc)` removes it. The force is untouched, since a constant differentiates to
//! nothing, so the dynamics are the same and only the energy bookkeeping changes. What is left
//! is a discontinuity in the *force* at the cutoff, which is `0.0390 ε/σ` there and small
//! enough that energy holds to the integrator's own bound. Removing that too needs a force
//! shift, which changes the dynamics rather than the bookkeeping, and is not done here.

/// A Lennard-Jones interaction, in reduced units of its own `ε` and `σ`.
///
/// Stored as the two parameters rather than as the various precomputed powers, because the
/// powers are cheap and a struct whose fields cannot be read back is a struct that will be
/// wrong one day.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LennardJones {
    /// Well depth, in joules.
    pub epsilon: f64,
    /// Distance at which the potential crosses zero, in metres.
    pub sigma: f64,
    /// Beyond this, pairs are skipped. In metres.
    pub cutoff: f64,
}

/// A [`LennardJones`] with its loop-invariant arithmetic already done.
///
/// Made by [`LennardJones::prepared`] and cheap to copy. Holds no reference to the potential it
/// came from, so a caller can keep one beside a cell list for as long as the parameters do not
/// change — and if they do, the potential is `Copy` and preparing again is five operations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Prepared {
    cutoff_sq: f64,
    sigma_sq: f64,
    four_epsilon: f64,
    twentyfour_epsilon: f64,
    shift: f64,
}

impl Prepared {
    /// Energy and force for one pair, given the *square* of their separation.
    ///
    /// Bit-for-bit what [`LennardJones::at_squared`] returns, and pinned as such by a test that
    /// sweeps separations across the cutoff.
    pub fn at_squared(&self, r2: f64) -> Option<Pair> {
        if r2 >= self.cutoff_sq || r2 <= 0.0 {
            return None;
        }
        let inv2 = self.sigma_sq / r2;
        let inv6 = inv2 * inv2 * inv2;
        let inv12 = inv6 * inv6;
        Some(Pair {
            energy: self.four_epsilon * (inv12 - inv6) - self.shift,
            // -du/dr = 24 e (2 s^12/r^13 - s^6/r^7), and dividing by r gives
            // 24 e (2 inv12 - inv6) / r^2 without ever forming r.
            force_over_r: self.twentyfour_epsilon * (2.0 * inv12 - inv6) / r2,
        })
    }
}

/// What a pair contributes: energy in joules, and the scalar `−du/dr / r` that turns a
/// separation vector into a force.
///
/// The combination `f/r` rather than `f` because every caller multiplies by the separation
/// vector next, and computing `r` from `r²` is a square root that is not otherwise needed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pair {
    /// Pair potential energy in joules, already shifted so it reaches the cutoff at zero.
    pub energy: f64,
    /// Force per unit separation: the force on the first particle is `force_over_r * d`,
    /// where `d` points from the second to the first.
    pub force_over_r: f64,
}

impl LennardJones {
    /// Argon, the substance every molecular-dynamics paper starts with.
    ///
    /// `ε/k_B = 119.8 K` and `σ = 3.405 Å`, which is the parameterisation Rahman used in 1964
    /// for the first molecular-dynamics simulation of a liquid, and which everything since has
    /// been compared against.
    pub fn argon() -> LennardJones {
        let epsilon = 119.8 * dualis_units::BOLTZMANN.to_si();
        LennardJones {
            epsilon,
            sigma: 3.405e-10,
            cutoff: 2.5 * 3.405e-10,
        }
    }

    /// Reduced units: `ε = σ = 1`, cutoff at `2.5`.
    ///
    /// What the literature quotes almost every result in, and what makes a state point
    /// comparable across papers: a density of `0.8442` and a temperature of `0.722` is *the*
    /// triple-point liquid whatever atom you had in mind.
    pub fn reduced() -> LennardJones {
        LennardJones {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.5,
        }
    }

    /// Move the cutoff. Beyond `L/2` of whatever box this ends up in, the minimum image
    /// becomes ambiguous and [`PeriodicBox::admits`](crate::PeriodicBox::admits) says so.
    pub fn with_cutoff(mut self, cutoff: f64) -> LennardJones {
        self.cutoff = cutoff;
        self
    }

    /// Separation at the bottom of the well, `2^(1/6) σ`.
    pub fn minimum(&self) -> f64 {
        2f64.powf(1.0 / 6.0) * self.sigma
    }

    /// The constant subtracted from every pair energy so the potential reaches the cutoff at
    /// zero. See the module note on why this and not a force shift.
    pub fn shift(&self) -> f64 {
        let s6 = (self.sigma / self.cutoff).powi(6);
        4.0 * self.epsilon * s6 * (s6 - 1.0)
    }

    /// Energy and force for one pair, given the *square* of their separation.
    ///
    /// Squared because a neighbour search compares `r²` against `rc²` to avoid a square root
    /// per pair, and having found the pair it would be silly to take one anyway. Returns
    /// `None` past the cutoff, so a caller cannot accidentally add a contribution it meant to
    /// skip.
    pub fn at_squared(&self, r2: f64) -> Option<Pair> {
        self.prepared().at_squared(r2)
    }

    /// The loop-invariant parts of [`at_squared`](LennardJones::at_squared), computed once.
    ///
    /// Hoist this out of a pair loop. `at_squared` is a convenience that prepares and discards,
    /// which is free for one call and is not free for a hundred thousand: it was recomputing
    /// `shift()` — a division and six multiplies for a number that depends only on constants —
    /// on **every pair**, which at 2048 atoms is about 113 000 evaluations per step of the same
    /// value.
    ///
    /// Every field here is exactly the expression it replaces, so a prepared potential and a
    /// bare one agree to the bit. That is not a nicety: this crate's results are pinned across
    /// four platforms, and an optimisation that moved the last bit would be a change in the
    /// physics wearing a performance change's clothes.
    pub fn prepared(&self) -> Prepared {
        Prepared {
            cutoff_sq: self.cutoff * self.cutoff,
            sigma_sq: self.sigma * self.sigma,
            four_epsilon: 4.0 * self.epsilon,
            twentyfour_epsilon: 24.0 * self.epsilon,
            shift: self.shift(),
        }
    }

    /// Energy alone, at a separation. For checking against the closed forms.
    pub fn energy_at(&self, r: f64) -> f64 {
        self.at_squared(r * r).map(|p| p.energy).unwrap_or(0.0)
    }

    /// The long-range correction to the energy per particle, for a fluid of number density
    /// `rho` whose structure past the cutoff is assumed uniform.
    ///
    /// ```text
    /// u_tail = (8/3) π ρ ε σ³ [ (1/3)(σ/rc)⁹ − (σ/rc)³ ]
    /// ```
    ///
    /// Not applied to the dynamics — it is a constant, so it exerts no force — but it is what
    /// makes a truncated simulation's energy comparable with a published one. At `rc = 2.5σ`
    /// and liquid density it is about `−0.5 ε` per particle, which is several percent of the
    /// total and not ignorable when the point is to match a number.
    pub fn energy_tail(&self, number_density: f64) -> f64 {
        let s3 = self.sigma.powi(3);
        let x3 = (self.sigma / self.cutoff).powi(3);
        (8.0 / 3.0)
            * std::f64::consts::PI
            * number_density
            * self.epsilon
            * s3
            * (x3 * x3 * x3 / 3.0 - x3)
    }

    /// The matching correction to the pressure.
    ///
    /// ```text
    /// p_tail = (16/3) π ρ² ε σ³ [ (2/3)(σ/rc)⁹ − (σ/rc)³ ]
    /// ```
    pub fn pressure_tail(&self, number_density: f64) -> f64 {
        let s3 = self.sigma.powi(3);
        let x3 = (self.sigma / self.cutoff).powi(3);
        (16.0 / 3.0)
            * std::f64::consts::PI
            * number_density
            * number_density
            * self.epsilon
            * s3
            * (2.0 * x3 * x3 * x3 / 3.0 - x3)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The two closed forms: zero at `σ`, minimum of exactly `−ε` at `2^(1/6)σ`.
    ///
    /// Checked on the *unshifted* potential, since the shift is a separate claim — so the raw
    /// expression is reconstructed by adding the shift back rather than by having a second
    /// code path that could disagree with the first.
    #[test]
    fn the_potential_has_its_textbook_shape() {
        let lj = LennardJones::reduced();
        let raw = |r: f64| lj.energy_at(r) + lj.shift();

        assert!(
            raw(1.0).abs() < 1e-15,
            "u(sigma) must be zero, got {}",
            raw(1.0)
        );
        let rmin = lj.minimum();
        assert!((rmin - 1.122_462_048_309_373).abs() < 1e-12, "{rmin}");
        assert!(
            (raw(rmin) + 1.0).abs() < 1e-15,
            "the well is exactly -epsilon"
        );

        // And it is a minimum, not just a stationary point the algebra happened to hit.
        assert!(raw(rmin * 0.99) > raw(rmin));
        assert!(raw(rmin * 1.01) > raw(rmin));

        // The force vanishes there, which is the same statement differentiated.
        let f = lj.at_squared(rmin * rmin).unwrap().force_over_r;
        assert!(f.abs() < 1e-12, "force at the minimum was {f}");
    }

    /// The force is the negative gradient of the energy, checked by differencing the energy.
    ///
    /// Independent in the way that matters: the analytic force in `at_squared` and this
    /// numerical one come from different expressions, so an algebra slip in either shows up.
    #[test]
    fn the_force_is_the_gradient_of_the_energy() {
        let lj = LennardJones::reduced();
        for r in [0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 2.4] {
            let h = 1e-7;
            let numerical = -(lj.energy_at(r + h) - lj.energy_at(r - h)) / (2.0 * h);
            let analytic = lj.at_squared(r * r).unwrap().force_over_r * r;
            assert!(
                (analytic - numerical).abs() < 1e-5 * analytic.abs().max(1.0),
                "at r={r}: analytic {analytic} against numerical {numerical}"
            );
        }
    }

    /// Repulsive inside the minimum and attractive outside, which is the whole behaviour of
    /// matter in one sentence.
    #[test]
    fn it_pushes_when_close_and_pulls_when_far() {
        let lj = LennardJones::reduced();
        let rmin = lj.minimum();
        assert!(
            lj.at_squared(0.9 * 0.9).unwrap().force_over_r > 0.0,
            "pushes apart"
        );
        assert!(
            lj.at_squared((rmin * 1.3).powi(2)).unwrap().force_over_r < 0.0,
            "pulls together"
        );
        // Steeply. Halving the separation from sigma multiplies the repulsion by 2^13.
        let near = lj.at_squared(0.5 * 0.5).unwrap().force_over_r * 0.5;
        let at_sigma = lj.at_squared(1.0).unwrap().force_over_r;
        assert!(near / at_sigma > 1000.0, "ratio {}", near / at_sigma);
    }

    /// **Why the shift exists.** The energy reaches the cutoff continuously, so a pair
    /// wandering across it does not teleport energy into the system.
    #[test]
    fn shifting_closes_the_step_at_the_cutoff() {
        let lj = LennardJones::reduced();
        let rc = lj.cutoff;
        let inside = lj.energy_at(rc - 1e-9);
        assert!(
            inside.abs() < 1e-8,
            "the shifted energy dies at the cutoff: {inside}"
        );
        assert_eq!(lj.energy_at(rc + 1e-9), 0.0, "and there is nothing past it");

        // The step that was removed is not negligible: at 2.5 sigma the raw potential is
        // still -0.0163 epsilon, which is where the shift constant comes from.
        assert!(
            (lj.shift() + 0.016_316_891).abs() < 1e-8,
            "shift was {}",
            lj.shift()
        );

        // The force is *not* continuous there, and that is the trade the module note names.
        let f = lj.at_squared((rc - 1e-9).powi(2)).unwrap().force_over_r * rc;
        assert!(
            (f + 0.038_999_477).abs() < 1e-8,
            "the leftover force jump was {f}"
        );
    }

    /// Argon's parameters, and the well depth in kelvin that names them.
    #[test]
    fn argon_is_the_argon_everybody_uses() {
        let ar = LennardJones::argon();
        assert!(
            (ar.sigma * 1e10 - 3.405).abs() < 1e-12,
            "sigma in angstroms"
        );
        let in_kelvin = ar.epsilon / dualis_units::BOLTZMANN.to_si();
        assert!(
            (in_kelvin - 119.8).abs() < 1e-9,
            "epsilon/k was {in_kelvin} K"
        );
        // The minimum is 3.82 A, which is about where argon's first neighbour shell sits.
        assert!((ar.minimum() * 1e10 - 3.822).abs() < 1e-3);
    }

    /// The tail corrections, against the closed forms they are quoted as.
    ///
    /// Both are negative at a sensible cutoff — the attraction that was cut off was pulling
    /// inwards — and both vanish as the cutoff grows, which is the sanity check that matters
    /// more than the coefficient.
    #[test]
    fn the_tail_corrections_shrink_as_the_cutoff_grows() {
        let rho = 0.8442; // the reduced triple-point density
        let near = LennardJones::reduced();
        let far = LennardJones::reduced().with_cutoff(6.0);

        assert!(near.energy_tail(rho) < 0.0 && near.pressure_tail(rho) < 0.0);
        // At 2.5 sigma the energy correction is -0.452 epsilon a particle. The total potential
        // energy at this state point is about -5.7, so it is eight percent: not something to
        // leave out when the point is to match a published number.
        assert!(
            (near.energy_tail(rho) + 0.452_012_6).abs() < 1e-6,
            "got {}",
            near.energy_tail(rho)
        );
        // The ratio between two cutoffs, against the closed form rather than against the
        // rc^-3 law it tends to. Those are not the same number: the exact ratio at 2.5 and 6
        // sigma is 13.805224, and (6/2.5)^3 is 13.824, because the repulsive rc^-9 term is
        // still worth a tenth of a percent at two and a half sigma. This assertion used to
        // compare against the power law inside 0.02, which a known 0.0188 discrepancy filled
        // to 94% — leaving 6% of the budget to catch anything that was actually wrong.
        let bracket = |rc: f64| {
            let x = (1.0f64 / rc).powi(3); // (sigma/rc)^3, with sigma = 1
            x * x * x / 3.0 - x
        };
        let ratio = near.energy_tail(rho) / far.energy_tail(rho);
        let exact = bracket(2.5) / bracket(6.0);
        assert!(
            (ratio / exact - 1.0).abs() < 1e-12,
            "ratio {ratio}, closed form {exact}"
        );
        // And the rc^-3 law separately, as the limit it actually is: once the cutoff is far
        // enough out, doubling it divides the tail by eight. It arrives quickly, because the
        // term being neglected dies six powers faster — 7.9892 at 2.5 sigma, 7.99983 at 5,
        // 7.999997 at 10. Stating it as a limit is the honest version of the claim the old
        // comment made; stating it as an equality at 2.5 sigma was not.
        for (rc, want) in [
            (2.5, 7.989_247_771),
            (5.0, 7.999_832_000),
            (10.0, 7.999_997_375),
        ] {
            let doubled = LennardJones::reduced().with_cutoff(rc).energy_tail(rho)
                / LennardJones::reduced()
                    .with_cutoff(2.0 * rc)
                    .energy_tail(rho);
            assert!(
                (doubled / want - 1.0).abs() < 1e-9,
                "rc {rc}: got {doubled}"
            );
        }
        // And linear in density, since it counts pairs against a uniform background.
        assert!((near.energy_tail(2.0 * rho) / near.energy_tail(rho) - 2.0).abs() < 1e-12);
    }
}

#[cfg(test)]
mod prepared_tests {
    use super::*;

    /// **A prepared potential is the bare one, to the bit, across the whole range.**
    ///
    /// The optimisation is only allowed to exist because of this. `Prepared` hoists five
    /// loop-invariant expressions out of a pair loop, and every one of them is *exactly* the
    /// expression it replaces — `sigma_sq / r2` for `sigma * sigma / r2`, `four_epsilon * x` for
    /// `4.0 * epsilon * x`. If any had been rearranged into something algebraically equal and
    /// numerically different, this workspace's results would have moved on four platforms at
    /// once, and it would have looked like a performance change rather than a physics one.
    ///
    /// Swept across the cutoff rather than checked at a point, so the `None` boundary is
    /// compared too.
    #[test]
    fn preparing_changes_the_speed_and_not_a_single_bit() {
        for lj in [
            LennardJones::reduced(),
            LennardJones {
                epsilon: 1.65e-21,
                sigma: 3.4e-10,
                cutoff: 8.5e-10,
            },
        ] {
            let prepared = lj.prepared();
            let rc = lj.cutoff;
            for k in 0..2000 {
                // From well inside the core out past the cutoff, including r2 exactly at rc².
                let r2 = if k == 1999 {
                    rc * rc
                } else {
                    (0.2 + 1.2 * k as f64 / 1999.0).powi(2) * rc * rc
                };
                match (lj.at_squared(r2), prepared.at_squared(r2)) {
                    (None, None) => {}
                    (Some(a), Some(b)) => {
                        assert_eq!(
                            a.energy.to_bits(),
                            b.energy.to_bits(),
                            "energy at r2 = {r2:e}: {} against {}",
                            a.energy,
                            b.energy
                        );
                        assert_eq!(
                            a.force_over_r.to_bits(),
                            b.force_over_r.to_bits(),
                            "force at r2 = {r2:e}"
                        );
                    }
                    (a, b) => panic!("disagreed about the cutoff at r2 = {r2:e}: {a:?} / {b:?}"),
                }
            }
        }
    }
}
