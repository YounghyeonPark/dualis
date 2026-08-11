#![deny(missing_docs)]

//! Electromagnetic fields in three dimensions, on the grid Yee built for them.
//!
//! ```text
//!   ∂D/∂t =  ∇×H − J        ∇·D = ρ
//!   ∂B/∂t = −∇×E            ∇·B = 0
//! ```
//!
//! Four equations, and the interesting thing about them is that the two on the right are not
//! independent: take the divergence of either curl equation and the divergence constraint follows
//! from the one beside it, because `∇·(∇×F) = 0` for any `F`. So a scheme that gets the curls right
//! gets the divergences free — *if* its discrete operators satisfy the same identity.
//!
//! # The Yee grid, and the reason it is the grid
//!
//! `E` lives on cell **edges** and `H` on cell **faces**, each component offset half a cell from
//! the others:
//!
//! ```text
//!   Ex at (i+½, j,   k  )      Hx at (i,   j+½, k+½)
//!   Ey at (i,   j+½, k  )      Hy at (i+½, j,   k+½)
//!   Ez at (i,   j,   k+½)      Hz at (i+½, j+½, k  )
//! ```
//!
//! Every curl is then a difference of quantities that already sit where the curl belongs, with no
//! interpolation anywhere. And the discrete divergence of the discrete curl is **identically zero**
//! — every term appears twice with opposite signs and cancels in exact arithmetic, not to within a
//! tolerance.
//!
//! That is why [`Cavity::magnetic_divergence`] is `0` and stays `0` after ten thousand steps rather
//! than drifting, and it is the strongest single statement that can be made about a
//! discretisation: not that it conserves something to `1e-12`, but that the conservation is an
//! algebraic identity of the update itself.
//!
//! A collocated grid has no such identity. It runs, it looks like electromagnetism, and it
//! accumulates magnetic charge.
//!
//! # What is checkable here
//!
//! ```text
//!   ∇·B = 0                exactly, and it is an identity rather than a tolerance
//!   f_mnp = (c/2)√(…)      the cavity's resonances, approached at second order in dx
//!   ½(εE² + μH²)           constant in a lossless box with conducting walls
//!   dt ≤ dx/(c√3)          the Courant limit, refused past rather than run past
//!   ⟨U_E⟩ = ⟨U_H⟩          equipartition in a standing mode
//! ```
//!
//! # What is deliberately not here
//!
//! No dispersive or anisotropic media, no nonlinearity, no moving charges.
//!
//! The boundary is either a perfect conductor or [`Boundary::Open`], which is Mur's first-order
//! condition — exact for a wave arriving along the normal and progressively worse away from it. It
//! is **not** a perfectly matched layer. A PML would be right to four digits at every angle and is
//! a larger piece of machinery than this crate has earned; what is here is measured instead of
//! claimed: a line source in an open 24³ box is down to **0.149%** of its energy after two and a
//! half crossings, where a conducting box still has all of it.

use dualis_units::{Frequency, Velocity};

mod cavity;

pub use cavity::{Boundary, Cavity, Wall, COURANT_3D};

/// Permittivity of free space, F/m.
pub const EPSILON_0: f64 = 8.854_187_812_8e-12;
/// Permeability of free space, H/m.
pub const MU_0: f64 = 1.256_637_062_12e-6;

/// A linear, isotropic, non-dispersive medium.
#[derive(Clone, Copy, Debug)]
pub struct Medium {
    /// Relative permittivity.
    pub relative_permittivity: f64,
    /// Relative permeability.
    pub relative_permeability: f64,
    /// Conductivity, S/m. Zero is lossless, which is the only case with an exact energy statement.
    pub conductivity: f64,
}

impl Medium {
    /// Vacuum.
    pub fn vacuum() -> Medium {
        Medium {
            relative_permittivity: 1.0,
            relative_permeability: 1.0,
            conductivity: 0.0,
        }
    }

    /// A lossless dielectric.
    pub fn dielectric(relative_permittivity: f64) -> Medium {
        Medium {
            relative_permittivity: relative_permittivity.max(1e-12),
            ..Medium::vacuum()
        }
    }

    /// Permittivity, F/m.
    pub fn permittivity(&self) -> f64 {
        EPSILON_0 * self.relative_permittivity
    }

    /// Permeability, H/m.
    pub fn permeability(&self) -> f64 {
        MU_0 * self.relative_permeability
    }

    /// The speed of light in it, `1/√(εμ)`.
    pub fn wave_speed(&self) -> Velocity {
        Velocity::from_si(1.0 / (self.permittivity() * self.permeability()).sqrt())
    }

    /// Whether anything in it dissipates.
    pub fn is_lossless(&self) -> bool {
        self.conductivity <= 0.0
    }
}

/// The resonant frequency of a rectangular cavity's `(m, n, p)` mode.
///
/// ```text
///   f = (c/2) √( (m/a)² + (n/b)² + (p/d)² )
/// ```
///
/// A property of the box and the medium and nothing else. At least two of the three indices must
/// be non-zero for a mode to exist at all, which this does not police — a caller asking for
/// `(1,0,0)` gets the frequency of something the cavity cannot hold.
pub fn cavity_frequency(size: [f64; 3], mode: [u32; 3], medium: &Medium) -> Frequency {
    let s: f64 = (0..3)
        .map(|a| (mode[a] as f64 / size[a]).powi(2))
        .sum::<f64>()
        .sqrt();
    Frequency::from_si(0.5 * medium.wave_speed().to_si() * s)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The speed of light comes out of the two constants**, rather than being one of them.
    #[test]
    fn the_wave_speed_is_one_over_root_epsilon_mu() {
        let c = Medium::vacuum().wave_speed().to_si();
        assert!(
            (c / 299_792_458.0 - 1.0).abs() < 1e-9,
            "1/sqrt(eps0 mu0) is c: {c:.6e}"
        );
        // And a dielectric slows it by exactly the index.
        let glass = Medium::dielectric(2.25);
        assert!(
            (glass.wave_speed().to_si() * 1.5 / c - 1.0).abs() < 1e-12,
            "n = sqrt(eps_r): {:.9}",
            c / glass.wave_speed().to_si()
        );
    }

    /// **A cube's degenerate modes are degenerate**, which the formula gives and nothing else has
    /// to arrange.
    #[test]
    fn a_cubes_modes_come_in_the_degeneracies_the_symmetry_demands() {
        let cube = [0.1, 0.1, 0.1];
        let v = Medium::vacuum();
        let a = cavity_frequency(cube, [1, 1, 0], &v).to_si();
        let b = cavity_frequency(cube, [1, 0, 1], &v).to_si();
        let c = cavity_frequency(cube, [0, 1, 1], &v).to_si();
        assert!(
            (a - b).abs() < 1e-6 && (b - c).abs() < 1e-6,
            "the three (1,1,0) permutations are one frequency in a cube: {a} {b} {c}"
        );
        // And the closed form, from the geometry alone: c√2/(2a).
        let closed = 299_792_458.0 * 2.0f64.sqrt() / (2.0 * 0.1);
        assert!(
            (a / closed - 1.0).abs() < 1e-9,
            "(1,1,0) in a 100 mm cube is c√2/2a: {a:.6e} against {closed:.6e}"
        );
    }
}
