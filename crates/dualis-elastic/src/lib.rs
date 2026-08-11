#![deny(missing_docs)]

//! Linear elasticity as a field: what a shape does under load, solved rather than stated.
//!
//! `dualis-mechanics` moves bodies that do not change shape. This is the other half — a body that
//! does, and does not move:
//!
//! ```text
//!   ∇·σ = 0,   σ = λ tr(ε) I + 2μ ε,   ε = ½(∇u + ∇uᵀ)
//! ```
//!
//! Elliptic, like [`Conductor`](https://docs.rs/dualis-electrical) and
//! [`Puck`](https://docs.rs/dualis-porous), and solved the same way — conjugate gradients on a
//! symmetric positive-definite system — except that the unknown is a **vector** at every node
//! rather than a scalar. That is the whole of what is new here, and it is why the operator is
//! assembled from an energy rather than differenced, which the element module's own docs
//! set out.
//!
//! # What comes out exactly, and what only converges
//!
//! Trilinear elements reproduce any linear displacement field exactly, so every case where the
//! strain is uniform is exact at any mesh size. That is four independent moduli, each a different
//! combination of `λ` and `μ`, and getting all four right is a much stronger statement than
//! getting one:
//!
//! ```text
//!   uniaxial stress    σ/ε  =  E                          sides free
//!   uniaxial strain    σ/ε  =  E(1−ν) / ((1+ν)(1−2ν))     sides held
//!   hydrostatic        p/ΔV/V = K = E / (3(1−2ν))         all six faces pressed
//!   simple shear       τ/γ  =  G = E / (2(1+ν))
//! ```
//!
//! A body that had `λ` and `μ` transposed reproduces none of them. A body that had the shear rows
//! of `D` carrying `2μ` instead of `μ` passes the first three and fails the fourth.
//!
//! **Bending only converges**, and from the stiff side. A fully integrated trilinear element
//! develops shear strain where it should be flexing, so a cantilever comes out too stiff and
//! approaches `PL³/3EI` from below as the mesh refines. That is stated rather than tuned away.
//!
//! # The energy identity, which is this domain's Tellegen
//!
//! At equilibrium the strain energy is half the work the boundary loads did:
//!
//! ```text
//!   2U = Σ f·u
//! ```
//!
//! Clapeyron's theorem. Both sides are computed independently — one from the stiffness and the
//! displacement, the other from the loads and the displacement — so the two agreeing is a check
//! and not a tautology. It is the same shape of statement `Conductor` makes when its field power
//! matches its terminal power, and it is the sharpest single thing that can be said about whether
//! a discretisation is self-consistent.
//!
//! # What is deliberately not here
//!
//! Small strain and linear material, which is the regime the closed forms above live in. No
//! plasticity, no large rotation, no contact, no fracture, no dynamics — this is
//! [`Kind::QuasiStatic`](dualis_core::Kind), a solve rather than a march.

use dualis_units::{Density, Pressure};

mod block;
mod element;

pub use block::{Block, Face};

/// A linear elastic material.
#[derive(Clone, Copy, Debug)]
pub struct Elastic {
    /// Young's modulus.
    pub youngs_modulus: Pressure,
    /// Poisson's ratio. Between −1 and ½ for a stable isotropic solid.
    pub poisson_ratio: f64,
    /// Density. Unused by the static solve and carried so a mass can be stated once.
    pub density: Density,
}

impl Elastic {
    /// A material from the pair an engineer quotes.
    ///
    /// # Both bounds are refused rather than clamped
    ///
    /// `ν = ½` is incompressible and `λ` is infinite there; `ν = −1` is the other end of
    /// thermodynamic stability and `G` diverges. Neither is a material this solver can represent,
    /// and clamping to "nearly" either one produces a stiffness that is finite, plausible and
    /// enormous — which is worse than a refusal, because it runs.
    ///
    /// The practical limit arrives before the mathematical one: at `ν = 0.49` a fully integrated
    /// trilinear element **volumetrically locks**, coming out many times too stiff, and no warning
    /// here can tell that from a genuinely stiff answer. Rubber wants a different element.
    pub fn new(youngs_modulus: Pressure, poisson_ratio: f64, density: Density) -> Option<Elastic> {
        if !(-0.999..0.4999).contains(&poisson_ratio) || youngs_modulus.to_si() <= 0.0 {
            return None;
        }
        Some(Elastic {
            youngs_modulus,
            poisson_ratio,
            density,
        })
    }

    /// Aluminium 6061-T6.
    pub fn aluminium_6061() -> Elastic {
        Elastic::new(Pressure::from_si(68.9e9), 0.33, Density::from_si(2700.0))
            .expect("6061 is a representable material")
    }

    /// Structural steel.
    pub fn steel() -> Elastic {
        Elastic::new(Pressure::from_si(200.0e9), 0.30, Density::from_si(7850.0))
            .expect("steel is a representable material")
    }

    /// The shear modulus, `E / (2(1+ν))`.
    pub fn shear_modulus(&self) -> Pressure {
        Pressure::from_si(self.youngs_modulus.to_si() / (2.0 * (1.0 + self.poisson_ratio)))
    }

    /// The bulk modulus, `E / (3(1−2ν))`.
    pub fn bulk_modulus(&self) -> Pressure {
        Pressure::from_si(self.youngs_modulus.to_si() / (3.0 * (1.0 - 2.0 * self.poisson_ratio)))
    }

    /// The **constrained** modulus, `E(1−ν) / ((1+ν)(1−2ν))`, which is `λ + 2μ`.
    ///
    /// What a block resists compression with when its sides cannot move — a puck in a basket, a
    /// core in a bore, soil under a footing. It is larger than `E`, by 1.5× at `ν = 0.3`, and
    /// mistaking one for the other is the commonest way to get a confined stiffness wrong.
    pub fn constrained_modulus(&self) -> Pressure {
        let (lambda, mu) = element::lame(self.youngs_modulus.to_si(), self.poisson_ratio);
        Pressure::from_si(lambda + 2.0 * mu)
    }

    pub(crate) fn lame(&self) -> (f64, f64) {
        element::lame(self.youngs_modulus.to_si(), self.poisson_ratio)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The four moduli agree with each other**, which is a closed-form identity and not a
    /// property of any solver.
    ///
    /// `E`, `G`, `K` and `M` are four combinations of two constants, so two identities hold among
    /// them for every material. Checking those here means the solver tests below are checking the
    /// discretisation rather than the arithmetic in this file.
    #[test]
    fn the_moduli_are_consistent_with_each_other() {
        for m in [Elastic::aluminium_6061(), Elastic::steel()] {
            let (e, nu) = (m.youngs_modulus.to_si(), m.poisson_ratio);
            let (g, k, mm) = (
                m.shear_modulus().to_si(),
                m.bulk_modulus().to_si(),
                m.constrained_modulus().to_si(),
            );
            // E = 9KG/(3K+G), the standard identity.
            let from_kg = 9.0 * k * g / (3.0 * k + g);
            assert!(
                (from_kg / e - 1.0).abs() < 1e-12,
                "E from K and G: {from_kg:.6e} against {e:.6e}"
            );
            // M = K + 4G/3.
            assert!(
                ((k + 4.0 * g / 3.0) / mm - 1.0).abs() < 1e-12,
                "M = K + 4G/3: {:.6e} against {mm:.6e}",
                k + 4.0 * g / 3.0
            );
            // And M is meaningfully above E, or the distinction this crate makes is not one.
            assert!(
                mm / e > 1.2,
                "the constrained modulus should be well above E at nu={nu}: {:.3}x",
                mm / e
            );
        }
    }

    /// **A material outside the stable range is refused rather than clamped.**
    #[test]
    fn an_impossible_poisson_ratio_is_refused() {
        assert!(Elastic::new(Pressure::from_si(1e9), 0.5, Density::from_si(1.0)).is_none());
        assert!(Elastic::new(Pressure::from_si(1e9), 0.6, Density::from_si(1.0)).is_none());
        assert!(Elastic::new(Pressure::from_si(1e9), -1.5, Density::from_si(1.0)).is_none());
        assert!(Elastic::new(Pressure::from_si(0.0), 0.3, Density::from_si(1.0)).is_none());
        // And a negative ratio, which is unusual but real, is allowed.
        assert!(Elastic::new(Pressure::from_si(1e9), -0.2, Density::from_si(1.0)).is_some());
    }
}
