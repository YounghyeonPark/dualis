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
//! plasticity, no large rotation, no contact, no fracture.
//!
//! **Dynamics is here now**, as [`Waves`] rather than as a mode on [`Block`]: `ρü = ∇·σ` on the same
//! element, marched with central differences. `Block` stays [`Kind::QuasiStatic`](dualis_core::Kind)
//! and is right about itself — a body with velocity in it has a different lifecycle, a different
//! stability limit and a different thing to conserve. What the two share is the element, which is the
//! part worth sharing: the static tests check that operator against four exact moduli and the
//! dynamic ones check it against two exact speeds.

use dualis_core::Substance;
use dualis_units::{Density, Pressure, Velocity};

mod block;
mod element;
mod waves;

pub use block::{Block, Face};
pub use waves::{Axis, Waves};

/// A linear elastic material.
#[derive(Clone, Copy, Debug)]
pub struct Elastic {
    /// Young's modulus.
    pub youngs_modulus: Pressure,
    /// Poisson's ratio. Between −1 and ½ for a stable isotropic solid.
    pub poisson_ratio: f64,
    /// Density. Unused by the static solve — it sets the **wave speeds**, which is the only place
    /// a mass enters a problem that has no inertia in it.
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

    /// A material from a [`Substance`], or `None` if that substance has no mechanical description.
    ///
    /// The kernel's catalogue carries `youngs_modulus`, `poisson_ratio` and `density`, which is
    /// exactly what this type needs — so the conversion is arithmetic-free and the only question is
    /// which crate it lives in. **Here**, because the kernel must not know that elasticity exists;
    /// a `Substance::elastic()` would have made it know.
    ///
    /// `None` for a fluid: [`Substance::water`](dualis_core::Substance::water) has no `mechanical`
    /// entry at all, because a liquid has no shear modulus to make one out of. That is the same
    /// asymmetry [`s_wave_speed`](Elastic::s_wave_speed) is about.
    ///
    /// # It drops the yield strength, and that is the one thing to know
    ///
    /// `Substance` says where a material **stops coming back**; this type has no yield and no
    /// plasticity, so it cannot represent that and does not pretend to. A solve past yield returns a
    /// displacement that is arithmetically correct and physically meaningless, and nothing in the
    /// answer says which.
    ///
    /// So keep the `Substance`. The strain at which the linear model stops applying is
    /// `yield_strength / E`, and across this catalogue it spans **130×**:
    ///
    /// ```text
    ///   ice        0.011%     brittle: 1 MPa against 9.1 GPa, and no plastic region beyond it
    ///   Cu ETP     0.060%
    ///   N-BK7      0.073%
    ///   Al 6061    0.401%
    ///   PLA        1.429%     a polymer has twenty times a metal's elastic room
    /// ```
    ///
    /// So "small strain" is not one number, and a bound that would be absurdly conservative for PLA
    /// is already past yield for ice. A load case that looked reasonable can be past several of these.
    pub fn from_substance(substance: &Substance) -> Option<Elastic> {
        let m = substance.mechanical?;
        Elastic::new(m.youngs_modulus, m.poisson_ratio, substance.density)
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

    /// The **pressure** wave speed, `√((λ+2μ)/ρ)` — the fast one, and the first arrival.
    ///
    /// The constrained modulus over the density, so it is the speed of a compression that cannot
    /// bulge sideways: a wave in the bulk, where the material a wavelength away is holding the sides
    /// still. It is **not** `√(E/ρ)`, which is the speed along a thin rod free to bulge, and for
    /// aluminium those are 6149 and 5051 m/s — a 22% difference from the same two constants.
    pub fn p_wave_speed(&self) -> Velocity {
        Velocity::from_si((self.constrained_modulus().to_si() / self.density.to_si()).sqrt())
    }

    /// The **shear** wave speed, `√(μ/ρ)` — the slow one, and the one a fluid does not have.
    ///
    /// A shear wave needs something to resist a change of shape at constant volume, which is what
    /// `μ` is. A fluid has none, so it carries no shear wave at all — which is why
    /// [`AcousticProps`](dualis_core::substance::AcousticProps) has room for one speed and a solid
    /// has two.
    pub fn s_wave_speed(&self) -> Velocity {
        Velocity::from_si((self.shear_modulus().to_si() / self.density.to_si()).sqrt())
    }

    /// `c_p / c_s = √(2(1−ν)/(1−2ν))` — a function of Poisson's ratio and **nothing else**.
    ///
    /// Both `E` and `ρ` cancel, which makes this the sharpest check a wave solver can be given: a
    /// scheme with the wrong stiffness or the wrong mass still has to land on it, and only a scheme
    /// with the wrong *operator* fails. It runs from `√2` at `ν = 0` to infinity as `ν → ½`, where
    /// the material becomes incompressible and the compression wave has nothing left to compress.
    pub fn speed_ratio(&self) -> f64 {
        let v = self.poisson_ratio;
        (2.0 * (1.0 - v) / (1.0 - 2.0 * v)).sqrt()
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
