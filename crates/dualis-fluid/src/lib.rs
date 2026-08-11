#![deny(missing_docs)]

//! Incompressible flow in three dimensions, by projection on a staggered grid.
//!
//! ```text
//!   ∂u/∂t + ∇·(uu) = −∇p/ρ + ν∇²u + g
//!   ∇·u = 0
//! ```
//!
//! # Why this is the hardest of this workspace's domains to trust
//!
//! Every other physics here has closed forms lying around. Fluids has few, its schemes trade
//! stability against numerical diffusion, and **"it looks like a fluid" is the easiest wrong
//! answer in computational physics to accept**: a scheme with the wrong viscosity still makes
//! plausible vortices, and a scheme that quietly loses momentum still makes a pretty picture.
//!
//! So this crate is built around the three exact solutions that exist, and each one is chosen to
//! be blind to a different mistake:
//!
//! ```text
//!   Poiseuille    u(y) = (g/2ν)·y(h−y)        exact — a quadratic, and a second difference of one
//!   Couette       u(y) = U·y/h                 exact — linear, and blind to advection entirely
//!   Taylor–Green  e^{−2νk²t}                   the full nonlinear equations, decay rate and all
//! ```
//!
//! The first two are unidirectional and steady, so the advection term is identically zero in
//! both — they cannot check it at all, and saying so is the point. Taylor–Green can: it is an
//! exact solution of the *complete* equations, in which the nonlinear term is balanced by the
//! pressure gradient rather than absent. Beside it sit two statements that hold at machine
//! precision and catch what a decay rate is too coarse to see: a uniform flow must remain exactly
//! uniform, and total momentum in a periodic box must not move at all.
//!
//! # The staggered grid, and its one guarantee
//!
//! Velocities on cell **faces**, pressure at cell **centres** — the same arrangement Yee uses for
//! electromagnetism, and for the same reason: the divergence of a face velocity lands naturally at
//! a cell centre, and the gradient of a centred pressure lands naturally on a face. No
//! interpolation, and no checkerboard pressure mode.
//!
//! After the projection, `∇·u` is **the residual of the pressure solve** and nothing else — see
//! [`Channel::divergence`]. That is weaker than electromagnetism's identity, which holds exactly:
//! here it holds to whatever the conjugate-gradient solve was asked for.
//!
//! # Two limits, and one of them is on the grid rather than on the step
//!
//! ```text
//!   dt ≤ dx²/(6ν)      viscous, the same Fourier limit conduction has
//!   dt ≤ dx/|u|max     advective, the Courant limit
//!   |u|dx/ν ≤ 2        the cell Reynolds number — a property of the *mesh*
//! ```
//!
//! The third is the one that surprises people. Central differences on the advection term go
//! unstable when a cell is too coarse for the viscosity to smooth what advection sharpens, and no
//! amount of shortening the step fixes it: the mesh is wrong. [`Channel::cell_reynolds`] reports
//! it and [`Channel::step`](dualis_core::Domain::step) refuses above two, rather than producing
//! the wiggles that a reader would take for turbulence.
//!
//! # What is deliberately not here
//!
//! No turbulence model, no compressibility, no free surface, no immersed geometry, no adaptive
//! mesh. A box with periodic sides and optional walls, which is exactly what the three exact
//! solutions live in.

use dualis_units::{Density, Diffusivity};

mod channel;

pub use channel::{Channel, Walls, CELL_REYNOLDS_LIMIT};

/// A Newtonian fluid.
#[derive(Clone, Copy, Debug)]
pub struct Fluid {
    /// Density.
    pub density: Density,
    /// **Kinematic** viscosity, `μ/ρ`, in m²/s.
    ///
    /// The one that appears in the equations above and in every closed form below. Tables quote
    /// the dynamic viscosity as often as this one and the two differ by a factor of a thousand for
    /// water, which is the sort of error dimensions cannot catch and a name can.
    pub kinematic_viscosity: Diffusivity,
}

impl Fluid {
    /// Water at 20 °C.
    pub fn water() -> Fluid {
        Fluid {
            density: Density::from_si(998.2),
            kinematic_viscosity: Diffusivity::from_si(1.004e-6),
        }
    }

    /// Air at 20 °C and one atmosphere.
    pub fn air() -> Fluid {
        Fluid {
            density: Density::from_si(1.204),
            kinematic_viscosity: Diffusivity::from_si(1.511e-5),
        }
    }

    /// A fluid with whatever properties, for a test that wants a convenient Reynolds number.
    pub fn new(density: Density, kinematic_viscosity: Diffusivity) -> Fluid {
        Fluid {
            density,
            kinematic_viscosity,
        }
    }

    /// Dynamic viscosity, `ρν`, in Pa·s.
    pub fn dynamic_viscosity(&self) -> f64 {
        self.density.to_si() * self.kinematic_viscosity.to_si()
    }
}

/// The mean speed of plane Poiseuille flow driven by a body force.
///
/// ```text
///   u(y) = (g/2ν)·y(h−y)      ⇒      ū = g h² / (12 ν)
/// ```
///
/// A closed form of the gap, the force and the viscosity, and of nothing else. The profile is a
/// **quadratic**, which is why a second-order scheme reproduces it exactly rather than nearly: the
/// second difference of a quadratic is that quadratic's second derivative, with no truncation
/// error at all.
pub fn poiseuille_mean_speed(force_per_mass: f64, gap: f64, kinematic_viscosity: f64) -> f64 {
    force_per_mass * gap * gap / (12.0 * kinematic_viscosity)
}

/// The rate at which a Taylor–Green vortex's velocity decays, `2νk²`.
///
/// The kinetic energy decays at twice this, because energy goes as the square. Confusing the two
/// is a factor of two in the viscosity, and it is the reason both are named here rather than one.
pub fn taylor_green_rate(wavenumber: f64, kinematic_viscosity: f64) -> f64 {
    2.0 * kinematic_viscosity * wavenumber * wavenumber
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The two viscosities differ by the density**, and both are named so neither is guessed.
    #[test]
    fn the_kinematic_and_dynamic_viscosities_are_a_density_apart() {
        let w = Fluid::water();
        assert!(
            (w.dynamic_viscosity() / 1.002e-3 - 1.0).abs() < 0.01,
            "water is about 1.0 mPa.s: {:.4e}",
            w.dynamic_viscosity()
        );
        // Air is fifteen times *more* viscous kinematically than water and fifty times less
        // dynamically, which is the whole reason the distinction is worth a name.
        let a = Fluid::air();
        assert!(
            a.kinematic_viscosity.to_si() > 10.0 * w.kinematic_viscosity.to_si()
                && a.dynamic_viscosity() < 0.1 * w.dynamic_viscosity(),
            "air: nu {:.3e} against water's {:.3e}, mu {:.3e} against {:.3e}",
            a.kinematic_viscosity.to_si(),
            w.kinematic_viscosity.to_si(),
            a.dynamic_viscosity(),
            w.dynamic_viscosity()
        );
    }

    /// **The Poiseuille mean is the profile's own average**, which is an integral and not a fit.
    #[test]
    fn the_poiseuille_mean_is_the_integral_of_its_profile() {
        let (g, h, nu) = (0.5, 0.02, 1e-5);
        let closed = poiseuille_mean_speed(g, h, nu);
        // Numerically integrate `(g/2nu) y (h-y)` over the gap, with enough points that the
        // trapezium rule's own error is below the comparison.
        let n = 100_000;
        let mut sum = 0.0;
        for i in 0..=n {
            let y = h * i as f64 / n as f64;
            let w = if i == 0 || i == n { 0.5 } else { 1.0 };
            sum += w * (g / (2.0 * nu)) * y * (h - y);
        }
        let mean = sum * (h / n as f64) / h;
        assert!(
            (mean / closed - 1.0).abs() < 1e-9,
            "gh^2/12nu is the mean of the parabola: {mean:.6e} against {closed:.6e}"
        );
    }

    /// **Energy decays at twice the velocity's rate**, which is the factor of two this crate names.
    #[test]
    fn the_energy_rate_is_twice_the_velocity_rate() {
        let (k, nu) = (100.0, 1e-4);
        let rate = taylor_green_rate(k, nu);
        assert!((rate - 2.0 * nu * k * k).abs() < 1e-12);
        // A velocity going as `e^{-rt}` gives an energy going as `e^{-2rt}`.
        let t = 0.37;
        let u = (-rate * t).exp();
        assert!(
            ((u * u) / (-2.0 * rate * t).exp() - 1.0).abs() < 1e-12,
            "energy is the square of the velocity"
        );
    }
}
