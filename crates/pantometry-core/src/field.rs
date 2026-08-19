//! Quantities that vary over space and time.
//!
//! Continuum physics is fields most of the way down: a temperature field, a
//! velocity field, a stress field, an electromagnetic field. They differ in what
//! they hold and in the equation that governs them, but the three operators used
//! to write those equations — gradient, divergence, Laplacian — are the same
//! three every time. Implementing them once is most of what a kernel owes a
//! continuum domain.
//!
//! # Values are raw SI numbers, and that is a real limitation
//!
//! A [`ScalarField`] returns `f64` and a [`VectorField`] returns `DVec3`, both in
//! SI base units, rather than the dimensioned types from `pantometry-units`. The
//! dimension varies per field — kelvin here, pascals there — and expressing "the
//! gradient of this field has the field's dimension divided by a length" needs
//! arithmetic on const generic parameters, which is unstable. So the dimension
//! lives on the concrete type that implements the trait, and is documented rather
//! than checked, at exactly this one boundary.
//!
//! # The derivatives are finite differences, and the step matters
//!
//! Central differences, second-order accurate, with the step passed in. There is
//! no good default for it: too large and the truncation error dominates, too small
//! and cancellation in the subtraction does. For a field with curvature scale `Lc`
//! and values near 1, `h ≈ Lc * 1e-5` is a reasonable start, and a field that
//! knows its own analytic derivative should override these methods and skip the
//! question entirely.

use glam::DVec3;
use pantometry_units::{Length, LengthVec, Time};

/// A scalar field: temperature, pressure, concentration, potential.
pub trait ScalarField {
    /// The value at a place and time, in this field's SI base unit.
    fn at(&self, p: LengthVec, t: Time) -> f64;

    /// What this field is measured in — `"Pa"`, `"C"`, `"V"`.
    ///
    /// The fifth thing a layer above had to know a domain by name to find out. A legend, an axis
    /// and a CSV header all need it, and none of them should have to match on domain types to
    /// get a two-character string.
    ///
    /// Defaults to `""`, which renders as a quantity with no unit — visibly odd rather than
    /// quietly wrong, which is the best a default can do here.
    fn unit(&self) -> &'static str {
        ""
    }

    /// ∇f, by central differences over `h`. Units are the field's per metre.
    fn gradient(&self, p: LengthVec, t: Time, h: Length) -> DVec3 {
        let step = h.to_si();
        let mut out = DVec3::ZERO;
        for axis in 0..3 {
            let mut d = DVec3::ZERO;
            d[axis] = step;
            let plus = self.at(p + LengthVec::from_si(d), t);
            let minus = self.at(p - LengthVec::from_si(d), t);
            out[axis] = (plus - minus) / (2.0 * step);
        }
        out
    }

    /// ∇²f — the operator that makes diffusion diffuse. A point hotter than the
    /// average of its neighbours has a negative Laplacian and cools; that is the
    /// whole content of the heat equation.
    fn laplacian(&self, p: LengthVec, t: Time, h: Length) -> f64 {
        let step = h.to_si();
        let centre = self.at(p, t);
        let mut sum = 0.0;
        for axis in 0..3 {
            let mut d = DVec3::ZERO;
            d[axis] = step;
            sum += self.at(p + LengthVec::from_si(d), t) + self.at(p - LengthVec::from_si(d), t)
                - 2.0 * centre;
        }
        sum / (step * step)
    }

    /// ∂f/∂t, by a central difference in time.
    fn rate(&self, p: LengthVec, t: Time, dt: Time) -> f64 {
        (self.at(p, t + dt) - self.at(p, t - dt)) / (2.0 * dt.to_si())
    }
}

/// A vector field: velocity, force per volume, electric field, heat flux.
pub trait VectorField {
    /// The value at a place and time, in this field's SI base unit.
    fn at(&self, p: LengthVec, t: Time) -> DVec3;

    /// ∇·F — net outflow per unit volume. Zero everywhere means nothing is being
    /// created or destroyed, which is how a conservation law is written locally.
    fn divergence(&self, p: LengthVec, t: Time, h: Length) -> f64 {
        let step = h.to_si();
        let mut sum = 0.0;
        for axis in 0..3 {
            let mut d = DVec3::ZERO;
            d[axis] = step;
            sum += (self.at(p + LengthVec::from_si(d), t)[axis]
                - self.at(p - LengthVec::from_si(d), t)[axis])
                / (2.0 * step);
        }
        sum
    }

    /// ∇×F — local rotation. Zero everywhere means the field is a gradient of
    /// something, which is why an electrostatic field has a potential and a
    /// magnetic one does not.
    fn curl(&self, p: LengthVec, t: Time, h: Length) -> DVec3 {
        let step = h.to_si();
        let d = |axis: usize| {
            let mut e = DVec3::ZERO;
            e[axis] = step;
            let plus = self.at(p + LengthVec::from_si(e), t);
            let minus = self.at(p - LengthVec::from_si(e), t);
            (plus - minus) / (2.0 * step)
        };
        let (dx, dy, dz) = (d(0), d(1), d(2));
        DVec3::new(dy.z - dz.y, dz.x - dx.z, dx.y - dy.x)
    }
}

/// The same value everywhere and always. Useful as an ambient condition, and as
/// the thing a finite-difference test should return zero derivatives for.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Uniform(pub f64);

impl ScalarField for Uniform {
    fn at(&self, _p: LengthVec, _t: Time) -> f64 {
        self.0
    }
}

/// A field from a closure, for the common case where the physics is a formula.
pub struct Analytic<F>(pub F);

impl<F: Fn(LengthVec, Time) -> f64> ScalarField for Analytic<F> {
    fn at(&self, p: LengthVec, t: Time) -> f64 {
        (self.0)(p, t)
    }
}

impl<F: Fn(LengthVec, Time) -> DVec3> VectorField for Analytic<F> {
    fn at(&self, p: LengthVec, t: Time) -> DVec3 {
        (self.0)(p, t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: Length = Length::from_si(1e-4);

    /// Checked against closed forms, not against another discretisation: for
    /// f = x² + 2y² + 3z², ∇f = (2x, 4y, 6z) and ∇²f = 2 + 4 + 6 = 12, exactly,
    /// everywhere.
    #[test]
    fn derivatives_match_the_closed_form() {
        let f = Analytic(|p: LengthVec, _t: Time| {
            let v = p.to_si();
            v.x * v.x + 2.0 * v.y * v.y + 3.0 * v.z * v.z
        });
        let p = LengthVec::m(0.3, -0.2, 0.5);
        let t = Time::ZERO;

        let grad = f.gradient(p, t, H);
        let expected = DVec3::new(0.6, -0.8, 3.0);
        assert!((grad - expected).length() < 1e-9, "got {grad}");

        // A quadratic's second difference is exact, so this is a tight check.
        let lap = f.laplacian(p, t, H);
        assert!((lap - 12.0).abs() < 1e-6, "got {lap}");
    }

    /// Divergence and curl, against a field whose answers are known: a rigid
    /// rotation about z has zero divergence and a curl of exactly 2ω ẑ.
    #[test]
    fn divergence_and_curl_match_the_closed_form() {
        let omega = 3.0;
        let rotation = Analytic(|p: LengthVec, _t: Time| {
            let v = p.to_si();
            DVec3::new(-3.0 * v.y, 3.0 * v.x, 0.0)
        });
        let p = LengthVec::m(0.4, 0.1, -0.2);
        let t = Time::ZERO;

        assert!(
            VectorField::divergence(&rotation, p, t, H).abs() < 1e-9,
            "a rotation moves fluid around, not outwards"
        );
        let curl = rotation.curl(p, t, H);
        assert!(
            (curl - DVec3::new(0.0, 0.0, 2.0 * omega)).length() < 1e-9,
            "got {curl}"
        );
    }

    /// A radial outflow has divergence 3 and no curl — the complement of the
    /// previous test, and the pair distinguishes a sign error from a correct
    /// implementation.
    #[test]
    fn a_radial_field_diverges_without_rotating() {
        let outflow = Analytic(|p: LengthVec, _t: Time| p.to_si());
        let p = LengthVec::m(0.2, -0.4, 0.7);
        let t = Time::ZERO;
        assert!((VectorField::divergence(&outflow, p, t, H) - 3.0).abs() < 1e-9);
        assert!(outflow.curl(p, t, H).length() < 1e-9);
    }

    /// Time derivatives too, and a field that does not move has none.
    #[test]
    fn a_uniform_field_has_no_derivatives_of_any_kind() {
        let u = Uniform(4.2);
        let p = LengthVec::m(1.0, 2.0, 3.0);
        assert_eq!(u.at(p, Time::s(9.0)), 4.2);
        assert_eq!(u.gradient(p, Time::ZERO, H), DVec3::ZERO);
        assert_eq!(u.laplacian(p, Time::ZERO, H), 0.0);
        assert_eq!(u.rate(p, Time::ZERO, Time::s(0.1)), 0.0);
    }

    #[test]
    fn time_derivatives_match_the_closed_form() {
        // f = t², so df/dt = 2t.
        let f = Analytic(|_p: LengthVec, t: Time| t.to_si() * t.to_si());
        let rate = f.rate(LengthVec::ZERO, Time::s(3.0), Time::s(1e-4));
        assert!((rate - 6.0).abs() < 1e-9, "got {rate}");
    }
}
