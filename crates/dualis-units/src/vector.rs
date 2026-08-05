//! Vector quantities: three components sharing one dimension.
//!
//! A position, a velocity, a force and a field are all `DVec3` to the compiler,
//! and adding two of them is a bug the compiler cannot see. [`QVec3`] carries the
//! same seven exponents as [`Qty`], so a displacement and a velocity stop being
//! interchangeable.
//!
//! Two of the operations are worth noticing, because they fall out of the
//! dimensions rather than being decided:
//!
//! - [`QVec3::normalize`] returns a bare `DVec3`. A direction has no dimension —
//!   dividing a length by a length leaves a pure number — so a unit vector is
//!   exactly the right type for "which way", and a ray direction cannot be
//!   mistaken for a displacement.
//! - [`QVec3::length`] returns the scalar of the *same* dimension, which needs no
//!   exponent arithmetic and so works for every dimension at once.
//!
//! `dot` and `cross` are missing on purpose: both change the dimension, and there
//! is no way to express "the square of L" in a const generic parameter without
//! unstable features. [`QVec3::along`] covers the case that actually comes up —
//! projecting onto a unit direction, which preserves the dimension.

use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use glam::DVec3;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::{Force, Length, Mass, Qty, Time, Velocity};

/// Three components of one dimension, stored in SI base units.
#[derive(Clone, Copy, PartialEq, Default)]
pub struct QVec3<
    const L: i8,
    const M: i8,
    const T: i8,
    const I: i8,
    const K: i8,
    const N: i8,
    const J: i8,
>(DVec3);

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    QVec3<L, M, T, I, K, N, J>
{
    pub const ZERO: Self = QVec3(DVec3::ZERO);

    /// Wrap a vector already in SI base units.
    pub fn from_si(v: DVec3) -> Self {
        QVec3(v)
    }

    /// The components in SI base units.
    pub fn to_si(self) -> DVec3 {
        self.0
    }

    pub fn new(
        x: Qty<L, M, T, I, K, N, J>,
        y: Qty<L, M, T, I, K, N, J>,
        z: Qty<L, M, T, I, K, N, J>,
    ) -> Self {
        QVec3(DVec3::new(x.to_si(), y.to_si(), z.to_si()))
    }

    pub fn splat(v: Qty<L, M, T, I, K, N, J>) -> Self {
        QVec3(DVec3::splat(v.to_si()))
    }

    pub fn x(self) -> Qty<L, M, T, I, K, N, J> {
        Qty::from_si(self.0.x)
    }

    pub fn y(self) -> Qty<L, M, T, I, K, N, J> {
        Qty::from_si(self.0.y)
    }

    pub fn z(self) -> Qty<L, M, T, I, K, N, J> {
        Qty::from_si(self.0.z)
    }

    /// Magnitude, which keeps the dimension.
    pub fn length(self) -> Qty<L, M, T, I, K, N, J> {
        Qty::from_si(self.0.length())
    }

    /// Which way it points — a pure number, because a direction is a length over
    /// a length. Zero-length vectors give zero rather than a NaN.
    pub fn normalize(self) -> DVec3 {
        self.0.normalize_or_zero()
    }

    /// The component along a unit direction. Projection does not change the
    /// dimension, which is why this one is expressible and `dot` is not.
    pub fn along(self, direction: DVec3) -> Qty<L, M, T, I, K, N, J> {
        Qty::from_si(self.0.dot(direction))
    }

    /// The part of this vector perpendicular to a unit direction.
    pub fn perpendicular_to(self, direction: DVec3) -> Self {
        QVec3(self.0 - direction * self.0.dot(direction))
    }

    pub fn is_finite(self) -> bool {
        self.0.is_finite()
    }

    pub fn lerp(self, other: Self, t: f64) -> Self {
        QVec3(self.0 + (other.0 - self.0) * t)
    }
}

macro_rules! generic_vec_op {
    ($trait:ident, $method:ident, $op:tt) => {
        impl<
                const L: i8,
                const M: i8,
                const T: i8,
                const I: i8,
                const K: i8,
                const N: i8,
                const J: i8,
            > $trait for QVec3<L, M, T, I, K, N, J>
        {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self {
                QVec3(self.0 $op rhs.0)
            }
        }
    };
}

generic_vec_op!(Add, add, +);
generic_vec_op!(Sub, sub, -);

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    AddAssign for QVec3<L, M, T, I, K, N, J>
{
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    SubAssign for QVec3<L, M, T, I, K, N, J>
{
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8> Neg
    for QVec3<L, M, T, I, K, N, J>
{
    type Output = Self;
    fn neg(self) -> Self {
        QVec3(-self.0)
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    Mul<f64> for QVec3<L, M, T, I, K, N, J>
{
    type Output = Self;
    fn mul(self, k: f64) -> Self {
        QVec3(self.0 * k)
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    Div<f64> for QVec3<L, M, T, I, K, N, J>
{
    type Output = Self;
    fn div(self, k: f64) -> Self {
        QVec3(self.0 / k)
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    Mul<QVec3<L, M, T, I, K, N, J>> for f64
{
    type Output = QVec3<L, M, T, I, K, N, J>;
    fn mul(self, v: QVec3<L, M, T, I, K, N, J>) -> QVec3<L, M, T, I, K, N, J> {
        QVec3(v.0 * self)
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    core::fmt::Debug for QVec3<L, M, T, I, K, N, J>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "[{}, {}, {}", self.0.x, self.0.y, self.0.z)?;
        for (symbol, exponent) in [
            ("m", L),
            ("kg", M),
            ("s", T),
            ("A", I),
            ("K", K),
            ("mol", N),
            ("cd", J),
        ] {
            match exponent {
                0 => {}
                1 => write!(f, "·{symbol}")?,
                e => write!(f, "·{symbol}^{e}")?,
            }
        }
        write!(f, "]")
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    Serialize for QVec3<L, M, T, I, K, N, J>
{
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        [self.0.x, self.0.y, self.0.z].serialize(s)
    }
}

impl<
        'de,
        const L: i8,
        const M: i8,
        const T: i8,
        const I: i8,
        const K: i8,
        const N: i8,
        const J: i8,
    > Deserialize<'de> for QVec3<L, M, T, I, K, N, J>
{
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        <[f64; 3]>::deserialize(d).map(|[x, y, z]| QVec3(DVec3::new(x, y, z)))
    }
}

/// A position or a displacement, m. The same dimension, and deliberately the same
/// type: the difference between them is a choice of origin, not of physics.
pub type LengthVec = QVec3<1, 0, 0, 0, 0, 0, 0>;
pub type VelocityVec = QVec3<1, 0, -1, 0, 0, 0, 0>;
pub type AccelerationVec = QVec3<1, 0, -2, 0, 0, 0, 0>;
pub type ForceVec = QVec3<1, 1, -2, 0, 0, 0, 0>;
pub type MomentumVec = QVec3<1, 1, -1, 0, 0, 0, 0>;

/// A scalar times a vector, and the division that undoes it.
macro_rules! scaled_by {
    ($vec:ty, $scalar:ty => $out:ty) => {
        impl Mul<$scalar> for $vec {
            type Output = $out;
            fn mul(self, k: $scalar) -> $out {
                QVec3(self.0 * k.to_si())
            }
        }
        impl Mul<$vec> for $scalar {
            type Output = $out;
            fn mul(self, v: $vec) -> $out {
                QVec3(v.0 * self.to_si())
            }
        }
        impl Div<$scalar> for $out {
            type Output = $vec;
            fn div(self, k: $scalar) -> $vec {
                QVec3(self.0 / k.to_si())
            }
        }
    };
}

scaled_by!(VelocityVec, Time => LengthVec);
scaled_by!(AccelerationVec, Time => VelocityVec);
scaled_by!(ForceVec, Time => MomentumVec);
scaled_by!(VelocityVec, Mass => MomentumVec);
scaled_by!(AccelerationVec, Mass => ForceVec);

impl LengthVec {
    pub fn mm(x: f64, y: f64, z: f64) -> LengthVec {
        QVec3(DVec3::new(x, y, z) * 1e-3)
    }
    pub fn m(x: f64, y: f64, z: f64) -> LengthVec {
        QVec3(DVec3::new(x, y, z))
    }
    pub fn in_mm(self) -> DVec3 {
        self.0 * 1e3
    }
}

impl VelocityVec {
    pub fn mm_per_s(x: f64, y: f64, z: f64) -> VelocityVec {
        QVec3(DVec3::new(x, y, z) * 1e-3)
    }
    pub fn m_per_s(x: f64, y: f64, z: f64) -> VelocityVec {
        QVec3(DVec3::new(x, y, z))
    }
}

/// Distance between two points, which is what a length actually measures.
pub fn distance(a: LengthVec, b: LengthVec) -> Length {
    (a - b).length()
}

/// Newton's second law, with the dimensions doing the checking.
pub fn newton_second(mass: Mass, acceleration: AccelerationVec) -> ForceVec {
    mass * acceleration
}

/// Momentum of a moving mass.
pub fn momentum(mass: Mass, velocity: VelocityVec) -> MomentumVec {
    mass * velocity
}

/// Kinetic energy, ½mv². Needs the squared magnitude, so it is written out here
/// rather than falling out of an operator.
pub fn kinetic_energy(mass: Mass, velocity: VelocityVec) -> crate::Energy {
    let v = velocity.to_si().length();
    Qty::from_si(0.5 * mass.to_si() * v * v)
}

/// Speed acquired, and distance covered, under a constant acceleration.
pub fn free_travel(v0: VelocityVec, a: AccelerationVec, t: Time) -> (VelocityVec, LengthVec) {
    let v = v0 + a * t;
    let x = v0 * t + (a * t) * t * 0.5;
    (v, x)
}

/// Force needed to hold `mass` in a circle — the check that a rotating stage's
/// bearing can take what a scan rate asks of it.
pub fn centripetal(mass: Mass, speed: Velocity, radius: Length) -> Force {
    Qty::from_si(mass.to_si() * speed.to_si() * speed.to_si() / radius.to_si())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Energy, Temperature};

    /// A direction is dimensionless, which means a ray direction and a
    /// displacement are different types and cannot be swapped.
    #[test]
    fn normalising_strips_the_dimension() {
        let d = LengthVec::mm(3.0, 4.0, 0.0);
        assert!((d.length().in_mm() - 5.0).abs() < 1e-12);
        let dir: DVec3 = d.normalize();
        assert!((dir.length() - 1.0).abs() < 1e-15);
        assert!((dir - DVec3::new(0.6, 0.8, 0.0)).length() < 1e-15);
        // A zero vector has no direction, and says so rather than producing NaN.
        assert_eq!(LengthVec::ZERO.normalize(), DVec3::ZERO);
    }

    /// Projection preserves the dimension, which is what makes it expressible
    /// where `dot` is not.
    #[test]
    fn projection_keeps_the_dimension() {
        let v = VelocityVec::mm_per_s(120.0, 0.0, -5.0);
        let axis = DVec3::X;
        let along: Velocity = v.along(axis);
        assert!((along.to_si() - 0.12).abs() < 1e-15);
        let across = v.perpendicular_to(axis);
        assert!((across.along(axis).to_si()).abs() < 1e-15);
        // The two parts add back up to the whole.
        let rebuilt = across + VelocityVec::from_si(axis * along.to_si());
        assert!((rebuilt - v).length().to_si() < 1e-15);
    }

    /// Kinematics with the dimensions doing the bookkeeping: 1 g for 2 s gives
    /// 19.6 m/s and 19.6 m, which are different numbers of different dimensions
    /// that happen to share digits.
    #[test]
    fn constant_acceleration_is_dimensionally_checked() {
        let a = AccelerationVec::from_si(DVec3::new(0.0, -crate::G0.to_si(), 0.0));
        let (v, x) = free_travel(VelocityVec::ZERO, a, Time::s(2.0));
        assert!((v.length().to_si() - 19.6133).abs() < 1e-3, "{v:?}");
        assert!((x.length().to_si() - 19.6133).abs() < 1e-3, "{x:?}");
        // Downwards, both of them.
        assert!(v.y().to_si() < 0.0 && x.y().to_si() < 0.0);
    }

    /// Newton's second law and the energy it does, checked against the closed
    /// form: work done equals the kinetic energy gained.
    #[test]
    fn work_equals_the_kinetic_energy_it_bought() {
        let m = Mass::kg(2.0);
        let a = AccelerationVec::from_si(DVec3::X * 3.0);
        let f: ForceVec = newton_second(m, a);
        assert!((f.length().to_si() - 6.0).abs() < 1e-12);

        let t = Time::s(4.0);
        let (v, x) = free_travel(VelocityVec::ZERO, a, t);
        let work: Energy = Qty::from_si(f.along(DVec3::X).to_si() * x.along(DVec3::X).to_si());
        let ke = kinetic_energy(m, v);
        assert!(
            (work - ke).abs().to_si() < 1e-9,
            "work {work:?} should equal kinetic energy {ke:?}"
        );
    }

    /// Momentum is conserved in a collision, and the type system will not let a
    /// velocity be added to a momentum on the way there.
    #[test]
    fn momentum_adds_across_a_collision() {
        let p1 = momentum(Mass::kg(2.0), VelocityVec::m_per_s(3.0, 0.0, 0.0));
        let p2 = momentum(Mass::kg(1.0), VelocityVec::m_per_s(-4.0, 0.0, 0.0));
        let total: MomentumVec = p1 + p2;
        assert!((total.along(DVec3::X).to_si() - 2.0).abs() < 1e-12);
        // The combined mass therefore moves at 2/3 m/s.
        let after: VelocityVec = total / Mass::kg(3.0);
        assert!((after.along(DVec3::X).to_si() - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn vectors_round_trip_through_json() {
        let v = LengthVec::mm(1.0, 2.0, 3.0);
        let json = serde_json::to_string(&v).unwrap();
        assert_eq!(json, "[0.001,0.002,0.003]");
        assert_eq!(serde_json::from_str::<LengthVec>(&json).unwrap(), v);
    }

    #[test]
    fn debug_shows_the_dimension() {
        assert_eq!(
            format!("{:?}", ForceVec::from_si(DVec3::new(1.0, 0.0, 0.0))),
            "[1, 0, 0·m·kg·s^-2]"
        );
    }

    /// The whole point, stated as a compile-time fact rather than a runtime one:
    /// these lines do not compile, and the comments are the test.
    #[test]
    fn wrong_dimensions_do_not_compile() {
        let _ = LengthVec::mm(1.0, 0.0, 0.0);
        let _ = VelocityVec::mm_per_s(1.0, 0.0, 0.0);
        let _ = Temperature::kelvin(300.0);
        // let _ = _position + _velocity;        // mismatched types
        // let _ = _position.along(_velocity);   // expected DVec3
        // let _: Length = _temperature;         // mismatched types
    }
}
