//! Dimensional analysis: physical quantities that refuse to be added wrongly.
//!
//! ```
//! use dualis_units::{Area, Energy, Length, Mass, Power, SpecificHeat, Temperature, Time};
//!
//! // A unit-bearing constructor is the only place a factor of a thousand may appear.
//! let side = Length::mm(10.0);
//! let area: Area = side * side;                       // the dimension follows the product
//! assert!((area.to_si() - 1e-4).abs() < 1e-18);
//!
//! // Absorbed power over a time is an energy, and the type says so without being told.
//! let absorbed = Power::mw(96.0);
//! let heat: Energy = absorbed * Time::s(1.0);
//!
//! // Divide it by a heat capacity and a temperature comes out.
//! let capacity = Mass::g(2.0) * SpecificHeat::j_per_kg_k(858.0);
//! let rise: Temperature = heat / capacity;
//! assert!((rise.to_si() - 0.05594).abs() < 1e-4);
//! ```
//!
//! And the mistake the whole crate exists to prevent does not compile:
//!
//! ```compile_fail
//! use dualis_units::{Length, Time};
//! let nonsense = Length::mm(3.0) + Time::s(1.0);
//! ```
//!
//! One domain can get away with a convention. `dualis-core` began as optics and
//! said "millimetres, nanometres and seconds, everywhere" in a doc comment, and
//! that held because every number in the crate was a length, a wavelength or a
//! fraction. It stops holding the moment a second domain arrives: a kelvin, a
//! newton and a watt are all `f64`, they all add, and the compiler and the tests
//! both stay green while the physics goes wrong.
//!
//! So dimension lives in the type. [`Qty`] carries the seven SI base exponents as
//! const generic parameters, which makes `Length + Time` a compile error and
//! `Force * Length` an [`Energy`] — and costs nothing at runtime, since a `Qty`
//! is an `f64` and every operation on it is the `f64` operation.
//!
//! # Storage is always SI base units
//!
//! A `Qty` holds metres, kilograms, seconds, amperes, kelvin, moles, candela —
//! never millimetres, never nanometres. Those are *entry and exit* forms:
//!
//! ```
//! use dualis_units::{Length, Time, Velocity};
//!
//! let d = Length::mm(120.0);
//! let t = Time::ms(4.0);
//! let v: Velocity = d / t;
//! assert!((v.to_si() - 30.0).abs() < 1e-12);   // 30 m/s
//! assert!((d.in_nm() - 1.2e8).abs() < 1.0);
//! ```
//!
//! That way there is exactly one representation to reason about, and the
//! unit-bearing constructors are the only place a factor of 1000 can hide.
//!
//! # What this cannot do
//!
//! **Angles are dimensionless**, so [`Frequency`] and an angular velocity are the
//! same type — SI says radians are m/m, and no dimensional system can separate
//! them. Same for torque and energy. Where that distinction matters, it has to be
//! carried by a newtype in the domain crate, not here.
//!
//! **Only declared products compose.** `Length * Length` is an [`Area`] because
//! that pair is written down below. Deriving arbitrary products would need
//! arithmetic on const generic parameters, which is unstable, so the alternative
//! to a declared list is a dependency on `uom`. The list is cheap to extend, and
//! anything undeclared can always go through [`Qty::from_si`].

// Every public item carries a doc comment. Denied rather than warned: a public physics API
// whose `Length::mm` shows a blank summary in rustdoc is documented in the sense that a
// paragraph exists somewhere, and not in the sense a reader needs.
#![deny(missing_docs)]
#![forbid(unsafe_code)]

use core::fmt;
use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use serde::{Deserialize, Deserializer, Serialize, Serializer};

pub mod vector;
pub use vector::{AccelerationVec, ForceVec, LengthVec, MomentumVec, QVec3, VelocityVec};

/// A quantity, with the seven SI base dimensions in its type.
///
/// The parameters are the exponents of metre, kilogram, second, ampere, kelvin,
/// mole and candela, in that order, so a velocity (m·s⁻¹) is `Qty<1,0,-1,0,0,0,0>`
/// — which is what [`Velocity`] names.
///
/// Addition, subtraction, negation, comparison and scaling by a plain `f64` work
/// for every dimension. Multiplication and division between two quantities work
/// for the pairs declared in this module.
#[derive(Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct Qty<
    const L: i8,
    const M: i8,
    const T: i8,
    const I: i8,
    const K: i8,
    const N: i8,
    const J: i8,
>(f64);

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    Qty<L, M, T, I, K, N, J>
{
    /// Zero, which is the one value every dimension shares.
    pub const ZERO: Self = Qty(0.0);

    /// Wrap a number already in SI base units. The escape hatch: use it when a
    /// dimension has no name here, and name it if you use it twice.
    ///
    /// `const`, so a dimensioned constant can be written without a lazy static.
    pub const fn from_si(value: f64) -> Self {
        Qty(value)
    }

    /// The value in SI base units.
    pub const fn to_si(self) -> f64 {
        self.0
    }

    /// The seven exponents, for diagnostics and for a runtime dimension check at
    /// a boundary the type system does not cross (deserialisation, FFI).
    pub const fn dimension() -> [i8; 7] {
        [L, M, T, I, K, N, J]
    }

    /// Magnitude without its sign, in the same dimension.
    pub fn abs(self) -> Self {
        Qty(self.0.abs())
    }

    /// The smaller of two quantities of the same dimension.
    pub fn min(self, other: Self) -> Self {
        Qty(self.0.min(other.0))
    }

    /// The larger of two quantities of the same dimension.
    pub fn max(self, other: Self) -> Self {
        Qty(self.0.max(other.0))
    }

    /// Whether the magnitude is neither infinite nor NaN.
    ///
    /// Worth checking where a limit is reported rather than computed: several methods here
    /// return an infinity to mean "no limit", which is honest but arithmetic on it is not.
    pub fn is_finite(self) -> bool {
        self.0.is_finite()
    }

    /// Sign of the magnitude, as a plain number — a sign has no dimension.
    pub fn signum(self) -> f64 {
        self.0.signum()
    }

    /// Linear interpolation, which stays within the dimension.
    pub fn lerp(self, other: Self, t: f64) -> Self {
        Qty(self.0 + (other.0 - self.0) * t)
    }
}

// ---------------------------------------------------------------------------
// Dimension-preserving arithmetic: works for every dimension at once, because
// none of it changes the exponents.
// ---------------------------------------------------------------------------

macro_rules! generic_op {
    ($trait:ident, $method:ident, $op:tt) => {
        impl<
                const L: i8,
                const M: i8,
                const T: i8,
                const I: i8,
                const K: i8,
                const N: i8,
                const J: i8,
            > $trait for Qty<L, M, T, I, K, N, J>
        {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self {
                Qty(self.0 $op rhs.0)
            }
        }
    };
}

generic_op!(Add, add, +);
generic_op!(Sub, sub, -);

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    AddAssign for Qty<L, M, T, I, K, N, J>
{
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    SubAssign for Qty<L, M, T, I, K, N, J>
{
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8> Neg
    for Qty<L, M, T, I, K, N, J>
{
    type Output = Self;
    fn neg(self) -> Self {
        Qty(-self.0)
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    Mul<f64> for Qty<L, M, T, I, K, N, J>
{
    type Output = Self;
    fn mul(self, k: f64) -> Self {
        Qty(self.0 * k)
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    Div<f64> for Qty<L, M, T, I, K, N, J>
{
    type Output = Self;
    fn div(self, k: f64) -> Self {
        Qty(self.0 / k)
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    Mul<Qty<L, M, T, I, K, N, J>> for f64
{
    type Output = Qty<L, M, T, I, K, N, J>;
    fn mul(self, q: Qty<L, M, T, I, K, N, J>) -> Qty<L, M, T, I, K, N, J> {
        Qty(self * q.0)
    }
}

/// Dividing two quantities of the *same* dimension gives a plain number — which
/// is the one product rule that needs no exponent arithmetic, and the one every
/// tolerance check uses.
impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8> Div
    for Qty<L, M, T, I, K, N, J>
{
    type Output = f64;
    fn div(self, rhs: Self) -> f64 {
        self.0 / rhs.0
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    fmt::Debug for Qty<L, M, T, I, K, N, J>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)?;
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
        Ok(())
    }
}

impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    fmt::Display for Qty<L, M, T, I, K, N, J>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

// Serialised as the bare SI number: a scene file stays readable, and the
// dimension is carried by the field's type rather than repeated in the data.
impl<const L: i8, const M: i8, const T: i8, const I: i8, const K: i8, const N: i8, const J: i8>
    Serialize for Qty<L, M, T, I, K, N, J>
{
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.0.serialize(s)
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
    > Deserialize<'de> for Qty<L, M, T, I, K, N, J>
{
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        f64::deserialize(d).map(Qty)
    }
}

// ---------------------------------------------------------------------------
// The dimensions themselves.
// ---------------------------------------------------------------------------

/// A pure ratio: reflectance, duty cycle, refractive index, Strehl.
pub type Dimensionless = Qty<0, 0, 0, 0, 0, 0, 0>;

/// Metres.
pub type Length = Qty<1, 0, 0, 0, 0, 0, 0>;
/// Kilograms.
pub type Mass = Qty<0, 1, 0, 0, 0, 0, 0>;
/// Seconds.
pub type Time = Qty<0, 0, 1, 0, 0, 0, 0>;
/// Amperes.
pub type Current = Qty<0, 0, 0, 1, 0, 0, 0>;
/// Absolute temperature. Kelvin only — see [`Temperature::celsius`].
pub type Temperature = Qty<0, 0, 0, 0, 1, 0, 0>;
/// Moles.
pub type Amount = Qty<0, 0, 0, 0, 0, 1, 0>;
/// Candelas.
pub type LuminousIntensity = Qty<0, 0, 0, 0, 0, 0, 1>;

/// Square metres.
pub type Area = Qty<2, 0, 0, 0, 0, 0, 0>;
/// Cubic metres.
pub type Volume = Qty<3, 0, 0, 0, 0, 0, 0>;
/// Metres per second.
pub type Velocity = Qty<1, 0, -1, 0, 0, 0, 0>;
/// Metres per second squared.
pub type Acceleration = Qty<1, 0, -2, 0, 0, 0, 0>;
/// kg·m·s⁻¹ — mass times velocity, and the thing a closed system conserves
/// exactly rather than nearly.
pub type Momentum = Qty<1, 1, -1, 0, 0, 0, 0>;
/// Newtons.
pub type Force = Qty<1, 1, -2, 0, 0, 0, 0>;
/// Pascals. Also the unit of an energy density and of a stress, which are the same
/// dimension and not a coincidence.
pub type Pressure = Qty<-1, 1, -2, 0, 0, 0, 0>;
/// Joules.
pub type Energy = Qty<2, 1, -2, 0, 0, 0, 0>;
/// Watts.
pub type Power = Qty<2, 1, -3, 0, 0, 0, 0>;
/// kg·m⁻³. Note that a glass catalogue quotes g/cm³, a factor of a thousand away —
/// see [`Density::g_per_cm3`].
pub type Density = Qty<-3, 1, 0, 0, 0, 0, 0>;
/// Cycles per second. Dimensionally identical to an angular velocity, since a
/// radian is m/m — the type system cannot and should not pretend otherwise.
pub type Frequency = Qty<0, 0, -1, 0, 0, 0, 0>;
/// Power per unit area, W·m⁻². What a detector face actually receives.
pub type Irradiance = Qty<0, 1, -3, 0, 0, 0, 0>;
/// W·m⁻¹·K⁻¹ — the `k` of Fourier's law.
pub type ThermalConductivity = Qty<1, 1, -3, 0, -1, 0, 0>;
/// W·K⁻¹ — how fast heat crosses a joint, `UA`.
///
/// Dimensionally [`Power`] per [`Temperature`], and equivalently [`ThermalConductivity`] times
/// a [`Length`], which is the physically meaningful reading: `kA/L`. It is what a *contact*
/// resistance is measured in — a bolted joint, a winding pressed into a stator — and those have
/// no bulk conductivity to be derived from, which is why the quantity exists in its own right.
pub type Conductance = Qty<2, 1, -3, 0, -1, 0, 0>;
/// J·kg⁻¹·K⁻¹ — the `c_p` that says how much heat a gram of glass can hide.
pub type SpecificHeat = Qty<2, 0, -2, 0, -1, 0, 0>;
/// m²·s⁻¹ — thermal diffusivity `α = k/(ρ c_p)`, and also mass diffusivity.
pub type Diffusivity = Qty<2, 0, -1, 0, 0, 0, 0>;
/// K⁻¹ — the coefficient that turns absorbed light into a focus shift.
pub type ThermalExpansion = Qty<0, 0, 0, 0, -1, 0, 0>;
/// kg·m² — how hard a body is to spin up about an axis.
///
/// The rotational counterpart of mass, and unlike mass it depends on the axis: a
/// pencil is trivial to spin about its length and awkward about its middle. That
/// direction-dependence is why it is a tensor and why a free body's rotation is
/// interesting rather than uniform.
pub type MomentOfInertia = Qty<2, 1, 0, 0, 0, 0, 0>;
/// kg·m²·s⁻¹ — the rotational counterpart of momentum, and conserved for the same
/// reason.
pub type AngularMomentum = Qty<2, 1, -1, 0, 0, 0, 0>;
/// N·m⁻¹ — a spring's `k`, and the penalty stiffness a contact is modelled with.
///
/// This is what sets a mechanical solver's stability limit: a mass on a spring
/// oscillates with period `2π√(m/k)`, and an explicit integrator has to resolve that
/// period whether or not anyone cares about it. Stiff contact is expensive for
/// exactly this reason.
pub type Stiffness = Qty<0, 1, -2, 0, 0, 0, 0>;
/// N·s·m⁻¹ — a dashpot's `c`. Force proportional to velocity, and the only place a
/// mechanical simulation loses energy on purpose.
pub type Damping = Qty<0, 1, -1, 0, 0, 0, 0>;
/// Coulombs.
pub type Charge = Qty<0, 0, 1, 1, 0, 0, 0>;
/// Volts.
pub type Voltage = Qty<2, 1, -3, -1, 0, 0, 0>;
/// J·K⁻¹ — mass times specific heat. How much heat a thing can hide before it
/// shows up as a temperature.
pub type HeatCapacity = Qty<2, 1, -2, 0, -1, 0, 0>;

// ---------------------------------------------------------------------------
// Declared products. Each line also gives the two divisions that undo it.
// ---------------------------------------------------------------------------

macro_rules! product {
    ($a:ty, $b:ty => $c:ty) => {
        impl Mul<$b> for $a {
            type Output = $c;
            fn mul(self, rhs: $b) -> $c {
                Qty(self.0 * rhs.0)
            }
        }
        impl Mul<$a> for $b {
            type Output = $c;
            fn mul(self, rhs: $a) -> $c {
                Qty(self.0 * rhs.0)
            }
        }
        impl Div<$b> for $c {
            type Output = $a;
            fn div(self, rhs: $b) -> $a {
                Qty(self.0 / rhs.0)
            }
        }
        impl Div<$a> for $c {
            type Output = $b;
            fn div(self, rhs: $a) -> $b {
                Qty(self.0 / rhs.0)
            }
        }
    };
}

macro_rules! square {
    ($a:ty => $c:ty) => {
        impl Mul<$a> for $a {
            type Output = $c;
            fn mul(self, rhs: $a) -> $c {
                Qty(self.0 * rhs.0)
            }
        }
        impl Div<$a> for $c {
            type Output = $a;
            fn div(self, rhs: $a) -> $a {
                Qty(self.0 / rhs.0)
            }
        }
    };
}

square!(Length => Area);
product!(Area, Length => Volume);
product!(Velocity, Time => Length);
product!(Acceleration, Time => Velocity);
product!(Mass, Acceleration => Force);
product!(Mass, Velocity => Momentum);
product!(Force, Length => Energy);
product!(Force, Time => Momentum);
product!(Pressure, Area => Force);
product!(Power, Time => Energy);
product!(Irradiance, Area => Power);
product!(Density, Volume => Mass);
product!(Current, Time => Charge);
product!(Voltage, Current => Power);
product!(Mass, Area => MomentOfInertia);
product!(MomentOfInertia, Frequency => AngularMomentum);
product!(Stiffness, Length => Force);
product!(Damping, Velocity => Force);
product!(Mass, SpecificHeat => HeatCapacity);
// UA·ΔT is watts, and C/UA is a time — the two identities a thermal network is built out of,
// so the type system checks them rather than a comment claiming them.
product!(Conductance, Temperature => Power);
product!(Conductance, Time => HeatCapacity);
product!(HeatCapacity, Temperature => Energy);
product!(Frequency, Time => Dimensionless);

impl Area {
    /// The side of a square of this area. The one root worth naming, because it
    /// is how a beam radius comes back out of a spot area.
    pub fn sqrt(self) -> Length {
        Qty(self.0.sqrt())
    }
}

// ---------------------------------------------------------------------------
// Unit-bearing entry and exit. The only place a factor of 1000 may appear.
// ---------------------------------------------------------------------------

impl Length {
    /// Metres.
    pub fn m(v: f64) -> Length {
        Qty(v)
    }
    /// Millimetres.
    pub fn mm(v: f64) -> Length {
        Qty(v * 1e-3)
    }
    /// Micrometres.
    pub fn um(v: f64) -> Length {
        Qty(v * 1e-6)
    }
    /// Nanometres. The wavelength unit, and why every `Spectrum` field is named `_nm`.
    pub fn nm(v: f64) -> Length {
        Qty(v * 1e-9)
    }
    /// As millimetres.
    pub fn in_mm(self) -> f64 {
        self.0 * 1e3
    }
    /// As micrometres.
    pub fn in_um(self) -> f64 {
        self.0 * 1e6
    }
    /// As nanometres.
    pub fn in_nm(self) -> f64 {
        self.0 * 1e9
    }
}

impl Time {
    /// Seconds.
    pub fn s(v: f64) -> Time {
        Qty(v)
    }
    /// Milliseconds.
    pub fn ms(v: f64) -> Time {
        Qty(v * 1e-3)
    }
    /// Microseconds.
    pub fn us(v: f64) -> Time {
        Qty(v * 1e-6)
    }
    /// Nanoseconds.
    pub fn ns(v: f64) -> Time {
        Qty(v * 1e-9)
    }
    /// As milliseconds.
    pub fn in_ms(self) -> f64 {
        self.0 * 1e3
    }
    /// As microseconds.
    pub fn in_us(self) -> f64 {
        self.0 * 1e6
    }
}

impl Temperature {
    /// Kelvin, which is what is stored.
    pub fn kelvin(v: f64) -> Temperature {
        Qty(v)
    }
    /// Celsius is an *offset* scale, not a scaled one, which is why it gets a
    /// named constructor rather than a factor: 20 °C is 293.15 K, and a
    /// temperature *difference* of 20 K is a different thing entirely.
    pub fn celsius(v: f64) -> Temperature {
        Qty(v + 273.15)
    }
    /// As degrees Celsius. Subtracts the offset; see [`Temperature::celsius`].
    pub fn in_celsius(self) -> f64 {
        self.0 - 273.15
    }
}

impl Mass {
    /// Kilograms.
    pub fn kg(v: f64) -> Mass {
        Qty(v)
    }
    /// Grams.
    pub fn g(v: f64) -> Mass {
        Qty(v * 1e-3)
    }
}

impl Density {
    /// The way a glass catalogue quotes it: N-BK7 is 2.51 g/cm³.
    pub fn g_per_cm3(v: f64) -> Density {
        Qty(v * 1e3)
    }
    /// Kilograms per cubic metre, which is what is stored.
    pub fn kg_per_m3(v: f64) -> Density {
        Qty(v)
    }
}

impl Power {
    /// Watts.
    pub fn w(v: f64) -> Power {
        Qty(v)
    }
    /// Milliwatts.
    pub fn mw(v: f64) -> Power {
        Qty(v * 1e-3)
    }
    /// Microwatts.
    pub fn uw(v: f64) -> Power {
        Qty(v * 1e-6)
    }
    /// As milliwatts.
    pub fn in_mw(self) -> f64 {
        self.0 * 1e3
    }
}

impl Energy {
    /// Joules.
    pub fn j(v: f64) -> Energy {
        Qty(v)
    }
    /// Millijoules.
    pub fn mj(v: f64) -> Energy {
        Qty(v * 1e-3)
    }
}

impl Frequency {
    /// Hertz.
    pub fn hz(v: f64) -> Frequency {
        Qty(v)
    }
    /// Kilohertz.
    pub fn khz(v: f64) -> Frequency {
        Qty(v * 1e3)
    }
    /// Megahertz.
    pub fn mhz(v: f64) -> Frequency {
        Qty(v * 1e6)
    }
    /// Period: one over the frequency. Named because `1.0 / f` cannot typecheck.
    pub fn period(self) -> Time {
        Qty(1.0 / self.0)
    }
}

impl Velocity {
    /// Metres per second.
    pub fn m_per_s(v: f64) -> Velocity {
        Qty(v)
    }
    /// Millimetres per second.
    pub fn mm_per_s(v: f64) -> Velocity {
        Qty(v * 1e-3)
    }
}

impl Irradiance {
    /// Watts per square metre, which is what is stored.
    pub fn w_per_m2(v: f64) -> Irradiance {
        Qty(v)
    }
    /// How an illumination spec is usually written: mW/cm².
    pub fn mw_per_cm2(v: f64) -> Irradiance {
        Qty(v * 10.0)
    }
}

impl ThermalConductivity {
    /// W·m⁻¹·K⁻¹, the unit a materials table uses.
    pub fn w_per_m_k(v: f64) -> ThermalConductivity {
        Qty(v)
    }
}

impl SpecificHeat {
    /// J·kg⁻¹·K⁻¹, the unit a materials table uses.
    pub fn j_per_kg_k(v: f64) -> SpecificHeat {
        Qty(v)
    }
}

impl Conductance {
    /// Watts per kelvin.
    pub fn w_per_k(v: f64) -> Conductance {
        Qty(v)
    }
}

impl HeatCapacity {
    /// Joules per kelvin. The companion to [`Conductance::w_per_k`]: their ratio is a time
    /// constant, and the type system says so.
    pub fn j_per_k(v: f64) -> HeatCapacity {
        Qty(v)
    }
}

impl ThermalExpansion {
    /// Catalogues quote it in parts per million per kelvin: N-BK7 is 7.1.
    pub fn ppm_per_k(v: f64) -> ThermalExpansion {
        Qty(v * 1e-6)
    }
}

impl Dimensionless {
    /// A bare ratio, for the one case where a number genuinely has no dimension:
    /// a reflectance, a duty cycle, a refractive index.
    pub fn ratio(v: f64) -> Dimensionless {
        Qty(v)
    }
}

// ---------------------------------------------------------------------------
// Physical constants, in SI base units, so that a formula written with them
// carries its own dimensional proof.
// ---------------------------------------------------------------------------

/// Speed of light in vacuum, m·s⁻¹ (exact by definition).
pub const C: Velocity = Qty(299_792_458.0);
/// Planck constant, J·s (exact by definition).
pub const PLANCK: Qty<2, 1, -1, 0, 0, 0, 0> = Qty(6.626_070_15e-34);
/// Boltzmann constant, J·K⁻¹ (exact by definition).
pub const BOLTZMANN: HeatCapacity = Qty(1.380_649e-23);
/// Stefan-Boltzmann constant, W·m⁻²·K⁻⁴ — radiative exchange lives on this.
pub const STEFAN_BOLTZMANN: Qty<0, 1, -3, 0, -4, 0, 0> = Qty(5.670_374_419e-8);
/// Standard gravity, m·s⁻².
pub const G0: Acceleration = Qty(9.806_65);

/// Energy of one photon at a vacuum wavelength: `E = hc/λ`.
///
/// The bridge between a spectrum and a photon count, and the reason a detector's
/// response is not the same shape as a lamp's output.
pub fn photon_energy(wavelength: Length) -> Energy {
    Qty(PLANCK.0 * C.0 / wavelength.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The point of the crate: a product of dimensions lands on the type that
    /// names it, whichever route it took there.
    #[test]
    fn products_land_on_the_named_dimension() {
        let m = Mass::kg(2.0);
        let a = Acceleration::from_si(3.0);
        let f: Force = m * a;
        assert!((f.to_si() - 6.0).abs() < 1e-12);

        // Two different routes to the same energy, and they unify.
        let by_work: Energy = f * Length::m(4.0);
        let by_power: Energy = Power::w(24.0) * Time::s(1.0);
        assert!((by_work - by_power).abs().to_si() < 1e-12);

        // And multiplication commutes, as it must.
        let swapped: Force = a * m;
        assert_eq!(f, swapped);
    }

    /// Millimetres and nanometres are entry forms only; storage is metres. A
    /// wavelength and a lens diameter therefore compare correctly without anyone
    /// remembering which convention each was written in.
    #[test]
    fn unit_prefixes_are_only_a_doorway() {
        let lens = Length::mm(25.4);
        let green = Length::nm(550.0);
        assert!(lens > green);
        assert!((lens.to_si() - 0.0254).abs() < 1e-15);
        assert!((green.in_nm() - 550.0).abs() < 1e-9);
        // 25.4 mm is 46181.8... wavelengths of green light.
        let waves = lens / green;
        assert!((waves - 46_181.8).abs() < 0.1, "got {waves}");
    }

    /// Dividing like by like gives a plain number, which is what every
    /// tolerance and every reflectance is.
    #[test]
    fn like_over_like_is_a_bare_number() {
        let reflected = Power::mw(0.42);
        let incident = Power::mw(10.0);
        let r: f64 = reflected / incident;
        assert!((r - 0.042).abs() < 1e-12);
    }

    /// Celsius offsets rather than scales, and getting that wrong is a 273 K
    /// error that no dimensional check would ever catch.
    #[test]
    fn celsius_is_an_offset_not_a_factor() {
        assert!((Temperature::celsius(20.0).to_si() - 293.15).abs() < 1e-12);
        assert!((Temperature::kelvin(293.15).in_celsius() - 20.0).abs() < 1e-12);
        // A *difference* of 20 K is not 293.15 K, and only one of these is a
        // temperature you can put in Stefan-Boltzmann.
        let rise = Temperature::kelvin(313.15) - Temperature::kelvin(293.15);
        assert!((rise.to_si() - 20.0).abs() < 1e-12);
    }

    /// The optics-to-thermal chain this whole crate exists to make safe: a
    /// surface absorbs a fraction of an irradiance over an area, and the watts
    /// that result heat a mass with a known specific heat.
    #[test]
    fn absorbed_light_becomes_a_temperature_rise() {
        let irradiance = Irradiance::mw_per_cm2(50.0); // 500 W/m²
        let area: Area = Length::mm(10.0) * Length::mm(10.0); // 1e-4 m²
        let absorptance = 0.02; // what SurfaceOptics::absorptance returns
        let absorbed: Power = irradiance * area * absorptance;
        assert!((absorbed.to_si() - 0.001).abs() < 1e-12, "{absorbed:?}");

        // 1 mW into a 2 g piece of glass for 1 s.
        let glass = Mass::g(2.0);
        let c_p = SpecificHeat::j_per_kg_k(858.0); // N-BK7
        let capacity: HeatCapacity = glass * c_p;
        let heat: Energy = absorbed * Time::s(1.0);
        let rise: Temperature = heat / capacity;
        assert!(
            (rise.to_si() - 0.000_582_7).abs() < 1e-7,
            "expected about 0.58 mK, got {rise:?}"
        );
    }

    /// A photon at 550 nm carries 3.6e-19 J, and the count per watt follows.
    /// This is the number that separates radiometry from photon counting.
    #[test]
    fn photon_energy_matches_the_textbook_figure() {
        let e = photon_energy(Length::nm(550.0));
        assert!((e.to_si() - 3.612e-19).abs() < 1e-21, "{e:?}");
        // 2.26 eV, and about 2.77e18 photons in a joule.
        let per_joule = Energy::j(1.0) / e;
        assert!((per_joule - 2.768e18).abs() < 1e15, "got {per_joule:e}");
    }

    /// Debug prints the dimension, so a mismatch found at a boundary can be
    /// reported in a form a human recognises.
    #[test]
    fn debug_shows_the_dimension() {
        assert_eq!(format!("{:?}", Force::from_si(6.0)), "6·m·kg·s^-2");
        assert_eq!(format!("{:?}", Dimensionless::ratio(0.5)), "0.5");
        assert_eq!(Force::dimension(), [1, 1, -2, 0, 0, 0, 0]);
    }

    /// Serialised as the bare SI number: the dimension is in the field's type,
    /// not repeated in every scene file.
    #[test]
    fn serialises_as_a_bare_si_number() {
        let json = serde_json::to_string(&Length::mm(25.4)).unwrap();
        assert_eq!(json, "0.0254");
        let back: Length = serde_json::from_str(&json).unwrap();
        assert_eq!(back, Length::mm(25.4));
    }
}
