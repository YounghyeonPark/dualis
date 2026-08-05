//! How a scene changes with time, where it changes in closed form.
//!
//! Everything here is a function of a time: ask for the world at `t` and you get
//! it, with no state carried between calls. That is what lets a video frame be
//! sampled several times across its own exposure — which is where motion blur
//! comes from — and it keeps a recording reproducible, since frame 7 does not
//! depend on having rendered frame 6.
//!
//! A scene as authored is the scene at `t = 0`.
//!
//! # Where this stops
//!
//! Drift, oscillation and spin have closed forms. Three bodies under gravity do
//! not, and neither do contact, heat or a stiff reaction. Those go through
//! [`Integrator`](crate::integrator::Integrator), which gives up the
//! frame-independence in exchange for being able to express them at all. This
//! module is the fast path, not the general one — and a domain built on it can
//! declare itself [`Kind::QuasiStatic`](crate::sim::Kind::QuasiStatic), because
//! there is no state for a scheduler to march.

use dualis_units::{Frequency, LengthVec, Time, VelocityVec};
use glam::{DQuat, DVec3};
use serde::{Deserialize, Serialize};

/// Rigid motion of one element, assembly or source.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Motion {
    /// Constant velocity: a conveyor belt, a translation stage, or the drift of a
    /// star field across an untracked telescope.
    Drift { velocity: VelocityVec },
    /// Sinusoidal displacement about the authored position, which is the centre of
    /// the swing. `amplitude` carries both size and direction, so a floor
    /// vibration and an axial focus dither are the same variant.
    Oscillate {
        amplitude: LengthVec,
        frequency: Frequency,
        #[serde(default)]
        phase_deg: f64,
    },
    /// Rotation at `rate_deg_per_s` about `axis` through `pivot`. A rotating stage
    /// or a scanning mirror: it turns what it carries as well as moving it, so an
    /// element's own axis tilts too.
    Spin {
        axis: DVec3,
        #[serde(default)]
        pivot: LengthVec,
        rate_deg_per_s: f64,
    },
}

impl Motion {
    /// The rotation this motion has accumulated by time `t`. Identity for the
    /// purely translational variants.
    pub fn rotation_at(&self, t: Time) -> DQuat {
        match self {
            Motion::Spin {
                axis,
                rate_deg_per_s,
                ..
            } => DQuat::from_axis_angle(
                axis.normalize_or(DVec3::Z),
                (rate_deg_per_s * t.to_si()).to_radians(),
            ),
            _ => DQuat::IDENTITY,
        }
    }

    /// Where a point authored at `p` sits at time `t`.
    pub fn move_point(&self, p: LengthVec, t: Time) -> LengthVec {
        match self {
            Motion::Drift { velocity } => p + *velocity * t,
            Motion::Oscillate {
                amplitude,
                frequency,
                phase_deg,
            } => {
                let phase =
                    std::f64::consts::TAU * frequency.to_si() * t.to_si() + phase_deg.to_radians();
                p + *amplitude * phase.sin()
            }
            Motion::Spin { pivot, .. } => {
                let offset = (p - *pivot).to_si();
                *pivot + LengthVec::from_si(self.rotation_at(t) * offset)
            }
        }
    }

    /// Where a direction authored as `d` points at time `t`. A direction is
    /// dimensionless, and translation does not change it.
    pub fn turn(&self, d: DVec3, t: Time) -> DVec3 {
        match self {
            Motion::Spin { .. } => self.rotation_at(t) * d,
            _ => d,
        }
    }

    /// Velocity of the point authored at `p`, at time `t`.
    ///
    /// Differentiated in closed form rather than by finite difference, so it is
    /// exact — which matters because this is what sets a motion-blur streak length
    /// and what a Doppler shift would be computed from.
    pub fn velocity_at(&self, p: LengthVec, t: Time) -> VelocityVec {
        match self {
            Motion::Drift { velocity } => *velocity,
            Motion::Oscillate {
                amplitude,
                frequency,
                phase_deg,
            } => {
                let omega = std::f64::consts::TAU * frequency.to_si();
                let phase = omega * t.to_si() + phase_deg.to_radians();
                VelocityVec::from_si(amplitude.to_si() * omega * phase.cos())
            }
            Motion::Spin {
                axis,
                pivot,
                rate_deg_per_s,
            } => {
                // v = ω × r, with ω along the axis.
                let omega = axis.normalize_or(DVec3::Z) * rate_deg_per_s.to_radians();
                let r = (self.move_point(p, t) - *pivot).to_si();
                VelocityVec::from_si(omega.cross(r))
            }
        }
    }
}

/// How a light is gated in time.
///
/// This is how machine vision freezes a moving part: fire the light for a small
/// fraction of the frame and the object barely moves while it is lit. The trade is
/// brightness — gating away nine tenths of the time delivers a tenth of the energy
/// — and nothing here hides that. Watch out for auto-exposure, which will happily
/// lengthen the integration to win the light back and undo the whole point; a
/// strobed station wants a fixed exposure.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Strobe {
    /// Time between pulses. Match it to the frame period to fire once per frame.
    pub period: Time,
    /// Fraction of each period the light is on, 0..1.
    pub duty: f64,
    /// Where the pulse sits inside the period.
    #[serde(default)]
    pub phase: Time,
}

impl Strobe {
    pub fn new(period: Time, duty: f64) -> Strobe {
        Strobe {
            period,
            duty,
            phase: Time::ZERO,
        }
    }

    /// Whether the light is on at time `t`.
    pub fn is_on(&self, t: Time) -> bool {
        if self.period.to_si() <= 0.0 {
            return true;
        }
        let duty = self.duty.clamp(0.0, 1.0);
        if duty >= 1.0 {
            return true;
        }
        (((t - self.phase).to_si()) / self.period.to_si()).rem_euclid(1.0) < duty
    }

    /// Total on-time from the first pulse edge up to `x`. Whole periods each
    /// contribute `duty`, and the part-period at the end contributes whatever of
    /// the pulse it reaches.
    fn on_time_to(&self, x: Time) -> f64 {
        let duty = self.duty.clamp(0.0, 1.0);
        let u = (x - self.phase).to_si() / self.period.to_si();
        let whole = u.floor();
        (whole * duty + (u - whole).min(duty)) * self.period.to_si()
    }

    /// Fraction of a window of `length` opening at `start` for which the light is
    /// on. Exact, so it can be checked against the duty cycle — and it is what
    /// says how much of an exposure's light a strobe actually delivers.
    pub fn on_fraction(&self, start: Time, length: Time) -> f64 {
        if self.period.to_si() <= 0.0 || self.duty >= 1.0 {
            return 1.0;
        }
        if length.to_si() <= 0.0 {
            return f64::from(self.is_on(start));
        }
        ((self.on_time_to(start + length) - self.on_time_to(start)) / length.to_si())
            .clamp(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dualis_units::Length;

    /// Drift is exactly velocity x time, and it does not turn anything.
    #[test]
    fn drift_moves_at_its_velocity() {
        let m = Motion::Drift {
            velocity: VelocityVec::mm_per_s(120.0, 0.0, -5.0),
        };
        assert_eq!(m.move_point(LengthVec::ZERO, Time::ZERO), LengthVec::ZERO);
        let after = m.move_point(LengthVec::mm(1.0, 2.0, 3.0), Time::s(0.25));
        assert!((after - LengthVec::mm(31.0, 2.0, 1.75)).length().to_si() < 1e-12);
        assert_eq!(m.turn(DVec3::Z, Time::s(10.0)), DVec3::Z);
        // A drift's velocity is the same everywhere and always.
        assert_eq!(
            m.velocity_at(LengthVec::mm(5.0, 0.0, 0.0), Time::s(3.0)),
            VelocityVec::mm_per_s(120.0, 0.0, -5.0)
        );
    }

    /// A sine wave, checked at the four points where its value is exact.
    #[test]
    fn oscillation_swings_about_where_it_was_authored() {
        let amplitude = LengthVec::mm(0.0, 0.02, 0.0);
        let m = Motion::Oscillate {
            amplitude,
            frequency: Frequency::hz(50.0),
            phase_deg: 0.0,
        };
        let p = LengthVec::mm(1.0, 1.0, 1.0);
        let period = Time::s(1.0 / 50.0);
        let at = |t: Time| m.move_point(p, t);
        assert!((at(Time::ZERO) - p).length().to_si() < 1e-12);
        assert!((at(period / 4.0) - (p + amplitude)).length().to_si() < 1e-12);
        assert!((at(period / 2.0) - p).length().to_si() < 1e-12);
        assert!((at(period * 0.75) - (p - amplitude)).length().to_si() < 1e-12);
        // And it is periodic, which a drift is not.
        assert!((at(period) - p).length().to_si() < 1e-12);
    }

    /// The oscillator's velocity is fastest through the centre and zero at the
    /// extremes, with a peak of exactly `2 pi f A` — a closed form, so this checks
    /// the derivative rather than assuming it.
    #[test]
    fn oscillation_is_fastest_through_the_centre() {
        let amplitude = LengthVec::mm(0.0, 0.02, 0.0);
        let f = 50.0;
        let m = Motion::Oscillate {
            amplitude,
            frequency: Frequency::hz(f),
            phase_deg: 0.0,
        };
        let p = LengthVec::ZERO;
        let period = Time::s(1.0 / f);
        let peak = std::f64::consts::TAU * f * amplitude.length().to_si();
        assert!((m.velocity_at(p, Time::ZERO).length().to_si() - peak).abs() < 1e-12);
        // At the top of the swing it has stopped.
        assert!(m.velocity_at(p, period / 4.0).length().to_si() < 1e-12);
        // 6283 micrometres per second for a 20 micrometre swing at 50 Hz — a small
        // motion moving fast enough to smear a millisecond exposure by six
        // micrometres, which is what decides whether the image is sharp.
        let peak_um_per_s = peak * 1e6;
        assert!((peak_um_per_s - 6283.2).abs() < 0.1, "got {peak_um_per_s}");
    }

    /// A quarter turn about +z takes +x to +y, and keeps every point at its own
    /// radius from the pivot.
    #[test]
    fn a_spin_turns_both_position_and_axis() {
        let pivot = LengthVec::mm(5.0, 0.0, 0.0);
        let m = Motion::Spin {
            axis: DVec3::Z,
            pivot,
            rate_deg_per_s: 90.0,
        };
        let p = LengthVec::mm(7.0, 0.0, 0.0);
        let after = m.move_point(p, Time::s(1.0));
        assert!(
            (after - LengthVec::mm(5.0, 2.0, 0.0)).length().to_si() < 1e-12,
            "got {after:?}"
        );
        assert!(
            ((after - pivot).length() - (p - pivot).length())
                .abs()
                .to_si()
                < 1e-12
        );
        let axis = m.turn(DVec3::X, Time::s(1.0));
        assert!((axis - DVec3::Y).length() < 1e-12, "got {axis}");
    }

    /// Rotational velocity is `omega x r`: perpendicular to the radius, and
    /// proportional to it. A point on the pivot does not move at all.
    #[test]
    fn spin_velocity_grows_with_the_radius() {
        let pivot = LengthVec::ZERO;
        let m = Motion::Spin {
            axis: DVec3::Z,
            pivot,
            rate_deg_per_s: 90.0,
        };
        let omega = 90f64.to_radians();
        for radius_mm in [1.0, 7.0, 100.0] {
            let p = LengthVec::mm(radius_mm, 0.0, 0.0);
            let v = m.velocity_at(p, Time::ZERO);
            let expected = omega * Length::mm(radius_mm).to_si();
            assert!((v.length().to_si() - expected).abs() < 1e-12, "{v:?}");
            // Perpendicular to the radius, as circular motion must be.
            assert!(v.along(p.normalize()).to_si().abs() < 1e-12);
        }
        assert!(m.velocity_at(pivot, Time::ZERO).length().to_si() < 1e-12);
    }

    /// A strobe is on for its duty cycle and nothing more, and averaged over a
    /// whole period the fraction it is on *is* the duty cycle.
    #[test]
    fn a_strobe_is_on_for_its_duty_cycle() {
        let s = Strobe::new(Time::s(1.0 / 60.0), 0.1);
        assert!(s.is_on(Time::ZERO));
        assert!(s.is_on(Time::s(0.9 * 0.1 / 60.0)));
        assert!(!s.is_on(Time::s(0.5 / 60.0)));
        // Periodic: the same instant one period later.
        assert!(s.is_on(Time::s(1.0 / 60.0)));

        for window in [1.0 / 60.0, 3.0 / 60.0, 1.0] {
            let f = s.on_fraction(Time::ZERO, Time::s(window));
            assert!(
                (f - 0.1).abs() < 1e-9,
                "over {window} s the light should be on a tenth of the time, got {f}"
            );
        }
        // A window that opens inside the pulse and closes before the next one.
        let f = s.on_fraction(Time::ZERO, Time::s(0.05 / 60.0));
        assert!((f - 1.0).abs() < 1e-9, "still inside the pulse, got {f}");
        // Duty 1 and a zero period are both "always on".
        assert_eq!(
            Strobe::new(Time::s(1.0), 1.0).on_fraction(Time::s(0.3), Time::s(0.4)),
            1.0
        );
        assert_eq!(
            Strobe::new(Time::ZERO, 0.1).on_fraction(Time::s(0.3), Time::s(0.4)),
            1.0
        );
    }
}
