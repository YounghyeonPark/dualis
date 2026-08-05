//! How a scene changes with time.
//!
//! Everything here is a closed-form function of a time in seconds: ask for the
//! world at `t` and you get it, with no state carried between calls. That is
//! what lets a video frame be sampled several times across its own exposure —
//! which is where motion blur comes from — and it keeps a recording
//! reproducible, since frame 7 does not depend on having rendered frame 6.
//!
//! A scene as authored is the scene at `t = 0`.

use glam::{DQuat, DVec3};
use serde::{Deserialize, Serialize};

/// Rigid motion of one element, assembly or source.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Motion {
    /// Constant velocity, mm/s: a conveyor belt, a translation stage, or the
    /// drift of a star field across an untracked telescope.
    Drift { velocity: DVec3 },
    /// Sinusoidal displacement about the authored position, which is the centre
    /// of the swing. `amplitude` carries both size (mm) and direction, so a
    /// floor vibration and an axial focus dither are the same variant.
    Oscillate {
        amplitude: DVec3,
        frequency_hz: f64,
        #[serde(default)]
        phase_deg: f64,
    },
    /// Rotation at `rate_deg_per_s` about `axis` through `pivot`. A rotating
    /// stage or a scanning mirror: it turns what it carries as well as moving
    /// it, so an element's own axis tilts too.
    Spin {
        axis: DVec3,
        #[serde(default)]
        pivot: DVec3,
        rate_deg_per_s: f64,
    },
}

impl Motion {
    /// The rotation this motion has accumulated by time `t`. Identity for the
    /// purely translational variants.
    pub fn rotation_at(&self, t: f64) -> DQuat {
        match self {
            Motion::Spin {
                axis,
                rate_deg_per_s,
                ..
            } => DQuat::from_axis_angle(
                axis.normalize_or(DVec3::Z),
                (rate_deg_per_s * t).to_radians(),
            ),
            _ => DQuat::IDENTITY,
        }
    }

    /// Where a point authored at `p` sits at time `t`.
    pub fn move_point(&self, p: DVec3, t: f64) -> DVec3 {
        match self {
            Motion::Drift { velocity } => p + *velocity * t,
            Motion::Oscillate {
                amplitude,
                frequency_hz,
                phase_deg,
            } => {
                let phase = std::f64::consts::TAU * frequency_hz * t + phase_deg.to_radians();
                p + *amplitude * phase.sin()
            }
            Motion::Spin { pivot, .. } => *pivot + self.rotation_at(t) * (p - *pivot),
        }
    }

    /// Where a direction authored as `d` points at time `t`.
    pub fn turn(&self, d: DVec3, t: f64) -> DVec3 {
        match self {
            Motion::Spin { .. } => self.rotation_at(t) * d,
            _ => d,
        }
    }
}

/// How a light is gated in time.
///
/// This is how machine vision freezes a moving part: fire the light for a small
/// fraction of the frame and the object barely moves while it is lit. The trade
/// is brightness — gating away nine tenths of the time delivers a tenth of the
/// energy — and nothing here hides that. Watch out for auto-exposure, which
/// will happily lengthen the integration to win the light back and undo the
/// whole point; a strobed station wants a fixed exposure.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Strobe {
    /// Seconds between pulses. Match it to the frame period to fire once per
    /// frame.
    pub period_s: f64,
    /// Fraction of each period the light is on, 0..1.
    pub duty: f64,
    /// Where the pulse sits inside the period, in seconds.
    #[serde(default)]
    pub phase_s: f64,
}

impl Strobe {
    /// Whether the light is on at time `t`.
    pub fn is_on(&self, t: f64) -> bool {
        if self.period_s <= 0.0 {
            return true;
        }
        let duty = self.duty.clamp(0.0, 1.0);
        if duty >= 1.0 {
            return true;
        }
        ((t - self.phase_s) / self.period_s).rem_euclid(1.0) < duty
    }

    /// Total on-time from the first pulse edge up to `x`. Whole periods each
    /// contribute `duty`, and the part-period at the end contributes whatever of
    /// the pulse it reaches.
    fn on_time_to(&self, x: f64) -> f64 {
        let duty = self.duty.clamp(0.0, 1.0);
        let u = (x - self.phase_s) / self.period_s;
        let whole = u.floor();
        (whole * duty + (u - whole).min(duty)) * self.period_s
    }

    /// Fraction of a window of `length` seconds opening at `start` for which the
    /// light is on. Exact, so it can be checked against the duty cycle — and it
    /// is what says how much of an exposure's light a strobe actually delivers.
    pub fn on_fraction(&self, start: f64, length: f64) -> f64 {
        if self.period_s <= 0.0 || self.duty >= 1.0 {
            return 1.0;
        }
        if length <= 0.0 {
            return f64::from(self.is_on(start));
        }
        ((self.on_time_to(start + length) - self.on_time_to(start)) / length).clamp(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Drift is exactly velocity x time, and it does not turn anything.
    #[test]
    fn drift_moves_at_its_velocity() {
        let m = Motion::Drift {
            velocity: DVec3::new(120.0, 0.0, -5.0),
        };
        assert_eq!(m.move_point(DVec3::ZERO, 0.0), DVec3::ZERO);
        let after = m.move_point(DVec3::new(1.0, 2.0, 3.0), 0.25);
        assert!((after - DVec3::new(31.0, 2.0, 1.75)).length() < 1e-12);
        assert_eq!(m.turn(DVec3::Z, 10.0), DVec3::Z);
    }

    /// A sine wave, checked at the four points where its value is exact.
    #[test]
    fn oscillation_swings_about_where_it_was_authored() {
        let amplitude = DVec3::new(0.0, 0.02, 0.0);
        let m = Motion::Oscillate {
            amplitude,
            frequency_hz: 50.0,
            phase_deg: 0.0,
        };
        let p = DVec3::new(1.0, 1.0, 1.0);
        let period = 1.0 / 50.0;
        assert!((m.move_point(p, 0.0) - p).length() < 1e-12);
        assert!((m.move_point(p, period / 4.0) - (p + amplitude)).length() < 1e-12);
        assert!((m.move_point(p, period / 2.0) - p).length() < 1e-12);
        assert!((m.move_point(p, 3.0 * period / 4.0) - (p - amplitude)).length() < 1e-12);
        // And it is periodic, which a drift is not.
        assert!((m.move_point(p, period) - p).length() < 1e-12);
    }

    /// A quarter turn about +z takes +x to +y, and keeps every point at its own
    /// radius from the pivot.
    #[test]
    fn a_spin_turns_both_position_and_axis() {
        let pivot = DVec3::new(5.0, 0.0, 0.0);
        let m = Motion::Spin {
            axis: DVec3::Z,
            pivot,
            rate_deg_per_s: 90.0,
        };
        let p = DVec3::new(7.0, 0.0, 0.0);
        let after = m.move_point(p, 1.0);
        assert!(
            (after - DVec3::new(5.0, 2.0, 0.0)).length() < 1e-12,
            "got {after}"
        );
        assert!(((after - pivot).length() - (p - pivot).length()).abs() < 1e-12);
        let axis = m.turn(DVec3::X, 1.0);
        assert!((axis - DVec3::Y).length() < 1e-12, "got {axis}");
    }

    /// A strobe is on for its duty cycle and nothing more, and averaged over a
    /// whole period the fraction it is on *is* the duty cycle.
    #[test]
    fn a_strobe_is_on_for_its_duty_cycle() {
        let s = Strobe {
            period_s: 1.0 / 60.0,
            duty: 0.1,
            phase_s: 0.0,
        };
        assert!(s.is_on(0.0));
        assert!(s.is_on(0.9 * 0.1 / 60.0));
        assert!(!s.is_on(0.5 / 60.0));
        // Periodic: the same instant one period later.
        assert!(s.is_on(1.0 / 60.0));

        for window in [1.0 / 60.0, 3.0 / 60.0, 1.0] {
            let f = s.on_fraction(0.0, window);
            assert!(
                (f - 0.1).abs() < 1e-9,
                "over {window} s the light should be on a tenth of the time, got {f}"
            );
        }
        // A window that opens inside the pulse and closes before the next one.
        let f = s.on_fraction(0.0, 0.05 / 60.0);
        assert!((f - 1.0).abs() < 1e-9, "still inside the pulse, got {f}");
        // Duty 1 and a zero period are both "always on".
        assert_eq!(
            Strobe {
                period_s: 1.0,
                duty: 1.0,
                phase_s: 0.0
            }
            .on_fraction(0.3, 0.4),
            1.0
        );
    }
}
