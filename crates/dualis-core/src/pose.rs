//! Where something is, and which way it faces.
//!
//! The kernel had fields — functions of position — and [`Interface`](crate::Interface) for a
//! discretised boundary, and no way to say where either of them *was*. A domain's coordinates
//! were world coordinates, implicitly, so two grids could not be placed against each other at
//! all: not rotated, not offset, not stacked. Almost everything a scene wants sits behind this.
//!
//! # A rigid motion, and nothing more
//!
//! A [`Pose`] is a rotation and a translation. It is deliberately **not** a general transform:
//!
//! - **No scale.** A scaled metre is not a metre. Every quantity in this workspace carries its
//!   dimension in its type precisely so a factor of a thousand appears in exactly one place, and
//!   a transform that could stretch space would put one everywhere — silently, in a matrix.
//! - **No shear, no projection.** Both change lengths and angles, and a conservation law stated
//!   over a sheared volume is a different law.
//!
//! What is left is an isometry, which preserves every distance and angle exactly. That is the
//! only class of placement a physics can be moved by without its physics changing, and it is
//! why this type is this small.
//!
//! # Which placement this is
//!
//! **This is physical placement**: a pose changes what the physics computes. Two solids in
//! contact, a lens at a distance, a grid rotated against its neighbour.
//!
//! The other kind — a position handed to something that *has* no geometry, purely so a viewer
//! can draw it — is deliberately not here. A [`ThermalNetwork`] node has a capacity and not a
//! position, and a conductance is not a distance; giving one a coordinate is a statement about a
//! diagram, not about heat. That belongs to the scene layer, above this one, where the physics
//! cannot reach it. If the two shared a type, a drawing coordinate would eventually arrive in a
//! conductance and nothing would fail loudly.
//!
//! [`ThermalNetwork`]: https://docs.rs/dualis-thermal
//!
//! ```
//! use dualis_core::Pose;
//! use dualis_units::{Length, LengthVec};
//! use glam::{DQuat, DVec3};
//!
//! // A grid whose own origin sits 2 m along x, turned a quarter turn about z.
//! let placed = Pose::new(
//!     LengthVec::m(2.0, 0.0, 0.0),
//!     DQuat::from_rotation_z(std::f64::consts::FRAC_PI_2),
//! );
//!
//! // Its local +x axis points along world +y.
//! let facing = placed.direction_to_world(DVec3::X);
//! assert!((facing - DVec3::Y).length() < 1e-15);
//!
//! // And a point one metre out along that axis lands at (2, 1, 0).
//! let p = placed.point_to_world(LengthVec::m(1.0, 0.0, 0.0));
//! assert!((p.to_si() - DVec3::new(2.0, 1.0, 0.0)).length() < 1e-15);
//! ```

use dualis_units::LengthVec;
use glam::{DQuat, DVec3};

/// A rotation and a translation: where a domain's own coordinates sit in the world.
///
/// Cheap to copy. Composition is [`then`](Pose::then) and reads left to right.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pose {
    translation: LengthVec,
    rotation: DQuat,
}

impl Default for Pose {
    fn default() -> Pose {
        Pose::IDENTITY
    }
}

impl Pose {
    /// At the origin, unrotated. What a domain has until something places it.
    pub const IDENTITY: Pose = Pose {
        translation: LengthVec::ZERO,
        rotation: DQuat::IDENTITY,
    };

    /// A rotation and a translation.
    ///
    /// The rotation is normalised on the way in. A quaternion accumulated by repeated
    /// multiplication drifts off the unit sphere, and an unnormalised one stretches space —
    /// which is exactly the thing this type exists not to do.
    pub fn new(translation: LengthVec, rotation: DQuat) -> Pose {
        Pose {
            translation,
            rotation: rotation.normalize(),
        }
    }

    /// Moved, not turned.
    pub fn at(translation: LengthVec) -> Pose {
        Pose::new(translation, DQuat::IDENTITY)
    }

    /// Turned, not moved.
    pub fn turned(rotation: DQuat) -> Pose {
        Pose::new(LengthVec::ZERO, rotation)
    }

    /// Where the local origin sits in the world.
    pub fn translation(&self) -> LengthVec {
        self.translation
    }

    /// How the local axes are turned. Always a unit quaternion.
    pub fn rotation(&self) -> DQuat {
        self.rotation
    }

    /// A point in local coordinates, in the world.
    pub fn point_to_world(&self, local: LengthVec) -> LengthVec {
        LengthVec::from_si(self.rotation * local.to_si() + self.translation.to_si())
    }

    /// A point in world coordinates, in the local frame.
    ///
    /// The exact inverse of [`point_to_world`](Pose::point_to_world) up to rounding — a rigid
    /// motion has an exact inverse, unlike anything that scales.
    pub fn point_to_local(&self, world: LengthVec) -> LengthVec {
        LengthVec::from_si(
            self.rotation
                .conjugate()
                .mul_vec3(world.to_si() - self.translation.to_si()),
        )
    }

    /// A direction in local coordinates, in the world.
    ///
    /// Rotated and **not** translated, which is the whole difference between a direction and a
    /// point. A surface normal moved by a translation would stop being normal to anything.
    pub fn direction_to_world(&self, local: DVec3) -> DVec3 {
        self.rotation * local
    }

    /// A direction in world coordinates, in the local frame.
    pub fn direction_to_local(&self, world: DVec3) -> DVec3 {
        self.rotation.conjugate().mul_vec3(world)
    }

    /// This pose, then `outer`: the result of placing something by `self` inside a frame that is
    /// itself placed by `outer`.
    ///
    /// Reads left to right, so `a.then(b).then(c)` applies `a` first. Associative, and a test
    /// says so — matrix composition order is the classic place a sign or an order gets lost, and
    /// it is silent when it does.
    pub fn then(&self, outer: Pose) -> Pose {
        Pose {
            translation: LengthVec::from_si(
                outer.rotation * self.translation.to_si() + outer.translation.to_si(),
            ),
            rotation: (outer.rotation * self.rotation).normalize(),
        }
    }

    /// The pose that undoes this one.
    pub fn inverse(&self) -> Pose {
        let r = self.rotation.conjugate();
        Pose {
            translation: LengthVec::from_si(-(r.mul_vec3(self.translation.to_si()))),
            rotation: r,
        }
    }

    /// Whether this pose moves anything, within a tolerance.
    ///
    /// For a caller that can take a fast path when nothing is placed — sampling a field through
    /// an identity pose should cost nothing.
    pub fn is_identity(&self, tol: f64) -> bool {
        self.translation.to_si().length() <= tol && (self.rotation.w.abs() - 1.0).abs() <= tol
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dualis_units::Length;
    use std::f64::consts::{FRAC_PI_2, PI};

    fn wild() -> Pose {
        Pose::new(
            LengthVec::m(1.5, -0.25, 3.0),
            DQuat::from_euler(glam::EulerRot::XYZ, 0.3, -1.1, 2.2),
        )
    }

    /// **An isometry preserves every distance, which is the property that lets physics be moved.**
    ///
    /// If a pose could change a length, then a conservation law stated in local coordinates
    /// would not be the same law in world coordinates, and placing a domain would silently
    /// change its physics. That is the reason this type has no scale and the reason this is the
    /// first test.
    #[test]
    fn placing_something_cannot_change_a_distance() {
        let p = wild();
        let pairs = [
            (LengthVec::m(0.0, 0.0, 0.0), LengthVec::m(1.0, 0.0, 0.0)),
            (LengthVec::m(-2.0, 5.0, 0.5), LengthVec::m(3.0, -1.0, 4.0)),
            (LengthVec::m(1e-6, 0.0, 0.0), LengthVec::m(-1e-6, 0.0, 0.0)),
        ];
        for (a, b) in pairs {
            let before = (a.to_si() - b.to_si()).length();
            let after = (p.point_to_world(a).to_si() - p.point_to_world(b).to_si()).length();
            assert!(
                (after - before).abs() <= 1e-15 * before.max(1.0),
                "a {before} m separation became {after} m"
            );
        }
    }

    /// **A round trip returns what went in**, which is what "rigid" buys and what a scaling
    /// transform could not offer.
    #[test]
    fn local_to_world_and_back_is_the_identity() {
        let p = wild();
        for v in [
            LengthVec::m(0.0, 0.0, 0.0),
            LengthVec::m(1.0, 2.0, 3.0),
            LengthVec::m(-7.5, 0.0, 1e3),
        ] {
            let back = p.point_to_local(p.point_to_world(v));
            assert!(
                (back.to_si() - v.to_si()).length() < 1e-12,
                "{:?} came back as {:?}",
                v.to_si(),
                back.to_si()
            );
        }
        // And for directions, which do not translate.
        for d in [
            DVec3::X,
            DVec3::Y,
            DVec3::Z,
            DVec3::new(1.0, -2.0, 0.5).normalize(),
        ] {
            let back = p.direction_to_local(p.direction_to_world(d));
            assert!((back - d).length() < 1e-12);
        }
    }

    /// **A direction is rotated and not translated.**
    ///
    /// The one-line difference between a point and a direction, and the bug that follows from
    /// missing it is a surface normal that stops being normal to its surface as soon as anything
    /// is moved off the origin — while every length still checks out.
    #[test]
    fn a_direction_ignores_the_translation() {
        let far = Pose::at(LengthVec::m(1e4, -1e4, 1e4));
        for d in [DVec3::X, DVec3::Y, DVec3::Z] {
            assert_eq!(
                far.direction_to_world(d),
                d,
                "translation rotated a direction"
            );
        }
        // A normal stays perpendicular to the surface it came from, wherever the surface goes.
        let p = wild();
        let (u, v) = (DVec3::X, DVec3::Y);
        let n = u.cross(v);
        let (u2, v2, n2) = (
            p.direction_to_world(u),
            p.direction_to_world(v),
            p.direction_to_world(n),
        );
        assert!(n2.dot(u2).abs() < 1e-15 && n2.dot(v2).abs() < 1e-15);
        assert!((u2.cross(v2) - n2).length() < 1e-15, "handedness flipped");
    }

    /// **Composition is associative, and reads left to right.**
    ///
    /// Order and handedness are where a transform quietly goes wrong: the result is still a
    /// valid pose, still preserves lengths, and puts everything in the wrong place. Checked
    /// against a case whose answer is known by hand as well as against associativity.
    #[test]
    fn composing_is_associative_and_in_the_order_it_reads() {
        let (a, b, c) = (
            Pose::at(LengthVec::m(1.0, 0.0, 0.0)),
            Pose::turned(DQuat::from_rotation_z(FRAC_PI_2)),
            Pose::at(LengthVec::m(0.0, 0.0, 2.0)),
        );
        let left = a.then(b).then(c);
        let right = a.then(b.then(c));
        assert!((left.translation.to_si() - right.translation.to_si()).length() < 1e-14);
        assert!(left.rotation.abs_diff_eq(right.rotation, 1e-14));

        // By hand: translate 1 m along x, then turn a quarter turn about z — the point lands on
        // +y, not +x. If `then` had composed the other way it would still be at +x.
        let origin = a.then(b).point_to_world(LengthVec::ZERO);
        assert!(
            (origin.to_si() - DVec3::new(0.0, 1.0, 0.0)).length() < 1e-14,
            "a then b put the origin at {:?}",
            origin.to_si()
        );
        // And the other order does something different, which is what makes the check mean
        // anything.
        let swapped = b.then(a).point_to_world(LengthVec::ZERO);
        assert!((swapped.to_si() - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-14);
    }

    /// **A pose and its inverse cancel exactly enough**, which a scaling transform could not
    /// promise.
    #[test]
    fn a_pose_undoes_itself() {
        let p = wild();
        let both = p.then(p.inverse());
        assert!(both.is_identity(1e-12), "{both:?}");
        assert!(p.inverse().then(p).is_identity(1e-12));

        // The inverse of the identity is the identity, and a half turn is its own inverse.
        assert!(Pose::IDENTITY.inverse().is_identity(0.0));
        let half = Pose::turned(DQuat::from_rotation_y(PI));
        assert!(half.then(half).is_identity(1e-14));
    }

    /// The identity does nothing, and says it does nothing.
    #[test]
    fn the_identity_is_free() {
        let v = LengthVec::m(3.0, -1.0, 0.25);
        assert_eq!(Pose::IDENTITY.point_to_world(v).to_si(), v.to_si());
        assert_eq!(Pose::default(), Pose::IDENTITY);
        assert!(Pose::IDENTITY.is_identity(0.0));
        assert!(!Pose::at(LengthVec::m(1e-3, 0.0, 0.0)).is_identity(1e-9));
        // A metre is a metre: no constructor here can change one.
        assert_eq!(
            Pose::at(LengthVec::m(5.0, 0.0, 0.0))
                .point_to_world(LengthVec::m(1.0, 0.0, 0.0))
                .to_si()
                .x,
            Length::m(6.0).to_si()
        );
    }
}
