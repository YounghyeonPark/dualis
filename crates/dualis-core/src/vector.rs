//! Vector arithmetic that no single domain owns.
//!
//! An orthonormal basis about a direction is needed to sample a hemisphere, to
//! orient a scatter lobe and to build a local frame for a contact normal.
//! Reflection about a plane is a light ray off a mirror and a ball off a wall, in
//! the same three lines. Neither belongs to optics, so neither lives there.
//!
//! Refraction does belong to optics — Snell's law is about refractive indices —
//! and is in `dualis-optics` accordingly.

use glam::DVec3;

/// Two unit vectors completing an orthonormal frame with `dir`.
///
/// `dir` must be normalised. The helper axis is chosen by which component of
/// `dir` is small, which keeps the cross product well conditioned — picking a
/// fixed helper would give a near-zero cross product whenever `dir` happened to
/// be parallel to it, and the resulting basis would be numerical noise.
pub fn basis_for(dir: DVec3) -> (DVec3, DVec3) {
    let helper = if dir.x.abs() < 0.9 {
        DVec3::X
    } else {
        DVec3::Y
    };
    let u = helper.cross(dir).normalize();
    let v = dir.cross(u);
    (u, v)
}

/// Flip `n` so it opposes the incident direction (`d · n < 0`).
///
/// A surface has two sides and a stored normal points at one of them. Whether a
/// ray is arriving from the front or the back is a fact about the ray, not the
/// surface, so the normal is oriented per hit rather than per surface.
pub fn oriented_against(n: DVec3, d: DVec3) -> DVec3 {
    if d.dot(n) > 0.0 {
        -n
    } else {
        n
    }
}

/// Mirror `d` about the plane with normal `n`. Both should be normalised.
pub fn reflect(d: DVec3, n: DVec3) -> DVec3 {
    (d - n * (2.0 * d.dot(n))).normalize()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The basis is orthonormal and right-handed for every direction, including
    /// the axis-aligned ones where a badly chosen helper vector would collapse it.
    #[test]
    fn the_basis_is_orthonormal_everywhere() {
        let directions = [
            DVec3::X,
            DVec3::Y,
            DVec3::Z,
            -DVec3::X,
            -DVec3::Z,
            DVec3::new(1.0, 1.0, 1.0).normalize(),
            DVec3::new(0.9999, 0.01, 0.0).normalize(),
            DVec3::new(-0.7, 0.2, 0.68).normalize(),
        ];
        for d in directions {
            let (u, v) = basis_for(d);
            assert!((u.length() - 1.0).abs() < 1e-12, "u not unit for {d}");
            assert!((v.length() - 1.0).abs() < 1e-12, "v not unit for {d}");
            assert!(u.dot(d).abs() < 1e-12, "u not perpendicular to {d}");
            assert!(v.dot(d).abs() < 1e-12, "v not perpendicular to {d}");
            assert!(u.dot(v).abs() < 1e-12, "u and v not perpendicular for {d}");
            // Right-handed: u x v = d.
            assert!((u.cross(v) - d).length() < 1e-12, "left-handed for {d}");
        }
    }

    /// Reflection reverses the normal component and keeps the tangential one, so
    /// the angle in equals the angle out and the reflection of a reflection is the
    /// original.
    #[test]
    fn reflection_preserves_the_angle_and_is_its_own_inverse() {
        let n = DVec3::new(0.0, 1.0, 0.0);
        let d = DVec3::new(1.0, -1.0, 0.0).normalize();
        let r = reflect(d, n);
        assert!((r - DVec3::new(1.0, 1.0, 0.0).normalize()).length() < 1e-12);
        // Equal angles either side.
        assert!((d.dot(n).abs() - r.dot(n).abs()).abs() < 1e-12);
        // Reflecting the reflection comes back, up to direction.
        assert!((reflect(r, n) - d).length() < 1e-12);
        // Normal incidence turns straight around.
        assert!((reflect(-n, n) - n).length() < 1e-12);
    }

    #[test]
    fn a_normal_is_oriented_per_hit_not_per_surface() {
        let n = DVec3::Z;
        // Arriving from the front, the normal is left alone.
        assert_eq!(oriented_against(n, -DVec3::Z), n);
        // Arriving from behind, it is flipped to face the ray.
        assert_eq!(oriented_against(n, DVec3::Z), -n);
        // And either way the result opposes the ray, which is the contract.
        for d in [DVec3::Z, -DVec3::Z, DVec3::new(0.3, 0.1, 0.9).normalize()] {
            assert!(oriented_against(n, d).dot(d) <= 0.0);
        }
    }
}
