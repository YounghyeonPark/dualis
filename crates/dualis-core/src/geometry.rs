//! Ray geometry and the sampling patterns that feed it.
//!
//! Analytic intersections against the shapes optics is actually made of — a
//! spherical cap, a flat annulus, a cylinder — together with Snell's law and the
//! disc samplings that turn an aperture into a bundle of rays. None of it knows
//! what a lens is; this is the arithmetic underneath one.

use glam::DVec3;

use crate::rng::Rng;

/// Smallest hit distance worth accepting, in millimetres. Anything nearer is the
/// surface the ray just left, found again through floating-point noise.
pub const EPS: f64 = 1e-6;
/// Flip `n` so it opposes the incident direction (d . n < 0).
pub fn oriented_against(n: DVec3, d: DVec3) -> DVec3 {
    if d.dot(n) > 0.0 {
        -n
    } else {
        n
    }
}

pub fn reflect(d: DVec3, n: DVec3) -> DVec3 {
    (d - n * (2.0 * d.dot(n))).normalize()
}

/// Vector-form Snell's law. `n` must oppose `d`; `eta` = n1/n2.
/// Returns None on total internal reflection.
pub fn refract(d: DVec3, n: DVec3, eta: f64) -> Option<DVec3> {
    let cos_i = -d.dot(n);
    let sin2_t = eta * eta * (1.0 - cos_i * cos_i);
    if sin2_t > 1.0 {
        return None;
    }
    let cos_t = (1.0 - sin2_t).sqrt();
    Some((d * eta + n * (eta * cos_i - cos_t)).normalize())
}

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

/// Unit-disc hexapolar sample pattern: chief point plus `rings` rings with
/// 6*i points on ring i. Deterministic, no RNG (matters for WASM + tests).
pub fn hexapolar_unit(rings: u32) -> Vec<(f64, f64)> {
    let mut pts = vec![(0.0, 0.0)];
    for i in 1..=rings {
        let r = i as f64 / rings.max(1) as f64;
        let n = 6 * i;
        for k in 0..n {
            let phi = std::f64::consts::TAU * k as f64 / n as f64;
            pts.push((r * phi.cos(), r * phi.sin()));
        }
    }
    pts
}

/// As [`hexapolar_unit`], but with each sample jittered inside its own cell.
///
/// An extended emitter sampled on a bare hexapolar grid renders as spokes and
/// rings â€” a picture of the sampling, not of the light. Jittering breaks that
/// up; the seed is fixed, so a scene still traces identically every time.
pub fn jittered_disc(rings: u32) -> Vec<(f64, f64)> {
    if rings == 0 {
        return vec![(0.0, 0.0)];
    }
    let mut rng = Rng::new(0x5A17_7E3D);
    let dr = 1.0 / rings as f64;
    let mut pts = Vec::new();
    // The chief point wanders over the innermost cell rather than sitting
    // exactly on the axis, which would otherwise be the one unjittered sample.
    let (cx, cy) = rng.in_disc(dr / 2.0);
    pts.push((cx, cy));
    for i in 1..=rings {
        let n = 6 * i;
        let dphi = std::f64::consts::TAU / n as f64;
        for k in 0..n {
            let r = (i as f64 * dr + rng.range(-0.5, 0.5) * dr).clamp(0.0, 1.0);
            let phi = dphi * (k as f64 + rng.range(-0.5, 0.5));
            pts.push((r * phi.cos(), r * phi.sin()));
        }
    }
    pts
}

/// Surface sag: axial offset (along +axis from the vertex) of the surface at
/// radial height h. Zero for flat surfaces (R = 0).
pub fn sag(r: f64, h: f64) -> f64 {
    if r == 0.0 {
        return 0.0;
    }
    let h = h.min(r.abs());
    r - r.signum() * (r * r - h * h).sqrt()
}

/// Intersect a ray with a spherical cap (or flat disc when r == 0) of
/// semi-aperture `semi_ap`, vertex `v`, axis `a` (unit). Returns the nearest
/// valid (t, point, geometric unit normal).
pub fn cap_intersect(
    origin: DVec3,
    dir: DVec3,
    v: DVec3,
    a: DVec3,
    r: f64,
    semi_ap: f64,
) -> Option<(f64, DVec3, DVec3)> {
    if r == 0.0 {
        let t = plane_intersect(origin, dir, v, a)?;
        let p = origin + dir * t;
        let w = p - v;
        let rho2 = w.length_squared() - w.dot(a).powi(2);
        if rho2 <= semi_ap * semi_ap {
            return Some((t, p, a));
        }
        return None;
    }
    let semi_ap = semi_ap.min(r.abs());
    let center = v + a * r;
    let oc = origin - center;
    // |oc + t*dir|^2 = r^2 with |dir| = 1
    let b = oc.dot(dir);
    let c = oc.length_squared() - r * r;
    let disc = b * b - c;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    for t in [-b - sq, -b + sq] {
        if t <= EPS {
            continue;
        }
        let p = origin + dir * t;
        // Must be on the vertex-side hemisphere: (p - center) . a has the
        // opposite sign of r there ((v - center) . a = -r).
        if (p - center).dot(a) * r.signum() > 1e-9 {
            continue;
        }
        let w = p - v;
        let rho2 = w.length_squared() - w.dot(a).powi(2);
        if rho2 > semi_ap * semi_ap + 1e-9 {
            continue;
        }
        return Some((t, p, (p - center) / r.abs()));
    }
    None
}

pub fn plane_intersect(origin: DVec3, dir: DVec3, point: DVec3, normal: DVec3) -> Option<f64> {
    let denom = dir.dot(normal);
    if denom.abs() < 1e-12 {
        return None;
    }
    let t = (point - origin).dot(normal) / denom;
    (t > EPS).then_some(t)
}

pub fn annulus_intersect(
    origin: DVec3,
    dir: DVec3,
    center: DVec3,
    normal: DVec3,
    r_min: f64,
    r_max: f64,
) -> Option<(f64, DVec3)> {
    let t = plane_intersect(origin, dir, center, normal)?;
    let p = origin + dir * t;
    let w = p - center;
    let rho2 = w.length_squared() - w.dot(normal).powi(2);
    (rho2 >= r_min * r_min && rho2 <= r_max * r_max).then_some((t, p))
}

/// Intersect with an open cylinder of radius `radius` around the element
/// axis, restricted to axial coordinate (measured from `v` along `a`) in
/// [z0, z1]. Returns nearest (t, point, outward normal).
pub fn cylinder_intersect(
    origin: DVec3,
    dir: DVec3,
    v: DVec3,
    a: DVec3,
    radius: f64,
    z0: f64,
    z1: f64,
) -> Option<(f64, DVec3, DVec3)> {
    let oc = origin - v;
    let d_perp = dir - a * dir.dot(a);
    let oc_perp = oc - a * oc.dot(a);
    let qa = d_perp.length_squared();
    if qa < 1e-16 {
        return None; // ray parallel to axis
    }
    let qb = oc_perp.dot(d_perp);
    let qc = oc_perp.length_squared() - radius * radius;
    let disc = qb * qb - qa * qc;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    for t in [(-qb - sq) / qa, (-qb + sq) / qa] {
        if t <= EPS {
            continue;
        }
        let p = origin + dir * t;
        let z = (p - v).dot(a);
        if z < z0 - 1e-9 || z > z1 + 1e-9 {
            continue;
        }
        let n = (p - v - a * z).normalize();
        return Some((t, p, n));
    }
    None
}
