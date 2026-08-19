//! Ray geometry and the sampling patterns that feed it.
//!
//! Analytic intersections against the shapes optics is actually made of — a
//! spherical or conic cap, a flat annulus, a cylinder — together with Snell's law
//! and the disc samplings that turn an aperture into a bundle of rays. None of it
//! knows what a lens is; this is the arithmetic underneath one.
//!
//! # A ray, a hit, and a dimensionless direction
//!
//! [`Ray`] and [`Hit`] replace what used to be three different tuple shapes —
//! `Option<f64>`, `Option<(f64, DVec3)>` and `Option<(f64, DVec3, DVec3)>` — one
//! per intersection function.
//!
//! A ray's origin is a [`LengthVec`] and its direction is a bare `DVec3`, and that
//! is not an inconsistency: a direction is a length over a length, so it is
//! genuinely dimensionless, and the type therefore refuses to let a displacement
//! be passed as a direction or the reverse.
//!
//! # Metres inside, dimensions at the edge
//!
//! The intersection kernels solve their quadratics in raw SI numbers. Writing
//! `oc.length_squared() - r * r` needs the square of a length, and expressing that
//! in a const generic parameter needs unstable features — so the dimensions are
//! checked at the signature and dropped for the six lines in the middle. That is the
//! one place in this crate where the checking stops, and it is deliberate.

use glam::DVec3;
use pantometry_core::Rng;
use pantometry_units::{Length, LengthVec};

/// Smallest hit distance worth accepting. Anything nearer is the surface the ray
/// just left, found again through floating-point noise.
pub const EPS: Length = Length::from_si(1e-9);

/// A ray: where it starts and which way it goes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ray {
    /// Where it starts.
    pub origin: LengthVec,
    /// Unit direction. Dimensionless, because a direction is.
    pub dir: DVec3,
}

impl Ray {
    /// A ray with a normalised direction.
    pub fn new(origin: LengthVec, dir: DVec3) -> Ray {
        Ray {
            origin,
            dir: dir.normalize(),
        }
    }

    /// The point `t` along the ray.
    pub fn at(&self, t: Length) -> LengthVec {
        self.origin + LengthVec::from_si(self.dir * t.to_si())
    }

    /// The same ray started at `t`, nudged past the surface it is leaving so the
    /// next intersection does not find it again.
    pub fn advance(&self, t: Length) -> Ray {
        Ray {
            origin: self.at(t + EPS),
            dir: self.dir,
        }
    }

    /// The same ray going a new way from where it got to.
    pub fn redirect(&self, t: Length, dir: DVec3) -> Ray {
        Ray {
            origin: self.at(t),
            dir: dir.normalize(),
        }
    }
}

/// Where a ray met a surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    /// Distance along the ray.
    pub t: Length,
    /// Where it landed, in the same frame as the ray.
    pub point: LengthVec,
    /// Geometric unit normal, as the surface stores it — not yet oriented against
    /// the ray. Use [`oriented_against`](pantometry_core::oriented_against) for that.
    pub normal: DVec3,
}

/// Vector-form Snell's law. `n` must oppose `d`; `eta` is n1/n2.
///
/// Returns `None` on total internal reflection, which is not a failure — it is the
/// physics saying there is nowhere for a transmitted ray to go, and the caller
/// should reflect instead.
pub fn refract(d: DVec3, n: DVec3, eta: f64) -> Option<DVec3> {
    let cos_i = -d.dot(n);
    let sin2_t = eta * eta * (1.0 - cos_i * cos_i);
    if sin2_t > 1.0 {
        return None;
    }
    let cos_t = (1.0 - sin2_t).sqrt();
    Some((d * eta + n * (eta * cos_i - cos_t)).normalize())
}

/// Unit-disc hexapolar sample pattern: chief point plus `rings` rings with 6i points
/// on ring i. Deterministic, no RNG (matters for WASM and for tests).
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
/// An extended emitter sampled on a bare hexapolar grid renders as spokes and rings
/// — a picture of the sampling, not of the light. Jittering breaks that up, and
/// `seed` decides how: the same seed always gives the same pattern, so a scene still
/// traces identically every time, while a different seed gives an independent
/// pattern to average against the first.
pub fn jittered_disc(rings: u32, seed: u64) -> Vec<(f64, f64)> {
    if rings == 0 {
        return vec![(0.0, 0.0)];
    }
    let mut rng = Rng::new(seed);
    let dr = 1.0 / rings as f64;
    let mut pts = Vec::new();
    // The chief point wanders over the innermost cell rather than sitting exactly on
    // the axis, which would otherwise be the one unjittered sample.
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

/// Sag of a spherical surface: axial offset from the vertex at radial height `h`.
/// Zero for a flat surface (R = 0).
pub fn sag(r: Length, h: Length) -> Length {
    Length::from_si(conic_sag_si(r.to_si(), 0.0, h.to_si()))
}

/// Sag of a conic surface, with conic constant `k`.
///
/// `k = 0` is a sphere, `k = -1` a paraboloid, `-1 < k < 0` an ellipsoid and
/// `k < -1` a hyperboloid. The distinction is not academic: a paraboloid focuses an
/// on-axis collimated beam with no spherical aberration at all, which is why a
/// collector mirror is one and a sphere is a compromise.
pub fn conic_sag(r: Length, k: f64, h: Length) -> Length {
    Length::from_si(conic_sag_si(r.to_si(), k, h.to_si()))
}

/// `z = h² / (R (1 + sqrt(1 - (1+k) h²/R²)))`, in metres.
fn conic_sag_si(r: f64, k: f64, h: f64) -> f64 {
    if r == 0.0 {
        return 0.0;
    }
    let h2 = h * h;
    let inner = 1.0 - (1.0 + k) * h2 / (r * r);
    if inner < 0.0 {
        // Past where the conic exists; clamp to its edge rather than returning NaN.
        return r;
    }
    h2 / (r * (1.0 + inner.sqrt()))
}

/// Intersect a ray with a spherical cap (or flat disc when `r` is zero) of
/// semi-aperture `semi_ap`, vertex `v`, unit axis `a`.
pub fn cap_intersect(ray: Ray, v: LengthVec, a: DVec3, r: Length, semi_ap: Length) -> Option<Hit> {
    let (origin, dir) = (ray.origin.to_si(), ray.dir);
    let (v_si, r_si, semi_si) = (v.to_si(), r.to_si(), semi_ap.to_si());
    let eps = EPS.to_si();

    if r_si == 0.0 {
        let t = plane_intersect(ray, v, a)?;
        let p = ray.at(t).to_si();
        let w = p - v_si;
        let rho2 = w.length_squared() - w.dot(a).powi(2);
        if rho2 <= semi_si * semi_si {
            return Some(Hit {
                t,
                point: LengthVec::from_si(p),
                normal: a,
            });
        }
        return None;
    }

    let semi_si = semi_si.min(r_si.abs());
    let center = v_si + a * r_si;
    let oc = origin - center;
    // |oc + t*dir|² = r² with |dir| = 1
    let b = oc.dot(dir);
    let c = oc.length_squared() - r_si * r_si;
    let disc = b * b - c;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    for t in [-b - sq, -b + sq] {
        if t <= eps {
            continue;
        }
        let p = origin + dir * t;
        // Must be on the vertex-side hemisphere: (p - center) . a has the opposite
        // sign of r there, since (v - center) . a = -r.
        if (p - center).dot(a) * r_si.signum() > 1e-12 {
            continue;
        }
        let w = p - v_si;
        let rho2 = w.length_squared() - w.dot(a).powi(2);
        if rho2 > semi_si * semi_si + 1e-12 {
            continue;
        }
        return Some(Hit {
            t: Length::from_si(t),
            point: LengthVec::from_si(p),
            normal: (p - center) / r_si.abs(),
        });
    }
    None
}

/// Intersect a ray with a conic cap, by Newton iteration on the sag.
///
/// A conic has no closed-form ray intersection in general, so this starts from the
/// spherical solution and refines. Fixed iteration count rather than a convergence
/// loop, so the arithmetic path does not depend on the values — the same reason
/// [`Integrator`](pantometry_core::Integrator) refuses adaptive steps.
pub fn conic_intersect(
    ray: Ray,
    v: LengthVec,
    a: DVec3,
    r: Length,
    k: f64,
    semi_ap: Length,
) -> Option<Hit> {
    if k == 0.0 {
        return cap_intersect(ray, v, a, r, semi_ap);
    }
    // Start from the sphere of the same vertex curvature.
    let start = cap_intersect(ray, v, a, r, semi_ap.max(Length::from_si(r.to_si().abs())))?;
    let (origin, dir) = (ray.origin.to_si(), ray.dir);
    let (v_si, r_si, semi_si) = (v.to_si(), r.to_si(), semi_ap.to_si());

    // f(t) = axial position of the ray minus the surface's sag at that height.
    let f = |t: f64| {
        let p = origin + dir * t;
        let w = p - v_si;
        let z = w.dot(a);
        let h = (w.length_squared() - z * z).max(0.0).sqrt();
        z - conic_sag_si(r_si, k, h)
    };

    let mut t = start.t.to_si();
    const ITERATIONS: usize = 6;
    for _ in 0..ITERATIONS {
        let h = 1e-9_f64.max(t.abs() * 1e-9);
        let slope = (f(t + h) - f(t - h)) / (2.0 * h);
        if slope.abs() < 1e-18 {
            break;
        }
        t -= f(t) / slope;
    }
    // Newton can walk off to infinity or land behind the ray's origin; both mean
    // there is no hit rather than a hit at a strange place.
    if !t.is_finite() || t <= EPS.to_si() {
        return None;
    }

    let p = origin + dir * t;
    let w = p - v_si;
    let z = w.dot(a);
    let radial = w - a * z;
    let h = radial.length();
    if h > semi_si + 1e-12 {
        return None;
    }
    // Normal from the surface gradient: dz/dh in the meridional plane.
    let dz_dh = {
        let step = 1e-9;
        (conic_sag_si(r_si, k, h + step) - conic_sag_si(r_si, k, h - step)) / (2.0 * step)
    };
    let radial_dir = if h > 1e-15 { radial / h } else { DVec3::ZERO };
    let normal = (a - radial_dir * dz_dh).normalize();
    Some(Hit {
        t: Length::from_si(t),
        point: LengthVec::from_si(p),
        normal,
    })
}

/// Distance along the ray to a plane through `point` with unit `normal`.
pub fn plane_intersect(ray: Ray, point: LengthVec, normal: DVec3) -> Option<Length> {
    let denom = ray.dir.dot(normal);
    if denom.abs() < 1e-12 {
        return None;
    }
    let t = (point - ray.origin).to_si().dot(normal) / denom;
    (t > EPS.to_si()).then(|| Length::from_si(t))
}

/// Intersect a flat annulus — an aperture stop, a lens edge, a baffle ring.
pub fn annulus_intersect(
    ray: Ray,
    center: LengthVec,
    normal: DVec3,
    r_min: Length,
    r_max: Length,
) -> Option<Hit> {
    let t = plane_intersect(ray, center, normal)?;
    let p = ray.at(t).to_si();
    let w = p - center.to_si();
    let rho2 = w.length_squared() - w.dot(normal).powi(2);
    let (lo, hi) = (r_min.to_si(), r_max.to_si());
    (rho2 >= lo * lo && rho2 <= hi * hi).then(|| Hit {
        t,
        point: LengthVec::from_si(p),
        normal,
    })
}

/// Intersect an open cylinder of the given radius about the axis through `v`,
/// restricted to axial coordinate in `[z0, z1]` measured from `v` along `a`.
///
/// This is a lens barrel, a tube, the inside of a bore — the surfaces stray light
/// actually bounces off on its way to becoming glare.
pub fn cylinder_intersect(
    ray: Ray,
    v: LengthVec,
    a: DVec3,
    radius: Length,
    z0: Length,
    z1: Length,
) -> Option<Hit> {
    let (origin, dir) = (ray.origin.to_si(), ray.dir);
    let (v_si, radius_si) = (v.to_si(), radius.to_si());
    let (lo, hi) = (z0.to_si(), z1.to_si());
    let eps = EPS.to_si();

    let oc = origin - v_si;
    let d_perp = dir - a * dir.dot(a);
    let oc_perp = oc - a * oc.dot(a);
    let qa = d_perp.length_squared();
    if qa < 1e-16 {
        return None; // ray parallel to axis
    }
    let qb = oc_perp.dot(d_perp);
    let qc = oc_perp.length_squared() - radius_si * radius_si;
    let disc = qb * qb - qa * qc;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    for t in [(-qb - sq) / qa, (-qb + sq) / qa] {
        if t <= eps {
            continue;
        }
        let p = origin + dir * t;
        let z = (p - v_si).dot(a);
        if z < lo - 1e-12 || z > hi + 1e-12 {
            continue;
        }
        let n = (p - v_si - a * z).normalize();
        return Some(Hit {
            t: Length::from_si(t),
            point: LengthVec::from_si(p),
            normal: n,
        });
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use pantometry_core::oriented_against;

    fn ray_along(dir: DVec3, from: LengthVec) -> Ray {
        Ray::new(from, dir)
    }

    /// A flat disc is hit where the plane is, and missed outside the aperture.
    #[test]
    fn a_flat_cap_is_a_disc() {
        let ray = ray_along(DVec3::Z, LengthVec::mm(0.0, 0.0, -10.0));
        let hit = cap_intersect(
            ray,
            LengthVec::ZERO,
            DVec3::Z,
            Length::ZERO,
            Length::mm(5.0),
        )
        .expect("straight down the axis");
        assert!((hit.t.in_mm() - 10.0).abs() < 1e-9);
        assert!((hit.point.in_mm() - DVec3::ZERO).length() < 1e-9);
        assert_eq!(hit.normal, DVec3::Z);

        // Outside the semi-aperture, nothing.
        let off = ray_along(DVec3::Z, LengthVec::mm(7.0, 0.0, -10.0));
        assert!(cap_intersect(
            off,
            LengthVec::ZERO,
            DVec3::Z,
            Length::ZERO,
            Length::mm(5.0)
        )
        .is_none());
    }

    /// A spherical cap: an axial ray meets the vertex, and the normal there is the
    /// axis. The sag is what puts the rim behind the vertex.
    #[test]
    fn a_spherical_cap_meets_its_vertex_and_carries_its_sag() {
        let r = Length::mm(50.0);
        let semi = Length::mm(12.5);
        let ray = ray_along(DVec3::Z, LengthVec::mm(0.0, 0.0, -20.0));
        let hit = cap_intersect(ray, LengthVec::ZERO, DVec3::Z, r, semi).unwrap();
        assert!((hit.t.in_mm() - 20.0).abs() < 1e-9, "{:?}", hit.t);
        // The normal points outward from the sphere's centre, which sits at
        // `v + a*r` — so at the vertex of a positive-radius cap it is *against* the
        // axis. Orienting it for a particular ray is `oriented_against`'s job.
        assert!(
            (hit.normal + DVec3::Z).length() < 1e-9,
            "got {}",
            hit.normal
        );
        assert_eq!(oriented_against(hit.normal, ray.dir), -DVec3::Z);

        // Sag at the rim, against the closed form R - sqrt(R^2 - h^2).
        let expected = 50.0 - (50.0f64 * 50.0 - 12.5 * 12.5).sqrt();
        assert!((sag(r, semi).in_mm() - expected).abs() < 1e-9);
        assert_eq!(
            sag(Length::ZERO, semi),
            Length::ZERO,
            "a flat surface has no sag"
        );
    }

    /// A paraboloid is the shape that focuses a collimated beam perfectly, and it
    /// differs from the sphere of the same vertex curvature exactly by its conic
    /// constant. Near the axis they agree; away from it they do not, and the
    /// departure grows as the fourth power of the height.
    #[test]
    fn a_conic_departs_from_the_sphere_away_from_the_vertex() {
        let r = Length::mm(50.0);
        // Near the axis, every conic looks like the sphere: at 1 mm out on a 50 mm
        // radius, the sphere and the paraboloid are a nanometre apart.
        for h_mm in [0.0, 0.1, 1.0] {
            let h = Length::mm(h_mm);
            let gap = (conic_sag(r, -1.0, h) - sag(r, h)).abs();
            assert!(
                gap.in_nm() < 2.0,
                "at {h_mm} mm the conic and sphere differ by {} nm",
                gap.in_nm()
            );
        }
        // At the rim of an f/2 mirror they are 25 micrometres apart — about fifty
        // wavelengths, which is why a fast collector is a paraboloid and a sphere
        // of the same curvature is not a substitute for one.
        let rim = Length::mm(12.5);
        let departure = (conic_sag(r, -1.0, rim) - sag(r, rim)).abs();
        assert!(
            (departure.in_um() - 25.2).abs() < 0.1,
            "departure {} um",
            departure.in_um()
        );
        // The departure is quartic in height to leading order, so doubling the
        // height multiplies it by about sixteen. Only *about*: at h/R = 0.25 the
        // sixth-order term is already contributing a couple of percent, and the
        // measured ratio is 16.4. That excess is the aspheric series refusing to be
        // a single term, which is exactly why a real prescription carries more than
        // one coefficient.
        let half = (conic_sag(r, -1.0, Length::mm(6.25)) - sag(r, Length::mm(6.25))).abs();
        let ratio = departure / half;
        assert!((ratio - 16.4).abs() < 0.2, "ratio {ratio}");
        // A hyperboloid departs the other way from a paraboloid.
        assert!(conic_sag(r, -2.0, rim) < conic_sag(r, -1.0, rim));
        assert!(conic_sag(r, 0.0, rim) > conic_sag(r, -1.0, rim));
    }

    /// The conic intersection agrees with the spherical one when the conic constant
    /// says it is a sphere, which is the only case where an independent answer
    /// exists to check against.
    #[test]
    fn the_conic_solver_reproduces_the_sphere() {
        let r = Length::mm(50.0);
        let semi = Length::mm(12.5);
        for offset_mm in [0.0, 2.0, 6.0, 11.0] {
            let ray = ray_along(DVec3::Z, LengthVec::mm(offset_mm, 0.0, -20.0));
            let sphere = cap_intersect(ray, LengthVec::ZERO, DVec3::Z, r, semi).unwrap();
            let conic = conic_intersect(ray, LengthVec::ZERO, DVec3::Z, r, 0.0, semi).unwrap();
            assert_eq!(sphere, conic, "k = 0 must be the spherical path");

            // And Newton's method on a real conic lands on the surface: the residual
            // between the hit's axial position and the sag at its height is nil.
            let hit = conic_intersect(ray, LengthVec::ZERO, DVec3::Z, r, -1.0, semi).unwrap();
            let w = hit.point.to_si();
            let h = Length::from_si((w.x * w.x + w.y * w.y).sqrt());
            let residual = Length::from_si(w.z) - conic_sag(r, -1.0, h);
            assert!(
                residual.abs().in_nm() < 0.1,
                "at {offset_mm} mm the hit is {} nm off the surface",
                residual.in_nm()
            );
            assert!((hit.normal.length() - 1.0).abs() < 1e-12);
        }
    }

    /// Snell's law, against the closed form: entering glass bends the ray towards
    /// the normal by exactly the ratio of the sines.
    #[test]
    fn refraction_obeys_the_sine_ratio() {
        let n = DVec3::Z;
        let theta_i = 30f64.to_radians();
        let d = DVec3::new(theta_i.sin(), 0.0, -theta_i.cos());
        let inward = oriented_against(n, d);
        let out = refract(d, inward, 1.0 / 1.5).expect("air into glass always transmits");
        let theta_t = out.x.atan2(-out.z);
        // sin(30) = 1.5 sin(theta_t)
        assert!(
            (theta_i.sin() - 1.5 * theta_t.sin()).abs() < 1e-12,
            "theta_t {}",
            theta_t.to_degrees()
        );
        assert!((out.length() - 1.0).abs() < 1e-12);
        // Straight in comes straight out.
        let straight = refract(-n, n, 1.0 / 1.5).unwrap();
        assert!((straight + n).length() < 1e-12);
    }

    /// Past the critical angle there is nowhere for the ray to go, and Snell's law
    /// says so rather than returning a nonsense direction.
    #[test]
    fn refraction_reports_total_internal_reflection() {
        let n = DVec3::Z;
        // 45 degrees inside n = 1.5 is past the 41.8 degree critical angle.
        let theta = 45f64.to_radians();
        let d = DVec3::new(theta.sin(), 0.0, -theta.cos());
        assert!(refract(d, n, 1.5).is_none());
        // 30 degrees is not.
        let shallow = 30f64.to_radians();
        let d = DVec3::new(shallow.sin(), 0.0, -shallow.cos());
        assert!(refract(d, n, 1.5).is_some());
    }

    /// A cylinder is hit on its inside, within its length, and the normal points
    /// away from the axis.
    #[test]
    fn a_cylinder_is_hit_within_its_length() {
        let ray = ray_along(DVec3::X, LengthVec::mm(0.0, 0.0, 5.0));
        let hit = cylinder_intersect(
            ray,
            LengthVec::ZERO,
            DVec3::Z,
            Length::mm(10.0),
            Length::ZERO,
            Length::mm(20.0),
        )
        .expect("straight out to the wall");
        assert!((hit.t.in_mm() - 10.0).abs() < 1e-9);
        assert!((hit.normal - DVec3::X).length() < 1e-9);

        // Beyond the end of the tube, nothing.
        let past = ray_along(DVec3::X, LengthVec::mm(0.0, 0.0, 25.0));
        assert!(cylinder_intersect(
            past,
            LengthVec::ZERO,
            DVec3::Z,
            Length::mm(10.0),
            Length::ZERO,
            Length::mm(20.0)
        )
        .is_none());
        // A ray along the axis never meets the wall.
        let axial = ray_along(DVec3::Z, LengthVec::ZERO);
        assert!(cylinder_intersect(
            axial,
            LengthVec::ZERO,
            DVec3::Z,
            Length::mm(10.0),
            Length::ZERO,
            Length::mm(20.0)
        )
        .is_none());
    }

    /// An annulus catches what is between its radii and nothing else — which is what
    /// makes it an aperture rather than a disc.
    #[test]
    fn an_annulus_catches_only_its_ring() {
        let ring = |x_mm: f64| {
            annulus_intersect(
                ray_along(DVec3::Z, LengthVec::mm(x_mm, 0.0, -5.0)),
                LengthVec::ZERO,
                DVec3::Z,
                Length::mm(4.0),
                Length::mm(8.0),
            )
        };
        assert!(ring(6.0).is_some(), "inside the ring");
        assert!(ring(2.0).is_none(), "through the hole");
        assert!(ring(9.0).is_none(), "outside the rim");
    }

    /// A ray does not find the surface it just left, which is what EPS is for.
    #[test]
    fn advancing_escapes_the_surface_it_just_left() {
        let ray = ray_along(DVec3::Z, LengthVec::mm(0.0, 0.0, -10.0));
        let hit = cap_intersect(
            ray,
            LengthVec::ZERO,
            DVec3::Z,
            Length::ZERO,
            Length::mm(5.0),
        )
        .unwrap();
        let next = ray.advance(hit.t);
        assert!(
            cap_intersect(
                next,
                LengthVec::ZERO,
                DVec3::Z,
                Length::ZERO,
                Length::mm(5.0)
            )
            .is_none(),
            "the same plane must not be hit twice"
        );
        assert!((next.origin - hit.point).length().in_nm() > 0.0);
    }

    /// Hexapolar sampling is exact and deterministic; jittering is deterministic per
    /// seed and different between seeds.
    #[test]
    fn disc_sampling_is_deterministic_and_seed_dependent() {
        // 1 + 6 + 12 + 18 points for three rings.
        assert_eq!(hexapolar_unit(3).len(), 37);
        assert_eq!(hexapolar_unit(3), hexapolar_unit(3));
        assert_eq!(hexapolar_unit(0), vec![(0.0, 0.0)]);
        // Every point is inside the unit disc.
        for (x, y) in hexapolar_unit(4) {
            assert!(x * x + y * y <= 1.0 + 1e-12);
        }

        let a = jittered_disc(3, 0x5A17_7E3D);
        assert_eq!(a.len(), 37);
        assert_eq!(a, jittered_disc(3, 0x5A17_7E3D), "same seed, same pattern");
        assert_ne!(a, jittered_disc(3, 1), "a different seed is independent");
        for (x, y) in a {
            assert!(x * x + y * y <= 1.0 + 1e-12);
        }
    }

    /// A ray's direction is dimensionless and its origin is not, so the two cannot
    /// be confused — and moving along it a length gives a place.
    #[test]
    fn a_ray_carries_a_place_and_a_bare_direction() {
        let ray = Ray::new(LengthVec::mm(1.0, 2.0, 3.0), DVec3::new(0.0, 0.0, 2.0));
        assert!(
            (ray.dir.length() - 1.0).abs() < 1e-15,
            "normalised on build"
        );
        let there = ray.at(Length::mm(10.0));
        assert!((there.in_mm() - DVec3::new(1.0, 2.0, 13.0)).length() < 1e-9);
        let turned = ray.redirect(Length::mm(10.0), DVec3::X);
        assert_eq!(turned.dir, DVec3::X);
        assert_eq!(turned.origin, there);
    }
}
