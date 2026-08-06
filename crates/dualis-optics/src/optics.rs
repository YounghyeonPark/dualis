//! What a surface does to light: reflect it, transmit it, or absorb it.
//!
//! Every surface splits incident power three ways, and the three add to one. That
//! is not a modelling choice, it is conservation of energy, so this module stores
//! only two of them — [`SurfaceOptics::reflectance`] and
//! [`SurfaceOptics::transmittance`] — and computes absorptance as the remainder.
//! There is no way to write down a surface that gains energy.
//!
//! Glass is the exception, and an important one: its reflectance is *not* a free
//! parameter. It follows from the refractive indices either side and the angle of
//! incidence, by Fresnel's equations. A lens surface therefore carries a
//! [`SurfaceFinish`] — bare, or coated — rather than a reflectance curve, and the
//! number comes out of the physics.
//!
//! # Where the absorbed light goes
//!
//! [`SurfaceOptics::absorptance`] says what fraction of the light did not leave.
//! That energy is not gone: it is heat, and it is the seam where a thermal domain
//! attaches — publish it on the [`Exchange`](dualis_core::Exchange) and something
//! else can receive it. Until then this module computes the number and hands it
//! back, which is the honest thing for an optics crate to do.

use dualis_core::{basis_for, oriented_against, reflect, Rng, Violation};
use dualis_units::Length;
use glam::DVec3;
use serde::{Deserialize, Serialize};

use crate::spectrum::{Spectrum, VISIBLE_RANGE_NM};

/// Spectral reflectance and transmittance of a surface. Absorptance is whatever is
/// left, so the three always sum to one.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SurfaceOptics {
    /// Fraction sent back, per wavelength. `diffuse` decides how much of it leaves
    /// specularly and how much scatters.
    pub reflectance: Spectrum,
    /// Fraction let through, per wavelength.
    pub transmittance: Spectrum,
    /// How much of the reflected light leaves diffusely rather than as a mirror
    /// image, 0..1.
    ///
    /// This does not change the energy budget — the light was already counted in
    /// `reflectance` — only where it goes. And where it goes matters enormously: a
    /// polished surface sends its reflection somewhere definite, where you can
    /// design a baffle to catch it, while a matte one sprays it over the whole
    /// hemisphere. That spray is **veiling glare**, the haze that lifts the black
    /// level of an otherwise good image, and it is why the inside of a lens barrel
    /// is flocked rather than merely painted.
    #[serde(default)]
    pub diffuse: f64,
}

/// What happened to a photon at a surface, and which way it left.
///
/// The complement of [`SurfaceOptics::split`]: that says how the energy divides,
/// this says where each share goes. Both are needed, and until now only the first
/// existed — `diffuse` described veiling glare without ever producing a direction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Scatter {
    /// Left as a mirror image, in a direction a baffle can be designed around.
    Specular(DVec3),
    /// Left into the hemisphere, cosine-weighted. This is the glare.
    Diffuse(DVec3),
    /// Went through, in the same direction — refraction is
    /// [`refract`](crate::geometry::refract)'s job, since it needs the indices.
    Transmitted(DVec3),
    /// Became heat.
    Absorbed,
}

impl SurfaceOptics {
    /// Absorptance at a wavelength: everything neither reflected nor transmitted.
    /// Clamped at zero, so an over-specified surface loses the excess rather than
    /// manufacturing light.
    pub fn absorptance(&self, wavelength: Length) -> f64 {
        let (r, t) = self.split(wavelength);
        (1.0 - r - t).max(0.0)
    }

    /// Reflectance and transmittance at a wavelength, rescaled if they sum past one
    /// so that no surface ever returns more light than reached it.
    pub fn split(&self, wavelength: Length) -> (f64, f64) {
        self.split_nm(wavelength.in_nm())
    }

    fn split_nm(&self, wavelength_nm: f64) -> (f64, f64) {
        let r = self.reflectance.at_nm(wavelength_nm).clamp(0.0, 1.0);
        let t = self.transmittance.at_nm(wavelength_nm).clamp(0.0, 1.0);
        let total = r + t;
        if total > 1.0 {
            (r / total, t / total)
        } else {
            (r, t)
        }
    }

    /// Draw one outgoing event, with the weight to carry along it.
    ///
    /// The three shares are sampled in proportion, so the returned weight is 1 for
    /// every branch and the estimator needs no correction — a photon is not split
    /// into three fractional photons, it goes one way. Absorbed returns a weight of
    /// zero and the energy it took, which is what a thermal coupling wants.
    ///
    /// `incident` points *at* the surface and `normal` is the geometric normal;
    /// which side the ray arrived from is worked out here.
    pub fn sample(
        &self,
        incident: DVec3,
        normal: DVec3,
        wavelength: Length,
        rng: &mut Rng,
    ) -> (Scatter, f64) {
        let (r, t) = self.split(wavelength);
        let n = oriented_against(normal, incident);
        let u = rng.unit();
        if u < r {
            let diffuse = self.diffuse.clamp(0.0, 1.0);
            // Whether this particular reflected photon scatters is its own draw:
            // `diffuse` is the fraction that do, not a blur applied to all of them.
            if diffuse > 0.0 && rng.unit() < diffuse {
                (Scatter::Diffuse(rng.cosine_hemisphere(n)), 1.0)
            } else {
                (Scatter::Specular(reflect(incident, n)), 1.0)
            }
        } else if u < r + t {
            (Scatter::Transmitted(incident), 1.0)
        } else {
            (Scatter::Absorbed, 0.0)
        }
    }

    /// Whether reflectance and transmittance stay within the energy budget over the
    /// whole working range. The error names the wavelength that broke it.
    pub fn validate(&self) -> Result<(), Violation> {
        let (lo, hi) = VISIBLE_RANGE_NM;
        const STEPS: usize = 150;
        for i in 0..=STEPS {
            let w = lo + (hi - lo) * i as f64 / STEPS as f64;
            let (r, t) = (self.reflectance.at_nm(w), self.transmittance.at_nm(w));
            if r < 0.0 || t < 0.0 {
                return Err(Violation {
                    quantity: "energy".to_string(),
                    site: format!("surface at {w:.0} nm (negative R or T)"),
                    before: 1.0,
                    after: r + t,
                    // The budget is normalised to one, so that is the scale.
                    scale: 1.0,
                    tolerance: 0.0,
                });
            }
            if r + t > 1.0 + 1e-9 {
                return Err(Violation {
                    quantity: "energy".to_string(),
                    site: format!("surface at {w:.0} nm"),
                    before: 1.0,
                    after: r + t,
                    scale: 1.0,
                    tolerance: 1e-9,
                });
            }
        }
        Ok(())
    }

    /// Perfectly opaque and black: absorbs everything. An aperture blade, a bench,
    /// the inside of a tube.
    pub fn black() -> SurfaceOptics {
        SurfaceOptics {
            reflectance: Spectrum::constant(0.0),
            transmittance: Spectrum::constant(0.0),
            diffuse: 0.0,
        }
    }

    /// Protected aluminium, the ordinary mirror coating: 92% across the visible,
    /// falling off in the near-UV where the protective overcoat absorbs. Four
    /// bounces of that is a third of the light gone, which is why fold count
    /// matters in a real instrument.
    pub fn aluminium() -> SurfaceOptics {
        SurfaceOptics {
            reflectance: Spectrum::curve(vec![
                (350.0, 0.88),
                (400.0, 0.91),
                (500.0, 0.92),
                (650.0, 0.91),
                (900.0, 0.94),
                (1100.0, 0.95),
            ]),
            transmittance: Spectrum::constant(0.0),
            diffuse: 0.0,
        }
    }

    /// Enhanced silver: 98% through the visible, poor in the blue and UV.
    pub fn silver() -> SurfaceOptics {
        SurfaceOptics {
            reflectance: Spectrum::curve(vec![
                (350.0, 0.25),
                (400.0, 0.90),
                (450.0, 0.97),
                (550.0, 0.98),
                (900.0, 0.99),
                (1100.0, 0.99),
            ]),
            transmittance: Spectrum::constant(0.0),
            diffuse: 0.0,
        }
    }

    /// A mirror of uniform reflectance; the rest is absorbed in the coating.
    /// Protected aluminium is about 0.90, enhanced silver 0.98.
    pub fn mirror(reflectance: f64) -> SurfaceOptics {
        SurfaceOptics {
            reflectance: Spectrum::constant(reflectance),
            transmittance: Spectrum::constant(0.0),
            diffuse: 0.0,
        }
    }

    /// A beamsplitter: part through, part back, nothing lost.
    pub fn beamsplitter(reflectance: f64) -> SurfaceOptics {
        SurfaceOptics {
            reflectance: Spectrum::constant(reflectance),
            transmittance: Spectrum::constant(1.0 - reflectance),
            diffuse: 0.0,
        }
    }

    /// An absorbing filter: passes its bands, swallows the rest.
    pub fn filter(bands: Vec<[f64; 2]>, transmission: f64, edge_nm: f64) -> SurfaceOptics {
        SurfaceOptics {
            reflectance: Spectrum::constant(0.0),
            transmittance: Spectrum::interference_bands(bands, transmission, edge_nm),
            diffuse: 0.0,
        }
    }

    /// A dichroic beamsplitter: passes its bands and *reflects* what it blocks, so
    /// the blocked light leaves down another arm instead of being lost. This is what
    /// lets one objective carry excitation down and fluorescence back.
    pub fn dichroic(bands: Vec<[f64; 2]>, transmission: f64, edge_nm: f64) -> SurfaceOptics {
        let pass = Spectrum::interference_bands(bands.clone(), transmission, edge_nm);
        // Reflect the complement: high where the passband is not.
        SurfaceOptics {
            reflectance: Spectrum::Bands {
                bands,
                in_band: 1.0 - transmission,
                out_of_band: 0.98,
                edge_nm,
            },
            transmittance: pass,
            diffuse: 0.0,
        }
    }

    /// A detector face: absorbs nearly everything, and reflects the few percent that
    /// a real sensor's cover glass and silicon send back — which is where sensor
    /// ghosting comes from.
    pub fn detector() -> SurfaceOptics {
        SurfaceOptics {
            reflectance: Spectrum::constant(0.03),
            transmittance: Spectrum::constant(0.0),
            diffuse: 0.0,
        }
    }

    /// Flock paper, or anything else meant to kill a stray reflection: it still
    /// returns a few percent, and it sprays all of it.
    pub fn flocked() -> SurfaceOptics {
        SurfaceOptics {
            reflectance: Spectrum::constant(0.02),
            transmittance: Spectrum::constant(0.0),
            diffuse: 1.0,
        }
    }
}

/// How a glass surface is finished, which is the *only* freedom there is: the bare
/// reflectance itself follows from the refractive indices and the angle.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum SurfaceFinish {
    /// Bare polished glass. About 4% per surface at normal incidence for n = 1.5,
    /// rising steeply past 60°, and the reason an uncoated multi-element lens is dim
    /// and full of ghosts.
    Bare,
    /// An anti-reflection coating. `residual` is the fraction of the bare Fresnel
    /// reflectance left, per wavelength: a good broadband V-coat leaves about 0.05
    /// of it (4% becomes 0.2%) in its design band and degrades outside it, which is
    /// why coatings are specified spectrally.
    Coated {
        /// Fraction of the bare Fresnel reflectance that survives, per wavelength.
        residual: Spectrum,
    },
    /// A mirrored surface: reflect rather than refract, with this reflectance.
    Mirrored {
        /// Reflectance per wavelength. Whatever is not reflected is absorbed, since a mirror
        /// transmits nothing — which is what makes a mirror a heat source.
        reflectance: Spectrum,
    },
}

impl SurfaceFinish {
    /// A broadband AR coating: 0.25% residual reflection across the visible,
    /// worsening in the near-UV and near-IR where the layer stack detunes.
    pub fn broadband_ar() -> SurfaceFinish {
        SurfaceFinish::Coated {
            residual: Spectrum::curve(vec![
                (350.0, 0.60),
                (400.0, 0.15),
                (450.0, 0.06),
                (550.0, 0.04),
                (650.0, 0.06),
                (750.0, 0.15),
                (900.0, 0.45),
                (1100.0, 0.80),
            ]),
        }
    }

    /// Reflectance of this finish at a given index step and angle.
    pub fn reflectance_at(&self, n1: f64, n2: f64, cos_incident: f64, wavelength: Length) -> f64 {
        match self {
            SurfaceFinish::Bare => fresnel_reflectance(n1, n2, cos_incident),
            SurfaceFinish::Coated { residual } => {
                fresnel_reflectance(n1, n2, cos_incident) * residual.at(wavelength).max(0.0)
            }
            SurfaceFinish::Mirrored { reflectance } => reflectance.at(wavelength).clamp(0.0, 1.0),
        }
    }
}

impl Default for SurfaceFinish {
    /// Bare, because that is what a piece of glass is until someone coats it. A
    /// setup that says nothing about coatings gets the honest 4% per surface.
    fn default() -> Self {
        SurfaceFinish::Bare
    }
}

/// Unpolarised Fresnel reflectance at an interface, averaging the s and p
/// polarisations.
///
/// `cos_incident` is |cos θ| in the incident medium. Past the critical angle this
/// returns 1: total internal reflection, which is Fresnel's equations telling you
/// there is nowhere for the transmitted ray to go.
pub fn fresnel_reflectance(n1: f64, n2: f64, cos_incident: f64) -> f64 {
    let (rs, rp) = fresnel_split(n1, n2, cos_incident);
    ((rs + rp) / 2.0).clamp(0.0, 1.0)
}

/// The two polarisations separately, `(rs, rp)`.
///
/// Averaging them is what unpolarised light does, and [`fresnel_reflectance`] does
/// exactly that — but they are not equal, and the difference is the whole of
/// polarisation optics. `rp` falls to zero at Brewster's angle, which is why a
/// polarising filter can kill a reflection off water and why a laser cavity
/// window is tilted.
pub fn fresnel_split(n1: f64, n2: f64, cos_incident: f64) -> (f64, f64) {
    let cos_i = cos_incident.clamp(0.0, 1.0);
    let sin_i2 = 1.0 - cos_i * cos_i;
    let ratio = n1 / n2;
    let sin_t2 = ratio * ratio * sin_i2;
    if sin_t2 >= 1.0 {
        return (1.0, 1.0); // total internal reflection
    }
    let cos_t = (1.0 - sin_t2).sqrt();
    let rs = ((n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)).powi(2);
    let rp = ((n1 * cos_t - n2 * cos_i) / (n1 * cos_t + n2 * cos_i)).powi(2);
    (rs.clamp(0.0, 1.0), rp.clamp(0.0, 1.0))
}

/// Brewster's angle for an index step, in radians: where `rp` vanishes and the
/// reflection is completely s-polarised.
pub fn brewster_angle(n1: f64, n2: f64) -> f64 {
    (n2 / n1).atan()
}

/// Critical angle for total internal reflection, in radians, or `None` when going
/// into a denser medium — where there is no such angle.
pub fn critical_angle(n1: f64, n2: f64) -> Option<f64> {
    (n2 < n1).then(|| (n2 / n1).asin())
}

/// An orthonormal frame in which to express a scatter about `normal`. Re-exported
/// convenience so a consumer does not reach into the kernel for it.
pub fn scatter_frame(normal: DVec3) -> (DVec3, DVec3) {
    basis_for(normal)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nm(v: f64) -> Length {
        Length::nm(v)
    }

    /// Fresnel at normal incidence reduces to ((n2-n1)/(n2+n1))² — 4.0% for n = 1.5
    /// in air, the number every optics text quotes.
    #[test]
    fn fresnel_matches_the_closed_form_at_normal_incidence() {
        for n in [1.4f64, 1.5, 1.5168, 1.75, 2.0] {
            let expected = ((n - 1.0) / (n + 1.0)).powi(2);
            let got = fresnel_reflectance(1.0, n, 1.0);
            assert!(
                (got - expected).abs() < 1e-12,
                "n = {n}: expected {expected}, got {got}"
            );
        }
        // The textbook figure for crown glass.
        assert!((fresnel_reflectance(1.0, 1.5, 1.0) - 0.04).abs() < 0.0005);
    }

    /// Reflectance rises to 1 at grazing incidence, and is symmetric under swapping
    /// the media at normal incidence.
    #[test]
    fn fresnel_rises_to_grazing_and_is_reciprocal() {
        let mut previous = 0.0;
        for angle_deg in [0.0f64, 30.0, 45.0, 60.0, 75.0, 85.0, 89.9] {
            let r = fresnel_reflectance(1.0, 1.5168, angle_deg.to_radians().cos());
            assert!(
                r >= previous - 1e-12,
                "reflectance should not fall as the angle grows: {angle_deg}° gave {r}"
            );
            previous = r;
        }
        assert!(previous > 0.6, "grazing incidence should be a good mirror");
        assert!(
            (fresnel_reflectance(1.0, 1.5, 1.0) - fresnel_reflectance(1.5, 1.0, 1.0)).abs() < 1e-12
        );
    }

    /// Total internal reflection: past the critical angle nothing gets through, and
    /// the critical angle is where Snell's law runs out — 41.8° for n = 1.5.
    #[test]
    fn total_internal_reflection_starts_at_the_critical_angle() {
        let n = 1.5;
        let critical = critical_angle(n, 1.0).unwrap().to_degrees();
        assert!((critical - 41.81).abs() < 0.01, "critical angle {critical}");
        assert!(fresnel_reflectance(n, 1.0, (critical - 1.0).to_radians().cos()) < 1.0);
        assert_eq!(
            fresnel_reflectance(n, 1.0, (critical + 1.0).to_radians().cos()),
            1.0
        );
        assert_eq!(fresnel_reflectance(n, 1.0, 80f64.to_radians().cos()), 1.0);
        // Going into a denser medium there is no critical angle at all.
        assert_eq!(critical_angle(1.0, 1.5), None);
    }

    /// At Brewster's angle the p reflectance is exactly zero, so the reflection is
    /// purely s-polarised — 56.3° for crown glass in air. Averaging the two hides
    /// this completely, which is why both are exposed.
    #[test]
    fn p_polarisation_vanishes_at_brewsters_angle() {
        let n = 1.5;
        let brewster = brewster_angle(1.0, n);
        assert!((brewster.to_degrees() - 56.31).abs() < 0.01);
        let (rs, rp) = fresnel_split(1.0, n, brewster.cos());
        assert!(rp < 1e-15, "rp should vanish, got {rp}");
        assert!(rs > 0.14, "rs should not, got {rs}");
        // The unpolarised average is half of rs there, and nonzero.
        let unpolarised = fresnel_reflectance(1.0, n, brewster.cos());
        assert!((unpolarised - rs / 2.0).abs() < 1e-15);
        // Away from Brewster, rs always exceeds rp.
        for angle in [10.0f64, 30.0, 70.0, 85.0] {
            let (rs, rp) = fresnel_split(1.0, n, angle.to_radians().cos());
            assert!(rs >= rp - 1e-15, "{angle}°: rs {rs} < rp {rp}");
        }
    }

    /// An AR coating cuts the bare reflection by its residual, and being spectral,
    /// it is worse away from its design band.
    #[test]
    fn an_ar_coating_reduces_the_bare_reflection_spectrally() {
        let bare = SurfaceFinish::Bare;
        let coated = SurfaceFinish::broadband_ar();
        let n = 1.5168;
        let at = |f: &SurfaceFinish, w: f64| f.reflectance_at(1.0, n, 1.0, nm(w));
        assert!((at(&bare, 550.0) - 0.0421).abs() < 0.001);
        // A quarter of a percent in the middle of the band.
        assert!(at(&coated, 550.0) < 0.003, "got {}", at(&coated, 550.0));
        // And measurably worse in the near-IR, as a real V-coat is.
        assert!(at(&coated, 1000.0) > 4.0 * at(&coated, 550.0));
        // A coating never reflects more than the bare surface.
        for w in [350.0, 550.0, 1100.0] {
            assert!(at(&coated, w) <= at(&bare, w) + 1e-12);
        }
    }

    /// The three fractions sum to one at every wavelength, whatever is written down
    /// — that is the invariant the whole tracer relies on.
    #[test]
    fn reflectance_transmittance_and_absorptance_sum_to_one() {
        for optics in [
            SurfaceOptics::black(),
            SurfaceOptics::mirror(0.9),
            SurfaceOptics::beamsplitter(0.5),
            SurfaceOptics::detector(),
            SurfaceOptics::flocked(),
            SurfaceOptics::filter(vec![[500.0, 560.0]], 0.95, 8.0),
            SurfaceOptics::dichroic(vec![[500.0, 560.0]], 0.95, 8.0),
        ] {
            optics
                .validate()
                .expect("built-in optics must conserve energy");
            for w in [400.0, 488.0, 530.0, 610.0, 900.0] {
                let (r, t) = optics.split(nm(w));
                let a = optics.absorptance(nm(w));
                assert!(
                    (r + t + a - 1.0).abs() < 1e-9,
                    "R {r} + T {t} + A {a} at {w} nm"
                );
            }
        }
    }

    /// An over-specified surface is caught, not silently allowed to amplify.
    #[test]
    fn a_surface_that_would_gain_energy_is_rejected() {
        let impossible = SurfaceOptics {
            reflectance: Spectrum::constant(0.7),
            transmittance: Spectrum::constant(0.6),
            diffuse: 0.0,
        };
        let err = impossible.validate().expect_err("1.3 must not be allowed");
        assert_eq!(err.quantity, "energy");
        assert!((err.after - 1.3).abs() < 1e-12, "{err}");
        assert!(err.to_string().contains("created"), "{err}");
        // And if it is traced anyway, `split` renormalises rather than gaining.
        let (r, t) = impossible.split(nm(550.0));
        assert!((r + t - 1.0).abs() < 1e-9);
        assert!(r > t, "the larger share should stay larger");
    }

    /// `diffuse` moves light around without creating or destroying any: it says how
    /// the reflected share leaves, not how much of it there is.
    #[test]
    fn diffuse_scatter_does_not_change_the_energy_budget() {
        let polished = SurfaceOptics {
            reflectance: Spectrum::constant(0.04),
            transmittance: Spectrum::constant(0.0),
            diffuse: 0.0,
        };
        let matte = SurfaceOptics {
            diffuse: 0.9,
            ..polished.clone()
        };
        for w in [400.0, 550.0, 900.0] {
            assert_eq!(polished.split(nm(w)), matte.split(nm(w)));
            assert_eq!(polished.absorptance(nm(w)), matte.absorptance(nm(w)));
        }
        matte
            .validate()
            .expect("a matte surface still conserves energy");
    }

    /// A dichroic passes its band and reflects the rest — the two curves are
    /// complements, which is what makes it a beamsplitter and not a filter.
    #[test]
    fn a_dichroic_reflects_what_it_does_not_pass() {
        let d = SurfaceOptics::dichroic(vec![[495.0, 545.0]], 0.95, 10.0);
        let (r_in, t_in) = d.split(nm(520.0));
        assert!(t_in > 0.9 && r_in < 0.1, "in band: R {r_in}, T {t_in}");
        let (r_out, t_out) = d.split(nm(405.0));
        assert!(
            r_out > 0.9 && t_out < 0.01,
            "out of band: R {r_out}, T {t_out}"
        );
        // And almost nothing is lost either way: a dichroic redirects, it does not
        // absorb.
        assert!(d.absorptance(nm(520.0)) < 0.06);
        assert!(d.absorptance(nm(405.0)) < 0.03);
    }

    /// Sampling the three branches in proportion recovers the energy split, over
    /// enough draws. This is the half of the surface model that did not exist
    /// before: the budget said how much, and now something says where.
    #[test]
    fn sampling_recovers_the_energy_split() {
        let splitter = SurfaceOptics::beamsplitter(0.3);
        let mut rng = Rng::new(20260805);
        let incident = DVec3::new(0.0, 0.0, -1.0);
        let normal = DVec3::Z;
        const N: usize = 200_000;
        let (mut reflected, mut transmitted, mut absorbed) = (0, 0, 0);
        for _ in 0..N {
            match splitter.sample(incident, normal, nm(550.0), &mut rng).0 {
                Scatter::Specular(d) => {
                    // A mirror reflection off +z goes back the way it came.
                    assert!((d - DVec3::Z).length() < 1e-12, "got {d}");
                    reflected += 1;
                }
                Scatter::Diffuse(_) => reflected += 1,
                Scatter::Transmitted(d) => {
                    assert_eq!(d, incident);
                    transmitted += 1;
                }
                Scatter::Absorbed => absorbed += 1,
            }
        }
        let f = |n: usize| n as f64 / N as f64;
        assert!((f(reflected) - 0.3).abs() < 0.005, "R {}", f(reflected));
        assert!((f(transmitted) - 0.7).abs() < 0.005, "T {}", f(transmitted));
        assert_eq!(absorbed, 0, "a lossless splitter absorbs nothing");
    }

    /// A flocked surface scatters nearly all of what little it reflects, and every
    /// scattered ray leaves on the side it arrived from. That is what makes it
    /// glare rather than a ghost.
    #[test]
    fn a_flocked_surface_sprays_its_reflection_into_the_hemisphere() {
        let flock = SurfaceOptics::flocked();
        let mut rng = Rng::new(7);
        let normal = DVec3::Z;
        let incident = DVec3::new(0.4, 0.0, -0.916_515).normalize();
        let mut diffuse = 0;
        let mut specular = 0;
        for _ in 0..100_000 {
            match flock.sample(incident, normal, nm(550.0), &mut rng).0 {
                Scatter::Diffuse(d) => {
                    assert!(d.dot(normal) > -1e-9, "left through the surface: {d}");
                    diffuse += 1;
                }
                Scatter::Specular(_) => specular += 1,
                Scatter::Transmitted(_) => panic!("flock is opaque"),
                Scatter::Absorbed => {}
            }
        }
        assert!(diffuse > 1500, "expected about 2000 diffuse, got {diffuse}");
        assert_eq!(specular, 0, "diffuse 1.0 leaves nothing specular");
    }

    /// Absorbed light is counted, and it is the number a thermal domain receives.
    /// A black surface takes all of it; aluminium takes about 8%.
    #[test]
    fn absorption_is_the_number_a_thermal_domain_would_receive() {
        assert!((SurfaceOptics::black().absorptance(nm(550.0)) - 1.0).abs() < 1e-12);
        let al = SurfaceOptics::aluminium().absorptance(nm(500.0));
        assert!((al - 0.08).abs() < 0.005, "got {al}");
        // Four aluminium folds absorb a third of the beam, which is where the heat
        // in a folded relay comes from.
        let survived = (1.0f64 - al).powi(4);
        assert!(survived < 0.72 && survived > 0.70, "survived {survived}");
    }

    /// Sampling is reproducible, which is the invariant the whole crate rests on —
    /// and it is per-index reproducible, so a parallel trace gets the same answer.
    #[test]
    fn sampling_is_reproducible_per_index() {
        let surface = SurfaceOptics::flocked();
        let draw = |index: u64| {
            let mut rng = Rng::for_index(99, index);
            surface.sample(-DVec3::Z, DVec3::Z, nm(550.0), &mut rng).0
        };
        let forward: Vec<_> = (0..64).map(draw).collect();
        let shuffled: Vec<_> = (0..64).rev().map(draw).collect();
        assert_eq!(
            forward,
            shuffled.into_iter().rev().collect::<Vec<_>>(),
            "the order the rays were traced in must not matter"
        );
    }
}
