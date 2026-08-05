//! What a surface does to light: reflect it, transmit it, or absorb it.
//!
//! Every surface splits incident power three ways, and the three add to one.
//! That is not a modelling choice, it is conservation of energy, so this module
//! stores only two of them — [`SurfaceOptics::reflectance`] and
//! [`SurfaceOptics::transmittance`] — and computes absorptance as the remainder.
//! There is no way to write down a surface that gains energy.
//!
//! Glass is the exception, and an important one: its reflectance is *not* a free
//! parameter. It follows from the refractive indices either side and the angle
//! of incidence, by Fresnel's equations. A lens surface therefore carries a
//! [`SurfaceFinish`] — bare, or coated — rather than a reflectance curve, and
//! the number comes out of the physics.

use serde::{Deserialize, Serialize};

use crate::spectrum::{Spectrum, VISIBLE_RANGE_NM};

/// Spectral reflectance and transmittance of a surface. Absorptance is whatever
/// is left, so the three always sum to one.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SurfaceOptics {
    /// Fraction sent back, per wavelength. `diffuse` decides how much of it
    /// leaves specularly and how much scatters.
    pub reflectance: Spectrum,
    /// Fraction let through, per wavelength.
    pub transmittance: Spectrum,
    /// How much of the reflected light leaves diffusely rather than as a mirror
    /// image, 0..1.
    ///
    /// This does not change the energy budget — the light was already counted in
    /// `reflectance` — only where it goes. And where it goes matters enormously:
    /// a polished surface sends its reflection somewhere definite, where you can
    /// design a baffle to catch it, while a matte one sprays it over the whole
    /// hemisphere. That spray is **veiling glare**, the haze that lifts the black
    /// level of an otherwise good image, and it is why the inside of a lens
    /// barrel is flocked rather than merely painted.
    #[serde(default)]
    pub diffuse: f64,
}

impl SurfaceOptics {
    /// Absorptance at a wavelength: everything neither reflected nor
    /// transmitted. Clamped at zero, so an over-specified surface loses the
    /// excess rather than manufacturing light.
    pub fn absorptance(&self, wavelength_nm: f64) -> f64 {
        (1.0 - self.reflectance.at(wavelength_nm) - self.transmittance.at(wavelength_nm)).max(0.0)
    }

    /// Reflectance and transmittance at a wavelength, rescaled if they sum past
    /// one so that no surface ever returns more light than reached it.
    pub fn split(&self, wavelength_nm: f64) -> (f64, f64) {
        let r = self.reflectance.at(wavelength_nm).clamp(0.0, 1.0);
        let t = self.transmittance.at(wavelength_nm).clamp(0.0, 1.0);
        let total = r + t;
        if total > 1.0 {
            (r / total, t / total)
        } else {
            (r, t)
        }
    }

    /// Whether reflectance and transmittance stay within the energy budget over
    /// the whole working range. `Err` names the wavelength that broke it.
    pub fn validate(&self) -> Result<(), String> {
        let (lo, hi) = VISIBLE_RANGE_NM;
        const STEPS: usize = 150;
        for i in 0..=STEPS {
            let w = lo + (hi - lo) * i as f64 / STEPS as f64;
            let (r, t) = (self.reflectance.at(w), self.transmittance.at(w));
            if r < 0.0 || t < 0.0 {
                return Err(format!(
                    "negative reflectance or transmittance at {w:.0} nm (R {r:.3}, T {t:.3})"
                ));
            }
            if r + t > 1.0 + 1e-9 {
                return Err(format!(
                    "reflectance plus transmittance is {:.3} at {w:.0} nm, which would \
                     create light out of nothing",
                    r + t
                ));
            }
        }
        Ok(())
    }

    /// Perfectly opaque and black: absorbs everything. An aperture blade, a
    /// bench, the inside of a tube.
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

    /// A dichroic beamsplitter: passes its bands and *reflects* what it blocks,
    /// so the blocked light leaves down another arm instead of being lost. This
    /// is what lets one objective carry excitation down and fluorescence back.
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

    /// A detector face: absorbs nearly everything, and reflects the few percent
    /// that a real sensor's cover glass and silicon send back — which is where
    /// sensor ghosting comes from.
    pub fn detector() -> SurfaceOptics {
        SurfaceOptics {
            reflectance: Spectrum::constant(0.03),
            transmittance: Spectrum::constant(0.0),
            diffuse: 0.0,
        }
    }
}

/// How a glass surface is finished, which is the *only* freedom there is: the
/// bare reflectance itself follows from the refractive indices and the angle.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum SurfaceFinish {
    /// Bare polished glass. About 4% per surface at normal incidence for
    /// n = 1.5, rising steeply past 60°, and the reason an uncoated multi-element
    /// lens is dim and full of ghosts.
    Bare,
    /// An anti-reflection coating. `residual` is the fraction of the bare Fresnel
    /// reflectance left, per wavelength: a good broadband V-coat leaves about
    /// 0.05 of it (4% becomes 0.2%) in its design band and degrades outside it,
    /// which is why coatings are specified spectrally.
    Coated { residual: Spectrum },
    /// A mirrored surface: reflect rather than refract, with this reflectance.
    Mirrored { reflectance: Spectrum },
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
    pub fn reflectance_at(&self, n1: f64, n2: f64, cos_incident: f64, wavelength_nm: f64) -> f64 {
        match self {
            SurfaceFinish::Bare => fresnel_reflectance(n1, n2, cos_incident),
            SurfaceFinish::Coated { residual } => {
                fresnel_reflectance(n1, n2, cos_incident) * residual.at(wavelength_nm).max(0.0)
            }
            SurfaceFinish::Mirrored { reflectance } => {
                reflectance.at(wavelength_nm).clamp(0.0, 1.0)
            }
        }
    }
}

impl Default for SurfaceFinish {
    /// Bare, because that is what a piece of glass is until someone coats it.
    /// A setup that says nothing about coatings gets the honest 4% per surface.
    fn default() -> Self {
        SurfaceFinish::Bare
    }
}

/// Unpolarised Fresnel reflectance at an interface, averaging the s and p
/// polarisations.
///
/// `cos_incident` is |cos θ| in the incident medium. Past the critical angle
/// this returns 1: total internal reflection, which is Fresnel's equations
/// telling you there is nowhere for the transmitted ray to go.
pub fn fresnel_reflectance(n1: f64, n2: f64, cos_incident: f64) -> f64 {
    let cos_i = cos_incident.clamp(0.0, 1.0);
    let sin_i2 = 1.0 - cos_i * cos_i;
    let ratio = n1 / n2;
    let sin_t2 = ratio * ratio * sin_i2;
    if sin_t2 >= 1.0 {
        return 1.0; // total internal reflection
    }
    let cos_t = (1.0 - sin_t2).sqrt();
    let rs = ((n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)).powi(2);
    let rp = ((n1 * cos_t - n2 * cos_i) / (n1 * cos_t + n2 * cos_i)).powi(2);
    ((rs + rp) / 2.0).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Fresnel at normal incidence reduces to ((n2-n1)/(n2+n1))^2 — 4.0% for
    /// n = 1.5 in air, the number every optics text quotes.
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

    /// Reflectance rises to 1 at grazing incidence, and is symmetric under
    /// swapping the media at normal incidence.
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

    /// Total internal reflection: past the critical angle nothing gets through,
    /// and the critical angle is where Snell's law runs out — 41.8° for n = 1.5.
    #[test]
    fn total_internal_reflection_starts_at_the_critical_angle() {
        let n = 1.5;
        let critical = (1.0f64 / n).asin().to_degrees();
        assert!((critical - 41.81).abs() < 0.01, "critical angle {critical}");
        assert!(fresnel_reflectance(n, 1.0, (critical - 1.0).to_radians().cos()) < 1.0);
        assert_eq!(
            fresnel_reflectance(n, 1.0, (critical + 1.0).to_radians().cos()),
            1.0
        );
        assert_eq!(fresnel_reflectance(n, 1.0, 80f64.to_radians().cos()), 1.0);
    }

    /// An AR coating cuts the bare reflection by its residual, and being
    /// spectral, it is worse away from its design band.
    #[test]
    fn an_ar_coating_reduces_the_bare_reflection_spectrally() {
        let bare = SurfaceFinish::Bare;
        let coated = SurfaceFinish::broadband_ar();
        let n = 1.5168;
        let at = |f: &SurfaceFinish, w: f64| f.reflectance_at(1.0, n, 1.0, w);
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

    /// The three fractions sum to one at every wavelength, whatever is written
    /// down — that is the invariant the whole tracer relies on.
    #[test]
    fn reflectance_transmittance_and_absorptance_sum_to_one() {
        for optics in [
            SurfaceOptics::black(),
            SurfaceOptics::mirror(0.9),
            SurfaceOptics::beamsplitter(0.5),
            SurfaceOptics::detector(),
            SurfaceOptics::filter(vec![[500.0, 560.0]], 0.95, 8.0),
            SurfaceOptics::dichroic(vec![[500.0, 560.0]], 0.95, 8.0),
        ] {
            optics
                .validate()
                .expect("built-in optics must conserve energy");
            for w in [400.0, 488.0, 530.0, 610.0, 900.0] {
                let (r, t) = optics.split(w);
                let a = optics.absorptance(w);
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
        assert!(err.contains("out of nothing"), "{err}");
        // And if it is traced anyway, `split` renormalises rather than gaining.
        let (r, t) = impossible.split(550.0);
        assert!((r + t - 1.0).abs() < 1e-9);
        assert!(r > t, "the larger share should stay larger");
    }

    /// `diffuse` moves light around without creating or destroying any: it says
    /// how the reflected share leaves, not how much of it there is.
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
            assert_eq!(polished.split(w), matte.split(w));
            assert_eq!(polished.absorptance(w), matte.absorptance(w));
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
        let (r_in, t_in) = d.split(520.0);
        assert!(t_in > 0.9 && r_in < 0.1, "in band: R {r_in}, T {t_in}");
        let (r_out, t_out) = d.split(405.0);
        assert!(
            r_out > 0.9 && t_out < 0.01,
            "out of band: R {r_out}, T {t_out}"
        );
        // And almost nothing is lost either way: a dichroic redirects, it does
        // not absorb.
        assert!(d.absorptance(520.0) < 0.06);
        assert!(d.absorptance(405.0) < 0.03);
    }
}
