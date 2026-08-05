//! Refractive index against wavelength, and how much light survives the glass.
//!
//! The two are separate properties and both matter. Dispersion decides where the
//! colours land; internal transmittance decides which colours arrive at all. A
//! thick piece of dense flint is visibly yellow, and that is not a coating effect —
//! it is the bulk absorbing blue.
//!
//! This is the *optical* half of what a material is. The thermal, mechanical and
//! acoustic halves are [`Substance`](dualis_core::Substance) in the kernel, which
//! must not know that optics exists; a simulation that heats a lens pairs the two.

use dualis_units::Length;
use serde::{Deserialize, Serialize};

use crate::spectrum::Spectrum;

/// An optical material: how it bends light, and how much of it survives the
/// journey through the glass.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Material {
    /// Fraction surviving 10 mm of glass, per wavelength — the way glass catalogs
    /// quote it. `None` is a perfectly clear material, which no real glass is but
    /// which is often close enough to be the sensible default.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub internal_transmittance: Option<Spectrum>,
    #[serde(flatten)]
    pub dispersion: Dispersion,
}

/// How a material's refractive index varies with wavelength.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Dispersion {
    /// Wavelength-independent index (ideal glass).
    Constant { n: f64 },
    /// Three-term Sellmeier dispersion: n² = 1 + Σ Bᵢ λ² / (λ² - Cᵢ), wavelength in
    /// micrometres. Optional `name` keeps catalog identity.
    Sellmeier {
        #[serde(default)]
        name: String,
        b: [f64; 3],
        c: [f64; 3],
    },
    /// Cauchy dispersion: n = A + B/λ² + C/λ⁴, λ in micrometres.
    ///
    /// Less accurate than Sellmeier and valid only across the visible, but it is
    /// the form a great deal of published data for plastics, liquids and coatings
    /// comes in — and a two-term fit is often all that was ever measured.
    Cauchy { a: f64, b: f64, c: f64 },
}

impl From<Dispersion> for Material {
    fn from(dispersion: Dispersion) -> Material {
        Material {
            internal_transmittance: None,
            dispersion,
        }
    }
}

/// The Fraunhofer lines an optical prescription is written against: the helium d
/// line and the hydrogen F and C lines. A glass is specified at these three
/// wavelengths and nowhere else.
pub const D_LINE: Length = Length::from_si(587.5618e-9);
pub const F_LINE: Length = Length::from_si(486.1327e-9);
pub const C_LINE: Length = Length::from_si(656.2725e-9);

impl Material {
    /// A material of fixed index and no absorption.
    pub fn constant(n: f64) -> Material {
        Dispersion::Constant { n }.into()
    }

    /// Vacuum, and close enough to air for anything but a metrology simulation.
    pub fn air() -> Material {
        Material::constant(1.0)
    }

    /// Water at 20 °C, by a Cauchy fit — the immersion medium for a dipping
    /// objective, and the thing a sample sits in.
    pub fn water() -> Material {
        Material {
            internal_transmittance: Some(Spectrum::curve(vec![
                (350.0, 0.9999),
                (700.0, 0.9995),
                (900.0, 0.993),
                (1100.0, 0.94),
            ])),
            dispersion: Dispersion::Cauchy {
                a: 1.324,
                b: 0.00306,
                c: 0.0,
            },
        }
    }

    /// Microscope immersion oil, matched to N-BK7 at the d line — which is the
    /// whole point of it: an index step of nothing means no reflection and no
    /// refraction at the coverslip, so the objective can collect past NA 1.
    pub fn immersion_oil() -> Material {
        Material {
            internal_transmittance: Some(Spectrum::curve(vec![
                (400.0, 0.96),
                (450.0, 0.99),
                (550.0, 0.998),
                (900.0, 0.995),
            ])),
            dispersion: Dispersion::Cauchy {
                a: 1.5064,
                b: 0.00381,
                c: 0.0,
            },
        }
    }

    /// Refractive index at a vacuum wavelength.
    pub fn index(&self, wavelength: Length) -> f64 {
        self.index_nm(wavelength.in_nm())
    }

    /// Refractive index at a wavelength in nanometres — what the sampling loops in
    /// this crate call.
    pub fn index_nm(&self, wavelength_nm: f64) -> f64 {
        match &self.dispersion {
            Dispersion::Constant { n } => *n,
            Dispersion::Sellmeier { b, c, .. } => {
                let l = wavelength_nm / 1000.0; // micrometres
                let l2 = l * l;
                let mut n2 = 1.0;
                for i in 0..3 {
                    n2 += b[i] * l2 / (l2 - c[i]);
                }
                n2.max(1.0).sqrt()
            }
            Dispersion::Cauchy { a, b, c } => {
                let l2 = (wavelength_nm / 1000.0).powi(2);
                (a + b / l2 + c / (l2 * l2)).max(1.0)
            }
        }
    }

    /// Index at the d line, which is what "n = 1.5168" means when a catalogue says
    /// it without qualification.
    pub fn n_d(&self) -> f64 {
        self.index(D_LINE)
    }

    /// Abbe number, `v_d = (n_d - 1) / (n_F - n_C)` — how much the glass spreads
    /// colour relative to how much it bends it.
    ///
    /// This is the second number every glass is known by, and the one that decides
    /// what it is *for*: a crown near 64 and a flint near 25 cancel each other's
    /// chromatic aberration in an achromatic doublet, which no single glass can do.
    /// A dispersionless material returns infinity, honestly — it spreads nothing.
    pub fn abbe(&self) -> f64 {
        let spread = self.index(F_LINE) - self.index(C_LINE);
        if spread.abs() < 1e-12 {
            return f64::INFINITY;
        }
        (self.n_d() - 1.0) / spread
    }

    /// dn/dλ, per metre of wavelength, by a central difference over 1 nm.
    ///
    /// Negative for normal dispersion, since blue light is bent more than red. This
    /// is what a chromatic focal shift is computed from: a lens's power follows
    /// `n - 1`, so a spread in `n` is a spread in focus.
    pub fn dn_dlambda(&self, wavelength: Length) -> f64 {
        let w = wavelength.in_nm();
        let h = 1.0;
        // Per nm from the difference, then per metre for the return.
        (self.index_nm(w + h) - self.index_nm(w - h)) / (2.0 * h) * 1e9
    }

    /// Catalog name, if this came from one.
    pub fn name(&self) -> &str {
        match &self.dispersion {
            Dispersion::Sellmeier { name, .. } => name,
            _ => "",
        }
    }

    /// Fraction of light surviving `path` of this glass at this wavelength.
    ///
    /// Beer-Lambert from the catalog figure: transmittance is exponential in path
    /// length, so `t10` over 10 mm becomes `t10^(path/10mm)`. A clear material
    /// returns 1 and costs nothing.
    pub fn bulk_transmittance(&self, wavelength: Length, path: Length) -> f64 {
        let Some(t10) = &self.internal_transmittance else {
            return 1.0;
        };
        let path_mm = path.in_mm();
        if path_mm <= 0.0 {
            return 1.0; // no glass crossed, nothing absorbed
        }
        let t = t10.at(wavelength).clamp(0.0, 1.0);
        if t >= 1.0 {
            return 1.0;
        }
        if t <= 0.0 {
            return 0.0;
        }
        t.powf(path_mm / 10.0)
    }

    /// Look up a glass by catalog name (case-insensitive, '-' optional).
    pub fn from_catalog(name: &str) -> Option<Material> {
        let key: String = name
            .to_ascii_uppercase()
            .chars()
            .filter(|ch| ch.is_ascii_alphanumeric())
            .collect();
        let (b, c) = match key.as_str() {
            "NBK7" | "BK7" => (
                [1.039_612_12, 0.231_792_344, 1.010_469_45],
                [0.006_000_698_67, 0.020_017_914_4, 103.560_653],
            ),
            "NSF11" | "SF11" => (
                [1.737_596_95, 0.313_747_346, 1.898_781_01],
                [0.013_188_707, 0.062_306_814_2, 155.236_29],
            ),
            "F2" => (
                [1.345_333_59, 0.209_073_176, 0.937_357_162],
                [0.009_977_438_71, 0.047_045_076_7, 111.886_764],
            ),
            "NSF5" | "SF5" => (
                [1.521_463_48, 0.198_570_18, 1.659_180_6],
                [0.011_931_495_4, 0.054_690_775_3, 112.041_888],
            ),
            "NBAF10" | "BAF10" => (
                [1.585_149_5, 0.143_559_385, 1.085_212_69],
                [0.009_265_681_25, 0.042_412_386_2, 105.613_573],
            ),
            "FUSEDSILICA" | "SILICA" | "SIO2" | "FSILICA" => (
                [0.696_166_3, 0.407_942_6, 0.897_479_4],
                [0.004_679_148_26, 0.013_512_063_1, 97.934_002_5],
            ),
            _ => return None,
        };
        Some(Material {
            internal_transmittance: Some(internal_transmittance(&key)),
            dispersion: Dispersion::Sellmeier {
                name: name.to_string(),
                b,
                c,
            },
        })
    }

    /// Names of all built-in catalog glasses.
    pub fn catalog_names() -> &'static [&'static str] {
        &["N-BK7", "N-SF11", "F2", "N-SF5", "N-BAF10", "Fused Silica"]
    }
}

/// Internal transmittance per 10 mm, from the shape catalogs publish.
///
/// The pattern is always the same: near-perfect through the visible, falling off
/// into the near-UV. Where it falls off is what separates the glasses — fused
/// silica is still clear at 300 nm, while a dense flint like N-SF11 has lost most
/// of its 380 nm light and looks faintly yellow in thick pieces.
fn internal_transmittance(key: &str) -> Spectrum {
    let samples = match key {
        "FUSEDSILICA" | "SILICA" | "SIO2" | "FSILICA" => vec![
            (250.0, 0.99),
            (300.0, 0.999),
            (400.0, 0.9995),
            (700.0, 0.9995),
            (1100.0, 0.999),
        ],
        "NSF11" | "SF11" => vec![
            (350.0, 0.02),
            (380.0, 0.28),
            (400.0, 0.75),
            (420.0, 0.93),
            (450.0, 0.985),
            (500.0, 0.995),
            (700.0, 0.998),
            (1100.0, 0.997),
        ],
        "NSF5" | "SF5" => vec![
            (350.0, 0.10),
            (380.0, 0.62),
            (400.0, 0.88),
            (420.0, 0.965),
            (450.0, 0.99),
            (500.0, 0.996),
            (700.0, 0.998),
            (1100.0, 0.997),
        ],
        "F2" => vec![
            (350.0, 0.08),
            (380.0, 0.55),
            (400.0, 0.85),
            (420.0, 0.96),
            (450.0, 0.99),
            (500.0, 0.996),
            (700.0, 0.998),
            (1100.0, 0.997),
        ],
        "NBAF10" | "BAF10" => vec![
            (350.0, 0.55),
            (380.0, 0.90),
            (400.0, 0.965),
            (450.0, 0.993),
            (500.0, 0.997),
            (700.0, 0.998),
            (1100.0, 0.997),
        ],
        // N-BK7 and anything else crown-like: clear from the near-UV up.
        _ => vec![
            (300.0, 0.30),
            (334.0, 0.80),
            (365.0, 0.96),
            (400.0, 0.993),
            (450.0, 0.998),
            (700.0, 0.999),
            (1100.0, 0.998),
        ],
    };
    Spectrum::curve(samples)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nm(v: f64) -> Length {
        Length::nm(v)
    }

    #[test]
    fn nbk7_index_at_d_line() {
        let bk7 = Material::from_catalog("N-BK7").unwrap();
        // n_d (587.56 nm) for N-BK7 is 1.5168
        let n = bk7.n_d();
        assert!((n - 1.5168).abs() < 1e-3, "n_d = {n}");
        assert_eq!(n, bk7.index(D_LINE));
    }

    #[test]
    fn dispersion_direction() {
        let bk7 = Material::from_catalog("N-BK7").unwrap();
        // Normal dispersion: blue index > red index.
        assert!(bk7.index(F_LINE) > bk7.index(C_LINE));
        // Which means dn/dlambda is negative.
        assert!(bk7.dn_dlambda(nm(550.0)) < 0.0);
    }

    #[test]
    fn catalog_key_normalization() {
        assert!(Material::from_catalog("n-bk7").is_some());
        assert!(Material::from_catalog("BK7").is_some());
        assert!(Material::from_catalog("unobtainium").is_none());
    }

    /// The two numbers a glass is known by, against the catalogue: N-BK7 is
    /// 1.5168/64.2 and N-SF11 is 1.7847/25.7. A crown and a flint that far apart in
    /// Abbe number are what make an achromatic doublet possible.
    #[test]
    fn abbe_numbers_match_the_catalogue() {
        let crown = Material::from_catalog("N-BK7").unwrap();
        let flint = Material::from_catalog("N-SF11").unwrap();
        assert!(
            (crown.abbe() - 64.17).abs() < 0.3,
            "crown v_d {}",
            crown.abbe()
        );
        assert!(
            (flint.abbe() - 25.68).abs() < 0.3,
            "flint v_d {}",
            flint.abbe()
        );
        // The flint bends more and spreads far more, which is the trade.
        assert!(flint.n_d() > crown.n_d());
        assert!(flint.abbe() < crown.abbe() / 2.0);
        // A dispersionless material spreads nothing, and says so.
        assert_eq!(Material::constant(1.5).abbe(), f64::INFINITY);
        assert_eq!(Material::air().n_d(), 1.0);
    }

    /// Immersion oil exists to match the glass: at the d line the index step is
    /// small enough that the coverslip interface almost vanishes, which is what
    /// lets an objective work past NA 1.
    #[test]
    fn immersion_oil_matches_the_crown_it_is_used_with() {
        let oil = Material::immersion_oil();
        let glass = Material::from_catalog("N-BK7").unwrap();
        assert!(
            (oil.n_d() - glass.n_d()).abs() < 0.005,
            "oil {} vs glass {}",
            oil.n_d(),
            glass.n_d()
        );
        // Water does not match, which is why a water objective is a different
        // design and not just the same lens with a different drop on it.
        assert!((Material::water().n_d() - glass.n_d()).abs() > 0.15);
        assert!((Material::water().n_d() - 1.333).abs() < 0.003);
    }

    /// Bulk absorption is Beer-Lambert: doubling the path squares the survival
    /// fraction. Checked against the catalog figure it is built from.
    #[test]
    fn bulk_transmittance_follows_beer_lambert() {
        let flint = Material::from_catalog("N-SF11").unwrap();
        let t10 = flint.internal_transmittance.as_ref().unwrap().at(nm(400.0));
        let at = |mm: f64| flint.bulk_transmittance(nm(400.0), Length::mm(mm));
        assert!((at(10.0) - t10).abs() < 1e-9);
        // 20 mm is the 10 mm figure squared, and 5 mm its square root.
        assert!((at(20.0) - t10 * t10).abs() < 1e-9);
        assert!((at(5.0) - t10.sqrt()).abs() < 1e-9);
        assert_eq!(at(0.0), 1.0);
        // A material with no curve is perfectly clear, at no cost.
        assert_eq!(
            Material::constant(1.5).bulk_transmittance(nm(400.0), Length::mm(50.0)),
            1.0
        );
    }

    /// Dense flint really does eat blue, and crown really does not. This is why a
    /// thick flint element looks yellow, and it is a bulk effect, not a coating one.
    #[test]
    fn dense_flint_absorbs_blue_and_crown_does_not() {
        let flint = Material::from_catalog("N-SF11").unwrap();
        let crown = Material::from_catalog("N-BK7").unwrap();
        let silica = Material::from_catalog("Fused Silica").unwrap();
        let through = |m: &Material, w: f64| m.bulk_transmittance(nm(w), Length::mm(10.0));
        let (f, c, s) = (
            through(&flint, 380.0),
            through(&crown, 380.0),
            through(&silica, 380.0),
        );
        assert!(
            f < 0.4,
            "N-SF11 should lose most of its 380 nm light, got {f}"
        );
        assert!(c > 0.95, "N-BK7 should pass 380 nm, got {c}");
        assert!(s > 0.99, "fused silica should be clear at 380 nm, got {s}");
        // And all of them are effectively clear in the green.
        for m in [&flint, &crown, &silica] {
            assert!(through(m, 550.0) > 0.99);
        }
    }

    /// A Cauchy fit reproduces the shape a Sellmeier one does over the visible,
    /// which is the only claim it makes.
    #[test]
    fn cauchy_dispersion_behaves_like_dispersion() {
        let oil = Material::immersion_oil();
        assert!(oil.index(F_LINE) > oil.index(C_LINE), "normal dispersion");
        assert!(oil.dn_dlambda(nm(550.0)) < 0.0);
        assert!(oil.abbe() > 20.0 && oil.abbe() < 80.0, "v_d {}", oil.abbe());
        // Index never drops below vacuum, however the fit is extrapolated.
        assert!(oil.index(nm(50_000.0)) >= 1.0);
    }
}
