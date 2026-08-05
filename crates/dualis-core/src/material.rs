use serde::{Deserialize, Serialize};

use crate::spectrum::Spectrum;

/// An optical material: how it bends light, and how much of it survives the
/// journey through the glass.
///
/// The two are separate properties and both matter. Dispersion decides where the
/// colours land; internal transmittance decides which colours arrive at all. A
/// thick piece of dense flint is visibly yellow, and that is not a coating
/// effect — it is the bulk absorbing blue.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Material {
    /// Fraction surviving 10 mm of glass, per wavelength — the way glass
    /// catalogs quote it. `None` is a perfectly clear material, which no real
    /// glass is but which is often close enough to be the sensible default.
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
    /// Three-term Sellmeier dispersion: n^2 = 1 + sum(B_i * l^2 / (l^2 - C_i)),
    /// wavelength in micrometers. Optional `name` keeps catalog identity.
    Sellmeier {
        #[serde(default)]
        name: String,
        b: [f64; 3],
        c: [f64; 3],
    },
}

impl From<Dispersion> for Material {
    fn from(dispersion: Dispersion) -> Material {
        Material {
            internal_transmittance: None,
            dispersion,
        }
    }
}

impl Material {
    /// A material of fixed index and no absorption.
    pub fn constant(n: f64) -> Material {
        Dispersion::Constant { n }.into()
    }

    /// Refractive index at the given vacuum wavelength in nanometers.
    pub fn index(&self, wavelength_nm: f64) -> f64 {
        match &self.dispersion {
            Dispersion::Constant { n } => *n,
            Dispersion::Sellmeier { b, c, .. } => {
                let l = wavelength_nm / 1000.0; // micrometers
                let l2 = l * l;
                let mut n2 = 1.0;
                for i in 0..3 {
                    n2 += b[i] * l2 / (l2 - c[i]);
                }
                n2.max(1.0).sqrt()
            }
        }
    }

    /// Catalog name, if this came from one.
    pub fn name(&self) -> &str {
        match &self.dispersion {
            Dispersion::Sellmeier { name, .. } => name,
            Dispersion::Constant { .. } => "",
        }
    }

    /// Fraction of light surviving `path_mm` of this glass at this wavelength.
    ///
    /// Beer-Lambert from the catalog figure: transmittance is exponential in
    /// path length, so `t10` over 10 mm becomes `t10^(path/10)`. A clear
    /// material returns 1 and costs nothing.
    pub fn bulk_transmittance(&self, wavelength_nm: f64, path_mm: f64) -> f64 {
        let Some(t10) = &self.internal_transmittance else {
            return 1.0;
        };
        if path_mm <= 0.0 {
            return 1.0; // no glass crossed, nothing absorbed
        }
        let t = t10.at(wavelength_nm).clamp(0.0, 1.0);
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

    #[test]
    fn nbk7_index_at_d_line() {
        let bk7 = Material::from_catalog("N-BK7").unwrap();
        // n_d (587.56 nm) for N-BK7 is 1.5168
        let n = bk7.index(587.56);
        assert!((n - 1.5168).abs() < 1e-3, "n_d = {n}");
    }

    #[test]
    fn dispersion_direction() {
        let bk7 = Material::from_catalog("N-BK7").unwrap();
        // Normal dispersion: blue index > red index.
        assert!(bk7.index(486.13) > bk7.index(656.27));
    }

    #[test]
    fn catalog_key_normalization() {
        assert!(Material::from_catalog("n-bk7").is_some());
        assert!(Material::from_catalog("BK7").is_some());
        assert!(Material::from_catalog("unobtainium").is_none());
    }

    /// Bulk absorption is Beer-Lambert: doubling the path squares the survival
    /// fraction. Checked against the catalog figure it is built from.
    #[test]
    fn bulk_transmittance_follows_beer_lambert() {
        let flint = Material::from_catalog("N-SF11").unwrap();
        let t10 = flint.internal_transmittance.as_ref().unwrap().at(400.0);
        assert!((flint.bulk_transmittance(400.0, 10.0) - t10).abs() < 1e-9);
        // 20 mm is the 10 mm figure squared, and 5 mm its square root.
        assert!((flint.bulk_transmittance(400.0, 20.0) - t10 * t10).abs() < 1e-9);
        assert!((flint.bulk_transmittance(400.0, 5.0) - t10.sqrt()).abs() < 1e-9);
        assert_eq!(flint.bulk_transmittance(400.0, 0.0), 1.0);
        // A material with no curve is perfectly clear, at no cost.
        assert_eq!(Material::constant(1.5).bulk_transmittance(400.0, 50.0), 1.0);
    }

    /// Dense flint really does eat blue, and crown really does not. This is why
    /// a thick flint element looks yellow, and it is a bulk effect, not a
    /// coating one.
    #[test]
    fn dense_flint_absorbs_blue_and_crown_does_not() {
        let flint = Material::from_catalog("N-SF11").unwrap();
        let crown = Material::from_catalog("N-BK7").unwrap();
        let silica = Material::from_catalog("Fused Silica").unwrap();
        // Through 10 mm at 380 nm.
        let (f, c, s) = (
            flint.bulk_transmittance(380.0, 10.0),
            crown.bulk_transmittance(380.0, 10.0),
            silica.bulk_transmittance(380.0, 10.0),
        );
        assert!(
            f < 0.4,
            "N-SF11 should lose most of its 380 nm light, got {f}"
        );
        assert!(c > 0.95, "N-BK7 should pass 380 nm, got {c}");
        assert!(s > 0.99, "fused silica should be clear at 380 nm, got {s}");
        // And all of them are effectively clear in the green.
        for m in [&flint, &crown, &silica] {
            assert!(m.bulk_transmittance(550.0, 10.0) > 0.99);
        }
    }
}
