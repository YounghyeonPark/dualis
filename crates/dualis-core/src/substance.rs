//! What a piece of matter is, across every domain that cares.
//!
//! N-BK7 is not only a refractive index. It is 2.51 g/cm³, it conducts 1.11 W/m·K,
//! it holds 858 J/kg·K, it expands 7.1 ppm per kelvin and it fails at about
//! 60 MPa. A thermal solver needs the middle three, a mechanical one the last two,
//! and an optical one the index — but it is *one piece of glass*, and if each
//! domain carries its own idea of what it is made of then nothing can be coupled,
//! because there is no single object for a coupling to be about.
//!
//! So the properties live together and each domain reads the part it needs.
//! Everything is [`Option`]: a thermal-only simulation is not made to invent a
//! Young's modulus, and a property that is absent says so rather than defaulting
//! to a plausible lie.
//!
//! # Optics is deliberately missing
//!
//! There is no optical field here. This crate is the kernel and must not know that
//! optics exists — refractive index is `dualis-optics`'s
//! `Material`, and a consumer that needs both pairs them. Putting it here would
//! make the kernel depend on a domain, which is the one structural rule the split
//! exists to enforce.

use dualis_units::{
    Density, Diffusivity, HeatCapacity, Length, Mass, Pressure, SpecificHeat, Temperature,
    ThermalConductivity, ThermalExpansion, Velocity, Volume,
};
use serde::{Deserialize, Serialize};

/// A material, as much of it as is known.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Substance {
    /// What it is called. Free text: a catalogue designation, a common name, whatever the
    /// caller will recognise in a violation message.
    pub name: String,
    /// kg·m⁻³. The one property everything has, which is why it is not optional.
    pub density: Density,
    /// Conductivity, specific heat, emissivity and a service limit, if they are known.
    ///
    /// Optional because a substance is often only known as far as it needed to be. A domain
    /// asking for what is not here gets `None` rather than a plausible default, which is the
    /// difference between "unknown" and "zero".
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub thermal: Option<ThermalProps>,
    /// Stiffness, restitution and friction, if they are known.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mechanical: Option<MechanicalProps>,
    /// Sound speed and absorption, if they are known.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub acoustic: Option<AcousticProps>,
}

/// What it does with heat.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ThermalProps {
    /// Fourier's `k`: how fast heat moves through it.
    pub conductivity: ThermalConductivity,
    /// `c_p`: how much heat it takes to warm it.
    pub specific_heat: SpecificHeat,
    /// Linear expansion per kelvin — the property that turns absorbed light into
    /// a focus shift.
    pub expansion: ThermalExpansion,
    /// Emissivity, 0..1, for radiative exchange. 1 is a blackbody; polished metal
    /// is near 0.05, which is why a shiny shield works.
    pub emissivity: f64,
}

/// What it does under load.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct MechanicalProps {
    /// Young's modulus.
    pub youngs_modulus: Pressure,
    /// Poisson's ratio, dimensionless — how much it bulges sideways when squeezed.
    pub poisson_ratio: f64,
    /// Where it stops coming back. For a brittle material this is the fracture
    /// stress, and there is no plastic region before it.
    pub yield_strength: Pressure,
}

/// What it does with sound.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct AcousticProps {
    /// Longitudinal speed of sound.
    pub sound_speed: Velocity,
}

impl Substance {
    /// Thermal diffusivity, `α = k / (ρ c_p)` — the m²/s that decides how fast a
    /// temperature *front* moves, as opposed to how much heat flows.
    ///
    /// This is the number that sets an explicit heat solver's stability limit
    /// (`dt < dx²/2α`), so a thermal domain asks for it before it can say how big
    /// a step it can take.
    pub fn diffusivity(&self) -> Option<Diffusivity> {
        let t = self.thermal?;
        Some(Diffusivity::from_si(
            t.conductivity.to_si() / (self.density.to_si() * t.specific_heat.to_si()),
        ))
    }

    /// Heat capacity of a given volume of this substance.
    pub fn heat_capacity(&self, volume: Volume) -> Option<HeatCapacity> {
        let t = self.thermal?;
        Some(self.mass_of(volume) * t.specific_heat)
    }

    /// Mass of a given volume of it.
    pub fn mass_of(&self, volume: Volume) -> Mass {
        self.density * volume
    }

    /// How much a length of this substance grows for a temperature rise. Linear,
    /// which is a good approximation for the tens of kelvin an instrument sees and
    /// a poor one for hundreds.
    pub fn expansion_of(&self, length: Length, rise: Temperature) -> Option<Length> {
        let t = self.thermal?;
        Some(Length::from_si(
            length.to_si() * t.expansion.to_si() * rise.to_si(),
        ))
    }

    /// Stress produced by preventing that expansion — the reason a lens bonded
    /// rigidly into a metal mount cracks when it warms up.
    ///
    /// `σ = E α ΔT`, independent of size, which is why scaling the part down does
    /// not help.
    pub fn constrained_stress(&self, rise: Temperature) -> Option<Pressure> {
        let t = self.thermal?;
        let m = self.mechanical?;
        Some(Pressure::from_si(
            m.youngs_modulus.to_si() * t.expansion.to_si() * rise.to_si(),
        ))
    }

    /// Whether that stress would break it.
    pub fn survives(&self, rise: Temperature) -> Option<bool> {
        let stress = self.constrained_stress(rise)?;
        let limit = self.mechanical?.yield_strength;
        Some(stress < limit)
    }

    /// N-BK7, the borosilicate crown that most of an optical bench is made of.
    pub fn borosilicate_crown() -> Substance {
        Substance {
            name: "N-BK7".to_string(),
            density: Density::g_per_cm3(2.51),
            thermal: Some(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(1.114),
                specific_heat: SpecificHeat::j_per_kg_k(858.0),
                expansion: ThermalExpansion::ppm_per_k(7.1),
                emissivity: 0.90,
            }),
            mechanical: Some(MechanicalProps {
                youngs_modulus: Pressure::from_si(82.0e9),
                poisson_ratio: 0.206,
                // Brittle: this is a fracture stress, not a yield point.
                yield_strength: Pressure::from_si(60.0e6),
            }),
            acoustic: Some(AcousticProps {
                sound_speed: Velocity::m_per_s(5_680.0),
            }),
        }
    }

    /// 6061 aluminium: what the mount holding the glass is made of, and the
    /// reason a mount-and-lens pair moves when the room does — its expansion is
    /// three times the glass's.
    pub fn aluminium_6061() -> Substance {
        Substance {
            name: "Al 6061".to_string(),
            density: Density::g_per_cm3(2.70),
            thermal: Some(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(167.0),
                specific_heat: SpecificHeat::j_per_kg_k(896.0),
                expansion: ThermalExpansion::ppm_per_k(23.6),
                emissivity: 0.09,
            }),
            mechanical: Some(MechanicalProps {
                youngs_modulus: Pressure::from_si(68.9e9),
                poisson_ratio: 0.33,
                yield_strength: Pressure::from_si(276.0e6),
            }),
            acoustic: Some(AcousticProps {
                sound_speed: Velocity::m_per_s(6_320.0),
            }),
        }
    }

    /// Water at 20 °C.
    pub fn water() -> Substance {
        Substance {
            name: "water".to_string(),
            density: Density::g_per_cm3(0.998),
            thermal: Some(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(0.598),
                specific_heat: SpecificHeat::j_per_kg_k(4_182.0),
                expansion: ThermalExpansion::ppm_per_k(69.0),
                emissivity: 0.96,
            }),
            mechanical: None,
            acoustic: Some(AcousticProps {
                sound_speed: Velocity::m_per_s(1_482.0),
            }),
        }
    }

    /// A substance with nothing known but how heavy it is.
    pub fn bulk(name: &str, density: Density) -> Substance {
        Substance {
            name: name.to_string(),
            density,
            thermal: None,
            mechanical: None,
            acoustic: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dualis_units::Length;

    /// Diffusivity against the published figure: N-BK7 is about 5.2e-7 m²/s, and
    /// aluminium is 130 times faster, which is why one of them is a heat spreader
    /// and the other is not.
    #[test]
    fn diffusivity_matches_the_published_figures() {
        let glass = Substance::borosilicate_crown().diffusivity().unwrap();
        let metal = Substance::aluminium_6061().diffusivity().unwrap();
        assert!(
            (glass.to_si() - 5.17e-7).abs() < 1e-8,
            "N-BK7 diffusivity {glass:?}"
        );
        assert!(
            (metal.to_si() - 6.9e-5).abs() < 1e-6,
            "aluminium diffusivity {metal:?}"
        );
        assert!(metal.to_si() / glass.to_si() > 100.0);
    }

    /// The number an explicit heat solver needs: over a 1 mm cell, N-BK7 is stable
    /// to about a second and aluminium to about 7 ms. That two-orders-of-magnitude
    /// gap between two parts of the same instrument is exactly why
    /// `Schedule::Multirate` exists.
    #[test]
    fn stability_limits_differ_by_two_orders_of_magnitude() {
        let cell = Length::mm(1.0);
        let limit = |s: &Substance| {
            let a = s.diffusivity().unwrap().to_si();
            cell.to_si() * cell.to_si() / (2.0 * a)
        };
        let glass = limit(&Substance::borosilicate_crown());
        let metal = limit(&Substance::aluminium_6061());
        assert!((glass - 0.97).abs() < 0.1, "glass limit {glass} s");
        assert!((metal - 0.0072).abs() < 0.001, "metal limit {metal} s");
        assert!(glass / metal > 100.0);
    }

    /// Heat capacity of a real piece of glass: a 25 mm disc 5 mm thick is 6.2 g
    /// and holds 5.3 J per kelvin.
    #[test]
    fn a_lens_sized_piece_holds_a_few_joules_per_kelvin() {
        let glass = Substance::borosilicate_crown();
        let volume = Volume::from_si(std::f64::consts::PI * (0.0125f64).powi(2) * 0.005);
        let mass = glass.mass_of(volume);
        assert!((mass.to_si() * 1e3 - 6.16).abs() < 0.05, "{mass:?}");
        let capacity = glass.heat_capacity(volume).unwrap();
        assert!((capacity.to_si() - 5.28).abs() < 0.05, "{capacity:?}");
    }

    /// Thermal expansion, and the reason a bonded lens cracks: constrained
    /// stress is `E α ΔT` and does not depend on the size of the part, so a 60 K
    /// rise breaks N-BK7 whatever shape it is in.
    #[test]
    fn constrained_expansion_breaks_glass_before_metal() {
        let glass = Substance::borosilicate_crown();
        let metal = Substance::aluminium_6061();

        // 20 K over 100 mm of glass is 14 micrometres — small, and far more than
        // a wavelength.
        let growth = glass
            .expansion_of(Length::mm(100.0), Temperature::from_si(20.0))
            .unwrap();
        assert!((growth.in_um() - 14.2).abs() < 0.1, "{growth:?}");

        // Held rigidly, that same 20 K is 11.6 MPa: survivable.
        let stress = glass
            .constrained_stress(Temperature::from_si(20.0))
            .unwrap();
        assert!((stress.to_si() / 1e6 - 11.6).abs() < 0.2, "{stress:?}");
        assert_eq!(glass.survives(Temperature::from_si(20.0)), Some(true));
        // 120 K is not.
        assert_eq!(glass.survives(Temperature::from_si(120.0)), Some(false));
        // The aluminium mount takes it easily despite expanding three times more,
        // because it yields at 276 MPa rather than fracturing at 60.
        assert_eq!(metal.survives(Temperature::from_si(120.0)), Some(true));
    }

    /// A property that is not known reports that, rather than defaulting to a
    /// plausible number that would be silently wrong.
    #[test]
    fn unknown_properties_are_absent_not_guessed() {
        let unknown = Substance::bulk("unobtainium", Density::g_per_cm3(19.0));
        assert_eq!(unknown.diffusivity(), None);
        assert_eq!(unknown.heat_capacity(Volume::from_si(1e-6)), None);
        assert_eq!(unknown.survives(Temperature::from_si(50.0)), None);
        // But what *is* known still works.
        assert!((unknown.mass_of(Volume::from_si(1e-6)).to_si() - 0.019).abs() < 1e-9);
        // Water has no mechanical properties, and asking gives None rather than
        // an answer about a Young's modulus it does not have.
        assert_eq!(
            Substance::water().constrained_stress(Temperature::from_si(10.0)),
            None
        );
        assert!(Substance::water().diffusivity().is_some());
    }

    #[test]
    fn substances_round_trip_through_json() {
        let glass = Substance::borosilicate_crown();
        let json = serde_json::to_string(&glass).unwrap();
        assert_eq!(serde_json::from_str::<Substance>(&json).unwrap(), glass);
        // Absent properties are omitted rather than serialised as null.
        let plain = Substance::bulk("x", Density::kg_per_m3(1.0));
        let json = serde_json::to_string(&plain).unwrap();
        assert!(!json.contains("thermal"), "{json}");
    }
}
