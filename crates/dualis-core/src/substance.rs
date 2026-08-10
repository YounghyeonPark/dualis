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
///
/// # The two fields worth checking against the real part
///
/// `specific_heat` and `emissivity` are where a wrong number does the most damage, and they are
/// the two a user is most likely to carry over from something that looked close enough.
///
/// **A composite assembly is not a billet of its main metal.** A BLDC motor is copper, electrical
/// steel, magnets and air; its bulk `c_p` is nearer 450 J/kg/K than aluminium's 896. Reaching for
/// [`Substance::aluminium_6061`] because it is the metal in the catalogue **doubles the thermal
/// time constant** and changes the conclusion, with nothing to warn you — the answer stays
/// plausible, it is just for a different object. Use [`Substance::with_specific_heat`] on
/// whichever entry is closest and put the real figure in.
///
/// **Emissivity is a surface, not a substance.** The same 6061 is 0.09 polished and about 0.9
/// anodised, a factor of ten in the radiative path — which
/// [`Environment::loss_from`](../../dualis_thermal/struct.Environment.html) says is the same order
/// as still-air convection at room temperature. [`Substance::with_emissivity`] exists so a finish
/// does not have to become a new material.
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

    /// The same substance with a different surface finish.
    ///
    /// Emissivity is a property of the surface and not of the material, so anodised aluminium is
    /// not a new entry in the catalogue — it is `aluminium_6061().with_emissivity(0.9)`. The
    /// factor of ten between polished and anodised 6061 lands squarely on the radiative loss
    /// path, which is the same order as still-air convection at room temperature.
    ///
    /// Clamped to `0..=1`: a surface cannot radiate more than a blackbody, and a negative
    /// emissivity would make a body warm itself.
    ///
    /// Does nothing to a substance whose thermal properties are unknown, because `None` means
    /// unknown rather than zero and inventing three of the four fields to set the fourth would
    /// be worse than declining.
    pub fn with_emissivity(mut self, emissivity: f64) -> Substance {
        if let Some(t) = self.thermal.as_mut() {
            t.emissivity = emissivity.clamp(0.0, 1.0);
        }
        self
    }

    /// The same substance with a different heat capacity.
    ///
    /// For the assembly case: a motor, a populated board, a printed part with infill. The bulk
    /// `c_p` of a mixture is not the `c_p` of its main constituent, and this is the field where
    /// that difference is worth a factor of two.
    ///
    /// Does nothing to a substance whose thermal properties are unknown, for the reason in
    /// [`Substance::with_emissivity`].
    pub fn with_specific_heat(mut self, specific_heat: SpecificHeat) -> Substance {
        if let Some(t) = self.thermal.as_mut() {
            t.specific_heat = specific_heat;
        }
        self
    }

    /// Austenitic stainless, 304/18-8. What a portafilter basket, a boiler and most food-contact
    /// hardware is.
    ///
    /// # It is a poor conductor and that is the point
    ///
    /// 16.2 W/m/K against aluminium's 167 — a factor of ten — while holding **more** heat per unit
    /// volume, 4.0 MJ/m³/K against 2.4. So a steel part is a better *reservoir* and a worse
    /// *spreader* than an aluminium one of the same size, which is why a group head is brass and
    /// a basket is not.
    ///
    /// For an explicit conduction solver that combination is also the difference between a step of
    /// 2.4 ms and one of 41 ms on a millimetre grid, because the limit goes as the diffusivity and
    /// steel's is seventeen times lower. Reaching for aluminium because it is the metal already in
    /// the catalogue costs an order of magnitude in run time *and* understates the thermal mass.
    pub fn stainless_304() -> Substance {
        Substance {
            name: "304 stainless".to_string(),
            density: Density::g_per_cm3(8.00),
            thermal: Some(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(16.2),
                specific_heat: SpecificHeat::j_per_kg_k(500.0),
                expansion: ThermalExpansion::ppm_per_k(17.3),
                // Rolled and passivated rather than mirror-polished, which is what a basket is.
                emissivity: 0.28,
            }),
            mechanical: Some(MechanicalProps {
                youngs_modulus: Pressure::from_si(193.0e9),
                poisson_ratio: 0.29,
                yield_strength: Pressure::from_si(215.0e6),
            }),
            acoustic: Some(AcousticProps {
                sound_speed: Velocity::m_per_s(5_790.0),
            }),
        }
    }

    /// Electrolytic tough-pitch copper: windings, heat spreaders, planes.
    ///
    /// The values are uncontroversial to three figures. The emissivity is not: this is **bright
    /// polished** copper at 0.04, and copper oxidises — a tarnished surface runs 0.4 to 0.8, a
    /// factor of fifteen on the radiative path. If the part has been in air for a week, say so
    /// with [`Substance::with_emissivity`].
    pub fn copper() -> Substance {
        Substance {
            name: "Cu ETP".to_string(),
            density: Density::g_per_cm3(8.96),
            thermal: Some(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(401.0),
                specific_heat: SpecificHeat::j_per_kg_k(385.0),
                expansion: ThermalExpansion::ppm_per_k(16.5),
                emissivity: 0.04,
            }),
            mechanical: Some(MechanicalProps {
                youngs_modulus: Pressure::from_si(117.0e9),
                poisson_ratio: 0.34,
                yield_strength: Pressure::from_si(70.0e6),
            }),
            acoustic: Some(AcousticProps {
                sound_speed: Velocity::m_per_s(4_760.0),
            }),
        }
    }

    /// FR-4 glass-epoxy laminate: the board a driver sits on.
    ///
    /// **The conductivity is the through-plane one**, 0.3 W/m/K, and that is the number a
    /// designer wants because it is the one heat has to cross to reach the far side. In-plane it
    /// is nearer 0.8, because the copper-free glass weave conducts better along its fibres —
    /// a factor of about three, and `ThermalProps` carries one scalar, so the choice has to be
    /// stated rather than averaged. Any real board is dominated by its copper pour anyway, which
    /// is not laminate at all.
    ///
    /// The expansion is likewise in-plane, 14 ppm/K. Through-thickness FR-4 expands four to five
    /// times faster and goes higher again above its glass transition, which is what breaks
    /// plated through-holes; that regime is not modelled here.
    pub fn fr4() -> Substance {
        Substance {
            name: "FR-4".to_string(),
            density: Density::g_per_cm3(1.85),
            thermal: Some(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(0.30),
                specific_heat: SpecificHeat::j_per_kg_k(1_100.0),
                expansion: ThermalExpansion::ppm_per_k(14.0),
                emissivity: 0.90,
            }),
            mechanical: Some(MechanicalProps {
                youngs_modulus: Pressure::from_si(22.0e9),
                poisson_ratio: 0.16,
                yield_strength: Pressure::from_si(300.0e6),
            }),
            acoustic: None,
        }
    }

    /// Non-oriented silicon electrical steel: motor and transformer laminations.
    ///
    /// **Grade-dependent, and the spread is wide.** Silicon content trades core loss against
    /// conductivity: 25 W/m/K here is mid-range for non-oriented sheet, and grades run from about
    /// 20 to 30. Stacked laminations conduct far worse *across* the stack than the sheet does,
    /// because the interlaminar varnish dominates — a stack is not this substance at all, and
    /// treating it as one overstates the conduction out of a motor.
    ///
    /// Emissivity 0.3 is varnished sheet; bare mill finish is lower and rusty is much higher.
    pub fn electrical_steel() -> Substance {
        Substance {
            name: "electrical steel (non-oriented)".to_string(),
            density: Density::g_per_cm3(7.65),
            thermal: Some(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(25.0),
                specific_heat: SpecificHeat::j_per_kg_k(460.0),
                expansion: ThermalExpansion::ppm_per_k(12.0),
                emissivity: 0.30,
            }),
            mechanical: Some(MechanicalProps {
                youngs_modulus: Pressure::from_si(200.0e9),
                poisson_ratio: 0.29,
                yield_strength: Pressure::from_si(350.0e6),
            }),
            acoustic: Some(AcousticProps {
                sound_speed: Velocity::m_per_s(5_100.0),
            }),
        }
    }

    /// Solid cast PLA: printed structure, if it were solid, which it is not.
    ///
    /// **A printed part is not this substance.** Infill and layer adhesion move the effective
    /// conductivity and density more than the polymer chemistry does: at 20% infill the density
    /// is a fifth of this and the through-layer conductivity is lower again, because the path
    /// crosses voids and weld lines rather than bulk. Scale the density by the infill fraction
    /// at the very least, and treat the conductivity as an upper bound.
    ///
    /// That is not a caveat about precision. It is the difference between a part that survives
    /// and one that creeps: PLA softens around 60 °C, and [`Substance::survives`] is checking
    /// against a number the print may not reach in practice.
    pub fn pla() -> Substance {
        Substance {
            name: "PLA (solid)".to_string(),
            density: Density::g_per_cm3(1.24),
            thermal: Some(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(0.13),
                specific_heat: SpecificHeat::j_per_kg_k(1_800.0),
                expansion: ThermalExpansion::ppm_per_k(70.0),
                emissivity: 0.90,
            }),
            mechanical: Some(MechanicalProps {
                youngs_modulus: Pressure::from_si(3.5e9),
                poisson_ratio: 0.36,
                yield_strength: Pressure::from_si(50.0e6),
            }),
            acoustic: None,
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
    /// The builders change one field and leave the rest alone.
    ///
    /// Emissivity is a surface and not a substance, so anodised 6061 has to be reachable without
    /// a second catalogue entry — that was the reported friction, and the workaround was reaching
    /// into `thermal.as_mut()` by hand.
    #[test]
    fn a_finish_is_not_a_new_material() {
        let polished = Substance::aluminium_6061();
        let anodised = Substance::aluminium_6061().with_emissivity(0.9);
        let (p, a) = (polished.thermal.unwrap(), anodised.thermal.unwrap());

        assert_eq!(p.emissivity, 0.09);
        assert_eq!(a.emissivity, 0.9);
        // Everything else survives, including the name: it is the same alloy.
        assert_eq!(p.conductivity, a.conductivity);
        assert_eq!(p.specific_heat, a.specific_heat);
        assert_eq!(p.expansion, a.expansion);
        assert_eq!(polished.density, anodised.density);
        assert_eq!(polished.name, anodised.name);

        // A surface cannot out-radiate a blackbody, nor warm itself.
        assert_eq!(
            Substance::aluminium_6061()
                .with_emissivity(4.0)
                .thermal
                .unwrap()
                .emissivity,
            1.0
        );
        assert_eq!(
            Substance::aluminium_6061()
                .with_emissivity(-1.0)
                .thermal
                .unwrap()
                .emissivity,
            0.0
        );
    }

    /// **The trap the docs now warn about, as a number.**
    ///
    /// A lumped time constant is `C/(hA)` and `C` is `rho V c_p`, so reaching for aluminium's
    /// 896 J/kg/K to stand in for a motor's ~450 doubles it. That was the reported failure: the
    /// catalogue offered exactly one metal, reaching for it was the reasonable thing to do, and
    /// it changed a conclusion with nothing to say so.
    ///
    /// Asserted on the ratio rather than on either value, because the ratio is the claim.
    #[test]
    fn the_specific_heat_a_user_borrows_is_worth_a_factor_of_two() {
        let volume = Volume::from_si(3.456e-4);
        let billet = Substance::aluminium_6061();
        let assembly =
            Substance::aluminium_6061().with_specific_heat(SpecificHeat::j_per_kg_k(450.0));

        let c_billet = billet.heat_capacity(volume).unwrap().to_si();
        let c_assembly = assembly.heat_capacity(volume).unwrap().to_si();
        let ratio = c_billet / c_assembly;
        assert!(
            (ratio - 896.0 / 450.0).abs() < 1e-12,
            "the capacity ratio is the specific-heat ratio: {ratio}"
        );
        assert!(ratio > 1.9, "a borrowed c_p is worth about two: {ratio}");
    }

    /// The new entries carry heat capacity and expansion, and their ordering is the physics.
    ///
    /// Not asserting the values against themselves — that would check nothing. The orderings
    /// are the claims: copper conducts far better than steel and steel far better than laminate
    /// and plastic; a polymer expands several times faster than a metal; and copper stores less
    /// heat per kilogram than aluminium while storing more per unit volume, which is why a heat
    /// spreader is copper and a heatsink is aluminium.
    #[test]
    fn the_new_entries_are_ordered_the_way_the_physics_is() {
        let cu = Substance::copper().thermal.unwrap();
        let steel = Substance::electrical_steel().thermal.unwrap();
        let fr4 = Substance::fr4().thermal.unwrap();
        let pla = Substance::pla().thermal.unwrap();
        let al = Substance::aluminium_6061().thermal.unwrap();

        // Conduction, over four orders of magnitude.
        assert!(cu.conductivity > al.conductivity);
        assert!(al.conductivity > steel.conductivity);
        assert!(steel.conductivity.to_si() > 50.0 * fr4.conductivity.to_si());
        assert!(fr4.conductivity > pla.conductivity);

        // Expansion: a polymer moves about three times faster than the fastest metal here.
        // 70 ppm/K against 6061's 23.6 is 2.97, and the first version of this asserted 3.0 --
        // a claim written from the adjective rather than from the numbers.
        let ratio = pla.expansion.to_si() / al.expansion.to_si();
        assert!(
            (2.5..3.5).contains(&ratio),
            "PLA against 6061 is {ratio:.2}x"
        );
        assert!(al.expansion > cu.expansion && cu.expansion > steel.expansion);

        // Per kilogram copper stores less than aluminium; per unit volume it stores more.
        let v = Volume::from_si(1e-3);
        assert!(cu.specific_heat < al.specific_heat);
        assert!(
            Substance::copper().heat_capacity(v).unwrap()
                > Substance::aluminium_6061().heat_capacity(v).unwrap()
        );

        // The insulators are the emitters, which is why a black plastic case sheds heat a bare
        // metal one does not.
        assert!(fr4.emissivity > 0.8 && pla.emissivity > 0.8);
        assert!(cu.emissivity < 0.1);
    }
}
