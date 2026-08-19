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
//! optics exists — refractive index is `pantometry-optics`'s
//! `Material`, and a consumer that needs both pairs them. Putting it here would
//! make the kernel depend on a domain, which is the one structural rule the split
//! exists to enforce.

use pantometry_units::{
    Density, Diffusivity, Energy, HeatCapacity, LatentHeat, Length, Mass, Pressure, SpecificHeat,
    Temperature, ThermalConductivity, ThermalExpansion, Velocity, Volume,
};
use serde::{Deserialize, Serialize};

/// A material, as much of it as is known.
///
/// # A key this type does not know is refused, not dropped
///
/// `deny_unknown_fields`, here and on all four property blocks. `serde` discards unknown keys by
/// default, which is right for a wire protocol that must tolerate a newer peer and wrong for a
/// material somebody wrote down: a mistyped `"thermalz"` would leave the whole thermal block absent
/// and the substance would run as one whose conductivity is *unknown* rather than as one whose file
/// has a typo in it.
///
/// The same rule `pantometry-world`'s scene format has, for the same reason, and it was added after a
/// test asked whether a typo was caught and found that it was not.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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
    /// What it takes to melt it, if it is a substance that melts at a temperature.
    ///
    /// Absent for most of them, and absent is not zero — it means the substance is being modelled
    /// as never changing phase, which is the right model for a heat sink and the wrong one for ice.
    /// A domain that finds it absent does not change phase; one that finds it present must account
    /// for the latent heat or its books will not balance.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fusion: Option<FusionProps>,
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
/// [`Environment::loss_from`](../../pantometry_thermal/struct.Environment.html) says is the same order
/// as still-air convection at room temperature. [`Substance::with_emissivity`] exists so a finish
/// does not have to become a new material.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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
#[serde(deny_unknown_fields)]
pub struct MechanicalProps {
    /// Young's modulus.
    pub youngs_modulus: Pressure,
    /// Poisson's ratio, dimensionless — how much it bulges sideways when squeezed.
    pub poisson_ratio: f64,
    /// Where it stops coming back. For a brittle material this is the fracture
    /// stress, and there is no plastic region before it.
    pub yield_strength: Pressure,
}

/// What it takes to melt it.
///
/// # One temperature, and the substances that do not have one
///
/// A pure substance melts at a temperature; an alloy, a polymer and a rock melt over a *range*, and
/// this cannot say so. That is a real restriction rather than a simplification to be embarrassed
/// about — the sharp-interface problem is the one with an exact solution to check against, and a
/// mushy range is a different model with a different closed form.
///
/// So this is right for water, for a pure metal and for a paraffin phase-change material sold on
/// its plateau. It is wrong for solder, and a domain given it for solder will put the whole latent
/// heat on one temperature instead of spreading it over the twenty kelvin it really occupies.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FusionProps {
    /// The temperature at which it changes phase, and holds there while it does.
    pub melting_point: Temperature,
    /// The heat one kilogram absorbs melting, at no change in temperature.
    pub latent_heat: LatentHeat,
    /// What the **liquid** conducts and holds, if it differs from the solid.
    ///
    /// `None` is the **one-phase** model: the liquid is taken to have the solid's conductivity and
    /// specific heat. That is not a simplification to apologise for — it is exact whenever the liquid
    /// sits at the melting point, because then no heat flows through it whatever its properties are,
    /// and it is the case Stefan's original problem and Neumann's solution are about.
    ///
    /// It is wrong the moment the liquid is **superheated**, and wrong by a lot. Water conducts a
    /// quarter of what ice does and holds twice as much, and a liquid 20 K above freezing slows a
    /// front by **16%** — from 15.85 mm to 13.33 mm at 900 s. That is far more than any
    /// discretisation error, so the one-phase answer is not a slightly worse two-phase answer.
    ///
    /// The same type as the solid's, because a phase is a thing that conducts and holds heat and there
    /// is no reason to describe it differently. `expansion` and `emissivity` are carried and unused by
    /// conduction; give the liquid's if they are known.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub liquid: Option<ThermalProps>,
}

impl FusionProps {
    /// A one-phase description: a melting point and a latent heat, and no separate liquid.
    ///
    /// # Because adding `liquid` broke every literal, twice over
    ///
    /// `Substance` gained `fusion` and every struct literal outside this crate stopped compiling;
    /// builders were added so the next field would not do it again. Then `FusionProps` gained `liquid`
    /// and did exactly that one level down, to the tests written the week before.
    ///
    /// So this pair exists for the same reason `Substance::with_*` does. A field added below here
    /// costs nothing to a caller who went through `new` and [`with_liquid`](FusionProps::with_liquid).
    pub fn new(melting_point: Temperature, latent_heat: LatentHeat) -> FusionProps {
        FusionProps {
            melting_point,
            latent_heat,
            liquid: None,
        }
    }

    /// Name the liquid phase's conductivity and specific heat, making a block **two-phase**.
    ///
    /// Read [`liquid`](FusionProps::liquid) before reaching for this: it is the right model for a
    /// superheated liquid and it costs first-order accuracy at the interface, so it is not a strictly
    /// better version of the one-phase model.
    pub fn with_liquid(mut self, liquid: ThermalProps) -> FusionProps {
        self.liquid = Some(liquid);
        self
    }

    /// How many kelvin of sensible heat the phase change is worth: `L / c_p`.
    ///
    /// The reciprocal of the Stefan number, and the number that says whether latent heat matters at
    /// all in a given problem. For ice it is **163 K**, so a freezing front driven by a 10 K
    /// undercooling is overwhelmingly a latent-heat problem and only incidentally a conduction one.
    pub fn sensible_equivalent(&self, specific_heat: SpecificHeat) -> Temperature {
        Temperature::from_si(self.latent_heat.to_si() / specific_heat.to_si())
    }
}

/// What it does with sound.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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

    /// The joules a given volume absorbs changing phase, or `None` if it does not.
    ///
    /// The companion to [`heat_capacity`](Substance::heat_capacity), and the pair is what a domain
    /// needs to keep books that balance across a melting front: one buys kelvin and the other buys
    /// none. For a cubic millimetre of ice they are 1.88 mJ/K and **306 mJ** — so the phase change
    /// is worth 163 K of warming, and a scheme that dropped it would run the front 163 times too
    /// fast rather than slightly wrong.
    pub fn latent_energy(&self, volume: Volume) -> Option<Energy> {
        let f = self.fusion?;
        Some(self.mass_of(volume) * f.latent_heat)
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

    /// Every substance this crate ships, by the short name a *file* writes.
    ///
    /// The constructors below are the Rust door and are enough for a caller who knows at compile
    /// time what the thing is made of. This is the other door: a name that arrived as text, from
    /// JSON, a command line, or a spreadsheet column.
    ///
    /// # Why a slug and not [`Substance::name`]
    ///
    /// `name` is free text meant for a human reading a violation message — `"Al 6061"`,
    /// `"N-BK7"`. A key that a file is matched against has to be stable, lowercase and
    /// unpunctuated, and the two are different jobs: renaming `"Al 6061"` to `"Aluminium 6061-T6"`
    /// improves one and breaks every file that used the other.
    ///
    /// # This list and [`Substance::from_name`] are checked against each other
    ///
    /// Both directions, in `any_material.rs`: every slug here resolves, and every constructor is
    /// reachable through some slug here. A pair of hand-written lists that agree until they do not
    /// is how a catalogue grows an entry that exists and cannot be named — which is what had
    /// happened to `water`, present in this crate since 0.1.0 and unreachable from a scene file for
    /// **eleven** releases — every version in which a scene could name a material at all — because
    /// the scene format kept its own eight-name copy of this.
    pub const CATALOGUE: [&'static str; 9] = [
        "aluminium",
        "borosilicate",
        "copper",
        "electrical_steel",
        "fr4",
        "ice",
        "pla",
        "stainless_304",
        "water",
    ];

    /// Look one up by the name in [`Substance::CATALOGUE`], or `None`.
    ///
    /// `None` rather than a panic or a default: a name that arrived as text is a name that can be
    /// wrong, and the caller is the one who knows what to say about it and where the text came from.
    /// Substituting a plausible material for an unrecognised name is the failure this whole file is
    /// arranged against.
    ///
    /// **A catalogue is not the answer to "any material".** Nine entries cannot be, and adding a
    /// tenth does not change that. The answer is that a `Substance` is data — `Deserialize` and
    /// [`Substance::check`] — so anything with a datasheet can be declared without this crate
    /// learning it exists. This function is for the common case where the material is one of nine
    /// and typing out its properties would be worse.
    pub fn from_name(name: &str) -> Option<Substance> {
        Some(match name {
            "aluminium" => Substance::aluminium_6061(),
            "borosilicate" => Substance::borosilicate_crown(),
            "copper" => Substance::copper(),
            "electrical_steel" => Substance::electrical_steel(),
            "fr4" => Substance::fr4(),
            // `ice` is the only entry that changes phase, and so the only one for which a domain
            // reports a melted volume. `water` is the *liquid* — no `fusion`, because a substance
            // whose fusion is present is one being modelled as freezable, and a coolant mass in a
            // network is not. The pair is deliberate and they are not interchangeable.
            "ice" => Substance::ice(),
            "water" => Substance::water(),
            "pla" => Substance::pla(),
            "stainless_304" => Substance::stainless_304(),
            _ => return None,
        })
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
            fusion: None,
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
            fusion: None,
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
            fusion: None,
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
            fusion: None,
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
            fusion: None,
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
            fusion: None,
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
            fusion: None,
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
            fusion: None,
            mechanical: None,
            acoustic: Some(AcousticProps {
                sound_speed: Velocity::m_per_s(1_482.0),
            }),
        }
    }

    /// Ice at 0 °C, and the only entry in this catalogue that changes phase.
    ///
    /// The canonical Stefan material, and the numbers are the ones the closed-form tests need.
    /// 2.22 W/m·K is **four times** liquid water's 0.598, which is the thing about ice that surprises
    /// people and the reason a lake freezes downward at all.
    ///
    /// # This is the solid, and the one-phase model uses it for both sides
    ///
    /// A domain given this for a melting problem is taking the liquid's conductivity and specific
    /// heat to be the solid's, which they are not — water conducts a quarter as well and holds twice
    /// as much. That is **Stefan's original one-phase problem**, and it is exact when the liquid is
    /// already at the melting point so no heat flows through it: a lake freezing from a cold sky,
    /// where all the resistance is in the ice.
    ///
    /// It is not right for melting a block of ice into water that then warms up. Use
    /// [`Substance::water`] for the liquid and note that a cell cannot currently be both.
    pub fn ice() -> Substance {
        Substance {
            name: "ice".to_string(),
            density: Density::g_per_cm3(0.917),
            thermal: Some(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(2.22),
                specific_heat: SpecificHeat::j_per_kg_k(2_050.0),
                // Ice's expansion is anisotropic and this is the polycrystalline mean.
                expansion: ThermalExpansion::ppm_per_k(51.0),
                emissivity: 0.97,
            }),
            fusion: Some(FusionProps {
                melting_point: Temperature::celsius(0.0),
                latent_heat: LatentHeat::kj_per_kg(333.55),
                // **One-phase, deliberately, and this was measured before being decided.**
                //
                // Giving this entry water as its liquid made the one-phase answer *worse*: the front
                // in `a_freezing_front.rs` went from 0.43% out to 6.9% at forty cells, sixteen times
                // worse for a problem whose physics had not changed. The cause is the mushy cell,
                // whose conductivity is a mixture — so the cell holding the interface conducts partly
                // like water, and the heat reaching the interface has less conductance than the exact
                // solution gives it. That is a first-order error at the front and it is the price of
                // two phases.
                //
                // A default that silently costs a caller a factor of sixteen is the wrong default. Two
                // phases are opt-in: `Substance::ice().with_fusion(FusionProps { liquid: Some(..), .. })`
                // and `crates/pantometry-thermal/tests/two_phase_stefan.rs` shows it against the two-phase
                // Neumann solution, where it is the *right* answer and the one-phase model is 16% out.
                liquid: None,
            }),
            mechanical: Some(MechanicalProps {
                youngs_modulus: Pressure::from_si(9.1e9),
                poisson_ratio: 0.33,
                // Tensile strength, and ice is brittle: there is no yield before it.
                yield_strength: Pressure::from_si(1.0e6),
            }),
            acoustic: Some(AcousticProps {
                sound_speed: Velocity::m_per_s(3_840.0),
            }),
        }
    }

    /// Give it thermal properties, or replace the ones it has.
    ///
    /// # Why builders exist, when the fields are already public
    ///
    /// Because a struct literal names **every** field, so it breaks the moment this type learns one.
    /// `fusion` was added for latent heat and every literal outside this crate stopped compiling —
    /// allowed in `0.x`, and still a cost paid by exactly the callers this catalogue is least able to
    /// help: the ones whose material is not in it.
    ///
    /// A chain of `with_*` on [`bulk`](Substance::bulk) is immune to that, and it is how any real
    /// material becomes expressible without waiting for it to be added here:
    ///
    /// ```
    /// # use pantometry_core::substance::{MechanicalProps, Substance, ThermalProps};
    /// # use pantometry_core::units::*;
    /// // Ti-6Al-4V, from a datasheet rather than from this crate.
    /// let titanium = Substance::bulk("Ti-6Al-4V", Density::g_per_cm3(4.43))
    ///     .with_thermal(ThermalProps {
    ///         conductivity: ThermalConductivity::w_per_m_k(6.7),
    ///         specific_heat: SpecificHeat::j_per_kg_k(526.0),
    ///         expansion: ThermalExpansion::ppm_per_k(8.6),
    ///         emissivity: 0.30,
    ///     })
    ///     .with_mechanical(MechanicalProps {
    ///         youngs_modulus: Pressure::from_si(113.8e9),
    ///         poisson_ratio: 0.342,
    ///         yield_strength: Pressure::from_si(880.0e6),
    ///     });
    /// assert!(titanium.check().is_ok());
    /// ```
    ///
    /// **Enumeration does not reach "every material" and data does.** This catalogue holds nine
    /// entries because each is a set of numbers somebody has to be answerable for; a caller with a
    /// datasheet is answerable for theirs. [`check`](Substance::check) is what the library can still
    /// do for them.
    pub fn with_thermal(mut self, thermal: ThermalProps) -> Substance {
        self.thermal = Some(thermal);
        self
    }

    /// Give it mechanical properties, or replace the ones it has.
    pub fn with_mechanical(mut self, mechanical: MechanicalProps) -> Substance {
        self.mechanical = Some(mechanical);
        self
    }

    /// Give it acoustic properties, or replace the ones it has.
    pub fn with_acoustic(mut self, acoustic: AcousticProps) -> Substance {
        self.acoustic = Some(acoustic);
        self
    }

    /// Give it a phase change, or replace the one it has.
    pub fn with_fusion(mut self, fusion: FusionProps) -> Substance {
        self.fusion = Some(fusion);
        self
    }

    /// Every problem with this substance's numbers, or `Ok` if there are none.
    ///
    /// For a material that came from outside this crate — a datasheet, a JSON file, a builder chain —
    /// where nobody has checked the numbers against anything. It cannot tell whether a conductivity is
    /// *right*; it can tell whether it is **possible**, and an impossible one otherwise produces an
    /// answer that is plausible and wrong.
    ///
    /// Reports all of them at once rather than the first, because a material transcribed from the
    /// wrong column is usually wrong in several places.
    ///
    /// # The one check that is not a bound on a single field
    ///
    /// If a substance states both a sound speed and elastic constants, those are **three independent
    /// numbers describing one thing**, and they have to agree. A longitudinal wave in a solid is
    /// bounded below by the rod speed `sqrt(E/rho)` — free to bulge sideways — and above by the bulk
    /// speed `sqrt((lambda+2mu)/rho)`, fully constrained, so a stated speed must sit near one of them.
    ///
    /// **15%**, and it is measured rather than chosen: across this catalogue every entry is within
    /// 6.2% of whichever it means, and the gap is there because a tensile test and an ultrasonic
    /// measurement are not the same measurement — read as a bulk wave, copper's stated speed implies
    /// 132 GPa against the 117 in its own entry. So the bound cannot be tighter than that, and at 15%
    /// it still catches a **shear** speed transcribed by mistake, which sits 45% below the rod speed.
    pub fn check(&self) -> Result<(), String> {
        let mut wrong: Vec<String> = Vec::new();
        let positive = |what: &str, v: f64, out: &mut Vec<String>| {
            if !(v.is_finite() && v > 0.0) {
                out.push(format!("{what} must be finite and positive, is {v}"));
            }
        };
        positive("density", self.density.to_si(), &mut wrong);
        if let Some(t) = self.thermal {
            positive("conductivity", t.conductivity.to_si(), &mut wrong);
            positive("specific_heat", t.specific_heat.to_si(), &mut wrong);
            if !(0.0..=1.0).contains(&t.emissivity) {
                wrong.push(format!(
                    "emissivity is a fraction of a blackbody's and must be in 0..=1, is {}",
                    t.emissivity
                ));
            }
            if !t.expansion.to_si().is_finite() {
                wrong.push("expansion must be finite".to_string());
            }
        }
        if let Some(m) = self.mechanical {
            positive("youngs_modulus", m.youngs_modulus.to_si(), &mut wrong);
            positive("yield_strength", m.yield_strength.to_si(), &mut wrong);
            // The same range `pantometry-elastic` refuses outside of: at one half the material is
            // incompressible and lambda is infinite, at minus one the shear modulus diverges.
            if !(-1.0 < m.poisson_ratio && m.poisson_ratio < 0.5) {
                wrong.push(format!(
                    "poisson_ratio must be in (-1, 0.5) for a stable isotropic solid, is {}",
                    m.poisson_ratio
                ));
            }
        }
        if let Some(a) = self.acoustic {
            positive("sound_speed", a.sound_speed.to_si(), &mut wrong);
        }
        if let Some(f) = self.fusion {
            positive("latent_heat", f.latent_heat.to_si(), &mut wrong);
            positive("melting_point", f.melting_point.to_si(), &mut wrong);
        }
        // The cross-check, and only when there is something to cross.
        if let (Some(m), Some(a)) = (self.mechanical, self.acoustic) {
            let (e, nu, rho) = (
                m.youngs_modulus.to_si(),
                m.poisson_ratio,
                self.density.to_si(),
            );
            if e > 0.0 && rho > 0.0 && -1.0 < nu && nu < 0.5 {
                let rod = (e / rho).sqrt();
                let bulk = (e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu) * rho)).sqrt();
                let c = a.sound_speed.to_si();
                let gap = (c / rod - 1.0).abs().min((c / bulk - 1.0).abs());
                if gap > 0.15 {
                    wrong.push(format!(
                        "sound_speed {c:.0} m/s is {:.0}% from the nearer of the rod speed {rod:.0} \
                         and the bulk speed {bulk:.0} that its own E, nu and density give — so one of \
                         the four is not this material's, or the speed is a shear wave",
                        gap * 100.0
                    ));
                }
            }
        }
        if wrong.is_empty() {
            Ok(())
        } else {
            Err(format!("{}: {}", self.name, wrong.join("; ")))
        }
    }

    /// A substance with nothing known but how heavy it is.
    pub fn bulk(name: &str, density: Density) -> Substance {
        Substance {
            name: name.to_string(),
            density,
            thermal: None,
            fusion: None,
            mechanical: None,
            acoustic: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pantometry_units::Length;

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
