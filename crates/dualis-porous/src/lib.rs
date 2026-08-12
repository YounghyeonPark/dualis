#![deny(missing_docs)]

//! Flow through a packed bed, the heat it carries, and the dissolution it drives.
//!
//! A bed of ground solid with liquid forced through it under pressure. That is a coffee puck, and
//! it is also a chromatography column, a sand filter, a catalyst bed, a leaching heap and an
//! aquifer. What they share is the reason this is a domain rather than an example: on any scale
//! larger than a pore, momentum in a packed bed is not Navier–Stokes but **Darcy's law**, and
//! Darcy's law is elliptic.
//!
//! ```text
//!   u = −(k/μ) ∇p              Darcy: velocity is proportional to the pressure gradient
//!   ∇·u = 0                    incompressible
//!   ⇒  ∇·((k/μ) ∇p) = 0        an elliptic problem, solved rather than marched
//! ```
//!
//! That is the same operator [`Conductor`](https://docs.rs/dualis-electrical) solves for electric
//! potential, with mobility `k/μ` in place of conductivity. The reuse is not a coincidence or a
//! trick — a resistor network and a packed bed are the same mathematics, and getting the same
//! answer out of the same discretisation is one of the checks this crate makes.
//!
//! On top of the flow ride two things that are *not* elliptic and are marched:
//!
//! - **heat**, advected by the liquid and conducted through the bed, and
//! - **solute**, dissolving out of the particles and carried away by the liquid.
//!
//! # The three numbers a barista turns, and what each one is
//!
//! | knob | what it changes | how |
//! | --- | --- | --- |
//! | grind | permeability **and** extraction rate | `k ∝ d²`, and `K ∝ 1/d²` |
//! | temperature | extraction rate **and** viscosity | Arrhenius, and `μ(T)` |
//! | pressure | flow rate, linearly | `Q ∝ Δp` |
//!
//! The grind row is the interesting one and it is why espresso is hard. A finer grind makes the
//! bed **less** permeable, so the shot runs longer, *and* it makes each particle extract faster,
//! because a solute molecule has less distance to diffuse out. Both effects push yield the same
//! way, and they are separate physics with the same `d²` in them pointing in opposite directions.
//! Nothing here codes that relationship: it falls out of [`Grind::permeability`] and
//! [`Grind::extraction_rate`] being derived independently.
//!
//! # What is deliberately not modelled
//!
//! Stated because each of these is a real effect that a reader will otherwise assume is in here:
//!
//! - **Wetting.** The bed starts saturated. Pre-infusion is a two-phase problem — Richards'
//!   equation, a capillary pressure curve, a moving front — and a single-phase solver that
//!   pretends otherwise would report a first second of the shot that is fiction.
//! - **Compaction and swelling.** Porosity is fixed. A real puck compresses under 9 bar and the
//!   grounds swell as they wet, both of which change `k` during the shot.
//! - **Crema.** No dissolved CO₂, no gas phase, no bubbles.
//! - **Fines migration.** Particles do not move, so a puck cannot clog during its own shot.
//! - **The solute's effect on the liquid.** Viscosity and density are pure water's. At 10% TDS
//!   that is a real approximation, not a negligible one.
//!
//! # Determinism
//!
//! Nothing here is random and nothing consults a clock. The pressure solve is conjugate gradients
//! with a fixed iteration order from a fixed start; the marched parts are explicit and swept in
//! index order.

use dualis_units::{
    Area, Density, Diffusivity, DynamicViscosity, Length, SpecificHeat, Temperature,
    ThermalConductivity,
};

mod puck;

pub use puck::{Basket, Observable, Puck, PuckField, Shot};

/// The gas constant, J·mol⁻¹·K⁻¹.
const R_GAS: f64 = 8.314_462_618;

/// The liquid being pushed through the bed.
///
/// Temperature-dependent where it matters. Viscosity is the one that surprises people: water at
/// 70 °C is a third more viscous than at 100 °C, so a brew temperature change moves the flow rate
/// on its own, before any effect on extraction.
#[derive(Clone, Copy, Debug)]
pub struct Liquid {
    /// Density, treated as constant. Water changes 4% over the brewing range and the flow is
    /// pressure-driven, so this is the one property where the temperature dependence does not
    /// earn its complexity.
    pub density: Density,
    /// Specific heat.
    pub specific_heat: SpecificHeat,
    /// Thermal conductivity of the liquid itself.
    pub conductivity: ThermalConductivity,
    /// Andrade prefactor `A` in `μ = A·exp(B/T)`, Pa·s.
    viscosity_a: f64,
    /// Andrade exponent `B`, K.
    viscosity_b: f64,
}

impl Liquid {
    /// Water, as it is between about 60 °C and boiling.
    ///
    /// The viscosity is an Andrade fit, `μ = A·exp(B/T)`, with `A` and `B` set from the tabulated
    /// values at 70 °C and 93 °C. It is checked against 100 °C, which the fit was **not** given:
    /// it predicts 2.83×10⁻⁴ Pa·s against a tabulated 2.82×10⁻⁴, so the two-point fit extrapolates
    /// to within half a percent across the range anyone brews in.
    pub fn water() -> Liquid {
        Liquid {
            density: Density::from_si(965.0),
            specific_heat: SpecificHeat::from_si(4205.0),
            conductivity: ThermalConductivity::from_si(0.675),
            viscosity_a: 4.861e-6,
            viscosity_b: 1516.7,
        }
    }

    /// Dynamic viscosity at a temperature.
    pub fn viscosity(&self, t: Temperature) -> DynamicViscosity {
        DynamicViscosity::from_si(self.viscosity_a * (self.viscosity_b / t.to_si()).exp())
    }

    /// Volumetric heat capacity, `ρ·c`, J·m⁻³·K⁻¹.
    pub(crate) fn rho_c(&self) -> f64 {
        self.density.to_si() * self.specific_heat.to_si()
    }
}

/// The solid the bed is made of.
///
/// Coffee is the case this was written for, but nothing here is coffee-specific: a `Bed` is a
/// packed granular solid with a soluble fraction.
#[derive(Clone, Copy, Debug)]
pub struct Bed {
    /// **Apparent** density of a particle: its mass over the volume it occupies, including the
    /// pores inside it.
    ///
    /// Three densities get confused here and picking the wrong one is a factor of two in the dose:
    ///
    /// | | roasted coffee | what it is |
    /// | --- | --- | --- |
    /// | skeletal | ≈1400 kg/m³ | the cell-wall material, with no voids at all |
    /// | **apparent** | **≈600 kg/m³** | a particle, including the voids inside it |
    /// | bulk | ≈330 kg/m³ | a tamped bed, including the gaps between particles |
    ///
    /// This is the middle one, because [`Puck`]'s porosity is the *inter*-particle
    /// void and the two multiply to give the bulk. Using the skeletal 1400 with an inter-particle
    /// porosity of 0.45 puts 40 g in an 18 g basket — which is what the first version of this
    /// crate did, and it was visible only because a dose is a number a person recognises.
    pub solid_density: Density,
    /// Specific heat of the solid.
    pub specific_heat: SpecificHeat,
    /// Thermal conductivity of the solid **grain**, not of the bed.
    ///
    /// The bed's is [`Puck::bed_conductivity`], which combines this with the liquid's by Maxwell–Eucken
    /// with the liquid as the continuous phase. Read that before changing this one: it used to be an
    /// arithmetic mean, which is the Voigt bound and was 11.0% high, and the reason nobody noticed is
    /// that a bed under flow is isothermal so conduction carries nothing at all.
    pub conductivity: ThermalConductivity,
    /// The fraction of the dry solid mass that can dissolve at all.
    ///
    /// For roasted coffee this is about 0.30, and it is the ceiling every yield figure is
    /// measured against: a 20% extraction yield means two-thirds of what *could* dissolve did.
    pub soluble_fraction: f64,
    /// Effective diffusivity of the solute inside a particle at [`Bed::reference_temperature`].
    ///
    /// Not the bulk diffusivity of the solute in free water — an order of magnitude smaller,
    /// because the path out of a particle is tortuous and partly blocked.
    pub solute_diffusivity: Diffusivity,
    /// Activation energy for that diffusivity, J/mol.
    pub activation_energy: f64,
    /// The pore-liquid concentration in equilibrium with **fully loaded** grounds.
    ///
    /// # This is the term that makes channelling exist
    ///
    /// Without it, dissolution is `dm/dt = −K·m` and depends on nothing outside the particle: every
    /// cell at the same temperature extracts at the same rate no matter how much liquid passes it.
    /// A bed with a channel through it then extracts perfectly evenly, which is the opposite of
    /// what a channel does and is what the first version of this crate reported.
    ///
    /// With it, the driving force is a *difference*:
    ///
    /// ```text
    ///   dm/dt = −K·(m − m₀·C/C_sat)
    /// ```
    ///
    /// A cell the flow rushes past is kept at low `C` and keeps extracting; a cell the flow avoids
    /// fills up its own pore liquid and stalls. That is the whole mechanism, and it is why a
    /// channelled shot is simultaneously over- and under-extracted.
    ///
    /// 300 kg/m³ — about 30% by mass, which is where a coffee concentrate stops taking up more
    /// solids. It is a calibration in the same sense as the diffusivity: the *form* is a linear
    /// isotherm and the number is one measurement.
    pub saturation_concentration: Density,
    /// The temperature `solute_diffusivity` is quoted at.
    pub reference_temperature: Temperature,
}

impl Bed {
    /// Roasted coffee.
    ///
    /// The diffusivity is a **calibration**, and saying so is the point. It was set so that one
    /// conventional shot lands where a barista aims it — 17.6 g of 250 µm grind at `ε = 0.45` in a
    /// 58 mm basket, 9 bar at 93 °C, giving **38.6 g in 25 s at 19.8% yield and 8.3% TDS**. That
    /// is one measurement fixing one number.
    ///
    /// What is *not* calibrated is the form. The rate goes as `1/d²` because the solute diffuses
    /// out of a sphere, and it slows as the surrounding liquid fills up because the driving force
    /// is a difference. Both are checked against their closed forms rather than fitted.
    ///
    /// 2.0×10⁻¹⁰ m²/s is about a seventh of the bulk diffusivity of a small sugar in water at this
    /// temperature — a tortuosity of seven, which is ordinary for a porous particle. That is why
    /// the number is believable and not merely fitted to.
    pub fn roasted_coffee() -> Bed {
        Bed {
            solid_density: Density::from_si(600.0),
            specific_heat: SpecificHeat::from_si(1670.0),
            conductivity: ThermalConductivity::from_si(0.15),
            soluble_fraction: 0.30,
            solute_diffusivity: Diffusivity::from_si(2.0e-10),
            activation_energy: 30_000.0,
            saturation_concentration: Density::from_si(300.0),
            reference_temperature: Temperature::celsius(93.0),
        }
    }

    /// Volumetric heat capacity of the solid, J·m⁻³·K⁻¹.
    pub(crate) fn rho_c(&self) -> f64 {
        self.solid_density.to_si() * self.specific_heat.to_si()
    }

    /// The solute diffusivity at a temperature, by Arrhenius.
    pub fn diffusivity_at(&self, t: Temperature) -> Diffusivity {
        let (t, t_ref) = (t.to_si(), self.reference_temperature.to_si());
        let factor = (-self.activation_energy / R_GAS * (1.0 / t - 1.0 / t_ref)).exp();
        Diffusivity::from_si(self.solute_diffusivity.to_si() * factor)
    }
}

/// How finely the solid is ground, and the two quite different lengths that follow from it.
///
/// # Why one grind setting needs two diameters
///
/// **Extraction** is diffusion out of a particle, and the length that governs it is the particle's
/// own radius. The sieve diameter is the right number.
///
/// **Flow** is not. A real grind is not one size — it is a coarse mode plus a tail of fines, and
/// the fines lodge in the gaps between the coarse particles and carry most of the pressure drop. A
/// Kozeny–Carman permeability computed from the sieve diameter is wrong by four orders of
/// magnitude in the direction that matters: it predicts an espresso puck passes twenty litres a
/// second.
///
/// So the flow uses a **hydraulic diameter**, which is smaller. [`Grind::sieved`] takes the ratio
/// as [`Grind::FINES_RATIO`] — a calibration, stated as one, obtained by matching one measured
/// shot and used for nothing else. The `d²` scaling on either side is physics and is tested; the
/// single constant between them is not, and is why [`Grind::hydraulic`] exists for somebody who
/// has measured their own bed.
#[derive(Clone, Copy, Debug)]
pub struct Grind {
    sieve: f64,
    hydraulic: f64,
}

impl Grind {
    /// Hydraulic diameter as a fraction of the sieve diameter, for a typical espresso grind.
    ///
    /// One hundred and sixtieth. Fitted to a single measurement — a 58 mm basket 20 mm deep at
    /// `ε = 0.45`, which holds 17.4 g, with 9 bar across the bed at 93 °C, running at 1.5 g/s.
    /// That corresponds to a bed permeability of 4.1×10⁻¹⁵ m², which is the number actually being
    /// asserted; the diameter is how it is carried to other grinds.
    ///
    /// It is large, and that is the finding rather than an embarrassment: a Kozeny–Carman
    /// permeability computed from the sieve diameter says an espresso puck passes twenty litres a
    /// second. Two orders of magnitude of the resistance are in the fines and in the constrictions
    /// inside the particles, neither of which a sieve measures.
    ///
    /// It is a property of the *grinder's particle size distribution*, not of coffee: a grinder
    /// that makes fewer fines has a larger ratio and a puck that runs faster at the same setting,
    /// which is most of what people are buying when they buy a better burr set.
    pub const FINES_RATIO: f64 = 1.0 / 160.0;

    /// A grind quoted the way a sieve or a laser sizer quotes it.
    pub fn sieved(diameter: Length) -> Grind {
        Grind {
            sieve: diameter.to_si().max(f64::MIN_POSITIVE),
            hydraulic: diameter.to_si().max(f64::MIN_POSITIVE) * Grind::FINES_RATIO,
        }
    }

    /// A grind whose two diameters are both known — because the bed's permeability was measured
    /// rather than predicted.
    pub fn hydraulic(sieve: Length, hydraulic: Length) -> Grind {
        Grind {
            sieve: sieve.to_si().max(f64::MIN_POSITIVE),
            hydraulic: hydraulic.to_si().max(f64::MIN_POSITIVE),
        }
    }

    /// A conventional espresso grind: 250 µm.
    pub fn espresso() -> Grind {
        Grind::sieved(Length::from_si(250e-6))
    }

    /// The same distribution scaled by a factor. Below one is finer.
    pub fn scaled(&self, factor: f64) -> Grind {
        Grind {
            sieve: self.sieve * factor,
            hydraulic: self.hydraulic * factor,
        }
    }

    /// The sieve diameter — the particle size a grinder is set by.
    pub fn sieve_diameter(&self) -> Length {
        Length::from_si(self.sieve)
    }

    /// The hydraulic diameter the permeability is computed from.
    pub fn hydraulic_diameter(&self) -> Length {
        Length::from_si(self.hydraulic)
    }

    /// Kozeny–Carman permeability at a porosity.
    ///
    /// ```text
    ///   k = d² ε³ / (180 (1−ε)²)
    /// ```
    ///
    /// The `ε³/(1−ε)²` is the part people get wrong by remembering it as `ε³/(1−ε)`. Both rise
    /// with porosity, so a bed that is packed looser flows faster either way and no ordinary test
    /// distinguishes them — the difference only shows in *how much*, which is what
    /// `the_porosity_dependence_is_kozeny_carman` measures over a range wide enough to tell.
    pub fn permeability(&self, porosity: f64) -> Area {
        let e = porosity.clamp(0.0, 0.999);
        Area::from_si(self.hydraulic * self.hydraulic * e * e * e / (180.0 * (1.0 - e) * (1.0 - e)))
    }

    /// The first-order extraction rate constant at a temperature, s⁻¹.
    ///
    /// Diffusion out of a sphere of radius `a` decays as a series whose slowest term is
    /// `exp(−π²Dt/a²)`. After the first fraction of a second the faster terms are gone and what is
    /// left is a single exponential, so the whole of the particle-scale physics collapses to one
    /// rate:
    ///
    /// ```text
    ///   K = π² D(T) / a²,   a = d/2   ⇒   K = 4π² D(T) / d²
    /// ```
    ///
    /// This is the *sieve* diameter, not the hydraulic one — see the type's documentation for why
    /// those are different numbers.
    pub fn extraction_rate(&self, bed: &Bed, t: Temperature) -> f64 {
        let a = self.sieve / 2.0;
        std::f64::consts::PI.powi(2) * bed.diffusivity_at(t).to_si() / (a * a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The viscosity fit extrapolates to a point it was not given.**
    ///
    /// Two tabulated values set `A` and `B`; a third, at 100 °C, checks them. A fit that only
    /// reproduces its own two points has been assumed rather than tested.
    #[test]
    fn the_viscosity_fit_predicts_a_point_it_was_not_fitted_to() {
        let w = Liquid::water();
        let at_100 = w.viscosity(Temperature::celsius(100.0)).to_si();
        let tabulated = 2.82e-4;
        let off = (at_100 - tabulated).abs() / tabulated;
        assert!(
            off < 0.01,
            "the Andrade fit should hold to a percent outside its two points: {at_100:.4e} \
             against {tabulated:.4e}, off by {:.2}%",
            off * 100.0
        );
        // And it is doing real work: the range it spans is not flat.
        let at_70 = w.viscosity(Temperature::celsius(70.0)).to_si();
        assert!(
            at_70 / at_100 > 1.3,
            "water is a third more viscous at 70 C than at 100; ratio was {:.3}",
            at_70 / at_100
        );
    }

    /// **Permeability goes as the square of the grind, and the extraction rate as its inverse.**
    ///
    /// The two are derived from separate physics — Kozeny–Carman for one, sphere diffusion for the
    /// other — and this is the check that they carry the exponent each of them should, in the
    /// direction each should. A sign error either way would still produce a shot.
    #[test]
    fn halving_the_grind_quarters_the_permeability_and_quadruples_the_rate() {
        let bed = Bed::roasted_coffee();
        let coarse = Grind::espresso();
        let fine = coarse.scaled(0.5);
        let porosity = 0.30;

        let kr = coarse.permeability(porosity).to_si() / fine.permeability(porosity).to_si();
        let er = fine.extraction_rate(&bed, Temperature::celsius(93.0))
            / coarse.extraction_rate(&bed, Temperature::celsius(93.0));
        assert!(
            (kr - 4.0).abs() < 1e-12,
            "k goes as d^2: ratio was {kr:.6} for a halving"
        );
        assert!(
            (er - 4.0).abs() < 1e-12,
            "K goes as 1/d^2: ratio was {er:.6} for a halving"
        );
    }

    /// **The porosity dependence is `ε³/(1−ε)²` and not `ε³/(1−ε)`.**
    ///
    /// Both are increasing in `ε`, so no test that only checks the direction can tell them apart.
    /// Measured over a range wide enough that they part company: from 0.25 to 0.45 the correct
    /// form gives 10.84× and the misremembered one 7.95×, which is 36% apart. Anything narrower
    /// than that range and the test passes for either.
    #[test]
    fn the_porosity_dependence_is_kozeny_carman() {
        let g = Grind::espresso();
        let ratio = g.permeability(0.45).to_si() / g.permeability(0.25).to_si();
        let correct = (0.45f64.powi(3) / 0.55f64.powi(2)) / (0.25f64.powi(3) / 0.75f64.powi(2));
        let misremembered = (0.45f64.powi(3) / 0.55) / (0.25f64.powi(3) / 0.75);
        assert!(
            (ratio - correct).abs() / correct < 1e-12,
            "measured {ratio:.4}, the law gives {correct:.4}"
        );
        assert!(
            (correct / misremembered - 1.0).abs() > 0.3,
            "this test only means something if the two forms differ here: {correct:.3} against \
             {misremembered:.3}"
        );
    }

    /// **The Arrhenius factor is the right size for the activation energy quoted.**
    ///
    /// Checked against the closed form rather than against itself, at a step a barista actually
    /// makes: five degrees. 30 kJ/mol over 5 K near 93 °C is a 13% change in rate, which is why
    /// brew temperature is a knob at all and also why it is a gentle one.
    #[test]
    fn five_degrees_moves_the_extraction_rate_thirteen_percent() {
        let bed = Bed::roasted_coffee();
        let g = Grind::espresso();
        let hot = g.extraction_rate(&bed, Temperature::celsius(93.0));
        let cool = g.extraction_rate(&bed, Temperature::celsius(88.0));
        let predicted =
            (-30_000.0 / R_GAS * (1.0 / Temperature::celsius(88.0).to_si() - 1.0 / 366.15)).exp();
        assert!(
            (cool / hot - predicted).abs() < 1e-12,
            "measured {:.6}, Arrhenius gives {predicted:.6}",
            cool / hot
        );
        assert!(
            (0.85..0.89).contains(&(cool / hot)),
            "five degrees should be about 13%: was {:.1}%",
            (1.0 - cool / hot) * 100.0
        );
    }

    /// **The calibration lands where it was aimed, in the two places it was aimed at.**
    ///
    /// `FINES_RATIO` and `solute_diffusivity` were each set from one measurement, and this names
    /// them — so that a change to either shows up as a change to a shot rather than as a number
    /// moving in a table.
    ///
    /// The flow is checkable here in closed form. The yield is not, because it now depends on how
    /// the liquid loads up as it crosses the bed, and that is a transport problem rather than an
    /// arithmetic one — it is checked end to end in `a_conventional_shot_lands_where_a_barista_aims`.
    #[test]
    fn the_calibration_reproduces_a_conventional_shot() {
        let bed = Bed::roasted_coffee();
        let g = Grind::espresso();
        let (porosity, depth, radius, dp) = (0.45, 0.020, 0.029, 9.0e5);
        let area = std::f64::consts::PI * radius * radius;
        let mu = Liquid::water()
            .viscosity(Temperature::celsius(93.0))
            .to_si();

        // The dose is a consequence of the packing, and it is the number a person recognises.
        let dose = (1.0 - porosity) * bed.solid_density.to_si() * area * depth * 1000.0;
        assert!(
            (16.0..20.0).contains(&dose),
            "a 58 mm basket 20 mm deep holds about 18 g; this packing gives {dose:.2} g"
        );

        let q = g.permeability(porosity).to_si() * area * dp / (mu * depth);
        let g_per_s = q * Liquid::water().density.to_si() * 1000.0;
        assert!(
            (1.3..1.8).contains(&g_per_s),
            "a conventional shot runs at about 1.5 g/s; this one at {g_per_s:.3}"
        );

        // And the diffusivity is a plausible fraction of a bulk one rather than whatever made the
        // fit work. Sucrose in water at 93 C is near 1.9e-9 m2/s.
        let tortuosity = 1.9e-9 / bed.solute_diffusivity.to_si();
        assert!(
            (3.0..20.0).contains(&tortuosity),
            "an effective diffusivity should be a few-fold below bulk, not orders: {tortuosity:.1}x"
        );
    }
}
