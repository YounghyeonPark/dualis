//! Two substances made into one, and the honest answer about what that one is.
//!
//! A motor is copper, electrical steel, magnets and air. A populated board is FR-4, copper and
//! solder. A printed part is PLA and voids. A thermal buffer is a wax in an aluminium matrix. Every
//! one of those wants to be **one** [`Substance`] so a lumped model or a coarse grid can hold it, and
//! the properties of that one substance are not the properties of its main constituent —
//! [`Substance::with_specific_heat`] says so and is worth a factor of two on a motor.
//!
//! What this module refuses to do is make the numbers up. The properties of a mixture divide into
//! three kinds, and conflating them is the whole failure mode:
//!
//! | | what a mixture rule can say |
//! | --- | --- |
//! | density, volumetric heat capacity, latent heat | **exact**, from conservation alone |
//! | conductivity, stiffness | **bounded**. No single value exists without knowing the microstructure |
//! | emissivity | **nothing**. It is a property of the surface, and a mixture has no surface |
//!
//! The middle row is the one that matters, and it is not a small effect. A 50/50 aluminium and
//! borosilicate composite conducts somewhere between 2.21 and 84.06 W/m·K — a factor of **38** — and
//! which end depends entirely on whether the glass is in plates across the flux or in fibres along it.
//! Both extremes are realisable. A library that answered `0.5·167 + 0.5·1.114` would be handing back
//! the upper bound of a 38-fold range as if it were a measurement. On a half-copper, half-FR-4 board it
//! is a **335-fold** range.
//!
//! So [`Mix`] reports the exact properties as values and the bounded ones as **bounds**, and
//! [`Mix::as_substance`] makes the caller supply the conductivity and refuses one outside the bounds.
//! Choosing is the caller's job; the library's job is to refuse an impossible choice.
//!
//! # The bounds, and why there are two pairs
//!
//! **Voigt and Reuss** are the arithmetic and harmonic means by volume fraction. They hold for any
//! microstructure whatever and they are *attained* — a laminate loaded along its layers is exactly
//! Voigt, and the same laminate across its layers is exactly Reuss. `a_composite.rs` measures both on
//! one block, which is also the demonstration that no single number can exist: the same composite has
//! two different conductivities depending on the direction of the flux.
//!
//! **Hashin–Shtrikman** are tighter and buy that with an assumption: the microstructure is
//! statistically *isotropic*. A laminate is not, which is why HS does not contain the laminate's own
//! answers and is not a replacement for the outer pair. For a foam, a filled polymer or a packed
//! powder it is the pair to use, and it is narrower — for the aluminium and glass above, 4.33 to 67.6
//! against 2.21 to 84.06, which takes the range from 38-fold to 16-fold. Narrower and still wide: with
//! a 150-fold contrast between the phases no bound is going to be comfortable, and that is the honest
//! state of the problem rather than a deficiency of the bound.
//!
//! # What is not here, and why not
//!
//! No elastic or acoustic bounds. Voigt–Reuss on the bulk and shear moduli is the same theorem and
//! the algebra is no harder — but `dualis-elastic` takes one substance per body and has no per-cell
//! material, so there is nothing in this workspace a bound on stiffness could be *checked* against.
//! A bound nothing can falsify is a comment, and this workspace does not ship those as API. It becomes
//! available the moment `Elastic` grows a `fill`.
//!
//! No yield strength either, and that one is not a missing feature. A composite's yield is governed by
//! the weaker phase and by the interface between them, so it is not a mixture of the two yields in any
//! ordering — Voigt–Reuss does not bound it and a rule of mixtures for it would be wrong rather than
//! imprecise.

use crate::substance::{FusionProps, Substance, ThermalProps};
use dualis_units::{
    Density, HeatCapacity, LatentHeat, Mass, SpecificHeat, Temperature, ThermalConductivity,
    ThermalExpansion, Volume,
};

/// A composite: substances and the fraction of the **volume** each occupies.
///
/// Volume and not mass, because volume is what the geometry gives and what every bound is written in.
/// The mass fractions are derived — see [`Mix::mass_fraction`] — and getting the two the wrong way
/// round is the classic mistake this type exists to prevent.
#[derive(Clone, Debug, PartialEq)]
pub struct Mix {
    parts: Vec<(Substance, f64)>,
}

impl Mix {
    /// A mixture from substances and their volume fractions.
    ///
    /// The fractions must be positive and sum to one to within `1e-9`. Refused rather than normalised,
    /// and that is deliberate: fractions that do not sum to one are a transcription mistake, and
    /// normalising them silently would turn 45% and 50% into 47.4% and 52.6% and answer a question
    /// nobody asked. `1e-9` and not exact equality because `0.3 + 0.3 + 0.4` is not `1.0` in binary.
    ///
    /// A single part is allowed and is not a mistake — a mixture of one is how a caller writes "this
    /// is not a composite after all" without changing the shape of the code around it, and every bound
    /// below collapses to the substance's own value.
    pub fn of(parts: &[(Substance, f64)]) -> Result<Mix, String> {
        if parts.is_empty() {
            return Err("a mixture needs at least one substance".to_string());
        }
        let mut total = 0.0;
        for (s, f) in parts {
            if !(f.is_finite() && *f > 0.0) {
                return Err(format!(
                    "{}: a volume fraction must be finite and positive, is {f}",
                    s.name
                ));
            }
            total += f;
        }
        if (total - 1.0).abs() > 1e-9 {
            return Err(format!(
                "volume fractions must sum to 1, sum to {total} — they are not normalised for you, \
                 because 45% and 50% is a transcription mistake and not a request for 47.4% and 52.6%"
            ));
        }
        Ok(Mix {
            parts: parts.to_vec(),
        })
    }

    /// The substances and their volume fractions, in the order given.
    pub fn parts(&self) -> &[(Substance, f64)] {
        &self.parts
    }

    /// **Exact.** `ρ = Σ φᵢ ρᵢ`, which is mass conservation and not a model.
    pub fn density(&self) -> Density {
        Density::from_si(
            self.parts
                .iter()
                .map(|(s, f)| f * s.density.to_si())
                .sum::<f64>(),
        )
    }

    /// What fraction of the **mass** the `i`th part is, or `None` if there is no such part.
    ///
    /// `wᵢ = φᵢ ρᵢ / ρ`. The conversion that makes the difference between a correct specific heat and
    /// one that is wrong by the density ratio.
    pub fn mass_fraction(&self, i: usize) -> Option<f64> {
        let (s, f) = self.parts.get(i)?;
        Some(f * s.density.to_si() / self.density().to_si())
    }

    /// **Exact**, and the field where the volume-and-mass confusion costs the most.
    ///
    /// Volumetric heat capacity `ρc` is volume-additive — `ρc = Σ φᵢ ρᵢ cᵢ` — because a joule stored in
    /// a cubic metre of composite is the joules stored in each part's share of that cubic metre. Divide
    /// by the mixture's own density and the per-kilogram figure is **mass**-weighted:
    /// `c = Σ wᵢ cᵢ`.
    ///
    /// Volume-weighting `c_p` directly is the mistake, and how large it is depends entirely on which
    /// pair — which is what makes it dangerous, because it is invisible on the example somebody checks
    /// with. Measured, at half and half by volume:
    ///
    /// ```text
    ///   aluminium + borosilicate    877.7 correct    877.0 volume-weighted     0.08% out
    ///   aluminium + copper          503.3            640.5                    27.25%
    ///   copper + FR-4               507.4            742.5                    46.34%
    /// ```
    ///
    /// The first pair have nearly the same density, so the two rules agree to a tenth of a percent and a
    /// caller who tried it there would conclude the distinction does not matter. On a copper and FR-4
    /// board — the case this module's opening paragraph names — it is 46%. The error always runs toward
    /// the **lighter** constituent, because volume weighting over-counts a light phase whose `c_p` per
    /// kilogram is high.
    ///
    /// `None` if any part does not state a specific heat, because a mixture containing an unknown is
    /// unknown and not the average of what happens to be known.
    pub fn specific_heat(&self) -> Option<SpecificHeat> {
        let mut volumetric = 0.0;
        for (s, f) in &self.parts {
            let t = s.thermal?;
            volumetric += f * s.density.to_si() * t.specific_heat.to_si();
        }
        Some(SpecificHeat::from_si(volumetric / self.density().to_si()))
    }

    /// What a volume of the mixture holds per kelvin. Exact, for the reason above.
    pub fn heat_capacity(&self, volume: Volume) -> Option<HeatCapacity> {
        Some(HeatCapacity::from_si(
            volume.to_si() * self.density().to_si() * self.specific_heat()?.to_si(),
        ))
    }

    /// The mass of a volume of the mixture.
    pub fn mass_of(&self, volume: Volume) -> Mass {
        Mass::from_si(volume.to_si() * self.density().to_si())
    }

    /// **Voigt and Reuss**, in that order: the arithmetic and harmonic means by volume fraction.
    ///
    /// `(k_low, k_high)` — Reuss first, because a returned pair should read low to high.
    ///
    /// These hold for **any** microstructure and both are attained, so they are the widest correct
    /// answer and also the tightest one that assumes nothing. A laminate is the witness for both:
    /// flux along the layers gives Voigt exactly and across them gives Reuss exactly, which
    /// `a_composite.rs` measures on one block to machine precision.
    ///
    /// `None` if any part's conductivity is unknown.
    pub fn conductivity_bounds(&self) -> Option<(ThermalConductivity, ThermalConductivity)> {
        let mut voigt = 0.0;
        let mut reciprocal = 0.0;
        for (s, f) in &self.parts {
            let k = s.thermal?.conductivity.to_si();
            voigt += f * k;
            reciprocal += f / k;
        }
        Some((
            ThermalConductivity::from_si(1.0 / reciprocal),
            ThermalConductivity::from_si(voigt),
        ))
    }

    /// **Hashin–Shtrikman**, for a mixture of exactly two substances whose microstructure is
    /// statistically isotropic — a foam, a filled polymer, a packed powder.
    ///
    /// `(k_low, k_high)`, and both lie inside [`Mix::conductivity_bounds`]. The price of the tighter
    /// pair is the isotropy assumption, so this is **not** a strictly better answer: a laminate is
    /// anisotropic and its own exact conductivities fall outside these. Reaching for HS on a layered
    /// material is the one way to be wrong with it.
    ///
    /// ```text
    ///   k± = k_a + φ_b / ( 1/(k_b − k_a) + φ_a/(3 k_a) )
    /// ```
    ///
    /// with `a` the more conductive phase for the upper bound and the less for the lower. The bounds
    /// are attained by coated-sphere assemblages, which is why they are bounds and not a fit.
    ///
    /// `None` unless there are exactly two parts with known conductivities. Two parts of *equal*
    /// conductivity give a degenerate pair, which is correct — both bounds are that conductivity.
    pub fn hashin_shtrikman(&self) -> Option<(ThermalConductivity, ThermalConductivity)> {
        if self.parts.len() != 2 {
            return None;
        }
        let (k0, f0) = (
            self.parts[0].0.thermal?.conductivity.to_si(),
            self.parts[0].1,
        );
        let (k1, f1) = (
            self.parts[1].0.thermal?.conductivity.to_si(),
            self.parts[1].1,
        );
        if k0 == k1 {
            return Some((
                ThermalConductivity::from_si(k0),
                ThermalConductivity::from_si(k0),
            ));
        }
        // `host` is the phase the assemblage is built around; taking the stiffer one gives the upper
        // bound and the softer one the lower.
        let bound = |host: f64, host_f: f64, guest: f64, guest_f: f64| {
            host + guest_f / (1.0 / (guest - host) + host_f / (3.0 * host))
        };
        let (hi_host, hi_f, hi_guest, hi_gf) = if k0 > k1 {
            (k0, f0, k1, f1)
        } else {
            (k1, f1, k0, f0)
        };
        Some((
            ThermalConductivity::from_si(bound(hi_guest, hi_gf, hi_host, hi_f)),
            ThermalConductivity::from_si(bound(hi_host, hi_f, hi_guest, hi_gf)),
        ))
    }

    /// **Exact.** The latent heat per kilogram *of the mixture*, when exactly one part melts.
    ///
    /// `L = w L_part` at that part's melting point, by mass conservation. This is the phase-change
    /// composite — a wax in an aluminium matrix, a salt hydrate in a foam — and it is the case where a
    /// mixture rule is not an approximation at all: the joules are the joules, and diluting the PCM
    /// dilutes them in exact proportion to its mass fraction.
    ///
    /// `w` and not the volume fraction, and here the difference is the whole answer rather than a
    /// refinement. Wax at 814 kg/m³ filling **80% of the volume** of an aluminium matrix is only
    /// **54.7% of the mass**, so a composite whose wax stores 244 kJ/kg stores 133.4 — and using the
    /// volume fraction would claim 195, **46% high**, on the one property a buffer is bought for.
    ///
    /// `None` if nothing melts. **`Err` if two or more parts melt**, because two melting points is not
    /// one melting point and a composite with two plateaux cannot be described by a type that has one
    /// `melting_point` field. Refused rather than answered with the larger, the first, or an average.
    #[allow(clippy::type_complexity)]
    pub fn fusion(&self) -> Result<Option<(Temperature, LatentHeat)>, String> {
        let melting: Vec<usize> = (0..self.parts.len())
            .filter(|i| self.parts[*i].0.fusion.is_some())
            .collect();
        match melting.as_slice() {
            [] => Ok(None),
            [i] => {
                let f = self.parts[*i].0.fusion.expect("filtered on it");
                let w = self.mass_fraction(*i).expect("index came from the list");
                Ok(Some((
                    f.melting_point,
                    LatentHeat::from_si(w * f.latent_heat.to_si()),
                )))
            }
            many => Err(format!(
                "{} of the parts melt — {} — and a composite with two plateaux is not describable by \
                 one melting point. Mix the non-melting parts and handle the phase change as its own \
                 region, or model the second one as inert",
                many.len(),
                many.iter()
                    .map(|i| format!("{:?}", self.parts[*i].0.name))
                    .collect::<Vec<_>>()
                    .join(" and ")
            )),
        }
    }

    /// The mixture as one [`Substance`], with the conductivity and emissivity the caller chose.
    ///
    /// # Why those two are arguments
    ///
    /// **The conductivity, because no single value exists.** The bounds are as far as physics goes
    /// without the microstructure, and picking a point inside them is a modelling decision with a
    /// reason attached — the midpoint of the Hashin–Shtrikman pair for a foam, the Voigt bound for
    /// unidirectional fibres along the flux, a measured value if there is one. A value outside
    /// [`Mix::conductivity_bounds`] is **refused**: it is not conservative or approximate, it is
    /// impossible, and no microstructure realises it.
    ///
    /// **The emissivity, because a mixture has no surface.** It is not a bulk property and it does not
    /// mix. A half-copper, half-FR-4 board is 0.05 if it is bare metal on the outside and 0.9 if it is
    /// green solder mask, and the volume fractions say nothing about which.
    ///
    /// Expansion is volume-weighted, and that one *is* an approximation rather than a bound — the Voigt
    /// rule for CTE, which overestimates when the stiff phase is the one that expands less, because the
    /// stiff phase restrains the compliant one. Turner's and Kerner's rules weight by stiffness and are
    /// what a thermal-stress calculation wants; this is the number for a rough length change, and it
    /// is documented here rather than silently returned as if it were the exact ones above.
    ///
    /// `Err` if any part's thermal properties are unknown, if the conductivity is outside the bounds,
    /// if the emissivity is not a fraction, or if more than one part melts.
    pub fn as_substance(
        &self,
        name: &str,
        conductivity: ThermalConductivity,
        emissivity: f64,
    ) -> Result<Substance, String> {
        let (low, high) = self
            .conductivity_bounds()
            .ok_or_else(|| format!("{name}: a part does not state its conductivity"))?;
        let k = conductivity.to_si();
        if !(k.is_finite() && k >= low.to_si() * (1.0 - 1e-12) && k <= high.to_si() * (1.0 + 1e-12))
        {
            return Err(format!(
                "{name}: no microstructure of these parts conducts {k} W/m/K — the Voigt and Reuss \
                 bounds are {} to {}, and they are attained, so this is impossible rather than \
                 merely unlikely",
                low.to_si(),
                high.to_si()
            ));
        }
        if !(0.0..=1.0).contains(&emissivity) {
            return Err(format!(
                "{name}: emissivity is a fraction of a blackbody's and must be in 0..=1, is \
                 {emissivity}"
            ));
        }
        let specific_heat = self
            .specific_heat()
            .ok_or_else(|| format!("{name}: a part does not state its specific heat"))?;
        let expansion = ThermalExpansion::from_si(
            self.parts
                .iter()
                .map(|(s, f)| f * s.thermal.map_or(f64::NAN, |t| t.expansion.to_si()))
                .sum::<f64>(),
        );
        let mut out = Substance::bulk(name, self.density()).with_thermal(ThermalProps {
            conductivity,
            specific_heat,
            expansion,
            emissivity,
        });
        if let Some((point, latent)) = self.fusion().map_err(|e| format!("{name}: {e}"))? {
            out = out.with_fusion(FusionProps::new(point, latent));
        }
        Ok(out)
    }
}
