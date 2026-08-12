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
//! There is an HS pair for stiffness too — [`Mix::shear_hashin_shtrikman`] and
//! [`Mix::bulk_hashin_shtrikman`] — and how *tight* it is turns out to depend on which modulus you ask
//! about. For a three-dimensional checkerboard, measured, the upper bound on the shear modulus is tight to
//! within 0.5% and the one on the bulk modulus is at least 2.8% loose. That is not something the algebra
//! says, and `a_checkerboard.rs` is where it is measured.
//!
//! # Stiffness, and the reason there is no effective `(E, ν)`
//!
//! [`Mix::shear_bounds`] and [`Mix::p_wave_modulus_bounds`] arrived once `Waves::fill` existed, because
//! until there was per-element material in `dualis-elastic` there was nothing in this workspace a bound
//! on stiffness could be *checked* against, and a bound nothing can falsify is a comment rather than an
//! API. Both are now checked against **Backus averaging** — the exact long-wavelength moduli of a layered
//! elastic medium — in `crates/dualis-elastic/tests/a_layered_wave.rs`, and the harmonic end again
//! statically in `a_layered_block.rs`, which gets it nine orders sharper because an elliptic solve has no
//! time in it: `4.8e-13` against the wave's `3.5e-4`.
//!
//! What is deliberately absent is a Young's modulus and a Poisson ratio for the mixture, and that
//! absence is the physics rather than a gap. **A composite of two isotropic materials is generally
//! anisotropic.** A laminate has a different stiffness along its layers than across them — measured, the
//! shear modulus differs by a factor of **5.5** for aluminium against PLA — so there is no single pair
//! `(E, ν)` that describes it, and a function returning one would be inventing an isotropy the material
//! does not have. `Mix` therefore does not produce an [`crate::substance::MechanicalProps`] and
//! [`Mix::as_substance`] leaves the mechanical block absent.
//!
//! No yield strength either, and that one is not a missing feature. A composite's yield is governed by
//! the weaker phase and by the interface between them, so it is not a mixture of the two yields in any
//! ordering — Voigt–Reuss does not bound it and a rule of mixtures for it would be wrong rather than
//! imprecise.

use crate::substance::{FusionProps, Substance, ThermalProps};
use dualis_units::{
    Density, HeatCapacity, LatentHeat, Mass, Pressure, SpecificHeat, Temperature,
    ThermalConductivity, ThermalExpansion, Volume,
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

    /// **Voigt and Reuss on the shear modulus**, low first: `(⟨1/G⟩⁻¹, ⟨G⟩)`.
    ///
    /// `G = E/(2(1+ν))` for each part, weighted by volume fraction. Both ends are **attained**, and by
    /// the same witness in two directions — which is what makes this a range of achievable values rather
    /// than a hedge:
    ///
    /// - a laminate sheared **in** its layer planes carries a uniform shear strain, so the stresses add
    ///   and the effective modulus is `⟨G⟩` exactly. That is `C66` in the layered-medium literature;
    /// - the same laminate sheared **across** its layers carries a uniform shear stress, so the strains
    ///   add and it is `⟨1/G⟩⁻¹` exactly. That is `C44`.
    ///
    /// Both are Backus's 1962 results for a finely layered elastic medium, and both are measured against
    /// marched wave speeds in `a_layered_wave.rs` — **5.5× apart** for aluminium against PLA, from one
    /// block, with each end converging at second order to 0.06% or better.
    ///
    /// `None` if any part does not state its mechanical properties.
    pub fn shear_bounds(&self) -> Option<(Pressure, Pressure)> {
        self.moduli_bounds(|e, nu| e / (2.0 * (1.0 + nu)))
    }

    /// **Voigt and Reuss on the bulk modulus**, low first: `(⟨1/K⟩⁻¹, ⟨K⟩)`.
    ///
    /// `K = E/(3(1−2ν))`. Neither end is attained by a laminate, and no witness in this workspace attains
    /// either — a laminate under hydrostatic stress has lateral constraint between its layers, so it
    /// reaches neither the uniform-stress nor the uniform-strain state. They are correct bounds all the
    /// same, and [`Mix::bulk_hashin_shtrikman`] is the pair to prefer for anything isotropic.
    ///
    /// Stated because it is the difference between this pair and [`Mix::shear_bounds`], whose ends both
    /// *are* attained and measured. A bound that is reachable and a bound that merely holds are different
    /// things to a caller choosing a number inside them.
    pub fn bulk_bounds(&self) -> Option<(Pressure, Pressure)> {
        self.moduli_bounds(|e, nu| e / (3.0 * (1.0 - 2.0 * nu)))
    }

    /// **Hashin–Shtrikman on the shear modulus** for exactly two phases, low first.
    ///
    /// ```text
    ///   G± = G_r + f_o / [ 1/(G_o − G_r) + 6 f_r (K_r + 2G_r) / (5 G_r (3K_r + 4G_r)) ]
    /// ```
    ///
    /// with `r` the reference phase — the stiffer one for the upper bound, the softer for the lower — and
    /// `o` the other. Hashin and Shtrikman 1963, and the same trade the conductivity pair makes: tighter
    /// than Voigt–Reuss, at the price of assuming the microstructure is statistically **isotropic**. For
    /// aluminium against PLA at half and half it takes the range from 5.5-fold to 2.8-fold.
    ///
    /// # The check that says this is the right algebra
    ///
    /// Taken with the **matrix** as reference it is identically the **Mori–Tanaka** estimate for spherical
    /// inclusions — a separately derived result written as a different rational function — and
    /// `a_mixture.rs` measures the two agreeing to `2.2e-16` across two decades of inclusion fraction.
    /// That equivalence is a theorem rather than a coincidence: the bound is attained by a coated-sphere
    /// assemblage, which is what Mori–Tanaka describes.
    ///
    /// # What no witness here attains
    ///
    /// Unlike the conductivity pair, the elastic HS bounds are **not** bracketed from below by a
    /// measurement in this workspace, and the reason is worth knowing. A resolved isotropic geometry has
    /// to be driven by something, and an affine displacement on the boundary is a *kinematically
    /// admissible* field — so its energy is an **upper** estimate of the effective modulus, above the true
    /// value however fine the mesh. `a_checkerboard.rs` measures that estimate converging down to
    /// **1.005×** the upper bound for a well-resolved board and states plainly that it cannot cross it.
    /// That is evidence the bound is nearly tight and it is not the same as bracketing.
    ///
    /// `None` unless there are exactly two parts with mechanical properties.
    pub fn shear_hashin_shtrikman(&self) -> Option<(Pressure, Pressure)> {
        self.hashin_shtrikman_pair(false)
    }

    /// **Hashin–Shtrikman on the bulk modulus** for exactly two phases, low first.
    ///
    /// ```text
    ///   K± = K_r + f_o / [ 1/(K_o − K_r) + 3 f_r / (3K_r + 4G_r) ]
    /// ```
    ///
    /// Everything [`Mix::shear_hashin_shtrikman`] says applies, including the Mori–Tanaka equivalence,
    /// which for `K` is the more familiar of the two. Note the shear modulus of the *reference* phase
    /// appears in a bound on `K`: a stiff inclusion resists the hydrostatic compression of its
    /// surroundings partly in shear, so the two moduli do not separate.
    pub fn bulk_hashin_shtrikman(&self) -> Option<(Pressure, Pressure)> {
        self.hashin_shtrikman_pair(true)
    }

    /// Both elastic Hashin–Shtrikman pairs, since they differ only in one term.
    ///
    /// # The reference phase is not chosen, it is tried both ways
    ///
    /// The textbook prescription is "put the stiffest phase in the reference position for the upper bound
    /// and the softest for the lower", and that instruction only means something for a **well-ordered**
    /// pair — one phase larger in both `K` and `G`. Aluminium against borosilicate is not: aluminium has
    /// the larger bulk modulus, 67.5 GPa against 46.5, and the *smaller* shear modulus, 25.9 against 34.0.
    ///
    /// A first version tested which phase had the larger value of the modulus being bounded and used that
    /// one as the upper reference. For that pair at a tenth aluminium it returned a lower bound of 48.2312
    /// GPa above an upper bound of 48.1922 — **the pair inverted**, by 0.08%, which is small enough that
    /// only a test sweeping fractions and pairs would see it.
    ///
    /// So both evaluations are computed and then **ordered**. That is what "interchange which phase is
    /// subscripted one" actually prescribes, it needs no notion of stiffer, and it is right for a
    /// well-ordered pair and for the other kind alike.
    fn hashin_shtrikman_pair(&self, bulk: bool) -> Option<(Pressure, Pressure)> {
        if self.parts.len() != 2 {
            return None;
        }
        let of = |i: usize| -> Option<(f64, f64, f64)> {
            let m = self.parts[i].0.mechanical?;
            let (e, nu) = (m.youngs_modulus.to_si(), m.poisson_ratio);
            Some((
                e / (3.0 * (1.0 - 2.0 * nu)),
                e / (2.0 * (1.0 + nu)),
                self.parts[i].1,
            ))
        };
        let (a, b) = (of(0)?, of(1)?);
        let bound = |r: (f64, f64, f64), o: (f64, f64, f64)| {
            let (kr, gr, fr) = r;
            let (ko, go, fo) = o;
            if bulk {
                kr + fo / (1.0 / (ko - kr) + 3.0 * fr / (3.0 * kr + 4.0 * gr))
            } else {
                gr + fo
                    / (1.0 / (go - gr)
                        + 6.0 * fr * (kr + 2.0 * gr) / (5.0 * gr * (3.0 * kr + 4.0 * gr)))
            }
        };
        // Equal moduli make the pair degenerate, which is correct — a mixture of two things with the same
        // `G` has that `G` — and the expression divides by their difference, so it has to come first.
        let (ma, mb) = if bulk { (a.0, b.0) } else { (a.1, b.1) };
        if ma == mb {
            return Some((Pressure::from_si(ma), Pressure::from_si(ma)));
        }
        let (one, other) = (bound(a, b), bound(b, a));
        Some((
            Pressure::from_si(one.min(other)),
            Pressure::from_si(one.max(other)),
        ))
    }

    /// **Voigt and Reuss on the P-wave modulus** `M = λ + 2μ`, low first: `(⟨1/M⟩⁻¹, ⟨M⟩)`.
    ///
    /// `M = E(1−ν)/((1+ν)(1−2ν))`, the modulus that relates stress to strain when the lateral strain is
    /// held at zero — so it is the one a compression wave travels on and the one a thin layer bonded
    /// between stiff neighbours actually feels.
    ///
    /// **Both ends are attained**, and a first draft of this documentation said the high end was not —
    /// the measurement is what corrected it.
    ///
    /// - a laminate compressed **across** its layers carries a uniform normal stress, so the compliances
    ///   add and `⟨1/M⟩⁻¹` is exact. Backus's `C33`, and the speed of a compression wave through the
    ///   stack;
    /// - the same laminate compressed **along** its layers, with the lateral strain held at zero
    ///   *pointwise*, carries a uniform strain, so the stresses add and `⟨M⟩` is exact.
    ///
    /// The second needs that constraint said out loud, because it is what the Voigt bound *is*. A laminate
    /// whose lateral contraction is free gives neither bound: the layers each want to contract differently
    /// and the ones beside them prevent it, and the answer is Backus's `C11`, which carries a correction
    /// term in `⟨λ/M⟩` and for aluminium against PLA is 43.77 GPa against `⟨M⟩`'s 53.98 — **18.9% below**.
    /// `a_layered_wave.rs` says why it does not measure that one: it needs the lateral strain zero on
    /// average but free locally, and `Waves::hold` holds a component everywhere or nowhere.
    ///
    /// `None` if any part does not state its mechanical properties.
    pub fn p_wave_modulus_bounds(&self) -> Option<(Pressure, Pressure)> {
        self.moduli_bounds(|e, nu| e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)))
    }

    /// The Voigt and Reuss pair for any modulus derivable from `(E, ν)`, low first.
    ///
    /// One helper because the two public pairs differ only in which modulus, and writing the weighting
    /// twice is how two bounds come to disagree about what a volume fraction is.
    fn moduli_bounds(&self, modulus: impl Fn(f64, f64) -> f64) -> Option<(Pressure, Pressure)> {
        let mut voigt = 0.0;
        let mut reciprocal = 0.0;
        for (s, f) in &self.parts {
            let m = s.mechanical?;
            let value = modulus(m.youngs_modulus.to_si(), m.poisson_ratio);
            if !(value.is_finite() && value > 0.0) {
                return None;
            }
            voigt += f * value;
            reciprocal += f / value;
        }
        Some((
            Pressure::from_si(1.0 / reciprocal),
            Pressure::from_si(voigt),
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
