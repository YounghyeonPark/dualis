//! A lamp shining on a coated surface, as a domain.
//!
//! The fifth field, and the one the library cannot supply itself: `pantometry-optics` has no
//! [`Domain`] in it at all. It has spectra, Fresnel coefficients, coatings and dispersion —
//! all of which answer a question rather than march a state — so nothing in it is a thing
//! that steps. The workspace's own coupling tests define an absorbing surface locally for the
//! same reason.
//!
//! That makes this the fourth `Domain` written outside the library, and the one that most
//! looks like what a consumer would actually write: real optics on the front, a heat channel
//! on the back, and no knowledge whatever of who takes the joules.
//!
//! # What the optics actually decides
//!
//! A tungsten lamp is not white and a coating is not flat, so how much of a hundred watts
//! ends up as heat depends on both curves and on where they overlap. Here:
//!
//! - the source is a blackbody at its colour temperature, integrated over the visible range;
//! - the surface is a mirror whose reflectance *varies with wavelength* — aluminium-like,
//!   worse in the blue — and whatever it does not reflect it absorbs, which is what makes a
//!   mirror a heat source;
//! - the absorptance is sampled every 5 nm and handed to the radiometry as a spectrum, and
//!   the integral of the product is the absorbed power.
//!
//! Change the colour temperature and the absorbed fraction changes, because a hotter lamp
//! puts more of its output into the blue where the mirror is worse. A flat reflectance would
//! make that dependence vanish, and the whole spectral apparatus with it — which is why the
//! mirror here is not flat and why `tests/scene.rs` checks the difference between 2800 K and
//! 6500 K rather than checking a single number.

use pantometry::core::Reading;
use pantometry::optics::spectrum::Spectrum as Spec;
use pantometry::prelude::*;

/// An aluminium-like mirror: worse in the blue, which is what makes the lamp's colour matter.
///
/// Real numbers, near enough — evaporated aluminium runs about 92% at 400 nm and 97% by
/// 700 nm, with the dip near 800 nm that gives it its faint yellow cast.
pub fn aluminium_mirror() -> Spectrum {
    Spectrum::curve(vec![
        (380.0, 0.905),
        (450.0, 0.918),
        (500.0, 0.925),
        (550.0, 0.935),
        (600.0, 0.947),
        (650.0, 0.960),
        (700.0, 0.968),
    ])
}

/// A lamp, a coating, and the heat that comes of the two meeting.
pub struct Light {
    name: String,
    lamp: SpectralPower,
    absorptance: Spectrum,
    /// Joules not yet spent, so the books close. See [`crate::heater::Heater`].
    reserve: f64,
    saved: Option<f64>,
}

impl Light {
    /// A blackbody lamp of `watts` at `colour_k`, falling on a mirror of this reflectance.
    pub fn new(name: impl Into<String>, watts: f64, colour_k: f64, finish: Spectrum) -> Light {
        let lamp = SpectralPower::new(
            Spectrum::blackbody(colour_k),
            Power::w(watts.max(0.0)),
            VISIBLE_RANGE,
        );
        // The coating sampled onto a spectrum the radiometry can integrate. This is the
        // discretisation step where a careless coupling loses energy, and the bus audit is
        // what would catch it.
        let optics = SurfaceOptics {
            reflectance: finish,
            transmittance: Spectrum::constant(0.0),
            diffuse: 0.0,
        };
        let absorptance = Spec::curve(
            (0..=60)
                .map(|i| {
                    let nm = 400.0 + i as f64 * 5.0;
                    (nm, optics.absorptance(Length::nm(nm)))
                })
                .collect(),
        );
        let reserve = f64::INFINITY;
        Light {
            name: name.into(),
            lamp,
            absorptance,
            reserve,
            saved: None,
        }
    }

    /// Give it a finite budget, so a run has an end and the ledger has something to say.
    pub fn with_reserve(mut self, joules: f64) -> Light {
        self.reserve = joules.max(0.0);
        self
    }

    /// What the surface takes out of the beam, in watts.
    pub fn absorbed_power(&self) -> Power {
        self.lamp.absorbed_by(&self.absorptance)
    }

    /// The fraction of the lamp's output that ends up as heat.
    ///
    /// The number the optics exists to produce: it depends on the whole overlap of the
    /// source's spectrum with the coating's, not on either alone.
    pub fn absorbed_fraction(&self) -> f64 {
        let total = self.lamp.total().to_si();
        if total > 0.0 {
            self.absorbed_power().to_si() / total
        } else {
            0.0
        }
    }

    /// Joules left to spend.
    pub fn reserve(&self) -> Energy {
        Energy::from_si(self.reserve)
    }
}

impl Domain for Light {
    fn name(&self) -> &str {
        &self.name
    }

    /// Light crosses the instrument in nanoseconds. Against a thermal timescale it is solved
    /// instantly and never subcycled — a solve is a solve.
    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let joules = (self.absorbed_power().to_si() * dt.to_si()).min(self.reserve);
        self.reserve -= joules;
        bus.publish(HEAT, joules);
        Ok(())
    }

    fn ledger(&self) -> Ledger {
        // An infinite reserve is not a number the audit can work with, and saying so is
        // better than reporting a finite lie. A scene wanting the books to close gives it
        // a budget.
        if self.reserve.is_finite() {
            Ledger::new().with(quantity::ENERGY, self.reserve)
        } else {
            Ledger::new()
        }
    }

    fn checkpoint(&mut self) {
        self.saved = Some(self.reserve);
    }

    fn restore(&mut self) {
        if let Some(r) = self.saved {
            self.reserve = r;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    /// What is left to spend.
    ///
    /// Not what arrives: a lamp on a mirror spends far more than it delivers, and the coating
    /// decides the difference. The consumer reports what it absorbed.
    fn readings(&self) -> Vec<Reading> {
        vec![Reading::new(
            &self.name,
            "reserve",
            self.reserve().to_si(),
            "J",
        )]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }
}
