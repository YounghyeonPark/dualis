//! A source that heats *where it lands*, over a shared boundary.
//!
//! The spatial counterpart of [`Heater`](crate::heater::Heater), and the reason this file
//! exists separately: a plain channel carries an amount, and an [`Interface`] carries an
//! amount **and a place**. That distinction is the newest part of the kernel and had no
//! consumer outside the workspace until this crate grew one.
//!
//! # What a consumer has to get right, and what the kernel will not let it get wrong
//!
//! Both sides of a spatial coupling must agree on the boundary's *discretisation* — the same
//! name and the same number of faces. Nothing derives one from the other: the bar builds its
//! interface from its own cell count and this builds a matching one from the scene, so the
//! agreement is stated twice and could be stated inconsistently.
//!
//! It cannot be stated inconsistently *and get away with it*. [`Exchange::publish_on`] refuses
//! a flux whose face count differs from the interface's, and says both numbers. Silently
//! padding or truncating would put energy on the wrong part of the boundary, which is worse
//! than losing it — losing it the audit would catch on its own.

use dualis::core::Reading;
use dualis::prelude::*;

/// A beam that pays joules onto a boundary with a Gaussian profile across it.
pub struct Beam {
    name: String,
    boundary: Interface,
    watts: f64,
    reserve: f64,
    /// Waist as a fraction of the boundary's span, so the shape is scale-free and the scene
    /// does not have to know the bar's length in two places.
    waist_fraction: f64,
    saved: Option<f64>,
}

impl Beam {
    /// A beam of `watts` with `reserve` joules to spend, landing on `boundary`.
    pub fn new(
        name: impl Into<String>,
        boundary: Interface,
        watts: f64,
        reserve: f64,
        waist_fraction: f64,
    ) -> Beam {
        Beam {
            name: name.into(),
            boundary,
            watts: watts.max(0.0),
            reserve: reserve.max(0.0),
            waist_fraction: waist_fraction.max(1e-6),
            saved: None,
        }
    }

    /// Joules left to spend.
    pub fn reserve(&self) -> Energy {
        Energy::from_si(self.reserve)
    }

    /// The boundary this publishes onto. Both sides have to agree on its face count.
    pub fn boundary(&self) -> &Interface {
        &self.boundary
    }
}

impl Domain for Beam {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let joules = (self.watts * dt.to_si()).min(self.reserve).max(0.0);
        // `exp(-2(u - ½)²/w²)` in the boundary's own cumulative-area coordinate, which is
        // how the physics writes a Gaussian and not how a grid indexes one. The kernel
        // normalises the shape so the total is exactly `joules` — the profile decides where,
        // never how much, which is what keeps the audit an equality.
        let w = self.waist_fraction;
        let flux = Flux::profiled(joules, &self.boundary, |u| {
            (-2.0 * ((u - 0.5) / w).powi(2)).exp()
        });
        bus.publish_on(&self.boundary, HEAT, &flux)?;
        self.reserve -= joules;
        Ok(())
    }

    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.reserve)
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
