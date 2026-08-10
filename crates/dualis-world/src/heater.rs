//! A heat source, written here rather than taken from the library.
//!
//! Deliberately here. The point of the coupled scene is to drive `Exchange` from *outside*
//! the workspace, and a publisher that lives in this crate is the only way to find out
//! whether [`Domain`] can be implemented by a consumer without reaching for anything the
//! library keeps to itself. Every other domain in a scene so far was one the library already
//! provided, which tested the constructors and nothing else.

use dualis::core::Reading;
use dualis::prelude::*;

/// An element that pays joules onto the bus out of a finite tank.
///
/// The tank is what makes the books closable. A source with an unlimited supply creates
/// energy from nothing every step, and the kernel is right to refuse that — so the choice is
/// between modelling where the energy comes from and turning the audit off. This models it.
pub struct Heater {
    name: String,
    watts: f64,
    /// Joules not yet spent. This *is* the ledger entry, which is what closes the books:
    /// what leaves here arrives somewhere else and is reported there.
    reserve: f64,
    saved: Option<f64>,
}

impl Heater {
    /// An element of `watts` with `reserve` joules to spend.
    pub fn new(name: impl Into<String>, watts: f64, reserve: f64) -> Heater {
        Heater {
            name: name.into(),
            watts: watts.max(0.0),
            reserve: reserve.max(0.0),
            saved: None,
        }
    }

    /// Joules left in the tank.
    pub fn reserve(&self) -> Energy {
        Energy::from_si(self.reserve)
    }
}

impl Domain for Heater {
    fn name(&self) -> &str {
        &self.name
    }

    /// Quasi-static: the element has no state that needs marching, only a budget that
    /// depletes. Subcycling it would divide the same payment into smaller ones for nothing.
    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let joules = (self.watts * dt.to_si()).min(self.reserve).max(0.0);
        self.reserve -= joules;
        bus.publish(HEAT, joules);
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

    /// What is left in the tank.
    ///
    /// The whole state of a source, and the number that says whether a run outlasted its energy
    /// — a scene that ends with a full tank never got going, and one that empties halfway has a
    /// flat second half nobody asked for.
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
