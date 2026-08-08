//! Conservation, as a thing a process must answer for rather than a property it
//! is trusted to have.
//!
//! `SurfaceOptics` stores reflectance and transmittance and computes absorptance
//! as the remainder, so a surface cannot return more light than reached it. That
//! is the right idea and the wrong scope: it protects one quantity at one kind of
//! boundary. Momentum in a collision, charge across a junction, mass through a
//! pipe and energy across a coupling interface are all the same problem, and all
//! of them are places a simulation can quietly manufacture something.
//!
//! So a process reports what it holds, before and after, and the difference is
//! checked. A [`Ledger`] is that report; [`audit`] is that check; a [`Violation`]
//! names what went missing and where.
//!
//! # Relative to what, exactly
//!
//! Floating-point arithmetic loses the low bits of every sum, so no real
//! integrator conserves anything exactly and a test for exact equality would fail
//! on correct code. The loss is bounded and relative — but relative to the wrong
//! thing if one is not careful, and this is the trap:
//!
//! A well-formed system's conserved total is often **exactly zero**. One domain
//! holds a debt of 28.9 J and another holds 28.9 J of asset, and the sum is nothing.
//! Comparing the residual against that sum makes every rounding error a 100%
//! relative error, and the audit fires on correct code.
//!
//! So a [`Ledger`] entry records the largest magnitude that went into it as well as
//! the total, and [`audit`] judges the change against *that*. Rounding error scales
//! with the size of the numbers being added, not with the size of their sum, and
//! this is the version of the tolerance that says so.

use std::collections::BTreeMap;
use std::fmt;

/// One quantity's books: the net total, and the size of the entries it came from.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Entry {
    total: f64,
    /// Largest single contribution, which is the scale rounding error lives on.
    scale: f64,
}

/// What a process claims to be holding, by quantity name, in SI base units.
///
/// A `BTreeMap` rather than a `Vec` or a `HashMap`: names must come out in one
/// order for the audit report to be reproducible, and a hash map's order is not
/// one order.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Ledger(BTreeMap<&'static str, Entry>);

/// The quantities worth naming. Strings rather than an enum, so a domain crate
/// can add its own without editing the kernel — but these spellings are the ones
/// [`audit`] will match across domains, so use them.
pub mod quantity {
    /// Joules. The channel four of the five domains publish and consume on.
    pub const ENERGY: &str = "energy";
    /// kg·m·s⁻¹. Audited component by component, which makes the smallest component the
    /// binding one — see [`audit`](super::audit).
    pub const MOMENTUM: &str = "momentum";
    /// Kilograms.
    pub const MASS: &str = "mass";
    /// Coulombs.
    pub const CHARGE: &str = "charge";
    /// A count, not an energy. A photon budget and a joule budget are different books, and
    /// a detector is where the two stop being interchangeable.
    pub const PHOTONS: &str = "photons";
}

impl Ledger {
    /// An empty ledger, holding nothing.
    pub fn new() -> Ledger {
        Ledger(BTreeMap::new())
    }

    /// Record a total. Repeating a name adds to it, since a domain made of parts
    /// reports the sum of its parts.
    pub fn with(mut self, quantity: &'static str, si_total: f64) -> Ledger {
        self.add(quantity, si_total);
        self
    }

    /// Add to a quantity's total, in SI base units.
    ///
    /// Also raises that entry's `scale` to the largest contribution seen, which is what makes
    /// a relative tolerance mean anything when the net total is near zero.
    pub fn add(&mut self, quantity: &'static str, si_total: f64) {
        let entry = self.0.entry(quantity).or_default();
        entry.total += si_total;
        entry.scale = entry.scale.max(si_total.abs());
    }

    /// The net total for a quantity.
    pub fn get(&self, quantity: &str) -> Option<f64> {
        self.0.get(quantity).map(|e| e.total)
    }

    /// The largest single entry that went into a quantity — the scale on which
    /// rounding error in its total should be judged.
    pub fn scale_of(&self, quantity: &str) -> Option<f64> {
        self.0.get(quantity).map(|e| e.scale)
    }

    /// Whether anything at all has been recorded.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Names and net totals, in a fixed order.
    pub fn quantities(&self) -> impl Iterator<Item = (&'static str, f64)> + '_ {
        self.0.iter().map(|(k, e)| (*k, e.total))
    }

    /// Sum of two ledgers — how a simulation totals its domains. The scales carry
    /// over as the larger of the two, so a big domain's rounding budget is not
    /// shrunk by being added to a small one.
    pub fn merged(mut self, other: &Ledger) -> Ledger {
        for (name, entry) in other.0.iter() {
            let mine = self.0.entry(name).or_default();
            mine.total += entry.total;
            mine.scale = mine.scale.max(entry.scale);
        }
        self
    }
}

/// A conservation law that did not hold.
#[derive(Clone, Debug, PartialEq)]
pub struct Violation {
    /// Which law: one of [`quantity`], or a domain's own name for it.
    pub quantity: String,
    /// Where it broke — a domain name, a coupling name, a wavelength.
    pub site: String,
    /// What the quantity was, in SI base units.
    pub before: f64,
    /// What it became.
    pub after: f64,
    /// What the discrepancy was measured against — the largest entry that went into
    /// the books, not the net total, since a correct system's net is often zero. Zero
    /// means "use the totals", which is what a non-conservation error does.
    pub scale: f64,
    /// The tolerance that was being applied, so the report says how badly.
    pub tolerance: f64,
}

impl Violation {
    /// For the cases that are not a before/after comparison at all: a surface
    /// specified to reflect more than it receives, an iteration that never
    /// converged.
    pub fn at(site: impl Into<String>, quantity: impl Into<String>, detail: f64) -> Violation {
        Violation {
            quantity: quantity.into(),
            site: site.into(),
            before: detail,
            after: detail,
            scale: detail.abs(),
            tolerance: 0.0,
        }
    }

    /// Absolute size of the discrepancy.
    pub fn error(&self) -> f64 {
        (self.after - self.before).abs()
    }

    /// Discrepancy as a fraction of the scale it was judged against.
    pub fn relative_error(&self) -> f64 {
        let scale = if self.scale > 0.0 {
            self.scale
        } else {
            self.before.abs().max(self.after.abs())
        };
        if scale == 0.0 {
            0.0
        } else {
            self.error() / scale
        }
    }
}

impl fmt::Display for Violation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // `Violation::at` builds the cases that are not a before/after comparison — a surface
        // specified to reflect more than it receives, a substance with no heat capacity, an
        // iteration that never converged. Those carry a *message* in `quantity`, not a
        // quantity, and reading it as one produced the first error a consumer ever saw from
        // this library: "substance has no heat capacity is not conserved at plate: inf".
        if self.tolerance == 0.0 && self.before == self.after {
            return write!(f, "at {}: {} ({})", self.site, self.quantity, self.before);
        }
        if self.before == self.after {
            return write!(
                f,
                "{} is not conserved at {}: {}",
                self.quantity, self.site, self.before
            );
        }
        let verb = if self.after > self.before {
            "created"
        } else {
            "destroyed"
        };
        write!(
            f,
            "{} {} at {}: {:.6e} became {:.6e}, a relative change of {:.3e} against a \
             tolerance of {:.3e}",
            self.quantity,
            verb,
            self.site,
            self.before,
            self.after,
            self.relative_error(),
            self.tolerance
        )
    }
}

impl std::error::Error for Violation {}

/// Compare two ledgers and name the first quantity that moved by more than
/// `rel_tol`, in the fixed order the ledger keeps its names.
///
/// The change is measured against the largest of the two totals *and* the largest
/// single entry either ledger recorded. See the module docs for why the totals alone
/// are not enough: a correct system whose books cancel to zero would otherwise turn
/// every rounding error into a 100% relative error.
///
/// Quantities present in only one of the two are treated as having been zero in
/// the other, so a process that starts reporting momentum halfway through gets
/// caught rather than excused.
pub fn audit(site: &str, before: &Ledger, after: &Ledger, rel_tol: f64) -> Result<(), Violation> {
    let mut names: Vec<&'static str> = before.0.keys().copied().collect();
    for name in after.0.keys() {
        if !before.0.contains_key(name) {
            names.push(name);
        }
    }
    names.sort_unstable();

    for name in names {
        let b = before.get(name).unwrap_or(0.0);
        let a = after.get(name).unwrap_or(0.0);
        let scale = b
            .abs()
            .max(a.abs())
            .max(before.scale_of(name).unwrap_or(0.0))
            .max(after.scale_of(name).unwrap_or(0.0));
        // Two numbers that are both denormal are equal for every purpose a
        // simulation has.
        if scale < 1e-300 {
            continue;
        }
        if (a - b).abs() / scale > rel_tol {
            return Err(Violation {
                quantity: name.to_string(),
                site: site.to_string(),
                before: b,
                after: a,
                scale,
                tolerance: rel_tol,
            });
        }
    }
    Ok(())
}

/// Something that can say what it is holding. Implemented by domains, and by
/// anything else whose books are worth checking.
pub trait Conserves {
    /// What this is currently holding.
    fn ledger(&self) -> Ledger;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_ledger_that_did_not_move_passes() {
        let before = Ledger::new()
            .with(quantity::ENERGY, 3.7)
            .with(quantity::MASS, 2.0);
        // Losing the last bits of a double is arithmetic, not a leak.
        let after = Ledger::new()
            .with(quantity::ENERGY, 3.7 + 4e-16)
            .with(quantity::MASS, 2.0);
        assert!(audit("test", &before, &after, 1e-12).is_ok());
    }

    /// The failure is named, sited and quantified, because "conservation failed"
    /// is not a debuggable message.
    #[test]
    fn a_leak_is_named_and_sited() {
        let before = Ledger::new().with(quantity::ENERGY, 1.0);
        let after = Ledger::new().with(quantity::ENERGY, 0.6);
        let err = audit("thermal", &before, &after, 1e-9).expect_err("40% is not arithmetic");
        assert_eq!(err.quantity, "energy");
        assert_eq!(err.site, "thermal");
        assert!((err.relative_error() - 0.4).abs() < 1e-12);
        let text = err.to_string();
        assert!(text.contains("destroyed"), "{text}");
        assert!(text.contains("thermal"), "{text}");
    }

    #[test]
    fn creating_something_reads_differently_from_losing_it() {
        let before = Ledger::new().with(quantity::PHOTONS, 1e6);
        let after = Ledger::new().with(quantity::PHOTONS, 1.5e6);
        let err = audit("optics", &before, &after, 1e-6).unwrap_err();
        assert!(err.to_string().contains("created"), "{err}");
    }

    /// A quantity that appears out of nowhere is a violation, not an exemption —
    /// this is the case a naive "compare the keys they share" audit would miss.
    #[test]
    fn a_quantity_absent_before_is_still_audited() {
        let before = Ledger::new().with(quantity::ENERGY, 1.0);
        let after = Ledger::new()
            .with(quantity::ENERGY, 1.0)
            .with(quantity::MOMENTUM, 5.0);
        let err = audit("contact", &before, &after, 1e-9).expect_err("momentum from nowhere");
        assert_eq!(err.quantity, "momentum");
        assert_eq!(err.before, 0.0);
    }

    /// Reports come out in one order, so a failing run names the same quantity
    /// every time rather than whichever one the hash landed on first.
    #[test]
    fn the_audit_order_is_fixed() {
        let before = Ledger::new()
            .with(quantity::MOMENTUM, 1.0)
            .with(quantity::CHARGE, 1.0)
            .with(quantity::ENERGY, 1.0);
        let after = Ledger::new()
            .with(quantity::MOMENTUM, 2.0)
            .with(quantity::CHARGE, 2.0)
            .with(quantity::ENERGY, 2.0);
        // Three laws broken at once; the alphabetically first is reported, every
        // time, on every platform.
        for _ in 0..8 {
            let err = audit("s", &before, &after, 1e-9).unwrap_err();
            assert_eq!(err.quantity, "charge");
        }
    }

    #[test]
    fn ledgers_merge_by_summing() {
        let a = Ledger::new().with(quantity::ENERGY, 1.5);
        let b = Ledger::new()
            .with(quantity::ENERGY, 2.5)
            .with(quantity::MASS, 1.0);
        let total = a.merged(&b);
        assert_eq!(total.get(quantity::ENERGY), Some(4.0));
        assert_eq!(total.get(quantity::MASS), Some(1.0));
    }

    /// Zero against zero is not a hundred-percent error.
    #[test]
    fn nothing_compared_to_nothing_is_fine() {
        let z = Ledger::new().with(quantity::ENERGY, 0.0);
        assert!(audit("s", &z, &z, 0.0).is_ok());
    }
}
