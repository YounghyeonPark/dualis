//! Shared plumbing for the examples.
//!
//! Cargo builds `examples/*.rs` and `examples/*/main.rs` as binaries; a directory with no
//! `main.rs` is not a target, so this one is invisible to the build and reachable only
//! through `mod common;` in each example. That is the standard way to share code between
//! examples without inventing a crate for it.

#![allow(dead_code)]

pub mod svg;

use std::path::Path;

/// Where an example writes, taken from the command line.
///
/// With no argument an example prints its numbers and asserts them; with a path it also
/// writes a picture. CI runs them without a path, so the assertions are checked on every
/// commit and no generated file is ever committed — a picture in the repository would go
/// stale the first time the physics changed underneath it.
pub fn output_path() -> Option<String> {
    std::env::args().nth(1)
}

pub fn write(path: &str, contents: &str) {
    if let Some(parent) = Path::new(path).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).expect("could not create the output directory");
        }
    }
    std::fs::write(path, contents).unwrap_or_else(|e| panic!("could not write {path}: {e}"));
    println!("wrote {path}");
}

/// Report a computed value against what it should be, and fail loudly if it is not.
///
/// Every example is a test as well as a demonstration. An example that quietly produced
/// nonsense would be worse than no example, because it reads as a claim that the library
/// works — so each number printed here has been checked against a closed form or an
/// independent calculation, and the check runs whether or not anyone is looking at the
/// picture.
pub fn check(label: &str, value: f64, expected: f64, rel_tol: f64, unit: &str) {
    let scale = expected.abs().max(value.abs()).max(f64::MIN_POSITIVE);
    let error = (value - expected).abs() / scale;
    println!("  {label:<44} {value:>12.4} {unit:<8} (expected {expected:.4}, off by {error:.2e})");
    assert!(
        error <= rel_tol,
        "{label}: got {value} {unit}, expected {expected} — relative error {error:.3e} exceeds {rel_tol:.3e}"
    );
}

/// Report a value that should be zero, judged against a scale supplied by the caller.
///
/// [`check`] cannot do this: with an expected value of zero every nonzero result is a 100%
/// relative error, so the tolerance has nothing to mean. It is the same trap the kernel's
/// own [`audit`](pantometry_core::audit) fell into — a correct system's net is often exactly
/// zero, which is why [`Ledger`](pantometry_core::Ledger) records the largest contribution
/// alongside the total. The scale here plays that part: for a conservation residual it is
/// the energy that actually crossed, not the fraction that failed to.
pub fn check_zero(label: &str, value: f64, scale: f64, rel_tol: f64, unit: &str) {
    let error = value.abs() / scale.abs().max(f64::MIN_POSITIVE);
    println!("  {label:<44} {value:>12.3e} {unit:<8} (against a scale of {scale:.4}, off by {error:.2e})");
    assert!(
        error <= rel_tol,
        "{label}: {value} {unit} against a scale of {scale} — relative {error:.3e} exceeds {rel_tol:.3e}"
    );
}

/// Report a value that has no closed form to check, only a range it must lie in.
pub fn check_between(label: &str, value: f64, lo: f64, hi: f64, unit: &str) {
    println!("  {label:<44} {value:>12.4} {unit:<8} (between {lo} and {hi})");
    assert!(
        value >= lo && value <= hi,
        "{label}: {value} {unit} is outside [{lo}, {hi}]"
    );
}

pub fn heading(s: &str) {
    println!("\n{s}");
    println!("{}", "-".repeat(s.len()));
}
