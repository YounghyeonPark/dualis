# Contributing

Issues and pull requests are welcome. This file is short on process and long on the two or
three conventions that are unusual enough that a contributor would otherwise discover them by
having a change sent back.

## Run what CI runs

```sh
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
cargo test --locked --workspace --release
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --no-deps
cargo deny check                       # needs `cargo install cargo-deny`
```

And the examples, which are tests that print:

```sh
for e in beam_hot_spot airy_pattern detector_snr room_modes melting agents_quickstart readme_check; do
  cargo run --locked --release --example "$e" || break
done
```

CI additionally builds on Rust 1.78, builds for `wasm32-unknown-unknown`, and runs the whole
suite under `wasm32-wasip1` with wasmtime. Those three catch different things and none of them
is decoration — see below.

## The conventions that are not obvious

### Check against a closed form, never against another implementation

A test that compares a function to a second implementation of the same idea passes when both
are wrong in the same way, which is the usual way to be wrong. So every claim here is checked
against something structurally different: an analytic result, an exactly known limit, a
quantity computed by a different route, or a convergence *rate*.

That last one is worth naming, because it caught a real defect. The acoustic domain read every
room mode 1.4% low, which looks exactly like discretisation error and invites a loose
tolerance. What showed it was not discretisation was that refining the grid **halved** the
error instead of quartering it — a second-order scheme converging at first order means the
boundary condition is wrong. No single run could have said that.

### A tolerance has to be earned

`1e-9` because the arithmetic earns it, or `2e-2` because a penalty contact is a non-smooth
potential and the symplectic bound does not apply. Never a number chosen because it made the
test pass. If a tolerance is loose, the comment says which physical or numerical effect is
using up the budget.

Judge a residual against a **scale**, not against itself. A correct system's net conserved
quantity is often exactly zero, so a relative tolerance on it is meaningless; `Ledger` records
the largest contribution beside the total for this reason, and `Violation` carries it.

### Nothing is random, and nothing consults a clock

Fixed steps, no wall clock, no unordered reduction. Randomness comes from
`Rng::for_index(seed, index)`, which hashes a work item into its own stateless stream — so
results are identical whatever order the work is done in, and a parallel run agrees with a
serial one bit for bit. `rng::tests::the_stream_is_pinned` fixes the generator's output as a
constant, and **changing that constant is never the fix**.

This is also why the `wasm32-wasip1` job exists. It runs the pinned digest and every
closed-form comparison under a different target, so a platform that rounded differently would
fail there rather than quietly producing different physics.

### The kernel must never depend on a domain

`dualis-core` knows about conservation, integration, scheduling and boundaries. It knows
nothing about light, heat, motion, sound or matter. If a new physics needs the kernel changed,
the kernel was wrong — with one narrow exception, which is that the *coupling mechanism* itself
can turn out to be under-specified. That happened once and the README explains why it was not a
violation of the rule.

Domain crates do not depend on each other either. Rustdoc enforces this in a way worth
knowing: an intra-doc link from `dualis-acoustic` to `Bar1D` does not resolve, because the
dependency is not there, and `-D warnings` turns that into a build failure.

### Every public item is documented

`#![deny(missing_docs)]` in all ten crates. A one-line summary that names the unit is enough
for a constructor; anything with a trap in it should say what the trap is.

### Say what was wrong, not only what changed

Commit messages here record the mistake as well as the fix, because the mistake is usually the
more useful half. Several of them document a wrong assumption of the author's that the tests
caught — which grid a boundary condition belongs to, why a statistical test passed on one seed
and not on three others, why an obvious-looking neighbour shell turned out to discriminate
nothing. If you find something like that, write it down.

## What a good pull request looks like

- One idea. A fix and a refactor in the same change are two changes.
- Tests that would fail without it, checked against something independent.
- No new dependency without saying what it buys. The workspace has twelve, three of which
  reach a *published* artifact — the unpublished application links four more — and `deny.toml`
  gates the licences.
- `cargo fmt` clean and `clippy -D warnings` clean.

## Reporting something

An issue with a failing case is worth more than a description. If it is a physics
disagreement, the most useful form is: what the code produced, what the closed form says, and
how the two diverge as something is refined.

## Licence

By contributing you agree that your work is dual-licensed under MIT and Apache-2.0, as the
rest of the workspace is. The README's licence section states this in the standard wording.
