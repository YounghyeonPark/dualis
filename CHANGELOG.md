# Changelog

Notable changes, in the format of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This workspace follows [semantic versioning](https://semver.org/). It is `0.x`, so the API is
explicitly not stable and a minor bump may break you. The first consumer exists now, and it
has already found three places the API should probably change — see `dualis-world` below.

Entries record what was *found* as well as what was added, because several of the more useful
changes here were corrections to a mistaken assumption rather than new features. The commit
messages carry the full account.

## [Unreleased]

Version bumped to 0.2.0 in the tree; **not yet published**. crates.io still carries 0.1.0, so
`dualis = "0.1"` is what a consumer gets and none of the changes below are in it.

### Changed — breaking

The first consumer went looking for the API's shape and found the same decision three times:
`&'static str` and `impl Domain` are free when every name is a literal in a test and every
domain type is known when you compile. They are not free for an application. All three are
reversed, and the cost was one order of magnitude smaller than the argument for keeping them —
**no existing call site changed**, because `&str: Into<String>`.

- **`Domain::name` returns `&str`.** Domains store a `String`; constructors take
  `impl Into<String>`, so `Bar1D::new("bar", ..)` and a name read out of a scene file both
  work. `Interface` followed, and `Exchange`'s spatial channel map is keyed by an owned
  interface name.
- **`Report::substeps` is `Vec<(String, u32)>`.** The only breakage across 349 tests: five
  comparisons against string literals. A report outlives the borrow it would otherwise hold.
- **`Simulation::with_boxed(Box<dyn Domain>)`**, plus `impl Domain for Box<dyn Domain>`. The
  simulation always stored boxes; now a caller who built one at run time can hand it over.
- **`Domain::as_field() -> Option<&dyn ScalarField>`** and `Simulation::field(name)`, both
  opt-in and `None` by default, in the style of `as_any`. `ScalarField` was written as the
  interface a visualiser reads a simulation through and was unreachable from `&dyn Domain`, so
  the workspace's own visualiser downcast to concrete types instead — precisely what the
  interface existed to avoid. It no longer names a single domain type.
- **`Room` is in the prelude**, which it should always have been; `Tube` already was.
- **`Bar1D::exposing` takes `impl Into<String>`.** The sweep behind `Domain::name` had matched
  on the parameter being called `name` and missed this one, so a boundary name — the thing two
  domains have to agree on, and therefore exactly the kind of thing that comes from a file —
  was still compile-time only. The kernel's own worked examples were also still teaching
  `fn name(&self) -> &'static str`.

### Added

- **`dualis-world`** — the first consumer, and not published. Scenes described as JSON, built
  into a `Simulation`, run, and drawn as an SVG filmstrip with no dependency. It exists to use
  the SDK from outside rather than to be a good application, and it reports what that was like
  in `crates/dualis-world/FRICTION.md`: eight findings, six fixed. Two are declined in
  writing rather than actioned — sharing the examples' SVG plotting, which would mean
  committing to a public drawing API in a workspace whose scope excludes rendering, and
  validating a scene's schedule against its domains at build time, which an application can
  already do with `Domain::max_stable_dt`. The report also says what it does *not* cover.

  The scene format now couples: a heater defined **in the application** publishes joules and
  a bar takes them, with the kernel auditing the crossing at 1e-9 and the bar's rise checked
  against `6 J / ρVc_p` computed outside the library. That closes the gap the report itself
  had flagged — until then no `publish`, `take` or `Exchange` appeared anywhere in the
  consumer, so the part of the API this workspace exists for had only ever been driven from
  inside. Writing a `Domain` from outside needed nothing beyond `dualis::prelude`.
  Excluded from the wasm, determinism and 1.78 jobs, which are promises the *library* makes to
  the people who depend on it.

### Found and fixed

- **A first-order startup error in `Room`, and in `Tube`.** A mode released from rest follows
  `|cos(2πft)|`, but the gap converged at first order against grid resolution where the
  scheme's interior is second. `released_from` left the velocity at `t = 0`; a staggered
  leapfrog carries it at `t = −h/2`, so the first velocity update travelled a whole step where
  it was owed half. `O(h)`, permanent, and `h` follows `dx` through the CFL condition.

  Fixed: the first velocity update takes half a step. Second order now, and the worst
  departure over 20 ms fell from 0.0528 to 0.00238 at 31 cells. Found by the workspace's own
  application checking itself against a closed form, because nothing inside the library was
  checking a rate.

  Two things came out of the fix. **A test had turned the bug into the specification** — one
  step from rest was asserted to move the pressure by `h²c²∇²p`, where Taylor gives `½h²c²∇²p`
  since `ṗ(0) = 0`; the test was missing the half because the scheme was, and `Tube` had the
  same pair. And **the old startup conserved energy exactly while the correct one does not**:
  with `v = 0` read as the half-step value, `Σ∇·(p∇p) = 0` at a rigid wall makes the first
  step's energy change cancel to the last bit. Starting correctly breaks that by `O(h²)` at
  the first step only — 0.42% at 31 cells, quartering on refinement, and 1e-15 thereafter. The
  old code had bought exact bookkeeping by making the scheme first order, which is this
  workspace's own documented trap appearing a second time in the same crate. Energy is now
  reported against the released state as its datum, with the difference available from
  `Room::startup_adjustment` and bounded at 25% so a real first-step bug cannot hide there.

- **The API is comfortable only when the set of domains is known at compile time.**
  `Simulation::with` takes `impl Domain` and there is no `impl Domain for Box<dyn Domain>`;
  domain names are `&'static str`, so a name from a file has to be leaked; and a renderer
  cannot get a `&dyn ScalarField` from a `&dyn Domain`, so it downcasts and knows every domain
  by name — which is what `ScalarField` existed to avoid. Three symptoms of one position. It
  may be the right position, but nobody chose it.

## [0.1.0] — 2026-08-07

First release. All eight crates published to crates.io together.

Everything below was in this release. Two notes on the publish itself, since they cost time
and are not obvious from the outside:

- crates.io requires a **verified** email address, not merely a registered one, and reports
  its absence as a `400` at the first upload rather than at login.
- New crates are rate limited to a burst of five, then roughly one every ten minutes. A
  workspace of eight publishes five, stops, and has to be resumed — so `cargo publish
  --workspace` is not atomic and a partial publish is the normal outcome, not a fault.

### Added

- **`dualis-units`** — dimensional analysis with the SI exponents in the type, so `Length +
  Time` does not compile. Const-generic `Qty` and `QVec3`, macro-generated products, and
  unit-bearing constructors as the only place a factor of a thousand may appear.
- **`dualis-core`** — the kernel. Conservation as an audit (`Ledger`, `audit`, `Violation`),
  fixed-step integrators including `velocity_verlet`, multi-domain scheduling with
  quasi-static, multirate and iterative coupling, deterministic sampling through
  `Rng::for_index`, an accurate discrete Fourier transform, closed-form rigid motion, scalar
  and vector fields, and shared boundaries (`Interface`, `Flux`).
- **`dualis-optics`** — spectral radiometry, surface optics with Fresnel and coatings,
  Sellmeier and Cauchy dispersion, ray geometry, Airy diffraction and the ideal MTF, Zernike
  wavefronts and aberrated PSFs, angular-spectrum propagation, partial coherence, and a
  detector with the four noises that come with counting photons.
- **`dualis-thermal`** — lumped masses and explicit one-dimensional conduction, with radiative
  and convective loss.
- **`dualis-mechanics`** — exact N-body, Barnes-Hut with a quadrupole, penalty contact with
  Coulomb friction, rolling, and rigid-body rotation with Euler's equations.
- **`dualis-acoustic`** — the linear wave equation on a staggered grid, in a tube and in a
  room, with impedance boundaries.
- **`dualis-molecular`** — Lennard-Jones fluids in periodic boxes, cell lists, a Langevin
  thermostat, virial pressure, and radial distribution functions.
- **Five examples**, each of which asserts its numbers and is run by CI. Give one a path and
  it writes an SVG; the plotting has no dependency. Two further examples are checks rather
  than showcases: `agents_quickstart` and `readme_check`.
- **`AGENTS.md` and `CLAUDE.md`** — the API on one page for a consumer, and the gate and
  conventions for a contributor. Written after an AI agent looked for dualis on `PATH`, as a
  Python package, and in a consuming repository, found it in none of the three, and used
  MuJoCo instead. That is a distribution failure and not a documentation one, but the
  quickstart it now lands on is `examples/agents_quickstart.rs`, which CI runs, so it cannot
  drift from the library the way a hand-written snippet does.
- **CI** across Linux, macOS, Windows, two WebAssembly targets and Rust 1.78, with formatting,
  clippy at `-D warnings`, rustdoc at `-D warnings`, licence and advisory checks, and the
  examples.

### Fixed

- **The acoustic wall weighting.** A pressure sample sits *on* a wall, so it owns half a cell,
  and both `Tube` and `Room` divided its divergence by the whole `dx`. Every mode read low —
  1.4% on an 89-cell room — and the scheme converged at first order despite a second-order
  interior. Found by the *rate* of convergence rather than by its size.
- **`Tube`'s absorbing ends** consequently needed their own step limit, `Z·dx/2ρc²`. At the
  full CFL limit the corrected boundary inverts a wave instead of absorbing it: stable, silent
  and wrong.
- **`Bar1D`'s enthalpy reference.** Measured from absolute zero, the bar in these tests holds
  1.42 kJ, so a millijoule arriving is a change in the seventh significant figure and the
  audit's relative check asked for precision the arithmetic had thrown away — a floor that
  grew from 1.6e-12 J at 41 cells to 7.3e-12 J at 161. Measured from the initial temperature,
  the number being summed *is* the change.
- **A statistical test that passed on one seed.** The ideal-gas check asserted that doubling
  the density doubles the departure from `PV = Nk_BT`; across four seeds the ratio came out
  1.35, 1.61, 2.34 and 2.92, averaging to 2.06. It now averages over seeds, because the fix for
  a noisy statistical test is more samples and not a wider tolerance.
- **Licence texts** now ship inside each crate. The packages declared `MIT OR Apache-2.0` and
  contained neither, which `cargo package` does not warn about.
- **Five tolerances that were loose rather than wrong.** A tail-correction ratio checked
  against the `rc⁻³` power law inside 0.02, where the exact ratio differs from that law by
  0.0188 — a known discrepancy filling 94% of the budget. It is now checked against the closed
  form to 1e-12, with the power law asserted separately as the limit it actually is. The
  Langevin settling test averaged one seed and allowed 5%; the seed-to-seed spread is 0.96%,
  so it now averages four seeds and allows 2% — more samples, and a *tighter* tolerance. Two
  gradients expected to be zero were checked with a relative comparison against zero. A
  crystal's negative pressure was asserted as `pressure.max(0.0)`, which an identically zero
  virial satisfies; it is now a band that excludes zero. And the bouncing ball's energy
  handover compared the heat that crossed against `start_energy - 0.0f64.max(0.0)` — the final
  energy was never computed, so the equality its comment promised was a one-sided bound.
  `ContactSystem` now answers `as_any` so the equality can be checked; publishing 5% extra
  heat fails it.

### Documented

- Every public item, with `#![deny(missing_docs)]` in all eight crates so it stays that way.
  244 items had no doc comment; the concentration was in `dualis-units`, which is the API
  nobody can avoid.
- A `compile_fail` doctest proving `Length + Time` does not build — the workspace's reason for
  existing, previously asserted only in prose.

[Unreleased]: https://github.com/YounghyeonPark/dualis/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.1.0
