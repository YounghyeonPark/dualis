# Changelog

Notable changes, in the format of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This workspace follows [semantic versioning](https://semver.org/). It is `0.x`, so the API is
explicitly not stable and a minor bump may break you. The first consumer exists now, and it
has already found sixteen places it is awkward, eleven of which have been changed — see the
0.2.0 entry below.

Entries record what was *found* as well as what was added, because several of the more useful
changes here were corrections to a mistaken assumption rather than new features. The commit
messages carry the full account.

## [Unreleased]

### Added

- **Python bindings**, in `bindings/python`, as their own cargo workspace. `pip install` the
  wheel and `import dualis`: a `Simulation`, the library's heater, bar and lumped-mass domains,
  and the conservation audit as a `dualis.Violation` carrying `quantity`, `site`, `before`,
  `after`, `scale` and `tolerance` — addressable rather than a sentence to parse. A refused step
  does not move the clock.

  SI floats at the boundary with the unit in the parameter name, because the dimensional types
  are a compile-time thing Python cannot have, and a runtime wrapper would cost per operation to
  catch an error a Python caller does not make. What crosses instead is the audit.

  A domain cannot be written in Python yet, and the reasons are in its README. Enough to *run and
  audit* coupled physics, not enough to *extend* it.

  Separate workspace because pyo3 brings about fifteen crates and links libpython, and the
  library's twelve external dependencies, its `deny.toml` allow-list and its WebAssembly jobs are
  promises that should not have to accommodate a Python extension. Verified rather than assumed:
  the library workspace still resolves to exactly twelve external crates. An abi3 wheel, so one
  build serves 3.10 upward; CI builds it, installs it and runs its six tests, each against a
  number computed in the test file rather than read off the simulation.

## [0.2.0] — 2026-08-08

Breaking, and almost nothing broke: `&str: Into<String>` meant not one existing call site
changed. Everything in it came from the workspace acquiring its first consumer, and then from
two subagents built out of what that consumer taught.

The headline is not the ergonomics. It is that a first-order accuracy defect in the kernel's own
scheduler — in the schedule chosen *for* accuracy — was found by an application comparing a
coupled run against the closed form of its own recursion, having survived every test the library
had while the conservation audit reported clean to 1e-12.

### Fixed — the kernel scheduler

- **`Schedule::Multirate` does not refine a coupled quantity.** `sweep` steps one domain to
  completion before the next, so a quasi-static publisher puts a whole outer step's joules on the
  bus at once and a subcycling consumer takes all of them on its *first* substep. Refining the
  substep therefore does not move the answer: the error is first order in the **outer** step and
  independent of the substep. Measured on a lumped plate under a lamp — 26.2% low at a 300 s
  outer step, 13.8% at 150 s, 7.1% at 75 s, whatever the substep count — and at 300 s the
  schedule chosen *for* accuracy is worse than `Staggered`, with the errors on opposite sides.

  Every one of those runs passes the conservation audit at ~1e-12. The total that crossed is
  exactly right; only its distribution in time is wrong, and a `Ledger` has no representation for
  *when*. The time-domain twin of the reason `audit_transfers` became a per-face check.

  **Fixed** by `Exchange::take_share(channel, dt)`, new in the kernel: `advance` tells the bus
  what interval the sweep covers and a subcycling consumer asks for its substep's share instead
  of the lot. `Bar1D` and `LumpedMass` use it. Error at a 300 s outer step went from 1.89 K to
  0.304 K — from the worse of the two schedules to fourteen times better than the alternative.

  The share is apportioned against the time *remaining* rather than the whole interval, which is
  what leaves the channel exactly empty: `A·dt/T` with both reduced keeps `A/T`, so the last
  substep takes the remainder. Against the whole interval, `n` shares strand `O(n·ε·A)` where
  `audit_transfers` uses an absolute tolerance.

  Built rather than deferred because the recommendation to wait for a second consumer did not
  survive checking: `04-heater-and-bar` already pairs a quasi-static heater with a bar that
  subcycles hard, so a shipped scene had the defect.

### Fixed — what the two new subagents found

- **`Violation::at`'s cases printed an ungrammatical sentence.** They carry a *message* in
  `quantity` and `Display` had no branch for them, so the first error a consumer ever saw read
  "substance has no heat capacity is not conserved at plate: inf". A third branch.
- **`Report` could not be named without a module path**, though it is what `advance` returns.
  Added to `dualis-core`'s root and the prelude.
- **`Substance` was in the prelude and unbuildable from it**: `bulk` leaves `thermal: None`,
  which `LumpedMass` refuses to step, and the three types needed to supply one were not exported.
  `ThermalProps`, `MechanicalProps`, `AcousticProps`, `ThermalConductivity` and
  `ThermalExpansion` now are.
- **The scene format discarded unknown keys.** `serde` does that by default, which is right for
  a wire protocol and wrong for a saved document: `main.rs`'s own built-in scene kept the
  pre-`release` spelling for two commits and nothing failed, because the keys were dropped and
  the field fell back to its `Default`. Editing them was a no-op that reported success.
  `deny_unknown_fields` on `Scene`, `DomainSpec`, `Release`, `ScheduleSpec` and `Boundary`.
- **The round-trip test could not have caught that.** Both sides of its byte assertion were
  serialiser output, so the hand-written spelling never entered any comparison. It now parses the
  text a person would type and requires each stated value to survive.
- **An unrecognised `finish` produced a silent zero-watt lamp.** The early return also skipped
  `with_reserve`, so the reserve stayed infinite, so `Light::ledger` reported nothing, so the
  audit had nothing to compare — and the scene ran green at `conservation_tolerance(0.0)`, the
  strictest setting expressible, with the lamp doing nothing. One character, `aluminium` against
  `aluminum`. `DomainSpec::build` is fallible now and names the finishes it knows.
- **Two domains could share a name.** `Simulation::domain` takes the first match, so the second
  was never sampled and the first was drawn twice under the second's label and geometry — a
  500 °C bar reported as 20 °C, twice. `World::build` refuses it.
- **A scene whose every domain lacked a field wrote a zero-byte SVG and exited 0**, and "0 KiB"
  could not distinguish that from a legitimate 937-byte strip. The report is a row per *domain*
  now, naming the ones with nothing to draw, with the run-wide extremum beside the final value —
  because a ball that bounced half a metre and one that never moved both end at zero. An empty
  picture is refused rather than written.
- **One colour scale spanned panels of different units**, so a 1 Pa room beside a 7546 m/s orbit
  rendered as an empty bordered square while the numbers beside it looked fine. One extent per
  panel now, still shared across frames.

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
- **`NBody`, `TreeNBody`, `RigidBody` and `Rolling` answer `as_any`.** All four returned the
  default `None`, so `Simulation::domain_as` could not reach any of them and a renderer got
  nothing back. The failure was a picture with no bodies in it — not an error, not a violation,
  just an empty frame, which is the least debuggable outcome there is. Optics, thermal,
  acoustic and molecular had all opted in because tests inside the workspace had reached for
  them; mechanics had not, because none had.
- **`Bar1D::exposing` takes `impl Into<String>`.** The sweep behind `Domain::name` had matched
  on the parameter being called `name` and missed this one, so a boundary name — the thing two
  domains have to agree on, and therefore exactly the kind of thing that comes from a file —
  was still compile-time only. The kernel's own worked examples were also still teaching
  `fn name(&self) -> &'static str`.

### Added

- **`dualis-world`** — the first consumer, and not published. Scenes described as JSON, built
  into a `Simulation`, run, and drawn as an SVG filmstrip with no dependency. It exists to use
  the SDK from outside rather than to be a good application, and it reports what that was like
  in `crates/dualis-world/FRICTION.md`: twelve findings, seven fixed. Five are recorded rather
  than actioned, and the reasons differ: sharing the examples' SVG plotting would mean
  committing to a public drawing API in a workspace whose scope excludes rendering; validating
  a scene's schedule against its domains at build time is something an application can already
  do with `Domain::max_stable_dt`; and a duplicated face count needs no format change because
  `Exchange::publish_on` already refuses the mismatch by name. The report also says what it
  does *not* cover.

  **Ten scenes ship**, in `crates/dualis-world/scenes/`, covering all five domains,
  and CI runs every one through the real binary as well as the test harness. Three acoustic —
  two room modes and a clap that reflects off all four walls; two thermal, which are the same
  heat told over a plain channel and over a shared boundary and are the argument for
  `Interface` in one picture; two mechanical — four satellites and a bouncing ball whose
  dashpot heat a thermal lump takes; and two molecular, the same 108 atoms at `T* = 0.15` and
  `T* = 1.4`, which side by side are melting. Each has one number asserted, chosen to be a
  property of the physics rather than of the file, and a scene that ships without a claim
  fails the test rather than passing quietly.

  A room can be released as a Gaussian pulse now and not only as a mode, and a panel can hold
  bodies as well as a sampled field — an orbit is a countable number of things at places, and
  rasterising one would invent a continuum it does not have.

  **Bodies are drawn in three dimensions.** The physics always had them: `NBody`,
  `ContactSystem` and `Fluid` all carry `DVec3`, and flattening to a plane was the renderer's
  simplification. They are projected axonometrically now, sorted back to front, with radius
  growing toward the viewer and colour mixed toward the plate for distance — all three are
  needed or the picture is flat however true the coordinates are. Periodic cells get a
  wireframe, because they are a real boundary. The orbit scene tilts its satellites out of one
  plane, which is what makes the third axis carry anything.

  **Optics has a scene**, and it is the fifth `Domain` written outside the library: a
  blackbody lamp on an aluminium mirror whose reflectance falls off in the blue, so the
  colour temperature decides how much of a hundred watts becomes heat. Asserted as the
  difference between 2800 K and 6500 K rather than as one number, because a flat reflectance
  would make `Spectrum` and `SurfaceOptics::absorptance` an expensive way to multiply by a
  constant.

  The scene format couples both ways. On a plain channel, a heater defined **in the
  application** publishes joules and a bar takes them. Over a shared boundary, a beam
  publishes a Gaussian `Flux` onto an `Interface` the bar exposes, and the bar ends up hotter
  in the middle than at the ends by a ratio the scene never states. Both are audited at 1e-9,
  and the totals are checked against `Q / ρVc_p` computed outside the library.

  That closes the gap the report itself had flagged — until then no `publish`, `take`,
  `Exchange`, `Interface` or `Flux` appeared anywhere in the consumer, so the part of the API
  this workspace exists for had only ever been driven from inside. Writing both domains from
  outside needed nothing beyond `dualis::prelude`.
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

[Unreleased]: https://github.com/YounghyeonPark/dualis/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.2.0
[0.1.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.1.0
