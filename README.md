# dualis

[![CI](https://github.com/YounghyeonPark/dualis-core/actions/workflows/ci.yml/badge.svg)](https://github.com/YounghyeonPark/dualis-core/actions/workflows/ci.yml)

Physics for simulated worlds — a kernel that knows nothing about any particular
physics, and domains built on it that do.

Extracted from [AiryTrace](https://github.com/YounghyeonPark/airytrace), which is
its first consumer.

## The crates

| Crate | |
| --- | --- |
| `dualis-units` | Dimensional analysis. SI quantities and vectors whose dimension lives in the type, so `Length + Time` does not compile |
| `dualis-core` | The kernel: conservation audits, fixed-step integrators, fields, multi-domain scheduling, deterministic sampling, closed-form rigid motion |
| `dualis-optics` | Light: spectral radiometry, surface optics, dispersion, ray geometry, diffraction |
| `dualis-thermal` | Heat: lumped masses, explicit conduction, radiative and convective loss |
| `dualis` | A facade over the four, and where the cross-domain integration tests live |

```text
dualis-units    no dependencies but glam and serde
dualis-core     depends on units
dualis-optics   depends on core        ─┐  these two do not know
dualis-thermal  depends on core        ─┘  about each other
dualis          depends on all of them
```

**The kernel must never depend on a domain.** If a new physics needs the kernel
changed, the kernel was wrong — that rule is what makes "add contact, add sound" a
matter of writing a crate rather than editing this one. Optics and thermal are the
proof: neither names the other, neither contains a line about coupling, and they
exchange energy anyway.

Storage is SI base units everywhere: metres, kilograms, seconds, kelvin. Millimetres
and nanometres are entry and exit forms, and the unit-bearing constructors are the
only place a factor of a thousand can hide.

## Two invariants

Both are now enforced rather than promised, and both survived being generalised past
optics.

**Nothing is created or destroyed without being noticed.** A `Ledger` is what a
process claims to hold, `audit` is the check, and a `Violation` names what went
missing and where. Energy crossing between domains goes over an `Exchange`, which
compares what was published against what was consumed — because each domain
conserving energy internally says nothing about the interface between two
discretisations of the same surface, and that interface is where it actually leaks.
`SurfaceOptics` still cannot be written down in a form that returns more light than
reached it.

**Nothing is random.** `Rng::for_index(seed, index)` hashes a work item into its own
stateless stream, so ray 10 000 can be drawn before ray 3 and neither result
changes. That is what lets a simulation be parallel *and* bit-reproducible — a
single shared generator loses reproducibility precisely when the run gets big enough
to need it. `rng::tests::the_stream_is_pinned` fixes the generator's output as a
constant, and changing that constant is never the fix.

## Running several domains at once

Domains do not agree on how big a step is: an explicit FDTD solver on a nanometre
grid is stable to about 10⁻¹⁷ s, heat conduction to 10⁻⁹ s, rigid contact to 10⁻⁴ s,
and a thermal drift that defocuses an instrument plays out over seconds. Stepping
everything at the smallest limit integrates the slow domains ten billion times for
nothing. Three mechanisms deal with that:

- `Kind::QuasiStatic` — a domain with no state to roll forward, re-solved on demand
  rather than stepped. Light crosses an instrument in nanoseconds; against a thermal
  timescale that is zero, so optics is never integrated. The largest single saving
  available.
- `Schedule::Multirate` — each evolving domain takes as many equal substeps as its
  own stability limit needs. Integer counts from a fixed limit, never adaptive, so
  two runs follow the same arithmetic path.
- `Schedule::Iterative` — repeats the sweep until the residuals settle. Necessary
  because a staggered coupling can be unstable no matter how small the step is; the
  standard example is fluid-structure interaction at comparable densities. Failing
  to converge is a `Violation`, not a result — an unconverged coupling produces
  numbers that look like physics.

## Closed form where there is one

`Motion` is a function of `t`: ask for the world at 0.7 s and you get it, without
having computed 0.6 s first. An exposure can therefore be sampled at seven instants
for motion blur, and frame 7 of a recording does not depend on frame 6.

That is not available in general — three bodies under gravity have no closed form,
and neither do contact or heat — so `Integrator` rolls those forward instead, under
three rules that keep the reproducibility: fixed steps, no wall clock, ordered
reduction. And for anything conservative it is `velocity_verlet` rather than `Rk4`:
RK4 is more accurate per step and steadily *dissipates*, while the symplectic method
is second order and holds its energy within a bound. The test module proves that on a
harmonic oscillator against the closed-form energy.

## A taste

```rust
use dualis_optics::{Material, Spectrum, SurfaceFinish, fresnel_reflectance};
use dualis_units::Length;

// Reflectance is not a setting. It follows from the refractive indices.
let bk7 = Material::from_catalog("N-BK7").unwrap();
let n = bk7.index(Length::nm(587.56));                 // 1.5168
let bare = fresnel_reflectance(1.0, n, 1.0);           // 0.0421 — the textbook 4%

// A coating can only scale that down, and it does so spectrally.
let coated = SurfaceFinish::broadband_ar()
    .reflectance_at(1.0, n, 1.0, Length::nm(550.0));
assert!(coated < bare / 10.0);

// A lamp has a temperature, and Planck decides what colour it is.
let tungsten = Spectrum::blackbody(3200.0);
assert!(tungsten.at(Length::nm(450.0)) < 0.45 * tungsten.at(Length::nm(650.0)));
```

The dimensions carry the coupling between domains too, and check it:

```rust
use dualis_units::{Area, HeatCapacity, Irradiance, Length, Mass, Power, SpecificHeat, Temperature};

let absorbed: Power = Irradiance::mw_per_cm2(50.0)
    * (Length::mm(10.0) * Length::mm(10.0))   // an Area
    * 0.02;                                   // what SurfaceOptics::absorptance returns

let capacity: HeatCapacity = Mass::g(2.0) * SpecificHeat::j_per_kg_k(858.0);
let rise: Temperature = (absorbed * dualis_units::Time::s(1.0)) / capacity;
// 0.58 mK, and every step of that chain landed on the dimension that names it.
```

## Watts, photons, and the chain that closes

A spectrum is a shape until it is integrated against something. `SpectralPower`
carries a shape and a total wattage, and answers two different questions from the
same distribution: `through()` gives watts and `photon_rate()` gives photons per
second. They are not proportional — a photon at 450 nm carries 1.44 times the energy
of one at 650 nm, so a milliwatt of blue is *fewer* photons than a milliwatt of red,
and every silicon detector responds to the count. Filtering moves the mean wavelength
as well as the total, which is why scaling an unfiltered photon rate by a power
fraction is wrong by nearly three times for a tungsten lamp.

That number is also the seam between the two domains, and the integration test
follows it end to end: a dichroic under a 5 W lamp absorbs 96 mW, which warms a 25 mm
lens 10.0 K above ambient, which grows its 100 mm mount by 7.10 µm — 81% of the depth
of focus at NA 0.25, and more than all of it at NA 0.30. Radiometry comes from the
optics crate, heat capacity and expansion from the kernel's `Substance`, depth of
focus from diffraction; the dimensions are what let them compose, and the `Exchange`
audit is what proves nothing leaked on the way.

## What is not here

No scene graph, no renderer.

Wave optics is present but analytic only: Airy patterns, encircled energy, the ideal
MTF, Rayleigh and Abbe limits, depth of focus, Strehl. There is no wavefront
propagation and no pupil FFT, so an aberrated PSF cannot be computed — only the
diffraction-limited ceiling a real system is measured against.

Meshes and grids are no longer excluded. They were, on the grounds that adding them
before a second consumer would be guessing at an interface; finite elements and
finite differences need them, so that decision has been reversed deliberately rather
than drifted away from.

And some things stay out for reasons that will not change: general relativity and
quantum field theory are research subjects rather than simulation targets, turbulent
DNS at world scale is a question about supercomputer budgets rather than about API
design, and no single `f64` state vector spans the fifteen decades from a nucleus to
a galaxy — a simulation has to declare which regime it is in.

## Development

Everything CI runs, in the order it runs it:

```sh
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --no-deps
cargo test --locked --workspace
cargo test --locked --workspace --release
cargo build --locked --workspace --target wasm32-unknown-unknown
cargo test  --locked --workspace --target wasm32-wasip1   # needs wasmtime on PATH
```

Two of those exist to enforce claims rather than to catch typos. The test suite runs
on Linux, macOS and Windows because `rng::tests::the_stream_is_pinned` asserts a
hardcoded digest of ten thousand draws — a platform that rounded differently would
fail there rather than quietly rendering a different image. And the suite runs again
under `wasm32-wasip1` in wasmtime, because "and in WebAssembly" is a claim about
results, which compiling for the target does not establish.

`--locked` throughout: a stale `Cargo.lock` should fail the build rather than be
silently updated, so CI compiles what a contributor compiled.

`-D warnings` is passed to clippy and rustdoc rather than set in `RUSTFLAGS`, since
`RUSTFLAGS` reaches dependencies too and would break the build on somebody else's
warning. Clippy runs the rustc lints as well, so our own warnings are still errors.

**There is no declared MSRV yet.** The newest language feature in the workspace is
`let ... else`, which landed in 1.65, so the floor is likely lower than the 1.75 the
CI job pins — but an MSRV that has never been compiled is worse than none, because a
consumer's toolchain check would be enforcing a guess. When that job passes, the
number moves into `[workspace.package] rust-version`.

## Licence

Dual-licensed under either [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE) at your
option — the Rust ecosystem's convention. MIT is short and universally understood;
Apache-2.0 adds an explicit patent grant for anyone who needs one.
