# dualis

[![CI](https://github.com/YounghyeonPark/dualis/actions/workflows/ci.yml/badge.svg)](https://github.com/YounghyeonPark/dualis/actions/workflows/ci.yml)

Physics for simulated worlds — a kernel that knows nothing about any particular
physics, and domains built on it that do.

Extracted from AiryTrace, a non-sequential ray tracer, now archived. Nothing consumes
this workspace yet — so every claim in it is checked against a closed form or an
independent computation rather than against an application that might be wrong in the
same direction.

## The crates

| Crate | |
| --- | --- |
| `dualis-units` | Dimensional analysis. SI quantities and vectors whose dimension lives in the type, so `Length + Time` does not compile |
| `dualis-core` | The kernel: conservation audits, fixed-step integrators, fields, shared boundaries, multi-domain scheduling, deterministic sampling, closed-form rigid motion |
| `dualis-optics` | Light: spectral radiometry, surface optics, dispersion, ray geometry, diffraction |
| `dualis-thermal` | Heat: lumped masses, explicit conduction, radiative and convective loss |
| `dualis-mechanics` | Motion under force: N-body, Barnes-Hut, penalty contact, rigid rotation |
| `dualis-acoustic` | Sound: the wave equation on a staggered grid, impedance boundaries |
| `dualis` | A facade over the six, and where the cross-domain integration tests live |

```text
dualis-units       no dependencies but glam and serde
dualis-core        depends on units
dualis-optics      depends on core     ─┐
dualis-thermal     depends on core      ├─  none of these knows about
dualis-mechanics   depends on core      │   any of the others
dualis-acoustic    depends on core     ─┘
dualis             depends on all of them
```

**The kernel must never depend on a domain.** If a new physics needs the kernel
changed, the kernel was wrong — that rule is what makes "add sound, add fluids" a
matter of writing a crate rather than editing this one.

Four domains are now the proof rather than an assertion. Optics publishes absorbed
light as heat and thermal consumes it; mechanics publishes a dashpot's dissipation on
the same channel and the same thermal domain consumes that too, with nothing changed
on either side; acoustics publishes what an absorbing duct end radiates onto the same
channel again. None of the four names another and none of them needed the kernel
changed.

One thing did need the kernel changed, and it is worth being precise about why that is
not a violation of the rule. The rule is that a *domain* must never force a kernel edit.
What forced this one was the coupling mechanism itself being under-specified from the
start: a bus carrying one number per channel could not say *where* on a surface a
quantity crossed, so no domain could ask for that and none of the four ever did. Adding
`Interface` and `Flux` did not teach the kernel any physics — a discretised boundary is
not optics or heat — and no existing domain had to change to keep working.

They also brought the audit three conserved quantities instead of one, and the three
hold to wildly different tolerances for structural reasons rather than through
differences in effort:

| Quantity | Holds to | Because |
| --- | --- | --- |
| Linear momentum (`NBody`) | 1e-13 | Exact by construction: equal and opposite forces cancel bit for bit |
| Linear momentum (`TreeNBody`) | θ-dependent | Each body sees its own approximation of the rest, so nothing cancels |
| Angular momentum (`RigidBody`) | 1e-9 | Nothing cancels here either; it is only as good as RK4 plus a quaternion renormalisation |
| Energy across a coupling | 1e-9 | Both sides are closed-form evaluations, nothing is integrated |
| Energy through a contact | 2e-2 | A penalty contact is a non-smooth potential, so the symplectic bound does not apply |

Auditing a vector component by component also makes the *smallest* component the
binding constraint, since the absolute error is set by the whole vector while the scale
it is judged against is only that component. That is a property of the kernel's
per-quantity audit, and it is written down where it will be met.

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

The parallel half of that is executed rather than argued. `TreeNBody::with_threads`
changes how many threads evaluate forces, and
`tree::tests::parallel_and_sequential_agree_bit_for_bit` asserts the answers are
identical across one, two, four, eight and sixteen of them — through a whole
integration, not just one evaluation. It holds because each thread owns a disjoint
range of the output: there is no reduction, so there is no summation order to vary.

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

## Where two domains meet

Sharing a clock is not enough; they also have to share a *place*. `Exchange::publish`
carries an amount and nothing else, so "a 1 mm beam on a 20 mm mirror" and "one watt
spread over the whole mirror" are the same message — and a coating fails at its hot
spot, not at its average.

An `Interface` is a boundary cut into faces that both sides address, and a `Flux` is a
quantity spread over them. Both sides share one discretisation on purpose: interpolating
between two meshes is where energy quietly goes missing, so a face-count disagreement is
refused rather than papered over, and a caller who genuinely needs to cross grids says so
with `Flux::resample`, which conserves the total by construction. Spatial channels are
audited **face by face** — a redistribution that keeps the total but moves it to the wrong
end of a mirror is exactly the bug a total-only check cannot see.

Three things fell out of building it that are worth knowing.

An enthalpy's reference point is arbitrary, so it should be chosen for precision, and the
obvious choice is the bad one. Measured from absolute zero a warm aluminium bar holds
228 kJ, so a millijoule arriving is a change in the ninth significant figure and
differencing two such numbers leaves a rounding floor of about 10⁻¹¹ J whatever the
transfer was — worse as the grid is refined, because there are more absolute temperatures
to add up. Measured from the initial temperature, the number being summed *is* the change.

An insulated bar under a continuous beam does not flatten out. It settles into a fixed
shape and rides upward on a mean that climbs forever, so the hot spot is permanent and the
lumped model's error in kelvin never shrinks; what shrinks is the *fraction*. The
intuition that says it evens out eventually is about a bar heated once.

And how often the domains talk is itself an error source — here the largest one, and the
only one no conservation check can catch. A 10 ms coupling window on a bar whose own step
is 0.44 ms delivers exactly the right joules and reads the peak temperature 12% low,
because the bar spends most of each window relaxing with nothing arriving. Nothing is
lost, nothing is created, both domains stay inside their stability limits, and the answer
is still wrong. The only thing that finds it is a solution that never went through the
coupling: `tests/beam_heats_where_it_lands.rs` computes one by quadrature from the
steady-state energy balance, checks *that* against the exact `ṫL²/12α` for a point
source, and then watches the coupled answer converge on it as the window shrinks.

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

## Numerical against analytic, on purpose

Two ways of computing the same optics, kept side by side so each checks the other.
`diffraction` answers what a perfect system does, in closed form from Bessel
functions. `wavefront` answers what an imperfect one does, by transforming a pupil.
Set the aberrations to zero and the second must reproduce the first — and it does, in
three independent places that share no code:

| | |
| --- | --- |
| Airy profile | The transformed pupil agrees with `[2J₁(v)/v]²` to under 0.5% of the peak, out to 3 λ/D |
| Ideal MTF | The transform of the PSF agrees with `(2/π)(arccos s − s√(1−s²))`, the closed-form autocorrelation of a disc |
| Strehl ratio | Small aberrations follow `exp(−(2πσ)²)`, the Maréchal approximation, whichever mode produced the error |

The comparison is what makes either side trustworthy. A pupil transform has a great
deal of room to be subtly wrong — a sampling factor, a sign, a normalisation — and
none of it shows up as anything but a plausible picture.

It also catches mistakes in the *questions*. The first Airy zero is at 1.22 λ/D and
at 0.61 λ/NA; those are the same zero and the factor of two is the numerical aperture.
Asking for encircled energy at the wrong one of them gives 57% instead of 84%, which
is exactly what happened the first time.

## What the audit tolerance is really measuring

A conservation audit needs a tolerance, and the right one is a property of the
*integrator*, not a concession to sloppiness. Three cases in this workspace, all
different:

- **Momentum under N-body gravity**: 1e-11 and could be tighter. Newton's third law
  is structural here, so the only loss is turning each force into an acceleration and
  back.
- **Energy across an optics-to-thermal coupling**: 1e-9. Both sides are closed-form
  evaluations, and nothing is being integrated.
- **Energy through a penalty contact**: 2e-2. Semi-implicit Euler is symplectic, so
  its energy error is meant to stay bounded — but that guarantee needs a smooth
  potential, and a contact that switches on the moment two things touch is not one.
  Each transition shifts the shadow Hamiltonian the method actually conserves, so the
  true energy takes an `O(dt)` step at every bounce and those accumulate. Resolving
  the contact more finely shrinks each step, which is why `ContactSystem` reports a
  stability limit a hundred times smaller than bare stability needs, but no step size
  recovers the bound.

The third is the interesting one, because a tolerance chosen without knowing that
would either fail on correct code or hide a real leak.

## What is not here

No scene graph, no renderer.

Rigid bodies are spheres where they collide. `RigidBody` rotates with an applied torque
and an arbitrary inertia tensor, but `Sphere` and `Rolling` are the only things that
touch anything, and a sphere is chosen because its inertia is the same about every axis
— so a contact does not have to carry the orientation through. Boxes hitting boxes is a
different module.

Acoustics is linear and up to two dimensions: a tube or a room, not a hall. No 3D grid,
no scattering geometry, and no nonlinearity, so nothing here shocks up or distorts.

**There is no fluid domain, and that is deliberate.** Sound *is* the fluid domain here:
it is what a fluid does when the variations are small enough to linearise, which is
exactly the regime where every answer has a closed form to check against. Full
Navier-Stokes has none, and a solver that could not be validated against anything would
be decoration. See the note below on turbulence.

Gravity comes both ways and the choice is a real one. `NBody` sums every pair: exact,
momentum conserved to the last bit, `O(n²)`, and awkward to parallelise precisely
because the `i < j` pairing that makes it exact has two threads writing to one body.
`TreeNBody` is Barnes-Hut: `O(n log n)`, embarrassingly parallel, and it **gives up
exact momentum** — each body sees its own approximation of the rest, so their mutual
forces no longer cancel. The drift is a knob rather than a defect, closing with the
opening angle and vanishing at `θ = 0`, and the audit tolerance has to be the one the
angle earns. The expansion carries the quadrupole as well as the monopole, which buys
back most of that accuracy at close angles — a factor of six at `θ = 0.3` — but nothing
past it, so `θ` above about 1 is still asking a centre of mass to stand in for a group
that is not far enough away.

Fields propagate between planes by the angular spectrum, but only through free space:
there is nothing to put in the beam's way except an aperture, and the grid's reach is
`NΔ²/λ`, past which the propagator refuses rather than aliasing.

Partial coherence is a transfer function rather than a simulation. The two exact limits
are there, and so is the coherence a source has at a distance, but computing an
arbitrary object's partially coherent image needs the transmission cross-coefficients —
a four-dimensional integral, and not here.

The tree stops at the quadrupole. Higher multipoles would buy another order in the
opening angle, and a proper Fast Multipole Method would change the complexity rather
than the constant.

Meshes and grids are no longer excluded. They were, on the grounds that adding them
before a second consumer would be guessing at an interface; finite elements and
finite differences need them, so that decision has been reversed deliberately rather
than drifted away from.

An `Interface` is one-dimensional: a boundary is a sequence of faces, and `Flux::resample`
remaps between two of them by overlap in cumulative area. That is enough for a mirror
face, a bar's side, a row of pixels, and any boundary whose faces have a natural order —
and it is not enough for a triangulated surface or an arbitrary mesh-to-mesh projection,
where the overlaps are not an interval intersection and there is no order to walk. The
conservation argument generalises; the implementation does not, and the doc comment says
so where a caller will meet it rather than only here.

An `Interface` also carries areas and an order and nothing else — no coordinates, no
normals, no connectivity. A domain that needs to know where a face *is* in space still
has nowhere to put that, so a beam profile has to be handed over in the boundary's own
coordinate rather than computed from geometry. That is the next thing this layer is short
of, and it is deliberately not guessed at before something needs it.

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

**The MSRV is set by the lockfile, not by the code.** CI builds on 1.78, and that
number came out of a failure worth recording: the job was first pinned at 1.75 and
died with `failed to parse lock file` before compiling anything, because `Cargo.lock`
is format version 4 and cargo could not read that until 1.78. The newest language
feature in the workspace is `let ... else` from 1.65, so the source would go lower.

Which floor applies depends on who is asking. A consumer depending on `dualis-optics`
never receives this lockfile, so their constraint is the source and its dependencies.
CI passes `--locked` deliberately, so its constraint is the lockfile format. The
declared `rust-version` follows CI, because it is the stronger of the two and a
declared MSRV should be a promise about what has been compiled.

## Licence

Dual-licensed under either [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE) at your
option — the Rust ecosystem's convention. MIT is short and universally understood;
Apache-2.0 adds an explicit patent grant for anyone who needs one.
