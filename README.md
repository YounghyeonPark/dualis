# dualis

[![CI](https://github.com/YounghyeonPark/dualis/actions/workflows/ci.yml/badge.svg)](https://github.com/YounghyeonPark/dualis/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/dualis.svg)](https://crates.io/crates/dualis)
[![docs.rs](https://docs.rs/dualis/badge.svg)](https://docs.rs/dualis)

Physics for simulated worlds — a kernel that knows nothing about any particular
physics, and five domains built on it that do: **light, heat, motion, sound, and matter
one atom at a time.**

Dimensions live in the type system, so `Length + Time` does not compile. Conservation is
audited rather than assumed, and a `Violation` names what went missing and where. Every
result is reproducible bit for bit — across platforms, across optimisation levels, in
WebAssembly, and under sixteen threads as readily as one.

```sh
cargo add dualis                             # one dependency, all eight crates
```

Or, in a clone of this repository:

```sh
cargo run --release --example melting        # a crystal melting, read off its own structure
cargo run --release --example beam_hot_spot  # a laser on a mirror, and the hot spot a lumped model misses
cargo test --workspace                       # 357 tests, all against closed forms
```

Add `out.svg` to either example and it draws the result. There are five of them; the table
is [further down](#examples).

> **Reading this as an AI agent, or in a hurry?** [AGENTS.md](AGENTS.md) is the whole API on
> one page, and `cargo run --example agents_quickstart` is a runnable version of it —
> including a deliberate 10% energy leak, so you can see what the audit says when a model is
> wrong. Working *on* dualis rather than with it: [CLAUDE.md](CLAUDE.md).

Every claim in this workspace is checked against a closed form or against an independent
computation rather than against an application that might be wrong in the same direction.
Where no closed form exists, the README says so.

There is now one consumer, `dualis-world`, and its first job was not to be a good
application but to use the SDK the way a stranger would. It came back with five places the
API is awkward and one real defect in the acoustic startup, all in
[`crates/dualis-world/FRICTION.md`](crates/dualis-world/FRICTION.md) — twelve of them now,
seven fixed and five argued down in writing, one of those because the kernel already refuses
the mistake it describes. Not
one of the library's own tests could have found them — the ergonomic ones because a test is
written by somebody who already knows the shape, and the defect because nothing here was
checking a *rate*. There are two tests for that rate now.

## The crates

| Crate | |
| --- | --- |
| `dualis-units` | Dimensional analysis. SI quantities and vectors whose dimension lives in the type, so `Length + Time` does not compile |
| `dualis-core` | The kernel: conservation audits, fixed-step integrators, fields, shared boundaries, multi-domain scheduling, deterministic sampling, closed-form rigid motion |
| `dualis-optics` | Light: spectral radiometry, surface optics, dispersion, ray geometry, diffraction |
| `dualis-thermal` | Heat: lumped masses, explicit conduction, radiative and convective loss |
| `dualis-mechanics` | Motion under force: N-body, Barnes-Hut, penalty contact, rigid rotation |
| `dualis-acoustic` | Sound: the wave equation on a staggered grid, impedance boundaries |
| `dualis-molecular` | Matter atom by atom: Lennard-Jones fluids in periodic boxes, cell lists, a Langevin bath, radial distributions |
| `dualis` | A facade over the other seven, and where the cross-domain integration tests live |
| `dualis-world` | The first consumer, and not published. Worlds described as data: built, coupled over the bus, run and drawn, with ten scenes across all five domains that CI runs. It exists to use the SDK from outside and write down where that is awkward |

```text
dualis-units       no dependencies but glam and serde
dualis-core        depends on units
dualis-optics      depends on core     ─┐
dualis-thermal     depends on core      ├─  none of these knows about
dualis-mechanics   depends on core      │   any of the others
dualis-acoustic    depends on core      │
dualis-molecular   depends on core     ─┘
dualis             depends on all of them
dualis-world       depends on the facade, and nothing depends on it
```

**The kernel must never depend on a domain.** If a new physics needs the kernel
changed, the kernel was wrong — that rule is what makes "add sound, add fluids" a
matter of writing a crate rather than editing this one.

Five domains are now the proof rather than an assertion. Optics publishes absorbed
light as heat and thermal consumes it; mechanics publishes a dashpot's dissipation on
the same channel and the same thermal domain consumes that too, with nothing changed
on either side; acoustics publishes what an absorbing duct end radiates onto the same
channel again; and molecular dynamics arrived last without asking for anything. None of
the five names another and none of them needed the kernel changed.

One thing did need the kernel changed, and it is worth being precise about why that is
not a violation of the rule. The rule is that a *domain* must never force a kernel edit.
What forced this one was the coupling mechanism itself being under-specified from the
start: a bus carrying one number per channel could not say *where* on a surface a
quantity crossed, so no domain could ask for that and none of them ever did. Adding
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
obvious choice is the bad one. Measured from absolute zero the aluminium bar in these
tests holds 1.42 kJ, so a millijoule arriving is a change in the seventh significant figure,
and differencing two such numbers leaves a rounding floor of a few times 10⁻¹² J whatever
the transfer was — worse as the grid is refined, because there are more absolute
temperatures to add up: 1.6×10⁻¹² J at 41 cells against 7.3×10⁻¹² J at 161. Measured from the initial temperature, the number being summed *is* the change.

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

## Examples

```
cargo run --example beam_hot_spot            # numbers, checked
cargo run --example beam_hot_spot out.svg    # and a picture
```

| | |
| --- | --- |
| `beam_hot_spot` | A 100 W laser on a mirror. Optics and heat meeting over a shared boundary, and the peak temperature a lumped model cannot see |
| `airy_pattern` | What a *perfect* lens does to a point: the Airy pattern, its encircled energy, and the MTF that follows from both |
| `detector_snr` | Why read noise is the first number on a datasheet, with the Poisson statistics sampled rather than asserted |
| `room_modes` | Why a small room booms at one note, and why that note is not a note. Four mode shapes and a corner trace |
| `melting` | A Lennard-Jones crystal melting, read off its own radial distribution rather than declared |

Two more are run by CI without being in the table, because they are checks rather than
showcases: `agents_quickstart`, the runnable form of [AGENTS.md](AGENTS.md), and
`readme_check`, which re-runs this file's own code so the snippets above cannot rot.

Each one prints its numbers and **asserts** them, so CI runs all of them on every commit.
An example is a claim that the library works, which makes a quietly broken one worse than
no example at all — every value printed has been checked against a closed form or against
a calculation that did not go through the same code. Give a path and it also writes an SVG;
give none and it just checks. Nothing generated is committed.

Plotting has no dependency. SVG is text, so it is a `format!` and a file write — no
encoder, no fonts, and it opens by double-click. `examples/common/svg.rs` is about three
hundred and fifty lines and is the right size for this job; when it stops being, the answer is
`dualis-world` rather than a bigger version of it. That turned out to be a prediction rather
than a plan: `dualis-world` had to write its own renderer, because this one lives under
`examples/` where no other crate can reach it — `FRICTION.md` finding 4.

The examples also exist to keep the library honest in a way tests cannot. `ScalarField` was
written as the interface a visualiser would read a simulation through and then sat with no
implementor at all, which meant "is it the right interface" was a guess. There are two now —
a one-dimensional bar governed by diffusion and a two-dimensional room governed by a wave —
and the pair has said more about the interface than either could alone. See the next section.

An actual visualiser has since said something neither could: `ScalarField` is the right
*shape* and is unreachable through the thing you have. A renderer holds `&dyn Domain` and
there is no way to ask it for a `&dyn ScalarField`, so `dualis-world` downcasts to concrete
types instead and knows every domain by name — which is exactly what the interface existed to
avoid. `FRICTION.md` finding 3.

## What implementing the field interface found

The kind of thing that only surfaces when something uses a trait, and more of it once there
were two implementors that disagreed.

**`ScalarField::at` takes a time, and a marched domain has only "now".** A closed-form field
like `Motion` answers for any instant; `Bar1D` holds one state and no history, so the
argument is ignored. That would have made the default `rate` — a difference across two
times — read zero, which is worse than unavailable, since the bar is visibly heating. The
fix is not a workaround: a diffusive field's time derivative *is* `α∇²T`, so the governing
equation supplies from the present state exactly what the finite difference wanted history
for. And it is not approximately right. The explicit update is `T += α·dt·∇²T` on the same
stencil, so the rate the field reports is bit-for-bit the step the domain is about to take.

**A field's gradient cannot always be the derivative of its own values.** `at` interpolates
linearly between cell centres, so its exact second derivative is zero between nodes and
infinite at them — a Laplacian read off the interpolant is useless, and must come from the
cell stencil instead. The cost is that `gradient` is the derivative the *scheme* uses rather
than of `at` exactly. That is a property of sampling a discrete field, not a rough edge that
could be polished out, and it is documented where a caller meets it.

**`rate` is not always about the next step.** `Bar1D` marches forward, so the rate it reports
is exactly the change the next step will make. `Room` is a leapfrog and stores its velocities
half a step behind its pressures, so the velocity for the next update does not exist yet —
what it can report exactly is `(pⁿ − pⁿ⁻¹)/h`, the step just taken. Being half a step behind
is the better trade rather than a compromise: a centred difference at the midpoint is second
order where a forward one at the present instant is first.

**And a claim of mine that was simply wrong.** I wrote that the mirrored stencil gives a zero
gradient at an insulated wall. For `Bar1D` it does not: its grid is cell-centred, so the first
sample sits half a cell inside the wall where the temperature really is still changing.
`Room`'s grid is node-centred, a node sits *on* the wall, and there the gradient is zero to
the last bit. Same physics, different sampling — and only the second implementor could have
shown that the first one's behaviour was about the grid rather than about the boundary.

Two smaller things a caller will meet. A grid-backed field interpolates in `at` but snaps to
the nearest node for its derivatives, because the derivatives are the scheme's own stencils
and the scheme has values only at nodes — so combining them at an arbitrary point divides a
number computed here by one computed slightly over there. And a signed field needs a
different colour map from a positive one: a temperature climbs from a floor, a pressure swings
about zero, and a ramp built for the first draws a room at rest as though half of it were
cold.

## A boundary defect, and the thing that found it

Both `Tube` and `Room` put their pressure samples *on* the walls rather than half a cell
inside, so a wall sample owns half a cell. Both divided its divergence by the whole `dx`. The
walls came out twice as heavy as they are, and every mode read low: 1.4% on an 89-cell room,
5.4% on a coarse one.

A percent and a half looks like discretisation error and invites a loose tolerance. What
showed it was not: refining the grid **halved** the error instead of quartering it. A scheme
second order in the interior converging at first order overall means the boundary is first
order, and that is a wrong condition rather than a coarse one. The order of convergence found
it; no single run could have.

Two independent confirmations that the factor of two was the right one. The mirrored
five-point stencil — the standard second-order treatment of a zero-gradient boundary — gives
`2(p₁ − p₀)/dx²` at a wall, and the scheme was producing exactly half of that; with the fix
they agree, which is now asserted directly by stepping once from rest and comparing against
`ScalarField::laplacian` at every node. And the conservation audit still holds at 1e-15,
because the energy carries the same half weights the update does. Those two had to move
together: changing the update alone breaks the audit, and changing the audit alone hides the
bug.

The fix cost something real, which is the interesting part. `Tube`'s absorbing ends had been
built against the unweighted update, and at the full CFL limit the corrected boundary drains
a half-width cell with a factor of exactly **−1** — it inverts the wave instead of swallowing
it. Not a divergence: a perfectly stable run that quietly reflects. So `max_stable_dt` now
reports the impedance's own limit as well as the wave's, `Z·dx/2ρc²`, which halves the step
for a matched end and leaves a closed one alone.

### The same crate, the same shape, a second time

There was another one, in the boundary in *time*, and it was still there after all of the
above. `released_from` set the velocity to zero at `t = 0`, but a staggered leapfrog carries
velocity at `t = −h/2`, so the first velocity update travelled a whole step where the initial
condition entitled it to half. `O(h)`, permanent, and `h` follows `dx` — first order again,
from a scheme whose interior and whose walls were now both second order.

Found the same way, by the rate. Not by anything in the library: it took `dualis-world`, the
first consumer, checking a released mode against `|cos(2πft)|`. The worst departure over 20 ms
went from 0.0528 to 0.00238 at 31 cells once the first update took half a step.

Two things fell out of it that are worth more than the fix.

**A test had turned the bug into the specification.** One step from rest was asserted to move
the pressure by `h²c²∇²p` — but from rest `ṗ(0) = 0`, so Taylor gives `½h²c²∇²p`. The test was
missing the half because the scheme was, and `Tube` carried the identical pair. Both were
written by reading the implementation. A test written from the closed form cannot do that,
which is why the conventions say to write it that way.

**And the old startup conserved energy exactly, while the correct one does not.** With `v = 0`
read as the half-step value, `Σ∇·(p∇p) = 0` at a rigid wall makes the first step's energy
change cancel to the last bit. Starting correctly breaks that cancellation by `O(h²)` — 0.42%
at 31 cells, quartering on refinement, and only at the first step; from there the invariant
holds to 1e-15. So the old code had bought exact bookkeeping by making the scheme first order.
That is the trap this repository already had written down — *the functional and the update
were consistent with each other and both wrong* — turning up a second time in the same crate,
and the audit that was supposed to catch it was the thing keeping it in place. The energy is
now measured against the released state as its datum, with the difference reported by
`Room::startup_adjustment` and bounded so a real first-step bug cannot hide in it.

## The domain whose answers are distributions

Every domain but one is checked against a value: a Fresnel coefficient, a mode frequency, a
temperature rise. `dualis-molecular` mostly cannot be. A hundred atoms in a box are chaotic by
construction — perturb one coordinate in its last bit and the trajectories separate completely
inside a few hundred steps — so "close enough" is not available as a notion, and the physics
lives in what the trajectory averages to.

There turns out to be plenty of that, and all of it exact:

| Claim | Exact form |
| --- | --- |
| Equipartition | `⟨KE⟩ = (3N − 3) k_BT / 2` |
| Ideal gas, dilute | `PV = N k_BT`, departing as `1 + B₂(T)ρ` |
| Virial pressure, any density | `PV = N k_BT + ⟨Σ f·r⟩/3` |
| Lennard-Jones well | `−ε` exactly, at `2^(1/6)σ` |
| Momentum | To the last bit, by Newton's third law |
| Energy, unthermostatted | Bounded rather than drifting, because Verlet is symplectic |

The `− 3` in the first row is not a rounding. Momentum is conserved and the drift was removed,
so the centre of mass is frozen and three degrees of freedom are gone; counting all `3N` makes
every reported temperature low by `1/N`, which is a percent at a hundred atoms and looks
exactly like statistics.

Two things this domain leans on harder than any other. `Rng::for_index` keys the Langevin noise
on `(seed, step, particle)`, which is the *only* route back to a reproducible run in a chaotic
system — there is no "close enough" to fall back on. And the symplectic argument behind
`velocity_verlet`, which the kernel proves on a harmonic oscillator and which is relied on here
for a many-body potential with a truncated force.

The `melting` example is where this pays off. Three state points are equilibrated and the
radial distribution function decides which phase each one is, rather than anyone declaring it:
`g(r)` peaks at 5.1 for the crystal, 2.9 for the liquid and 1.4 for the gas, while the *number*
of first neighbours barely moves between solid and liquid — 12.4 against 12.8. Melting costs the
order and keeps the packing, which is why a liquid is nearly as dense as its solid and a gas is
a thousand times thinner.

The crystal panel is checked exactly before anything harder is attempted: an fcc lattice's
neighbour shells are at `1 : √2 : √3 : 2` holding 12, 6, 24 and 12 atoms, and the histogram
reproduces all four to 1e-12. Combinatorics rather than measurement.

And one thing the data corrected. Long-range order should be checked at the *third* shell, not
the fourth. The fourth was the obvious choice and reads 1.045 for the crystal against 0.999 for
the liquid — no discrimination at all — because thermal broadening washes out a thin shell at a
large radius long before a fat one nearer in. The third holds 24 atoms, twice as many as the
fourth, and separates them cleanly at 2.81 against 1.27.

One test in the crate is worth describing because of how it was wrong first. The ideal-gas check
asserted that doubling the density doubles the departure, ran on one seed, and passed. Across
four seeds the ratio came out 1.35, 1.61, 2.34 and 2.92 — averaging to 2.06, but landing
anywhere. A hundred atoms is a small sample and two thousand correlated snapshots are far fewer
independent ones than they look. The test now averages over seeds, and the fix for a noisy
statistical test is more samples rather than a wider tolerance.

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
Asking for encircled energy at the wrong one of them gives 59% instead of 84%, which
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

No scene graph, and no renderer *in the library*. `dualis-world` has both — a JSON scene
format and an SVG filmstrip — which is why it is a separate crate and is not published.

Rigid bodies are spheres where they collide. `RigidBody` rotates with an applied torque
and an arbitrary inertia tensor, but `Sphere` and `Rolling` are the only things that
touch anything, and a sphere is chosen because its inertia is the same about every axis
— so a contact does not have to carry the orientation through. Boxes hitting boxes is a
different module.

Acoustics is linear and up to two dimensions: a tube or a room, not a hall. No 3D grid,
no scattering geometry, and no nonlinearity, so nothing here shocks up or distorts.

Molecular dynamics is monatomic, so the crate's name is aspirational by one step: no bonds,
angles or torsions, and therefore no molecules. And no electrostatics, which for a charged
system is the hard part rather than a missing feature — Coulomb falls off as `1/r` and cannot
be cut off at all, so it needs Ewald summation or a particle-mesh method, and that is a larger
piece of work than everything currently in the crate. No constraints, no barostat, no
free-energy machinery.

Both acoustic domains had a boundary defect until recently, and how it was found is worth
more than the fix. See [the section above](#a-boundary-defect-and-the-thing-that-found-it).

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
back most of that accuracy at close angles — a factor of 6.5 at `θ = 0.3` — but nothing
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

[`CONTRIBUTING.md`](CONTRIBUTING.md) has the conventions that are unusual enough to be worth
stating — check against a closed form and never against another implementation, earn every
tolerance, judge a residual against a scale, and keep the kernel ignorant of every domain.
[`CHANGELOG.md`](CHANGELOG.md) records what was found as well as what was added.

Everything CI runs, in the order it runs it:

```sh
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --no-deps
cargo test --locked --workspace
cargo test --locked --workspace --release
cargo build --locked --workspace --exclude dualis-world --target wasm32-unknown-unknown
cargo test  --locked --workspace --exclude dualis-world --target wasm32-wasip1   # wasmtime
cargo +1.78 build --locked --workspace --exclude dualis-world
cargo deny check
```

The three `--exclude`s are deliberate: WebAssembly, determinism and the 1.78 floor are
promises the *library* makes to the people who depend on it, and `dualis-world` is an
unpublished application with no dependents. CI additionally *runs* every example and every
scene rather than only compiling them — see [Examples](#examples) and
`crates/dualis-world/scenes/`.

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
option — the Rust ecosystem's convention, and not indecision. Each covers a gap the other
has. Apache-2.0 carries an explicit patent grant and MIT does not, which is the first thing a
corporate legal review asks about. MIT is compatible with GPLv2 and Apache-2.0 is not, which
matters here because the scientific-computing world has plenty of GPLv2 code. Offering both
lets a consumer take whichever they need, so `OR` is strictly *less* restrictive than either
alone.

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion
in the work by you, as defined in the Apache-2.0 licence, shall be dual-licensed as above,
without any additional terms or conditions.

Every dependency is permissive too, and that is checked rather than remembered: `deny.toml`
holds an allow-list and CI fails on anything outside it. Twelve external crates at the time of
writing, of which three reach a *published* artifact — `glam`, `serde` and `serde_core`, all under
the same `MIT OR Apache-2.0`. `dualis-world` also links `serde_json` and its three
transitive crates, but it is not published. The rest are compile-time or test-only.
