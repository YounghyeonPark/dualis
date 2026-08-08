# What the first consumer found

`dualis-world` exists to use the SDK from outside and write down where that is awkward. A
library with no consumers is a library whose ergonomics nobody has measured, and none of the
345 tests inside the workspace can answer this question about themselves — they are written by
someone who already knows the shape.

Everything below was hit while building the smallest thing that loads a scene, runs it, couples
two domains over a plain channel and two more over a shared boundary, and draws the result. None of it is a bug in the physics except finding 6, which is — and which no test inside the
library could have found, because none of them was checking a rate.

**Six of the ten are fixed**, and four are recorded rather than actioned — three of
those because the kernel already refuses the mistake they describe. The entries are
kept rather than deleted, because what the API used to be is the argument for what it is — and because the next consumer should be able
to see that the answer to "this is awkward" was to change the library rather than to work
around it. Each fixed entry says what was done.

---

## 1. A domain cannot be built behind a `dyn`

`Simulation::with` takes `impl Domain + 'static` by value, and there is no
`impl Domain for Box<dyn Domain>`. Internally the simulation already stores
`Vec<Box<dyn Domain>>`, so the boxing happens either way — it just cannot happen on the
caller's side.

The consequence is that a data-driven builder must be a `match` with one arm per domain type,
and each arm has to call `with` separately:

```rust
sim = match spec {
    DomainSpec::Room { .. } => sim.with(AcousticRoom::of_air(..)),
    DomainSpec::Bar  { .. } => sim.with(Bar1D::new(..)),
};
```

That works and is what `World::build` does. What it forecloses is a *registry*: a third party
cannot add a domain type to the scene format without editing this match. For a physics SDK
whose central claim is that domains are pluggable, the plug is only available at compile time.

**Fixed.** The kernel gained `Simulation::with_boxed(Box<dyn Domain>)` and an
`impl Domain for Box<dyn Domain>` that delegates every method. `DomainSpec::build` now returns
a box and `World::build` is a three-line loop. The `match` over domain types still exists, but
it is confined to one function and is the *scene format's* business rather than the kernel's —
an out-of-tree domain can be boxed and added without this crate knowing.

## 2. Domain names are `&'static str`, so they cannot come from data

Every constructor takes `name: &'static str`, and `Domain::name` returns one. A name read out
of a JSON file is a `String`. `World::build` therefore calls `Box::leak`.

The leak is bounded by the number of domains in a scene, so it is survivable rather than
dangerous. But it is the API stating that names are compile-time things, and for an
application they are exactly the opposite: they are what the user typed.

This is the friction that felt worst in practice, because it is unavoidable and it appears at
the very first thing an application does.

**Fixed, and it cost less than expected.** `Domain::name` returns `&str`, every domain stores
a `String`, and every constructor takes `impl Into<String>`. Because `&str: Into<String>`,
**not one existing call site changed** — `Bar1D::new("bar", ..)` still compiles. The only
breakage in 349 tests was five comparisons against `Report::substeps`, which had to become
owned for the same reason. `Interface` followed, and `Exchange`'s spatial map is keyed by an
owned interface name now.

The `Copy` that was lost was never in a hot path: a name is read to report a violation and to
look a domain up, a handful of times per step.

## 3. Reading state back needs the concrete type

`Simulation::domain_as::<T>` needs `T`, so the renderer knows every domain type just as the
builder does. `ScalarField` is exactly the abstraction that would avoid this — sample a field,
draw it, never ask what it is — but there is no way to get a `&dyn ScalarField` from a
`&dyn Domain`.

The result is a second `match` over the same enum, in `World::capture`, for no reason other
than downcasting.

**Fixed.** `Domain::as_field` returns `Option<&dyn ScalarField>` and defaults to `None`;
`Bar1D` and `Room` implement it in one line each; `Simulation::field(name)` returns one.
`World::capture` no longer mentions `Room` or `Bar1D` at all — it asks each domain for a field
and samples it. That is what `ScalarField` was written for.

One thing the fix does not give away: a `ScalarField` is a function of position and does not
know where it stops, so the *extent* to sample over still comes from the caller. That is the
right division — a field that knew its own bounds would be a mesh — and the scene has the
bounds already.

## 4. The examples' plotting is not reachable

`crates/dualis/examples/common/svg.rs` is about three hundred and fifty lines of dependency-
free SVG plotting, and it lives under `examples/`, so no other crate can use it. This crate
has its own smaller renderer that overlaps with it substantially.

Not obviously wrong — the examples are meant to be self-contained, and a plotting API is a
commitment. But it means the first thing a consumer wants to do after running a simulation is
something the workspace already solved and cannot share.

**Not fixed, and the only one of the six that is not.** Sharing it means either a `dualis-plot`
crate or a feature-gated module in the facade, and either way it is a public API for drawing
that would have to be supported, versioned and documented — for a workspace whose stated
scope excludes rendering. Two overlapping private renderers is the cheaper mistake for now.

Revisit when there is a second consumer. One application writing its own hundred lines of SVG
is not evidence; two would be.

## 5. `Room` is not in the prelude

`dualis::prelude` re-exports `Tube` but not `Room`, though they are the two headline types of
the same crate. Reached through `dualis::acoustic::Room` instead.

**Fixed.** One line. It was an oversight, as suspected.

## 6. `Room` has a first-order startup error — a real defect

This one is physics, not ergonomics, and it was found by the app checking itself against a
closed form rather than by any test in the library.

A room released in its `(1, 1)` mode should follow `|cos(2 pi f t)|` at every point. It does,
but the gap converges at **first** order against grid resolution, where the scheme's interior
is second:

```text
  31 cells   0.0528       241 cells  0.0076
  61 cells   0.0265       481 cells  0.0039
 121 cells   0.0151
```

Halving on refinement, not quartering. The cause looks like the leapfrog's startup:
`Room::released_from` sets the velocity array to zero at `t = 0`, but a staggered scheme
carries velocity at half steps, so what is wanted is `v(-dt/2)`. For a mode released from rest
that is `-sin(pi f dt)`, not zero — an `O(dt)` error, and `dt` follows `dx` through the CFL
condition, giving exactly the first order observed.

This is the same shape as the wall-weighting defect the workspace already found and fixed: a
second-order interior dragged to first order by how the boundary — here the boundary in
*time* — is handled. It was found the same way too, by the rate rather than the size.

**Fixed, and `Tube` had it too.** The first velocity update now travels half a step; every
one after it travels a whole one. Second order, and the error at 31 cells fell by a factor of
22 — 0.0528 to 0.00238. `tests/scene.rs` and a new pair in `dualis-acoustic` pin the rate.

Two things the fix turned up that were not visible from the outside:

- **A test had turned the bug into the specification.** `one_step_from_rest_is_the_laplacian_the_field_reports`
  asserted that one step from rest moves the pressure by `h²c²∇²p`. From rest `ṗ(0) = 0`, so
  Taylor gives `½h²c²∇²p` — the test was missing the half, and it passed because the scheme
  was missing it too. `Tube` had the matching test with the matching error. Both were written
  by reading the implementation, which is the failure mode a test written from the closed form
  does not have.
- **The old startup conserved energy *exactly*, and the fix does not.** Not a regression: with
  `v = 0` treated as the half-step value, `Σ∇·(p∇p) = 0` at a rigid wall makes the first step's
  energy change cancel to the last bit. Starting correctly breaks that cancellation by
  `−h²Σ(∇p)²/8ρ` — 0.42% of the total at 31 cells, quartering on refinement, and *only at the
  first step*; from there the invariant holds to 1e-15.

  So the old code bought exact bookkeeping by making the scheme first order. That is the
  workspace's own documented trap — "the energy functional and the update were consistent with
  each other and both wrong" — appearing a second time in the same crate.

  The energy is now reported against the released state as its datum, with the one-off
  difference kept in `Room::startup_adjustment` where it can be asked for, and bounded at 25%
  so a real first-step bug cannot hide in it.

---

## 7. The name change had missed `Bar1D::exposing`

Found by going looking, after the coupling scene needed a boundary. `Domain::name` and every
constructor took `impl Into<String>` after finding 2, but `exposing(boundary: &'static str, ..)`
did not — the sweep had matched on the parameter being called `name`. A boundary name is data
in exactly the same way and for exactly the same reason: two domains agree on it, and what they
agree on can come from a file.

**Fixed.** One signature. Worth its own entry because it is what an incomplete refactor looks
like from outside: the API is *mostly* consistent, and the one place it is not is the place
nobody had reached yet.

The crate-level documentation in `dualis-core` was also still teaching
`fn name(&self) -> &'static str` in both of its worked examples. They compiled, so nothing
failed; they were simply showing the reader the idiom that had just been removed.

## 8. Nothing checks a schedule against the domains until the first step

A scene picks its schedule by name. `staggered` with a half-second frame is thirty-eight times
the bar's explicit-diffusion limit, and the run is refused — correctly, by name, with the
limit and the value, which is the whole argument for this library and it works.

But it is refused *when the step is taken*, not when the scene is built. `Domain::max_stable_dt`
is public, so an application can ask every domain what it can survive and refuse at build time
where the message can name the file and the line. This one does not yet.

Not a library defect. A note about where the natural seam is, and the sort of thing only
somebody loading scenes from disk would think to want.

## 9. A spatial coupling makes both sides state the discretisation, twice

`Bar1D::exposing(name, face_area)` builds its own `Interface` with one face per cell. A
publisher has to build a matching one, and there is no way to derive it from the bar: by the
time a `Box<dyn Domain>` exists, its interface is behind the trait, and a per-spec builder
cannot see another spec's product anyway.

So a scene says the face count twice — once as the bar's `cells` and once as the beam's
`faces` — and can say it inconsistently.

**Not fixed, and the kernel is the reason it does not need to be.** `publish_on` refuses a
flux whose face count differs from the interface's and reports both numbers. That is the right
place for the check: silently padding or truncating would put energy on the wrong part of the
boundary while keeping the total exactly right, which is the one failure a conservation audit
cannot see and the whole reason the spatial channel exists. `a_boundary_the_two_sides_cut_differently_is_refused`
asserts it.

What would remove the duplication is `exposing` taking an `Interface` rather than building
one — then a scene constructs a single boundary and hands clones to both sides. Worth doing
when a second spatial consumer exists; with one, the duplication is two integers in a file
and the kernel already refuses the mistake.

## 10. A sampled field is not the state, and averaging it is not averaging the state

The renderer samples `ScalarField` at evenly spaced points including both ends. `Bar1D`'s grid
is cell-centred, so the two end samples sit half a cell outside the outermost cell centres.
Averaging the samples therefore comes out about `1/2n` low against averaging the cells — 1.2%
at 41 cells.

Found by an assertion failing, not by reasoning: a test checked that the bar held every joule
the beam paid, computed the mean from the render panel, and missed by 1.2%.

**Not a defect anywhere.** `ScalarField` is a function of position and is behaving exactly as
documented; the renderer is sampling it exactly as it should. But an application that reported
a mean temperature from its own render buffer would be wrong by that much with nothing to tell
it, and the two numbers look interchangeable right up until they are compared. The test now
reads the total from the domain and the shape from the panel, and pins the gap between them so
it stays understood rather than rediscovered.

---

## What this says about the exercise

Ten findings from about six hundred lines of application code: six ergonomic, one a real defect
in the physics, and three notes about where a check belongs or why one is already in the right
place. Six are fixed.

Findings 1, 2, 3 and 7 were the same shape. **The API was comfortable when the set of domains
was known at compile time and awkward the moment it was not** — and that was never a decision
anybody made. It is the shape that falls out of writing a library with no consumer, where
`&'static str` costs nothing because every name is a literal in a test.

It has been made deliberately now, in the other direction, and the cost was a tenth of the
argument for keeping it: **no existing call site changed**, five test comparisons did, and the
application lost its leak, both of its downcast matches and about forty lines.

Finding 6 is the one that matters most, and the one nothing inside could have produced. A
first-order startup error survived a second-order interior, a second-order wall fix, exact
energy conservation and 345 passing tests — two of which had turned the bug into the
specification by asserting what the implementation did. It took an outside program comparing a
released mode against `|cos(2 pi f t)|` at four grid resolutions. Nothing in the library was
checking a *rate*.

That is the case for building a consumer early, and it is stronger than the ergonomic half.
None of this was visible from inside.

## What this report does not cover

**Three of the five domains.** Optics, mechanics and molecular have no `DomainSpec` variant, so
their constructors have not been driven from data at all. `Fluid::lattice` alone takes five
arguments of four different kinds; whether that survives contact with a scene file is untested.
