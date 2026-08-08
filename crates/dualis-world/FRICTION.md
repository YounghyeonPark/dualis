# What the first consumer found

`dualis-world` exists to use the SDK from outside and write down where that is awkward. A
library with no consumers is a library whose ergonomics nobody has measured, and none of the
345 tests inside the workspace can answer this question about themselves — they are written by
someone who already knows the shape.

Everything below was hit while building the smallest thing that loads a scene, runs it and
draws it. None of it is a bug in the physics except finding 6, which is.

**Five of the six are now fixed.** The entries are kept rather than deleted, because what the
API used to be is the argument for what it is — and because the next consumer should be able
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

**Not fixed here.** A half-step kick at release would change the numbers in every acoustic
test, and `Tube` very likely has the same startup, so it deserves its own pass rather than
being smuggled in with an application. `tests/scene.rs` pins the rate at first order, so
fixing it will fail that test — which is correct, since the table above would then be wrong.

---

## What this says about the exercise

Five ergonomic frictions and one real defect, from about three hundred lines of application
code. Findings 1, 2 and 3 were the same shape: **the API was comfortable when the set of
domains was known at compile time and awkward the moment it was not.** That was never a
decision anybody made — it was the shape that falls out of writing a library with no consumer,
where `&'static str` is free because every name is a literal in a test.

It has been made deliberately now, in the other direction, and the cost was smaller than the
argument for keeping it: no existing call site changed, five test comparisons did, and the
application lost its leak, both of its downcast matches and about forty lines.

That is the case for building a consumer early. None of this was visible from inside.
