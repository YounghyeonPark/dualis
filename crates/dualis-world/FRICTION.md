# What the first consumer found

`dualis-world` exists to use the SDK from outside and write down where that is awkward. A
library with no consumers is a library whose ergonomics nobody has measured, and none of the
345 tests inside the workspace can answer this question about themselves — they are written by
someone who already knows the shape.

Everything below was hit while building the smallest thing that loads a scene, runs it and
draws it. None of it is a bug in the physics except finding 6, which is.

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

**Fix, if wanted:** `impl Domain for Box<dyn Domain>` in the kernel, delegating every method.
Six lines, no change to any existing caller, and it makes `with(boxed)` work.

## 2. Domain names are `&'static str`, so they cannot come from data

Every constructor takes `name: &'static str`, and `Domain::name` returns one. A name read out
of a JSON file is a `String`. `World::build` therefore calls `Box::leak`.

The leak is bounded by the number of domains in a scene, so it is survivable rather than
dangerous. But it is the API stating that names are compile-time things, and for an
application they are exactly the opposite: they are what the user typed.

This is the friction that felt worst in practice, because it is unavoidable and it appears at
the very first thing an application does.

**Fix, if wanted:** `name(&self) -> &str` and `name: impl Into<String>` on the constructors.
It is a breaking change to eight crates and it removes a `Copy` from a hot path that does not
need one — the name is only read for reporting and lookup.

## 3. Reading state back needs the concrete type

`Simulation::domain_as::<T>` needs `T`, so the renderer knows every domain type just as the
builder does. `ScalarField` is exactly the abstraction that would avoid this — sample a field,
draw it, never ask what it is — but there is no way to get a `&dyn ScalarField` from a
`&dyn Domain`.

The result is a second `match` over the same enum, in `World::capture`, for no reason other
than downcasting.

**Fix, if wanted:** an optional `fn as_field(&self) -> Option<&dyn ScalarField>` on `Domain`,
defaulting to `None`, in the same style as the existing `as_any`. Domains that have a field
implement it in one line and a renderer becomes domain-agnostic.

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

Almost certainly an oversight rather than a decision. A one-line change.

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
code. Findings 1, 2 and 3 are the same shape: **the API is comfortable when the set of domains
is known at compile time and awkward the moment it is not.** That is a coherent position for a
library to take, and it may even be the right one — but it was not a decision anybody made,
and it is worth making deliberately now rather than discovering it again in a year.
