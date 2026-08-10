# The shape of dualis

The goal is to reproduce physical law in three dimensions, and to do it in a structure that can
accept physics nobody has written yet without the parts already written having to change.

That second clause is the whole design. "All of physics" is unbounded, so no architecture can
contain it by enumeration — the only thing an architecture can do is stay open. Everything below
is in service of one property: **a new physics is a new crate, and nothing existing moves.**

This document is the map. It says what each layer owns, what is built, what is missing, and which
rules cannot be broken without losing the property above. It is not a roadmap and it does not
promise dates.

---

## The three layers

```
    ┌─────────────────────────────────────────────────────────────┐
    │  ANALYSIS      what a person or an agent asks of a run       │
    │                cameras · 3D · 2D · graphs · measurements     │
    │                                              `dualis-view`   │
    └───────────────────────────▲─────────────────────────────────┘
                                │  reads
    ┌───────────────────────────┴─────────────────────────────────┐
    │  SCENE         where things are, and how they meet           │
    │                placement · interfaces · one clock            │
    │                                             `dualis-scene`   │
    └───────────────────────────▲─────────────────────────────────┘
                                │  reads
    ┌───────────────────────────┴─────────────────────────────────┐
    │  PHYSICS       what evolves, and what it conserves           │
    │                the kernel, and one crate per physics         │
    │                            `dualis-core` + six domain crates │
    └─────────────────────────────────────────────────────────────┘

As of 0.9.0 all three are crates on crates.io. Before that the upper two lived inside an
unpublished application, so a consumer could run a simulation and had no way to see it.
```

**The arrows point one way and that is load-bearing.** Analysis reads the scene; the scene reads
physics; physics reads neither. A domain that could see the scene could see another domain
through it, and the moment that is possible the crate split stops meaning anything.

### Layer 1 — physics

A domain owns state and advances it. It declares what it conserves, and it meets other domains
only on the `Exchange`, which carries **amounts** and not state.

The kernel — `dualis-core` — knows conservation, integration, scheduling, sampling and boundaries.
It knows no physics. That is the rule that makes the goal reachable at all: light, heat, motion,
sound, matter and electricity have each arrived as a crate, and none of them required the kernel
to learn anything about them.

### Layer 2 — scene

Where things are, and how they meet. `dualis-scene`.

It owns three things:

- **Placement.** A domain's state lives in its own coordinates; the scene says where those sit in
  the world. A grid gets a pose, a body set gets a pose, and a lumped model gets a *presentational*
  position — see below, because that distinction is the whole difficulty.
- **Interfaces.** Already in the kernel: `Interface` and `Flux` carry an amount *and where on a
  surface it crossed*, audited face by face rather than on the total.
- **The clock.** One `Simulation`, one schedule, one audit across everything placed in it.

`capture` turns all of that into a `Frame` at an instant, by asking each domain what it *offers*
— `as_field` for a continuum, `as_bodies` for a countable set, `readings` for scalars. It names
no domain, and `knows_no_physics.rs` demonstrates rather than asserts that: it defines a physics
inside the test file and captures it whole.

### Layer 3 — analysis

What a person asks of a run. `dualis-view`: a filmstrip as SVG, a self-contained HTML report,
a CSV of every domain's scalars, and the frames as JSON.

The rule here is proven and should be extended rather than replaced: **views dispatch on the
shape of the data, not on the name of the domain.** Scalars over time become a chart, a 1D field
a profile, a 2D field a heatmap, points a 3D scene. A new domain gets a correct picture without
this layer learning it exists.

Its tests are driven by frames written out by hand rather than by a simulation, which is the only
way to check that rule at all: a test that ran a real scene could not tell *a heatmap because the
data is a 2D grid* apart from *a heatmap because that domain was a room*.

One more rule holds throughout, and it is the easiest to break by accident: **the scale is fixed
across a run**. A picture that renormalises per frame makes a decay look like a steady state, and
per-frame normalisation is what you get if you do not think about it.

---

## Where the work is: 3D is not the default yet

The physics layer is six crates deep and dimensionally uneven. This is the honest state.

| crate | space it lives in | for the goal |
| --- | --- | --- |
| `dualis-mechanics` | **3D** — bodies, contacts, rigid rotation | done |
| `dualis-molecular` | **3D** — atoms in a periodic box | done |
| `dualis-optics` | **3D rays**; no volumetric field | rays done, fields missing |
| `dualis-acoustic` | **3D** `Hall`; **2D** `Room`; **1D** `Tube` | done |
| `dualis-thermal` | **3D** `Solid3D`; **1D** `Bar1D`; `ThermalNetwork` is a graph with no space | conduction done |
| `dualis-electrical` | **3D** `Conductor`; `Winding` is a lumped `I²R` | done |

Two observations follow, and they point in opposite directions.

**The reductions were deliberate.** `Room` is 2D because a third dimension costs √3 in the
stability limit on top of the obvious factor in cells, and it said so before anyone asked.
`ThermalNetwork` refuses to have positions because a conductance is not a distance. These are not
unfinished work. They are the cheap models that answer engineering questions, and a platform that
deleted them would be worse at its job.

**And they are not the goal.** A 3D platform needs 3D conduction, a 3D wave and a field
formulation of current, and those are three new crates or three additions, none of which requires
the kernel to change.

Both are true. The resolution is in the next section.

---

## Lumped models are reductions, not exceptions

A lumped thermal mass is a 3D conduction problem with the interior collapsed to one temperature.
That is valid when the Biot number is below about 0.1, and `LumpedMass::biot_number` exists so a
caller can find out rather than assume.

So in a platform whose goal is 3D, a lumped model is **a reduction with a stated validity**, and
it earns its place by being fast where the reduction holds. The alternative is not a purer
platform; it is an unusable one. Heat conduction through a motor housing, resolved as particles,
is picosecond steps against a two-thousand-second time constant — about 10¹⁵ steps for one answer
a graph of four nodes gives immediately.

What the platform owes such a model is not deletion but **honesty about what it is**: which 3D
problem it reduces, under what condition, and what it cannot show. `ThermalNetwork` cannot show a
hot spot, and its documentation says so.

---

## Placement: built, in two halves that must not touch

The kernel has `ScalarField` and `VectorField` — functions of position — `Interface` and `Flux`
for discretised boundaries, and now `Pose`. Until `Pose` a domain's coordinates *were* world
coordinates, implicitly, and two grids could not be placed against each other at all.

`Pose` is a rigid motion and deliberately nothing more: a rotation and a translation, no scale,
no shear, no projection. An isometry preserves every distance and angle exactly, and that is the
only class of placement a physics can be moved by without its physics changing — a conservation
law stated over a sheared volume is a different law, and a scaled metre is not a metre.

Placement has two uses, and they must not share a type:

- **Physical placement** — changes what the physics computes. Two solids in contact, a lens at a
  distance, a grid rotated against another. This is `Pose`, and it is **built**. The scene layer
  assigns it; a domain reads only its own coordinates.
- **Presentational placement** — a position given to something that has none, purely so a viewer
  can draw it. A thermal network node on a diagram. This is `Placement::marker`, and it is
  **built**, in `dualis-scene`: above the kernel, above every domain, where no physics can reach
  it.

Keeping them apart is structural rather than a naming convention. They are separated by *which
crate they live in* — if they shared a type, someone would eventually feed a drawing coordinate
into a conductance and nothing would fail loudly.

`Placement` also carries an `Extent`, and that third field is the one nobody predicted. A
`ScalarField` is a function of position and does not stop anywhere; a field that knew its own
bounds would be a mesh. So the region to sample has to come from above, and the scene is where
the size was written down in the first place.

---

## Rules that hold the whole thing up

These are not style. Each one is what makes some part of the goal reachable.

1. **The kernel must never depend on a domain.** Without this, "add a physics" means "edit the
   kernel", and the goal is a rewrite each time.
2. **No domain may depend on another.** They meet on the bus. Six domains have now been added
   without this breaking, which is the evidence that the split is real.
3. **The arrows point one way.** Analysis → scene → physics. A domain that can see the scene can
   see another domain through it. This is enforced by cargo rather than by discipline now that
   each layer is a crate: `dualis-scene` does not appear in any domain's manifest, so a domain
   reaching upward does not compile.
4. **Conservation is audited, not assumed.** The audit is what makes an unfamiliar coupling
   trustworthy, which matters more as the number of domains grows, not less.
5. **Results are bit-for-bit across platforms.** A simulation that gives a different answer on a
   different machine is not a measurement.
6. **A number in prose is a number under test.** This repository has shipped stale figures
   repeatedly; the ones that stopped recurring are the ones a test now checks.

Rule 4 had a known limit and it is closed. The tolerance used to be one number for the whole
simulation; it is now **one per quantity**, and a domain can additionally claim
`books_balance` and be checked **on its own scale** rather than against the sum of every ledger.
Both halves matter and neither substitutes for the other: the first separates schemes carrying
different quantities, the second separates domains carrying the same one.

---

## What is missing, in the order it blocks things

1. ~~**Placement in the kernel.**~~ Done: `Pose` is a rigid motion — rotation and translation,
   no scale, no shear — so two things can be positioned relative to each other. The
   *presentational* half is still missing on purpose and waits for the scene crate, where the
   physics cannot reach it.
2. ~~**Scene and analysis as libraries.**~~ Done. Both were inside `dualis-world`, which is
   `publish = false`, so a consumer who could state a simulation and run it could reach neither
   the shape of the answer nor any view of it.

   `dualis-scene` is layer 2: `Placement`, `Extent`, `Frame`, `Panel`, `capture`,
   `settle_framing`. `dualis-view` is layer 3: a filmstrip, a self-contained HTML report, CSV and
   JSON, with the view chosen by the shape of the data.

   **Neither names a domain, and both prove it by construction.** `dualis-scene`'s test defines a
   physics inside the test file and captures it whole; `dualis-view`'s tests are driven by frames
   written out by hand, which is the only way to tell "a heatmap because the data is a 2D grid"
   apart from "a heatmap because that domain was a room".

   Getting there needed **five** things the layer had been doing by matching on domain types:
   `Domain::readings`, `Domain::as_bodies`, `Simulation::domains`, `ScalarField::unit`, and
   `Placement::extent` for the region a field occupies. Each was invisible while one crate did
   everything, and each is a thing the *kernel* was missing rather than an invention of the
   layers above it.

   What is left in `dualis-world` is what an application actually is: a file format, the domain
   types that format names, and one place saying how far each field extends.
3. ~~**3D field domains.**~~ Done. `Solid3D` is conduction through a block — a seven-point
   stencil, insulated faces, `dx²/6α`, checked against the exact eigenvalue of its own discrete
   operator. `Hall` is the wave equation with a ceiling — a staggered grid, rigid surfaces,
   `dx/(c√3)`, checked against the rigid-wall mode frequencies and a second-order convergence
   rate measured across three doublings.

   Four of the six domains are three-dimensional now. What is left is `dualis-electrical`, which
   is gap 4 below, and `dualis-optics`, whose rays are already 3D and whose *fields* are not.

   Building the first one **found a gap in the layer above it**, which is what a first
   three-dimensional anything is for. The second one found nothing, which is the evidence that
   the first fix was the right shape: `Hall` needed no change anywhere outside its own crate and
   the scene format. `Extent::samples` was a pair and the sampler built its
   position as `(u, v, 0)`, so a solid would have been captured as its `z = 0` face — silently,
   because a slice of a block is a perfectly plausible picture of a block. Nothing in
   `dualis-scene` could have noticed on its own; every field it had ever been handed was flat.
   `samples` is a triple now, `PanelData::Field` carries `nz`, and the type system made all three
   view sites decide what to do about it.
4. ~~**A field formulation of electricity.**~~ Done. `Conductor` solves `∇·(σ∇φ) = 0` by
   conjugate gradients and reads `J = −σ∇φ` off it, so a resistance is a property of a shape.
   `ρL/A` comes out exactly for a uniform block; a notch comes out as whatever the notch gives,
   which is the point — spreading resistance has no closed form for an arbitrary geometry.

   It is the first **elliptic** domain here, and the first whose failure mode is a solver rather
   than a stability limit. An iterative solve stopped at its iteration cap returns a field that is
   smooth, bounded and shaped exactly like an answer, so `step` refuses one that did not converge
   and the residual is a *reading* rather than an internal number.
5. **Per-quantity tolerances.** Done. `Tolerances` and `Simulation::conservation_tolerance_for`
   give each conserved quantity its own number, with a default for the rest. A Barnes-Hut tree
   gives up exact momentum by construction while energy in a rigid room is exact to `1e-15`, and
   under one number either the momentum check refuses a correct run or the energy check stops
   being able to see anything. Both failures are demonstrated in
   `per_quantity_tolerances.rs` rather than asserted.

   **The other half is done too**, and it turned out easier than this document predicted.
   `Domain::books_balance` is an opt-in claim that a domain's ledger changes by exactly what it
   took from the bus minus what it published, and a domain that makes it is checked **on its own
   scale** every step. A domain holding a microjoule beside one holding a kilojoule can lose a
   fifth of itself without moving the sum by more than `2e-10`; on its own scale a fifth is a
   fifth.

   This document predicted the blocker would be ledgers that are honest approximations, and named
   `Room::startup_adjustment` as the example. **That was wrong.** `Room::energy` reports the
   released state's datum plus an offset chosen so the two agree, so its books balance exactly
   from the first step — the `O(h²)` correction is *inside* the number it reports and does not
   move it. Every domain in the workspace except `LumpedMass` takes the claim and passes.

   `LumpedMass` correctly declines, and it is the reason the check is opt-in rather than
   automatic: it loses heat to an environment that is not on the bus, so its ledger does not
   balance against bus traffic alone. That is a boundary being modelled, not a leak, and a check
   that accused it would be the wrong check.
6. **A renderer with depth.** The analysis layer draws with painter's algorithm on a 2D canvas.
   Real 3D content deserves real depth buffering — but content first: a better renderer of four
   points is still four points.

Beyond these, the physics itself is open-ended: electromagnetism, fluid dynamics, elasticity.
Each is a crate on the kernel. None is a change to the layers.

---

## How to judge a proposed change

Against the property at the top. A change is good if a new physics still costs one crate and
nothing existing moves. A change is suspect if it makes the kernel know a domain, lets a domain
see another, points an arrow upward, or replaces a checked number with a claimed one.
