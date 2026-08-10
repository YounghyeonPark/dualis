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
    └───────────────────────────▲─────────────────────────────────┘
                                │  reads
    ┌───────────────────────────┴─────────────────────────────────┐
    │  SCENE         where things are, and how they meet           │
    │                placement · interfaces · one clock            │
    └───────────────────────────▲─────────────────────────────────┘
                                │  reads
    ┌───────────────────────────┴─────────────────────────────────┐
    │  PHYSICS       what evolves, and what it conserves           │
    │                the kernel, and one crate per physics         │
    └─────────────────────────────────────────────────────────────┘
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

Where things are, and how they meet. This layer is the one the goal most needs and the one that
barely exists today.

It owns three things:

- **Placement.** A domain's state lives in its own coordinates; the scene says where those sit in
  the world. A grid gets a pose, a body set gets a pose, and a lumped model gets a *presentational*
  position — see below, because that distinction is the whole difficulty.
- **Interfaces.** Already in the kernel: `Interface` and `Flux` carry an amount *and where on a
  surface it crossed*, audited face by face rather than on the total.
- **The clock.** One `Simulation`, one schedule, one audit across everything placed in it.

### Layer 3 — analysis

What a person asks of a run. Views — a camera, a 3D scene, a 2D section, a graph — and derived
quantities.

The rule here is already proven and should be extended rather than replaced: **views dispatch on
the shape of the data, not on the name of the domain.** Scalars over time become a chart, a 1D
field a profile, a 2D field a heatmap, points a 3D scene. A new domain gets a correct picture
without this layer learning it exists.

---

## Where the work is: 3D is not the default yet

The physics layer is six crates deep and dimensionally uneven. This is the honest state.

| crate | space it lives in | for the goal |
| --- | --- | --- |
| `dualis-mechanics` | **3D** — bodies, contacts, rigid rotation | done |
| `dualis-molecular` | **3D** — atoms in a periodic box | done |
| `dualis-optics` | **3D rays**; no volumetric field | rays done, fields missing |
| `dualis-acoustic` | **2D** — `Room` says so in four places | needs a 3D wave |
| `dualis-thermal` | **1D** `Bar1D`; `ThermalNetwork` is a graph with no space | needs 3D conduction |
| `dualis-electrical` | **none** — scalars only | needs a field formulation |

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

## Placement: half built

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
  can draw it. A thermal network node on a diagram. This is **not built**, on purpose: it belongs
  to the scene crate, above the kernel, where the physics cannot reach it.

Keeping them apart is structural rather than a naming convention. If they shared a type, someone
would eventually feed a drawing coordinate into a conductance and nothing would fail loudly — so
the second one lives in a crate the first cannot see.

---

## Rules that hold the whole thing up

These are not style. Each one is what makes some part of the goal reachable.

1. **The kernel must never depend on a domain.** Without this, "add a physics" means "edit the
   kernel", and the goal is a rewrite each time.
2. **No domain may depend on another.** They meet on the bus. Six domains have now been added
   without this breaking, which is the evidence that the split is real.
3. **The arrows point one way.** Analysis → scene → physics. A domain that can see the scene can
   see another domain.
4. **Conservation is audited, not assumed.** The audit is what makes an unfamiliar coupling
   trustworthy, which matters more as the number of domains grows, not less.
5. **Results are bit-for-bit across platforms.** A simulation that gives a different answer on a
   different machine is not a measurement.
6. **A number in prose is a number under test.** This repository has shipped stale figures
   repeatedly; the ones that stopped recurring are the ones a test now checks.

Rule 4 has a known limit worth stating here: **the tolerance is one number for the whole
simulation.** Mix a molecular fluid at 5e-2 with a room at 1e-9 and the loose domain sets what the
strict one is checked against. A platform of many domains will need per-domain or per-quantity
tolerances, and that is a kernel change nobody has needed yet.

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
3. **3D field domains.** Conduction through a solid, and the wave equation with a ceiling. Both
   are the existing 1D and 2D schemes with one more index, and both have strong closed forms.
4. **A field formulation of electricity.** Current density and potential on a grid, so `I²R`
   becomes a consequence rather than a parameter.
5. **Per-quantity tolerances**, when a scene mixes domains whose achievable accuracies differ by
   orders of magnitude.
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
