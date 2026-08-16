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
    │                            `dualis-core` + ten domain crates │
    └───────────────────────────▲─────────────────────────────────┘
                                │  fills
    ┌───────────────────────────┴─────────────────────────────────┐
    │  INPUT         what a person designed, made into cells       │
    │                a mesh · a rasterisation · what it lost       │
    │                                             `dualis-shape`   │
    └─────────────────────────────────────────────────────────────┘

As of 0.9.0 the upper three are crates on crates.io. Before that the top two lived inside an
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

### Layer 0 — input

`dualis-shape`, and it is under the physics rather than beside it. Every domain
here takes a **structured grid**; every file a person designs is a **surface**. This is the bridge,
and it is the answer to a question the other three layers cannot be asked: *where does the geometry
come from?* Until now, from a closure written by hand, which is a fine answer for a test and no
answer at all for someone who has a part.

It depends on `dualis-units` and nothing else in the workspace, and **no domain depends on it**. It
produces a predicate, `|i, j, k| voxels.contains(i, j, k)`; a domain's `fill` consumes one. That is
the entire coupling, and it is why adding geometry cost the domains nothing — `Solid3D::fill`,
`Block::fill` and `Waves::fill` already had exactly that signature for other reasons.

**The layer exists because the cell size is a physics decision disguised as a performance one.**
Picking `cell_mm` picks three things at once: the stability limit and step count as `dx²`, the
discretisation error, and — the one with no symptom — *how much of the shape survives*. A 0.5 mm rib
voxelised at 2 mm is not a thin rib, it is gone, and the simulation runs perfectly well without it,
conserves energy, audits clean and answers a question about a different object. So the crate's
output is not only cells but a `Loss`: volume error, the fraction of the volume in boundary cells,
runs one or two cells thick, and rows the rasteriser could not decide.

Two things it deliberately does not do. It does not **choose** the cell size — that would be the
library guessing, and refusing to guess is the property the rest of this document is about. And it
does not repair, refine or simplify a mesh; it reads what was exported, measures it, and says
whether the answer is trustworthy.

What is still missing is named honestly in *What is not built* below: a grid has no **void**, so the
cells a part does not occupy are some other substance rather than nothing, and two parts voxelised
separately have no way to touch.

---

## Where the work is: 3D is not the default yet

The physics layer is eleven crates deep and dimensionally uneven. This is the honest state.

| crate | space it lives in | for the goal |
| --- | --- | --- |
| `dualis-mechanics` | **3D** — bodies, contacts, rigid rotation | done |
| `dualis-molecular` | **3D** — atoms in a periodic box | done |
| `dualis-optics` | **3D rays**; no volumetric field | rays done, fields missing |
| `dualis-acoustic` | **3D** `Hall`; **2D** `Room`; **1D** `Tube` | done |
| `dualis-thermal` | **3D** `Solid3D`; **1D** `Bar1D`; `ThermalNetwork` is a graph with no space | conduction done |
| `dualis-electrical` | **3D** `Conductor`; `Winding` is a lumped `I²R` | done |
| `dualis-porous` | **3D** — Darcy flow, advected heat and dissolution in a packed bed | done |
| `dualis-elastic` | **3D** — linear elasticity on trilinear elements, statics and waves | done |
| `dualis-em` | **3D** — Maxwell on a Yee grid, conducting walls | done |
| `dualis-fluid` | **3D** — incompressible Navier–Stokes by projection | done |
| `dualis-quantum` | **1D** `Well` — a wavefunction between walls | arrived one-dimensional, as every wave here did; higher dimensions are cells and cost, not new physics |

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
2. **No domain may depend on another.** They meet on the bus. Eleven domains have now been added
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

   Nine of the eleven domains are three-dimensional now. What is left is `dualis-optics`, whose rays
   are already 3D and whose *fields* are not — and gap 4 below, which closed `dualis-electrical`,
   is why this sentence used to say eight.

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
6. **A renderer with depth.** Partly answered, and the answer turned out not to be a depth
   buffer. Content came first, as this list said it should: three domains produce volumes now, and
   what a volume wants is not occlusion but **integration along a ray**. `dualis-view` raycasts a
   3D field — trilinear sampling, front-to-back compositing, rotatable — beside the slice montage,
   because a render shows shape and cannot be read for values while a montage is the reverse.

   Depth buffering is still absent and still not obviously needed. Bodies are points, which
   painter's algorithm sorts correctly, and a field is composited rather than occluded. It becomes
   worth doing when something here has *surfaces* — a mesh, an isosurface — and nothing does.
   `dualis-shape` now reads meshes but does not hand one to the view, so this is unchanged.
7. **Geometry from a designed file.** Done as far as one object goes: `dualis-shape` reads an STL,
   measures it, and rasterises it into the predicate a domain's `fill` already took. Nothing in the
   physics changed.

   The rest of it is not done, and the three that remain are the ones that block a real assembly:

   - **A grid has no void.** Every cell is some substance, so a part in a box is surrounded by
     another material rather than by nothing. Insulating it is a substance with a low conductivity,
     which is not the same thing and is wrong for radiation, for contact and for anything with a
     free surface.
   - **Two parts have no way to touch.** Each rasterises onto its own grid from its own bounding
     box. Placing them in one grid needs a shared origin and a rule for a cell both of them claim,
     and neither exists. This is the one that arrives first in practice.
   - **Nothing chooses the cell size.** The `Loss` report says what a choice cost, which is the
     honest half; picking a cell that resolves the smallest feature and reporting what was traded
     is a layer above this one, and would be the first place in the workspace that *guesses* — so
     it has to guess visibly.

Beyond these, the physics itself is open-ended, and two of the three that list named are done.

**Elasticity** is `dualis-elastic`: the third elliptic domain and the first whose unknown is a
vector. **Electromagnetism** is `dualis-em`: Maxwell on a Yee grid, which is the second hyperbolic
domain and the first to carry a constraint — `∇·B = 0` — that the update preserves as an *identity*
rather than to a tolerance. Neither needed the kernel or either layer to change at all, which is
the claim rule 4 makes and the ninth time it has held.

**Fluid dynamics** is `dualis-fluid`, and it was the hardest of the three to make checkable for
the reason predicted: few exact solutions, schemes that trade stability against diffusion, and "it
looks like a fluid" as the easiest wrong answer in computational physics to accept. It is built
around the three that exist — Poiseuille, Couette and Taylor–Green — each chosen to be blind to a
different mistake, with two machine-precision statements beside them that a decay rate is too
coarse to see.

That closes the list this section opened with. **The physics layer is where it was aimed**: ten
crates, each one added without the kernel or either layer above it changing, which is this
document's rule 1 held ten times. `CLAUDE.md` numbers the same rule 4, which is why a reader who has
both open should trust the wording over the number.

What was open was not a list of missing physics but a list of missing *depth* in what is here, and
that list is closed too. All three entries were the second half of a crate that already existed, and
none of them cost a new one.

**Small-strain** was the last, and the answer turned out not to be finite strain. `dualis-elastic` was
static — `∇·σ = 0`, with a `density` field its own documentation called unused — so the missing half
was **inertia**, not large deformation. `Waves` is `ρü = ∇·σ` on the same trilinear element, marched
with central differences and a lumped mass, and it is checked against two exact speeds where the
static solver is checked against four exact moduli.

The sharp statement is `c_p/c_s = √(2(1−ν)/(1−2ν))`, which holds to about `1e-8` — four orders tighter
than either speed alone, because `E` and `ρ` cancel out of it algebraically *and* the two modes share
a shape so the mesh error cancels out of it numerically. Each speed on its own is second order:
0.161%, 0.040%, 0.010% over 16, 32 and 64 elements.

Finite strain remains absent and is now a genuine *new* thing rather than a missing half. It needs
Newton iteration, and more importantly its closed forms are weak — a neo-Hookean bar's uniaxial
extension is an answer that depends on which material model was chosen, which is the closest this
workspace gets to the thing it wrote down as verifying nothing.

**Single-phase** was the second of the three, and `Solid3D` accounts for latent heat now — in **two
phases** as of the release after, where the liquid has its own conductivity and specific heat and the
check is the two-phase Neumann solution. The one-phase model is not a worse version of it: it is exact
while the liquid sits at the melting point, and switching two phases on costs first-order accuracy at
the interface, so `Substance::ice` keeps `liquid: None` and two phases are opt-in. The reduction between
them is the check on the algebra — zero superheat gives back the one-phase condition, proportionally, at
`9.4465e-3` per kelvin. `Substance`
carries `FusionProps`, `Substance::ice` is the catalogue's only entry that melts, and the check is
**Neumann's exact solution** of Stefan's problem — which is why this problem and not another: almost
nothing with a moving boundary has one.

```text
  X(t) = 2λ√(αt)      where  λ e^{λ²} erf(λ) = St/√π,   St = c(T_m − T_s)/L
```

A *position*, not a rate and not a limit. Measured at 1 mm for ice under 20 K of undercooling: 5.311
against 5.288 mm at 100 s, 10.570 against 10.567 at 400, 15.843 against 15.848 at 900. And nine times
the time gives 2.9984× the depth against `√9 = 3`, which is the `√t` law with no closed form in it at
all.

The scheme is **enthalpy**, bookkept as a temperature and a melted fraction: a cell's state is one
monotone number — `T − T_m` below, `φ·ℓ` inside, `ℓ + T − T_m` above — so the sweep adds energy and
inverts. Nothing tracks the front, nothing iterates, and energy is conserved because energy *is* the
state. It is not **apparent heat capacity**, which smears `L` over a temperature interval and lets a
step big enough to cross the interval skip the latent heat entirely; here an overshoot's remainder
lands on the far side because the inverse map says where that much energy goes, verified with ten
times the latent heat in one delivery.

Two things the tests found that are worth carrying: the profile's convergence order is **not** clean
between adjacent resolutions (0.78 then 1.84) because a fixed grid makes the front advance in a
staircase and the field depends on where it sits between cell centres — over the full fourfold
refinement it is 6.2×, between first and second order. And a checkpoint of temperatures alone is not
a checkpoint: 0 °C is ice, water or any mixture, so `checkpoint` carries the phase, which
`Schedule::Iterative` and the audit's retry both depend on.

**Single-material** was the third entry and `Solid3D::fill` is the first answer to it — and it
demonstrates the shape the rest of this list wants, because it changed nothing outside the crate
and it turned on a closed form rather than on a feature.

The number is the conductivity **on a face**, which is the harmonic mean of the two cells' and not
the arithmetic one. That is not a coarser convention: for aluminium against borosilicate the two
are 2.21 and 84.1 W/m/K, a factor of 38, and the arithmetic one short-circuits the interface. The
harmonic mean earns an *equality* — with the material interface on a cell face, the discrete chain
of face resistances is exactly `Σ Lᵢ/(kᵢA)` at every resolution, so a layered wall's resistance has
no discretisation error at all. The arithmetic mean's is first order, which is the dangerous kind:
8.9% at twelve cells, 1.0% at ninety-six, and about a thousand cells to reach in 0.1% what harmonic
has at three. A single-resolution check would have read that as a tolerance.

It also moved the stability limit off being a material property, and in the direction nobody
expects. `max_stable_dt` sums the actual face conductances now, so it is `minᵢ Cᵢ/(dx·Σ_f k_f)` —
which is `dx²/(2α)`, `dx²/(4α)` or `dx²/(6α)` according to how many axes have more than one cell
(this domain charged every shape the 3D rate, and a bar-shaped block three times the steps it
needed), and for a filled block is usually far *looser* than `dx²/(6·α_max)`. One aluminium cell
inside borosilicate is stable at exactly `k/k_face` = **75×** aluminium's own limit, because heat
cannot reach a cell faster than its worst face delivers it.

The consumer half is `regions` in the scene format and `19-a-coating-stops-the-heat`, where the
whole temperature drop sits on one face. A region that selects **no cells** is refused rather than
ignored: a mistyped bound otherwise gives a block of one material that runs, audits, renders and
answers the wrong question with nothing saying the coating was never applied.

### Every material, and the honest answer about a mixture of two

Resolving a material *per cell* answers "this part is two things". It does not answer "this part is
made of a thing that is itself two things", which is what a motor, a populated board, a printed part
and a phase-change buffer all are — and answering that turned out to need something the kernel did not
have.

**A catalogue is not the answer to "any material".** Nine entries against hundreds of thousands, and
adding a tenth does not change the shape of that. The answer is that a `Substance` is *data* —
`Deserialize` and a validator — so anything with a datasheet can be declared without the library
learning it exists. `Substance::from_name` puts the catalogue's spelling in one place so a consumer
writing a data-driven front end reads it rather than retyping it; that is the same defect three times
now, and the third copy was in the Python bindings holding five of the nine.

**`Mix` is the kernel's answer to a composite, and its shape is a refusal.** A mixture's properties
divide into three kinds and conflating them is the whole failure mode:

| | what a mixture rule can say |
| --- | --- |
| density, volumetric heat capacity, latent heat | **exact**, from conservation alone |
| conductivity, stiffness | **bounded**. No single value exists without the microstructure |
| emissivity | **nothing**. It is a property of the surface, and a mixture has no surface |

The middle row is why this belongs in the kernel rather than in a helper somebody writes once. Half
aluminium and half borosilicate conducts anywhere between 2.21 and 84.06 W/m·K — a factor of **38** —
and *both ends are attained by real arrangements of those two materials*. A library that returned
`0.5·167 + 0.5·1.114` would hand back the top of a 38-fold range as though it were a measurement. So
`Mix` returns the pair, and `as_substance` makes the caller choose inside it and refuses a choice
outside.

That this belongs to the kernel and not to a domain is worth stating as an instance of the rule: `Mix`
knows conservation and algebra and no physics. Every check on it is a domain's — a laminate in
`Solid3D` attaining Reuss across its layers and Voigt along them, a laminate in `Waves` and `Block`
attaining the shear pair, a checkerboard landing inside Hashin–Shtrikman — and none of those required
the kernel to know they existed.

**The bounds are checked against closed forms that were derived elsewhere**, which is rule 1 applied to
algebra rather than to a solver: Hashin–Shtrikman equals Maxwell–Garnett for conductivity to `6.7e-16`
and Mori–Tanaka for stiffness to `2.2e-16`, and both of those are separately derived rational functions
that agree with it at every fraction rather than in a limit.

The consumer half is `materials` and `composites` in the scene format, and scenes `21` and `22`. The
pair is the demonstration: the same wax buffer, then four fifths of it with an aluminium matrix, melting
at exactly `1/0.8` the rate — a ratio with the density and the latent heat cancelled out of it.

Two entries that used to be on this list are answered, and both the same way — not by a domain
absorbing another's job, but by a test that fails if two of them ever stop being limits of one
physics. There are three of those now:

- `fields_and_rays.rs` — a Yee grid against Fresnel's algebra, to 0.5% at eighty cells per
  wavelength, converging at second order, with a quarter-wave coating taking 14.4% to 0.23%.
- `loss_and_lumps.rs` — a field's decay in a conductor against the exact `α(ω)`, then against its
  own famous `√(2/ωμσ)` limit, beside a `Winding` answering `0.1724 Ω` at every frequency.
- `a_slit.rs` — a diffraction pattern against `sinc²`, which is **not** a comparison the field
  solution should win or lose. Scalar theory assumes Kirchhoff's boundary condition and Maxwell's
  equations do not, so the two agree in a limit: measured, the difference falls from 0.277 at one
  wavelength to 0.0057 at twelve, a factor of 48. The value is knowing *where* a closed form is
  valid, which is worth as much as the closed form.

**`dualis-optics` having rays and no fields is the one already answered**, and not by adding fields
to it. `dualis-em` *is* the field formulation, and `crates/dualis/tests/fields_and_rays.rs` is where
the two have to agree: a Yee grid and Fresnel's algebra, sharing no code and not even depending on
each other, land on the same reflectance to 0.5% at eighty cells per wavelength, converging at
second order. That is the shape the rest of this list wants — not a domain absorbing another's job,
but a test that fails if they ever stop being two limits of one physics.

---

## Rendering, and the decision not to compete on it

RTX rendering, GPU physics, USD and an interactive editor were asked about directly. Three of the
four are additive and only heavy; the fourth is not.

**GPU physics conflicts with rule 5**, and `runtime/gpu` is the first answer to that. Atomics
and warp scheduling reorder floating-point reductions, and addition is not associative — so
bit-for-bit, which `rng::tests::the_stream_is_pinned` holds to a digest, would be given up.

The arrangement that survives it: **the CPU domain is the reference and the GPU is an
accelerator**, with a test that measures how far apart they land. The stencil itself is safe —
each cell reads six neighbours and writes itself, with no reduction and so no order to depend on
— and reductions stay on the CPU, summed in index order after a readback.

The measurement is 191× on a 64³ grid, and it comes with a cost that is not about speed: WGSL has
no `f64`, so an accelerated domain is single precision against the library's double. It conserves
to `5e-11` where the CPU holds `9e-15`, which is below `Simulation`'s default `1e-9` audit — so a
scene using one has to loosen `conservation_tolerance_for(ENERGY, ..)`, and choosing that number
is choosing what the run may lose.

That crate also found the thing worth knowing about single precision here: it was never the
problem. Storing absolute kelvin made `sum − 6·centre` subtract two numbers agreeing to five
digits, keeping less than one digit of the answer. Storing the **deviation** from the reference —
which the linear stencil commutes with exactly — improved the divergence 1660× at no cost.

The other three belong **above the layers**, in workspaces of their own, for the reason
`bindings/python` already established. Measured: the library resolves 12 external crates, the
python bindings 15, and a wgpu stack **86**. `runtime/viewer` is the first of those.

And on the rest, the recommendation is not to build them. Omniverse and Isaac Sim exist and are
enormous; what this workspace has that they do not is physics that is audited, deterministic and
checked against closed forms. So the move is to be *reachable from* them — `dualis_view::gltf`
does that in a few hundred lines and no dependency, and reaches Blender, three.js and every USD
pipeline. A USD writer is the next rung, and only worth it if the readers really are Omniverse.

### Where each of the four stands

| | |
| --- | --- |
| a viewer | `runtime/viewer` — a wgpu window that reads the run **file** and does not link `dualis` |
| export | `dualis_view::gltf` — no dependency, reaches Blender, three.js and USD tools |
| GPU physics | `runtime/gpu` — 191× at 64³, single precision, CPU as the reference |
| an editor | a skeleton, at `runtime/editor`: the scene's JSON checked as you type beside a wireframe of every placed extent, run and verify as buttons, `viewer-core`'s camera. The platform section below is what it grows into, under whose rules |

**The editor was last for a reason that is now gone.** An editor writes files, and the scene
format is `dualis-world`'s, which is `publish = false` precisely because a file format is a
compatibility promise this one was not ready to make — it had already changed under its own users
once, silently, when `mode` became `release`.

It carries a `format` number now, with the narrow promise attached: **within one version, a file
that loads today loads tomorrow.** Absence means 1, which is what every scene written before the
field existed is, so old files are readable by construction rather than by special case. A version
this build does not know is refused rather than half-run — `deny_unknown_fields` catches a key that
was *added* but cannot catch one whose meaning changed, and that is exactly the gap a version
number closes.

`--check` is the other half: parse and build without running, reporting a parse failure as
`file:line:column` with the keys that were expected. That is what an editor needs while somebody is
typing, and CI runs it over every scene so it is not the one entry point nothing exercises.

What is left for an editor is a GUI, and that is a product decision rather than an architectural
one. The format is ready to be edited — and `runtime/editor` is the first GUI over it: a third
workspace, linking `dualis` where the viewer deliberately does not, reusing `viewer-core`'s
camera rather than writing that arithmetic a third time, and split into a GUI-free
`editor-core` (checked, placed, run, verified — tested headlessly) and a shell that paints only
shapes. The platform section below is the set of rules it was built under.

---

## A platform of its own, and the rules that keep it open

The section above says what not to build and where to be reachable instead, and that answer is
complete for everything an external tool can express. **Most of what this workspace holds, no
external tool can express.** There is no USD schema for a `Violation` — nor for a `Ledger`, an
impedance boundary, a multirate schedule, a radial distribution function, a mode shape, a
spectrum, a Biot number or a `Loss` report. Subtract geometry and appearance from what dualis
models and nearly everything left has no authoring or inspection surface anywhere. Interop cannot
reach that; only a native platform can — authoring the scene format, running it, inspecting the
result, and verifying it.

So the direction is both, with the boundary stated once: **what an external tool can express is
reached through files** — glTF today, a usda writer as the next rung — **and what none of them
can express is the native platform's reason to exist.**

### Verification is the platform's identity

The most dangerous errors this workspace has found were all passed by the audit, and that list is
the platform's specification rather than a feature backlog:

- **The coupling window.** A 10 ms window delivered exactly the right joules and read the peak
  temperature 12% low — the one error class no conservation check catches.
- **A collapsed convergence order.** Refining the grid halved the error instead of quartering it,
  which is what found both acoustic boundary defects; no single run could have.
- **Rasterisation loss.** A rib finer than the grid does not fail, it disappears. `Loss` measures
  that, and nothing yet turns the measurement into a verdict a person reads.
- **Margins, not verdicts.** How far each dt sits from `max_stable_dt`, and how many digits each
  audit residual has left against its tolerance, per quantity — a pass with no margin is a
  different fact from a pass.

Every one of these needs *reruns* — a second resolution, a halved window — so a tool that only
reads a run file cannot do it. The platform links `dualis` where the viewer deliberately does
not, and that is the structural reason the panel is native: an incumbent without the audit and
the determinism underneath cannot build it. Determinism adds one thing free of charge: a run
digest, so "this result reproduces bit for bit on any machine" is a checkable fact on a report
rather than a claim. And where a run stands beside a measured system — the digital-twin case —
the measurement takes the closed form's seat, and the same battery is the comparison.

Of the platform's verbs, all three exist inside `dualis-world`: `--check` parses and builds
without running, reporting `file:line:column`; a scene runs to a file; and `verify` runs the
battery above — built before any GUI, because a CLI that earns trust is the platform with the
smallest possible surface. Its maiden run did what the battery is for: the default scene's peak
converged at first order where the scheme is second, and the cause measured out to be the room's
**height**, quantised to whole cells and converging on the stated 3.1 m at first order — a
contaminant the battery's own documentation names, found in the first file pointed at it.

### The role boundary: the format graduates, the stranger remains

`dualis-world` holds three roles today, and they pull apart. It is the stranger — the consumer
whose naivety produces `FRICTION.md`, which requires it to stay unpublished. It owns the scene
format. And it is the runtime that builds, runs and draws what the format states. The moment a
platform exists, the second and third are product, and a product cannot be played by an actor
whose job is not knowing things.

The resolution is the move already made at 0.9.0, when scene and view left this same crate:
**the format and its builder graduate into a published crate, and the platform and the CLI both
consume it.** One schema, two front ends, no forks. That crate is the one place above the domains
allowed to name them — it is the composition root, which is already true today inside
`dualis-world` rather than a permission being invented.

### Three rules, so the platform does not undo the library

1. **Its own workspace.** Measured three times now: the library resolves 12 external crates, the
   python bindings 15, the viewer's wgpu stack 86 — and a GUI stack is heavier than a viewer's.
   The library's lockfile, licence gate, WebAssembly and MSRV promises cannot carry it. Unlike
   the viewer it links `dualis`, because it has to run things; the rendering half is
   `viewer-core`, reused rather than rewritten.

2. **The inspection half must not enumerate domains.** A panel hard-coded per domain means the
   eleventh physics costs a platform edit, and the property this document exists to protect is
   lost one layer up — the bottleneck has moved, not gone. The mechanism is the one `dualis-view`
   already proves: dispatch on the shape of the data. A new domain offering an existing shape is
   drawn for free; a genuinely new shape costs one panel *kind*. The authoring half is different —
   it is the composition root and names domains — and the difference stays structural for the
   same reason `Pose` and `Placement::marker` do not share a type: the two halves do not share a
   crate.

3. **The stranger's method is inherited, not retired.** `FRICTION.md`'s findings have come from
   five sources; building the platform is the sixth — direct manipulation. Dragging two parts
   together arrives at "two parts have no way to touch"; a cell-size control arrives at "nothing
   chooses the cell size"; drawing an enclosure arrives at "a grid has no void" — the gaps this
   document already names, found again from outside, in the order a user meets them. The
   obligation that transfers is the writing-down: every place the format or the library cannot
   say what a user just tried to do becomes a finding before it becomes a feature.

---

## How to judge a proposed change

Against the property at the top. A change is good if a new physics still costs one crate and
nothing existing moves. A change is suspect if it makes the kernel know a domain, lets a domain
see another, points an arrow upward, or replaces a checked number with a claimed one.
