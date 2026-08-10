# What the first consumer found

`dualis-world` exists to use the SDK from outside and write down where that is awkward. A
library with no consumers is a library whose ergonomics nobody has measured, and none of the
434 tests inside the workspace can answer this question about themselves — they are written by
someone who already knows the shape.

Everything below was hit while building the smallest thing that loads a scene, runs it, couples
two domains over a plain channel and two more over a shared boundary, and draws the result. None of it is a bug in the physics except finding 6, which is — and which no test inside the
library could have found, because none of them was checking a rate.

**Sixteen of the twenty-two are fixed**, and six are recorded rather than actioned. The reasons
differ and are given in each: one because the kernel already refuses the mistake it describes,
one because it is documented rather than changed, and the rest on scope. The entries are
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

**Not fixed, and declined on scope rather than deferred.** Sharing it means either a `dualis-plot`
crate or a feature-gated module in the facade, and either way it is a public API for drawing
that would have to be supported, versioned and documented — for a workspace whose stated
scope excludes rendering. Two overlapping private renderers is the cheaper mistake for now.

Revisit when there is a second consumer. One application writing its own hundred lines of SVG
is not evidence; two would be.

**There is one now, and it declined to be evidence.** The sizing tool behind findings 19–21 draws
nothing at all — it prints five numbers and exits, because a sizing question wants a settled
answer rather than a picture. So the second application did not want the plotting, which leaves
this finding exactly where it was rather than settling it. That is a real answer and not a
dodge: the case for a `dualis-plot` needs two consumers that want to plot, and there is still
one.

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

Measured through this crate's own harness — `World::run` over 0.02 s in forty advances, each
multirate-subcycled. Stepping the room directly at its stability limit gives the same second
order and smaller absolute numbers, because the substep pattern differs. The *rate* is the
claim; the column is only reproducible where it was taken.

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

## 11. `as_field` covers half the domains, and there is no counterpart for the other half

Finding 3 gave `Domain` an `as_field`, and a renderer stopped needing to know what a room or a
bar was. Then the scene format grew orbits, a bouncing ball and a box of atoms, and the
renderer went straight back to `domain_as::<NBody>`, `domain_as::<ContactSystem>`,
`domain_as::<Fluid>`.

Not because the fix was wrong. Those three genuinely are not fields: they are a countable
number of bodies at places, and rasterising them would invent a continuum they do not have.
`as_field` returning `None` for them is the honest answer.

**Fixed.** `Domain::as_bodies` returns `Option<&dyn Bodies>`: count, position, a value to colour
by, and a *real* wall or `None`. Sixty lines of downcasting in this crate became one call.

It sat here for months and was paid the moment the layers were split apart -- a scene layer that
must name three physics to find out where anything *is* needs editing every time a fourth
arrives, and that is the one thing the structure exists to prevent. The trait draws a line the
old code could not: a periodic cell is a boundary condition and the domain reports it, while an
orbit's box is a property of the picture and nothing physical sits at its edge.

## 12. Four mechanics domains had never opted into `as_any`

`NBody`, `TreeNBody`, `RigidBody` and `Rolling` all returned the default `None`, so
`Simulation::domain_as` could not reach any of them. The orbit scene ran, conserved, and drew
nothing at all — `bodies()` returned `None` and the panel was silently dropped.

**Fixed.** Four one-line impls. The same shape as finding 7: an opt-in that everything with
state to show is expected to take, not taken in the places nobody had needed yet. Optics,
thermal, acoustic and molecular had all taken it, because tests inside the workspace had
reached for them; mechanics had not, because none had.

Worth noticing how it failed. Not a compile error and not a violation — a picture with nothing
in it. A renderer that skips what it cannot read is reasonable on its own and produces the
least debuggable outcome there is.

## 13. `Schedule::Multirate` front-loads a coupled quantity, and the audit cannot see it

**The most serious finding in this file, and it is in the kernel.** Found by building a lumped
plate under a lamp against the *published* 0.1.0 and comparing it to the closed form of its own
scheme — not by reading the code.

`Simulation::sweep` steps one domain to completion before the next. A quasi-static publisher is
never subcycled, so it puts a whole outer step's joules on the bus once; a subcycling consumer
then calls `Exchange::take` on its **first** substep and takes all of them. Every joule of the
interval is deposited at its beginning and decays for the rest of it.

So subcycling does not refine the answer. The limit of `u ← u·gⁿ + (P·dt/C)·g^(n−1)` with
`g = 1 − h/τ` as `n → ∞` is `u·e^(−dt/τ) + (P·dt/C)·e^(−dt/τ)`, which is not the solution: the
error is first order in the **outer** step and independent of the substep entirely.

```text
  outer dt   staggered   multirate    analytic    stag err   multi err
    300 s    303.670     300.033      301.920      1.75       1.89
    150 s    302.678     301.257      301.920      0.758      0.663
     75 s    302.276     301.758      301.920      0.356      0.163
```

At 300 s the schedule chosen *for* accuracy is the worse of the two, with the errors on opposite
sides. And every one of those runs passed the conservation audit at around 1e-12: the total that
crossed is exactly right and only its distribution in time is wrong, and a `Ledger` has no
representation for *when*. This is the time-domain twin of the reason `audit_transfers` had to
become a per-face check in space.

**Fixed.** The recommendation was to document it and wait for a second consumer, and the
argument changed on inspection: a *shipped* scene already had it. `04-heater-and-bar` runs a
quasi-static heater beside a bar that subcycles hard, so this was not a speculative API.

`Exchange::take_share(channel, dt)` is the fix. `Simulation::advance` tells the bus what interval
the sweep covers, and a subcycling consumer asks for its substep's share instead of the lot.
`Bar1D` and `LumpedMass` use it. At a 300 s outer step multirate went from 1.89 K of error to
0.304 K — from the worse of the two schedules to fourteen times better than the alternative.

The share is apportioned against the time **remaining**, not against the whole interval, and that
is what makes it exact: handing out `A·dt/T` and reducing both leaves `A/T` unchanged, so the last
substep receives the remainder and the channel ends empty. Against the whole interval instead,
`n` shares leave `O(n·ε·A)` stranded, and `audit_transfers` uses an *absolute* tolerance that
would eventually refuse a run that was arithmetically fine. Even so the comparison needs a slack
of `1e-12` of the interval, because three substeps of a third do not sum to one in binary and an
exact test misses the final share.

Two things the fix exposed, both mine:

- **My first reference was wrong, and refining the step made the disagreement worse.** I compared
  against `T_a + (P/hA)(1 − e^(−t/τ))`, which is the closed form of *linear* loss, on a plate whose
  `Environment` also radiates. So the real equilibrium was lower, the run sat below the reference,
  and finer steps moved *away* from it. Explicit Euler must overshoot, so a scheme sitting below
  a reference and diverging from it on refinement is a reference that is wrong. Fixed by setting
  the emissivity to zero — which is buildable from the prelude only because of finding 15.
- **A test of mine was averaging the render panel again.** `a_heater_and_a_bar_meet_on_the_bus`
  measured the bar's mean from the sampled field rather than from the cells, and passed at 1e-6
  only because the field was nearly uniform by the end. Changing *when* the heat arrives changed
  the profile enough to expose it at 4.1e-6. That is finding 10, in a test written after finding
  10 was written down.

## 14. `Report` cannot be named without going through a module

`Report` is the return type of `Simulation::advance`, the most-called method in the library, and
it was reachable only as `dualis::core::sim::Report`. You could use it inferred; you could not
write a function signature over it. Found by a consumer wanting a helper that takes one.

**Fixed.** Added to `dualis-core`'s root re-export and to the prelude. Two lines.

## 15. `Substance` was in the prelude and could not be built from it

`Substance::bulk` leaves `thermal: None`, which `LumpedMass` rightly refuses to step. Supplying
one needs `ThermalProps` and two unit types, and the prelude exported twenty-five unit types
without those two. So the material set reachable from the prelude was the three catalogue entries
and one that cannot be used.

**Fixed.** `ThermalProps`, `MechanicalProps`, `AcousticProps`, `ThermalConductivity` and
`ThermalExpansion` are in the prelude.

## 16. The first error a consumer ever saw was ungrammatical

`Violation::at` builds the cases that are not a before/after comparison — a substance with no
heat capacity, an iteration that never converged — and carries a *message* in `quantity`. `Display`
had no branch for them, so it read the message as a quantity name:

```text
substance has no heat capacity is not conserved at plate: inf
```

Correct use of the constructor, correct field, unreadable sentence.

**Fixed.** A third `Display` branch for the `tolerance == 0.0 && before == after` case:
`at {site}: {quantity} ({before})`.

---

## 17. `ThermalNetwork::nodes()` returned a count that could not be turned into anything

A `Node` handle can only come from `node` or `node_losing_to`. That is deliberate and it is the
reason a link naming a node that does not exist is unrepresentable — the case the conservation
audit is structurally blind to, since a link's `+q` and `−q` cancel identically.

But it also meant a caller holding a network it did not build had no way in. `nodes()` gave a
count, `label` and `temperature` both needed a handle, and there was no way to obtain one except
`node_named`, which needs a name you do not have. The count was information you could not act on.

Found the first time this application tried to print a network's node temperatures — the one
thing the domain exists to produce.

**Fixed.** `handles() -> impl Iterator<Item = (Node, &str)>`. A dozen lines, and it makes `nodes()`
mean something.

The same shape as findings 4, 8, 9 and 10: **the API is comfortable when the parts are known at
compile time and awkward the moment they are not.** A caller writing `let winding = net.node(…)`
never noticed, because they were holding the handles already. That is now five of seventeen
findings with one underlying cause, and the count is the argument.

---

## 18. A caller can read a domain and not write one, so a feedback loop is unclosable

`Simulation::domain_as::<T>` hands back a `&T`. There was no `&mut T`.

That is fine for the reader it was built for — a test asserting a profile, a renderer sampling a
field. It is not fine for the one thing an application genuinely has to do that the library
cannot: **close a feedback loop the bus cannot carry.**

Copper's resistance rises 0.393% per kelvin, so a winding that heats up dissipates more. The
temperature lives in `dualis-thermal`, the resistance in `dualis-electrical`, and neither can see
the other's state — correctly, since domains meeting only on `Exchange` is the property the crate
split exists to hold. The caller between frames can see both. It could read the temperature and
had no way to write the resistance, so the loop was closable from *nowhere at all*.

**Fixed.** `Domain::as_any_mut` and `Simulation::domain_as_mut`, mirroring the existing pair.
This does not weaken the rule: it is about what happens inside `step`, and this runs between
frames in code holding `&mut Simulation` that could drop the domain and rebuild it — so denying
it a write was never protecting anything.

The same opt-in hazard as findings 7 and 12, and this time it was handled in the same change:
`as_any_mut` defaults to `None`, so a domain that forgets it is silently unwritable rather than
broken. Every domain that implements `as_any` got the counterpart beside it — twelve at the time, and
fifteen sites now including the two the application defines — along with the `Box<dyn Domain>`
forwarding impl and the kernel's own front-page example.

Scene 13 is the loop closed, and it measures what the application-level version costs.

## 19. Dimensioned constructors have no unit a person types

Building a three-node network was six lines of `Volume::from_si(x * 1e-6)`. `Length` has `mm`,
`m` and `cm`; `Volume` had only `from_si`, and `Area` had none at all.

The type system's whole promise is that a factor of a thousand appears in exactly one place — a
unit-bearing constructor. Where there is no such constructor the factor moves to the call site,
which is precisely where the promise said it would not be.

**Fixed.** `Volume::cm3`/`mm3`/`m3`/`litres`, `Area::cm2`/`mm2`/`m2`.

## 20. `runaway_current` wanted a number the network held, and the hand version was wrong

`Winding::runaway_current(g)` takes the conductance of the whole path to ambient. A caller with a
`ThermalNetwork` had to assemble it — `1/(1/K₁ + 1/K₂ + 1/(h·A))` — out of the links and
environment the network is already holding.

Tedious, and worse than tedious. That formula is convection-only: the housing also radiates at
its operating temperature, so it gave 0.203 W/K where the truth is 0.220, and a threshold of
4.11 A where the truth is 4.28. **The library's own documentation was quoting the wrong number**,
and a network with one more joint, or an environment on an interior node, would have been wrong
by more with nothing to say so.

**Fixed.** `ThermalNetwork::path_conductance(node, at)` takes the slope of its own solved
balance, so every path out is in it by construction rather than by the caller remembering.
Shipped as 0.7.0 on its own rather than batched, because a wrong number in published
documentation is something a reader copies.

## 21. The electro-thermal fixed point is eight lines every consumer writes

Dissipation depends on winding temperature, which depends on dissipation. Closing it is a loop
over `steady_state` and `dissipation_at` — eight lines, and the same eight lines for anyone who
wants a settled coupled answer.

**Not fixed, and declined rather than deferred.** The library cannot close it: a domain has no
way to read another's state, and that is the property the crate split defends. A helper would
have to depend on both crates, which only the facade may do, and a physics-specific solver in
the facade is a worse precedent than eight lines in a consumer.

The general form is a state channel on `Exchange`, which stays undecided — and this entry is
evidence *against* it. Twice the hand-written loop has agreed with the stepped answer: 99.0 °C
from the sizing tool's fixed point, 99.02 °C from scene 13's marching. Nothing has yet needed
the kernel to carry state.

---

## 22. A view wants a unit the domain does not store, and there is nowhere to say so

`Bar1D` holds kelvin. Every picture of a bar in this application is in celsius. While this crate
sampled fields itself that was four characters — an offset carried beside the extent — and the
panel came back labelled `"C"`.

Splitting the scene layer out took the sampling with it, and the offset had nowhere to go.
`dualis-scene` cannot carry it: a conversion is a choice about presentation, and a layer that
knows no domain cannot know that this one field wants 273.15 subtracted and the pressure beside
it does not.

**Not fixed, and the shape of the answer is not obvious.** `ScalarField::unit` at least makes it
*detectable*: a panel now says `K` and a reading says `C`, and a test compares them by converting
from the two declared units rather than assuming they match. That is the honest state and it is
better than the old one, which was a silent relabelling — the offset and the label were applied
in the same expression, so nothing could have disagreed with anything.

What it wants is a view-level unit conversion with the *dimension* known, so that asking for
celsius from a pressure field fails rather than subtracting 273.15 from a pascal. The library has
that machinery — `Temperature` is a dimensioned type — and what it does not have is a way for a
`&dyn ScalarField` to say which dimension it returns. Recorded rather than guessed at.

## What this says about the exercise

Twenty-one findings, and the source has shifted four times. The first twelve came from writing the
application; four came from **running the two subagents built out of what the first twelve
taught** — one hunting outcomes that come out empty, one building against the *published* 0.1.0
rather than the working tree. The seventeenth came from a third source again: adding a domain the
library did not have, and finding that the *new* API had the old shape.

Sixteen are fixed. That line said "ten" until a test counted them, which is the failure
`prose-auditor` exists for and the second time this file has been the one carrying it — and the
count is now checked by `friction_counts.rs` against the headings, because the author evidently
cannot do it reliably and a reader cannot do it at a glance.

Findings 11 and 22 are the same event seen twice. Pulling the scene layer into its own crate
**paid** 11, which had sat unfixed for months, and **created** 22 in the same edit: a conversion
that was invisible while one crate did everything became a thing somebody has to declare. That is
what a layer boundary does — it turns assumptions into statements, and some of the statements
turn out to be missing.

Finding 13 is the one that changes the ledger on this exercise. Every earlier finding was
ergonomic or a defect in a domain; that one is a first-order accuracy defect in the *kernel's*
scheduler, in the schedule whose whole purpose is accuracy, invisible to the audit, and it was
found by comparing a coupled run against the closed form of its own recursion. Not by reading
the code — the code is doing exactly what it documents.

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

All six domains have scenes now, so the gap this section used to name is closed. What is left
is smaller and more specific.

**`TreeNBody`, `RigidBody` and the rest of mechanics.** Four types took `as_any` in this pass
but only `NBody` and `ContactSystem` have scene variants, so Barnes-Hut and rigid rotation are
still driven only from inside.

**A room in three dimensions.** `Room` is two-dimensional by construction and the crate says
why, so this is a limit of the physics rather than of the scene format — but it does mean the
only fields anyone can draw are flat, while every *body* domain is now drawn in space.
