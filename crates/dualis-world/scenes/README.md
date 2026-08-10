# Scenes

Eighteen worlds described as data, covering all seven of the library's domains — seventeen of them
one physics at a time, and one that is actually a world. Nothing here is Rust: the
physics, the resolution, the coupling and the run length are all in the file, and the same
binary runs all of them.

```sh
cargo run --release -p dualis-world -- scenes/03-room-pulse.json
cargo run --release -p dualis-world -- scenes/03-room-pulse.json out.svg
cargo run --release -p dualis-world -- --check scenes/03-room-pulse.json
```

`--check` parses and builds without running: the format version, the domain names, a `tracks`
pointing at a node the scene defines. It reports a parse failure as `file:line:column` with the
keys that were expected, which is what an editor puts a squiggle under. CI runs it over every
scene, because it would otherwise be the one entry point nothing exercises.

Every file carries a `format` number, and **absence means 1** — which is what all seventeen here
are, having been written before the field existed. A version this build cannot read is refused
rather than half-run: `deny_unknown_fields` catches a key that was *added*, but not one whose
meaning changed, and that is what the number is for.

No second argument prints the numbers and checks them. A second argument writes an asset, and
**the extension chooses which** — because a run has several shapes and only one is a picture:

| | |
| --- | --- |
| `out.html` | A report that **picks a view per domain from its shape**, and opens in a browser. A 3D field gets two: a rotatable render and every slice |
| `out.svg` | A filmstrip: every frame on one page, one colour scale throughout so frames compare |
| `out.csv` | Every domain's scalars over time, one row per frame, units in the header |
| `out.json` | The frames themselves — fields as grids, bodies as positions in space, readings beside them |
| `out.gltf` | The geometry of the last frame, for Blender, three.js, Omniverse or any USD tool |

`.csv` is the one that reaches the domains a picture cannot. Eight of these eighteen scenes have
a domain with no field and no bodies, and for several the scalar *is* the result: `13` is about a
winding whose resistance follows its own temperature, and it drew nothing at all. As a table it
shows the feedback directly — 12.46 W at 25 °C rising to 16.01 W at 99 °C, with the resistance
going 3.11 to 4.00 Ω beside it.

`.html` is for someone who can state a simulation and does not want to decide how to draw it.
The shape of the data makes the choice, not the domain's name: scalars over time become a line
chart, a 1D field a profile over a ghost of the whole run, a 2D field an animated heatmap, and
bodies a rotatable depth-sorted scene. A new domain gets a sensible picture without the reporter
learning about it — which is the same reason `Domain::as_field` exists.

One scale per panel, fixed across every frame, in all four views. A frame that rescales makes a
quantity *look* constant while it changes by orders of magnitude, and that is the one thing a
picture of a simulation must never do.

**All four writers are `dualis-view`, and none of them is this application's.** This crate is
`publish = false`, so while they lived here they were unreachable — a consumer who could state a
simulation could not draw it. `dualis_view::{html, svg, readings_csv, to_json}` takes the frames
`dualis_scene::capture` produces, and everything the table above describes is available to any
program without going near a scene file.

`.gltf` is the one that leaves this workspace. Eight of the eighteen scenes have geometry to
export — bodies, ray paths, a 3D field as its cell centres — and the other nine are **refused with
a reason** rather than written as an empty scene: a 1D or 2D field is a graph, not something to put
in a 3D viewer, and the message says which panel and why.

The recommendation on record is to export into the rendering tools rather than rebuild them. This
is that, and it cost no dependency: glTF is JSON with the binary base64'd inside it, which is the
same reason SVG was chosen over a raster format.

Nothing generated is committed.

`cargo run -p dualis-world -- --emit-default scene.json` writes a starting point.

## Sound — `dualis-acoustic`

| Scene | What it shows |
| --- | --- |
| `01-room-mode` | The (1,1) mode of a 4.4 × 3.1 m room: the whole field breathing in and out together, one nodal line each way |
| `02-room-higher-mode` | The (3,2) mode of the same room — more nodal lines, a higher note, and a shorter run to catch it |
| `03-room-pulse` | A clap near a corner. No standing shape, so it travels, reflects off all four walls and interferes with itself |
| `16-a-room-with-a-ceiling` | The same room with a **ceiling**, released in its oblique (1,1,1) mode. 97.46 Hz, and the peak rides `\|cos(2πft)\|` to within 0.3% at 23 nodes across |

`16` is the scene `01`–`03` cannot be. A floor plan does not have the floor-to-ceiling mode
inaccurately; it does not have it. At 2.4 m that mode is 71 Hz, and the number of modes below a
given frequency grows as `f³` rather than `f²` — which is why a real room's resonances merge into
a hiss where a two-dimensional model keeps them separable much further up.

It costs what the third dimension costs: 23 × 17 × 13 nodes rather than 23 × 17, and `√3` in the
Courant limit rather than `√2`. That is the trade, stated rather than hidden, and `Room` remains
the right model for a floor plan.

## Heat — `dualis-thermal`, and the two ways domains meet

| Scene | What it shows |
| --- | --- |
| `04-heater-and-bar` | A plain channel. The heater's joules have no *place*, so they land in the bar's first cell and conduct along it |
| `05-beam-on-bar` | A shared boundary. The beam's joules do have a place, so a spot forms where it points and then spreads |

Those two are the same physics told twice, and the difference between them is the whole
argument for `Interface` and `Flux`. A conservation audit passes either way; only the picture
tells them apart.

| Scene | What it shows |
| --- | --- |
| `15-a-hot-spot-in-a-block` | One cell of a 9×9×9 aluminium block starting 60 K hot, and the spot spreading in **three** dimensions. The scene that a one-dimensional bar cannot express: heat going sideways is the whole job of a spreader plate |

`15` is the first scene whose field is a volume, and the report draws it as **every z-slice at
once** rather than as one plane with the rest behind a slider — because a viewer who never
touches a slider would see a picture of a solid that was really a picture of one plane through
it. The filmstrip has no room for that, so it draws the middle slice and *says* `z-slice 5/9` in
the label.

Its claim in `tests/scene.rs` is the one only three dimensions can make: the neighbour one cell
away along **z** is exactly as warm as the one along x. A model that resolved a plane and stacked
it, or that used the wrong spacing on one axis, fails there and passes everything else.

| Scene | What it shows |
| --- | --- |
| `11-motor-thermal-network` | 12 W into a copper winding, out through electrical steel and an aluminium housing. Three materials, two joints, and the drop across each |

The one scene with **nothing to draw**. A network's nodes have capacities, not positions, so
`as_field` declines to invent a continuum and the renderer has no panel to make — the numbers
are the output. That is also the point: a `lump` would report the motor as one temperature, and
the winding runs 13 K above the housing at half a time constant in. The thing that fails is not
the thing you can measure.

The scene test carries a matching list, so an undrawable scene has to earn its place with an
explicit check rather than passing by having nothing to check.

## Electricity — `dualis-electrical`

| Scene | What it shows |
| --- | --- |
| `17-a-busbar-with-a-notch` | A 12 × 5 × 5 mm copper busbar with a notch three cells deep, driven at 1 mV. The resistance is **solved** rather than stated: 12.39 µΩ, against 8.28 µΩ for the full section |
| `18-an-espresso-shot` | Two 30 mm baskets, identical but for the ring against the wall. Nothing states a flow rate — Darcy's law is solved on the permeability that 250 µm at `ε = 0.45` gives. The gap makes that ring 2.9× more permeable, so it runs **79% faster** and delivers **4.54% TDS against 6.48%**: more liquid, less coffee in it |

The scene that shows why a field formulation is worth the solve. `ρL/A` is a statement about a
uniform bar, and a bar with a notch is not one — so the file states a *shape* and a material, and
the resistance comes out. The test asserts two bounds rather than the measured value: above
`ρL/A` for the full section, because removing conductor cannot help, and above a naive series
estimate that treats the notched slice as a shorter bar, because the current also has to spread
back out. The excess over the second **is** the spreading resistance, and it has no closed form
for this shape.

The potential is the panel, in volts. What a picture of it shows is where the current is going.

| Scene | What it shows |
| --- | --- |
| `12-winding-heats-a-motor` | The same motor as `11`, with the watts **computed** instead of stated: 62 m of 0.35 mm² copper at 1.75 A |

`11` says 12 W. This one derives 11.93 W from a length of wire, and the two settle within a
fifth of a kelvin of each other — so the guess was a good one, and this scene is what would have
caught it if it had not been. A stated number cannot be wrong, which is another way of saying it
is not a model.

Evaluated at 90 °C, which is worth 27.5%: copper gains 0.393% per kelvin, so the same coil on a
cold bench dissipates 9.35 W. The temperature is a parameter rather than an omission, and it is
one the *simulation* does not set — a domain cannot read another's state inside the step loop,
so closing that feedback is the caller's job and `dualis-electrical` says why.

| Scene | What it shows |
| --- | --- |
| `13-winding-that-heats-itself` | The same coil with `tracks` set: its resistance follows its own temperature, so it settles 16 K hotter than one held at ambient |

**The only place in this repository where two domains are coupled by hand**, and it is worth
saying why that is allowed. Domains never read each other *inside* the step loop — they meet on
the bus, which carries amounts and not state. This runs between frames, in the code that owns
the simulation. It needed `Simulation::domain_as_mut`, which did not exist: a caller could read a
domain and not write one, so this loop was closable from nowhere at all. FRICTION 18.

The amplification is `1/(1 − g)` with `g = I²R₂₀·α·R_th`, and the test checks it as a ratio
against the same scene with `tracks` removed. Measured 1.281. Convection alone predicts 1.310;
the housing's linearised radiative conductance at 74.6 °C is 0.036 W/K against 0.294 for
convection, and including it gives 1.280 — so the 2.2% gap is radiation stiffening the heat path,
not error.

## A world — four crates at once

| Scene | What it shows |
| --- | --- |
| `14-a-world` | A laser and a lamp both heating one bar, sound crossing a room, three planets orbiting. Five domains, one clock, one bus, one audit at 1e-9 |

Every other scene here is one physics, or two meeting. This is the first with more than two, and
it exists because the crate split's claim is that domains *compose* — a claim verified in pairs
and never beyond. Building it found two things pairwise coupling cannot reach.

**A second consumer of one channel silently got nothing.** `Exchange::take` empties a channel, so
the second domain to take gets zero while every total agrees, because everything published was
consumed. Two plates under one lamp warm at the rate of one plate. That is refused now, and it
had been in every released version.

**A world's tolerance is set by its loosest domain.** The first attempt included the atom box,
which runs at 5e-2 over six picoseconds; the rest of this scene holds 1e-9 over 0.2 seconds.
Those are not reconcilable, and that is physics rather than a defect — a Lennard-Jones fluid and
a planetary orbit do not share a clock. "Physics for simulated worlds" means *a* world, not all
of them at once.

Two consequences that shape any scene like this. Declaration order is execution order, so a
producer declared after its consumer publishes into a step that has already taken — the audit
catches it as *published but not consumed*, which is how the lamp got moved above the bar. And
the plain channel is a single global resource: many producers may write to it, exactly one domain
may consume it.

## Motion — `dualis-mechanics`

| Scene | What it shows |
| --- | --- |
| `06-orbits` | Four satellites on circular orbits round one mass, each tilted out of the reference plane. The inner ones lap the outer ones, which is Kepler's third law drawn rather than asserted |
| `07-bouncing-ball` | A penalty contact with a dashpot, losing height each bounce — and a thermal lump taking the joules it dissipates |

`07` needs that lump. Without a consumer the contact publishes heat that arrives nowhere and
the kernel refuses the step, which is correct and is worth seeing once.

## Matter — `dualis-molecular`

| Scene | What it shows |
| --- | --- |
| `08-atoms-crystal` | 108 atoms at `T* = 0.15`: still on their lattice sites, rattling in place |
| `09-atoms-liquid` | The same 108 at `T* = 1.4`: the lattice is gone and they wander |

Same seed, same density, same box — only the temperature differs, so the two pictures side by
side are melting. The periodic cell is drawn as a wireframe, because it is a real boundary
rather than the edge of a picture.

## Light — `dualis-optics`

| Scene | What it shows |
| --- | --- |
| `10-lamp-on-a-mirror` | A 100 W tungsten lamp on an aluminium mirror. The mirror is worse in the blue, so the lamp's *colour* decides how much becomes heat |

`dualis-optics` has no `Domain` in it — spectra and Fresnel coefficients answer questions
rather than march a state — so this one is written in the application, like the heater and the
beam. A flat reflectance would make the colour temperature irrelevant and the whole spectral
apparatus an expensive way to multiply by a constant, which is why the test compares 2800 K
against 6500 K rather than checking one number.

## Every one of them is run by CI

A scene in this repository is a claim, and one that parses and then produces nonsense is worse
than none at all. `tests/scene.rs` runs all seventeen on every commit and asserts one number each —
chosen to be a property of the physics rather than of the file, so it would change if the
library broke and not merely if the scene were edited. Adding a scene without a claim fails
the test rather than passing quietly. CI also runs the real binary on the real files, which is
what a reader will type.

Running a scene is not a weak check on its own: the conservation audit is live for the whole
run, so one that leaked energy, created it, or left it unclaimed on a channel fails before it
draws anything.

## Three numbers that look wrong and are not

**A room is not quite the size you asked for.** `Room::of_air` quantises the height to a whole
number of cells so they stay square and the stability limit stays isotropic. A 3.1 m room at 81
cells across is 3.08 m tall, and its (1,1) mode is 67.97 Hz rather than the 67.67 Hz the
nominal dimensions predict.

**The bar's peak is not its mean.** In `04`, six joules raise 20 mm of aluminium by 1.240 K on
average, but the peak reads 1.302 K after four seconds. Heat arriving on a plain channel has no
place, so `Bar1D` puts it in the first cell and conduction has not finished levelling it.

**A bar's panel is in kelvin and its readings are in celsius**, and both say which. The field
returns what the cells hold; celsius is a conversion a view chooses, and no view here makes it.
So `04` prints `|K| 294.4821` for a bar whose hottest cell is at 21.33 °C. It said `|C| 21.3321` while this
application sampled fields itself and applied the offset in the same expression as the label —
which meant nothing could disagree with anything. `FRICTION.md` 22 is that gap, still open.

**The bouncing ball ends at zero.** `07` prints the last frame, and by one second the ball has
stopped. The bounces are in the middle of the strip, not at its end.

## Three dimensions

The physics always had them: `NBody`, `ContactSystem` and `Fluid` all carry `DVec3`, and
flattening to a plane was the renderer's simplification rather than the simulation's. Bodies
are drawn in an axonometric projection now, sorted back to front, with radius growing toward
the viewer and colour mixed toward the plate for distance. Without all three the picture is
flat however true the coordinates are.

Not isometric, deliberately: a true isometric view puts the axes at 120 degrees and makes a
cube ambiguous, which is a bad way to read a periodic box. Rooms and bars stay as they are —
`Room` is two-dimensional by construction, and the crate says why: a third dimension costs a
factor of √3 in the stability limit as well as the obvious one in cells.
