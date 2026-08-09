# Scenes

Fourteen worlds described as data, covering all six of the library's domains — thirteen of them
one physics at a time, and one that is actually a world. Nothing here is Rust: the
physics, the resolution, the coupling and the run length are all in the file, and the same
binary runs all of them.

```sh
cargo run --release -p dualis-world -- scenes/03-room-pulse.json
cargo run --release -p dualis-world -- scenes/03-room-pulse.json out.svg
```

No second argument prints the numbers and checks them. A second argument also writes an SVG
filmstrip — one panel per captured frame, drawn left to right, on one colour scale so the
frames are comparable with each other. Nothing generated is committed.

`cargo run -p dualis-world -- --emit-default scene.json` writes a starting point.

## Sound — `dualis-acoustic`

| Scene | What it shows |
| --- | --- |
| `01-room-mode` | The (1,1) mode of a 4.4 × 3.1 m room: the whole field breathing in and out together, one nodal line each way |
| `02-room-higher-mode` | The (3,2) mode of the same room — more nodal lines, a higher note, and a shorter run to catch it |
| `03-room-pulse` | A clap near a corner. No standing shape, so it travels, reflects off all four walls and interferes with itself |

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
than none at all. `tests/scene.rs` runs all fourteen on every commit and asserts one number each —
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
