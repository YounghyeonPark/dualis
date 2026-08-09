# Scenes

Eleven worlds described as data, covering five of the library's six domains — electricity
has none yet, which is the obvious next scene. Nothing here is Rust: the
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
than none at all. `tests/scene.rs` runs all eleven on every commit and asserts one number each —
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
