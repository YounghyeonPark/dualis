# Scenes

Five worlds described as data. Nothing here is Rust: the physics, the resolution, the coupling
and the run length are all in the file, and the same binary runs all of them.

```sh
cargo run --release -p dualis-world -- scenes/03-room-pulse.json
cargo run --release -p dualis-world -- scenes/03-room-pulse.json out.svg
```

No second argument prints the numbers and checks them. A second argument also writes an SVG
filmstrip — one panel per captured frame, drawn left to right, on one colour scale so the
frames are comparable with each other. Nothing generated is committed.

`cargo run -p dualis-world -- --emit-default scene.json` writes a starting point.

| Scene | What it shows |
| --- | --- |
| `01-room-mode` | The (1,1) mode of a 4.4 × 3.1 m room: the whole field breathing in and out together, with one nodal line each way |
| `02-room-higher-mode` | The (3,2) mode of the same room — more nodal lines, a higher note, and a shorter run to catch it |
| `03-room-pulse` | A clap near a corner. No standing shape, so it travels, reflects off all four walls and interferes with itself |
| `04-heater-and-bar` | Two domains meeting on a plain channel. The heater's joules have no *place*, so they land in the bar's first cell and conduct along it |
| `05-beam-on-bar` | Two domains meeting on a shared boundary. The beam's joules do have a place, so a hot spot forms where it points and then spreads |

## Every one of them is run by CI

A scene in this repository is a claim, and one that parses and then produces nonsense is worse
than none at all. `tests/scene.rs` runs all five on every commit and asserts one number each —
chosen to be a property of the physics rather than of the file, so that it would change if the
library broke and not merely if the scene were edited.

Running a scene is not a weak check on its own: the conservation audit is live for the whole
run, so one that leaked energy, created it, or left it unclaimed on a channel fails before it
draws anything.

## Two numbers that look wrong and are not

**A room is not quite the size you asked for.** `Room::of_air` quantises the height to a whole
number of cells so they stay square and the stability limit stays isotropic. A 3.1 m room at 81
cells across is 3.08 m tall, and its (1,1) mode is 67.97 Hz rather than the 67.67 Hz the
nominal dimensions predict. The scheme is right; the room is 2 cm shorter than the file says.

**The bar's peak is not its mean.** In `04`, six joules raise 20 mm of aluminium by 1.240 K on
average, but the peak reads 1.302 K after four seconds. Heat arriving on a plain channel has no
place, so `Bar1D` puts it in the first cell, and conduction has not finished levelling it. That
gradient is the difference between `04` and `05` — and the reason the spatial channel exists.
