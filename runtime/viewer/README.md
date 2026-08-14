# A native viewer for a dualis run

```sh
cd runtime/viewer
cargo run --release -- ../../out.json                       # a window
cargo run --release -- ../../out.json --snapshot out.ppm             # one frame, no window
cargo run --release -- ../../out.json --snapshot out.ppm --frame 30  # a frame worth looking at
```

Drag to rotate, scroll to zoom, space to play, left and right to scrub.

The input is what any run writes: `dualis-world <scene> out.json`, or anything that calls
`dualis_view::to_json`. Two of the shipped scenes and the `optical_bench` example are committed
under `viewer-core/tests/runs/` as fixtures, trimmed to a couple of frames each.

## Why it is a separate workspace

The same reason `bindings/python` is, and the numbers are measured rather than argued:

| | external crates |
| --- | --- |
| the library, all sixteen published crates | **12** |
| `bindings/python`, which was split out for this | 15 |
| here | **86** |

`deny.toml` gates every one of the library's twelve by licence, the lockfile is what CI builds
with `--locked`, and the same crates compile to `wasm32` and to Rust 1.78. A GPU stack cannot be
carried by any of that, and nothing in `crates/` depends on this. The arrow points one way: this
reads what a run wrote.

## Why it does not depend on `dualis`

Deliberately, and it is the load-bearing decision here. If a viewer can be written against the
**file** and nothing else, the wire format is complete. If it had needed to link the library for
something the file did not carry, the format would have been the thing to fix — and finding that
out is worth more than the convenience.

It reads all three panel shapes a run can contain, and an unknown `kind` is an **error** rather
than a panel quietly skipped. A viewer written against an older library than the run should say
so, not open a window with something missing.

## Where the arithmetic lives

`viewer-core` has no GPU dependency and holds everything that could be got wrong the same way
twice:

- **one colour scale across the whole run** — a frame normalised to itself makes a decay look like
  a steady state;
- **one framing across the whole run** — otherwise a moving body looks still and the camera moves;
- **the projection**, including the clamp that stops a point behind the eye becoming a streak
  across the window;
- **the fit** — the focal length that puts the run's bounding box in the frame.

## Two things this file already claimed and the renderer did not do

Both found by pointing it at a run shaped unlike the fixtures — `portafilter_flow`, which is tall,
thin, and spans two orders of magnitude in value from the shower screen to the spout.

**The colour scale was per frame.** `Run::scale_of` computes the range over the whole run, was
tested, and **nothing called it**: `segments` measured the values it had in hand. So water leaving
a screen clean and arriving at a spout at 83 kg/m³ rendered mid-ramp the whole way down, because
at every instant it sat halfway between that instant's lightest and darkest — the exact failure
the bullet above says this crate exists to prevent. `segments` takes the span now, and
`the_shading_uses_the_runs_scale_and_not_the_frames` checks the *consumption* rather than the
accessor.

**The camera was set up for a cube.** `Framing` normalises by the longest side, so a tall thin run
fills one axis and a fraction of the others; at a fixed distance the portafilter came out at 15% of
the frame height and **0.29% of its pixels**. Backing in does not fix it — at this field of view a
unit subject fills the frame at a distance of 0.35, which is inside its own bounding box. Distance
is the wrong knob; it sets how strong the perspective is, and that was never the problem.

`Camera::fit` sets the **focal length**, and the projection is linear in it, so one pass is exact
and there is no convergence test to get wrong. The same portafilter now covers 2.5%, and a cube, a
plate and a portafilter all land their furthest corner at 0.850 of the half-frame.

Seven tests, against the real fixtures. `viewer` is then a surface, a line pipeline and an event
loop, and draws only what it is handed.

## `--snapshot`, and what it is for

A window nobody can photograph proves the program did not panic, which is a much weaker claim than
it looks: a wrong projection, an empty vertex buffer and a silently-failed pipeline all fail to
panic too.

`--snapshot` renders one frame to a texture with no window, reads it back and writes a PPM, using
the **same** camera, the same `viewer-core` segments and the same pipeline the window uses. On the
optical bench it reports

```text
1046 of 792000 pixels carry a line (0.13%), in 3 shades
```

and three shades is the check worth having — one per field angle, so the value-to-colour mapping
arrived intact.

The lit count is measured against **the corner pixel**, not a constant. The target is sRGB, so the
clear colour is stored far brighter than the linear number the pass was given — 56,66,77 rather
than 10,14,19 — and a fixed threshold called every pixel in the image a line and reported 100%.

## What is not here

Only paths are drawn. Fields and point clouds are different pipelines and the viewer says so and
exits rather than opening a blank window.

No RTX, no materials, no shadows, no editing, no USD. Those are a different product, and the
recommendation on record is to **export into** the tools that already do them rather than rebuild
them — `dualis-scene`'s frames are a few hundred lines away from glTF, which reaches Blender,
three.js and every USD viewer.
