# The editor, in a browser

```sh
cd runtime/editor
cargo build --release -p editor-wasm --target wasm32-unknown-unknown
cp target/wasm32-unknown-unknown/release/editor_wasm.wasm site/
cd site && python -m http.server 8000      # any static server; wasm needs http, not file://
```

Then open `http://localhost:8000`. Edit the scene on the left, **run**, scrub the frames, and
**verify**. The physics is running in the tab: there is no backend, and nothing is uploaded.

## Why there is no server

The whole library compiles to `wasm32-unknown-unknown` — kernel, eleven domains, the scene
format, the builder and the view layer — so every question this page asks is answered locally by
the same code the CLI runs. A backend would only be a slower copy of what is already in the tab.

## Two shells, one core

`editor-core` holds checking, placement geometry, running, verifying and the rule that decides
what colour a field's cells are. The native window and this page are both shells over it, and
the camera is `viewer-core`'s in both — the same argument that kept one camera keeps one colour
rule. What is only in `editor-wasm` is marshalling and one projection loop; what is only here is
a canvas, a text box and an event loop.

## Checked without a browser

`node editor-wasm/selftest.mjs` instantiates the same `.wasm` the page fetches and exercises
every export: a scene checks, a bad scene reports `line:column`, an unchecked scene is refused,
a run produces its frames, the field becomes cells with **Planck's colours** and the canvas says
so, a cool field says *false colour* instead, and the battery returns its report. The module
imports nothing, so any host can run it — and a page nobody clicks proves as little as a window
nobody photographs.

That test earned its place immediately: it found that `run` would happily execute the last text
handed to `check` **even when the check had failed**, because the page's disabled button was the
only thing stopping it. The module refuses now, with the check's own message.
