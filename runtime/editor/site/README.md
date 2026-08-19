# The editor, in a browser

```sh
cd runtime/editor
cargo build --release -p editor-wasm --target wasm32-unknown-unknown
cp target/wasm32-unknown-unknown/release/editor_wasm.wasm site/
cd site && python -m http.server 8000      # any static server; wasm needs http, not file://
```

Then open `http://localhost:8000`. Edit the scene on the left, **run**, scrub the frames, and
**verify**. The physics is running in the tab: there is no backend, and nothing is uploaded.

## CAD

**add CAD**, or drop the files on the window. A `.json` dropped here is a scene and anything else
is a part, because guessing wrong about which is which is worse than having two ways in.

The bytes go into the wasm module and stay there. A scene's `parts` names a file:

```json
"parts": [
  { "stl": "bracket.stl", "material": "aluminium" },
  { "stl": "insert.stl",  "material": "copper" }
]
```

On a machine that string is a path. Here there is no filesystem, so it is the **name the file was
dropped under** — and the two are the same scene, which is the point rather than a convenience.
`World::build` reads from a disk, `World::build_with` reads from wherever it is given, and a test
asserts the same STL voxelises to the same block cell by cell from either source. If those could
diverge this page would be a demo: a thing that looks like the product and answers a slightly
different question.

Naming a file that is not here is refused, and the refusal lists the names that **are** — somebody
with three files in a tab and one misspelling should be told which spellings exist rather than
"not found".

### fit grid

The step between dropping a file and having a scene. A `parts` entry needs `cells` and `cell_mm`,
and nothing about a CAD file says what those should be — while getting them wrong does not fail:
a part finer than the grid rasterises to **nothing**, and the run is perfectly well behaved about a
different object.

**fit grid** rasterises the uploaded parts at every candidate cell size and prints what each one
cost. The ladder is a statement about the assembly rather than about millimetres — the thinnest
dimension of any part gets 1, 2, 4, 8 … cells — and every number in the table was measured rather
than predicted, because `Loss::volume_error` is a lattice-point count that can get *worse* under
refinement and change sign. The rule for the recommended row is one sentence, so it can be
disagreed with by reading the table: the coarsest grid on which every part is present, nothing was
undecidable, and the worst boundary fraction is under a half. The fragment it prints goes straight
into a `block` domain.

The same thing from a terminal is `pantometry-world fit part.stl [more.stl …] --cells N`.

### what each part is, and assembling them

Dropping CAD lists each file with a **material menu**, and the menu is asked for rather than typed
out here: it is `Substance::CATALOGUE` itself, fetched across the boundary, so a page cannot offer
a name the builder then refuses. A part nothing was chosen for is aluminium, which is what a
`block` defaults to.

**assemble** writes the scene. It takes the grid the fitter measured and the material each file was
tagged with, and composes the `block` domain that puts them all on one mesh — everything outside
the parts becoming void, which is what an assembly in air *is*. That is the last of the four steps
the browser was for: upload, say what each thing is, put them together, watch it run.

Everything a part costs is reported the way the CLI reports it: filled cells, volume error, how
much of it is in boundary cells, thin runs, triangles under a cell. A rib finer than the grid does
not fail, it *disappears*, and the run is perfectly well behaved about a different object.

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

The CAD path is in there too, end to end in a host that is not a browser: the test writes a
**binary STL from the format's own layout**, hands it across the boundary, and checks that a
20 mm box on a 20 mm block fills all 64 cells — which is the claim that says the bytes were read
as millimetres and landed where the scene put them, rather than merely parsing. Then it forgets
the file and asserts the same scene is refused again, because a page that lets somebody remove a
part and keeps running the old bytes shows a stale answer as a fresh one.

That test earned its place immediately: it found that `run` would happily execute the last text
handed to `check` **even when the check had failed**, because the page's disabled button was the
only thing stopping it. The module refuses now, with the check's own message.
