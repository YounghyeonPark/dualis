# A scene editor beside a 3D view

```sh
cd runtime/editor
cargo run --release                 # opens on the built-in room
cargo run --release -- scene.json   # opens on a file
```

The left pane is the scene's JSON, checked **as you type** with the same two steps
`dualis-world --check` runs — parse, then build — with parse errors carried as `line:column`,
which is what that error format was designed for. The viewport draws every placed extent as a
wireframe, live from the text before anything runs. **Run** executes the scene off the UI
thread and overlays what came back, scrubbable by frame; **Verify** runs the battery from
`dualis-world verify` and shows the report the CLI prints, with the findings count in the
window title. Drag to rotate, scroll to zoom.

## Why this is a third workspace

The measured table, one row longer:

| | external crates |
| --- | --- |
| the library, all sixteen published crates | **12** |
| `bindings/python` | 15 |
| `runtime/viewer` | 86 |
| here | **371** |

`deny.toml`, `--locked`, `wasm32` and Rust 1.78 are the library's promises and none of them
can carry a window toolkit. Nothing in `crates/` depends on this; the arrow points one way.

**Unlike the viewer, this workspace links `dualis`** — deliberately, and the difference is the
point. The viewer proves the wire format is complete by never linking the library. The editor
exists to build, run and verify scenes, which cannot be done from a file alone; it is the
first consumer of the platform verbs from a GUI. Where the two overlap, they share: the
camera, the framing, the fit and the projection are `viewer-core`'s, imported rather than
rewritten, because that arithmetic has been wrong here before and twice is enough.

## The two halves, kept apart by crate

`ARCHITECTURE.md`'s platform rules: the **authoring half** is the composition root and may
name domains; the **inspection half** must dispatch on the shape of the data, so a new physics
costs no editor edit.

- `editor-core` is the authoring half's machinery, GUI-free and tested headlessly: checking
  (parse + build), placement geometry (posed corners of every extent, through
  `Pose::point_to_world` even though no scene can state a pose yet), and thin passes over run
  and verify. It knows domains only through `dualis-world`'s public API, where that knowledge
  already legitimately lives.
- `editor` is the shell: a text box, a canvas and an event loop. Everything it paints, it
  paints from a shape — a box, points, paths, a reading — with the domain name used only as a
  label.

## What the first version does not do

No structured editing — the text *is* the model, and `deny_unknown_fields` plus the version
check are what stand between a typo and a silently different scene. No field rendering in the
viewport — a field's box carries a note pointing at the HTML report, which is stated on the
canvas rather than quietly drawn as nothing. No file dialogs, no undo beyond the text box's
own. Each of those is worth adding only after using this one reports back what actually
chafes; that method has a name here, and `FRICTION.md` is where its findings go.
