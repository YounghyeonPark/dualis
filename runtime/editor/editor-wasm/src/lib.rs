//! The same editor machinery, compiled for a browser — no server, no bindings generator.
//!
//! The whole library compiles to `wasm32-unknown-unknown`: kernel, eleven domains, the scene
//! format, the builder and the view layer. So a scene can be checked, run, verified and drawn
//! entirely inside a page, and nothing here talks to a backend because there is nothing for a
//! backend to do. The physics that answers the question is the physics that shipped.
//!
//! # A second shell, not a second implementation
//!
//! Checking, placement geometry, running, verifying and the colour of a field's cells are
//! [`editor_core`]'s — the same machinery the native window uses and the same tests cover. The
//! camera and framing are `viewer-core`'s, as they are in the native shell and the HTML report.
//! What is only here is marshalling and one projection loop. That split is the reason this
//! crate is a workspace member beside `editor` rather than a workspace of its own: two shells
//! over one core is a structure, and two copies of a colour rule is a defect waiting.
//!
//! # Why raw exports rather than `wasm-bindgen`
//!
//! Because the whole surface is **text in, text out**: a scene is JSON, a run is JSON, a report
//! is a string. That needs one allocator pair and a length prefix, which is the code below, and
//! it costs the workspace no new dependency and the reader no new toolchain — `cargo build
//! --release -p editor-wasm --target wasm32-unknown-unknown` is the entire build. A bindings
//! generator earns its place when the boundary carries structs; this one carries bytes.
//!
//! # The contract with the page
//!
//! Every export returns a pointer to a **length-prefixed** buffer: four little-endian bytes of
//! length, then that many bytes of UTF-8. The page reads the length, decodes the body, then
//! calls [`dualis_free`] with the pointer and `4 + length`. Every return is JSON, errors
//! included, so a caller has one shape to parse rather than two.

/// What the page keeps between calls, held here so a megabyte of frames does not cross the
/// boundary once per repaint.
struct State {
    checked: editor_core::Checked,
    run: Option<viewer_core::Run>,
    text: String,
    /// CAD the page has handed over, by the name the scene will use for it. There is no
    /// filesystem here, so this is the only place a part's bytes can live.
    files: editor_core::Uploaded,
}

static mut STATE: Option<State> = None;

/// The state, created on first use.
///
/// `unsafe` because a `static mut` is, and sound because this is compiled for a
/// single-threaded target with no re-entrancy: the page calls one export, it returns, the page
/// calls the next. A `thread_local!` would say the same thing with more ceremony and no more
/// safety on a target that has one thread.
#[allow(static_mut_refs)]
fn state() -> &'static mut State {
    unsafe {
        if STATE.is_none() {
            STATE = Some(State {
                checked: editor_core::Checked::default(),
                run: None,
                text: String::new(),
                files: editor_core::Uploaded::new(),
            });
        }
        STATE.as_mut().expect("just set")
    }
}

/// Give the page a buffer to write a scene into.
#[no_mangle]
pub extern "C" fn dualis_alloc(len: usize) -> *mut u8 {
    let mut buf = Vec::<u8>::with_capacity(len);
    let ptr = buf.as_mut_ptr();
    std::mem::forget(buf);
    ptr
}

/// Take a buffer back.
///
/// # Safety
///
/// `ptr` and `len` must be a pair this module handed out and that has not been freed.
#[no_mangle]
pub unsafe extern "C" fn dualis_free(ptr: *mut u8, len: usize) {
    if !ptr.is_null() && len > 0 {
        drop(Vec::from_raw_parts(ptr, len, len));
    }
}

/// The uploaded set, cloned so a `&mut` state borrow is not held across the call.
fn s_files() -> editor_core::Uploaded {
    state().files.clone()
}

/// Take one uploaded file: a name and its bytes.
///
/// This is the whole reason the browser can do CAD at all. A scene's `parts` name a file, and on
/// a machine that means a path — in a page there is no path and no filesystem, only bytes that
/// arrived because somebody dropped a file on the window. The name is whatever the page calls it,
/// and a scene that says `"stl": "bracket.stl"` finds it if the page uploaded it under that name.
///
/// Replaces silently on a repeat, which is what re-dropping an edited export should do.
///
/// # Safety
///
/// `name` must be valid for `name_len` bytes of UTF-8 and `data` for `data_len` bytes.
#[no_mangle]
pub unsafe extern "C" fn dualis_part(
    name: *const u8,
    name_len: usize,
    data: *const u8,
    data_len: usize,
) -> *mut u8 {
    let label = borrow(name, name_len);
    if label.is_empty() {
        return give(serde_json::json!({ "error": "a part needs a name" }).to_string());
    }
    let bytes = if data.is_null() || data_len == 0 {
        Vec::new()
    } else {
        std::slice::from_raw_parts(data, data_len).to_vec()
    };
    // An empty upload is refused here rather than at build time. A zero-byte file is what a
    // failed read in the page looks like, and `Mesh::from_stl` would report it as a malformed
    // mesh — true, and about the wrong thing.
    if bytes.is_empty() {
        return give(
            serde_json::json!({ "error": format!("{label}: the upload was empty") }).to_string(),
        );
    }
    let s = state();
    s.files.insert(label, bytes);
    give(parts_summary(&s.files))
}

/// Measure a ladder of grids for the uploaded parts.
///
/// The step between "I dropped a CAD file on this" and "here is a scene": a user has no way to
/// know what `cells` and `cell_mm` hold their part, and getting it wrong does not fail — a part
/// finer than the grid rasterises to **nothing** and the run is well behaved about a different
/// object. Every row in the table was rasterised rather than predicted, because
/// `Loss::volume_error` is erratic under refinement and extrapolating it would be confidently
/// wrong.
///
/// Returns `{table, fragment, cell_mm, cells}` or `{error}`. `fragment` is null when no grid holds
/// every part under half its volume in boundary cells, which is an answer and not a failure.
///
/// # Safety
///
/// `material` must be valid for `material_len` bytes of UTF-8.
#[no_mangle]
pub unsafe extern "C" fn dualis_fit(
    budget_cells: u32,
    material: *const u8,
    material_len: usize,
) -> *mut u8 {
    let name = borrow(material, material_len);
    let name = if name.is_empty() {
        "aluminium".to_string()
    } else {
        name
    };
    let budget = if budget_cells == 0 {
        2_000_000
    } else {
        budget_cells as usize
    };
    match editor_core::fit(&s_files(), budget, &name) {
        Ok(json) => give(json),
        Err(e) => give(serde_json::json!({ "error": e }).to_string()),
    }
}

/// The material names a scene can use without declaring them.
///
/// Asked for rather than hardcoded in the page: a list of substances copied into a dropdown is a
/// list that drifts, and the failure is a menu offering something the builder then refuses. This is
/// `Substance::CATALOGUE` itself, so the two cannot disagree.
#[no_mangle]
pub extern "C" fn dualis_materials() -> *mut u8 {
    give(serde_json::json!({ "materials": editor_core::MATERIALS }).to_string())
}

/// Forget one uploaded file, or all of them when given an empty name.
///
/// # Safety
///
/// `name` must be valid for `name_len` bytes of UTF-8.
#[no_mangle]
pub unsafe extern "C" fn dualis_forget_part(name: *const u8, name_len: usize) -> *mut u8 {
    let label = borrow(name, name_len);
    let s = state();
    if label.is_empty() {
        s.files = editor_core::Uploaded::new();
    } else {
        s.files.remove(&label);
    }
    give(parts_summary(&s.files))
}

/// What the page shows in its upload list: the names, in order, and the total size.
fn parts_summary(files: &editor_core::Uploaded) -> String {
    serde_json::json!({
        "error": serde_json::Value::Null,
        "names": files.names(),
        "bytes": files.total_bytes(),
    })
    .to_string()
}

/// Read a caller's buffer as a string, replacing invalid UTF-8 rather than refusing it: a
/// mangled byte in a scene should reach the parser and be reported as bad JSON, not vanish
/// here.
///
/// # Safety
///
/// `ptr` must be valid for `len` bytes.
unsafe fn borrow(ptr: *const u8, len: usize) -> String {
    if ptr.is_null() || len == 0 {
        return String::new();
    }
    String::from_utf8_lossy(std::slice::from_raw_parts(ptr, len)).into_owned()
}

/// Hand a string back, length-prefixed, and forget it — the page frees it.
fn give(s: String) -> *mut u8 {
    let body = s.into_bytes();
    let mut out = Vec::with_capacity(4 + body.len());
    out.extend_from_slice(&(body.len() as u32).to_le_bytes());
    out.extend_from_slice(&body);
    let ptr = out.as_mut_ptr();
    std::mem::forget(out);
    ptr
}

/// Check a scene and keep it: the same parse-then-build the CLI's `--check` runs.
///
/// Returns `{error, summary, notes, boxes, bounds, edges}` — every field present on success and
/// failure alike, because a page redraws from whatever it has and an error must not blank the
/// viewport the last good parse filled.
///
/// # Safety
///
/// `ptr` must be valid for `len` bytes of UTF-8.
#[no_mangle]
pub unsafe extern "C" fn dualis_check(ptr: *const u8, len: usize) -> *mut u8 {
    let text = borrow(ptr, len);
    let checked = editor_core::check(&text, &s_files());
    let boxes: Vec<serde_json::Value> = checked
        .boxes
        .iter()
        .map(|b| serde_json::json!({ "name": b.name, "corners": b.corners }))
        .collect();
    let out = serde_json::json!({
        "error": checked.error,
        "summary": checked.summary,
        "notes": checked.notes,
        "boxes": boxes,
        "bounds": checked.bounds,
    });
    let s = state();
    s.text = text;
    s.checked = checked;
    give(out.to_string())
}

/// Run the scene last handed to [`dualis_check`], and keep the frames.
///
/// Returns `{frames}` or `{error}`. It is the CLI's run — same builder, same audit — so a
/// violation arrives worded by the kernel rather than paraphrased here.
#[no_mangle]
pub extern "C" fn dualis_run() -> *mut u8 {
    let s = state();
    // The text is stored whether or not it checked, because a shell re-checks from it — but a
    // scene that did not check must not be run off the back of that. The page disables the
    // button; a module that relied on the page to remember would run a scene the reader has
    // already been told is broken, and report it in a second, worse wording than the first.
    if let Some(why) = &s.checked.error {
        return give(serde_json::json!({ "error": why }).to_string());
    }
    let out = match editor_core::run(&s.text, &s.files.clone()) {
        Ok(json) => match viewer_core::Run::from_json(&json) {
            Ok(run) => {
                let n = run.frames.len();
                s.run = Some(run);
                serde_json::json!({ "frames": n })
            }
            // This module wrote that JSON one call ago, so failing to read it back is a wire
            // format defect rather than a user's mistake, and it says so.
            Err(e) => {
                serde_json::json!({ "error": format!("the run's own JSON did not read back: {e}") })
            }
        },
        Err(e) => {
            s.run = None;
            serde_json::json!({ "error": e })
        }
    };
    give(out.to_string())
}

/// Run the verification battery on the current scene.
///
/// Returns `{report, findings}` or `{error}` — the same report the CLI prints, because it is
/// the same battery.
#[no_mangle]
pub extern "C" fn dualis_verify(deep: u32) -> *mut u8 {
    let s = state();
    // The same refusal as `dualis_run`, and for the same reason.
    if let Some(why) = &s.checked.error {
        return give(serde_json::json!({ "error": why }).to_string());
    }
    let out = match editor_core::verify(&s.text, deep != 0, &s.files.clone()) {
        Ok((report, findings)) => serde_json::json!({ "report": report, "findings": findings }),
        Err(e) => serde_json::json!({ "error": e }),
    };
    give(out.to_string())
}

/// Project one frame for the page to paint.
///
/// Takes the camera as JSON — `{azimuth, elevation, distance, scale, aspect, frame, fit}` — and
/// returns primitives with every colour resolved: `lines`, `dots`, `labels`, `readings`. **The
/// page paints and decides nothing**, which is what keeps one camera and one colour rule across
/// the native window and the browser.
///
/// `x` and `y` come back in the projection's own normalised units, `-1..1` across the shorter
/// side, so the page maps them exactly as the native shell does: `px = cx + x·w/2`,
/// `py = cy − y·h/2`.
///
/// # Safety
///
/// `ptr` must be valid for `len` bytes of UTF-8 JSON.
#[no_mangle]
pub unsafe extern "C" fn dualis_draw(ptr: *const u8, len: usize) -> *mut u8 {
    let req: serde_json::Value = match serde_json::from_str(&borrow(ptr, len)) {
        Ok(v) => v,
        Err(e) => return give(serde_json::json!({ "error": e.to_string() }).to_string()),
    };
    let num = |k: &str, d: f64| req.get(k).and_then(|v| v.as_f64()).unwrap_or(d);
    let aspect = num("aspect", 1.0).max(1e-6);
    let want_fit = req.get("fit").and_then(|v| v.as_bool()).unwrap_or(false);
    let frame_index = num("frame", 0.0).max(0.0) as usize;

    let s = state();

    // The world this paint frames: scene geometry, widened by whatever the run reached — an
    // orbit's bodies live far outside any placed extent.
    let mut bounds = s.checked.bounds;
    if let Some(run) = &s.run {
        if let Some(f) = run.frames.first() {
            for p in &f.panels {
                bounds = Some(union(bounds, p.bounds()));
            }
        }
    }
    let Some(bounds) = bounds else {
        return give(
            serde_json::json!({
                "empty": "nothing in this scene has geometry — sources, lumps and networks are \
                          readings, not places; run to see what they report"
            })
            .to_string(),
        );
    };

    let framing = viewer_core::Framing::of(bounds);
    let mut camera = viewer_core::Camera {
        azimuth: num("azimuth", 0.7),
        elevation: num("elevation", 0.4),
        distance: num("distance", 2.5),
        scale: num("scale", 1.0),
    };
    if want_fit {
        camera.fit(bounds, &framing, aspect, 0.85);
    }
    let project = |p: [f64; 3]| camera.project(p, &framing, aspect);

    let mut lines: Vec<serde_json::Value> = Vec::new();
    let mut dots: Vec<serde_json::Value> = Vec::new();
    let mut labels: Vec<serde_json::Value> = Vec::new();
    let mut notes: Vec<&str> = Vec::new();
    let mut readings: Vec<serde_json::Value> = Vec::new();
    let mut t = 0.0;

    // The scene's own geometry, drawn from the text rather than the run, so a layout is visible
    // while editing and before anything is computed.
    for b in &s.checked.boxes {
        for (i, j) in editor_core::EDGES {
            let (a, c) = (project(b.corners[i]), project(b.corners[j]));
            lines.push(serde_json::json!([a.x, a.y, c.x, c.y, "#6b7683", 1.0]));
        }
        let at = project(b.corners[0]);
        labels.push(serde_json::json!([at.x, at.y, b.name]));
    }

    if let Some(run) = &s.run {
        let index = frame_index.min(run.frames.len().saturating_sub(1));
        if let Some(frame) = run.frames.get(index) {
            t = frame.t;
            for panel in &frame.panels {
                let scale = run.scale_of(panel.name());
                match panel {
                    viewer_core::Panel::Points {
                        positions, values, ..
                    } => {
                        for (i, v) in values.iter().enumerate() {
                            let q = project([
                                positions[3 * i],
                                positions[3 * i + 1],
                                positions[3 * i + 2],
                            ]);
                            dots.push(serde_json::json!([q.x, q.y, 0.012, shade(*v, scale)]));
                        }
                    }
                    viewer_core::Panel::Paths {
                        starts,
                        vertices,
                        values,
                        ..
                    } => {
                        for (r, v) in values.iter().enumerate() {
                            let lo = starts[r] as usize;
                            let hi = starts
                                .get(r + 1)
                                .map_or(vertices.len() / 3, |x| *x as usize);
                            for w in lo..hi.saturating_sub(1) {
                                let a = project([
                                    vertices[3 * w],
                                    vertices[3 * w + 1],
                                    vertices[3 * w + 2],
                                ]);
                                let c = project([
                                    vertices[3 * w + 3],
                                    vertices[3 * w + 4],
                                    vertices[3 * w + 5],
                                ]);
                                lines.push(serde_json::json!([
                                    a.x,
                                    a.y,
                                    c.x,
                                    c.y,
                                    shade(*v, scale),
                                    1.5
                                ]));
                            }
                        }
                    }
                    viewer_core::Panel::Field {
                        name,
                        unit,
                        nx,
                        ny,
                        nz,
                        values,
                    } => {
                        let Some(b) = s.checked.boxes.iter().find(|b| &b.name == name) else {
                            continue;
                        };
                        let out = editor_core::field_splats(
                            &b.corners,
                            (*nx, *ny, *nz),
                            values,
                            unit,
                            scale,
                        );
                        notes.push(out.note);
                        // The splat's screen size from the box's own projected span, so a coarse
                        // field draws fat cells and a fine one small ones.
                        let (p0, p7) = (project(b.corners[0]), project(b.corners[7]));
                        let span = ((p0.x - p7.x).powi(2) + (p0.y - p7.y).powi(2)).sqrt();
                        let radius = (span / (*nx).max(*ny).max(*nz) as f64 * 0.35
                            / out.stride.max(1) as f64)
                            .clamp(0.002, 0.1);
                        // Painter's algorithm: far to near, so a translucent volume composites
                        // in the right order. `depth` grows away from the eye.
                        let mut splats: Vec<(f64, viewer_core::Projected, [u8; 4])> = out
                            .splats
                            .iter()
                            .map(|sp| {
                                let q = project(sp.at);
                                (q.depth, q, sp.rgba)
                            })
                            .collect();
                        splats.sort_by(|a, b| {
                            b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal)
                        });
                        for (_, q, c) in splats {
                            dots.push(serde_json::json!([
                                q.x,
                                q.y,
                                radius,
                                format!(
                                    "rgba({},{},{},{:.3})",
                                    c[0],
                                    c[1],
                                    c[2],
                                    c[3] as f64 / 255.0
                                )
                            ]));
                        }
                    }
                }
            }
            readings = frame
                .readings
                .iter()
                .map(|r| {
                    serde_json::json!({
                        "domain": r.domain, "label": r.label, "value": r.value, "unit": r.unit
                    })
                })
                .collect();
        }
    }

    give(
        serde_json::json!({
            "lines": lines, "dots": dots, "labels": labels, "notes": notes,
            "readings": readings, "t": t,
            "camera": { "azimuth": camera.azimuth, "elevation": camera.elevation,
                        "distance": camera.distance, "scale": camera.scale },
        })
        .to_string(),
    )
}

/// A value on the run-wide scale as a CSS colour, for the shapes that are not fields.
fn shade(value: f64, scale: Option<(f64, f64)>) -> String {
    let t = match scale {
        Some((lo, hi)) if hi > lo => ((value - lo) / (hi - lo)).clamp(0.0, 1.0),
        _ => 0.5,
    };
    format!(
        "rgb({},{},{})",
        (60.0 + 195.0 * t) as u8,
        (90.0 + 40.0 * (1.0 - (2.0 * t - 1.0).abs())) as u8,
        (230.0 - 170.0 * t) as u8
    )
}

/// Union of an optional box with another box.
fn union(a: Option<[f64; 6]>, b: [f64; 6]) -> [f64; 6] {
    match a {
        None => b,
        Some(a) => [
            a[0].min(b[0]),
            a[1].min(b[1]),
            a[2].min(b[2]),
            a[3].max(b[3]),
            a[4].max(b[4]),
            a[5].max(b[5]),
        ],
    }
}
