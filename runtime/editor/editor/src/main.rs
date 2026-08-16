//! Edit a dualis scene beside a 3D view of it — and leave it open while a script does the
//! editing.
//!
//! ```text
//! cargo run --release                 # opens on the built-in room
//! cargo run --release -- scene.json   # opens on a file
//! ```
//!
//! The left pane is the scene's JSON, checked as you type with the same two steps
//! `dualis-world --check` runs; a parse error is shown with its `line:column`. The viewport
//! draws every placed extent as a wireframe — live, from the text, before anything runs — and
//! **Run** streams the run in as it computes: each frame appears when it is captured, the
//! slider grows, and **stop** ends a long run between frames. **Verify** runs the battery
//! from `dualis-world verify` and shows the same report the CLI prints.
//!
//! # The live loop
//!
//! With **watch file** on (it is on by default), the editor polls the file's modified time
//! and reloads when something else writes it — a script, an agent, another editor. With
//! **run on change** on as well, the reload runs the scene, so the loop
//! `script writes → editor rechecks → runs → draws` closes with no hand on the window.
//! An in-flight run is stopped and superseded when the file changes again, so the picture
//! converges on the latest text rather than queueing history.
//!
//! One rule keeps that honest: **unsaved edits in the pane are never clobbered.** If the pane
//! is dirty and the disk changes, the status line says so, loudly, and the disk's version
//! waits for an explicit `load`. The alternative silently discards somebody's typing, and
//! which of the two writers meant it is not the editor's call to make.
//!
//! Drag to rotate, scroll to zoom, and the camera is `viewer-core`'s — the same fit, the same
//! projection, the same clamps, because that arithmetic has been wrong here before and is not
//! being written a third time.
//!
//! # The two halves, kept apart
//!
//! ARCHITECTURE.md's platform rules: the authoring half is the composition root and may name
//! domains; the inspection half dispatches on the shape of the data. This file is the shell
//! over both and holds neither — authoring machinery lives in `editor-core`, and everything
//! painted below is painted from a *shape* (a box, points, paths, a reading) with the domain
//! name used only as a label. A domain added next year is drawn by the code below unchanged.

use eframe::egui;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{mpsc, Arc};
use std::time::{Duration, Instant, SystemTime};

fn main() -> eframe::Result {
    let path = std::env::args().nth(1);
    eframe::run_native(
        "dualis editor",
        eframe::NativeOptions::default(),
        Box::new(|_cc| Ok(Box::new(App::new(path)))),
    )
}

/// What a background job sent back.
enum Job {
    /// The run so far, as JSON — one of these per captured frame, plus the settled final.
    Frames(String),
    /// How the streaming run ended: finished, stopped, or the violation that refused it.
    RunEnded(Result<editor_core::RunEnd, String>),
    /// The battery's rendered report and its findings count, or why it could not start.
    Verified(Result<(String, usize), String>),
}

/// A parsed run being scrubbed through.
struct RunView {
    run: viewer_core::Run,
    frame: usize,
    /// Whether what is on screen is a prefix of a stopped run rather than the run. Set from
    /// [`editor_core::RunEnd`], drawn on the canvas, because a prefix that looks complete is
    /// a picture of something that did not happen.
    partial: bool,
}

struct App {
    /// The scene text being edited, and where it loads from and saves to.
    text: String,
    path: String,
    /// The result of checking `text`, refreshed whenever the text changes.
    checked: editor_core::Checked,
    /// The last run — possibly still growing, streamed frame by frame.
    run: Option<RunView>,
    /// The verify report, its findings count, and whether its window is open.
    verify: Option<(String, usize)>,
    verify_open: bool,
    deep: bool,
    /// The in-flight background job, if any. One at a time: a second run pressed mid-run
    /// would race two worlds for one pane.
    busy: Option<(&'static str, mpsc::Receiver<Job>)>,
    /// Raised to end a streaming run between frames; replaced on every spawn.
    stop: Arc<AtomicBool>,
    status: String,
    /// The viewport camera — `viewer-core`'s, shared arithmetic with the viewer and the HTML
    /// report so all three open on the same picture.
    camera: viewer_core::Camera,
    /// Fit the camera on the next paint. Set when geometry appears or is replaced, not every
    /// frame — a camera that re-fits while you drag is a camera you cannot aim.
    needs_fit: bool,

    // The live loop.
    /// Poll the file for outside writes.
    watch: bool,
    /// Run automatically after an outside write is loaded.
    auto_run: bool,
    /// The pane has edits the disk does not. While true, outside writes are announced and
    /// never applied.
    dirty: bool,
    /// The modified time last loaded or saved, so an outside write is a *different* mtime
    /// rather than any mtime.
    known_mtime: Option<SystemTime>,
    /// When the file was last polled; polling is cheap but not free, and sixty times a
    /// second buys nothing over twice a second.
    last_poll: Option<Instant>,
    /// An auto-run is owed as soon as the current job ends — set when the file changed while
    /// something was in flight, so the picture converges on the latest text.
    rerun_owed: bool,
}

impl App {
    fn new(path: Option<String>) -> App {
        let (text, path, status) = match path {
            Some(p) => match std::fs::read_to_string(&p) {
                Ok(t) => (t, p.clone(), format!("loaded {p}")),
                Err(e) => (
                    default_scene(),
                    p.clone(),
                    format!("{p}: {e}; opened the built-in scene"),
                ),
            },
            None => (
                default_scene(),
                String::from("scene.json"),
                String::from("the built-in scene"),
            ),
        };
        let checked = editor_core::check(&text);
        let known_mtime = mtime_of(&path);
        App {
            text,
            path,
            checked,
            run: None,
            verify: None,
            verify_open: false,
            deep: false,
            busy: None,
            stop: Arc::new(AtomicBool::new(false)),
            status,
            camera: viewer_core::Camera::default(),
            needs_fit: true,
            watch: true,
            auto_run: false,
            dirty: false,
            known_mtime,
            last_poll: None,
            rerun_owed: false,
        }
    }

    fn recheck(&mut self) {
        self.checked = editor_core::check(&self.text);
        self.run = None;
    }

    /// Start a background job, refusing a second while one is in flight.
    fn spawn(&mut self, label: &'static str, job: impl FnOnce(mpsc::Sender<Job>) + Send + 'static) {
        if self.busy.is_some() {
            self.status = format!("still busy with the last job; {label} not started");
            return;
        }
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || job(tx));
        self.busy = Some((label, rx));
        self.status = format!("{label}…");
    }

    fn start_run(&mut self) {
        let text = self.text.clone();
        self.stop = Arc::new(AtomicBool::new(false));
        let stop = self.stop.clone();
        self.spawn("running", move |tx| {
            let end = editor_core::run_streaming(&text, &stop, |json| {
                let _ = tx.send(Job::Frames(json));
            });
            let _ = tx.send(Job::RunEnded(end));
        });
    }

    fn start_verify(&mut self) {
        let text = self.text.clone();
        let deep = self.deep;
        self.spawn("verifying", move |tx| {
            let _ = tx.send(Job::Verified(editor_core::verify(&text, deep)));
        });
    }

    /// Drain everything the background job has sent. Drained, not sampled: a fast run can
    /// produce several frames per paint, and showing only the first of them would make the
    /// stream look slower than the simulation.
    fn poll_jobs(&mut self) {
        let Some((_, rx)) = &self.busy else { return };
        let mut ended = None;
        let mut latest_frames = None;
        while let Ok(job) = rx.try_recv() {
            match job {
                Job::Frames(json) => latest_frames = Some(json),
                other => {
                    ended = Some(other);
                    break;
                }
            }
        }
        if let Some(json) = latest_frames {
            match viewer_core::Run::from_json(&json) {
                Ok(run) => {
                    let last = run.frames.len().saturating_sub(1);
                    // Follow the tail while the run grows, unless the person has scrubbed
                    // back — a slider that snatches itself out of a hand is worse than one
                    // that lags.
                    let follow = self
                        .run
                        .as_ref()
                        .is_none_or(|v| v.frame + 1 >= v.run.frames.len());
                    let frame = if follow {
                        last
                    } else {
                        self.run.as_ref().map_or(last, |v| v.frame.min(last))
                    };
                    if self.run.is_none() {
                        self.needs_fit = true;
                    }
                    self.run = Some(RunView {
                        run,
                        frame,
                        partial: true,
                    });
                    self.status = format!("running: {} frame(s) so far", last + 1);
                }
                // The editor wrote this JSON one call ago, so the viewer failing to read it
                // back is a wire-format defect, not a user mistake — say so.
                Err(e) => self.status = format!("the run's own JSON did not read back: {e}"),
            }
        }
        match ended {
            None => {}
            Some(Job::RunEnded(end)) => {
                self.busy = None;
                match end {
                    Ok(editor_core::RunEnd::Finished) => {
                        if let Some(v) = &mut self.run {
                            v.partial = false;
                            self.status = format!("ran: {} frames", v.run.frames.len());
                        }
                    }
                    Ok(editor_core::RunEnd::Stopped) => {
                        self.status =
                            String::from("stopped — what is on screen is a prefix, not the run");
                    }
                    Err(e) => {
                        // The frames already streamed stay on screen beside the reason the
                        // run ended, which is the view somebody debugging a violation wants.
                        self.status = e;
                    }
                }
                if std::mem::take(&mut self.rerun_owed) {
                    self.start_run();
                }
            }
            Some(Job::Verified(result)) => {
                self.busy = None;
                match result {
                    Ok((report, findings)) => {
                        self.verify = Some((report, findings));
                        self.verify_open = true;
                        self.status = match findings {
                            0 => String::from("verified: no structural findings"),
                            n => format!("verify: {n} finding(s) — see the report"),
                        };
                    }
                    Err(e) => self.status = e,
                }
                if std::mem::take(&mut self.rerun_owed) {
                    self.start_run();
                }
            }
            Some(Job::Frames(_)) => unreachable!("frames are drained above"),
        }
    }

    /// Notice an outside write, and apply it only when the pane has nothing to lose.
    fn poll_disk(&mut self) {
        if !self.watch {
            return;
        }
        let due = self
            .last_poll
            .is_none_or(|t| t.elapsed() > Duration::from_millis(400));
        if !due {
            return;
        }
        self.last_poll = Some(Instant::now());
        let Some(disk) = mtime_of(&self.path) else {
            return;
        };
        if Some(disk) == self.known_mtime {
            return;
        }
        if self.dirty {
            // Announced, sticky, and not applied. Which of two writers meant it is not the
            // editor's call; the disk's version waits for an explicit `load`.
            self.status = format!(
                "{} changed on disk while the pane has unsaved edits — press load to take \
                 the disk's version",
                self.path
            );
            self.known_mtime = Some(disk);
            return;
        }
        match std::fs::read_to_string(&self.path) {
            Ok(t) => {
                self.known_mtime = Some(disk);
                self.text = t;
                self.recheck();
                self.status = format!("reloaded {} (changed outside)", self.path);
                if self.auto_run && self.checked.error.is_none() {
                    if self.busy.is_some() {
                        // Converge on the latest text: end the in-flight run at the next
                        // frame boundary and owe a fresh one.
                        self.stop.store(true, Ordering::Relaxed);
                        self.rerun_owed = true;
                    } else {
                        self.start_run();
                    }
                }
            }
            Err(e) => self.status = format!("{}: {e}", self.path),
        }
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_jobs();
        self.poll_disk();
        if self.busy.is_some() {
            // A stream draws itself; without this, frames wait for a mouse wiggle.
            ctx.request_repaint_after(Duration::from_millis(60));
        } else if self.watch {
            ctx.request_repaint_after(Duration::from_millis(400));
        }

        egui::TopBottomPanel::top("bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("file");
                ui.add(egui::TextEdit::singleline(&mut self.path).desired_width(220.0));
                if ui.button("load").clicked() {
                    match std::fs::read_to_string(&self.path) {
                        Ok(t) => {
                            self.text = t;
                            self.dirty = false;
                            self.known_mtime = mtime_of(&self.path);
                            self.recheck();
                            self.needs_fit = true;
                            self.status = format!("loaded {}", self.path);
                        }
                        Err(e) => self.status = format!("{}: {e}", self.path),
                    }
                }
                if ui.button("save").clicked() {
                    match std::fs::write(&self.path, &self.text) {
                        Ok(()) => {
                            self.dirty = false;
                            self.known_mtime = mtime_of(&self.path);
                            self.status = format!("saved {}", self.path);
                        }
                        Err(e) => self.status = format!("{}: {e}", self.path),
                    }
                }
                ui.separator();
                let runnable = self.checked.error.is_none() && self.busy.is_none();
                if ui.add_enabled(runnable, egui::Button::new("run")).clicked() {
                    self.start_run();
                }
                if let Some(("running", _)) = self.busy {
                    if ui.button("stop").clicked() {
                        self.stop.store(true, Ordering::Relaxed);
                    }
                }
                if ui
                    .add_enabled(runnable, egui::Button::new("verify"))
                    .clicked()
                {
                    self.start_verify();
                }
                ui.checkbox(&mut self.deep, "deep");
                ui.separator();
                ui.checkbox(&mut self.watch, "watch file");
                ui.checkbox(&mut self.auto_run, "run on change");
                if ui.button("fit view").clicked() {
                    self.needs_fit = true;
                }
            });
        });

        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            match self.checked.error.as_deref() {
                Some(e) => {
                    ui.colored_label(egui::Color32::from_rgb(220, 80, 80), e);
                }
                None => {
                    ui.label(&self.status);
                }
            }
        });

        egui::SidePanel::left("text")
            .resizable(true)
            .default_width(430.0)
            .show(ctx, |ui| {
                if let Some(summary) = &self.checked.summary {
                    ui.label(summary.clone());
                    for note in &self.checked.notes {
                        ui.colored_label(
                            egui::Color32::from_rgb(230, 180, 60),
                            format!("note: {note}"),
                        );
                    }
                    ui.separator();
                }
                egui::ScrollArea::vertical().show(ui, |ui| {
                    let edit = ui.add(
                        egui::TextEdit::multiline(&mut self.text)
                            .code_editor()
                            .desired_width(f32::INFINITY)
                            .desired_rows(34),
                    );
                    if edit.changed() {
                        self.dirty = true;
                        self.recheck();
                    }
                });
            });

        if self.verify_open {
            let mut open = true;
            if let Some((report, findings)) = &self.verify {
                let title = match findings {
                    0 => "verify — no structural findings".to_string(),
                    n => format!("verify — {n} FINDING(S)"),
                };
                egui::Window::new(title).open(&mut open).show(ctx, |ui| {
                    egui::ScrollArea::both().show(ui, |ui| {
                        ui.label(egui::RichText::new(report.clone()).monospace());
                    });
                });
            }
            self.verify_open = open;
        }

        egui::CentralPanel::default().show(ctx, |ui| {
            self.viewport(ui);
        });
    }
}

impl App {
    /// The 3D viewport: wireframe extents from the text, run panels over them by shape.
    fn viewport(&mut self, ui: &mut egui::Ui) {
        let (response, painter) =
            ui.allocate_painter(ui.available_size(), egui::Sense::click_and_drag());
        let rect = response.rect;
        let aspect = (rect.width() / rect.height().max(1.0)) as f64;

        if response.dragged() {
            let d = response.drag_delta();
            self.camera.turn(d.x as f64 * 0.01, d.y as f64 * 0.01);
        }
        if response.hovered() {
            let scroll = ui.input(|i| i.raw_scroll_delta.y) as f64;
            if scroll != 0.0 {
                self.camera.zoom(0.999_f64.powf(scroll * 3.0));
            }
        }

        // The world this paint frames: scene geometry, widened by whatever the run reached —
        // an orbit's bodies live far outside any placed extent.
        let mut bounds = self.checked.bounds;
        if let Some(view) = &self.run {
            if let Some(frame) = view.run.frames.first() {
                for panel in &frame.panels {
                    bounds = Some(union(bounds, panel.bounds()));
                }
            }
        }
        let Some(bounds) = bounds else {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "nothing in this scene has geometry — sources, lumps and networks are \
                 readings, not places; run to see what they report",
                egui::FontId::proportional(14.0),
                ui.visuals().weak_text_color(),
            );
            return;
        };
        let framing = viewer_core::Framing::of(bounds);
        if self.needs_fit {
            self.camera.fit(bounds, &framing, aspect, 0.85);
            self.needs_fit = false;
        }

        let to_screen = |p: [f64; 3]| -> egui::Pos2 {
            let q = self.camera.project(p, &framing, aspect);
            egui::pos2(
                rect.center().x + (q.x as f32) * rect.width() * 0.5,
                rect.center().y - (q.y as f32) * rect.height() * 0.5,
            )
        };

        // The scene's own geometry: every placed extent, wireframed. Drawn from the text, not
        // the run, so layout is visible while editing and before anything is computed.
        let wire = ui.visuals().weak_text_color();
        for placed in &self.checked.boxes {
            for (a, b) in editor_core::EDGES {
                painter.line_segment(
                    [to_screen(placed.corners[a]), to_screen(placed.corners[b])],
                    egui::Stroke::new(1.0, wire),
                );
            }
            painter.text(
                to_screen(placed.corners[0]),
                egui::Align2::LEFT_BOTTOM,
                &placed.name,
                egui::FontId::proportional(12.0),
                ui.visuals().text_color(),
            );
        }

        // The run, by shape. Points are circles, paths are polylines, both coloured by value
        // on the run-wide scale — one scale across the run, never per frame, for the reason
        // viewer-core states; while a run streams, "the run" is the run so far, and the
        // colours settle when it does. A field's box is already on screen; its values stay
        // the HTML report's job for now, and saying so beats quietly drawing nothing.
        let Some(view) = &mut self.run else {
            return;
        };
        if view.partial {
            painter.text(
                egui::pos2(rect.left() + 8.0, rect.top() + 8.0),
                egui::Align2::LEFT_TOP,
                "streaming — a prefix of the run",
                egui::FontId::proportional(12.0),
                egui::Color32::from_rgb(230, 180, 60),
            );
        }
        let frames = view.run.frames.len();
        if frames > 1 {
            ui.scope_builder(egui::UiBuilder::new().max_rect(rect.shrink(8.0)), |ui| {
                ui.with_layout(egui::Layout::bottom_up(egui::Align::LEFT), |ui| {
                    ui.horizontal(|ui| {
                        ui.add(egui::Slider::new(&mut view.frame, 0..=frames - 1).text("frame"));
                        ui.label(format!("t = {:.4} s", view.run.frames[view.frame].t));
                    });
                });
            });
        }
        let frame = &view.run.frames[view.frame.min(frames - 1)];
        for panel in &frame.panels {
            let scale = view.run.scale_of(panel.name());
            match panel {
                viewer_core::Panel::Points {
                    positions, values, ..
                } => {
                    for (i, value) in values.iter().enumerate() {
                        let p = [positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]];
                        painter.circle_filled(to_screen(p), 3.5, shade(*value, scale));
                    }
                }
                viewer_core::Panel::Paths {
                    starts,
                    vertices,
                    values,
                    ..
                } => {
                    for (r, value) in values.iter().enumerate() {
                        let lo = starts[r] as usize;
                        let hi = starts
                            .get(r + 1)
                            .map_or(vertices.len() / 3, |s| *s as usize);
                        for w in lo..hi.saturating_sub(1) {
                            let a = [vertices[3 * w], vertices[3 * w + 1], vertices[3 * w + 2]];
                            let b = [
                                vertices[3 * w + 3],
                                vertices[3 * w + 4],
                                vertices[3 * w + 5],
                            ];
                            painter.line_segment(
                                [to_screen(a), to_screen(b)],
                                egui::Stroke::new(1.5, shade(*value, scale)),
                            );
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
                    let Some(b) = self.checked.boxes.iter().find(|b| &b.name == name) else {
                        // A field with no placed extent has nowhere to be drawn, and saying so
                        // beats an absence: the scene stated a field and the picture has none.
                        continue;
                    };
                    if let Some(note) = draw_field(
                        &painter,
                        &to_screen,
                        b,
                        (*nx, *ny, *nz),
                        values,
                        unit,
                        scale,
                    ) {
                        painter.text(
                            to_screen(b.corners[0]),
                            egui::Align2::LEFT_TOP,
                            note,
                            egui::FontId::proportional(10.0),
                            ui.visuals().weak_text_color(),
                        );
                    }
                }
            }
        }

        // The frame's readings, top-right: the numbers for everything that has no picture.
        let readings = &frame.readings;
        if !readings.is_empty() {
            let mut y = rect.top() + 8.0;
            for r in readings {
                painter.text(
                    egui::pos2(rect.right() - 8.0, y),
                    egui::Align2::RIGHT_TOP,
                    format!("{} {} {:.4} {}", r.domain, r.label, r.value, r.unit),
                    egui::FontId::monospace(11.0),
                    ui.visuals().text_color(),
                );
                y += 14.0;
            }
        }
    }
}

/// The file's modified time, or `None` for a path that does not resolve to one.
fn mtime_of(path: &str) -> Option<SystemTime> {
    std::fs::metadata(path).and_then(|m| m.modified()).ok()
}

/// The most cells one field may draw in a frame.
///
/// A 100³ field is a million splats and a painter that stops painting. Past this the field is
/// **subsampled with a stride**, and the canvas says so — a silently decimated picture is a
/// picture of a coarser simulation than the one that ran, which is exactly the kind of quiet
/// substitution this workspace refuses elsewhere.
const MAX_VOXELS: usize = 8000;

/// Draw a field as depth-sorted voxels, and return the note the canvas should carry.
///
/// # Two colourings, and the physics decides which
///
/// A temperature field is drawn in the colour a body at that temperature **actually is** —
/// Planck's law through the CIE matching functions, from `dualis::view::colour` — whenever
/// anything in the field is hot enough to emit visible light at all. That is not a palette: a
/// melting block glows the orange a melting block glows, and nothing here picked it.
///
/// Below that, physics gives no colour — a body at 300 K emits nothing visible, and this
/// workspace holds no visible reflectance for it to have instead — so the field falls back to
/// the conventional ramp, which says *more* and *less* and does not pretend to say *looks
/// like*. The returned note states which of the two a reader is looking at, because a false
/// colour mistaken for a real one is a wrong answer that looks right.
fn draw_field(
    painter: &egui::Painter,
    to_screen: &impl Fn([f64; 3]) -> egui::Pos2,
    placed: &editor_core::PlacedBox,
    counts: (usize, usize, usize),
    values: &[f64],
    unit: &str,
    scale: Option<(f64, f64)>,
) -> Option<&'static str> {
    let (nx, ny, nz) = counts;
    if nx == 0 || ny == 0 || nz == 0 || values.len() < nx * ny * nz {
        return Some("field: the panel's values do not fill its grid — not drawn");
    }

    // The extent's axes, from the corners `editor-core` posed. Corner 0 is the low one and
    // bits 0, 1, 2 step one axis each, which is the order `EDGES` is written against.
    let o = placed.corners[0];
    let axis = |c: [f64; 3]| [c[0] - o[0], c[1] - o[1], c[2] - o[2]];
    let (ax, ay, az) = (
        axis(placed.corners[1]),
        axis(placed.corners[2]),
        axis(placed.corners[4]),
    );

    // A temperature in kelvin, or nothing. The unit is the panel's own — data, not a domain
    // name — which is the rule the view layer is held to and the reason this works for a
    // physics that does not exist yet.
    let to_kelvin = |v: f64| match unit {
        "K" => Some(v),
        "C" => Some(v + 273.15),
        _ => None,
    };
    // Does anything in this field glow? `glow_fraction` answers it with the visible share of
    // the radiated power rather than a threshold somebody chose.
    let hottest = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let peak_glow = to_kelvin(hottest).map_or(0.0, dualis::view::glow_fraction);
    let physical = peak_glow > 1e-6;

    // Subsample to a budget, and report it rather than quietly drawing a coarser field.
    let total = nx * ny * nz;
    let stride = if total > MAX_VOXELS {
        ((total as f64 / MAX_VOXELS as f64).cbrt().ceil() as usize).max(1)
    } else {
        1
    };

    // One splat's screen size: the box's own projected span divided by its grid, so a coarse
    // field draws fat cells and a fine one draws small ones instead of both drawing dots.
    let span_px = {
        let a = to_screen(placed.corners[0]);
        let b = to_screen(placed.corners[7]);
        ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt()
    };
    let radius =
        (span_px / (nx.max(ny).max(nz) as f32) * 0.7 / stride.max(1) as f32).clamp(1.0, 40.0);

    let frac = |i: usize, n: usize| {
        if n > 1 {
            i as f64 / (n - 1) as f64
        } else {
            0.5
        }
    };

    // Collected then depth-sorted: painter's algorithm, far to near, which is what the viewer
    // uses and what makes a translucent volume composite in the right order.
    let mut splats: Vec<(f64, egui::Pos2, egui::Color32)> = Vec::new();
    for k in (0..nz).step_by(stride) {
        for j in (0..ny).step_by(stride) {
            for i in (0..nx).step_by(stride) {
                let v = values[i + nx * (j + ny * k)];
                if !v.is_finite() {
                    continue;
                }
                let (u, w, t) = (frac(i, nx), frac(j, ny), frac(k, nz));
                let p = [
                    o[0] + ax[0] * u + ay[0] * w + az[0] * t,
                    o[1] + ax[1] * u + ay[1] * w + az[1] * t,
                    o[2] + ax[2] * u + ay[2] * w + az[2] * t,
                ];
                // Normalised position on the run-wide scale — one scale across the whole run,
                // never per frame, for the reason `viewer-core` states.
                let s = match scale {
                    Some((lo, hi)) if hi > lo => ((v - lo) / (hi - lo)).clamp(0.0, 1.0),
                    _ => 0.5,
                };
                let colour = if physical {
                    let kelvin = to_kelvin(v).unwrap_or(0.0);
                    let [r, g, b] = dualis::view::blackbody_srgb(kelvin);
                    // Brightness is the glow relative to the field's own hottest cell, so a
                    // cool corner of a glowing block is dark rather than merely bluer — which
                    // is what a photograph of it looks like.
                    let rel = (dualis::view::glow_fraction(kelvin) / peak_glow).clamp(0.0, 1.0);
                    let a = (rel.sqrt() * 235.0) as u8;
                    egui::Color32::from_rgba_unmultiplied(r, g, b, a.max(6))
                } else {
                    let c = shade(v, scale);
                    let a = (30.0 + 200.0 * s * s) as u8;
                    egui::Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), a)
                };
                let screen = to_screen(p);
                // Depth from the same projection the camera uses, so the sort agrees with the
                // picture rather than approximating it.
                splats.push((depth_of(to_screen, p), screen, colour));
            }
        }
    }
    splats.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    for (_, at, colour) in &splats {
        painter.circle_filled(*at, radius, *colour);
    }

    match (physical, stride) {
        (true, 1) => Some("field: colour is Planck's, not a palette"),
        (true, _) => Some("field: Planck colour, subsampled — see the report for every cell"),
        (false, 1) => Some("field: false colour — nothing here is hot enough to glow"),
        (false, _) => Some("field: false colour, subsampled — see the report for every cell"),
    }
}

/// A stand-in depth for the painter's sort: the screen position alone cannot order two points,
/// so the world point's distance along the view is recovered by projecting a point pushed
/// slightly away — cheaper than threading the camera through, and monotone in true depth,
/// which is all a sort needs.
fn depth_of(to_screen: &impl Fn([f64; 3]) -> egui::Pos2, p: [f64; 3]) -> f64 {
    // The projection is a perspective divide, so a fixed world offset shrinks with distance:
    // the screen separation between `p` and `p + eps` falls monotonically as `p` recedes.
    let a = to_screen(p);
    let b = to_screen([p[0], p[1], p[2] + 1e-3]);
    let sep = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt() as f64;
    // Larger separation means nearer, so invert for a far-to-near sort key.
    if sep > 0.0 {
        1.0 / sep
    } else {
        f64::MAX
    }
}

/// A value on the run-wide scale, as a colour. The two-ended ramp the HTML report uses in
/// spirit: low is blue, high is warm, and a missing scale (a panel of one constant) is drawn
/// mid-ramp rather than invisibly.
fn shade(value: f64, scale: Option<(f64, f64)>) -> egui::Color32 {
    let t = match scale {
        Some((lo, hi)) if hi > lo => ((value - lo) / (hi - lo)).clamp(0.0, 1.0),
        _ => 0.5,
    } as f32;
    egui::Color32::from_rgb(
        (60.0 + 195.0 * t) as u8,
        (90.0 + 40.0 * (1.0 - (2.0 * t - 1.0).abs())) as u8,
        (230.0 - 170.0 * t) as u8,
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

/// The scene the editor opens on with no file: the same built-in room `dualis-world` runs
/// with no arguments, so the two front ends agree about where "hello" is.
fn default_scene() -> String {
    String::from(
        r#"{
  "title": "a small room ringing in its (1,1) mode",
  "schedule": "multirate",
  "duration_s": 0.02,
  "frames": 11,
  "conservation_tolerance": 1e-6,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.4, "height_m": 3.1,
      "cells_across": 61,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } }
  ]
}
"#,
    )
}
