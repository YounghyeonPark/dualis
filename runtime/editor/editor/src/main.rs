//! Edit a dualis scene beside a 3D view of it.
//!
//! ```text
//! cargo run --release                 # opens on the built-in room
//! cargo run --release -- scene.json   # opens on a file
//! ```
//!
//! The left pane is the scene's JSON, checked as you type with the same two steps
//! `dualis-world --check` runs; a parse error is shown with its `line:column`. The viewport
//! draws every placed extent as a wireframe — live, from the text, before anything runs — and
//! after **Run** it overlays what the run produced, scrubbable by frame. **Verify** runs the
//! battery from `dualis-world verify` and shows the same report the CLI prints.
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
use std::sync::mpsc;

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
    /// The run's JSON, or the violation that stopped it.
    Ran(Result<String, String>),
    /// The battery's rendered report and its findings count, or why it could not start.
    Verified(Result<(String, usize), String>),
}

/// A parsed run being scrubbed through.
struct RunView {
    run: viewer_core::Run,
    frame: usize,
}

struct App {
    /// The scene text being edited, and where it loads from and saves to.
    text: String,
    path: String,
    /// The result of checking `text`, refreshed whenever the text changes.
    checked: editor_core::Checked,
    /// The last completed run, if any. Cleared when the text changes, because a picture of a
    /// run beside text that no longer produces it is a picture of something else.
    run: Option<RunView>,
    /// The verify report, its findings count, and whether its window is open.
    verify: Option<(String, usize)>,
    verify_open: bool,
    deep: bool,
    /// The in-flight background job, if any. One at a time: a second run pressed mid-run
    /// would race two worlds for one pane.
    busy: Option<(&'static str, mpsc::Receiver<Job>)>,
    status: String,
    /// The viewport camera — `viewer-core`'s, shared arithmetic with the viewer and the HTML
    /// report so all three open on the same picture.
    camera: viewer_core::Camera,
    /// Fit the camera on the next paint. Set when geometry appears or is replaced, not every
    /// frame — a camera that re-fits while you drag is a camera you cannot aim.
    needs_fit: bool,
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
        App {
            text,
            path,
            checked,
            run: None,
            verify: None,
            verify_open: false,
            deep: false,
            busy: None,
            status,
            camera: viewer_core::Camera::default(),
            needs_fit: true,
        }
    }

    fn recheck(&mut self) {
        self.checked = editor_core::check(&self.text);
        self.run = None;
    }

    /// Start a background job, refusing a second while one is in flight.
    fn spawn(&mut self, label: &'static str, job: impl FnOnce() -> Job + Send + 'static) {
        if self.busy.is_some() {
            self.status = format!("still busy with the last job; {label} not started");
            return;
        }
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let _ = tx.send(job());
        });
        self.busy = Some((label, rx));
        self.status = format!("{label}…");
    }

    fn poll_jobs(&mut self) {
        let Some((label, rx)) = &self.busy else {
            return;
        };
        let label = *label;
        match rx.try_recv() {
            Ok(Job::Ran(Ok(json))) => {
                self.busy = None;
                match viewer_core::Run::from_json(&json) {
                    Ok(run) => {
                        let frames = run.frames.len();
                        self.run = Some(RunView {
                            run,
                            frame: frames.saturating_sub(1),
                        });
                        self.needs_fit = true;
                        self.status = format!("ran: {frames} frames");
                    }
                    // The editor wrote this JSON one call ago, so the viewer failing to read
                    // it back is a wire-format defect, not a user mistake — say so.
                    Err(e) => self.status = format!("the run's own JSON did not read back: {e}"),
                }
            }
            Ok(Job::Ran(Err(e))) => {
                self.busy = None;
                self.status = e;
            }
            Ok(Job::Verified(Ok((report, findings)))) => {
                self.busy = None;
                self.verify = Some((report, findings));
                self.verify_open = true;
                self.status = match findings {
                    0 => String::from("verified: no structural findings"),
                    n => format!("verify: {n} finding(s) — see the report"),
                };
            }
            Ok(Job::Verified(Err(e))) => {
                self.busy = None;
                self.status = e;
            }
            Err(mpsc::TryRecvError::Empty) => {
                self.status = format!("{label}…");
            }
            Err(mpsc::TryRecvError::Disconnected) => {
                self.busy = None;
                self.status =
                    format!("{label} died without reporting — that is a bug worth a note");
            }
        }
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_jobs();
        if self.busy.is_some() {
            // A background job finishes whether or not the mouse moves; without this the
            // result waits for the next wiggle.
            ctx.request_repaint_after(std::time::Duration::from_millis(100));
        }

        egui::TopBottomPanel::top("bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("file");
                ui.add(egui::TextEdit::singleline(&mut self.path).desired_width(260.0));
                if ui.button("load").clicked() {
                    match std::fs::read_to_string(&self.path) {
                        Ok(t) => {
                            self.text = t;
                            self.recheck();
                            self.needs_fit = true;
                            self.status = format!("loaded {}", self.path);
                        }
                        Err(e) => self.status = format!("{}: {e}", self.path),
                    }
                }
                if ui.button("save").clicked() {
                    match std::fs::write(&self.path, &self.text) {
                        Ok(()) => self.status = format!("saved {}", self.path),
                        Err(e) => self.status = format!("{}: {e}", self.path),
                    }
                }
                ui.separator();
                let runnable = self.checked.error.is_none() && self.busy.is_none();
                if ui.add_enabled(runnable, egui::Button::new("run")).clicked() {
                    let text = self.text.clone();
                    self.spawn("running", move || Job::Ran(editor_core::run(&text)));
                }
                if ui
                    .add_enabled(runnable, egui::Button::new("verify"))
                    .clicked()
                {
                    let text = self.text.clone();
                    let deep = self.deep;
                    self.spawn("verifying", move || {
                        Job::Verified(editor_core::verify(&text, deep))
                    });
                }
                ui.checkbox(&mut self.deep, "deep");
                if ui.button("fit view").clicked() {
                    self.needs_fit = true;
                }
            });
        });

        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            let error = self.checked.error.as_deref();
            match error {
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
        // viewer-core states. A field's box is already on screen; its values stay the HTML
        // report's job for now, and saying so beats quietly drawing nothing.
        let Some(view) = &mut self.run else {
            return;
        };
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
        let frame = &view.run.frames[view.frame];
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
                viewer_core::Panel::Field { name, .. } => {
                    // Its extent is already wireframed above. Not silently: the box carries a
                    // note so a field is visibly "not drawn here" rather than absent.
                    if let Some(b) = self.checked.boxes.iter().find(|b| &b.name == name) {
                        painter.text(
                            to_screen(b.corners[0]),
                            egui::Align2::LEFT_TOP,
                            "field: values in the HTML report",
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
