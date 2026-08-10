//! A native window for a dualis run.
//!
//! ```text
//! cargo run --release -- run.json
//! ```
//!
//! Drag to rotate, scroll to zoom, space to play, left and right to scrub.
//!
//! ```text
//! cargo run --release -- run.json --snapshot out.ppm
//! ```
//!
//! renders one frame **without opening a window** and writes it out. That mode is not a
//! convenience: it is the only way anything can check that this file draws. A window nobody can
//! photograph proves it did not panic, which is a much weaker claim than it looks — every wrong
//! projection, every empty vertex buffer and every silently-failed pipeline also does not panic.
//!
//! # What is here and what is deliberately not
//!
//! This is the shell. Everything that could be got wrong twice — the colour scale across a run,
//! the framing, the projection — is in `viewer-core`, which has no GPU dependency and is tested
//! against real run files. What is left here is a surface, a line pipeline and an event loop.
//!
//! That split is not tidiness. A renderer is the one place where a wrong answer looks like a
//! picture, so the arithmetic lives where a test can reach it and this file only draws what it is
//! handed.
//!
//! # It does not depend on `dualis`
//!
//! Deliberately. It reads the JSON a run wrote and nothing else, so "the wire format carries
//! enough to draw a run" is demonstrated rather than asserted. If this needed to link the library
//! for something the file did not have, the format would be the thing to fix.

use std::sync::Arc;

use viewer_core::{segments, Camera, Framing, Panel, Run};
use winit::application::ApplicationHandler;
use winit::event::{ElementState, MouseScrollDelta, WindowEvent};
use winit::event_loop::{ActiveEventLoop, ControlFlow, EventLoop};
use winit::keyboard::{KeyCode, PhysicalKey};
use winit::window::{Window, WindowId};

/// One end of a line, as the shader wants it.
#[repr(C)]
#[derive(Clone, Copy)]
struct Vertex {
    position: [f32; 2],
    colour: [f32; 3],
}

/// A stand-in for `bytemuck`, which is one dependency for one trait.
mod bytemuck_lite {
    /// Marker for a type that can be reinterpreted as bytes.
    ///
    /// # Safety
    ///
    /// Implementors must be `#[repr(C)]`, contain no padding and no references. Both types here
    /// are arrays of `f32`, which satisfies all three.
    pub unsafe trait Pod: Copy {}

    /// The bytes of a slice.
    pub fn cast_slice<T: Pod>(v: &[T]) -> &[u8] {
        // SAFETY: `T: Pod` promises a plain-data layout, and the length is derived from the
        // slice's own, so the range is exactly the memory the slice owns.
        unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
    }
}

// SAFETY: `#[repr(C)]`, five `f32` with no padding and no references.
unsafe impl bytemuck_lite::Pod for Vertex {}

fn main() {
    let path = match std::env::args().nth(1) {
        Some(p) => p,
        None => {
            eprintln!("usage: viewer <run.json>");
            eprintln!("  produced by `dualis-world <scene> out.json`, or by any run that calls");
            eprintln!("  dualis_view::to_json");
            std::process::exit(2);
        }
    };
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("cannot read {path}: {e}");
            std::process::exit(1);
        }
    };
    let run = match Run::from_json(&text) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    };

    // Which panel to show. The first one with paths, because that is what this shell draws; a
    // field or a point cloud needs a different pipeline and saying so beats opening a blank
    // window and letting the reader wonder.
    let drawable = run.panels().into_iter().find(|name| {
        run.frames
            .iter()
            .flat_map(|f| &f.panels)
            .any(|p| p.name() == name && matches!(p, Panel::Paths { .. }))
    });
    let Some(panel) = drawable else {
        eprintln!(
            "{}: no panel of paths to draw. This shell draws paths; fields and point clouds are \
             different pipelines and are not built yet.",
            run.title
        );
        eprintln!("  panels in this run: {:?}", run.panels());
        std::process::exit(1);
    };

    println!(
        "{} — {} frames, showing {panel}",
        run.title,
        run.frames.len()
    );
    println!("  drag to rotate, scroll to zoom, space to play, left/right to scrub");

    // Headless: render one frame to a texture, read it back, write a PPM. No window, no display,
    // so a machine with neither can still check that the renderer puts something on the canvas.
    let mut args = std::env::args().skip(2);
    if args.next().as_deref() == Some("--snapshot") {
        let out = args.next().unwrap_or_else(|| "snapshot.ppm".into());
        // Which frame, because the first one is often the least interesting: a run that fills up
        // over its length — a shot, a spreading spot — has nothing in it at `t = 0`, and a
        // snapshot of that is a picture of an empty box that still counts as "the renderer works".
        let at: usize = match args.next().as_deref() {
            Some("--frame") => args
                .next()
                .and_then(|n| n.parse().ok())
                .unwrap_or(0)
                .min(run.frames.len().saturating_sub(1)),
            _ => 0,
        };
        let mut app = App::new(run, panel);
        app.frame = at;
        println!("  snapshot of frame {at} of {}", app.run.frames.len());
        match app.snapshot(1100, 720) {
            Ok(pixels) => {
                write_ppm(&out, 1100, 720, &pixels);
                // **Compared against the corner pixel, not against a constant.** The target is
                // sRGB, so the clear colour is stored far brighter than the linear number the
                // pass was given — 56,66,77 rather than 10,14,19. A fixed threshold called every
                // pixel in the image a line and reported 100%, which is exactly the kind of
                // "measurement" a renderer check exists to avoid.
                let background = [pixels[0], pixels[1], pixels[2]];
                let lit = pixels
                    .chunks_exact(4)
                    .filter(|p| p[..3] != background)
                    .count();
                let mut shades: Vec<[u8; 3]> = Vec::new();
                for p in pixels.chunks_exact(4) {
                    let c = [p[0], p[1], p[2]];
                    if c != background && !shades.contains(&c) {
                        shades.push(c);
                    }
                }
                println!(
                    "  wrote {out} — {lit} of {} pixels carry a line ({:.2}%), in {} shades",
                    1100 * 720,
                    100.0 * lit as f64 / (1100.0 * 720.0),
                    shades.len()
                );
                if lit == 0 {
                    eprintln!("  nothing was drawn");
                    std::process::exit(1);
                }
            }
            Err(e) => {
                eprintln!("  no GPU available for a snapshot: {e}");
                std::process::exit(3);
            }
        }
        return;
    }

    let event_loop = EventLoop::new().expect("an event loop");
    event_loop.set_control_flow(ControlFlow::Poll);
    let mut app = App::new(run, panel);
    event_loop.run_app(&mut app).expect("the window runs");
}

/// Everything the window holds between frames.
struct App {
    run: Run,
    panel: String,
    framing: Framing,
    camera: Camera,
    /// The run-wide range the shading is measured against.
    span: (f64, f64),
    frame: usize,
    playing: bool,
    dragging: Option<(f64, f64)>,
    gpu: Option<Gpu>,
}

struct Gpu {
    window: Arc<Window>,
    surface: wgpu::Surface<'static>,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    pipeline: wgpu::RenderPipeline,
}

impl App {
    fn new(run: Run, panel: String) -> App {
        let framing = Framing::of(
            run.framing_of(&panel)
                .unwrap_or([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]),
        );
        // Once, from the whole run. Re-fitting it per frame is what `Run::scale_of` exists to
        // stop, and for a while nothing called it.
        let span = run.scale_of(&panel).unwrap_or((0.0, 1.0));
        // Framed to what is actually there. The window can still be zoomed; `--snapshot` cannot,
        // and a fixed distance is a distance chosen for a cube.
        let mut camera = Camera::default();
        camera.fit(
            run.framing_of(&panel)
                .unwrap_or([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]),
            &framing,
            16.0 / 9.0,
            0.85,
        );
        App {
            run,
            panel,
            framing,
            span,
            camera,
            frame: 0,
            playing: true,
            dragging: None,
            gpu: None,
        }
    }

    /// The line vertices for the current frame.
    fn vertices(&self, aspect: f64) -> Vec<Vertex> {
        let Some(panel) = self
            .run
            .frames
            .get(self.frame)
            .and_then(|f| f.panels.iter().find(|p| p.name() == self.panel))
        else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for s in segments(panel, &self.camera, &self.framing, aspect, self.span) {
            let colour = ramp(s.shade);
            out.push(Vertex {
                position: [s.from.x as f32, s.from.y as f32],
                colour,
            });
            out.push(Vertex {
                position: [s.to.x as f32, s.to.y as f32],
                colour,
            });
        }
        out
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.gpu.is_some() {
            return;
        }
        let attrs = Window::default_attributes()
            .with_title(format!("dualis — {}", self.run.title))
            .with_inner_size(winit::dpi::LogicalSize::new(1100.0, 720.0));
        let window = Arc::new(event_loop.create_window(attrs).expect("a window"));

        let instance = wgpu::Instance::default();
        let surface = instance
            .create_surface(window.clone())
            .expect("a surface on that window");
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: Some(&surface),
            force_fallback_adapter: false,
        }))
        .expect("a GPU adapter");
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("dualis viewer"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::downlevel_defaults(),
                memory_hints: wgpu::MemoryHints::default(),
            },
            None,
        ))
        .expect("a device");

        let size = window.inner_size();
        let caps = surface.get_capabilities(&adapter);
        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: caps.formats[0],
            width: size.width.max(1),
            height: size.height.max(1),
            present_mode: wgpu::PresentMode::Fifo,
            alpha_mode: caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        let pipeline = line_pipeline(&device, config.format);

        self.gpu = Some(Gpu {
            window,
            surface,
            device,
            queue,
            config,
            pipeline,
        });
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, _id: WindowId, event: WindowEvent) {
        let Some(gpu) = self.gpu.as_mut() else {
            return;
        };
        match event {
            WindowEvent::CloseRequested => event_loop.exit(),
            WindowEvent::Resized(size) => {
                gpu.config.width = size.width.max(1);
                gpu.config.height = size.height.max(1);
                gpu.surface.configure(&gpu.device, &gpu.config);
            }
            WindowEvent::MouseInput { state, .. } => {
                self.dragging = matches!(state, ElementState::Pressed).then_some((0.0, 0.0));
            }
            WindowEvent::CursorMoved { position, .. } => {
                if let Some((lx, ly)) = self.dragging {
                    let (x, y) = (position.x, position.y);
                    if lx != 0.0 || ly != 0.0 {
                        self.camera.turn((x - lx) * 0.008, (y - ly) * 0.006);
                    }
                    self.dragging = Some((x, y));
                }
            }
            WindowEvent::MouseWheel { delta, .. } => {
                let step = match delta {
                    MouseScrollDelta::LineDelta(_, y) => y as f64,
                    MouseScrollDelta::PixelDelta(p) => p.y / 60.0,
                };
                self.camera.zoom(1.0 - step * 0.1);
            }
            WindowEvent::KeyboardInput { event, .. } if event.state == ElementState::Pressed => {
                let n = self.run.frames.len().max(1);
                match event.physical_key {
                    PhysicalKey::Code(KeyCode::Space) => self.playing = !self.playing,
                    PhysicalKey::Code(KeyCode::ArrowRight) => {
                        self.playing = false;
                        self.frame = (self.frame + 1) % n;
                    }
                    PhysicalKey::Code(KeyCode::ArrowLeft) => {
                        self.playing = false;
                        self.frame = (self.frame + n - 1) % n;
                    }
                    PhysicalKey::Code(KeyCode::Escape) => event_loop.exit(),
                    _ => {}
                }
            }
            WindowEvent::RedrawRequested => self.draw(),
            _ => {}
        }
    }

    fn about_to_wait(&mut self, _event_loop: &ActiveEventLoop) {
        if let Some(gpu) = &self.gpu {
            gpu.window.request_redraw();
        }
    }
}

impl App {
    fn draw(&mut self) {
        let Some(gpu) = self.gpu.as_ref() else {
            return;
        };
        let aspect = gpu.config.width as f64 / gpu.config.height.max(1) as f64;
        let verts = self.vertices(aspect);

        let Ok(surface_texture) = gpu.surface.get_current_texture() else {
            return;
        };
        let view = surface_texture
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());
        let buffer = if verts.is_empty() {
            None
        } else {
            Some(
                gpu.device
                    .create_buffer_init_lite(bytemuck_lite::cast_slice(&verts)),
            )
        };

        let mut encoder = gpu
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("lines"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(BACKGROUND),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            if let Some(b) = &buffer {
                pass.set_pipeline(&gpu.pipeline);
                pass.set_vertex_buffer(0, b.slice(..));
                pass.draw(0..verts.len() as u32, 0..1);
            }
        }
        gpu.queue.submit(Some(encoder.finish()));
        surface_texture.present();

        if self.playing && !self.run.frames.is_empty() {
            self.frame = (self.frame + 1) % self.run.frames.len();
        }
    }
}

impl App {
    /// Render one frame with no window at all, and hand back the pixels as RGBA.
    ///
    /// The same camera, the same `viewer-core` segments and the same pipeline the window uses —
    /// only the target differs. A snapshot path that built its own vertices would be checking
    /// itself.
    fn snapshot(&self, width: u32, height: u32) -> Result<Vec<u8>, String> {
        let instance = wgpu::Instance::default();
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .ok_or("no adapter")?;
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("snapshot"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::downlevel_defaults(),
                memory_hints: wgpu::MemoryHints::default(),
            },
            None,
        ))
        .map_err(|e| e.to_string())?;

        let format = wgpu::TextureFormat::Rgba8UnormSrgb;
        let target = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("snapshot"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let view = target.create_view(&wgpu::TextureViewDescriptor::default());
        let pipeline = line_pipeline(&device, format);

        let verts = self.vertices(width as f64 / height as f64);
        let buffer = (!verts.is_empty())
            .then(|| device.create_buffer_init_lite(bytemuck_lite::cast_slice(&verts)));

        // Copies out of a texture want rows padded to 256 bytes, so the readback is wider than the
        // image and the padding is dropped below rather than left in the file as a smear.
        let unpadded = width * 4;
        let padded = unpadded.div_ceil(256) * 256;
        let readback = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: (padded * height) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        let mut encoder =
            device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("snapshot"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(BACKGROUND),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            if let Some(b) = &buffer {
                pass.set_pipeline(&pipeline);
                pass.set_vertex_buffer(0, b.slice(..));
                pass.draw(0..verts.len() as u32, 0..1);
            }
        }
        encoder.copy_texture_to_buffer(
            wgpu::ImageCopyTexture {
                texture: &target,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::ImageCopyBuffer {
                buffer: &readback,
                layout: wgpu::ImageDataLayout {
                    offset: 0,
                    bytes_per_row: Some(padded),
                    rows_per_image: Some(height),
                },
            },
            wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
        );
        queue.submit(Some(encoder.finish()));

        let slice = readback.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        device.poll(wgpu::Maintain::Wait);
        rx.recv()
            .map_err(|e| e.to_string())?
            .map_err(|e| e.to_string())?;

        let mapped = slice.get_mapped_range();
        let mut out = Vec::with_capacity((unpadded * height) as usize);
        for row in 0..height {
            let from = (row * padded) as usize;
            out.extend_from_slice(&mapped[from..from + unpadded as usize]);
        }
        drop(mapped);
        readback.unmap();
        Ok(out)
    }
}

/// A plain PPM, because a picture nothing can open is not evidence and an encoder is a dependency.
fn write_ppm(path: &str, width: u32, height: u32, rgba: &[u8]) {
    let mut out = format!("P6\n{width} {height}\n255\n").into_bytes();
    for p in rgba.chunks_exact(4) {
        out.extend_from_slice(&p[..3]);
    }
    if let Err(e) = std::fs::write(path, out) {
        eprintln!("cannot write {path}: {e}");
        std::process::exit(1);
    }
}

/// The clear colour, shared by the window and the snapshot so the two agree.
const BACKGROUND: wgpu::Color = wgpu::Color {
    r: 0.039,
    g: 0.055,
    b: 0.075,
    a: 1.0,
};

/// The line pipeline, built once here so the window and the snapshot cannot drift apart.
fn line_pipeline(device: &wgpu::Device, format: wgpu::TextureFormat) -> wgpu::RenderPipeline {
    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("lines"),
        source: wgpu::ShaderSource::Wgsl(SHADER.into()),
    });
    let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: None,
        bind_group_layouts: &[],
        push_constant_ranges: &[],
    });
    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("lines"),
        layout: Some(&layout),
        vertex: wgpu::VertexState {
            module: &shader,
            entry_point: "vs",
            buffers: &[wgpu::VertexBufferLayout {
                array_stride: std::mem::size_of::<Vertex>() as u64,
                step_mode: wgpu::VertexStepMode::Vertex,
                attributes: &[
                    wgpu::VertexAttribute {
                        offset: 0,
                        shader_location: 0,
                        format: wgpu::VertexFormat::Float32x2,
                    },
                    wgpu::VertexAttribute {
                        offset: 8,
                        shader_location: 1,
                        format: wgpu::VertexFormat::Float32x3,
                    },
                ],
            }],
            compilation_options: Default::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: &shader,
            entry_point: "fs",
            targets: &[Some(format.into())],
            compilation_options: Default::default(),
        }),
        primitive: wgpu::PrimitiveState {
            topology: wgpu::PrimitiveTopology::LineList,
            ..Default::default()
        },
        depth_stencil: None,
        multisample: wgpu::MultisampleState::default(),
        multiview: None,
        cache: None,
    })
}

/// A vertex buffer from bytes, without pulling in `wgpu::util`.
trait BufferInit {
    fn create_buffer_init_lite(&self, bytes: &[u8]) -> wgpu::Buffer;
}

impl BufferInit for wgpu::Device {
    fn create_buffer_init_lite(&self, bytes: &[u8]) -> wgpu::Buffer {
        let buffer = self.create_buffer(&wgpu::BufferDescriptor {
            label: Some("lines"),
            size: bytes.len() as u64,
            usage: wgpu::BufferUsages::VERTEX,
            mapped_at_creation: true,
        });
        buffer
            .slice(..)
            .get_mapped_range_mut()
            .copy_from_slice(bytes);
        buffer.unmap();
        buffer
    }
}

/// The same ramp the HTML report uses, so the two views of one run agree about colour.
fn ramp(t: f64) -> [f32; 3] {
    let t = t.clamp(0.0, 1.0);
    let stops = [
        [24.0, 52.0, 110.0],
        [61.0, 139.0, 255.0],
        [255.0, 194.0, 71.0],
        [255.0, 82.0, 33.0],
    ];
    let (a, b, u) = if t < 0.34 {
        (stops[0], stops[1], t / 0.34)
    } else if t < 0.67 {
        (stops[1], stops[2], (t - 0.34) / 0.33)
    } else {
        (stops[2], stops[3], (t - 0.67) / 0.33)
    };
    [
        ((a[0] + (b[0] - a[0]) * u) / 255.0) as f32,
        ((a[1] + (b[1] - a[1]) * u) / 255.0) as f32,
        ((a[2] + (b[2] - a[2]) * u) / 255.0) as f32,
    ]
}

const SHADER: &str = r#"
struct In { @location(0) pos: vec2<f32>, @location(1) colour: vec3<f32> };
struct Out { @builtin(position) clip: vec4<f32>, @location(0) colour: vec3<f32> };

@vertex
fn vs(v: In) -> Out {
    var out: Out;
    out.clip = vec4<f32>(v.pos, 0.0, 1.0);
    out.colour = v.colour;
    return out;
}

@fragment
fn fs(v: Out) -> @location(0) vec4<f32> {
    return vec4<f32>(v.colour, 1.0);
}
"#;
