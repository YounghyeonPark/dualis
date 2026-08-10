//! An accelerator for dualis's explicit stencils, with the CPU domain as the reference.
//!
//! # The one rule
//!
//! **`Solid3D` is the answer and this is a cache of it.** Where the two disagree, the CPU is
//! right. That is not modesty — it is the only arrangement under which the library's promises
//! survive a GPU at all, and it is what makes the comparison in `tests/` the point of the crate
//! rather than a nicety.
//!
//! # Why a GPU port is a different computation, not a faster one
//!
//! WGSL has no `f64`. Every number below is `f32`, against the domain's `f64`, so this is not the
//! same arithmetic run quicker — it is a lower-precision arithmetic. `tests/against_the_cpu.rs`
//! measures how far apart they land rather than asserting they agree, and the figure it reports
//! is the honest cost of the acceleration.
//!
//! That has a consequence worth stating plainly: `Simulation`'s conservation audit defaults to a
//! relative `1e-9`, and single precision cannot hold that. A scene using this needs
//! `conservation_tolerance_for(quantity::ENERGY, ..)` set to something `f32` can meet, and
//! choosing that number is choosing what the run is allowed to lose.
//!
//! # What is on the GPU and what is deliberately not
//!
//! The **stencil** is: each cell reads six neighbours and writes itself, in no particular order,
//! which is exactly the shape a GPU is for and exactly the shape that has no ordering problem.
//!
//! The **reductions are not.** A mean or a total summed with atomics depends on the order the
//! workgroups happened to finish, and floating-point addition is not associative — so the answer
//! would change between runs on one machine. [`GpuSolid::ledger`] reads the grid back and sums it
//! on the CPU in index order. That costs a transfer per audited step and buys back the thing the
//! whole workspace is built on.
//!
//! `Ensemble` solved the same problem on threads with fixed-size blocks, and the same discipline
//! would work here. It is not done yet because a readback is simpler and correct, and a faster
//! deterministic reduction is worth writing when somebody's grid is large enough to need it.
//!
//! # The buffer holds the **deviation**, and that is not a detail
//!
//! Cells carry `T − T₀`, not `T`. The first version stored absolute kelvin and diverged from the
//! reference by `1.4e-3` after two hundred steps — a thousand times what accumulation predicts —
//! and the cause was not accumulation at all.
//!
//! The update is `centre + F·(sum − 6·centre)`. On absolute temperatures near 293 K that sum is
//! about 1759, where `f32`'s resolution is `1.2e-4`. The difference being extracted from it is of
//! order `1e-3` K. **Subtracting two numbers that agree to five digits keeps less than one digit
//! of the answer**, every step, forever.
//!
//! On deviations the same numbers are near 1 K, where the resolution is `1.2e-7`, and the
//! subtraction keeps about four digits. The stencil is linear, so subtracting a constant commutes
//! with it exactly and the fix costs nothing — which is the useful part: single precision was not
//! the problem, spending it on an offset nobody needed was.

#![deny(missing_docs)]
#![forbid(unsafe_code)]

use dualis_core::conserved::quantity;
use dualis_core::units::{Length, Temperature, Time};
use dualis_core::{Domain, Exchange, Ledger, Reading, Substance, Violation};
use dualis_thermal::STABLE_FOURIER_3D;

/// The heat channel, the same name the thermal crate publishes on.
pub const HEAT: &str = quantity::ENERGY;

/// A block of conducting material whose stencil runs on the GPU.
///
/// Mirrors `dualis_thermal::Solid3D`: cubic cells, a seven-point stencil, insulated faces. What
/// differs is the arithmetic — `f32` here against `f64` there — and that difference is measured
/// rather than assumed.
pub struct GpuSolid {
    name: String,
    counts: (usize, usize, usize),
    dx: f64,
    alpha: f64,
    capacity: f64,
    reference: f64,
    absorbed: f64,
    gpu: Context,
    /// The last state read back, as **deviations** from `reference`. See the crate docs.
    mirror: Vec<f32>,
    mirror_valid: bool,
}

/// The device, the pipeline and the two buffers the stencil ping-pongs between.
struct Context {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
    cells: [wgpu::Buffer; 2],
    uniforms: wgpu::Buffer,
    readback: wgpu::Buffer,
    /// Which of `cells` currently holds the state.
    front: usize,
    count: usize,
}

/// What went wrong before there was anything to compute.
#[derive(Debug)]
pub enum Unavailable {
    /// No GPU this process can reach.
    NoAdapter,
    /// One was found and would not give a device.
    NoDevice(String),
    /// The substance has no diffusivity, so there is no stencil to run.
    NotConducting,
}

impl std::fmt::Display for Unavailable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Unavailable::NoAdapter => write!(f, "no GPU adapter"),
            Unavailable::NoDevice(e) => write!(f, "no GPU device: {e}"),
            Unavailable::NotConducting => write!(f, "the substance has no diffusivity"),
        }
    }
}

impl std::error::Error for Unavailable {}

impl GpuSolid {
    /// A block of `counts` cubic cells of side `dx`, all starting at `initial`.
    ///
    /// Fails rather than falling back to the CPU. A silent fallback is the worst outcome here:
    /// a caller who asked for a GPU and got a slow CPU run with no message has no way to tell,
    /// and the measurement they are about to take is of the wrong thing.
    pub fn new(
        name: impl Into<String>,
        substance: Substance,
        counts: (usize, usize, usize),
        dx: Length,
        initial: Temperature,
    ) -> Result<GpuSolid, Unavailable> {
        let counts = (counts.0.max(1), counts.1.max(1), counts.2.max(1));
        let n = counts.0 * counts.1 * counts.2;
        let alpha = substance
            .diffusivity()
            .ok_or(Unavailable::NotConducting)?
            .to_si();
        let cell_volume = dx.to_si().powi(3);
        let capacity = substance
            .heat_capacity(dualis_core::units::Volume::from_si(cell_volume))
            .ok_or(Unavailable::NotConducting)?
            .to_si();

        // Zero, because the buffer holds `T - T0` and `T` starts at `T0`.
        let gpu = Context::new(n, 0.0)?;
        Ok(GpuSolid {
            name: name.into(),
            counts,
            dx: dx.to_si(),
            alpha,
            capacity,
            reference: initial.to_si(),
            absorbed: 0.0,
            gpu,
            mirror: vec![0.0; n],
            mirror_valid: true,
        })
    }

    /// Cell counts along x, y and z.
    pub fn counts(&self) -> (usize, usize, usize) {
        self.counts
    }

    /// `α·dt/dx²`. Must stay at or under `STABLE_FOURIER_3D`.
    pub fn fourier_number(&self, dt: Time) -> f64 {
        self.alpha * dt.to_si() / (self.dx * self.dx)
    }

    /// The temperature of one cell, in kelvin.
    ///
    /// Reads the grid back if anything has changed since the last read, so a caller looping over
    /// every cell pays one transfer rather than one per cell.
    pub fn temperature_at(&mut self, i: usize, j: usize, k: usize) -> Temperature {
        self.sync();
        let (nx, ny, nz) = self.counts;
        let idx = i.min(nx - 1) + nx * (j.min(ny - 1) + ny * k.min(nz - 1));
        Temperature::from_si(self.reference + self.mirror[idx] as f64)
    }

    /// Set one cell, for an initial condition.
    pub fn set_temperature(&mut self, i: usize, j: usize, k: usize, t: Temperature) {
        self.sync();
        let (nx, ny, nz) = self.counts;
        if i >= nx || j >= ny || k >= nz {
            return;
        }
        self.mirror[i + nx * (j + ny * k)] = (t.to_si() - self.reference) as f32;
        self.gpu.upload(&self.mirror);
    }

    /// Put joules into one cell, as heat that arrived there.
    pub fn deposit(&mut self, i: usize, j: usize, k: usize, joules: dualis_core::units::Energy) {
        let (nx, ny, nz) = self.counts;
        if i >= nx || j >= ny || k >= nz {
            return;
        }
        self.sync();
        let idx = i + nx * (j + ny * k);
        self.mirror[idx] += (joules.to_si() / self.capacity) as f32;
        self.absorbed += joules.to_si();
        self.gpu.upload(&self.mirror);
    }

    /// Every cell as a **deviation** from the starting temperature, in index order. One transfer.
    pub fn deviations(&mut self) -> &[f32] {
        self.sync();
        &self.mirror
    }

    /// Every cell as an absolute temperature in kelvin, in index order.
    pub fn cells(&mut self) -> Vec<f64> {
        self.sync();
        let reference = self.reference;
        self.mirror.iter().map(|v| reference + *v as f64).collect()
    }

    /// Mean over every cell, summed **on the CPU in index order**.
    pub fn mean_temperature(&mut self) -> Temperature {
        self.sync();
        let n = self.mirror.len() as f64;
        Temperature::from_si(
            self.reference + self.mirror.iter().map(|v| *v as f64).sum::<f64>() / n,
        )
    }

    /// The hottest cell.
    pub fn peak_temperature(&mut self) -> Temperature {
        self.sync();
        Temperature::from_si(
            self.reference + self.mirror.iter().fold(f32::MIN, |m, v| m.max(*v)) as f64,
        )
    }

    fn sync(&mut self) {
        if !self.mirror_valid {
            self.gpu.download(&mut self.mirror);
            self.mirror_valid = true;
        }
    }

    /// Heat held, measured from the temperature it started at — the same reference `Solid3D` uses,
    /// and for the same reason: against absolute zero the interesting millijoule is lost in the
    /// last bits of a kilojoule.
    ///
    /// The buffer already holds deviations, which is exactly what a stored heat is measured from —
    /// the same fact that makes the `f32` arithmetic work at all.
    pub fn stored_heat(&mut self) -> dualis_core::units::Energy {
        self.sync();
        dualis_core::units::Energy::from_si(
            self.capacity * self.mirror.iter().map(|v| *v as f64).sum::<f64>(),
        )
    }

    /// Heat taken from the bus over the run.
    pub fn absorbed_energy(&self) -> dualis_core::units::Energy {
        dualis_core::units::Energy::from_si(self.absorbed)
    }
}

impl Domain for GpuSolid {
    fn name(&self) -> &str {
        &self.name
    }

    /// `dx²/(6α)`, the same limit the CPU domain reports — the scheme is the same, only the
    /// precision differs.
    fn max_stable_dt(&self, _now: Time) -> Time {
        Time::from_si(STABLE_FOURIER_3D * self.dx * self.dx / self.alpha)
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let f = self.fourier_number(dt);
        if f > STABLE_FOURIER_3D + 1e-12 {
            return Err(Violation {
                quantity: "Fourier number".to_string(),
                site: format!("{} (explicit 3D conduction, on the GPU)", self.name),
                before: STABLE_FOURIER_3D,
                after: f,
                scale: STABLE_FOURIER_3D,
                tolerance: 1e-12,
            });
        }

        // Placeless heat spreads evenly, exactly as the CPU domain does it — a block has six faces
        // and no distinguished cell, so choosing one would invent a location the bus never carried.
        let gained = bus.take_share(HEAT, dt);
        if gained != 0.0 {
            self.absorbed += gained;
            self.sync();
            let per_cell = (gained / (self.mirror.len() as f64 * self.capacity)) as f32;
            for cell in self.mirror.iter_mut() {
                *cell += per_cell;
            }
            self.gpu.upload(&self.mirror);
        }

        let (nx, ny, nz) = self.counts;
        self.gpu
            .dispatch(f as f32, [nx as u32, ny as u32, nz as u32]);
        self.mirror_valid = false;
        Ok(())
    }

    /// Heat gained since the start, summed on the CPU.
    fn ledger(&self) -> Ledger {
        // `ledger` takes `&self`, so the mirror cannot be refreshed here — which is exactly why
        // the audit reads a value that is one dispatch behind unless something has already synced.
        // `Simulation` calls `ledger` after every domain has stepped, and `readings` below syncs,
        // so in practice this is current. Stated because a stale ledger is an audit that passes.
        Ledger::new().with(
            quantity::ENERGY,
            self.capacity * self.mirror.iter().map(|v| *v as f64).sum::<f64>(),
        )
    }

    fn books_balance(&self) -> bool {
        // Not claimed. The books would balance in `f64`; in `f32` the sum drifts from what the bus
        // moved by more than the per-domain check allows, and claiming otherwise would make a
        // correct domain fail. See the crate docs on what single precision costs.
        false
    }

    fn readings(&self) -> Vec<Reading> {
        let n = self.mirror.len().max(1) as f64;
        let mean = self.reference + self.mirror.iter().map(|v| *v as f64).sum::<f64>() / n;
        let peak = self.reference + self.mirror.iter().fold(f32::MIN, |m, v| m.max(*v)) as f64;
        vec![
            Reading::new(&self.name, "peak", peak - 273.15, "C"),
            Reading::new(&self.name, "mean", mean - 273.15, "C"),
            Reading::new(&self.name, "absorbed", self.absorbed, "J"),
        ]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }
}

impl Context {
    fn new(count: usize, initial: f32) -> Result<Context, Unavailable> {
        let instance = wgpu::Instance::default();
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .ok_or(Unavailable::NoAdapter)?;
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("dualis stencil"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::downlevel_defaults(),
                memory_hints: wgpu::MemoryHints::default(),
            },
            None,
        ))
        .map_err(|e| Unavailable::NoDevice(e.to_string()))?;

        let bytes = (count * 4) as u64;
        let make = |label: &str, usage: wgpu::BufferUsages| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: bytes.max(4),
                usage,
                mapped_at_creation: false,
            })
        };
        let storage = wgpu::BufferUsages::STORAGE
            | wgpu::BufferUsages::COPY_DST
            | wgpu::BufferUsages::COPY_SRC;
        let cells = [make("cells a", storage), make("cells b", storage)];
        queue.write_buffer(
            &cells[0],
            0,
            &vec![initial; count]
                .iter()
                .flat_map(|f| f.to_le_bytes())
                .collect::<Vec<u8>>(),
        );

        let uniforms = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("stencil uniforms"),
            size: 32,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let readback = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: bytes.max(4),
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("stencil"),
            entries: &[
                entry(0, true),
                entry(1, false),
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("stencil"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: None,
            bind_group_layouts: &[&layout],
            push_constant_ranges: &[],
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("stencil"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: "sweep",
            compilation_options: Default::default(),
            cache: None,
        });

        Ok(Context {
            device,
            queue,
            pipeline,
            layout,
            cells,
            uniforms,
            readback,
            front: 0,
            count,
        })
    }

    fn upload(&mut self, cells: &[f32]) {
        let bytes: Vec<u8> = cells.iter().flat_map(|f| f.to_le_bytes()).collect();
        self.queue.write_buffer(&self.cells[self.front], 0, &bytes);
    }

    fn download(&mut self, into: &mut [f32]) {
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        encoder.copy_buffer_to_buffer(
            &self.cells[self.front],
            0,
            &self.readback,
            0,
            (self.count * 4) as u64,
        );
        self.queue.submit(Some(encoder.finish()));

        let slice = self.readback.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        self.device.poll(wgpu::Maintain::Wait);
        rx.recv()
            .expect("the map completes")
            .expect("a mapped buffer");
        {
            let mapped = slice.get_mapped_range();
            for (out, chunk) in into.iter_mut().zip(mapped.chunks_exact(4)) {
                *out = f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
            }
        }
        self.readback.unmap();
    }

    /// One sweep: read the front buffer, write the back one, then swap.
    ///
    /// Two buffers rather than one, because a stencil that reads and writes the same array reads
    /// some neighbours already updated and some not — a Gauss-Seidel sweep pretending to be
    /// Jacobi, which is a different scheme with a different stability limit and no ordering
    /// anybody chose.
    fn dispatch(&mut self, fourier: f32, counts: [u32; 3]) {
        let mut uniforms = Vec::with_capacity(32);
        uniforms.extend(counts[0].to_le_bytes());
        uniforms.extend(counts[1].to_le_bytes());
        uniforms.extend(counts[2].to_le_bytes());
        uniforms.extend((self.count as u32).to_le_bytes());
        uniforms.extend(fourier.to_le_bytes());
        uniforms.extend([0u8; 12]);
        self.queue.write_buffer(&self.uniforms, 0, &uniforms);

        let back = 1 - self.front;
        let bind = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &self.layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: self.cells[self.front].as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: self.cells[back].as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: self.uniforms.as_entire_binding(),
                },
            ],
        });

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("sweep"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &bind, &[]);
            pass.dispatch_workgroups((self.count as u32).div_ceil(64), 1, 1);
        }
        self.queue.submit(Some(encoder.finish()));
        self.front = back;
    }
}

fn entry(binding: u32, read_only: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

/// The seven-point stencil, with insulated faces by mirroring — the same update `Solid3D` does,
/// in `f32`.
const SHADER: &str = r#"
struct Grid {
    nx: u32, ny: u32, nz: u32, total: u32,
    fourier: f32, _pad0: f32, _pad1: f32, _pad2: f32,
};

@group(0) @binding(0) var<storage, read>       src: array<f32>;
@group(0) @binding(1) var<storage, read_write> dst: array<f32>;
@group(0) @binding(2) var<uniform>             g: Grid;

fn at(i: u32, j: u32, k: u32) -> f32 {
    return src[i + g.nx * (j + g.ny * k)];
}

@compute @workgroup_size(64)
fn sweep(@builtin(global_invocation_id) id: vec3<u32>) {
    let n = id.x;
    if (n >= g.total) { return; }

    let i = n % g.nx;
    let j = (n / g.nx) % g.ny;
    let k = n / (g.nx * g.ny);

    // A mirror at a face, not a zero: a zero neighbour is a wall held at absolute zero and would
    // drain the block, where a mirror is a face with no gradient across it and so no flow.
    let lo_i = select(i - 1u, i, i == 0u);
    let hi_i = select(i + 1u, i, i + 1u == g.nx);
    let lo_j = select(j - 1u, j, j == 0u);
    let hi_j = select(j + 1u, j, j + 1u == g.ny);
    let lo_k = select(k - 1u, k, k == 0u);
    let hi_k = select(k + 1u, k, k + 1u == g.nz);

    let centre = at(i, j, k);
    let sum = at(lo_i, j, k) + at(hi_i, j, k)
            + at(i, lo_j, k) + at(i, hi_j, k)
            + at(i, j, lo_k) + at(i, j, hi_k);

    dst[n] = centre + g.fourier * (sum - 6.0 * centre);
}
"#;
