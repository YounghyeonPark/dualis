//! A portafilter, and the water going through it.
//!
//! ```text
//! cargo run --release --example portafilter_flow                # the numbers
//! cargo run --release --example portafilter_flow flow.html      # and the machine, turning
//! cargo run --release --example portafilter_flow flow.gltf      # or into Blender
//! cargo run --release --example portafilter_flow flow.json      # or the native window:
//! #   cd runtime/viewer && cargo run --release -- ../../flow.json
//! ```
//!
//! `espresso_shot` cuts the puck open and colours the inside. This one does the other thing: it
//! builds the hardware — shower screen, basket, body, spout — and lets you watch parcels of water
//! leave the screen, work down through the grounds, and come out the bottom darker than they went
//! in.
//!
//! # Where the water actually goes, and who decided
//!
//! Nothing here places a streamline by hand. `Puck` solves `∇·((k/μ)∇p) = 0` on the permeability
//! the grind and the packing give, and a parcel is advected by the **pore** velocity of that
//! solved field — `u/ε`, not `u`, which is a factor of 2.2 and the commonest error in porous
//! transport. Where a parcel goes is a consequence; the picture is a readout, not a drawing.
//!
//! The colour is what it is carrying. A parcel picks up dissolved coffee from every cell it
//! crosses, so it leaves the screen at zero and arrives at the spout loaded — which is why an
//! espresso *looks* like that, and it is the same field the yield is computed from.
//!
//! # What is solved and what is only drawn
//!
//! Stated plainly, because a picture does not distinguish them and a reader will assume the
//! flattering answer:
//!
//! | | |
//! | --- | --- |
//! | inside the basket | **solved.** Darcy, advection, dissolution |
//! | the headspace above the puck | drawn. A full plenum, so the mean speed is the superficial velocity — right, and not interesting |
//! | below the basket floor | drawn. A continuous jet — see below |
//! | the hardware | drawn. Geometry, no physics |
//!
//! The two baskets differ only in the ring against the wall — 0.60 against 0.45 — and everything
//! you can see between them follows from that one number.
//!
//! # The stream out of the spout is a stream, not parcels
//!
//! A parcel falls the 75 mm from the basket to the cup in about 0.12 s, and a frame is 0.7 s
//! apart, so sampling parcels there catches roughly none of them: the first version drew the
//! basket beautifully and had nothing at all between the spout and the cup.
//!
//! What is physically there is a **continuous jet**, so that is what is drawn — strands from the
//! spout to the cup, coloured by the concentration leaving the basket at that instant. It carries
//! the outlet's real number and none of its own; a parcel is a sample of something discrete, and
//! a stream is not discrete.
//!
//! # Read the two side by side as *equal time*, not equal weight
//!
//! Both run the same 25 s, so the channelled one is 42% heavier in the cup and its yield reads
//! *higher*. That is not the gap being an improvement — it is 42% more water having gone through.
//! `espresso_shot` makes the comparison a barista actually makes, which is to a weight, and there
//! the channelled shot loses. Here the comparison is to a clock, because two baskets on two
//! different clocks cannot share a frame.

mod common;

use common::{check, check_between, heading};
use pantometry::prelude::*;
use pantometry_view::{gltf, report};

/// A 58 mm basket in a 66 mm box on a 2 mm grid: a 20 mm bed with a 4 mm jacket.
const NX: usize = 33;
const NY: usize = 10;
const NZ: usize = 33;
const DX: f64 = 2e-3;
const RADIUS: f64 = 29e-3;
const POROSITY: f64 = 0.45;
/// The loose ring, as a puck that shrank from the wall.
const GAP_POROSITY: f64 = 0.60;

/// Headspace between the screen and the bed.
const HEADSPACE: f64 = 3e-3;
/// How many parcels are in flight at once.
const PARCELS: usize = 160;
/// How many positions of its own history a parcel draws behind it.
///
/// Eight frames at the pore velocity is about 7 mm, which reads as a streak on a 20 mm bed. Five
/// was 4 mm and read as a tick.
const TRAIL: usize = 8;
/// Frames in the animation, and how far apart.
const FRAMES: usize = 36;
const FRAME_S: f64 = 0.7;

/// Where a parcel is in its journey.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Phase {
    /// Between the screen and the top of the bed.
    Headspace,
    /// Inside the grounds, advected by the solved field.
    Bed,
    /// Out of the basket holes and falling to the spout.
    Jet,
    /// In the cup, and no longer drawn.
    Cup,
}

/// One parcel of water.
#[derive(Clone, Debug)]
struct Parcel {
    /// Position in the puck's own coordinates: `y` is the flow axis, increasing downward.
    at: [f64; 3],
    /// The last few positions, newest last.
    trail: Vec<[f64; 3]>,
    /// Dissolved coffee it is carrying, kg/m³.
    load: f64,
    phase: Phase,
    /// When it reached the top of the bed. `None` while still in the headspace.
    entered: Option<f64>,
    /// When it reached the basket floor. `None` until it does.
    through: Option<f64>,
    /// The radius it entered the bed at, for the ring statistics.
    entry_r: f64,
    /// Downward speed while falling free.
    fall: f64,
}

fn machine(channel: bool) -> Puck {
    let mut p = Puck::new(
        if channel { "wall gap" } else { "even" },
        Basket {
            counts: (NX, NY, NZ),
            cell: Length::from_si(DX),
            radius: Length::from_si(RADIUS),
            porosity: POROSITY,
            ..Basket::espresso()
        },
    );
    if channel {
        let ring = wall_ring(&p);
        p.repack(GAP_POROSITY, |i, j, k| ring[i + NX * (j + NY * k)]);
    }
    p
}

fn wall_ring(p: &Puck) -> Vec<bool> {
    let mut ring = vec![false; NX * NY * NZ];
    for k in 0..NZ {
        for j in 0..NY {
            for i in 0..NX {
                if !p.is_packed(i, j, k) {
                    continue;
                }
                ring[i + NX * (j + NY * k)] = i == 0
                    || k == 0
                    || i + 1 == NX
                    || k + 1 == NZ
                    || !p.is_packed(i - 1, j, k)
                    || !p.is_packed(i + 1, j, k)
                    || !p.is_packed(i, j, k - 1)
                    || !p.is_packed(i, j, k + 1);
            }
        }
    }
    ring
}

/// The centre of the basket in the puck's coordinates.
fn centre() -> (f64, f64) {
    (0.5 * NX as f64 * DX, 0.5 * NZ as f64 * DX)
}

/// Puck coordinates to drawing coordinates: centred on the axis, and `y` pointing **up**.
///
/// The puck's `y` increases along the flow, which is downward. A viewer draws `+y` up, so a
/// straight pass-through would hang the machine from the cup.
fn draw(p: [f64; 3]) -> [f64; 3] {
    let (cx, cz) = centre();
    [p[0] - cx, -p[1], p[2] - cz]
}

fn main() {
    heading("1. The two baskets, and the flow each one passes");

    let mut even = machine(false);
    let mut gap = machine(true);
    println!(
        "  {:<38} {:>10.2} g      identical dose; only the ring differs",
        "dose",
        even.dose().to_si() * 1000.0
    );

    // The flow split has a closed form, and it is an equality rather than a bound: every column
    // of cells from the screen to the holes is an independent series chain, so the basket is
    // columns in parallel and its conductance is their sum. Widening the ring changes only the
    // ring's term.
    let f = ring_fraction(&even);
    let mobility = |e: f64| e.powi(3) / (1.0 - e).powi(2);
    let predicted = (1.0 - f) + f * mobility(GAP_POROSITY) / mobility(POROSITY);
    let measured = gap.flow_rate().to_si() / even.flow_rate().to_si();
    println!(
        "  {:<38} {:>10.1} %      of the cross-section is the ring, counted rather than estimated",
        "ring",
        f * 100.0
    );
    check(
        "columns in parallel give the split",
        measured,
        predicted,
        1e-9,
        "x",
    );

    let pore = |p: &Puck| {
        p.flow_rate().to_si() / Liquid::water().density.to_si() / open_area(p) / POROSITY
    };
    let transit = NY as f64 * DX / pore(&even);
    println!(
        "  {:<38} {:>10.1} s      the pore transit, eps L / u -- not L / u, which is {:.1} s",
        "how long the water is in the grounds",
        transit,
        transit / POROSITY
    );

    // ========================================================================= 2. the parcels
    heading("2. Two hundred parcels, advected by the field rather than placed");

    let even_run = trace(&mut even, 1);
    let gap_run = trace(&mut gap, 2);

    for run in [&even_run, &gap_run] {
        println!(
            "  {:<20} {:>5} released, {:>4} still in the grounds, {:>4} in the cup",
            run.name, run.released, run.in_bed, run.arrived
        );
        assert_eq!(
            run.released,
            run.in_head + run.in_bed + run.in_jet + run.arrived,
            "{}: every parcel is somewhere",
            run.name
        );
        assert!(
            run.max_radius <= RADIUS + 1e-9,
            "{}: nothing leaves through the basket wall: {:.4} mm against {:.4} mm",
            run.name,
            run.max_radius * 1000.0,
            RADIUS * 1000.0
        );
    }

    // **The tracer reproduces the transit time the flow field implies**, which is the check that
    // matters: the advection is a separate piece of arithmetic from the pressure solve, and a
    // parcel moving at the Darcy velocity instead of the pore velocity would take 2.2x as long
    // and every picture would still look plausible.
    check_between(
        "mean transit against eps L / u",
        even_run.mean_transit / transit,
        0.85,
        1.15,
        "x",
    );
    assert!(
        even_run.mean_transit < 0.6 * transit / POROSITY,
        "the parcels are on the pore velocity, not the Darcy one: {:.2} s against {:.2} s",
        even_run.mean_transit,
        transit / POROSITY
    );

    // And they arrive carrying what the field says they crossed.
    println!(
        "  {:<38} {:>10.1} kg/m3  what a parcel has picked up by the spout",
        "load at the exit", even_run.mean_load
    );
    check_between(
        "the exit load against the outlet TDS",
        even_run.mean_load / (even.tds() * Liquid::water().density.to_si()),
        0.6,
        1.6,
        "x",
    );

    // ============================================================================== 3. the gap
    heading("3. What the ring does to where the water goes");

    println!(
        "  {:<38} {:>10.1} s      against {:.1} s through the even bed",
        "mean transit through the ring", gap_run.ring_transit, even_run.ring_transit
    );
    println!(
        "  {:<38} {:>10.2} x      the ring is faster than the core by this much",
        "ring against core, in speed",
        gap_run.core_transit / gap_run.ring_transit
    );
    assert!(
        gap_run.core_transit / gap_run.ring_transit > 1.3,
        "the gap must make the ring visibly quicker: {:.2}x",
        gap_run.core_transit / gap_run.ring_transit
    );
    assert!(
        (even_run.core_transit / even_run.ring_transit - 1.0).abs() < 0.15,
        "and an even bed must not: {:.3}x",
        even_run.core_transit / even_run.ring_transit
    );
    println!(
        "  {:<38} {:>10.1} kg/m3  against {:.1} through the even bed -- the same water, less coffee",
        "what the ring's water carries out", gap_run.ring_load, even_run.ring_load
    );
    // Said out loud, because the animation puts the two yields side by side at equal *time* and
    // the channelled one reads higher. It is higher because 42% more water went through it, not
    // because the gap helped -- `espresso_shot` compares to a weight, and there it loses.
    println!(
        "  {:<38} {:>10} both run the same clock, so the gap's cup is heavier and its yield reads",
        "note", ""
    );
    println!(
        "  {:<38} {:>10} higher. Compare to a weight, as espresso_shot does, and it loses.",
        "", ""
    );
    assert!(
        gap_run.ring_load < even_run.ring_load,
        "water that hurried through picked up less: {:.1} against {:.1}",
        gap_run.ring_load,
        even_run.ring_load
    );

    // ================================================================================== output
    match common::output_path() {
        Some(path) if path.ends_with(".gltf") => {
            let frame = Frame {
                time_s: even_run.frames[FRAMES - 1].time_s,
                panels: vec![
                    even_run.frames[FRAMES - 1].panels[0].clone(),
                    gap_run.frames[FRAMES - 1].panels[0].clone(),
                ],
                readings: Vec::new(),
            };
            let out = gltf::gltf("A portafilter, and the water through it", &frame);
            for why in &out.skipped {
                println!("  skipped: {why}");
            }
            common::write(&path, &out.document);
        }
        Some(path) => {
            let mut frames = Vec::with_capacity(FRAMES);
            for n in 0..FRAMES {
                frames.push(Frame {
                    time_s: even_run.frames[n].time_s,
                    panels: vec![
                        even_run.frames[n].panels[0].clone(),
                        gap_run.frames[n].panels[0].clone(),
                    ],
                    readings: even_run.frames[n]
                        .readings
                        .iter()
                        .chain(gap_run.frames[n].readings.iter())
                        .cloned()
                        .collect(),
                });
            }
            let title = "A portafilter, and the water through it";
            if path.ends_with(".json") {
                // The wire format, and `runtime/viewer` reads it. That viewer does not link this
                // library at all — it takes the file and nothing else — so a run that opens in
                // the native window is the format demonstrating it carries enough to draw a run.
                println!("\n  {FRAMES} frames. Open it in the native window:");
                println!("    cd runtime/viewer && cargo run --release -- ../../{path}");
                common::write(&path, &pantometry_view::to_json(title, &frames));
            } else {
                println!(
                    "\n  {FRAMES} frames of two baskets, {:.1} s apart. Drag to rotate, scroll \
                     to zoom.",
                    FRAME_S
                );
                common::write(&path, &report::html(title, &frames));
            }
        }
        None => println!(
            "\n  Pass a filename for the machine:\n    cargo run --release --example \
             portafilter_flow flow.html\n    cargo run --release --example portafilter_flow \
             flow.gltf\n    cargo run --release --example portafilter_flow flow.json"
        ),
    }
}

/// The concentration of what is leaving the basket, weighted by where it leaves.
///
/// The domain's own outlet, not an average over the parcels near it: the parcels are a sample and
/// this is the quantity the yield is computed from.
fn outlet_concentration(puck: &Puck) -> f64 {
    let (mut num, mut den) = (0.0, 0.0);
    for k in 0..NZ {
        for i in 0..NX {
            if !puck.is_packed(i, NY - 1, k) {
                continue;
            }
            let q = puck.pore_velocity_at(i, NY - 1, k).y.max(0.0);
            num += q * puck.concentration_at(i, NY - 1, k).to_si();
            den += q;
        }
    }
    if den > 0.0 {
        num / den
    } else {
        0.0
    }
}

/// The open cross-section the basket has on this grid.
fn open_area(p: &Puck) -> f64 {
    let mut n = 0;
    for k in 0..NZ {
        for i in 0..NX {
            if p.is_packed(i, 0, k) {
                n += 1;
            }
        }
    }
    n as f64 * DX * DX
}

/// The fraction of the cross-section the wall ring covers.
fn ring_fraction(p: &Puck) -> f64 {
    let ring = wall_ring(p);
    let (mut packed, mut edge) = (0usize, 0usize);
    for k in 0..NZ {
        for i in 0..NX {
            if p.is_packed(i, 0, k) {
                packed += 1;
                if ring[i + NX * (NY * k)] {
                    edge += 1;
                }
            }
        }
    }
    edge as f64 / packed as f64
}

/// Everything one basket's run produced.
struct Traced {
    name: &'static str,
    frames: Vec<Frame>,
    released: usize,
    in_head: usize,
    in_bed: usize,
    in_jet: usize,
    arrived: usize,
    /// Mean time in the grounds, over the parcels that got through.
    mean_transit: f64,
    /// The same, split by where the parcel entered.
    ring_transit: f64,
    core_transit: f64,
    /// Mean dissolved load at the spout, and the ring's share of it.
    mean_load: f64,
    ring_load: f64,
    /// The furthest any parcel got from the axis while in the bed.
    max_radius: f64,
}

/// Run the shot, advecting parcels through the field as it goes.
fn trace(puck: &mut Puck, seed: u64) -> Traced {
    let name: &'static str = if puck.name().starts_with("wall") {
        "wall gap"
    } else {
        "even"
    };
    let hardware = hardware_paths();
    let depth = NY as f64 * DX;
    let floor = depth;
    // Superficial velocity, which is the mean speed of water in the full plenum above the bed.
    let superficial = puck.flow_rate().to_si() / Liquid::water().density.to_si() / open_area(puck);

    let mut parcels: Vec<Parcel> = Vec::new();
    let mut bus = Exchange::new();
    let mut t = 0.0;
    let mut released = 0usize;
    let mut done: Vec<(f64, f64, f64)> = Vec::new(); // (transit, load, entry radius)
    let mut max_radius: f64 = 0.0;
    let mut frames = Vec::with_capacity(FRAMES);

    // Parcels are released steadily so the stream is continuous rather than a pulse. The rate is
    // set to keep `PARCELS` of them in flight once the pipeline has filled.
    let flight = HEADSPACE / superficial.max(1e-12) + depth / superficial.max(1e-12) * POROSITY;
    let per_second = PARCELS as f64 / flight.max(1e-9);
    let mut owed = 0.0;

    for n in 0..FRAMES {
        // Advance the physics to this frame, in steps the domain will accept.
        let target = n as f64 * FRAME_S;
        while t < target {
            let dt = (puck.max_stable_dt(Time::from_si(t)).to_si() * 0.5).min(target - t);
            if dt <= 0.0 {
                break;
            }
            puck.step(Time::from_si(t), Time::from_si(dt), &mut bus)
                .expect("stable");
            t += dt;
        }

        // Release what this frame is owed, on the screen's disc.
        owed += per_second * FRAME_S;
        while owed >= 1.0 {
            owed -= 1.0;
            let mut rng = Rng::for_index(seed, released as u64);
            // Uniform on a disc: the radius has to go as sqrt, or the middle is crowded.
            let r = RADIUS * rng.unit().sqrt();
            let a = rng.range(0.0, std::f64::consts::TAU);
            let (cx, cz) = centre();
            let at = [cx + r * a.cos(), -HEADSPACE, cz + r * a.sin()];
            parcels.push(Parcel {
                at,
                trail: vec![at],
                load: 0.0,
                phase: Phase::Headspace,
                entered: None,
                through: None,
                entry_r: r,
                fall: 0.0,
            });
            released += 1;
        }

        // Advect. Substeps because a parcel must not cross a cell in one go: the field is
        // piecewise constant, and a jump larger than a cell skips whatever was in between.
        let sub = 8;
        for _ in 0..sub {
            let h = FRAME_S / sub as f64;
            for p in parcels.iter_mut() {
                advance(p, puck, t, h, superficial, floor);
            }
        }
        for p in parcels.iter_mut() {
            if p.phase == Phase::Bed {
                let (cx, cz) = centre();
                let r = ((p.at[0] - cx).powi(2) + (p.at[2] - cz).powi(2)).sqrt();
                max_radius = max_radius.max(r);
            }
            p.trail.push(p.at);
            if p.trail.len() > TRAIL {
                p.trail.remove(0);
            }
        }
        for p in parcels.iter() {
            if p.phase == Phase::Cup {
                if let (Some(entered), Some(through)) = (p.entered, p.through) {
                    done.push((through - entered, p.load, p.entry_r));
                }
            }
        }
        parcels.retain(|p| p.phase != Phase::Cup);

        // The picture: the hardware and the water, in one panel, so they occlude each other.
        let mut runs: Vec<Vec<[f64; 3]>> = hardware.0.clone();
        let mut values: Vec<f64> = hardware.1.clone();
        for p in parcels.iter() {
            if p.trail.len() < 2 || p.phase == Phase::Jet {
                continue;
            }
            runs.push(p.trail.iter().map(|q| draw(*q)).collect());
            values.push(p.load);
        }
        // The jet, once anything has come through. Drawn as a stream because it is one — see the
        // module docs. Its colour is the flow-weighted concentration leaving the basket, which is
        // the domain's number and not the parcels'.
        if done.len() + parcels.iter().filter(|q| q.phase == Phase::Jet).count() > 0 {
            let out = outlet_concentration(puck);
            for strand in 0..3 {
                let a = std::f64::consts::TAU * strand as f64 / 3.0;
                let (wx, wz) = (1.1e-3 * a.cos(), 1.1e-3 * a.sin());
                runs.push(
                    (0..=10)
                        .map(|n| {
                            // Narrowing as it falls, the way a jet does once gravity has it.
                            let u = n as f64 / 10.0;
                            let taper = 1.0 - 0.55 * u;
                            [wx * taper, -(depth + 0.044 + u * 0.051), wz * taper]
                        })
                        .collect(),
                );
                values.push(out);
            }
        }
        frames.push(Frame {
            time_s: t,
            panels: vec![Panel {
                name: format!("{name} — the machine and the water"),
                unit: "kg/m3",
                data: PanelData::paths(runs, values),
            }],
            readings: vec![
                Reading::new(name, "in the cup", puck.delivered().to_si() * 1000.0, "g"),
                Reading::new(name, "yield", puck.yield_fraction() * 100.0, "%"),
                Reading::new(name, "parcels through", done.len() as f64, ""),
            ],
        });
    }

    let mean = |f: &dyn Fn(&(f64, f64, f64)) -> bool, g: &dyn Fn(&(f64, f64, f64)) -> f64| {
        let picked: Vec<f64> = done.iter().filter(|d| f(d)).map(g).collect();
        if picked.is_empty() {
            0.0
        } else {
            picked.iter().sum::<f64>() / picked.len() as f64
        }
    };
    // The ring is the outer annulus of the same area fraction the loose cells cover, so "ring"
    // means the same region here as it does in the packing.
    let ring_r = RADIUS * (1.0 - ring_fraction(puck)).sqrt();
    let all = |_: &(f64, f64, f64)| true;
    let in_ring = |d: &(f64, f64, f64)| d.2 >= ring_r;
    let in_core = |d: &(f64, f64, f64)| d.2 < ring_r;

    Traced {
        name,
        frames,
        released,
        in_head: parcels
            .iter()
            .filter(|p| p.phase == Phase::Headspace)
            .count(),
        in_bed: parcels.iter().filter(|p| p.phase == Phase::Bed).count(),
        in_jet: parcels.iter().filter(|p| p.phase == Phase::Jet).count(),
        arrived: done.len(),
        mean_transit: mean(&all, &|d| d.0),
        ring_transit: mean(&in_ring, &|d| d.0),
        core_transit: mean(&in_core, &|d| d.0),
        mean_load: mean(&all, &|d| d.1),
        ring_load: mean(&in_ring, &|d| d.1),
        max_radius,
    }
}

/// One substep for one parcel.
fn advance(p: &mut Parcel, puck: &Puck, now: f64, h: f64, superficial: f64, floor: f64) {
    match p.phase {
        Phase::Headspace => {
            p.at[1] += superficial * h;
            if p.at[1] >= 0.0 {
                p.at[1] = 0.0;
                p.phase = Phase::Bed;
                p.entered = Some(now);
            }
        }
        Phase::Bed => {
            // Midpoint, in the field the pressure solve produced. The field is piecewise
            // constant per cell, so the midpoint buys accuracy at a cell boundary and nothing
            // in the interior -- which is the honest description of what a first-order field
            // can support.
            let v0 = velocity(puck, p.at);
            let mid = [
                p.at[0] + 0.5 * h * v0[0],
                p.at[1] + 0.5 * h * v0[1],
                p.at[2] + 0.5 * h * v0[2],
            ];
            let v = velocity(puck, mid);
            for (at, vel) in p.at.iter_mut().zip(v) {
                *at += h * vel;
            }
            // What it picked up on the way.
            {
                let clamp = |v: f64, n: usize| ((v / DX).floor().max(0.0) as usize).min(n - 1);
                let c = puck
                    .concentration_at(clamp(p.at[0], NX), clamp(p.at[1], NY), clamp(p.at[2], NZ))
                    .to_si();
                // A parcel equilibrates with the pore liquid it is sitting in, over the time it
                // spends in a cell. One cell is `dx / |v|` seconds, so this is that relaxation.
                let tau = (DX / v_len(v).max(1e-9)).max(1e-6);
                p.load += (c - p.load) * (1.0 - (-h / tau).exp());
            }
            if p.at[1] >= floor {
                p.at[1] = floor;
                p.phase = Phase::Jet;
                p.through = Some(now);
            }
        }
        Phase::Jet => {
            // Drawn, not solved: a free jet under gravity, converging on the spout.
            let (cx, cz) = centre();
            p.fall += 9.81 * h;
            p.at[1] += p.fall * h;
            let pull = (0.6 * h * 9.81).min(1.0);
            p.at[0] += (cx - p.at[0]) * pull;
            p.at[2] += (cz - p.at[2]) * pull;
            if p.at[1] > floor + 0.075 {
                p.phase = Phase::Cup;
            }
        }
        Phase::Cup => {}
    }
}

fn v_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// The pore velocity at a point, in the puck's coordinates.
///
/// **Clamped to the grid, not zeroed outside it**, and that distinction cost an afternoon. The
/// midpoint of an RK2 step lands a little further along than the parcel is, so a parcel one cell
/// from the outlet probes *past* the outlet. Returning zero there makes the boundary an
/// attractor: the parcel stops dead a hair inside the last cell — measured, 19.979 mm of a 20 mm
/// bed — and stays there for the rest of the run.
///
/// Nothing about that looks like a bug. The stream is steady, the colours are right, the flow rate
/// is right, and the water simply never reaches the cup. It was visible only because the example
/// counts what arrives and the count was zero.
fn velocity(puck: &Puck, at: [f64; 3]) -> [f64; 3] {
    let clamp = |v: f64, n: usize| -> usize { ((v / DX).floor().max(0.0) as usize).min(n - 1) };
    let v = puck.pore_velocity_at(clamp(at[0], NX), clamp(at[1], NY), clamp(at[2], NZ));
    [v.x, v.y, v.z]
}

/// How many segments a full circle gets.
///
/// Forty-eight rather than seventy-two, and the reason is document size rather than looks: the
/// hardware is redrawn on every frame because the report indexes panels by position and every
/// frame must carry the same ones. At 36 frames of two baskets that is the wireframe seventy-two
/// times over, and it was two thirds of a nine-megabyte file.
const SEG: usize = 48;

/// The hardware, as paths in drawing coordinates, and a value for each.
///
/// Every one is zero, which puts the machine at the bottom of the colour ramp and the water above
/// it. That is a convention rather than a measurement and the legend says so — but it is the
/// convention that lets the geometry and the flow share one panel, and sharing one panel is what
/// makes them occlude each other properly.
fn hardware_paths() -> (Vec<Vec<[f64; 3]>>, Vec<f64>) {
    let mut runs: Vec<Vec<[f64; 3]>> = Vec::new();
    let depth = NY as f64 * DX;

    // A circle at a height, in drawing coordinates.
    let circle = |r: f64, y: f64, n: usize| -> Vec<[f64; 3]> {
        (0..=n)
            .map(|i| {
                let a = std::f64::consts::TAU * i as f64 / n as f64;
                [r * a.cos(), y, r * a.sin()]
            })
            .collect()
    };

    // --- the shower screen: a disc of holes, 3 mm above the bed -------------------------------
    let screen_y = HEADSPACE;
    runs.push(circle(RADIUS, screen_y, SEG));
    runs.push(circle(RADIUS * 0.97, screen_y, SEG));
    for ring in 1..=3 {
        let r = RADIUS * 0.25 * ring as f64;
        let holes = 6 * ring;
        for hole in 0..holes {
            let a = std::f64::consts::TAU * hole as f64 / holes as f64;
            let (hx, hz) = (r * a.cos(), r * a.sin());
            runs.push(
                (0..=5)
                    .map(|i| {
                        let b = std::f64::consts::TAU * i as f64 / 5.0;
                        [hx + 0.0012 * b.cos(), screen_y, hz + 0.0012 * b.sin()]
                    })
                    .collect(),
            );
        }
    }

    // --- the basket: tapered, with a rim and a perforated floor --------------------------------
    // Real baskets taper by about a millimetre over their depth, and the taper is what makes a
    // puck able to fall out.
    let bottom_r = RADIUS - 1.0e-3;
    for (r, y) in [
        (RADIUS, 0.0),
        (RADIUS - 0.25e-3, -depth * 0.25),
        (RADIUS - 0.5e-3, -depth * 0.5),
        (RADIUS - 0.75e-3, -depth * 0.75),
        (bottom_r, -depth),
    ] {
        runs.push(circle(r, y, SEG));
    }
    for rib in 0..24 {
        let a = std::f64::consts::TAU * rib as f64 / 24.0;
        runs.push(vec![
            [RADIUS * a.cos(), 0.0, RADIUS * a.sin()],
            [bottom_r * a.cos(), -depth, bottom_r * a.sin()],
        ]);
    }
    // The rim, which is what the basket hangs by.
    runs.push(circle(RADIUS + 1.5e-3, 0.0, SEG));
    for rib in 0..24 {
        let a = std::f64::consts::TAU * rib as f64 / 24.0;
        runs.push(vec![
            [RADIUS * a.cos(), 0.0, RADIUS * a.sin()],
            [
                (RADIUS + 1.5e-3) * a.cos(),
                0.0,
                (RADIUS + 1.5e-3) * a.sin(),
            ],
        ]);
    }
    // The floor's holes. A real 58 mm basket has several hundred; this draws a legible subset in
    // the right pattern and at the right size.
    for ring in 1..=3 {
        let r = bottom_r * 0.3 * ring as f64;
        let holes = 8 * ring;
        for hole in 0..holes {
            let a = std::f64::consts::TAU * hole as f64 / holes as f64;
            let (hx, hz) = (r * a.cos(), r * a.sin());
            runs.push(
                (0..=5)
                    .map(|i| {
                        let b = std::f64::consts::TAU * i as f64 / 5.0;
                        [hx + 0.0009 * b.cos(), -depth, hz + 0.0009 * b.sin()]
                    })
                    .collect(),
            );
        }
    }

    // --- the portafilter body and the spout ----------------------------------------------------
    let body_r = RADIUS + 4.5e-3;
    let (skirt, throat, spout) = (-depth - 12e-3, -depth - 30e-3, -depth - 44e-3);
    runs.push(circle(body_r, 0.0, SEG));
    runs.push(circle(body_r, skirt, SEG));
    runs.push(circle(7e-3, throat, SEG / 2));
    runs.push(circle(5.5e-3, spout, SEG / 2));
    for rib in 0..24 {
        let a = std::f64::consts::TAU * rib as f64 / 24.0;
        runs.push(vec![
            [body_r * a.cos(), 0.0, body_r * a.sin()],
            [body_r * a.cos(), skirt, body_r * a.sin()],
            [7e-3 * a.cos(), throat, 7e-3 * a.sin()],
            [5.5e-3 * a.cos(), spout, 5.5e-3 * a.sin()],
        ]);
    }

    // --- the cup ------------------------------------------------------------------------------
    let (cup_y, cup_r) = (-depth - 95e-3, 32e-3);
    runs.push(circle(cup_r, cup_y + 42e-3, SEG));
    runs.push(circle(cup_r * 0.72, cup_y, SEG));
    for rib in 0..20 {
        let a = std::f64::consts::TAU * rib as f64 / 20.0;
        runs.push(vec![
            [cup_r * a.cos(), cup_y + 42e-3, cup_r * a.sin()],
            [cup_r * 0.72 * a.cos(), cup_y, cup_r * 0.72 * a.sin()],
        ]);
    }

    let values = vec![0.0; runs.len()];
    (runs, values)
}
