//! The `verify` verb: the measurements a passing audit does not make.
//!
//! Every entry in this battery exists because of an error this workspace actually shipped and
//! the audit passed. A 10 ms coupling window delivered exactly the right joules and read a peak
//! temperature 12% low; a wall-weighting defect sat inside a 1.4% error that looked like
//! discretisation until refinement *halved* it instead of quartering it; and a run that
//! reproduces bit for bit is a claim no other tool can check because no other stack is
//! deterministic. None of those is visible in one run — which is the point of this module: it
//! runs the scene **again**, differently, and reports what moved.
//!
//! What it measures:
//!
//! - **Margins, not verdicts.** The worst conservation drift per quantity against its
//!   tolerance, and each evolving domain's stability limit against the coupling window. A pass
//!   with no margin left is a different fact from a pass.
//! - **Determinism.** The same scene run twice must produce byte-identical frames. A digest
//!   over the run is printed so "this reproduces" is checkable rather than claimed.
//! - **Window sensitivity.** The scene rerun with the coupling window halved. What moves is
//!   the error the window was hiding — the one error class no conservation check can see,
//!   because a ledger has no representation for *when*.
//! - **Resolution sensitivity.** The scene rerun with every discretised domain refined 2×.
//!   What moves is discretisation error, measured rather than assumed.
//!
//! `--deep` adds a third run to each sweep (window quartered, resolution 4×) and measures the
//! **order**: `p = log2(d1/d2)` where `d1` and `d2` are successive differences of the same
//! reading. A first-order coupling error reads `p ≈ 1`; a healthy second-order scheme reads
//! `p ≈ 2`; a second-order scheme with a first-order boundary reads between, falling toward 1
//! as the grid refines — which is exactly the signature that found both acoustic boundary
//! defects.
//!
//! # What the measured order actually is
//!
//! The order of the **reading**, sampling included — not of the scheme in the abstract. A peak
//! read at a cell centre moves with the cell, geometry quantised to a grid (a hall's ceiling, a
//! basket's radius) shifts by a fraction of a cell between resolutions, and both contaminate
//! the measurement toward first order. That is not a defect of the battery: the reading is what
//! a person believes, so its convergence is the one that matters. A contaminated order is a
//! true statement about the number being read.
//!
//! The **window** orders carry one more contaminant: a multirate domain's substep count is a
//! ceiling, so its inner step moves non-smoothly as the window halves, and a subcycled
//! domain's window order can read as almost anything small (a bar's peak measured 0.331 from
//! the ceiling alone). A window order means what it says for the quantities the *window*
//! governs — a coupling taken per frame, a feedback closed between frames — and the report
//! prints a reminder above those lines rather than leaving this paragraph to be found later.
//!
//! # What is refused rather than measured
//!
//! A refinement must state the **same physical problem** on a finer grid, and three things in
//! this format cannot: a [`HotSpot`](crate::HotSpot) is one cell, so refining halves its
//! physical size and divides its energy by eight; a [`DomainSpec::Conductor`]'s blocked cells
//! are one cell each, so a notch shrinks; and a channelled [`DomainSpec::Puck`]'s ring is one
//! cell wide by construction. Refining any of those would compare two different problems and
//! call the difference discretisation error, so the sweep is skipped and the report says which
//! domain and why — loudly, because a check that silently skipped is this workspace's oldest
//! failure shape.

use std::collections::BTreeMap;

use crate::{DomainSpec, Scene, World};
use dualis::core::Reading;
use dualis::prelude::*;

/// FNV-1a over the run's JSON, which is every panel, body and reading of every frame.
///
/// A digest and not a checksum-grade hash: its job is to make "these two runs are identical"
/// and "this run reproduced on another machine" checkable at a glance, not to resist an
/// adversary. Written out here in five lines rather than taken from a crate, because the
/// library's dependency count is a promise and this is not worth a dependency.
fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        h ^= *b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// One instrumented run: what [`World::run`] does, with the margins recorded.
#[derive(Debug)]
pub struct Measured {
    /// The final frame's readings — the scalars each domain chose to report, which are the
    /// observables the sweeps compare.
    pub readings: Vec<Reading>,
    /// FNV-1a digest over the frames' JSON. Two runs of one scene must agree on this exactly,
    /// on every platform, at every optimisation level.
    pub digest: u64,
    /// Per conserved quantity: the worst relative drift any single `advance` showed, and the
    /// tolerance it was judged against — measured with the **whole-simulation** audit's own
    /// scale rule, so these margins are against that check exactly.
    ///
    /// That is one of `advance`'s three refusal points, not all of them. The bus transfer
    /// audit has its own absolute tolerance, and a `books_balance` domain is additionally
    /// checked on its own scale — a small domain can sit an ulp from *that* refusal while
    /// these rows read comfortable, which is the very blindness that made the per-domain
    /// check exist. Neither margin is measured here yet, and the report heading says so
    /// rather than letting "conservation margins" claim more than it covers.
    pub drift: Vec<(String, f64, f64)>,
    /// Quantities that were on a ledger and **never rose above the audit's scale floor**, so no
    /// step ever audited them. The audit's own rule — two denormals are equal for every purpose
    /// a simulation has — and the right rule; but a margins table that simply omitted the row
    /// would read as "every quantity has margin", which is the shape of a check that turned
    /// itself off. A flat orbit's `momentum_z` is the live example.
    pub unaudited: Vec<String>,
    /// Ledger arithmetic that stopped being arithmetic: a non-finite total or drift. The audit
    /// itself cannot see these — `NaN > tol` is false — so a run can pass while poisoned, and a
    /// verification battery is where that blindness stops being mirrored. Each entry becomes a
    /// finding.
    pub anomalies: Vec<String>,
    /// Per evolving domain: the smallest stability limit it reported at any frame boundary, in
    /// seconds, and the largest substep count one advance asked of it.
    pub stability: Vec<(String, f64, u32)>,
}

/// Run a scene the way [`World::run`] does, measuring as it goes.
///
/// The loop is `World::run`'s — advance, close the feedback, capture — with the ledger read
/// before and after each advance and the stability limits read before each. A violation aborts
/// with the kernel's own message, because a run the audit refused has no margins to report.
fn run_measured(scene: &Scene) -> Result<Measured, String> {
    let mut world = World::build(scene.clone())?;
    let dt = Time::from_si(scene.duration_s / scene.frames as f64);
    let placed = world.placements();

    let mut drift: BTreeMap<String, (f64, f64)> = BTreeMap::new();
    let mut seen: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    let mut poisoned: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    let mut anomalies: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    let mut stability: BTreeMap<String, (f64, u32)> = BTreeMap::new();

    let mut frames = Vec::with_capacity(scene.frames + 1);
    frames.push(dualis::scene::capture(&world.sim, &placed));
    for _ in 0..scene.frames {
        for d in world.sim.domains() {
            if d.kind() == Kind::Evolving {
                let limit = d.max_stable_dt(world.sim.time()).to_si();
                let entry = stability
                    .entry(d.name().to_string())
                    .or_insert((f64::INFINITY, 0));
                entry.0 = entry.0.min(limit);
            }
        }

        let before = world.sim.ledger();
        let report = world.sim.advance(dt).map_err(|v| {
            format!(
                "the audit stopped the run at t = {:.4} s: {v}",
                world.sim.time().to_si()
            )
        })?;
        let after = world.sim.ledger();

        for (name, n) in &report.substeps {
            if let Some(entry) = stability.get_mut(name) {
                entry.1 = entry.1.max(*n);
            }
        }

        // The audit's own arithmetic, run on a step it passed: same scale rule, same
        // tolerance lookup, recording the worst instead of refusing over a line.
        let mut names: Vec<&str> = before.quantities().map(|(q, _)| q).collect();
        for (q, _) in after.quantities() {
            if !names.contains(&q) {
                names.push(q);
            }
        }
        names.sort_unstable();
        for name in names {
            seen.insert(name.to_string());
            let b = before.get(name).unwrap_or(0.0);
            let a = after.get(name).unwrap_or(0.0);
            let scale = b
                .abs()
                .max(a.abs())
                .max(before.scale_of(name).unwrap_or(0.0))
                .max(after.scale_of(name).unwrap_or(0.0));
            if scale < 1e-300 {
                continue;
            }
            let rel = (a - b).abs() / scale;
            // The audit itself is blind here: `NaN > tol` is false, so a poisoned ledger
            // passes every step. Recorded as an anomaly rather than folded into the maximum,
            // which `f64::max` would silently discard it from.
            if !rel.is_finite() {
                poisoned.insert(name.to_string());
                anomalies.insert(format!(
                    "the {name} ledger's step change is not a number ({b} to {a}) — arithmetic \
                     upstream is poisoned, and the audit cannot see it because NaN compares \
                     false against every tolerance"
                ));
                continue;
            }
            let tol = world.sim.tolerances().for_quantity(name);
            let entry = drift.entry(name.to_string()).or_insert((0.0, tol));
            entry.0 = entry.0.max(rel);
        }

        world.close_feedback();
        frames.push(dualis::scene::capture(&world.sim, &placed));
    }
    dualis::scene::settle_framing(&mut frames);

    let json = dualis::view::to_json(scene.title.as_str(), &frames);
    let readings = frames
        .last()
        .map(|f| f.readings.clone())
        .unwrap_or_default();
    let unaudited: Vec<String> = seen
        .iter()
        .filter(|q| !drift.contains_key(*q) && !poisoned.contains(*q))
        .cloned()
        .collect();
    Ok(Measured {
        readings,
        digest: fnv1a(json.as_bytes()),
        drift: drift.into_iter().map(|(q, (d, t))| (q, d, t)).collect(),
        unaudited,
        anomalies: anomalies.into_iter().collect(),
        stability: stability.into_iter().map(|(n, (l, s))| (n, l, s)).collect(),
    })
}

/// One reading compared across two runs of the same problem.
#[derive(Debug)]
pub struct Shift {
    /// Which domain reported it.
    pub domain: String,
    /// The reading's own label.
    pub label: String,
    /// Its unit, for the report.
    pub unit: &'static str,
    /// The value in the base run.
    pub base: f64,
    /// The value in the comparison run.
    pub other: f64,
    /// `|other − base|` judged against `max(|base|, |other|)` — zero when both are zero.
    ///
    /// A reading sitting near zero exaggerates this ratio; the absolute pair is beside it in
    /// the report for exactly that case.
    pub relative: f64,
}

/// A measured convergence order, or the honest reason there is none.
#[derive(Debug)]
pub enum Order {
    /// `log2(d1/d2)` of successive differences of one reading.
    Measured(f64),
    /// The differences are at the rounding floor, so the reading has converged below what this
    /// battery can measure — which is a good outcome stated as one, not a number invented from
    /// noise.
    BelowFloor,
    /// The **coarser** pair's difference is at the rounding floor while the finer pair's is
    /// not: the differences grew as the knob refined, which no asymptotic error does. Without
    /// this arm the ratio's logarithm printed `order -inf` — a number fabricated from rounding
    /// wearing the name of a measurement — and it is not [`Order::BelowFloor`], whose
    /// documentation calls it a good outcome; this one says the readings are not in any
    /// regime an order describes.
    NotAsymptotic,
}

/// The scale a reading's rounding actually lives on.
///
/// For most readings that is its magnitude. For a temperature it is not: the workspace's one
/// deliberate non-SI convention displays celsius, so the stored number is a kelvin quantity
/// displaced by 273.15 and its rounding sits at `~300·ε` **regardless of the displayed
/// value**. Judged on its displayed magnitude, a bar cooling through 0 °C reads a 3 mK
/// discretisation shift as 13.8% and an order fabricated from kelvin rounding as "does not
/// converge" — measured before this function existed. Keyed on the unit string because the
/// celsius convention is itself keyed there: `Reading::value`'s documentation states it, so
/// this is reading a contract rather than guessing from a name.
fn representation_scale(unit: &str, magnitude: f64) -> f64 {
    if unit == "C" {
        magnitude.max(273.15)
    } else {
        magnitude
    }
}

/// One sweep: the scene rerun with one knob moved, and what each reading did.
#[derive(Debug)]
pub struct Sweep {
    /// Every reading matched across the runs, in the base run's order.
    pub shifts: Vec<Shift>,
    /// The largest [`Shift::relative`] — the one-number summary of how much the knob was
    /// hiding.
    pub worst: f64,
    /// With `--deep`, the measured order per reading: `(domain, label, order)`. Empty without.
    pub orders: Vec<(String, String, Order)>,
    /// Readings present in one run and not the other — including, under `--deep`, the finest
    /// run, which is compared against the base for presence exactly so a reading cannot vanish
    /// between the middle and finest runs without a trace. Should be empty; every entry
    /// becomes a finding.
    pub unmatched: Vec<String>,
    /// Readings whose comparison stopped being arithmetic: a non-finite value, or two readings
    /// sharing one `(domain, label)` key so that `find` silently matched only the first. Both
    /// are unreachable with today's domains and both become findings if they ever happen,
    /// because each is a way for a row to lie while looking routine.
    pub broken: Vec<String>,
}

/// What became of one sweep — ran, was honestly out of scope, or **failed**.
///
/// The distinction `Skipped`/`Failed` is the whole point of this enum existing. A scene with
/// nothing to refine and a refined run the audit refused both used to fold into one error
/// string rendered as `SKIPPED:`, and the second is about the most alarming thing this battery
/// can learn — a scene that passes at its own settings and fails when a knob moves is a scene
/// one knob away from wrong. A skip is a fact about the scene; a failure is a finding, and it
/// carries the exit code.
#[derive(Debug)]
pub enum SweepOutcome {
    /// The sweep ran and compared.
    Ran(Sweep),
    /// The sweep does not apply to this scene, and this is why.
    Skipped(String),
    /// A run inside the sweep was refused, and this is what the kernel said.
    Failed(String),
}

fn compare(base: &[Reading], other: &[Reading]) -> Sweep {
    let mut shifts = Vec::new();
    let mut unmatched = Vec::new();
    let mut broken = Vec::new();

    // A duplicated key would make `find` silently match only the first holder, so the row
    // would look routine while comparing half the data. No domain emits one today; the check
    // is here because the failure would be invisible by construction if one ever did.
    for side in [base, other] {
        for (i, r) in side.iter().enumerate() {
            if side[..i]
                .iter()
                .any(|s| s.domain == r.domain && s.label == r.label)
            {
                broken.push(format!(
                    "{}/{}: two readings under one name — only the first was compared",
                    r.domain, r.label
                ));
            }
        }
    }

    for r in base {
        match other
            .iter()
            .find(|o| o.domain == r.domain && o.label == r.label)
        {
            Some(o) => {
                if !r.value.is_finite() || !o.value.is_finite() {
                    broken.push(format!(
                        "{}/{}: reads {} against {} — not arithmetic, and a percentage of it \
                         would be a lie",
                        r.domain, r.label, r.value, o.value
                    ));
                    continue;
                }
                let scale = representation_scale(r.unit, r.value.abs().max(o.value.abs()));
                let relative = if scale > 0.0 {
                    (o.value - r.value).abs() / scale
                } else {
                    0.0
                };
                shifts.push(Shift {
                    domain: r.domain.clone(),
                    label: r.label.clone(),
                    unit: r.unit,
                    base: r.value,
                    other: o.value,
                    relative,
                });
            }
            None => unmatched.push(format!("{}/{}", r.domain, r.label)),
        }
    }
    for o in other {
        if !base
            .iter()
            .any(|r| r.domain == o.domain && r.label == o.label)
        {
            unmatched.push(format!("{}/{}", o.domain, o.label));
        }
    }
    let worst = shifts.iter().fold(0.0f64, |m, s| m.max(s.relative));
    Sweep {
        shifts,
        worst,
        orders: Vec::new(),
        unmatched,
        broken,
    }
}

/// The order of each reading from three runs at successive halvings of one knob.
///
/// `d1 = f1 − f2` and `d2 = f2 − f4` are differences of the *same* reading, so the unknown
/// exact answer cancels and no closed form is needed — which is what lets this run on any
/// scene. The floor guards the division: differences below `1e-12` of the reading's
/// representation scale are rounding, and an order computed from rounding would be a number
/// that means nothing wearing the name of one that means a lot.
///
/// The `1e-12` is chosen, not derived, and here is what stands behind it: it is ~4.5e3 ulps,
/// which measured ~200× above the accumulated summation noise of the scenes probed
/// (`bar/mean` differences at 5.6e-15 relative) and ~10⁸ below the smallest real
/// discretisation differences the sweeps produce (~1e-4). A scene whose accumulated rounding
/// approached it from below — many more cells over many more substeps — would defeat it, and
/// the number that would fix it lives with the domain, not here; if that day comes, the floor
/// becomes something the domain reports, the way `Ledger` carries `scale` for the same
/// reason.
fn orders_of(f1: &[Reading], f2: &[Reading], f4: &[Reading]) -> Vec<(String, String, Order)> {
    let mut out = Vec::new();
    for r1 in f1 {
        let find = |rs: &[Reading]| {
            rs.iter()
                .find(|o| o.domain == r1.domain && o.label == r1.label)
                .map(|o| o.value)
        };
        let (Some(v2), Some(v4)) = (find(f2), find(f4)) else {
            continue;
        };
        let scale = representation_scale(r1.unit, r1.value.abs().max(v2.abs()).max(v4.abs()));
        let floor = 1e-12 * scale.max(1e-300);
        let d1 = (r1.value - v2).abs();
        let d2 = (v2 - v4).abs();
        let order = if d1 <= floor && d2 <= floor {
            Order::BelowFloor
        } else if d2 <= floor {
            // d1 above the floor with d2 below it is a reading that stopped moving between
            // the finer pair — converged, but the ratio would be infinite and say nothing.
            Order::BelowFloor
        } else if d1 <= floor {
            // The mirror case is not convergence at all: the coarser pair agreed to the bit
            // and the finer pair moved, which no asymptotic error does. The ratio's logarithm
            // here is -inf, and printing it was the fabrication this function's documentation
            // promises not to commit.
            Order::NotAsymptotic
        } else {
            Order::Measured((d1 / d2).log2())
        };
        out.push((r1.domain.clone(), r1.label.clone(), order));
    }
    out
}

/// The whole battery, ready to render.
#[derive(Debug)]
pub struct Battery {
    /// The base run's measurements.
    pub base: Measured,
    /// Whether a second identical run produced the identical digest.
    pub deterministic: bool,
    /// The coupling window halved (and quartered under `--deep`).
    pub window: SweepOutcome,
    /// Every discretised domain refined 2× (and 4× under `--deep`).
    pub resolution: SweepOutcome,
    /// Scene domains that report no readings at all, so the sweeps have nothing of theirs to
    /// compare. Rendered as a row inside each sweep — "no readings, not swept" — because a
    /// sweep section listing every *other* domain looks complete, and an absence a reader can
    /// see beats one they cannot. The same lesson `main.rs` learned about undrawable domains.
    pub unread: Vec<String>,
    /// Structural defects found. Any entry here is an exit-code failure: a run that is not
    /// deterministic, a staggered schedule past a stability limit, a sweep whose run the audit
    /// refused, a reading that vanished under a sweep, or arithmetic that stopped being
    /// arithmetic.
    pub findings: Vec<String>,
}

/// The evolving domains whose stability limit a non-subcycling schedule would overrun.
///
/// Checked on a **built, unrun** world at `t = 0`, because the hazard it names can keep the
/// base run from finishing at all: a staggered or one-way schedule takes the whole window in
/// one step, and a window past a wave domain's CFL limit amplifies every mode it has. Waiting
/// for the run to measure this would report the crash and not the cause.
fn stability_hazard(scene: &Scene) -> Result<Vec<String>, String> {
    if matches!(scene.schedule, crate::ScheduleSpec::Multirate) {
        return Ok(Vec::new());
    }
    let world = World::build(scene.clone())?;
    let window_s = scene.duration_s / scene.frames as f64;
    let mut hazards = Vec::new();
    for d in world.sim.domains() {
        if d.kind() != Kind::Evolving {
            continue;
        }
        let limit = d.max_stable_dt(Time::ZERO).to_si();
        if limit < window_s {
            hazards.push(format!(
                "{}: stability limit {limit:.3e} s is smaller than the {window_s:.3e} s window, \
                 and the {:?} schedule does not subcycle — this scene is silently unstable; use \
                 multirate or shrink the window",
                d.name(),
                scene.schedule
            ));
        }
    }
    Ok(hazards)
}

/// Run the battery.
///
/// The base scene runs twice (determinism), the window sweep once or twice more, and the
/// resolution sweep once or twice more — five to seven runs, the finest of which can cost
/// `2³ × 4` the base for a three-dimensional multirate domain. That cost is the price of the
/// measurements; a battery that only read the one run could only repeat what the audit said.
pub fn verify(scene: &Scene, deep: bool) -> Result<Battery, String> {
    let mut findings = stability_hazard(scene)?;

    let base = match run_measured(scene) {
        Ok(b) => b,
        // A base run the audit refused still gets the hazard reported, because the hazard is
        // usually *why*: an over-limit staggered window amplifies until the books cannot close.
        Err(e) if !findings.is_empty() => {
            return Err(format!("{e}\n  likely why: {}", findings.join("; ")))
        }
        Err(e) => return Err(e),
    };
    let again = run_measured(scene)?;
    let deterministic = base.digest == again.digest;

    if !deterministic {
        findings.push(format!(
            "two runs of this scene produced different bytes ({:016x} then {:016x}) — \
             something consults a clock, a shared generator, or an unordered reduction",
            base.digest, again.digest
        ));
    }

    findings.extend(base.anomalies.iter().cloned());

    // Which of the scene's domains contributed nothing for the sweeps to compare. A battery
    // that measured nothing must not look like one that measured everything and found it
    // still, so an empty reading set skips both sweeps *as* a finding rather than running
    // scenes to compare zero numbers.
    let unread: Vec<String> = scene
        .domains
        .iter()
        .map(|d| d.name().to_string())
        .filter(|n| !base.readings.iter().any(|r| &r.domain == n))
        .collect();
    let (window, resolution) = if base.readings.is_empty() {
        findings.push(
            "no domain in this scene reports a reading, so the sweeps have nothing to \
             measure — nothing here was verified, and running the sweeps would only have \
             compared zero numbers twice"
                .to_string(),
        );
        let why = "no domain reports a reading — nothing to compare".to_string();
        (
            SweepOutcome::Skipped(why.clone()),
            SweepOutcome::Skipped(why),
        )
    } else {
        (
            window_sweep(scene, &base, deep),
            resolution_sweep(scene, &base, deep),
        )
    };

    for (name, outcome) in [("window", &window), ("resolution", &resolution)] {
        match outcome {
            SweepOutcome::Ran(s) => {
                for missing in &s.unmatched {
                    findings.push(format!(
                        "{missing}: reported in one run of the {name} sweep and not another — \
                         a reading that appears or vanishes with a knob is itself a defect"
                    ));
                }
                for b in &s.broken {
                    findings.push(format!("{name} sweep: {b}"));
                }
            }
            SweepOutcome::Failed(why) => {
                findings.push(format!(
                    "the {name} sweep's run was refused: {why} — a scene that passes at its \
                     own settings and fails when a knob moves is one knob away from wrong, \
                     which is a finding and not a skip"
                ));
            }
            SweepOutcome::Skipped(_) => {}
        }
    }

    Ok(Battery {
        base,
        deterministic,
        window,
        resolution,
        unread,
        findings,
    })
}

/// Fold the finest run into a sweep: measure the orders, and presence-check the finest run
/// against the base so a reading cannot vanish between the middle and finest runs without a
/// trace — `orders_of` skips a reading it cannot find in all three, and a skip with nothing
/// counting it is how a list comes up one row short and looks complete.
fn deepen(sweep: &mut Sweep, base: &Measured, middle: &Measured, finest: &Measured) {
    sweep.orders = orders_of(&base.readings, &middle.readings, &finest.readings);
    let presence = compare(&base.readings, &finest.readings);
    for entry in presence.unmatched {
        if !sweep.unmatched.contains(&entry) {
            sweep.unmatched.push(entry);
        }
    }
    for entry in presence.broken {
        if !sweep.broken.contains(&entry) {
            sweep.broken.push(entry);
        }
    }
    // The finest run's own ledger anomalies. The middle run's are the caller's to carry —
    // splitting it this way rather than "everyone adds everything" is what keeps an anomaly
    // from being reported twice and read as two.
    for a in &finest.anomalies {
        sweep.broken.push(format!("in the finest run: {a}"));
    }
}

fn window_sweep(scene: &Scene, base: &Measured, deep: bool) -> SweepOutcome {
    let halve = |s: &Scene, factor: usize| {
        let mut finer = s.clone();
        finer.frames *= factor;
        finer
    };
    let half = match run_measured(&halve(scene, 2)) {
        Ok(m) => m,
        Err(e) => return SweepOutcome::Failed(format!("with the window halved, {e}")),
    };
    let mut sweep = compare(&base.readings, &half.readings);
    for a in &half.anomalies {
        sweep.broken.push(format!("in the halved-window run: {a}"));
    }
    if deep {
        match run_measured(&halve(scene, 4)) {
            Ok(quarter) => deepen(&mut sweep, base, &half, &quarter),
            Err(e) => return SweepOutcome::Failed(format!("with the window quartered, {e}")),
        }
    }
    SweepOutcome::Ran(sweep)
}

fn resolution_sweep(scene: &Scene, base: &Measured, deep: bool) -> SweepOutcome {
    let twice = match scene.refined() {
        Ok(s) => s,
        Err(why) => return SweepOutcome::Skipped(why),
    };
    let fine = match run_measured(&twice) {
        Ok(m) => m,
        Err(e) => return SweepOutcome::Failed(format!("refined 2x, {e}")),
    };
    let mut sweep = compare(&base.readings, &fine.readings);
    for a in &fine.anomalies {
        sweep.broken.push(format!("in the refined run: {a}"));
    }
    if deep {
        // A second refinement can only refuse if the first could, and the first was checked —
        // but "can only" is an argument, so the refusal path stays a skip rather than a panic.
        let four = match twice.refined() {
            Ok(s) => s,
            Err(why) => return SweepOutcome::Skipped(why),
        };
        match run_measured(&four) {
            Ok(finer) => deepen(&mut sweep, base, &fine, &finer),
            Err(e) => return SweepOutcome::Failed(format!("refined 4x, {e}")),
        }
    }
    SweepOutcome::Ran(sweep)
}

impl Scene {
    /// The same physical problem with every discretised domain refined 2×.
    ///
    /// Domains without a grid pass through unchanged — a heater's tank, an orbit's radii and a
    /// molecular fluid's atom count are the *system*, not a discretisation of one. Domains
    /// whose grid carries one-cell features refuse, naming the feature: see
    /// [`DomainSpec::refined`].
    pub fn refined(&self) -> Result<Scene, String> {
        let mut finer = self.clone();
        let mut refined_any = false;
        for spec in &mut finer.domains {
            if let Some(f) = spec.refined()? {
                *spec = f;
                refined_any = true;
            }
        }
        if !refined_any {
            return Err(
                "nothing in this scene is a discretisation — every domain *is* its own \
                 resolution, so there is no finer statement of the same problem to compare \
                 against"
                    .to_string(),
            );
        }
        Ok(finer)
    }
}

impl DomainSpec {
    /// This domain at doubled resolution, `Ok(None)` if it has no resolution, or the reason
    /// doubling would change the problem rather than refine it.
    ///
    /// The rule for what refines: the physical problem must be identical, only the grid finer.
    /// A bar doubles its cells and halves each exposed face's area, and the beam that lands on
    /// it doubles its faces to match — total area and total watts unchanged. A block doubles
    /// its counts, halves its cell, and doubles its region bounds, which lands every material
    /// boundary on the same physical plane.
    ///
    /// The rule for what refuses: any feature whose physical size is *defined* in cells. A
    /// [`HotSpot`](crate::HotSpot) is one cell of excess temperature, so at half the cell it
    /// holds an eighth of the joules; a conductor's blocked cell is a notch that shrinks; a
    /// puck's channel ring is one cell wide by construction. Each is a different problem at a
    /// different resolution, and a comparison across them would report the difference as
    /// discretisation error — a number that means nothing wearing the name of one that means a
    /// lot.
    pub fn refined(&self) -> Result<Option<DomainSpec>, String> {
        let spec = match self {
            // A room's samples sit on the walls — `width` is `(n − 1)·dx` — so halving the
            // spacing is `2n − 1`, exactly as for a hall below. The height is then quantised
            // to whole cells of the *new* spacing, and for a height not commensurate with the
            // width's cell that quantisation converges at first order; the module
            // documentation names it as the contaminant it is.
            DomainSpec::Room {
                name,
                width_m,
                height_m,
                cells_across,
                release,
            } => {
                // The domain clamps to 3 samples, so "refining" 2 gives 2·2−1 = 3 — the same
                // clamped grid, and the sweep would compare a run to itself and print the
                // cleanest convergence in the report. Refused for the same reason a scene
                // with nothing to refine is.
                if *cells_across < 3 {
                    return Err(format!(
                        "{name}: {cells_across} samples across is clamped to 3 by the domain, \
                         so a refinement would land on the same grid and compare the run to \
                         itself"
                    ));
                }
                Some(DomainSpec::Room {
                    name: name.clone(),
                    width_m: *width_m,
                    height_m: *height_m,
                    cells_across: cells_across * 2 - 1,
                    release: release.clone(),
                })
            }
            // `nodes_across` counts nodes, so halving the spacing is `2n − 1`, not `2n`:
            // the width holds `n − 1` intervals and must hold twice as many, each half the
            // size. `2n` would be a spacing of `(n−1)/(2n−1)` times the old one — close to
            // half, not half, and the residue would read as discretisation error.
            DomainSpec::Hall {
                name,
                width_m,
                height_m,
                depth_m,
                nodes_across,
                mode,
                amplitude_pa,
            } => {
                // The same clamp guard as the room's: the hall clamps to 2 nodes, and
                // 2·1 − 1 = 1 clamps straight back.
                if *nodes_across < 2 {
                    return Err(format!(
                        "{name}: {nodes_across} node(s) across is clamped to 2 by the domain, \
                         so a refinement would land on the same grid and compare the run to \
                         itself"
                    ));
                }
                Some(DomainSpec::Hall {
                    name: name.clone(),
                    width_m: *width_m,
                    height_m: *height_m,
                    depth_m: *depth_m,
                    nodes_across: nodes_across * 2 - 1,
                    mode: *mode,
                    amplitude_pa: *amplitude_pa,
                })
            }
            DomainSpec::Bar {
                name,
                length_mm,
                cells,
                area_mm2,
                initial_c,
                exposes,
            } => Some(DomainSpec::Bar {
                name: name.clone(),
                length_mm: *length_mm,
                cells: cells * 2,
                area_mm2: *area_mm2,
                initial_c: *initial_c,
                exposes: exposes.as_ref().map(|b| crate::Boundary {
                    name: b.name.clone(),
                    face_area_mm2: b.face_area_mm2 / 2.0,
                }),
            }),
            DomainSpec::Beam {
                name,
                onto,
                faces,
                face_area_mm2,
                watts,
                reserve_j,
                waist_fraction,
            } => Some(DomainSpec::Beam {
                name: name.clone(),
                onto: onto.clone(),
                faces: faces * 2,
                face_area_mm2: face_area_mm2 / 2.0,
                watts: *watts,
                reserve_j: *reserve_j,
                waist_fraction: *waist_fraction,
            }),
            DomainSpec::Block {
                name,
                cells,
                cell_mm,
                initial_c,
                material,
                regions,
                hot_spot,
                parts,
                cooling,
            } => {
                if hot_spot.is_some() {
                    return Err(format!(
                        "{name}: a hot spot is one cell, so refining halves its physical size \
                         and divides its energy by eight — a different problem, not a finer \
                         one. Verify this scene without the hot spot, or state the initial \
                         condition as a region"
                    ));
                }
                Some(DomainSpec::Block {
                    name: name.clone(),
                    cells: [cells[0] * 2, cells[1] * 2, cells[2] * 2],
                    cell_mm: cell_mm / 2.0,
                    initial_c: *initial_c,
                    material: material.clone(),
                    regions: regions
                        .iter()
                        .map(|r| crate::Region {
                            material: r.material.clone(),
                            from: [r.from[0] * 2, r.from[1] * 2, r.from[2] * 2],
                            to: [r.to[0] * 2, r.to[1] * 2, r.to[2] * 2],
                        })
                        .collect(),
                    hot_spot: None,
                    // The same meshes, re-voxelised at the finer cell — which is the whole
                    // reason a part is stated as geometry rather than as cells. A feature
                    // thinner than the coarse grid appears at the fine one, and the sweep
                    // reporting that as discretisation error is correct: it *is* the error,
                    // and it is the kind with no other symptom.
                    parts: parts.clone(),
                    // Carried through unchanged, and that is the point of stating a face's
                    // **whole** area rather than a cell's: the same part, exposed the same
                    // way, on a finer grid. A per-cell area would have to be halved here and
                    // the sweep would be comparing two different problems.
                    cooling: cooling.clone(),
                })
            }
            DomainSpec::Conductor {
                name,
                cells,
                cell_mm,
                resistivity_ohm_m,
                volts,
                blocked,
            } => {
                if !blocked.is_empty() {
                    return Err(format!(
                        "{name}: a blocked cell is a one-cell notch, so refining shrinks it — \
                         a different geometry, not a finer one"
                    ));
                }
                Some(DomainSpec::Conductor {
                    name: name.clone(),
                    cells: [cells[0] * 2, cells[1] * 2, cells[2] * 2],
                    cell_mm: cell_mm / 2.0,
                    resistivity_ohm_m: *resistivity_ohm_m,
                    volts: *volts,
                    blocked: Vec::new(),
                })
            }
            DomainSpec::Puck {
                name,
                cells,
                cell_mm,
                radius_mm,
                grind_um,
                porosity,
                bar,
                brew_c,
                wall_c,
                channel_porosity,
            } => {
                if channel_porosity.is_some() {
                    return Err(format!(
                        "{name}: the channel ring is one cell wide by construction, so refining \
                         halves the channel — a different fault, not a finer statement of it"
                    ));
                }
                Some(DomainSpec::Puck {
                    name: name.clone(),
                    cells: [cells[0] * 2, cells[1] * 2, cells[2] * 2],
                    cell_mm: cell_mm / 2.0,
                    radius_mm: *radius_mm,
                    grind_um: *grind_um,
                    porosity: *porosity,
                    bar: *bar,
                    brew_c: *brew_c,
                    wall_c: *wall_c,
                    channel_porosity: None,
                })
            }
            // The system, not a discretisation of one. A heater is a tank, an orbit is its
            // radii, a fluid is its atoms, a network is its nodes; twice as many atoms is a
            // different piece of matter.
            DomainSpec::Heater { .. }
            | DomainSpec::Lump { .. }
            | DomainSpec::Light { .. }
            | DomainSpec::Winding { .. }
            | DomainSpec::Network { .. }
            | DomainSpec::Orbit { .. }
            | DomainSpec::Bounce { .. }
            | DomainSpec::Atoms { .. } => None,
        };
        Ok(spec)
    }
}

impl Battery {
    /// The report, as deterministic text.
    ///
    /// Fixed orderings and fixed formats throughout, because this output is a promise the same
    /// scene prints the same report — and because the digest line is only worth printing from a
    /// tool whose own bytes reproduce.
    pub fn render(&self) -> String {
        use std::fmt::Write;
        let mut out = String::new();

        let _ = writeln!(out, "digest          {:016x}", self.base.digest);
        let _ = writeln!(
            out,
            "determinism     {}",
            if self.deterministic {
                "two runs, identical bytes"
            } else {
                "FAILED: two runs differed — see findings"
            }
        );

        let _ = writeln!(
            out,
            "\nconservation margins (whole-simulation audit, worst step; the bus and \
             per-domain checks have their own tolerances, not measured here)"
        );
        if self.base.drift.is_empty() {
            let _ = writeln!(out, "  nothing on any ledger — no margins to have");
        }
        for (q, worst, tol) in &self.base.drift {
            let headroom = if *worst > 0.0 {
                format!("{:.1} digits in hand", (tol / worst).log10())
            } else {
                "no drift measured".to_string()
            };
            let _ = writeln!(out, "  {q:<10} {worst:.3e} of {tol:.1e} — {headroom}");
        }
        // The row the audit's scale floor would otherwise erase: a quantity that never rose
        // above denormal was never audited, and a table that omitted it would read as "every
        // quantity has margin" — the shape of a check that turned itself off.
        for q in &self.base.unaudited {
            let _ = writeln!(out, "  {q:<10} never above the scale floor — not audited");
        }

        let _ = writeln!(out, "\nstability margins (limit vs the coupling window)");
        if self.base.stability.is_empty() {
            let _ = writeln!(out, "  no evolving domain — nothing has a limit");
        }
        for (name, limit, substeps) in &self.base.stability {
            let _ = writeln!(
                out,
                "  {name:<14} limit {limit:.3e} s, {substeps} substep(s) per window"
            );
        }

        let render_sweep = |out: &mut String, title: &str, sweep: &SweepOutcome, window: bool| {
            let _ = writeln!(out, "\n{title}");
            match sweep {
                SweepOutcome::Ran(s) => {
                    if window && !s.orders.is_empty() {
                        let _ = writeln!(
                            out,
                            "  (a subcycled domain's window order includes substep-ceiling \
                             rounding; it is meaningful for what the window itself governs)"
                        );
                    }
                    if s.shifts.is_empty() {
                        // Reachable only if every matched reading was broken, but a heading
                        // with nothing under it reads as "swept, nothing moved" — the exact
                        // misreading this line exists to prevent.
                        let _ = writeln!(out, "  nothing was compared — see findings");
                    }
                    for shift in &s.shifts {
                        let _ = writeln!(
                            out,
                            "  {:<14} {:<14} {:.6} -> {:.6} {} ({:.3}%)",
                            shift.domain,
                            shift.label,
                            shift.base,
                            shift.other,
                            shift.unit,
                            shift.relative * 100.0
                        );
                    }
                    for (domain, label, order) in &s.orders {
                        match order {
                            Order::Measured(p) => {
                                let _ = writeln!(
                                    out,
                                    "  {domain:<14} {label:<14} converges at order {p:.2}"
                                );
                            }
                            Order::BelowFloor => {
                                let _ = writeln!(
                                    out,
                                    "  {domain:<14} {label:<14} differences at the rounding \
                                     floor — converged below measurement"
                                );
                            }
                            Order::NotAsymptotic => {
                                let _ = writeln!(
                                    out,
                                    "  {domain:<14} {label:<14} differences grew as the knob \
                                     refined — not in a regime an order describes"
                                );
                            }
                        }
                    }
                    // The domains this sweep could not see. Without these rows a section
                    // listing every other domain looks complete, and a reader concludes the
                    // absent one was insensitive when in fact it was never measured.
                    for name in &self.unread {
                        let _ = writeln!(out, "  {name:<14} no readings — not swept");
                    }
                }
                SweepOutcome::Skipped(why) => {
                    let _ = writeln!(out, "  SKIPPED: {why}");
                }
                SweepOutcome::Failed(why) => {
                    let _ = writeln!(out, "  FAILED: {why}");
                }
            }
        };
        render_sweep(
            &mut out,
            "window sensitivity (the window halved; what moved was hiding in it)",
            &self.window,
            true,
        );
        render_sweep(
            &mut out,
            "resolution sensitivity (every grid refined 2x; what moved is discretisation)",
            &self.resolution,
            false,
        );

        if self.findings.is_empty() {
            let _ = writeln!(out, "\nno structural findings");
        } else {
            let _ = writeln!(out, "\nFINDINGS");
            for f in &self.findings {
                let _ = writeln!(out, "  ! {f}");
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The digest is pinned so a platform that rounded differently — or a formatting change
    /// that reordered the JSON — fails here rather than quietly re-stamping every report.
    #[test]
    fn the_digest_is_fnv1a() {
        assert_eq!(fnv1a(b""), 0xcbf29ce484222325);
        assert_eq!(fnv1a(b"a"), 0xaf63dc4c8601ec8c);
    }

    /// The floor guard: an order computed from rounding noise is refused in every arm — both
    /// differences below, the denominator below, and the **numerator** below, which is the arm
    /// that was missing: `order(1.0, 1.0+1e-16, 1.001)` printed "converges at order -inf", a
    /// number fabricated from rounding wearing the name of a measurement.
    #[test]
    fn an_order_is_not_invented_from_the_rounding_floor() {
        let r = |v: f64| Reading::new("d", "x", v, "");
        // Differences at 1e-16 of the scale: below the floor on both.
        let o = orders_of(&[r(1.0)], &[r(1.0 + 1e-16)], &[r(1.0 + 2e-16)]);
        assert!(matches!(o[0].2, Order::BelowFloor));
        // A real coarse difference over a rounding-level fine one: converged, not infinite.
        let o = orders_of(&[r(1.01)], &[r(1.0)], &[r(1.0 + 1e-16)]);
        assert!(matches!(o[0].2, Order::BelowFloor));
        // The mirror: a rounding-level coarse difference under a real fine one is not an
        // order of -inf, it is no order at all.
        let o = orders_of(&[r(1.0)], &[r(1.0 + 1e-16)], &[r(1.001)]);
        assert!(matches!(o[0].2, Order::NotAsymptotic));
        let o = orders_of(&[r(1.0)], &[r(1.0)], &[r(1.001)]);
        assert!(matches!(o[0].2, Order::NotAsymptotic));
        // Clean halving reads first order.
        let o = orders_of(&[r(1.4)], &[r(1.2)], &[r(1.1)]);
        match o[0].2 {
            Order::Measured(p) => assert!((p - 1.0).abs() < 1e-12),
            _ => panic!("a clean ratio was refused"),
        }
    }

    /// A temperature's rounding lives on the kelvin scale whatever celsius it displays: a
    /// 3 mK shift on a bar near 0 C must not read as 14%, and kelvin-rounding differences
    /// must not acquire an order. Measured before the fix: `bar/mean` at a −0.258 C datum
    /// reported a 13.8% resolution shift and `Measured(0.000)` where the 20 C datum reported
    /// 1.6e-4 and BelowFloor — identical physics, different offsets.
    #[test]
    fn a_celsius_reading_is_judged_on_the_kelvin_scale() {
        let c = |v: f64| Reading::new("bar", "mean", v, "C");
        // A 3.2 mK shift near 0 C: relative to the representation, not to the tiny display.
        let s = compare(&[c(0.00035)], &[c(0.00355)]);
        assert!(
            s.shifts[0].relative < 1e-4,
            "a millikelvin read as {:.3}",
            s.shifts[0].relative
        );
        // Kelvin-rounding differences (~6e-14 absolute) sit below the kelvin-scale floor
        // (~2.7e-10) and stay refused rather than becoming an order of zero.
        let o = orders_of(&[c(0.0)], &[c(6e-14)], &[c(1.2e-13)]);
        assert!(matches!(o[0].2, Order::BelowFloor));
    }

    /// The broken channel: a NaN and a duplicated key must each surface as a broken row
    /// rather than becoming a routine-looking shift or vanishing from the comparison.
    #[test]
    fn a_broken_reading_is_reported_rather_than_compared() {
        let r = |d: &str, l: &str, v: f64| Reading::new(d, l, v, "");
        let s = compare(&[r("d", "x", f64::NAN)], &[r("d", "x", 1.0)]);
        assert!(s.shifts.is_empty(), "a NaN was compared as arithmetic");
        assert_eq!(s.broken.len(), 1);

        let s = compare(&[r("d", "x", 1.0), r("d", "x", 2.0)], &[r("d", "x", 1.0)]);
        assert!(
            s.broken
                .iter()
                .any(|b| b.contains("two readings under one name")),
            "{:?}",
            s.broken
        );
    }
}
