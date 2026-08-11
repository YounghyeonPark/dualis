//! A block of more than one material, against the resistances its layers add up to.
//!
//! `Solid3D` was one material for ten domains' worth of history, and [`Solid3D::fill`] is the
//! second half of it rather than a new physics: the same `ρc ∂T/∂t = ∇·(k∇T)`, with `k` and `ρc`
//! allowed to vary from cell to cell. Nothing above the crate changed and the kernel did not.
//!
//! # What is worth checking, and it is one number
//!
//! Heat crosses a **face**, and a face has two materials touching it. The half cell either side is
//! a resistance in series, so the conductivity that governs the face is the harmonic mean
//! `2k_Lk_R/(k_L+k_R)` — and the arithmetic mean is not a coarser convention, it is a short circuit.
//! For aluminium against borosilicate the two are 2.21 and 84.1 W/m/K, a factor of 38.
//!
//! The harmonic mean earns an **equality**, not a tolerance. With the material interface on a cell
//! face, the discrete chain of face resistances is exactly `Σ Lᵢ/(kᵢA)` at every resolution: the
//! interface face contributes `1/(2k_L) + 1/(2k_R)` per unit length, which is precisely the two
//! half cells it stands for. So there is no discretisation error in a layered wall's resistance at
//! all, and the arithmetic mean's is first order — 8.9% at twelve cells, 1.0% at ninety-six, and
//! about a thousand cells to reach in 0.1% what harmonic has at three.
//!
//! # The stability limit stopped being a material property
//!
//! `max_stable_dt` sums the actual face conductances now instead of dividing by a diffusivity, and
//! that changed two answers. A block one cell thick across two axes reports `dx²/(2α)` — a bar's
//! limit, which is what that block is, where it used to report `dx²/(6α)` and cost three times the
//! steps. And a filled block is limited by its worst **cell**, which is usually far looser than its
//! fastest material: one aluminium cell inside borosilicate is stable at 75× aluminium's own
//! `dx²/(6α)`, because heat cannot arrive at a cell faster than its worst face delivers it.

use dualis_core::conserved::quantity;
use dualis_core::{
    units::{Energy, Length, SpecificHeat, Temperature, Time},
    Domain, Exchange, Substance,
};
use dualis_thermal::{Solid3D, STABLE_FOURIER_3D};

const ALU: f64 = 167.0;
const GLASS: f64 = 1.114;

/// A wall of `cells` one-millimetre cubes along z, aluminium then borosilicate, split in the middle.
///
/// One cell across x and y, so it is a bar: a layered wall has no transverse structure, and paying
/// for two more axes to resolve nothing would only make the test slower.
fn wall(cells: usize) -> Solid3D {
    let mut w = Solid3D::new(
        "wall",
        Substance::aluminium_6061(),
        (1, 1, cells),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    );
    w.fill(Substance::borosilicate_crown(), |_, _, k| k >= cells / 2);
    w
}

/// The chain of face resistances from the first cell centre to the last, in K/W.
fn chain(w: &Solid3D) -> f64 {
    let (_, _, nz) = w.counts();
    (1..nz)
        .map(|k| {
            1.0 / w
                .face_conductance((0, 0, k - 1), (0, 0, k))
                .expect("neighbours along z")
                .to_si()
        })
        .sum()
}

/// **A filled block that was filled with what it already held is the same block.**
///
/// The compatibility statement, and it has to be an equality: `fill` replaced the sweep, not just
/// extended it, so a uniform block now goes through the conductance form rather than through
/// `f·(ΣT − 6T)`. If that reformulation moved a uniform answer, every closed form the domain was
/// already checked against would have been quietly reinterpreted.
#[test]
fn filling_a_block_with_its_own_substance_changes_nothing() {
    let build = || {
        let mut b = Solid3D::new(
            "b",
            Substance::aluminium_6061(),
            (5, 4, 3),
            Length::mm(2.0),
            Temperature::celsius(20.0),
        );
        b.deposit(2, 2, 1, Energy::from_si(4.0));
        b
    };
    let (mut plain, mut filled) = (build(), build());
    filled.fill(Substance::aluminium_6061(), |_, _, _| true);
    // Filling with a substance already present must not grow the palette either — otherwise a
    // caller layering a wall in a loop would collect thousands of identical entries.
    assert_eq!(
        filled.substances(),
        1,
        "the same substance is the same entry"
    );

    let dt = Time::from_si(plain.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    for n in 0..40 {
        let t = Time::from_si(n as f64 * dt.to_si());
        plain.step(t, dt, &mut Exchange::new()).expect("stable");
        filled.step(t, dt, &mut Exchange::new()).expect("stable");
    }
    let (a, b) = (
        plain.peak_temperature().to_si(),
        filled.peak_temperature().to_si(),
    );
    assert!(
        (a - b).abs() < 1e-13,
        "a uniform fill is a no-op: {a:.15e} against {b:.15e}"
    );
    assert!((plain.mean_temperature().to_si() - filled.mean_temperature().to_si()).abs() < 1e-13);
}

/// **A layered wall's resistance is exactly what its layers add up to, at every resolution.**
///
/// The load-bearing test. `Σ 1/G_f` over the chain of faces against `Σ Lᵢ/(kᵢA)`, and it is an
/// equality rather than a tolerance because the harmonic mean makes it one: the interface face
/// stands for two half cells and contributes exactly their two resistances.
///
/// Four resolutions, because the whole claim is that this does not converge — it starts correct.
/// A single resolution would pass for a scheme that was accidentally right at that cell count.
#[test]
fn a_layered_wall_has_exactly_the_series_resistance_of_its_layers() {
    let dx = 1e-3;
    let area = dx * dx;
    for cells in [12, 24, 48, 96] {
        let w = wall(cells);
        let half = cells / 2;
        // From the first cell centre to the last. The material interface is the face between cells
        // `half-1` and `half`, so the aluminium path is `(half - 0.5)` cells long and the glass path
        // is `(cells - half - 0.5)`.
        let closed = (half as f64 - 0.5) * dx / (ALU * area)
            + (cells as f64 - half as f64 - 0.5) * dx / (GLASS * area);
        let got = chain(&w);
        let off = (got / closed - 1.0).abs();
        println!("  {cells:3} cells: {got:.6e} K/W against Σ L/kA {closed:.6e} — off {off:.2e}");
        assert!(
            off < 1e-14,
            "series resistance is exact, not approximate: {cells} cells off by {off:.3e}"
        );
    }
}

/// **The arithmetic mean converges, at first order, on what the harmonic mean has straight away.**
///
/// The reason the averaging is not a matter of taste. An arithmetic face at the interface
/// short-circuits it: 84.1 W/m/K where the series answer is 2.21, so the wall comes out with less
/// resistance than its own layers have.
///
/// The error is one face out of `n`, so it falls as `1/n` — which means it *does* vanish, and that
/// is exactly why it is dangerous. A single-resolution check would read as a discretisation error
/// and be given a tolerance. The claim here is the **rate**: halving with each refinement, and
/// still 1% at ninety-six cells where harmonic is at `1e-15`.
#[test]
fn arithmetic_averaging_is_first_order_where_harmonic_averaging_is_exact() {
    let dx = 1e-3;
    let area = dx * dx;
    let mut errors = Vec::new();
    for cells in [12, 24, 48, 96] {
        let half = cells / 2;
        let closed = (half as f64 - 0.5) * dx / (ALU * area)
            + (cells as f64 - half as f64 - 0.5) * dx / (GLASS * area);
        // The same chain of faces, with the one interface face averaged the other way. Everything
        // else is identical, so the whole difference below is that single face.
        let per = |k: f64| dx / (k * area);
        let arithmetic = (half as f64 - 1.0) * per(ALU)
            + (cells as f64 - half as f64 - 1.0) * per(GLASS)
            + per(0.5 * (ALU + GLASS));
        let off = (arithmetic / closed - 1.0).abs();
        println!(
            "  {cells:3} cells: arithmetic is {:.3}% short of the layers",
            off * 100.0
        );
        errors.push(off);
    }
    assert!(
        errors[0] > 0.08,
        "at twelve cells it should be badly wrong: {:.3}%",
        errors[0] * 100.0
    );
    for pair in errors.windows(2) {
        let rate = pair[0] / pair[1];
        assert!(
            (rate - 2.0).abs() < 0.1,
            "first order halves with each refinement: {rate:.3}"
        );
    }
    // And the cell count it would take to reach what harmonic already has.
    let needed = 96.0 * errors[3] / 0.001;
    println!("  0.1% would take about {needed:.0} cells; harmonic is there at twelve");
    assert!(
        errors[3] > 0.005,
        "still half a percent at 96: {:.4}",
        errors[3]
    );
}

/// **The marched steady state lands on that resistance and puts the interface where it belongs.**
///
/// The chain sum above is arithmetic on the conductances; this is the sweep actually realising it.
/// Both ends are clamped every step, which makes the boundary a Dirichlet one at the *cell centre*
/// — so the conducting length is `(n−1)dx` and not `n·dx`, and that is where the closed form is
/// evaluated.
///
/// The interface temperature is the second half and the one an engineer asks for. It is not halfway:
/// aluminium's 150× conductivity means it holds almost no gradient, so the interface sits within a
/// percent of the hot face and the whole drop is across the glass. A scheme that averaged
/// conductivities the wrong way would place it visibly further in.
#[test]
fn the_swept_steady_state_lands_on_the_layers_resistance() {
    let cells = 24;
    let mut w = wall(cells);
    let (hot, cold) = (Temperature::celsius(90.0), Temperature::celsius(20.0));
    let closed = chain(&w);

    let dt = Time::from_si(w.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
    // Long enough, and this took measuring rather than guessing. The slowest transient is the glass
    // half's, `L²/π²α` = 25.9 s over its 11.5 mm, and 130 s of it — five time constants — left the
    // flux **1.56%** high, which is what an unconverged march looks like when it looks converged.
    // Fifteen time constants puts the residual at `e⁻¹⁵`, and the tolerance below is then about the
    // discretisation and not about how long somebody was willing to wait.
    let steps = (400.0 / dt.to_si()) as usize;
    for n in 0..steps {
        w.set_temperature(0, 0, 0, hot);
        w.set_temperature(0, 0, cells - 1, cold);
        w.step(
            Time::from_si(n as f64 * dt.to_si()),
            dt,
            &mut Exchange::new(),
        )
        .expect("stable");
    }
    w.set_temperature(0, 0, 0, hot);
    w.set_temperature(0, 0, cells - 1, cold);

    // At steady state the flux is the same across every face, so measure it on one and compare with
    // the whole wall's ΔT over its whole resistance.
    let drop = hot.to_si() - cold.to_si();
    let face = w
        .face_conductance((0, 0, 5), (0, 0, 6))
        .expect("an aluminium face")
        .to_si();
    let flux = face * (w.temperature_at(0, 0, 5).to_si() - w.temperature_at(0, 0, 6).to_si());
    let want = drop / closed;
    println!(
        "  flux {flux:.6} W against ΔT/R {want:.6} W — off {:.3}%",
        (flux / want - 1.0).abs() * 100.0
    );
    assert!(
        (flux / want - 1.0).abs() < 0.01,
        "the swept wall carries ΔT/R: {flux:.6} against {want:.6}"
    );

    // Every cell centre, against the hot face less the drop over the chain of faces reaching it.
    // Stated for all of them rather than one, because the profile is what a two-layer wall *is*:
    // a kink at the interface, and a scheme that put the kink one cell out would still match at
    // both ends.
    let mut worst: f64 = 0.0;
    let mut running = 0.0;
    for k in 1..cells {
        running += 1.0
            / w.face_conductance((0, 0, k - 1), (0, 0, k))
                .expect("along z")
                .to_si();
        let off = (w.temperature_at(0, 0, k).to_si() - (hot.to_si() - want * running)).abs();
        worst = worst.max(off);
    }
    println!("  and every cell centre sits on ΔT/R to within {worst:.3e} K of 70");
    assert!(
        worst < 1e-3,
        "the whole profile is the chain of resistances: worst {worst:.3e} K"
    );

    // The **material interface** is the face between cells, not either cell's centre, and its
    // temperature is the hot face less the drop over the aluminium alone.
    let alu_path = (cells as f64 / 2.0 - 0.5) * 1e-3 / (ALU * 1e-6);
    let interface = hot.to_si() - want * alu_path;
    let held = (hot.to_si() - interface) / drop;
    println!(
        "  the metal holds {:.3} K of the {drop:.0} K — {:.2}% — so the interface is at {:.3} C",
        hot.to_si() - interface,
        held * 100.0,
        interface - 273.15
    );
    // Equal thicknesses, so the shares are `(1/k)` normalised — `k_g/(k_g + k_a)`, not `k_g/k_a`.
    let share = GLASS / (GLASS + ALU);
    assert!(
        (held - share).abs() < 1e-9,
        "equal thicknesses split the drop as 1/k does: {held:.6} against {share:.6}"
    );

    // Half a cell of the poor conductor holds six times what all of the good one does, which is the
    // whole argument for resolving the glass rather than the metal.
    let half_glass = 0.5e-3 / (GLASS * 1e-6) * want;
    println!(
        "  and half a glass cell holds {:.3} K — {:.1}x the entire aluminium half",
        half_glass,
        half_glass / (hot.to_si() - interface)
    );
    assert!(
        half_glass / (hot.to_si() - interface) > 6.0,
        "the gradient lives in the poor conductor: {:.2}x",
        half_glass / (hot.to_si() - interface)
    );
}

/// **Two materials, and the books still balance to the last bit.**
///
/// A conservation check that a uniform block cannot fail in the same way. The sweep computes a face
/// flux and takes it off one side and adds it to the other, so `Σ CᵢTᵢ` is invariant — but only if
/// the ledger weights each cell by *its own* capacity. A ledger using one capacity for the block
/// would report a leak that was not there, or hide one that was.
#[test]
fn a_two_material_block_conserves_exactly_and_its_capacity_is_a_sum() {
    let mut w = wall(24);
    let joules = 30.0;
    w.deposit(0, 0, 0, Energy::from_si(joules));
    let opening = w
        .ledger()
        .get(quantity::ENERGY)
        .expect("energy is on the books");
    assert!(
        (opening - joules).abs() < 1e-12 * joules,
        "the deposit is the opening balance: {opening:.12} against {joules}"
    );

    let dt = Time::from_si(w.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
    for n in 0..4_000 {
        w.step(
            Time::from_si(n as f64 * dt.to_si()),
            dt,
            &mut Exchange::new(),
        )
        .expect("stable");
    }
    let closing = w
        .ledger()
        .get(quantity::ENERGY)
        .expect("energy is on the books");
    println!("  {joules} J in, {closing:.12} J held after 4000 steps");
    assert!(
        (closing - joules).abs() < 1e-9 * joules,
        "an insulated filled block loses nothing: {closing:.12} against {joules}"
    );
    // The heat that spread out of the aluminium and into the glass is the point of the run — a
    // conservation check on a block that never moved would be a check on nothing.
    assert!(
        w.temperature_at(0, 0, 20).to_si() > w.temperature_at(0, 0, 23).to_si() + 1e-9,
        "the gradient should have reached the glass"
    );

    // And the total capacity is what the two halves hold, not twelve cells of either.
    let dx3 = 1e-9;
    let want = 12.0 * dx3 * (2700.0 * 896.0 + 2510.0 * 858.0);
    let got = w.heat_capacity().to_si();
    assert!(
        (got / want - 1.0).abs() < 1e-12,
        "capacity is a sum over cells: {got:.9e} against {want:.9e}"
    );
}

/// **Heat with no place still goes nowhere in particular, which now means capacity-weighted.**
///
/// A uniform *rise*, not equal joules. Equal joules would warm the glass more than the aluminium
/// and so would say where the heat landed — a claim the plain channel never carried, and one that
/// would look like a hot side.
#[test]
fn placeless_heat_leaves_a_two_material_block_uniform() {
    let mut w = wall(24);
    let mut bus = Exchange::new();
    let joules = 12.0;
    bus.publish(quantity::ENERGY, joules);
    w.step(Time::from_si(0.0), Time::from_si(1e-6), &mut bus)
        .expect("stable");

    let first = w.temperature_at(0, 0, 0).to_si();
    for k in 0..24 {
        let t = w.temperature_at(0, 0, k).to_si();
        assert!(
            (t - first).abs() < 1e-12,
            "cell {k} rose differently: {t} against {first}"
        );
    }
    let rise = first - Temperature::celsius(20.0).to_si();
    let want = joules / w.heat_capacity().to_si();
    assert!(
        (rise / want - 1.0).abs() < 1e-12,
        "the rise is the whole block's capacity: {rise:.12e} against {want:.12e}"
    );
}

/// **The limit is however many axes have cells, and the sharpest mode reaches exactly −1 at it.**
///
/// `1/2`, `1/4`, `1/6` for one, two and three active axes. This domain used to report `dx²/(6α)`
/// for every shape, which is safe and is three times too cautious for a block one cell thick — and
/// a bar-shaped `Solid3D` is exactly the shape the closed-form tests use to check the axes against
/// a one-dimensional answer.
///
/// An equality, because there is a closed form that makes it one: at the reported limit the
/// sharpest representable mode has amplification exactly `−1`, flipping sign every step and never
/// decaying. That is marginal stability, and it is the definition of the limit rather than a
/// symptom of it.
#[test]
fn the_stable_step_is_set_by_the_axes_that_have_cells() {
    let n = 8;
    let alpha = Substance::aluminium_6061().diffusivity().unwrap().to_si();
    let dx = 2e-3;
    for (counts, axes, mode) in [
        ((n, 1, 1), 1.0, (n, 0, 0)),
        ((n, n, 1), 2.0, (n, n, 0)),
        ((n, n, n), 3.0, (n, n, n)),
    ] {
        let b = Solid3D::new(
            "b",
            Substance::aluminium_6061(),
            counts,
            Length::mm(2.0),
            Temperature::celsius(20.0),
        );
        let limit = b.max_stable_dt(Time::from_si(0.0)).to_si();
        let closed = dx * dx / (2.0 * axes * alpha);
        let amp = b.mode_amplification(mode, Time::from_si(limit));
        println!(
            "  {axes:.0} active axes: limit {limit:.9e} s against dx²/{:.0}α {closed:.9e}, \
             sharpest mode amplifies by {amp:.15}",
            2.0 * axes
        );
        assert!(
            (limit / closed - 1.0).abs() < 1e-12,
            "the limit counts the axes that conduct: {limit:.9e} against {closed:.9e}"
        );
        assert!(
            (amp + 1.0).abs() < 1e-12,
            "marginal stability is amplification −1: {amp:.15}"
        );
        assert!(
            (b.stability_ratio(Time::from_si(limit)) - 1.0).abs() < 1e-12,
            "and the ratio the sweep guards on is one there"
        );
    }

    // Zero active axes is the end of the same sentence: a single insulated cell has no neighbour to
    // conduct to and therefore no limit at all. It used to report `dx²/6α` for a block that cannot
    // change, and `Time::INFINITY` is what `substeps_for` reads as "no limit" — one substep, which
    // is the right answer rather than an accident.
    let one = Solid3D::new(
        "one",
        Substance::aluminium_6061(),
        (1, 1, 1),
        Length::mm(2.0),
        Temperature::celsius(20.0),
    );
    assert!(one.max_stable_dt(Time::from_si(0.0)).to_si().is_infinite());
    assert_eq!(one.stability_ratio(Time::from_si(1.0)), 0.0);
}

/// **A fast material inside a slow one is stable at far more than its own diffusivity allows.**
///
/// The saving that comes from summing faces instead of dividing by `α_max`. One aluminium cell
/// surrounded by borosilicate cannot be heated at aluminium's rate, because every face into it is
/// the series mean 2.21 W/m/K rather than 167 — so `Σ k_f` is `6 × 2.21` where `dx²/(6α)` assumes
/// `6 × 167`, and the honest limit is seventy-five times larger.
///
/// The inclusion is **still the block's worst cell**; nothing moved to the glass. What changed is
/// only the conductivity on its faces, so the loosening is exactly `k/k_face` — 167/2.21 — and that
/// identity is the assertion rather than the measured 75.45, because a factor that only ever appears
/// as a decimal is a number nobody can check.
///
/// Then it is *used*: the block is marched a thousand steps at that limit and nothing grows. A
/// looser limit that was wrong would show up as an oscillation that doubles, which is what the
/// peak-against-bounds check is for.
#[test]
fn an_inclusion_is_limited_by_its_faces_and_not_by_its_own_diffusivity() {
    let n = 7;
    let mut b = Solid3D::new(
        "glass with a speck of metal in it",
        Substance::borosilicate_crown(),
        (n, n, n),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    );
    b.fill(Substance::aluminium_6061(), |i, j, k| {
        (i, j, k) == (n / 2, n / 2, n / 2)
    });

    let dx = 1e-3;
    let alpha_max = Substance::aluminium_6061().diffusivity().unwrap().to_si();
    let naive = dx * dx / (6.0 * alpha_max);
    let limit = b.max_stable_dt(Time::from_si(0.0)).to_si();
    println!(
        "  limit {limit:.6e} s against dx²/6α_max {naive:.6e} s — {:.1}x looser",
        limit / naive
    );
    // And the factor is not a magic number: the worst cell is still the aluminium one, and what
    // changed is only the conductivity on its faces. So the loosening is *exactly* `k / k_face`.
    let face = 2.0 * ALU * GLASS / (ALU + GLASS);
    println!(
        "  which is k/k_face = {:.4}/{:.4} = {:.3}",
        ALU,
        face,
        ALU / face
    );
    assert!(
        (limit / naive - ALU / face).abs() < 1e-9 * (ALU / face),
        "the loosening is the ratio of the metal's k to its faces': {:.6} against {:.6}",
        limit / naive,
        ALU / face
    );
    let want = 2700.0 * 896.0 * dx * dx / (6.0 * face);
    assert!(
        (limit / want - 1.0).abs() < 1e-12,
        "the worst cell is the inclusion against six series faces: {limit:.9e} against {want:.9e}"
    );

    // And it really is stable there. Sharp initial data, so the mode the limit is about is excited.
    b.deposit(n / 2, n / 2, n / 2, Energy::from_si(0.02));
    let start = b.peak_temperature().to_si();
    let dt = Time::from_si(limit);
    for m in 0..1_000 {
        b.step(
            Time::from_si(m as f64 * dt.to_si()),
            dt,
            &mut Exchange::new(),
        )
        .expect("at the limit is stable");
    }
    let end = b.peak_temperature().to_si();
    println!(
        "  peak {:.4} K over ambient at the start, {:.4} K after a thousand steps",
        start - 293.15,
        end - 293.15
    );
    assert!(
        end < start && end > 293.15,
        "it should relax, not grow and not undershoot: {start} then {end}"
    );

    // Past the limit is refused, and the violation names the number in the units it is stated in.
    let over = Time::from_si(limit * 1.05);
    let err = b
        .step(Time::from_si(0.0), over, &mut Exchange::new())
        .expect_err("5% past the limit must be refused");
    assert_eq!(err.quantity, "Fourier number");
    assert!(
        (err.after / STABLE_FOURIER_3D - 1.05).abs() < 1e-9,
        "by how much: {}",
        err.after
    );
}

/// **A substance that does not say what it conducts is refused, wherever it is.**
///
/// The silent failure this guards. A `fill` reaching for a substance whose thermal properties are
/// unknown gives that cell no conductivity and no capacity, and `NaN` loses every comparison — so a
/// limit computed with `max` would step straight past it, and the block would fill with `NaN` and
/// report a peak of `NaN` while every conservation check passed vacuously.
#[test]
fn a_cell_of_something_unknown_stops_the_sweep_rather_than_spreading_nan() {
    let mut b = Solid3D::new(
        "b",
        Substance::aluminium_6061(),
        (4, 4, 4),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    );
    assert!(b.max_stable_dt(Time::from_si(0.0)).to_si().is_finite());

    b.fill(
        Substance::bulk("mystery", dualis_core::units::Density::g_per_cm3(2.0)),
        |i, j, k| (i, j, k) == (1, 1, 1),
    );
    let err = b
        .step(
            Time::from_si(0.0),
            Time::from_si(1e-9),
            &mut Exchange::new(),
        )
        .expect_err("one unknown cell is enough to refuse");
    assert_eq!(err.quantity, "substance has no diffusivity");
    for k in 0..4 {
        let t = b.temperature_at(1, 1, k).to_si();
        assert!(t.is_finite(), "nothing was written: cell {k} is {t}");
    }
}

/// **A coating is thin, and the conductivity that matters is across it.**
///
/// The case the harmonic mean exists for, in the form somebody actually has: a thermal interface
/// one cell thick between two metals. Its own resistance is small, but it is in series with
/// everything, and the whole question is whether the model lets it be.
///
/// A hundred-micron layer of borosilicate between two aluminium halves adds `L/kA` — about ninety
/// times the aluminium it displaced. The wall's total resistance should rise by that, exactly, and
/// with an arithmetic mean the layer would be short-circuited into near-invisibility.
#[test]
fn a_one_cell_coating_adds_its_own_resistance_and_no_more() {
    let cells = 21;
    let dx = 1e-4;
    let area = dx * dx;
    let build = |coated: bool| {
        let mut w = Solid3D::new(
            "joint",
            Substance::aluminium_6061(),
            (1, 1, cells),
            Length::from_si(dx),
            Temperature::celsius(20.0),
        );
        if coated {
            w.fill(Substance::borosilicate_crown(), |_, _, k| k == cells / 2);
        }
        w
    };
    let bare = chain(&build(false));
    let coated = chain(&build(true));
    // The coating replaces one cell of aluminium with one of glass, so the change is the difference
    // of the two cells' own resistances — the interface faces either side each carry half of each.
    let want = dx / (GLASS * area) - dx / (ALU * area);
    let added = coated - bare;
    println!(
        "  bare {bare:.4} K/W, coated {coated:.4} K/W — the layer adds {added:.4}, want {want:.4}"
    );
    assert!(
        (added / want - 1.0).abs() < 1e-13,
        "a one-cell layer adds exactly its own resistance: {added:.9e} against {want:.9e}"
    );
    println!("  which is {:.1}x the whole bare joint", added / bare);
}

/// **The catalogue cannot reach the factor of two, and a doctored specific heat can.**
///
/// The bound in [`Solid3D::max_stable_dt`]'s other direction: `k_f ≤ 2·min(k_L, k_R)` means a cell
/// whose neighbours conduct *better* can see `Σ k_f` up to `12kᵢ` where its own diffusivity would
/// predict `6kᵢ`, making the limit twice as tight. Reaching it needs neighbours that conduct better
/// **and** store more, so that the middle cell is still the fastest — and no two real solids do
/// that, because volumetric heat capacity varies by one order of magnitude across all of them where
/// conductivity varies by four.
///
/// So this is a test about the *scheme* and the substances are not materials: borosilicate with a
/// specific heat of 1 J/kg/K and copper with a million. Said plainly, because a reader finding
/// those numbers should know they were chosen to sit on a bound rather than measured.
#[test]
fn the_factor_of_two_needs_a_material_that_does_not_exist() {
    let n = 5;
    let fast = Substance::borosilicate_crown().with_specific_heat(SpecificHeat::j_per_kg_k(1.0));
    let heavy = Substance::copper().with_specific_heat(SpecificHeat::j_per_kg_k(1.0e6));
    assert!(
        fast.diffusivity().unwrap() > heavy.diffusivity().unwrap(),
        "the middle cell has to be the fastest, or the bound is about the wrong cell"
    );

    let mut b = Solid3D::new(
        "bound",
        heavy,
        (n, n, n),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    );
    b.fill(fast.clone(), |i, j, k| (i, j, k) == (n / 2, n / 2, n / 2));

    let dx = 1e-3;
    let alpha_max = fast.diffusivity().unwrap().to_si();
    let naive = dx * dx / (6.0 * alpha_max);
    let limit = b.max_stable_dt(Time::from_si(0.0)).to_si();
    let tighter = naive / limit;
    println!("  limit {limit:.6e} s is {tighter:.4}x tighter than dx²/6α_max {naive:.6e} s");
    assert!(
        tighter > 1.98 && tighter < 2.0,
        "the bound is two and is approached, not passed: {tighter:.6}"
    );
}
