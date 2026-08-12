//! A substance that arrived as text, held to Neumann's exact solution.
//!
//! `any_material.rs` shows a datasheet material working in three domains and surviving a JSON round
//! trip. This file asks the harder question, which is the one somebody with a datasheet actually has:
//! **is it as accurate?**
//!
//! It could fail to be, and not because parsing is unreliable. A scheme's error depends on the
//! material — on the diffusivity through the step it forces, on the Stefan number through how much
//! latent heat sits in the cell holding the front. Every tolerance in this workspace was measured
//! against the catalogue's numbers, and "it works for ice" is not "it works".
//!
//! So: two substances this crate has never heard of, declared as JSON text, marched against the same
//! closed form ice is marched against — **and ice marched beside them through the identical harness**,
//! which is what turns a plausible result into a comparison. They are 82× apart in thermal diffusivity
//! and 21× apart in Stefan number across the sweep, which is the spread that makes the claim mean
//! something: a liquid metal and a wax are as far apart as two things that both freeze get.
//!
//! ```text
//!   X(t) = 2λ√(αt)      where  λ e^{λ²} erf(λ) = St/√π,   St = c(T_m − T_s)/L
//! ```
//!
//! Neumann's solution is the reason this check is possible at all: almost nothing with a moving
//! boundary has a closed form, and this one takes `k`, `ρ`, `c`, `L` and `T_m` as *inputs*. It does not
//! know which substance they belong to, so it cannot be kinder to one the library ships.
//!
//! # The tolerance, and the bound that was 200× too loose
//!
//! A first draft bounded the front error by `dx/X` — one cell of interface resolution, 8.3% with the
//! front at 12 mm and 1 mm cells. Measured, the error is **0.039% at worst**, and a tolerance two
//! hundred times looser than the effect it is meant to catch is not a check: it would have passed a
//! substance whose conductivity was 13% wrong.
//!
//! The reason it is so much better is worth keeping. The front is read from `melted_volume`, an
//! **integral of a conserved quantity**, and the enthalpy scheme conserves that exactly — so the
//! first-order error in the cell holding the interface does not reach the answer. A point sample of
//! temperature has no such protection, which is why the profile tests in `dualis-thermal` measure the
//! order they do and this one does not.
//!
//! The bound here is 0.06%, from a measured worst of 0.039% over two substances and three
//! undercoolings. What it buys is stated in `a_front_this_close_would_notice_a_property_that_is_wrong`
//! rather than assumed: it catches a conductivity that is 0.1% off.

use dualis::prelude::*;
use dualis::thermal::Solid3D;
use dualis::units::ThermalConductivity;
use dualis_core::substance::Substance;

/// Gallium, 99.99%, solid properties.
///
/// Chosen for being nothing like ice and for being unusually well pinned: its melting point is an
/// **ITS-90 fixed point**, 302.9146 K, a defined calibration temperature rather than a measured one.
/// Latent heat 80.16 kJ/kg (5.59 kJ/mol over 69.723 g/mol).
///
/// The conductivity is the polycrystalline figure near 293 K. Single-crystal gallium is strongly
/// anisotropic — about 88, 41 and 16 W/m/K along its three axes — and this domain is isotropic, so the
/// polycrystalline value is the one that belongs in it. Worth stating rather than hiding: the substance
/// modelled here is cast gallium, not a grown crystal.
///
/// No `liquid` block, so the melt sits at the melting point and takes no part. That is the one-phase
/// Stefan problem, which is what Neumann's solution solves and is *exact* rather than approximate for
/// it.
const GALLIUM: &str = r#"{
    "name": "gallium (99.99%), solid",
    "density": 5904.0,
    "thermal": { "conductivity": 40.6, "specific_heat": 371.0,
                 "expansion": 1.8e-5, "emissivity": 0.1 },
    "fusion": { "melting_point": 302.9146, "latent_heat": 80160.0 }
}"#;

/// n-Octadecane, solid properties. What a phase-change thermal buffer is actually made of.
///
/// Melting point 301.3 K (28.15 °C), latent heat 244 kJ/kg, `k` 0.358 W/m/K, `c_p` 1934 J/kg/K,
/// density 814 kg/m³ — the solid-phase column of the values quoted for it in the PCM literature.
///
/// It is here to be gallium's opposite. A wax stores about three times gallium's latent heat per
/// kilogram and conducts about a hundredth as well, which is exactly why it is used as a buffer and
/// exactly what makes it a different numerical problem: 82× the time constant, and a Stefan number
/// 1.7× higher at the same undercooling.
const OCTADECANE: &str = r#"{
    "name": "n-octadecane, solid",
    "density": 814.0,
    "thermal": { "conductivity": 0.358, "specific_heat": 1934.0,
                 "expansion": 8.0e-4, "emissivity": 0.9 },
    "fusion": { "melting_point": 301.3, "latent_heat": 244000.0 }
}"#;

/// `erf`, by its Maclaurin series. Forty terms, and the argument here never exceeds one.
fn erf(x: f64) -> f64 {
    let mut term = x;
    let mut sum = x;
    for n in 1..40 {
        term *= -x * x / n as f64;
        sum += term / (2 * n + 1) as f64;
    }
    sum * 2.0 / std::f64::consts::PI.sqrt()
}

/// Solve `λ e^{λ²} erf(λ) = St/√π` for `λ`, by bisection.
///
/// The left side is zero at zero and increases without bound, so bisection cannot miss it and needs no
/// derivative. Sixty halvings of `[0, 1]` is past what `f64` distinguishes.
fn lambda(stefan: f64) -> f64 {
    let target = stefan / std::f64::consts::PI.sqrt();
    let f = |l: f64| l * (l * l).exp() * erf(l) - target;
    let (mut lo, mut hi) = (0.0, 1.0);
    assert!(f(hi) > 0.0, "the root is not bracketed for St = {stefan}");
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Forty millimetres of column, so a front at eighteen never involves the far face and the column is
/// semi-infinite in the only sense that matters.
const DEPTH: f64 = 40e-3;
const DX: f64 = 1e-3;
/// The measured worst over two substances and three undercoolings is 0.0385%. See the module docs for
/// why this is not `dx/X`, and the teeth test for what it buys.
const BOUND: f64 = 6e-4;

/// What the closed form needs, read off the substance through its public accessors.
///
/// Through the accessors deliberately: a substance that came from a file and one that came from a
/// constructor are then indistinguishable here, which is the property under test.
fn stefan_and_alpha(s: &Substance, undercooling: f64) -> (f64, f64) {
    let volume = Volume::from_si(DX * DX * DX);
    let t = s
        .thermal
        .expect("the declaration states thermal properties");
    let alpha = s
        .diffusivity()
        .expect("so a diffusivity is available")
        .to_si();
    let latent =
        s.latent_energy(volume).expect("as is latent heat").to_si() / s.mass_of(volume).to_si();
    (t.specific_heat.to_si() * undercooling / latent, alpha)
}

/// Freeze a column of `s` inward from one clamped face until the front should have reached `target`,
/// and report where it actually is, where the closed form puts it, and how long it took.
///
/// **Marched to a target depth rather than to a target time**, which is what makes three substances
/// comparable: the front sits at the same place in every run, so the interface is resolved by the same
/// number of cells and the only thing that differs is the material. Equal times would compare a front
/// at 12 mm against one at 1 mm and report the difference as a property of the substance.
///
/// The column starts wholly liquid at exactly its melting point, which is Neumann's initial condition,
/// and is only sayable because `set_melted_fraction` exists: starting a hair above the melting point
/// instead would add `ε·c_p` of sensible heat the closed form does not have, and the front would arrive
/// early by an amount nothing explains.
fn freeze_until(s: &Substance, undercooling: f64, target: f64) -> (f64, f64, f64) {
    let (stefan, alpha) = stefan_and_alpha(s, undercooling);
    let l = lambda(stefan);
    let want_at = (target / (2.0 * l)).powi(2) / alpha;
    let (frozen, t) = freeze_for(s, undercooling, want_at);
    // The closed form at the time actually *reached*, not the time asked for. The march overshoots by
    // up to one step and at these step sizes that is 0.4% of `t` — comparable to everything being
    // measured, so it has to be taken out rather than tolerated.
    (frozen, 2.0 * l * (alpha * t).sqrt(), t)
}

/// Freeze for a given number of seconds, reporting the front and the time actually reached.
///
/// Split out from [`freeze_until`] because the two callers need opposite things, and conflating them
/// is a mistake this file made and caught. Asking for a *depth* makes three substances comparable;
/// asking for a *time* is the only way to see what a wrong property does, because a depth target
/// recomputes itself from the wrong property and the error cancels out exactly. The first draft of the
/// teeth test below reported a 1% conductivity error as moving the front by 0.0000%.
fn freeze_for(s: &Substance, undercooling: f64, want_at: f64) -> (f64, f64) {
    let point = s.fusion.expect("fusion properties").melting_point.to_si();
    let cells = (DEPTH / DX).round() as usize;
    let mut c = Solid3D::new(
        s.name.clone(),
        s.clone(),
        (1, 1, cells),
        Length::from_si(DX),
        Temperature::from_si(point),
    );
    for k in 0..cells {
        c.set_melted_fraction(0, 0, k, 1.0);
    }
    let cold = Temperature::from_si(point - undercooling);
    let dt = Time::from_si(c.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
    let mut t = 0.0;
    while t < want_at {
        c.set_temperature(0, 0, 0, cold);
        c.step(Time::from_si(t), dt, &mut Exchange::new())
            .expect("the step this domain sized is stable");
        t += dt.to_si();
    }
    c.set_temperature(0, 0, 0, cold);

    // From the melted *volume*, so the partial cell the front is inside counts for its fraction and the
    // answer is not quantised to `dx`. The half cell is the Dirichlet boundary: clamping cell zero puts
    // the surface at its centre, so half of that cell lies before `x = 0`.
    (DEPTH - c.melted_volume().to_si() / (DX * DX) - 0.5 * DX, t)
}

/// **A substance declared in a file reaches Neumann's front as closely as the one the crate ships.**
///
/// Six runs on declared substances, three on ice, all with the front at 12 mm so the comparison is
/// between materials and nothing else. Measured:
///
/// ```text
///   substance      5 K         20 K        60 K
///   gallium        0.0015%     0.0056%     0.0287%
///   octadecane     0.0013%     0.0170%     0.0385%
///   ice            0.0015%     0.0108%     0.0345%      <- the catalogue's own
/// ```
///
/// **Ice sits inside the range of the two, at every undercooling**, and that is the whole claim. The
/// error is a function of Stefan number — it grows 25× as `St` goes from 0.02 to 0.48 — and not of
/// where the substance came from. A material this library has never seen is not being treated as a
/// second-class one; there is no class.
///
/// Three undercoolings because one would not show that. The error at 5 K is twenty times smaller than
/// at 60 K, so a single measurement could be read as any accuracy you like by choosing the undercooling.
#[test]
fn a_declared_substance_freezes_where_the_closed_form_says() {
    let ice = Substance::ice();
    let declared: Vec<Substance> = [GALLIUM, OCTADECANE]
        .iter()
        .map(|text| {
            let s: Substance = serde_json::from_str(text).expect("the declaration parses");
            s.check()
                .unwrap_or_else(|e| panic!("and states a possible material: {e}"));
            s
        })
        .collect();

    let mut worst_declared = 0.0f64;
    let mut worst_shipped = 0.0f64;
    let mut stefans: Vec<f64> = Vec::new();
    for s in declared.iter().chain(std::iter::once(&ice)) {
        for undercooling in [5.0, 20.0, 60.0] {
            let (stefan, _) = stefan_and_alpha(s, undercooling);
            let (measured, exact, t) = freeze_until(s, undercooling, 12e-3);
            let error = (measured / exact - 1.0).abs();
            println!(
                "{:28} {undercooling:>4} K  St {stefan:.4}  t {t:8.1} s  {:.5} mm vs {:.5}  {:.4}%",
                s.name,
                measured * 1e3,
                exact * 1e3,
                error * 100.0
            );
            assert!(
                error < BOUND,
                "{}: front {:.5} mm against an exact {:.5} mm is {:.4}%, past the {:.3}% measured \
                 across this sweep — see the module docs for what this bound is and is not",
                s.name,
                measured * 1e3,
                exact * 1e3,
                error * 100.0,
                BOUND * 100.0
            );
            if std::ptr::eq(s, &ice) {
                worst_shipped = worst_shipped.max(error);
            } else {
                worst_declared = worst_declared.max(error);
                stefans.push(stefan);
            }
        }
    }

    // The claim stated as a comparison rather than as two numbers the reader has to compare. Ice is
    // *inside* the declared range, so this is generous to the possibility of failure: if a declared
    // substance were being handled worse, this is where it would show.
    assert!(
        worst_declared < 2.0 * worst_shipped,
        "the declared substances' worst error {:.4}% is more than twice the shipped material's \
         {:.4}%, so something is treating a material from a file differently",
        worst_declared * 100.0,
        worst_shipped * 100.0
    );
    let spread = stefans.iter().fold(0.0f64, |a, b| a.max(*b))
        / stefans.iter().fold(f64::MAX, |a, b| a.min(*b));
    assert!(
        spread > 20.0,
        "the sweep is meant to cover a wide range of Stefan number, covers {spread:.1}×"
    );
    println!(
        "declared worst {:.4}%, shipped worst {:.4}%, over a {spread:.0}× range of Stefan number",
        worst_declared * 100.0,
        worst_shipped * 100.0
    );
}

/// **The bound above is tight enough to catch a conductivity that is a fifth of a percent wrong.**
///
/// The question a tolerance has to answer is not "did it pass" but "what would it have caught", and
/// this is that question made executable. A property error scales the front by `√` of itself, so:
///
/// ```text
///   k off by  0.05%  ->  front moves 0.029%   (0.5× the bound: not caught)
///   k off by  1.00%  ->  front moves 0.520%   (8.7× the bound: caught)
///   k off by  5.00%  ->  front moves 2.602%   (43× the bound: caught)
/// ```
///
/// So the check sees a property error above about **0.1%** — measured at 0.1%, the front moves 0.058%,
/// which is 0.96× the bound and just misses. That case is deliberately not one of the three asserted:
/// an assertion sitting at 96% of its threshold passes today and fails on a rounding change, which
/// makes it a source of noise rather than a statement. The two either side of it are 0.5× and 8.7×.
///
/// The `dx/X` bound the first draft used would have seen 13%, which is to say it would have passed
/// almost any transcription mistake — and it would have *looked* earned, because it traces to a real
/// effect. A tolerance can be derived from the right physics and still be the wrong number.
///
/// Perturbed in the **solver** and not in the closed form, which is the direction that matters: the
/// file says one thing, the exact solution is asked the same question, and the disagreement is what a
/// wrong datasheet actually produces.
#[test]
fn a_front_this_close_would_notice_a_property_that_is_wrong() {
    let truth: Substance = serde_json::from_str(GALLIUM).expect("parses");
    let (reference, _, seconds) = freeze_until(&truth, 20.0, 12e-3);

    for (off, expect_caught) in [(0.0005, false), (0.01, true), (0.05, true)] {
        let mut wrong = truth.clone();
        let mut thermal = wrong.thermal.expect("thermal");
        thermal.conductivity = ThermalConductivity::w_per_m_k(40.6 * (1.0 + off));
        wrong.thermal = Some(thermal);
        // The same number of seconds as the true run, which is the whole point: asking for the same
        // *depth* would recompute the target from the wrong conductivity and the two errors would
        // cancel to the last digit. They did, in the first draft of this test.
        let (measured, _) = freeze_for(&wrong, 20.0, seconds);
        let moved = (measured / reference - 1.0).abs();
        println!(
            "k off by {:>5.2}%  ->  front moves {:.4}%  ({:.1}× the bound)",
            off * 100.0,
            moved * 100.0,
            moved / BOUND
        );
        assert_eq!(
            moved > BOUND,
            expect_caught,
            "a conductivity {:.2}% off moves the front {:.4}%, and the {:.3}% bound {} it",
            off * 100.0,
            moved * 100.0,
            BOUND * 100.0,
            if expect_caught { "misses" } else { "catches" }
        );
    }
}

/// **The `√t` law holds for a declared substance, and there is no material property in it.**
///
/// Nine times the time is three times the depth, whatever the substance and whatever `λ` came out as:
/// `α` and `λ` cancel out of the *ratio* of two fronts. So this is a statement about the scheme with
/// the material divided out, and it is the check a systematically wrong property cannot pass by
/// scaling the answer — scaling does not move a ratio.
///
/// Against `√(t₂/t₁)` and not against 3, because the march overshoots its target by up to one step and
/// the time ratio comes out 8.98 rather than 9.00. Comparing to 3 would charge the scheme 0.1% for the
/// harness's rounding, which is more than everything else measured here put together. Worst over six
/// runs: 0.043%.
#[test]
fn the_square_root_law_survives_a_material_the_crate_never_saw() {
    let mut worst = 0.0f64;
    for text in [GALLIUM, OCTADECANE] {
        let s: Substance = serde_json::from_str(text).expect("the declaration parses");
        for undercooling in [5.0, 20.0, 60.0] {
            let (near, _, t_near) = freeze_until(&s, undercooling, 6e-3);
            let (far, _, t_far) = freeze_until(&s, undercooling, 18e-3);
            let want = (t_far / t_near).sqrt();
            let off = (far / near / want - 1.0).abs();
            println!(
                "{:28} {undercooling:>4} K  {:.4} -> {:.4} mm  {:.5}× against √{:.4} = {want:.5}",
                s.name,
                near * 1e3,
                far * 1e3,
                far / near,
                t_far / t_near
            );
            assert!(
                off < 1e-3,
                "{}: {:.5}× the depth for {:.4}× the time is {:.4}% off √t",
                s.name,
                far / near,
                t_far / t_near,
                off * 100.0
            );
            worst = worst.max(off);
            // And the ratio really is near three, so the exponent is right and not merely consistent.
            assert!(
                (far / near - 3.0).abs() < 0.01,
                "{}: {:.5}× is not 3",
                s.name,
                far / near
            );
        }
    }
    println!("worst departure from √t: {:.4}%", worst * 100.0);
}

/// **A declared substance's stability limit is `dx²/2α` from its own three numbers, to machine
/// precision.**
///
/// The front check is a `1e-4` statement and cannot see a property wrong in the last digits. This one
/// can, and it is the reason to have both: the limit is `ρ c dx²/2k`, so it pins conductivity, density
/// and specific heat *individually* and no two errors in them cancel.
///
/// `2α` and not `6α` because the column has more than one cell along one axis only. This domain sums
/// the actual face conductances, so a bar-shaped block is not charged the three-dimensional rate.
#[test]
fn a_declared_substance_states_its_own_stability_limit() {
    for text in [GALLIUM, OCTADECANE] {
        let s: Substance = serde_json::from_str(text).expect("the declaration parses");
        let t = s.thermal.expect("thermal properties");
        let alpha = t.conductivity.to_si() / (s.density.to_si() * t.specific_heat.to_si());
        let point = s.fusion.expect("fusion").melting_point.to_si();
        let c = Solid3D::new(
            s.name.clone(),
            s.clone(),
            (1, 1, (DEPTH / DX).round() as usize),
            Length::from_si(DX),
            Temperature::from_si(point - 1.0),
        );
        let expected = DX * DX / (2.0 * alpha);
        let got = c.max_stable_dt(Time::from_si(0.0)).to_si();
        assert!(
            (got / expected - 1.0).abs() < 1e-12,
            "{}: limit {got:e} s against dx²/2α = {expected:e}",
            s.name
        );
        // And the substance's own accessor agrees with the arithmetic above, so the two routes to a
        // diffusivity cannot drift apart.
        assert!(
            (s.diffusivity().expect("diffusivity").to_si() / alpha - 1.0).abs() < 1e-15,
            "{}: `diffusivity` disagrees with k/(ρc)",
            s.name
        );
    }
}
