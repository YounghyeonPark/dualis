//! A laminate's three wave speeds, against Backus's exact averages.
//!
//! `Mix::shear_bounds` returns two numbers and claims they bracket every microstructure. As with
//! conductivity, the interesting property is not that they bound but that they are **attained** — and the
//! witness is again a laminate, again giving both ends from one block:
//!
//! ```text
//!   shear across the layers   ->  harmonic mean of mu    2.4517 GPa    (Backus C44)
//!   shear along  the layers   ->  arithmetic mean of mu  13.5945 GPa   (Backus C66)
//! ```
//!
//! Aluminium and PLA, half and half. A factor of **5.5**, decided by which way the shear goes. Same
//! block, same volume fractions.
//!
//! # Backus averaging is the closed form, and it is exact
//!
//! For a medium finely layered perpendicular to `z`, the long-wavelength effective moduli are exact
//! averages of the layers' — Backus 1962. Three of them are simple enough to be unmistakable:
//!
//! ```text
//!   C33 = <1/(lambda+2mu)>^-1     compression through the stack     harmonic M
//!   C44 = <1/mu>^-1               shear across the layers           harmonic mu
//!   C66 = <mu>                    shear in the layer planes         arithmetic mu
//!   rho = <rho>                                                     arithmetic
//! ```
//!
//! `Waves` can measure each separately because [`Waves::hold`] freezes two of the three displacement
//! components, which is what leaves a one-dimensional problem with one modulus in it. The frequency of a
//! clamped mode then gives the speed, and `c² ρ` gives the modulus back.
//!
//! # Why the arithmetic mean appears at all
//!
//! It is worth stating, because a first reading says it should not. For `C66` the wave travels **along**
//! the layers, so each layer could carry its own wave at its own speed and nothing would average. What
//! stops that is the shear stress *between* the layers: if two layers moved differently there would be a
//! `tau_yz` opposing it, and over a wavelength long compared with a layer that coupling wins and forces
//! equal strain. Equal strain is what makes stresses add, which is what makes the mean arithmetic.
//!
//! So `C66` is a **long-wavelength** result in a way `C44` is not, and the two errors in this file are
//! different animals: the element discretisation, which the homogeneous tests measure at second order,
//! and the ratio of layer thickness to wavelength, which only `C66` pays.
//!
//! # What this file costs
//!
//! **68 s in debug, 6.5 s in release**, which is 12% of the workspace's suite for five tests. Marching a
//! wave to a dozen periods at three resolutions is inherently that, and the number is written down here
//! because the first draft cost 5.4× more: the depth of the two propagate-along-the-layers blocks was
//! eight elements where four says the same thing, and the layer-thickness sweep spanned 64 elements where
//! 32 shows the same rate. `dualis-porous` records the same lesson about a jacket at 1 mm — a test that
//! costs more than the eleven beside it is a test somebody will stop running.

use dualis_core::mixture::Mix;
use dualis_core::substance::Substance;
use dualis_core::units::{Length, Time, Velocity};
use dualis_core::{Domain, Exchange};
use dualis_elastic::{Axis, Elastic, Waves};

/// The stiff phase and the compliant one. A 20-fold contrast in shear modulus, which is what makes the
/// two means 5.5 apart and the claim unmistakable — aluminium against borosilicate would be 1.9%, which a
/// mesh error could cover.
fn phases() -> (Elastic, Elastic) {
    (
        Elastic::aluminium_6061(),
        Elastic::from_substance(&Substance::pla()).expect("PLA states its mechanical properties"),
    )
}

/// The two as a `Mix`, half and half by volume, so the bounds come from the kernel rather than from
/// arithmetic written out again beside the assertion.
fn mix() -> Mix {
    Mix::of(&[(Substance::aluminium_6061(), 0.5), (Substance::pla(), 0.5)])
        .expect("fractions sum to one")
}

const DX: f64 = 1e-3;

/// March and return the frequency of the released mode, in hertz, with the leapfrog's dispersion taken
/// out.
///
/// The same method `two_speeds_marched.rs` uses on a homogeneous block, and deliberately unchanged: if
/// the measurement differed the comparison between a laminate and a uniform block would be between two
/// harnesses rather than two materials.
fn measured_frequency(w: &mut Waves, vary: Axis, along: Axis, dt: Time, steps: usize) -> f64 {
    let h = dt.to_si();
    let mut crossings = Vec::new();
    let mut previous = w.mode_amplitude(1, vary, along);
    assert!(
        previous.abs() > 0.0,
        "the mode has to be present before it can be measured"
    );
    for n in 0..steps {
        w.step(Time::from_si(n as f64 * h), dt, &mut Exchange::new())
            .expect("stable");
        let now = w.mode_amplitude(1, vary, along);
        if (previous <= 0.0) != (now <= 0.0) {
            let frac = previous / (previous - now);
            crossings.push((n as f64 + frac) * h);
        }
        previous = now;
    }
    assert!(
        crossings.len() >= 8,
        "not enough oscillation to measure: {} crossings",
        crossings.len()
    );
    let first = crossings[0];
    let last = crossings[crossings.len() - 1];
    let halves = (crossings.len() - 1) as f64;
    let period = 2.0 * (last - first) / halves;
    let omega_discrete = 2.0 * std::f64::consts::PI / period;
    let omega = 2.0 * (omega_discrete * h / 2.0).sin() / h;
    omega / (2.0 * std::f64::consts::PI)
}

/// Measure a speed from a laminate, and report it against `closed`.
///
/// `layers` is how many elements the stack has along its layering axis; the layers are one element thick
/// and alternate, so the volume fractions are exactly half. `span` is the propagation axis's element
/// count.
fn speed_of(
    counts: (usize, usize, usize),
    layers_along: Axis,
    thickness: usize,
    hold: [Axis; 2],
    vary: Axis,
    along: Axis,
    closed: f64,
) -> f64 {
    let (stiff, soft) = phases();
    let mut w = Waves::new("laminate", counts, Length::from_si(DX), stiff);
    let changed = w.fill(soft, |e_x, e_y, e_z| {
        let index = match layers_along {
            Axis::X => e_x,
            Axis::Y => e_y,
            Axis::Z => e_z,
        };
        (index / thickness) % 2 == 1
    });
    let elements = counts.0 * counts.1 * counts.2;
    assert_eq!(
        changed * 2,
        elements,
        "alternating layers should be exactly half of {elements} elements, are {changed}"
    );
    assert_eq!(w.materials().len(), 2, "two materials in the palette");

    for axis in hold {
        w.hold(axis);
    }
    w.clamp_ends(vary);
    w.release_mode(1, vary, along, Length::from_si(1e-9));

    let dt = Time::from_si(w.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    // Ten periods of the expected mode, so the crossing fit has plenty to work with. Sized from the
    // closed form rather than measured, because a step count derived from the answer would shrink
    // silently if the answer went wrong.
    let want = w.mode_frequency(1, vary, Velocity::from_si(closed)).to_si();
    let steps = (10.0 / (want * dt.to_si())).ceil() as usize;
    let f = measured_frequency(&mut w, vary, along, dt, steps);
    // Invert `f = n c / 2L`.
    let span = match vary {
        Axis::X => counts.0,
        Axis::Y => counts.1,
        Axis::Z => counts.2,
    };
    2.0 * f * span as f64 * DX
}

/// **Shear across the layers is the harmonic mean of `μ`; shear along them is the arithmetic mean. Same
/// block.**
///
/// The load-bearing test, and the two ends of `Mix::shear_bounds` measured rather than asserted. A factor
/// of 5.5 between them for a material whose volume fractions never changed — which is the whole reason
/// `Mix` returns a pair and refuses to pick a number.
///
/// Both are checked at three resolutions because the claim is a **rate**: a single mesh would pass for a
/// scheme that happened to be right at that spacing, and this workspace has shipped a first-order scheme
/// claiming second once.
#[test]
fn a_laminate_gives_both_shear_bounds_from_one_block() {
    let m = mix();
    let (reuss, voigt) = m.shear_bounds().expect("both state mechanical properties");
    let rho = m.density().to_si();
    println!(
        "  shear bounds: {:.4} to {:.4} GPa, a factor of {:.3}",
        reuss.to_si() / 1e9,
        voigt.to_si() / 1e9,
        voigt.to_si() / reuss.to_si()
    );
    assert!(
        voigt.to_si() / reuss.to_si() > 5.0,
        "the pair should be far apart for this contrast, is {:.3}x",
        voigt.to_si() / reuss.to_si()
    );

    // C44: propagate along z, polarised x, layers perpendicular to z. The shear stress is continuous
    // through the stack, so the strains add and the mean is harmonic.
    println!("  C44 — shear across the layers, against the harmonic mean:");
    let c44 = converging("C44", (reuss.to_si() / rho).sqrt(), [16, 32, 64], |n| {
        speed_of(
            (1, 1, n),
            Axis::Z,
            1,
            [Axis::Y, Axis::Z],
            Axis::Z,
            Axis::X,
            (reuss.to_si() / rho).sqrt(),
        )
    });

    // C66: propagate along x, polarised y, layers still perpendicular to z. The shear strain is uniform
    // across the layers, so the stresses add and the mean is arithmetic.
    //
    // Four elements of depth and not eight. The claim needs a stack, not a thick one, and the cost of
    // this file is `span × depth × steps` with steps rising with the span — so the depth is the one
    // factor that buys nothing and was halved once the numbers were in. `a_cold_basket…` in
    // `dualis-porous` records the same lesson: a test that costs more than the eleven beside it is a
    // test somebody will stop running.
    println!("  C66 — shear along the layers, against the arithmetic mean:");
    let c66 = converging("C66", (voigt.to_si() / rho).sqrt(), [16, 32, 64], |n| {
        speed_of(
            (n, 1, 4),
            Axis::Z,
            1,
            [Axis::X, Axis::Z],
            Axis::X,
            Axis::Y,
            (voigt.to_si() / rho).sqrt(),
        )
    });

    // The two ends, measured, from the same block. The ratio is the claim.
    let measured = c66 / c44;
    let want = (voigt.to_si() / reuss.to_si()).sqrt();
    println!(
        "  the two speeds are {measured:.5}x apart against sqrt(5.545) = {want:.5} — off {:.3}%",
        (measured / want - 1.0).abs() * 100.0
    );
    assert!(
        (measured / want - 1.0).abs() < 1e-3,
        "the ratio of the two shear speeds should be the square root of the ratio of the bounds: \
         {measured:.5} against {want:.5}"
    );
}

/// Measure `speed` at three spans against `closed`, print the table, and return the finest.
///
/// Two assertions, and they are different claims. The finest mesh is held to **0.1%**, measured at
/// 0.035% to 0.058% across the four moduli — a bound rather than the measurement, because it is the
/// number a caller can rely on. And the fourfold refinement has to improve by at least **12×**, measured
/// 15.4 to 17.2, where second order predicts 16: that is the claim a single mesh cannot make, and this
/// workspace has shipped a first-order scheme calling itself second once.
fn converging(
    what: &str,
    closed: f64,
    spans: [usize; 3],
    mut speed: impl FnMut(usize) -> f64,
) -> f64 {
    let mut errors = Vec::new();
    let mut finest = 0.0;
    for n in spans {
        let got = speed(n);
        let off = (got / closed - 1.0).abs();
        println!(
            "    {n:2} elements: {got:.3} m/s against {closed:.3} — off {:.3}%",
            off * 100.0
        );
        errors.push(off);
        finest = got;
    }
    assert!(
        errors[2] < 1e-3,
        "{what}: the finest mesh should be inside 0.1%, is {:.4}%",
        errors[2] * 100.0
    );
    let improvement = errors[0] / errors[2];
    println!("    fourfold refinement improved it {improvement:.1}x, second order predicts 16");
    assert!(
        improvement > 12.0,
        "{what}: a fourfold refinement should improve a second-order answer about sixteenfold, \
         improved it {improvement:.1}x"
    );
    finest
}

/// **The arithmetic mean is a long-wavelength result, and this is what a thick layer costs.**
///
/// `C66` pays an error the other two do not: it needs the layers to move together, which they only do
/// when a wavelength spans many of them. The sweep above cannot see that error on its own, because
/// refining the propagation span shrinks the element size *and* stretches the wavelength relative to the
/// layer at the same time — two second-order errors falling together, reported as one.
///
/// So this separates them. The span is fixed at 64 elements and the element size never changes; only the
/// layers get thicker, from one element to eight. Anything that moves is the Backus limit and nothing
/// else.
///
/// It matters, and second order in the ratio. Measured, at a span of 32 elements:
///
/// ```text
///   layers per wavelength    64      32      16       8
///   error                 0.234%  1.055%  4.877%  23.296%
/// ```
///
/// Quadrupling the thickness grew the error 22×, where second order predicts 16 — a little steeper,
/// because by eight layers per wavelength the expansion Backus averaging is the leading term of has
/// stopped being a good one. **The engineering number is the last column**: a caller laminating something
/// at eight layers per wavelength is 23% out, and that is the closed form's limit rather than the
/// solver's. Sixteen layers costs 4.9% and sixty-four costs 0.23%.
#[test]
fn the_arithmetic_mean_needs_the_layers_thin_against_the_wavelength() {
    let m = mix();
    let (_, voigt) = m.shear_bounds().expect("mechanical properties");
    let closed = (voigt.to_si() / m.density().to_si()).sqrt();
    // 32 elements of span, and 16 of depth so every thickness divides it evenly into an even number of
    // layers — which is what keeps the volume fractions exactly half at every thickness.
    let mut errors = Vec::new();
    for thickness in [1usize, 2, 4, 8] {
        let got = speed_of(
            (32, 1, 16),
            Axis::Z,
            thickness,
            [Axis::X, Axis::Z],
            Axis::X,
            Axis::Y,
            closed,
        );
        let off = (got / closed - 1.0).abs();
        let per_wavelength = 2.0 * 32.0 / thickness as f64;
        println!(
            "  layers {thickness} element(s) thick, {per_wavelength:.0} per wavelength: {got:.3} m/s \
             — off {:.3}%",
            off * 100.0
        );
        errors.push(off);
    }
    // Second order in layer thickness, and asserted between the second and fourth entries rather than
    // from the first: the first point is partly mesh error rather than Backus error, so including it
    // would flatter the rate.
    let growth = errors[3] / errors[1];
    println!(
        "  quadrupling the layer thickness grew the error {growth:.1}x, second order predicts 16"
    );
    assert!(
        growth > 12.0,
        "the Backus error is second order in layer thickness over wavelength; quadrupling the \
         thickness grew it {growth:.1}x"
    );
    // A baseline so the whole sweep cannot drift upward unnoticed. `0.5%` and not `0.1%`: at a span of
    // 32 elements one-element layers are 64 to a wavelength and the Backus error is **0.234%**, where the
    // 0.058% the C66 test reports is at a span of 64 and so twice as many layers to a wavelength. The
    // absolute accuracy claim belongs to that test, which sweeps resolution; this one is about the rate.
    assert!(
        errors[0] < 5e-3,
        "at one element per layer this should still be inside 0.5%, is {:.4}%",
        errors[0] * 100.0
    );
}

/// **Compression through the stack is the harmonic mean of `λ+2μ`.**
///
/// Backus's `C33`, and the low end of `Mix::p_wave_modulus_bounds`. Holding `x` and `y` at zero leaves
/// `ρü_z = M u_z''` with `M` the constrained modulus, and through a stack the normal stress is continuous
/// so the compliances add.
#[test]
fn compression_through_the_stack_is_the_harmonic_p_wave_modulus() {
    let m = mix();
    let (reuss, voigt) = m
        .p_wave_modulus_bounds()
        .expect("both state mechanical properties");
    let rho = m.density().to_si();
    println!(
        "  P-wave modulus bounds: {:.4} to {:.4} GPa",
        reuss.to_si() / 1e9,
        voigt.to_si() / 1e9
    );
    converging("C33", (reuss.to_si() / rho).sqrt(), [16, 32, 64], |n| {
        speed_of(
            (1, 1, n),
            Axis::Z,
            1,
            [Axis::X, Axis::Y],
            Axis::Z,
            Axis::Z,
            (reuss.to_si() / rho).sqrt(),
        )
    });
}

/// **A block filled with what it already held is the same block, bit for bit.**
///
/// The compatibility statement, and it can be an exact equality because an equal material resolves to the
/// palette entry it already has rather than adding a second one — so the assembly is the same assembly
/// and not merely an equivalent one.
#[test]
fn filling_with_what_it_already_held_changes_nothing() {
    let stiff = Elastic::aluminium_6061();
    let mut w = Waves::new("column", (1, 1, 16), Length::from_si(DX), stiff);
    w.hold(Axis::X);
    w.hold(Axis::Y);
    w.clamp_ends(Axis::Z);
    w.release_mode(1, Axis::Z, Axis::Z, Length::from_si(1e-9));
    let limit = w.max_stable_dt(Time::from_si(0.0));

    let changed = w.fill(stiff, |_, _, _| true);
    assert_eq!(
        changed, 0,
        "nothing changed, so nothing is reported changed"
    );
    assert_eq!(w.materials().len(), 1, "and the palette did not grow");
    assert_eq!(
        w.max_stable_dt(Time::from_si(0.0)).to_si(),
        limit.to_si(),
        "the stability limit is identical, not merely close"
    );

    // A partial fill with the same material is also a no-op, which is the case a `!=` check could get
    // wrong by counting every element it visited.
    assert_eq!(w.fill(stiff, |_, _, e_z| e_z % 2 == 0), 0);
}

/// **Compression *along* the layers, with the lateral strain held, is the arithmetic mean of `λ+2μ`.**
///
/// So **both** ends of [`Mix::p_wave_modulus_bounds`] are attained, and a first draft of this file
/// documented the high end as a bound that is not — which was wrong, and the measurement is what said
/// so. Holding `y` and `z` at zero on every node makes the lateral strain zero *pointwise* rather than
/// only on average, which is exactly the equal-strain condition Voigt assumes, and the layers then each
/// carry `Mᵢ ε` so the stresses add.
///
/// Longer spans than the shear tests need — 32 to 128 rather than 16 to 64 — because this configuration
/// starts further out: 0.817% at 32 elements against C66's 1.003% at 16, and the reason is the modulus
/// contrast the interlayer coupling has to overcome. Second order all the same — 17.1× over the fourfold
/// refinement — and 0.048% at the finest.
///
/// # What is *not* claimed, and why the API cannot pose it
///
/// A **free** laminate compressed along its layers is a different problem: the layers each want a
/// different lateral contraction, the ones beside them prevent it, and the answer is Backus's `C11` —
/// which for this pair is 43.77 GPa, 18.9% below the 53.98 of the Voigt bound. Measuring it needs the
/// lateral strain zero *macroscopically* while free locally, which means holding `u_z` on the two z faces
/// and nowhere else. [`Waves::hold`] holds a component on every node, so this API cannot express it.
///
/// Tried anyway, and the result is worth recording: releasing `u_z` everywhere gives 40.87 GPa, which is
/// neither `⟨M⟩` nor `C11`. A block four elements thick with traction-free faces carries **plate** modes,
/// whose modulus is neither — so the number is an answer to a third question. No claim is made about
/// `C11` here, and that is a limit of the boundary conditions rather than of the bound.
#[test]
fn compression_along_the_layers_is_the_arithmetic_p_wave_modulus() {
    let m = mix();
    let (_, voigt) = m
        .p_wave_modulus_bounds()
        .expect("both state mechanical properties");
    let rho = m.density().to_si();
    println!(
        "  against the arithmetic mean, {:.4} GPa:",
        voigt.to_si() / 1e9
    );
    converging(
        "C11-constrained",
        (voigt.to_si() / rho).sqrt(),
        [32, 64, 128],
        |n| {
            speed_of(
                (n, 1, 4),
                Axis::Z,
                1,
                [Axis::Y, Axis::Z],
                Axis::X,
                Axis::X,
                (voigt.to_si() / rho).sqrt(),
            )
        },
    );
}
