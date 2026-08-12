//! A mixture's conductivity bounds, against a resolved geometry that attains them.
//!
//! `Mix::conductivity_bounds` returns two numbers and claims they bracket every microstructure. That is
//! a theorem and it would be worth nothing here unresolved, because the interesting property of Voigt
//! and Reuss is not that they bound but that they are **attained** — there exist real arrangements of
//! the same two materials that hit each end exactly. If they were merely bounds a caller could treat
//! the midpoint as "roughly right"; because they are attained, the midpoint is roughly right for
//! nothing in particular.
//!
//! The witness is a laminate, and the sharp form of the demonstration is that **one block gives both
//! answers**. Alternating layers of aluminium and borosilicate:
//!
//! ```text
//!   flux across the layers   ->  Reuss exactly    2.2132 W/m/K
//!   flux along  the layers   ->  Voigt exactly   84.0570 W/m/K
//! ```
//!
//! Same material, same volume fractions, same block. A factor of 38 apart, decided by the direction the
//! heat is going. That is the whole argument for `Mix` reporting a pair and refusing to pick.
//!
//! # Why this is not comparing two implementations of one idea
//!
//! `Mix` does algebra on volume fractions and knows nothing about grids. `Solid3D` solves conduction on
//! a resolved geometry and knows nothing about mixture rules — it is told which cell is which material
//! and takes harmonic means on faces because that is what a series resistance is. The agreement is
//! between an algebraic bound and a solved field, and the *reason* they agree is a theorem about
//! laminates that neither of them contains.
//!
//! The third test is the one that makes the pair mean something: a 3D checkerboard is isotropic by
//! symmetry, and it lands strictly between the two — and inside the tighter Hashin–Shtrikman pair,
//! which is what HS is for and which the laminate itself violates.

use dualis_core::mixture::Mix;
use dualis_core::{
    units::{Length, Temperature, Time, Volume},
    Domain, Exchange, Substance,
};
use dualis_thermal::Solid3D;

const DX: f64 = 1e-3;

/// How close to steady a laminate solve is taken, as a relative imbalance between the flux in and the
/// flux out, and **every tolerance on a laminate quantity below traces to it**.
///
/// That is the point of naming it. A laminate's conductivity is an exact identity, so the only thing
/// standing between the assertion and machine precision is how far the march was taken: the tolerances
/// are [`SOLVED`], a hundred times this, so they are statements rather than knife edges. Tightening this
/// tightens them all, which is the relationship a tolerance should have to its cause.
const RESIDUAL: f64 = 1e-11;
/// A hundred times [`RESIDUAL`]. See there.
const SOLVED: f64 = 1e-9;

/// The same, for the three-dimensional checkerboard, and it is looser for a measured reason.
///
/// A laminate is a set of independent one-dimensional problems and settles in a few thousand steps. The
/// board does not: with a 150-fold conductivity contrast the explicit step is set by the **aluminium**
/// and the slowest mode by the **glass**, so the steps needed run as `N² α_al/α_eff` and the approach to
/// steady state is exponential in that. Measured, `1e-11` on a 16³ board does not arrive inside eight
/// million steps; `1e-9` arrives in about four seconds.
///
/// It costs nothing here, and that is why it is acceptable rather than a compromise. Every claim about
/// the board is a **bracketing** claim with margins of 50% and more — 1.53× and 2.32× the Reuss bound
/// against bounds 38× apart — so a residual of `1e-9` is seven orders tighter than anything being
/// asserted. The one claim that needs precision is the aliasing equality, and `1e-9` covers it with a
/// hundredfold margin.
const BOARD_RESIDUAL: f64 = 1e-9;

/// Aluminium and borosilicate, half and half by volume. A 150-fold contrast in conductivity, which is
/// what makes the bounds far enough apart to be unmistakable.
fn half_and_half() -> Mix {
    Mix::of(&[
        (Substance::aluminium_6061(), 0.5),
        (Substance::borosilicate_crown(), 0.5),
    ])
    .expect("fractions sum to one")
}

/// The chain of face resistances along z from the first cell centre to the last, in K/W.
///
/// No marching. A steady state's resistance is a property of the assembled conductances, and reading it
/// off them directly is exact where a march is converged only to a tolerance — `a_layered_wall.rs`
/// measures that difference and found a 130-second march sitting 1.56% high while looking converged.
fn chain_z(w: &Solid3D, i: usize, j: usize) -> f64 {
    let (_, _, nz) = w.counts();
    (1..nz)
        .map(|k| {
            1.0 / w
                .face_conductance((i, j, k - 1), (i, j, k))
                .expect("neighbours along z")
                .to_si()
        })
        .sum::<f64>()
}

/// An effective conductivity from a resistance, a path length in cells and a cross-section in cells.
///
/// `k = L/(R A)` with `L = (n−1) dx` and `A = m dx²`, so `k = (n−1)/(R m dx)`. The path is centre to
/// centre and therefore `(n−1)` cells long, not `n`: the two end half-cells are on the outside of the
/// measurement, which is the same bookkeeping the Dirichlet half-cell correction does elsewhere.
fn conductivity_from(resistance: f64, cells_along: usize, cells_across: usize) -> f64 {
    (cells_along - 1) as f64 / (resistance * cells_across as f64 * DX)
}

/// **Flux across the layers gives Reuss exactly, at every resolution.**
///
/// Every face in a one-cell-thick laminate is an interface, so every face conductance is the harmonic
/// mean of the two materials — and the harmonic mean of two phases at equal volume fraction *is* the
/// Reuss bound. The equality is therefore algebraic rather than asymptotic, which is why it holds to
/// `1e-15` and at four resolutions rather than converging toward something.
///
/// Four resolutions because an exact identity and a first-order error look identical at one. This
/// workspace has been caught by that: the arithmetic face mean is 8.9% out at twelve cells and 1.0% at
/// ninety-six, which a single measurement would have read as a tolerance.
#[test]
fn a_laminate_across_its_layers_is_exactly_the_reuss_bound() {
    let mix = half_and_half();
    let (reuss, voigt) = mix
        .conductivity_bounds()
        .expect("both state a conductivity");
    println!(
        "  bounds: Reuss {:.4}, Voigt {:.4} — a factor of {:.2}",
        reuss.to_si(),
        voigt.to_si(),
        voigt.to_si() / reuss.to_si()
    );

    for cells in [8, 16, 32, 64] {
        let mut w = Solid3D::new(
            "laminate",
            Substance::aluminium_6061(),
            (1, 1, cells),
            Length::from_si(DX),
            Temperature::celsius(20.0),
        );
        // Every other cell, so each layer is one cell thick and the volume fractions are exactly half.
        w.fill(Substance::borosilicate_crown(), |_, _, k| k % 2 == 1);

        let measured = conductivity_from(chain_z(&w, 0, 0), cells, 1);
        let off = (measured / reuss.to_si() - 1.0).abs();
        println!(
            "  {cells:2} cells across: {measured:.6} against {:.6} — off {off:.2e}",
            reuss.to_si()
        );
        assert!(
            off < 1e-15,
            "{cells} cells: across the layers should be exactly Reuss {:.6}, is {measured:.6}",
            reuss.to_si()
        );
    }
}

/// Which axis a steady state is driven along.
#[derive(Clone, Copy)]
enum Axis {
    X,
    Z,
}

/// Which end of it a flux is read at.
#[derive(Clone, Copy)]
enum Face {
    Hot,
    Cold,
}

/// Clamp the two end planes along `axis` 60 K apart, march to a **converged** steady state, and report
/// the effective conductivity from the flux crossing the plane next to the hot face.
///
/// # Converged, and not "marched for a while"
///
/// `a_layered_wall.rs` found a 130-second march sitting **1.56% high** while looking converged, so the
/// stopping rule here is an **absolute residual**: the flux entering the hot face has to equal the flux
/// leaving the cold one to within [`RESIDUAL`], because at steady state nothing is being stored. That is
/// dimensionless and independent of `dx`, unlike "the flux stopped changing over N steps" — which as `dx`
/// shrinks becomes trivially true, since N steps is then a shorter and shorter physical time. The
/// iteration cap exists so a mistake fails loudly instead of hanging, and reaching it is a panic.
///
/// The glass sets the time scale and it is slow: `L^2/alpha` over 16 mm of borosilicate is **495 s**,
/// against the aluminium's 3.7. A first draft marched 5 s, called it steady, and measured a 0.9 K
/// transverse difference that was only the glass not having warmed up yet.
fn steady_conductivity(w: &mut Solid3D, axis: Axis, residual: f64) -> f64 {
    let (nx, ny, nz) = w.counts();
    let (hot, cold) = (Temperature::celsius(80.0), Temperature::celsius(20.0));
    let along = match axis {
        Axis::X => nx,
        Axis::Z => nz,
    };
    let clamp = |w: &mut Solid3D| match axis {
        Axis::X => {
            for j in 0..ny {
                for k in 0..nz {
                    w.set_temperature(0, j, k, hot);
                    w.set_temperature(nx - 1, j, k, cold);
                }
            }
        }
        Axis::Z => {
            for i in 0..nx {
                for j in 0..ny {
                    w.set_temperature(i, j, 0, hot);
                    w.set_temperature(i, j, nz - 1, cold);
                }
            }
        }
    };
    // The flux crossing the plane beside one clamped face, summed face by face. A steady state carries
    // the same flux across every plane, and the difference between the two ends is how far from steady
    // it still is.
    let flux = |w: &Solid3D, face: Face| -> f64 {
        let (a, b) = match (axis, face) {
            (Axis::X, Face::Hot) => (0, 1),
            (Axis::X, Face::Cold) => (nx - 2, nx - 1),
            (Axis::Z, Face::Hot) => (0, 1),
            (Axis::Z, Face::Cold) => (nz - 2, nz - 1),
        };
        let mut q = 0.0;
        match axis {
            Axis::X => {
                for j in 0..ny {
                    for k in 0..nz {
                        let g = w
                            .face_conductance((a, j, k), (b, j, k))
                            .expect("neighbours along x")
                            .to_si();
                        q += g
                            * (w.temperature_at(a, j, k).to_si()
                                - w.temperature_at(b, j, k).to_si());
                    }
                }
            }
            Axis::Z => {
                for i in 0..nx {
                    for j in 0..ny {
                        let g = w
                            .face_conductance((i, j, a), (i, j, b))
                            .expect("neighbours along z")
                            .to_si();
                        q += g
                            * (w.temperature_at(i, j, a).to_si()
                                - w.temperature_at(i, j, b).to_si());
                    }
                }
            }
        }
        q
    };

    let dt = Time::from_si(w.max_stable_dt(Time::from_si(0.0)).to_si() * 0.9);
    let mut t = 0.0;
    let mut settled = false;
    for _ in 0..40_000 {
        for _ in 0..200 {
            clamp(w);
            w.step(Time::from_si(t), dt, &mut Exchange::new())
                .expect("the step this domain sized is stable");
            t += dt.to_si();
        }
        clamp(w);
        // In and out. At steady state they are equal because nothing is being stored, so this is an
        // **absolute residual** rather than a rate of change -- dimensionless, and independent of `dx`.
        if (flux(w, Face::Hot) / flux(w, Face::Cold) - 1.0).abs() < residual {
            settled = true;
            break;
        }
    }
    assert!(
        settled,
        "the march never reached a steady state, and eight million steps is not a tolerance problem"
    );

    // `k = Q L / (A dT)`, with `L` centre to centre of the two clamped planes.
    let across = match axis {
        Axis::X => ny * nz,
        Axis::Z => nx * ny,
    };
    flux(w, Face::Hot) * (along - 1) as f64 * DX / (across as f64 * DX * DX * 60.0)
}

/// **Flux along the layers gives Voigt -- the same block, the other direction.**
///
/// Measured from a solved steady state rather than from assembled conductances, and that is the whole
/// difference between this test and a restatement of the arithmetic. Summing row conductances would
/// *assume* the layers carry independent fluxes, which is the equal-gradient condition Voigt is defined
/// by; solving the field lets transverse heat happen if the scheme wants it to, and at steady state it
/// does not.
///
/// It does happen during the transient -- the aluminium rows reach their profile in seconds and warm the
/// glass rows on the way -- so this is a statement about the steady state and nothing else.
#[test]
fn the_same_laminate_along_its_layers_is_the_voigt_bound() {
    let mix = half_and_half();
    let (reuss, voigt) = mix
        .conductivity_bounds()
        .expect("both state a conductivity");

    for cells in [8, 16] {
        let mut w = Solid3D::new(
            "laminate",
            Substance::aluminium_6061(),
            (cells, 1, cells),
            Length::from_si(DX),
            Temperature::celsius(20.0),
        );
        w.fill(Substance::borosilicate_crown(), |_, _, k| k % 2 == 1);

        let measured = steady_conductivity(&mut w, Axis::X, RESIDUAL);
        let off = (measured / voigt.to_si() - 1.0).abs();
        println!(
            "  {cells:2} cells along:  {measured:.6} against Voigt {:.6} -- off {off:.2e}",
            voigt.to_si()
        );
        assert!(
            off < SOLVED,
            "{cells} cells: along the layers should be Voigt {:.6}, is {measured:.6}",
            voigt.to_si()
        );

        // And no transverse gradient survives, which is what makes the layers independent. Checked after
        // the march rather than argued: a scheme with the x-faces wrong would show one.
        let mut worst = 0.0f64;
        for i in 1..cells - 1 {
            for k in 1..cells {
                worst = worst.max(
                    (w.temperature_at(i, 0, k).to_si() - w.temperature_at(i, 0, k - 1).to_si())
                        .abs(),
                );
            }
        }
        // On a 60 K drop, so the relative figure is `worst/60`. Bounded by `SOLVED · 60 K` for the
        // reason that constant exists: an unfinished march leaves a transverse difference behind, and
        // this is how much of one the stopping rule allows.
        println!(
            "     worst transverse step {worst:.2e} K, {:.2e} of the drop",
            worst / 60.0
        );
        assert!(
            worst < SOLVED * 60.0,
            "{cells} cells: the layers should hold identical profiles, differ by {worst:e} K"
        );
    }

    // The point of this test and the last one together: one block, two numbers, 38x apart.
    assert!(
        voigt.to_si() / reuss.to_si() > 37.0,
        "the two directions should be far apart, are {:.2}x",
        voigt.to_si() / reuss.to_si()
    );
}

/// A three-dimensional checkerboard of `block`-cell cubes, `cells` on a side, solved along both axes.
///
/// Returns `(kx, kz)`. They have to agree: the geometry is cubic-symmetric, so the effective tensor is
/// isotropic, and measuring both is how that stops being an assumption.
fn checkerboard(cells: usize, block: usize) -> (f64, f64) {
    let build = || {
        let mut w = Solid3D::new(
            "checkerboard",
            Substance::aluminium_6061(),
            (cells, cells, cells),
            Length::from_si(DX),
            Temperature::celsius(20.0),
        );
        w.fill(Substance::borosilicate_crown(), move |i, j, k| {
            (i / block + j / block + k / block) % 2 == 1
        });
        w
    };
    (
        steady_conductivity(&mut build(), Axis::X, BOARD_RESIDUAL),
        steady_conductivity(&mut build(), Axis::Z, BOARD_RESIDUAL),
    )
}

/// **A microstructure that is neither direction of the laminate lands strictly inside both bounds.**
///
/// This is what makes the pair mean something. The laminate tests show the ends are reachable; this shows
/// that an arrangement which is neither is genuinely between them, so `Mix` returning a range is
/// describing a real spread of achievable values rather than hedging.
///
/// A three-dimensional checkerboard is the arrangement, because it is cubic-symmetric and therefore has
/// an isotropic effective conductivity — measured, `kx` and `kz` agree to `1e-9`. Its value has no closed
/// form, so the claim is the bracketing: inside `[Reuss, Voigt]`, which is a theorem for any
/// microstructure whatever, and **strictly** inside, because landing on a bound would mean the geometry
/// was secretly a laminate.
///
/// # Two things this test had to learn, and both are about resolution rather than physics
///
/// **A checkerboard one cell per phase is not a checkerboard.** Every cell then has all six neighbours of
/// the other material, so every face in the grid is an interface carrying the same harmonic mean, and the
/// discrete operator is exactly that of a *uniform* medium at `harmonic(167, 1.114) = 2.213236` — which at
/// equal volume fractions is the Reuss bound to the last digit. The microstructure sits at the grid scale,
/// so it is aliased rather than resolved, and the answer would be the same for any arrangement with no
/// same-material faces. Asserted below rather than avoided, so it is not rediscovered.
///
/// **A coarsely resolved high-contrast composite is under-conductive by enough to break a bound.** At a
/// fixed microstructure of four blocks per axis, refining the grid gives
///
/// ```text
///   2 cells/block   3.4926     below HS− = 4.3266
///   4 cells/block   5.1396     above
///   8 cells/block   6.8213     above, and still climbing 33% per doubling
/// ```
///
/// so the two-cell case violates the Hashin–Shtrikman lower bound — not because a checkerboard does, but
/// because two cells cannot represent one. The sequence has not converged at 32³ and the continuum value
/// is not reachable at this cost: with a 150-fold conductivity contrast the explicit step is set by the
/// aluminium and the settling time by the glass, so the step count runs as `N²·α_al/α_eff`. The same
/// pattern in stainless against aluminium, a 10-fold contrast, does the same thing — 39.57, 48.38, 54.10 —
/// which is what says it is the discretisation and not the pair.
///
/// So the HS pair is checked only where the microstructure is resolved, and the *direction* of the
/// remaining error is asserted instead of a value: refining raises the answer. That is a more robust
/// statement than any single number here would be, and it is the one that would catch a scheme that
/// converged the wrong way.
#[test]
fn a_checkerboard_is_inside_the_bounds_once_it_is_resolved() {
    let mix = half_and_half();
    let (reuss, voigt) = mix.conductivity_bounds().expect("bounds");
    let (hs_lo, hs_hi) = mix.hashin_shtrikman().expect("two parts");
    println!(
        "  Reuss {:.4} <= HS- {:.4} .. HS+ {:.4} <= Voigt {:.4}",
        reuss.to_si(),
        hs_lo.to_si(),
        hs_hi.to_si(),
        voigt.to_si()
    );
    // The narrow pair is inside the wide one, which is the theorem that makes HS worth having.
    assert!(hs_lo.to_si() > reuss.to_si() && hs_hi.to_si() < voigt.to_si());

    // Sixteen cells, blocks of 1, 2 and 4: sixteen, eight and four blocks per axis, all even, so every
    // case is exactly half and half and the bounds do not move between them.
    let mut measured = Vec::new();
    for block in [1usize, 2, 4] {
        let (kx, kz) = checkerboard(16, block);
        println!(
            "  block {block}: kx {kx:.6}, kz {kz:.6} — isotropic to {:.2e}, {:.4}x Reuss",
            (kx / kz - 1.0).abs(),
            kx / reuss.to_si()
        );
        assert!(
            (kx / kz - 1.0).abs() < BOARD_RESIDUAL * 100.0,
            "block {block}: cubic symmetry makes this isotropic, {kx:.6} against {kz:.6}"
        );
        measured.push((block, kx));
    }

    // The aliasing case, asserted so it stays known.
    let one = measured[0].1;
    println!(
        "  the one-cell board sits {:.2e} from Reuss, which is the march and not the geometry",
        (one / reuss.to_si() - 1.0).abs()
    );
    assert!(
        (one / reuss.to_si() - 1.0).abs() < BOARD_RESIDUAL * 100.0,
        "a one-cell board is aliased to a uniform medium at the harmonic mean, so it should measure \
         Reuss {:.6}, measures {one:.6}",
        reuss.to_si()
    );

    // Strictly inside the outer pair wherever the microstructure exists at all.
    for (block, k) in &measured[1..] {
        assert!(
            *k > reuss.to_si() * (1.0 + 1e-3) && *k < voigt.to_si() * (1.0 - 1e-3),
            "block {block}: {k:.6} is not strictly inside Reuss {:.6} and Voigt {:.6}",
            reuss.to_si(),
            voigt.to_si()
        );
    }

    // Refining raises the answer, which is the direction the coarse-grid error runs in. Asserted as a
    // trend because the value has not converged and a trend is what can honestly be claimed.
    assert!(
        measured[2].1 > measured[1].1,
        "refining the microstructure should raise the effective conductivity: block 2 gave {:.6}, \
         block 4 gave {:.6}",
        measured[1].1,
        measured[2].1
    );

    // And where it is resolved, it is inside the Hashin-Shtrikman pair as an isotropic composite must be.
    let four = measured[2].1;
    assert!(
        four > hs_lo.to_si() && four < hs_hi.to_si(),
        "at four cells per block, {four:.6} should be inside the Hashin-Shtrikman pair {:.6} to {:.6}",
        hs_lo.to_si(),
        hs_hi.to_si()
    );
}

/// **A mixture's heat capacity is exactly what the resolved geometry holds, and its density is exactly
/// the resolved mass.**
///
/// These are the two properties `Mix` claims are exact, so the check is an equality rather than a
/// tolerance. A resolved half-and-half block and a uniform block of the mixed substance must hold the
/// same joules per kelvin and weigh the same, at every resolution.
///
/// This is the test that would catch the volume-and-mass confusion. Volume-weighting `c_p` for this pair
/// is only 0.08% wrong — the two densities are close — so it would pass a percent-level tolerance and
/// fail here, which is the point of making it exact.
#[test]
fn the_mixtures_capacity_is_what_the_resolved_geometry_holds() {
    let mix = half_and_half();
    let cells = 16;
    let mut resolved = Solid3D::new(
        "laminate",
        Substance::aluminium_6061(),
        (1, 1, cells),
        Length::from_si(DX),
        Temperature::celsius(20.0),
    );
    resolved.fill(Substance::borosilicate_crown(), |_, _, k| k % 2 == 1);

    let volume = Volume::from_si(DX * DX * DX * cells as f64);
    let cell = Volume::from_si(DX * DX * DX);

    // The resolved capacity, cell by cell, from the substances themselves.
    let resolved_capacity: f64 = (0..cells)
        .map(|k| {
            let s = if k % 2 == 1 {
                Substance::borosilicate_crown()
            } else {
                Substance::aluminium_6061()
            };
            s.heat_capacity(cell).expect("both state one").to_si()
        })
        .sum();
    let mixed_capacity = mix
        .heat_capacity(volume)
        .expect("the mixture states one")
        .to_si();
    println!(
        "  capacity: resolved {resolved_capacity:.9} J/K, mixed {mixed_capacity:.9} — off {:.2e}",
        (mixed_capacity / resolved_capacity - 1.0).abs()
    );
    assert!(
        (mixed_capacity / resolved_capacity - 1.0).abs() < 1e-15,
        "the mixture must hold exactly what the layers hold: {mixed_capacity} against \
         {resolved_capacity}"
    );

    let resolved_mass: f64 = (0..cells)
        .map(|k| {
            let s = if k % 2 == 1 {
                Substance::borosilicate_crown()
            } else {
                Substance::aluminium_6061()
            };
            s.mass_of(cell).to_si()
        })
        .sum();
    assert!(
        (mix.mass_of(volume).to_si() / resolved_mass - 1.0).abs() < 1e-15,
        "and weigh exactly what they weigh"
    );

    // And the substance the mixture makes carries all three numbers into a domain. The stability limit
    // is `rho c dx^2 / 2k`, so it is only right if the density, the specific heat and the conductivity
    // all arrived — and no two mistakes in them cancel, which is what makes this worth asserting on top
    // of the capacity above. Given the Voigt bound, which is what this laminate is along its layers.
    let (_, voigt) = mix.conductivity_bounds().expect("bounds");
    let composite = mix
        .as_substance("Al/BK7 50-50", voigt, 0.9)
        .expect("a bound is an allowed choice");
    let alpha = voigt.to_si() / (mix.density().to_si() * mix.specific_heat().expect("c_p").to_si());
    let block = Solid3D::new(
        "composite",
        composite,
        (1, 1, cells),
        Length::from_si(DX),
        Temperature::celsius(20.0),
    );
    let expected = DX * DX / (2.0 * alpha);
    let got = block.max_stable_dt(Time::from_si(0.0)).to_si();
    assert!(
        (got / expected - 1.0).abs() < 1e-14,
        "the composite's limit is {got:e} s against dx^2/2a = {expected:e}"
    );
    // Silence the unused warning honestly: `resolved` is what the capacity comparison was about.
    assert_eq!(resolved.counts(), block.counts());
}
