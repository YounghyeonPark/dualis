# dualis-core

Physics for simulated worlds: spectral radiometry, surface optics, dispersion,
rigid motion and ray geometry.

This is the layer that does not know what it is being used to simulate. It knows
that light has a wavelength, that a surface divides incident power three ways,
that glass bends light by an amount its refractive index decides, that things
move, and that a random choice has to be reproducible. What you build on top —
a microscope, a camera, a room — is somebody else's problem.

Extracted from [AiryTrace](https://github.com/YounghyeonPark/airytrace), which
is its first consumer.

Units are **millimetres**, **nanometres** and **seconds**, everywhere.

## Two invariants

Everything here holds to two rules, and the tests exist to keep them holding.

**Energy is conserved at every surface.** Reflectance and transmittance are
stored; absorptance is whatever is left over. There is no way to write down a
surface that returns more light than reached it — `SurfaceOptics::validate`
rejects one, and `split` renormalises rather than amplifying if you trace it
anyway.

**Nothing is random.** Every stochastic choice comes from a seeded generator
with no global state and no entropy source, so two runs of the same scene agree
to the last bit, on every platform and under WebAssembly.

## What is in it

| Module | |
| --- | --- |
| `spectrum` | Wavelength-dependent quantities: Planck's law, Gaussians, measured curves, filter bands with real edges and finite blocking |
| `optics` | What a surface does to light — Fresnel from the refractive indices, spectral R/T/A, coatings, diffuse scatter |
| `material` | Refractive index against wavelength (Sellmeier, with a small glass catalogue) and how much light survives the glass (Beer-Lambert) |
| `motion` | Rigid motion and time gating: drift, oscillation, spin, strobe — all closed-form in `t`, so frame 7 does not depend on frame 6 |
| `geometry` | Ray intersections against caps, planes, annuli and cylinders; Snell's law; hexapolar and jittered disc sampling |
| `rng` | A deterministic xorshift64\* generator and the sampling built on it |

## What is deliberately not in it

No scene graph, no meshes, no acceleration structure, no renderer. Those belong
to whatever is built on top, and adding them here before there is a second
consumer would be guessing at an interface rather than discovering one.

## A taste

```rust
use dualis_core::{Material, Spectrum, SurfaceFinish, fresnel_reflectance};

// Reflectance is not a setting. It follows from the refractive indices.
let bk7 = Material::from_catalog("N-BK7").unwrap();
let n = bk7.index(587.56);                       // 1.5168
let bare = fresnel_reflectance(1.0, n, 1.0);     // 0.0421 — the textbook 4%

// A coating can only scale that down, and it does so spectrally.
let coated = SurfaceFinish::broadband_ar().reflectance_at(1.0, n, 1.0, 550.0);
assert!(coated < bare / 10.0);

// A lamp has a temperature, and Planck decides what colour it is.
let tungsten = Spectrum::blackbody(3200.0);
assert!(tungsten.at(450.0) < 0.45 * tungsten.at(650.0));
```

## Licence

Dual-licensed under either [MIT](LICENSE-MIT) or
[Apache-2.0](LICENSE-APACHE) at your option — the Rust ecosystem's convention.
MIT is short and universally understood; Apache-2.0 adds an explicit patent
grant for anyone who needs one.
