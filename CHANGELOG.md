# Changelog

Notable changes, in the format of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This workspace follows [semantic versioning](https://semver.org/). It is `0.x`, so the API is
explicitly not stable and a minor bump may break you. The first consumer exists now, and it
has already found twenty-four places it is awkward, eighteen of which have been changed — see
`crates/dualis-world/FRICTION.md`.

Entries record what was *found* as well as what was added, because several of the more useful
changes here were corrections to a mistaken assumption rather than new features. The commit
messages carry the full account.

`CLAUDE.md`'s gate is a `set -euo pipefail` script now, and carries a table of the five disguises in
which it has reported a pass it had not earned. Four of them are closed by that one line.

## [Unreleased]

### Added

- **Hashin–Shtrikman bounds for stiffness, and a witness that says how tight they are.**
  `Mix::bulk_bounds`, `Mix::shear_hashin_shtrikman`, `Mix::bulk_hashin_shtrikman`.
  `crates/dualis-elastic/tests/a_checkerboard.rs`, and two tests added to
  `crates/dualis-core/tests/a_mixture.rs`.

  The tighter elastic pair, assuming the microstructure is statistically isotropic. For aluminium against
  PLA at half and half it takes the shear range from 5.545-fold to 2.821-fold and the bulk range from
  4.57 to 2.37.

  **Checked against Mori–Tanaka**, which is the same equivalence the conductivity pair has with
  Maxwell–Garnett: HS with the matrix as reference *is* the Mori–Tanaka estimate for spherical inclusions,
  a separately derived result written as a different rational function, and the two agree to `2.2e-16`
  across two decades of inclusion fraction. Several plausible wrong versions of the HS expressions are also
  rational functions with the right limits at zero and one, so a dilute check would not have told them
  apart.

  Three things found, and two of them were mine:

  **The textbook prescription for which phase is the reference does not always mean anything.** "Stiffest
  for the upper bound" needs a well-ordered pair, and aluminium against borosilicate is not one — aluminium
  has the larger bulk modulus, 67.5 GPa against 46.5, and the *smaller* shear modulus, 25.9 against 34.0. A
  first version tested which phase had the larger value of the modulus being bounded, and for that pair at
  a tenth aluminium it returned a lower bound of 48.2312 GPa above an upper bound of 48.1922 — **the pair
  inverted**, by 0.08%, which only a sweep over fractions and pairs would see. Both evaluations are
  computed and then ordered now, which is what "interchange which phase is subscripted one" actually
  prescribes.

  **An affine boundary displacement bounds the effective modulus, not the bound.** The docs asserted that a
  kinematic estimate could not fall below HS+, reasoning that an upper estimate cannot cross an upper
  bound. It does not follow — the apparent stiffness bounds the *effective* stiffness, which is itself
  below HS+ — and the bulk modulus is the counterexample that fired the assertion. What is a theorem is
  `HS− ≤ C_effective ≤ C_apparent ≤ Voigt`, and only that is asserted.

  Read correctly, that makes the measurement a **ceiling on the truth**, and the two moduli answer
  differently: for shear the ceiling lands 0.505% above HS+, so the bound is tight to half a percent at
  worst; for bulk it lands 2.83% *below* HS+, so the bound is at least that loose. **The same pair is
  tight for one modulus of one geometry and loose for the other**, which the algebra does not say and a
  caller picking a number inside it would want to know.

  **A checkerboard one cell per phase is aliased to Voigt exactly**, arriving by a different route than the
  thermal one. There the cause was every face being an interface; here it is kinematic — with an affine
  boundary and a microstructure at the element scale a trilinear element has no freedom to relax into, so
  the affine field is the discrete solution and its energy is the volume average. Asserted, as it is in
  `a_composite.rs`.

  Two convergence parameters, separated deliberately: the period count is the fast knob and moves the
  answer 2.7% before plateauing, and cells-per-block is the slow one that moves *where* it plateaus — two
  cells sit 17.5% above HS+ and four get to 0.5%. A sweep changing both at once would report their sum as
  one rate, which is the mistake `a_layered_wave.rs` records making.

- **`Block::fill`, so the static solver takes per-element material too — and gets the same modulus nine
  orders sharper.** `Block::fill`, `Block::materials`, `Block::material_at`.
  `crates/dualis-elastic/tests/a_layered_block.rs`.

  `Waves` got per-element material first because the closed form that checks a composite's stiffness —
  Backus averaging — is about wave speeds, and there was somewhere to point it. That left one crate with
  per-element material in the wave solver and not the static one, which is an asymmetry with no reason
  behind it.

  Statics turns out to give the **sharper** measurement of the overlapping claim. A traction-driven column
  is an elliptic solve with no time in it — no dispersion relation, no period to fit, no second-order mesh
  error — so a laminate's harmonic constrained modulus comes out to solver tolerance:

  ```text
                                 worst error      debug cost
    a_layered_wave.rs   C33        3.5e-4            67.5 s
    a_layered_block.rs             4.8e-13            0.13 s
  ```

  Nine orders sharper and 520× cheaper. The error there is the conjugate gradient's tolerance accumulating
  over degrees of freedom — `5.9e-15` at eight elements against `4.8e-13` at sixty-four — and not a mesh
  error, which would not care how many iterations it took.

  **The layer thickness does not change the static answer**, and that contrast is a statement neither file
  could make alone. The wave measurement of the *arithmetic* mean is 23% out at eight layers per
  wavelength, because it needs the layers to move together and only a long wave makes them. A static
  uniform stress needs nothing of the kind — equilibrium makes the stress uniform whatever the layers look
  like — so the same column at one, two, four and eight elements per layer gives the same modulus to
  `4.8e-13` with no trend. Backus averaging is a limit for one of those quantities and an identity for the
  other.

  A `fill` on a solved block marks it **unconverged**. The displacement on record solves the previous
  assembly, and left alone it would still read as converged and still return a strain — an answer that is
  present and wrong rather than absent, which is the failure this workspace keeps finding.

  What statics cannot do is the arithmetic end, for the same reason the wave file cannot do `C11`:
  realising it needs the lateral strain zero on average while free locally, and `roller` constrains a face
  rather than a plane through the interior. Statics gets the sharp half; the wave gets both halves less
  sharply. Neither covers it alone.

- **Per-element material in `dualis-elastic`, and stiffness bounds that something can check.**
  `Waves::fill`, `Waves::materials`, `Waves::material_at`, `Elastic: PartialEq`, `Mix::shear_bounds`,
  `Mix::p_wave_modulus_bounds`. `crates/dualis-elastic/tests/a_layered_wave.rs`.

  `Mix` shipped without stiffness bounds a release ago, and the reason was stated at the time: `Waves`
  took one material per block, so there was nothing here a bound on stiffness could be *checked* against,
  and a bound nothing can falsify is a comment rather than an API. This closes that.

  The closed form is **Backus averaging** — the exact long-wavelength moduli of a layered elastic medium,
  Backus 1962 — and `Waves::hold` is what makes each one separately measurable, because freezing two
  displacement components leaves a one-dimensional problem with one modulus in it. Aluminium and PLA, half
  and half, a 20-fold contrast in shear modulus:

  ```text
                                            16/32/64 elements        improvement
    C44  harmonic mu   2.4517 GPa      0.558% / 0.143% / 0.036%        15.7x
    C66  arithmetic mu 13.5945 GPa     1.003% / 0.229% / 0.058%        17.2x
    C33  harmonic M    11.1237 GPa     0.533% / 0.143% / 0.035%        15.4x
    M    arithmetic    53.9839 GPa     0.817% / 0.189% / 0.048%        17.1x   (32/64/128)
  ```

  Second order throughout, where a fourfold refinement predicts 16. **Both ends of the shear pair come
  from one block**, 5.5× apart, decided by which way the shear goes — and the *ratio* of the two speeds
  lands within 0.023% of `√5.545`, tighter than either measurement because the mesh error cancels out of
  it.

  **A composite of two isotropic materials is generally anisotropic**, so `Mix` deliberately reports no
  effective `(E, ν)` and `as_substance` leaves the mechanical block absent. A laminate 5.5× stiffer in
  shear one way than the other has no single isotropic pair, and a function returning one would be
  inventing an isotropy the material does not have.

  Two things the measurements corrected:

  I documented the Voigt end of the P-wave modulus pair as a bound that is **not** attained. It is
  attained — holding the lateral strain at zero *pointwise* is exactly the equal-strain condition Voigt
  assumes, and it converges to `⟨M⟩` at second order. What is not measurable here is `C11`, the *free*
  laminate compressed along its layers, which needs the lateral strain zero on average but free locally;
  `Waves::hold` holds a component on every node or none, so this API cannot pose it. Tried anyway and got
  40.87 GPa, which is neither `⟨M⟩` nor `C11`'s 43.77 — a block four elements thick with free faces
  carries plate modes, so that number answers a third question. No claim is made about it.

  And the arithmetic mean is a **long-wavelength** result in a way the harmonic ones are not, which the
  resolution sweep alone cannot see: refining the span shrinks the element *and* stretches the wavelength
  against the layer at the same time, two second-order errors falling together and reported as one. Held
  the element size fixed and varied only the layer thickness: 0.234%, 1.055%, 4.877%, 23.296% at 64, 32,
  16 and 8 layers per wavelength. **At eight layers per wavelength a caller is 23% out, and that is the
  closed form's limit rather than the solver's.**

  The file costs 68 s in debug, and that is written into it, because the first draft cost 5.4× more for
  the same claims.


- **Composites: two substances made into one, with the bounded properties returned as bounds.**
  `dualis_core::mixture::Mix`, re-exported as `dualis::prelude::Mix`.
  `crates/dualis-core/tests/a_mixture.rs`, `crates/dualis-thermal/tests/a_composite.rs`.

  A motor is copper, steel, magnets and air; a board is FR-4, copper and solder; a buffer is wax in an
  aluminium matrix. Each wants to be **one** `Substance` so a lumped model can hold it, and
  `with_specific_heat` has been saying since 0.4.0 that the bulk `c_p` of such a thing is worth a factor
  of two — while leaving the caller to compute it by hand.

  The properties divide into three kinds and conflating them is the whole failure mode. Density,
  volumetric heat capacity and latent heat are **exact**, from conservation alone. Conductivity is
  **bounded** — Voigt and Reuss for any microstructure, Hashin–Shtrikman if it is isotropic. Emissivity is
  **nothing**: it is a property of the surface and a mixture has no surface, so `as_substance` takes it as
  an argument. A conductivity outside the bounds is **refused**, because no microstructure realises it.

  **One block gives both bounds, to machine precision.** A laminate of alternating aluminium and
  borosilicate measures Reuss `2.213236` across its layers and Voigt `84.057000` along them — same
  material, same fractions, a factor of 38 apart, decided by the direction of the flux. That is why `Mix`
  returns a pair and refuses to pick: a single number would be the upper bound of a 38-fold range
  presented as a measurement. On a half-copper, half-FR-4 board the range is 335-fold.

  Checked against a resolved geometry rather than against itself, and against Maxwell–Garnett for the
  algebra: the Hashin–Shtrikman lower bound with the matrix as host **is** Maxwell–Garnett, to `6.7e-16`
  across three decades of filler fraction.

  Three things found, all of them mine:

  A three-dimensional checkerboard one cell per phase is **not a checkerboard**. Every cell then has all
  six neighbours of the other material, so every face carries the same harmonic mean and the discrete
  operator is that of a *uniform* medium at `harmonic(167, 1.114) = 2.213236` — which at equal fractions
  is Reuss to the last digit. Two drafts reproduced the lower bound and read it as a striking result; the
  microstructure was aliased, not resolved. Asserted now, so it is not rediscovered a third time.

  A coarsely resolved high-contrast composite is under-conductive by enough to **break a bound**: at four
  blocks per axis the measured conductivity goes 3.4926, 5.1396, 6.8213 as the grid is refined, and the
  first of those is below the Hashin–Shtrikman lower bound of 4.3266. Not because a checkerboard violates
  HS, but because two cells cannot represent one. The same pattern in stainless against aluminium does the
  same thing, which is what says it is the discretisation and not the pair.

  A steady-state stopping rule of "the flux stopped changing over N steps" is **not resolution
  independent**: as `dx` shrinks, N steps is a shorter physical time, so the criterion becomes trivially
  true. It reported a 32³ block as converged after 7.5 microseconds of a 5-second process. The rule is an
  absolute residual now — flux in equals flux out, because at steady state nothing is stored — and every
  tolerance in that file traces to a named constant rather than to what made the run pass.

  Elastic and acoustic bounds are deliberately absent: Voigt–Reuss on the bulk and shear moduli is the
  same theorem, but `dualis-elastic` has no per-cell material, so there is nothing here a stiffness bound
  could be checked against. Yield strength is absent for a different reason — a composite's yield is
  governed by the weaker phase and the interface, so a rule of mixtures for it would be wrong rather than
  imprecise.

- **A scene can declare its own substances, so the catalogue stops being the limit.** `Scene.materials`,
  `Palette`, `Substance::CATALOGUE` and `Substance::from_name`.
  `crates/dualis/tests/substances_from_a_file.rs`, `crates/dualis-world/tests/declared_materials.rs`,
  `scenes/21-a-wax-thermal-buffer.json`.

  Nine catalogue entries against hundreds of thousands of materials: enumeration was never going to
  close that gap, and a format that can only *name* a substance can only describe nine kinds of thing.
  A `Substance` is data, so a scene writes one out and uses it exactly as it uses `"ice"`.

  **Checked against Neumann's exact solution, with ice marched beside it through the identical
  harness** — which is the only form in which "as accurate as one we ship" means anything. Gallium and
  n-octadecane, 82× apart in diffusivity, declared as JSON text: worst 0.039% over a 21× range of
  Stefan number, against the catalogue's own ice at 0.035%. Ice sits *inside* the declared range at
  every undercooling, so the error is a function of Stefan number and not of where the substance came
  from.

  Three things found on the way, two of them mistakes of mine:

  `MATERIALS` was a hand-written copy of the catalogue's spelling — eight names beside nine
  constructors — and the missing one was `water`. Unnameable from a scene for **eleven releases**, since
  the format learned to name a material at all. A name absent from a lookup is not a wrong answer, it is
  a substance that never appears, so nothing could have noticed. The lookup lives beside the catalogue
  now and a test checks both directions.

  The first tolerance was `dx/X` — one cell of interface resolution, 8.3%. Measured error: 0.039%,
  **two hundred times smaller**, because the front is read from an integral of a conserved quantity and
  the enthalpy scheme conserves that exactly, so the interface's first-order error never reaches the
  answer. A bound can trace to real physics and still be so loose it checks nothing: that one would have
  passed a conductivity 13% wrong. At 0.06% it catches 0.1%, and that is asserted rather than claimed —
  a test perturbs the conductivity and measures what the front does.

  The teeth test then failed to have teeth. It asked for the same front *depth* from the perturbed
  substance, which recomputes its own target from the wrong conductivity, so the two errors cancelled to
  the last digit and a 1% error reported as 0.0000%. The comment above it had named that exact trap.

  A declaration nothing uses is **refused**, and that is the non-obvious one: `material` on a block is
  optional and defaults to aluminium, so declaring a wax and forgetting to name it gives a block of metal
  that runs, audits and renders while answering about the wrong substance. Same shape as a region
  selecting no cells, one level up.

- **Two-phase conduction, against the two-phase Neumann solution.** `FusionProps::liquid`, `new` and
  `with_liquid`. `crates/dualis-thermal/tests/two_phase_stefan.rs`.

  The one-phase model is exact while the liquid sits at the melting point, because a face with no
  temperature difference carries no heat whatever its conductivity. It stops being right the moment the
  liquid is warmer, and not by a little: **20 K of superheat slows a freezing front 16%**, from 15.85 mm
  to 13.33 mm at 900 s.

  The closed form was established before any code changed, and its check is the **reduction**: setting
  the superheat to zero must give back the one-phase condition. It does, proportionally — the gap is
  `9.446e-12`, `9.447e-9`, `9.446e-6` at `1e-9`, `1e-6` and `1e-3` kelvin, which is `9.4465e-3` per
  kelvin over three decades. Three decades of exactly proportional error is a stronger statement than
  one equality.

  ```text
    superheat   front      two-phase says     one-phase says
        0 K   15.50 mm   15.85 ( 2.20%)     15.85 ( 2.20%)
        5 K   14.70      15.13 ( 2.89%)     15.85 ( 7.27%)
       20 K   12.83      13.33 ( 3.78%)     15.85 (19.07%)
  ```

  Worst 3.78% against a `dx/X` bound of 6.31%, and at 20 K the one-phase model is **19.1% out** — three
  times the bound. `dx/X` and not its square, because the cell holding the interface has a *mixed*
  conductivity and that is a first-order error.

  **Three things I had wrong, all caught by measurement:**

  - **The kelvin-normalised latent heat breaks with two capacities.** `L/c_p` multiplied by a cell's
    *current* capacity charged `4182/2050` = 2.04× to freeze water, and the front came out **27%
    short** with nothing pointing at the latent heat. The enthalpy map is in joules now, with both
    capacities named.
  - **"A liquid at the melting point cannot influence anything" is false for the scheme.** True of the
    continuum; a fixed grid puts a *wholly* liquid cell where the continuum has an interface, and its
    low conductivity throttles the front — 2.1% at twenty cells, closing 2.99× over a fourfold
    refinement. So `Substance::ice` keeps `liquid: None`: switching it on took the *one-phase* answer
    from 0.43% to 6.9%, sixteen times worse for physics that had not changed, and a default that costs
    a caller a factor of sixteen is the wrong default. Two phases are opt-in.
  - **Melting loosens the stability limit and freezing tightens it.** I checked one direction. A step
    sized on an all-liquid block is refused a few hundred steps later once there is ice in it, water's
    limit being 7.57× ice's. The guard caught it, which is how it was found.

  Cost measured before designing: a `resolve` is 1.6× a `step` at 40 cells and 4.2× at 4096, so a
  two-phase sweep runs 2.6× to 5.2× a one-phase one. Simple and whole-grid; an incremental update
  touching only the mush is the optimisation available if a problem needs it, and none does.

  And `FusionProps` gaining a field broke every struct literal — exactly what `Substance` gaining
  `fusion` did one level up, to tests written the week before. `FusionProps::new` and `with_liquid` are
  the same answer applied one level down.

### Changed

- **A packed bed's conductivity is Maxwell–Eucken now, not an arithmetic mean.** `Puck::bed_conductivity`
  and `Puck::conductivity_at`. `crates/dualis-porous/tests/the_beds_conductivity.rs`.

  It was `ε k_l + (1−ε) k_s`, which is the **Voigt bound** — exact only when the two phases lie in
  parallel with the flux, which in a bed of spheres nothing does. Nobody chose it; it is what you get if
  you do not think about it. The structural fact about a saturated bed is that the **liquid is the
  continuous phase** and the grains are dispersed in it, which is precisely what Maxwell–Eucken describes
  and is attained by a coated-sphere assemblage. For coffee at 45% porosity the old value was **11.0%
  high**: 0.38625 against 0.34811 W/m·K, and the honest range narrows from 1.674 to 1.184.

  Cross-checked rather than asserted: the puck's Maxwell–Eucken equals `Mix::hashin_shtrikman`'s bound
  built around the liquid to `2.1e-15` over 76 cases, and `Mix`'s form is itself checked against
  Maxwell–Garnett. Two independently written arrangements of one physics; the puck does not call `Mix`.

  **Why it went unexamined for four releases, measured both ways.** Under flow the bed is isothermal, so
  conduction carries no heat and `λ` is multiplied by zero — swinging it over a factor of eight leaves the
  extraction yield identical to `1e-14`, and scene `18` is bit-for-bit unchanged by this commit. It stops
  being free with a gradient: in a 20 °C basket the yield moves 4.9% per unit `ln λ`, so the range the old
  rule left open was worth about 2% in extraction. Both halves are now tests, so the choice cannot go back
  to being invisible.

  What moved: the cold-basket branch of `espresso_shot`, where the cup cools 3.0102 °C instead of 3.1829
  and the yield reads 19.43% instead of 19.44%. Both stay inside the capacity bound the example asserts,
  and no assertion in the workspace changed — every gradient claim here is an inequality or a structural
  statement rather than a recorded value, which is the only reason a model change was safe to make in one
  commit.

  One thing I chased that was not there: a 287% yield, above the 30% soluble ceiling. The domain saturates
  at exactly 30.00000% when driven to it; a percentage in the printout had been multiplied by a hundred
  twice.

### Fixed

- **`CITATION.cff`'s licence field, which failed the 0.13.0 Zenodo deposition.** `license:
  MIT OR Apache-2.0` is a valid SPDX *expression* and `Cargo.toml` is right to use it; CFF's schema
  takes an identifier or a **list** of them and an expression matches neither. Validated against the
  real schema rather than guessed, and now a two-element list.

  0.13.0 therefore has **no DOI** — it went to crates.io and PyPI and Zenodo rejected the metadata,
  reporting it only as a red *Failed* on its own web page. Nothing in the release, the tag or CI knew.
  0.14.0 will be the first version with one; re-depositing 0.13.0 would mean deleting and recreating
  its GitHub release, and it is not worth doing to a version that is already out.

  `crates/dualis-world/tests/citation_is_valid.rs` is four checks so the next one fails in the gate
  instead of on a web page: the licence is a list, every field a deposition is built from is present,
  the version matches the crate version, and the author block is a list of mappings rather than a bare
  string — that last being the silent half of this class, since a rejection is loud and a record with
  no creator is not.


## [0.13.0] — 2026-08-12

**Elastic waves, any material, and one behaviour change that will bite a deserialiser.**

`dualis-elastic` had a `density` field its own documentation called unused; it sets the wave speeds,
and `Waves` marches `ρü = ∇·σ` on the element `Block` solves statics with. That closed
`ARCHITECTURE.md`'s depth list — the last entry turned out not to be finite strain but **inertia**.

And `Substance` can describe a material this crate never chose. Four `with_*` builders and a `check`
that says whether a number is *possible*, because enumeration does not reach "every material" and data
does. `Substance::bulk` set all four property blocks to `None` and there was no way to fill any of them
except a struct literal — which 0.12.0 broke when it added `fusion`.

### Breaking

- **`Substance` and its four property blocks now refuse unknown keys.** `deny_unknown_fields`, so a
  JSON material with a mistyped `"thermalz"` is an error instead of a substance whose whole thermal
  block is silently absent. If you deserialise `Substance` from a file that carries extra keys, that
  file stops loading — which is the point: it was previously loading as a *different* material.

- **`Waves::new` takes `(name, counts, cell, material)`**, matching `Block::new`. It shipped in no
  release; noted because the two disagreed for a day on `main`.

### Added

- `dualis_elastic::Waves` — `ρü = ∇·σ`, central differences, lumped mass. `hold`, `clamp_ends`,
  `release_mode`, `mode_amplitude`, `mode_frequency`, `kinetic_energy`, `strain_energy`,
  `total_energy`, `displacement_at`. Second order in both speeds — 0.161%, 0.040%, 0.010% over 16, 32
  and 64 elements — and the **ratio** holds to `1e-8` because `E` and `ρ` cancel algebraically while
  the mesh error cancels numerically.
- `dualis_elastic::Axis`, replacing five `usize` axis arguments that disagreed about a fourth one.
- `Elastic::p_wave_speed`, `s_wave_speed`, `speed_ratio`, `from_substance`.
- `Block::as_field` and `Waves::as_field` — the displacement magnitude. Neither half of this crate
  could be drawn before, for as long as the crate has existed.
- `Face::on`, returning an `Axis` where `Face::axis` returns an index.
- `Substance::with_thermal`, `with_mechanical`, `with_acoustic`, `with_fusion`, `check`.

### Fixed

- `bindings/python/Cargo.lock` was stale at 0.10.0 through the whole of 0.12.0. Nothing in that job
  passes `--locked`, so cargo rewrote it silently on every build and the committed copy drifted.

### Documentation

- `RELEASING.md` is new: the release procedure was inside `CLAUDE.md`, which is loaded every session
  for a thing you do once per release. `CLAUDE.md` went from 246 lines to 115.
- `CONTRIBUTING.md` records **six** times the gate reported a pass it had not earned, one of which
  reached `main`, and the measured fact that `set -euo pipefail` does not prevent it in a pasted block.
- No `Co-Authored-By` trailer on commits from `05ed74e` onward.


### Added

- **Any material, not the nine in the catalogue.** `Substance::with_thermal`, `with_mechanical`,
  `with_acoustic`, `with_fusion` and `check`, plus `deny_unknown_fields` on `Substance` and all four
  property blocks. `crates/dualis/tests/any_material.rs`.

  Enumeration does not reach "every material" and data does. `Substance::bulk` set all four property
  blocks to `None` and there was **no way to fill any of them** except a struct literal naming every
  field — and that path broke earlier in this release when `fusion` was added, at the expense of
  exactly the callers the catalogue is least able to help. A chain of `with_*` is immune to a new
  field appearing.

  `check` is what the library can do for a number it did not choose: not whether it is *right*, but
  whether it is **possible**, reporting every problem at once because a material read off the wrong
  column is usually wrong in several places. Plus one check that is not a bound on a single field — a
  substance stating both a sound speed and elastic constants has three independent numbers describing
  one thing, and the speed must sit within **15%** of the rod or bulk speed those give. Measured
  rather than chosen: every catalogue entry is within 6.2% of whichever it means, and that gap exists
  because a tensile test and an ultrasonic measurement are not the same measurement.

  What it catches: a shear speed in the longitudinal slot (39% out), a modulus from the row below
  (27%), an emissivity written as a percentage, a negative conductivity, an incompressible Poisson
  ratio, and four wrong at once reported as four.

  **The test written to demonstrate the capability found the capability was unsafe.** `Substance` had
  no `deny_unknown_fields`, so a mistyped `"thermalz"` was silently dropped and the material ran as
  one whose conductivity is *unknown* rather than one whose file has a typo. The scene format already
  had that rule written down with its reasoning; `Substance` has it now, on all five types.

  Ti-6Al-4V and paraffin wax run through conduction, elastic waves and a melting plateau against the
  same closed forms a catalogue material does — `dx²/6α` to `1e-12`, `√(2(1−ν)/(1−2ν))` at `ν = 0.342`
  which no catalogue entry has, and `ρLV/P` to within one delivery.

- **`Elastic::from_substance`, and "small strain" spans 130× across the catalogue.** The conversion
  from `dualis-core`'s material catalogue to `dualis-elastic`'s material existed **in a test** —
  `two_wave_speeds.rs` built one by hand — which means every consumer wanting to solve an elastic
  problem with a catalogue material wrote the same four lines. It is in the library now, in
  `dualis-elastic`, because the kernel must not learn that elasticity exists.

  It **drops the yield strength**, which is the one thing to know: `Substance` says where a material
  stops coming back and this type has no yield, so a solve past it returns a displacement that is
  arithmetically correct and physically meaningless with nothing in the answer to say which. So the
  strain where the linear model ends is documented and measured, and the spread is the finding:

  ```text
    ice        0.011%     brittle: 1 MPa against 9.1 GPa
    Cu ETP     0.060%
    N-BK7      0.073%
    Al 6061    0.401%
    PLA        1.429%     a polymer has twenty times a metal's elastic room
  ```

  A first draft asserted every entry was under 1% and failed on PLA at 1.43% — correctly, because a
  polymer is not a metal. The bound was wrong, not the data, so the claim is now the **spread**: 130×,
  meaning a strain bound that would be absurdly conservative for PLA is already past yield for ice.

- **Neither half of `dualis-elastic` could be drawn, and its two constructors disagreed.** A
  stabilisation pass, continued: `Block::as_field`, `Face::on`, and `Waves::new`'s arguments in
  `Block::new`'s order.

  `Block` offered no field either, and had not since the crate existed — so an elastic run drew nothing
  and nothing said so, for statics as well as waves. Both now offer `|u|`, trilinearly interpolated,
  checked against the displacement the solver reports rather than against a picture.

  `Waves::new` took `(name, material, counts, cell)` where `Block::new` takes
  `(name, counts, cell, material)`. One crate, two orders, and only a caller writing both in one file
  would ever have noticed. Aligned while `Waves` is unreleased.

  `Face::axis` returns an index and predates `Axis`; `Face::on` returns the type. Two spellings of one
  idea inside one crate drift, so a test pins all six faces to both.

- **`dualis_elastic::Axis`, and `Waves` is drawable.** A stabilisation pass over what landed
  unreleased, which is the moment to change an API rather than after somebody depends on it.

  Five methods on `Waves` took an axis as a `usize` and **disagreed about a fourth one**: four returned
  silently — leaving a body with three free components where the caller believed it had one, and so a
  wave speed that is not the one being asked about — and `mode_frequency` clamped to 2, which is a
  different wrong answer to the same mistake. Three responses to one bad input in one type is not a
  policy, so the input is now unrepresentable.

  And `Waves` offers `as_field` — the displacement **magnitude**, trilinearly interpolated — so the
  analysis layer can draw it. Every other 3D field domain offers one, and a domain that offers none
  gets no picture from a layer whose whole rule is to dispatch on the shape of the data. `|u|` costs
  the sign, which is stated: a full standing wave shows both antinodes bright. `Block` still offers
  nothing here, which is a real gap and older than `Waves`; it is on a published type, so that is a
  separate decision.

  The per-step cost was **measured before deciding** anything about it: 2.47 ms/step at 16³ elements
  and 19.69 ms at 32³, linear in nodes, so the two `Vec` allocations a step makes are under 1% and are
  left alone. A wave crossing a 32 mm body is about 90 steps; a body left ringing for a millisecond is
  8,700 and three minutes.

- **Elastic waves: `dualis_elastic::Waves`, and `ARCHITECTURE.md`'s depth list is closed.**
  `ρü = ∇·σ` on the same trilinear element `Block` solves statics with, marched with central
  differences and a lumped mass. `Waves::hold`, `clamp_ends`, `release_mode`, `mode_amplitude`,
  `mode_frequency`, `kinetic_energy`, `strain_energy`, `total_energy`, `displacement_at`.
  `crates/dualis-elastic/tests/two_speeds_marched.rs`.

  A sibling type rather than a mode on `Block`, because `Block` is `Kind::QuasiStatic` and is right
  about itself — a body with velocity has a different lifecycle, a different stability limit and a
  different thing to conserve. What they share is the element, which is the part worth sharing: the
  static tests check that operator against four exact moduli and these check it against two exact
  speeds.

  ```text
    elements   c_p error   c_s error
          16     0.161%      0.161%
          32     0.040%      0.040%
          64     0.010%      0.010%
  ```

  Second order in both, a factor of four per halving to two figures. And the **ratio holds to about
  `1e-8`** at `ν` = 0.2, 0.33 and 0.45 — four orders tighter than either speed alone, because `E` and
  `ρ` cancel out of it algebraically *and* the two modes share a shape so the mesh error cancels out
  of it numerically. What is left being checked is the operator.

  The **leapfrog's own dispersion is removed rather than tolerated**: central differences on one
  eigenmode turn at `Ω = 2·arcsin(ω dt/2)` per step, so the measured `Ω` is inverted through
  `ω = 2·sin(Ω/2)/dt` before anything is compared. Comparing raw periods would add the time error to
  the space error and call the sum an accuracy.

  Energy swings by `2 sin(ωΔt/2)` — measured 0.0348 against 0.0348 predicted — and does **not** drift:
  the mean over the last eighth of a forty-period run matches the first to `2.5e-7`.

  The stability limit is `2/√λ_max(M⁻¹K)` by Gershgorin on the assembled rows, not a borrowed Courant
  number. Measured it is **1.229×** `dx/(c_p√3)`, so the Courant form would have been safe and 23%
  wasteful — the opposite direction from the one assumed before measuring it, which is why it is
  computed.

- **A solid has two wave speeds, and the catalogue's one field means both of them.**
  `Elastic::p_wave_speed`, `s_wave_speed` and `speed_ratio` in `dualis-elastic`, and
  `crates/dualis/tests/two_wave_speeds.rs` — the fourth cross-domain test of that shape, after
  `fields_and_rays`, `loss_and_lumps` and `a_slit`.

  `Elastic::density` was documented as "unused by the static solve and carried so a mass can be
  stated once". It sets the wave speeds, which is the only place a mass enters a problem with no
  inertia in it.

  The sharp statement is `c_p/c_s = √(2(1−ν)/(1−2ν))`: both `E` and `ρ` cancel, so it holds at
  **2.2e-16** for all six catalogue entries and a scheme with the wrong stiffness or the wrong mass
  still has to land on it.

  The finding is what `AcousticProps::sound_speed` turns out to mean. A fluid has one speed; a solid
  has two longitudinal ones — bulk `√((λ+2μ)/ρ)` and rod `√(E/ρ)` — and the six entries do not all
  mean the same one:

  ```text
                      rod   stated    bulk    vs bulk   vs rod
    ice              3150     3840    3834      +0.1%   +21.9%     bulk
    Al 6061          5052     6320    6149      +2.8%   +25.1%     bulk
    304 stainless    4912     5790    5623      +3.0%   +17.9%     bulk
    Cu ETP           3614     4760    4483      +6.2%   +31.7%     bulk
    N-BK7            5716     5680    6048      −6.1%    −0.6%     rod
    elec. steel      5113     5100    5853     −12.9%    −0.3%     rod
  ```

  Four bulk, two rod, every one within 6.2% of whichever it is. And read as a bulk wave the stated
  speed implies a modulus: 9.1 against 9.1 GPa for ice, but 72.8 against 68.9 for aluminium and 131.9
  against 117.0 for copper — because **a tensile test and an ultrasonic measurement are not the same
  measurement**, and the dynamic modulus of an annealed metal comes out higher. Ice agrees exactly
  because its elastic constants were back-calculated from velocity, so the one entry that agrees is
  the one where agreement was never independent.

  Two earlier drafts asserted `rod ≤ stated ≤ bulk` with 0.1% and then 2% of slack. Both were false —
  three entries sit above the bulk speed and two below the rod speed — and the second failed on
  aluminium, which is not an edge case. The assertion is now a classification: each entry is within 7%
  of one of the two, and that the split is *mixed* is itself asserted, so a future edit cannot make
  the point disappear by making them uniform.


## [0.12.0] — 2026-08-12

**Two of `ARCHITECTURE.md`'s three depth entries, and both the same shape.** Single-material and
single-phase are answered; small-strain is the one left. Neither cost a new crate, neither touched the
kernel or either layer above it, and each turned on a **closed form** rather than on a feature — the
harmonic mean making a layered wall's resistance exact, and Neumann's solution putting a freezing
front in a place rather than at a rate.

Also diffraction, which measures where scalar optics stops being true rather than assuming it is.

Twenty-two new public items across four crates. `Solid3D::max_stable_dt` reports a different number
for a thin or filled block than 0.11.0 did — read the entry below before pinning a step by hand.

### Added

- **Latent heat, against Neumann's exact solution of Stefan's problem.** `LatentHeat` in
  `dualis-units`, `FusionProps` and `Substance::fusion`, `Substance::latent_energy` and
  `Substance::ice` in `dualis-core`, `Solid3D::set_melted_fraction`, `melted_fraction_at` and
  `melted_volume` in `dualis-thermal`, `"ice"` in the scene format's materials.
  `crates/dualis-thermal/tests/a_freezing_front.rs` and `20-melting-a-block-of-ice`.

  The second of `ARCHITECTURE.md`'s three depth entries — every domain was single-phase — and it left
  one: small-strain. No new crate, and nothing outside `dualis-thermal` but the material data.

  **Why this problem:** it has an exact solution, which almost nothing with a moving boundary does. A
  semi-infinite liquid at its melting point, surface dropped to `T_s`, freezes to `X(t) = 2λ√(αt)`
  with `λ e^{λ²} erf(λ) = St/√π`. A *position*, not a rate and not a limit. Measured for ice under 20 K
  of undercooling:

  ```text
    dx = 1.0 mm    t = 100 s    5.311 mm against 5.288 mm     0.431%
                   t = 400 s   10.570 mm against 10.567 mm    0.032%
                   t = 900 s   15.843 mm against 15.848 mm    0.033%
    dx = 0.5 mm    t = 100 s    5.285 mm against 5.283 mm     0.032%
                   t = 900 s   15.847 mm against 15.848 mm    0.012%
  ```

  And nine times the time gives **2.9984×** the depth against `√9 = 3` — the `√t` law with no closed
  form in it at all, which a scheme with the right coefficient and the wrong power would fail.

  **The scheme is enthalpy**, bookkept as a temperature and a melted fraction: a cell's state is one
  monotone number — `T − T_m` below, `φ·ℓ` inside, `ℓ + T − T_m` above — so the sweep adds energy and
  inverts it. Nothing tracks the front, nothing iterates, and energy is conserved because energy *is*
  the state. Not **apparent heat capacity**, which smears `L` over a temperature interval and lets a
  step big enough to cross the interval skip the latent heat silently, running the front fast. Here an
  overshoot's remainder lands on the far side because the inverse map says where that much energy
  goes — verified at machine precision with ten times the latent heat in one delivery.

  A cell fed constant power holds at its melting point for exactly `ρLV/P`: measured 305.9 s against
  305.9 for a cubic millimetre of ice at a milliwatt. For ice `L/c_p` is **163 K**, so a scheme that
  dropped the latent heat would not be slightly wrong.

  Four things the tests found, all of them mine:

  - **The profile's convergence order is not clean** between adjacent resolutions — 0.78 then 1.84,
    and averaging five instants did not settle it. A fixed grid makes the front advance in a
    staircase, so the field at any moment depends on where the front sits between two cell centres and
    that phase is not a smooth function of `dx`. Over the full fourfold refinement it is 6.2×, between
    first and second order. Two earlier drafts asserted second order on no evidence and then first
    order on evidence that did not earn it.
  - **"Every cell behind the front" is the wrong metric.** It measures a different set of points at
    each resolution — the nearest one to the front moves inward as `dx` shrinks, into exactly the
    region a cell-centred profile represents worst. That alone read as order 0.66. Fixed depths.
  - **The front position is more accurate than the field it comes from** (0.03% against 0.8%), because
    `melted_volume` is an integral of a conserved quantity where a temperature is a point sample.
  - **A checkpoint of temperatures alone is not a checkpoint.** 0 °C is ice, water or any mixture, so
    `checkpoint`/`restore` carry the phase — a live path, since `Schedule::Iterative` restores every
    iteration and the conservation audit restores on a violation. Without it a half-frozen column
    comes back as a fully liquid one at the same temperatures, and the ledger balances, because the
    ledger reads the state that was corrupted.

- **A block can be made of more than one material.** `Solid3D::fill`, `substance_at`, `substances`,
  `face_conductance`, `heat_capacity` and `stability_ratio` in `dualis-thermal`, and `regions` plus
  `material` on the scene format's `block`. `crates/dualis-thermal/tests/a_layered_wall.rs` and
  `19-a-coating-stops-the-heat`.

  The first of `ARCHITECTURE.md`'s three *depth* entries — single-phase, single-material,
  small-strain — and it cost no new crate and nothing outside `dualis-thermal`.

  The number it turns on is the conductivity **on a face**, which is the harmonic mean
  `2k_Lk_R/(k_L+k_R)` and not the arithmetic one. For aluminium against borosilicate those are
  2.21 and 84.1 W/m/K, a factor of **38**, and the arithmetic mean short-circuits the interface.
  The harmonic mean earns an *equality*: with the material interface on a cell face the discrete
  chain of face resistances is exactly `Σ Lᵢ/(kᵢA)` at every resolution — measured at exactly zero
  for twelve cells and `2.2e-16` for twenty-four, forty-eight and ninety-six. A layered wall's
  resistance has no discretisation error at all.

  The arithmetic mean's is **first order**, which is the dangerous kind — 8.852%, 4.233%, 2.072%,
  1.025% over those four. It vanishes on refinement, so a single-resolution check would have read it
  as a discretisation error and given it a tolerance. Reaching 0.1% would take about **984** cells;
  harmonic is there at twelve.

  And that rate is itself an equality, which asserting `≈2` had obscured. The wrong resistance is one
  face's worth against `(n−1)/2` cells of the two layers, so the relative error is exactly
  `2(1/H − 1/A)/(1/k_a + 1/k_g) / (n−1)`. Multiplying it by `n−1` gives **0.973670** at all four
  resolutions and matches that closed form to `1e-12` — so the coefficient is pinned, not just the
  order. The first draft asserted the consecutive ratios were 2 within 0.1 and passed at 0.091,
  because the exact ratios are `23/11`, `47/23` and `95/47`: a tolerance absorbing a systematic
  offset it could have predicted.

- **Diffraction, and where scalar theory stops being true.** `single_slit_intensity` and `slit_zero`
  in `dualis-optics`, and `Cavity::obstruct` in `dualis-em` for a perfect conductor with a hole in
  it. `crates/dualis/tests/a_slit.rs` holds one against the other.

  `sinc²(π a sinθ/λ)` rests on **Kirchhoff's** boundary condition — the field in the opening is the
  incident field — and Maxwell's equations do not say that. The metal carries currents and the
  opening's field is perturbed over a strip of order `λ` at each edge, so the fraction of the
  aperture that is wrong goes as `λ/a`:

  ```text
    a = 12λ    0.0057     largest absolute difference in normalised intensity
    a =  6λ    0.0125
    a =  3λ    0.0311
    a =  1λ    0.2772     wrong by more than a quarter of the pattern it predicts
  ```

  A factor of **48**, monotone. And the first dark fringe lands at `sinθ = λ/a` to 0.66% at twelve
  wavelengths — the sharpest number in the pattern, and the one with no special constant in it,
  where a circular aperture's `1.22 λ/D` has the first zero of `J₁` in it.

  Below a wavelength a slit has **no** zeros at all, which is what `slit_zero` returning `None` is
  saying: a 0.75λ opening puts 79% of its axial intensity at 30° off axis where a 12λ one puts
  0.05%. That is the difference between an aperture and an antenna.

- **The far field out of a small box.** Not by making the box `a²/λ` deep — hundreds of wavelengths
  for a wide slit — but by recording the aperture plane one wavelength behind the screen and Fourier
  transforming it, which **is** the Fraunhofer limit exactly.

  Magnetic walls make the incident plane wave possible next to a screen, and their mirror images
  would make this a grating — except that the transform integrates over one box width and the
  neighbouring slits lie outside it.

- **Three things about the cost of asking.** A one-cell screen leaks, because the electric edges
  either side of it are shared with the vacuum beyond; two cells hold. A `pattern(width, samples)`
  that marched per call marched the same box eight times where five would do, at 120 s in a debug
  build — the mode four of CI's jobs run; separating the march from the sampling took it to 64 s.
  And the numbers in this entry are measured rather than remembered: the first draft of
  `single_slit_intensity`'s own tests claimed a sidelobe of `4/(2.25π²)` where it is `1/(2.25π²)`,
  and a grazing intensity of 0.8 where it is 0.2545.

### Changed

- **`Solid3D::max_stable_dt` sums the face conductances instead of dividing by a diffusivity**, so
  it is `minᵢ Cᵢ/(dx·Σ_f k_f)`. Two answers moved, and a caller stepping by hand will notice.

  It is now **shape-aware**: `dx²/(2α)`, `dx²/(4α)` or `dx²/(6α)` according to how many axes have
  more than one cell. This domain charged every shape the three-dimensional rate, so a bar-shaped
  block — the shape its own closed-form tests use to check the axes against a 1D answer — paid three
  times the steps it needed. At each of the three limits the sharpest representable mode amplifies
  by exactly `−1.000000000000000`, which is what makes the limit an equality rather than a caution.

  And a filled block is limited by its worst **cell**, which is usually far *looser* than
  `dx²/(6·α_max)`: one aluminium cell inside borosilicate is stable at exactly `k/k_face` = 167/2.21
  = **75.45×** aluminium's own limit, because heat cannot reach a cell faster than its worst face
  delivers it. It can go the other way by at most 2×, and reaching that needs neighbours that conduct
  better *and* store more — impossible for real solids, whose volumetric heat capacity spans one
  order of magnitude where conductivity spans four. Measured at 1.9945× with a doctored specific
  heat, which is the only way to get there.

- **`Solid3D`'s ledger and its placeless heat are weighted per cell.** `Σ CᵢTᵢ` is what a
  two-capacity block conserves, not the mean temperature, and heat off the plain channel spreads to
  a uniform *rise* — in proportion to capacity, not in equal joules, which would have warmed the
  poorer-storing material more and so claimed a location the bus never carried.

- **`ScalarField::rate` on a block is the conductance form**, `(1/Cᵢ)Σ_f G_f(T_f − Tᵢ)`, read off the
  same faces the sweep uses rather than reconstructed as `α·∇²T` from a diffusivity a filled block
  has no single value of. `laplacian` stays `∇²T`, which is what the trait asks for.

- **The scene format's material names live in one table**, `dualis_world::MATERIALS`, with
  `stainless_304`, `borosilicate` and `ice` added. It was a `match` inside the network builder; the
  moment a second domain wanted a material that would have been two lists that agree until they do not.

- **`Solid3D`'s ledger is in enthalpy**, so a melting front is on the books. A cell holding at its
  melting point while it absorbs 306 mJ per cubic millimetre has taken that heat in and its
  temperature says nothing about it; an audit reading temperature alone would call it a leak of exactly
  that. Every door heat comes in through — the sweep, `deposit`, the plain channel — goes through one
  place, so a phase change cannot be forgotten at one of them.

- **A `Solid3D` that can melt reports a `melted` reading** in mm³, and one that cannot does not. A
  column of zeros in every report tells a reader nothing; the condition is a property of the block's
  materials fixed at construction, so it is not a mode that can surprise anybody mid-run. Without it a
  phase change was the one thing the domain does that a report could not see.

## [0.11.0] — 2026-08-11

**Four new domains, and the list `ARCHITECTURE.md` opened with is closed.** Electromagnetism,
elasticity, fluids and flow through a packed bed are crates on the kernel now, and none of them
needed the kernel or either layer above it to change — rule 4 held ten times. Ten domains, sixteen
crates, 546 tests.

The other half of the release is two tests that make the *relationships* checkable:
`fields_and_rays.rs` puts a Yee grid against Fresnel's algebra, and `loss_and_lumps.rs` puts a
field's decay in a conductor against a `Winding` that answers the same number at every frequency.
Neither closes a gap by making one domain do another's job.

### Added

- **A waveguide, and its whole dispersion relation from one march.** A box with conducting side
  walls and open ends *is* a rectangular guide, so this needed nothing but a launcher that takes a
  transverse profile and a driven source.

  ```text
    f_c = c/2a                    exact, and a property of the cross-section alone
    β = √(k² − k_c²)              0.06% to 2.2% from 1.4 f_c to 3.1 f_c
    α = √(k_c² − k²)              0.4% to 5.5% below cutoff
    v_p v_g = c²                  1.0015, with v_p = 1.14c and v_g = 0.88c
  ```

  A pulse carries a band and the mode's `sin(πx/a)` profile is orthogonal to every other, so one
  march gives `β` at every frequency in it. `λ_g/λ` runs from **1.43** near cutoff to 1.05 at the
  top — a check at one high frequency would pass for a solver with no cutoff at all.

- **Five things about measuring a guided wave, four of them mistakes.**

  **A driven source is not a convenience below cutoff, it is the only way to ask.** An evanescent
  field is a near field: it does not travel, and with nothing driving it it decays to zero, so there
  is no steady state and no spatial profile to fit. The first attempt used an initial-value pulse and
  returned decay constants of 0.285, −0.891 and −0.102 against a closed form of 2.7 — a *negative*
  decay being the measurement saying it was fitting noise. `Cavity::impress` adds to `Ey` on a plane;
  it must be called **after** the step, because the electric update rewrites every interior face and
  erases a source added before it.

  **Mur's boundary is tuned to `c` and a guided mode is not.** Its phase velocity is
  `c/√(1−(f_c/f)²)` — 2c at `1.25 f_c` — so the absorbing faces are worst matched exactly where the
  measurement is most interesting. `β` came out 15% low there against 2.3% at `1.4 f_c`. Truncating
  the march before each reflection returns recovers the band above `1.4 f_c`; below it the near end's
  reflection and the direct pulse overlap and no window separates them. That is a limitation of the
  boundary rather than of the scheme, and which one it is is now written down.

  **A band that does not cover the frequency asked about returns a transform of rounding.** The first
  version read `β` correctly at its carrier, to 0.84%, and 87% wrong two octaves away.

  **A transform window that is not a whole number of cycles leaks between bins**, and with the
  amplitude falling by `e⁻²` across the stations that leakage was the difference between 13% off and
  0.4%.

  **A longer lever on a noisier point is not a better fit.** Widening the stations to hold a fixed
  number of decay lengths made the evanescent measurement worse — 12.4% where the tight span gives
  5.5% — because the far stations sit among the boundary's leakage.

- **What a lumped resistance cannot say, measured.** `Winding` states a resistance; `Conductor`
  solves for one from a shape but is quasi-static, so neither has a frequency in it. A field does.

  `crates/dualis/tests/loss_and_lumps.rs` measures the amplitude decay inside a conducting
  half-space against the closed form that holds for **any** conductivity, not just a good one:

  ```text
    p = σ/(ωε)                          the loss tangent
    α = ω√(εμ/2)·√(√(1+p²) − 1)        exactly
  ```

  A pulse carries a band, so one march gives the decay at every frequency in it: a discrete
  transform at each depth separates them, and three frequencies come out 0.04%, 0.89% and 1.15%
  from their own `α`. Measuring the *dependence* rather than one number is the point — a single
  driven sinusoid would hide whether the frequency enters correctly.

  Then the handbook `δ = √(2/(ωμσ))` as what it actually is: the `p ≫ 1` limit, off by **55.4%** at
  `p = 1`, 5.1% at 10, 0.50% at 100 and 0.050% at 1000. The 55% is the part worth knowing — a lossy
  dielectric is not a conductor, and using the skin depth on one is not a small error.

  Beside it, copper: 9.22 mm at 50 Hz, exactly a tenfold thinner per hundredfold in frequency, and
  a `Winding` reporting 0.1724 Ω at every one of them.

- **`dualis-em` carries structures now, and a Yee grid and Fresnel's algebra agree on a
  reflectance.**

  Per-cell media, so a slab or a coating goes into the box, and a *magnetic* boundary so a plane
  wave can exist in one at all: a wave along `z` polarised along `y` has `Ey` tangential to the `x`
  faces, and a conductor there turns it into a waveguide mode with a cutoff and the wrong phase
  velocity.

  `crates/dualis/tests/fields_and_rays.rs` is the check the two crates existed to make possible.
  `dualis-optics` gets a reflectance from algebra on two indices; `dualis-em` marches Maxwell and
  knows nothing about either. They share no code and do not depend on each other:

  ```text
    20 cells per wavelength    −10.543%
    40                          −2.556%     ratio 4.13
    80                          −0.634%     ratio 4.03

    n = 1.5   3.981% against 4.000%
    n = 2     11.041% against 11.111%
    n = 3.5   30.519% against 30.864%
  ```

  The claim is the **rate**, not a tolerance: a single resolution with a loose bound passes for a
  scheme converging at first order, and first order is what an interface placed half a cell wrong
  would give.

  Beside it, two things the algebra alone cannot say. A **quarter-wave coating** of the geometric
  mean index takes a 14.4% reflection to 0.23% — a factor of **61**, and it is interference between
  two surfaces rather than a property of either. And a slab **delays** a pulse by `(n−1)d/c` to
  2.2%, measured on the energy centroid.

- **Three things the cross-domain test found about how to measure.**

  **A monitor plane inside the pulse measures nothing.** The first geometry put a pulse of
  half-width 1.6λ at 4λ and the monitor at 5λ, so the "incident" window opened with the pulse
  already on it: the reflectance came out 34.7% where Fresnel says 14.8%.

  **The peak of an oscillating pulse is a peak of its carrier.** It moves in steps of half a period
  — 0.83 fs here, against delays of a few femtoseconds — so the slab delay read 9% high for that
  reason and no other. The energy centroid is immune to it.

  **A one-dimensional problem does not want a two-dimensional cross-section.** The testbed was
  4 × 4 cells across, every one holding the same field by construction, and the file took **486 s**
  in a debug build — the mode CI's OS matrix runs. Two by one is the smallest grid with an interior,
  gives the same numbers to the last digit, and takes 36 s.

- **An open boundary for `dualis-em`, and the source that showed why the first measurement of it
  was meaningless.**

  Mur's first-order condition: exact for a wave arriving along the normal, worse away from it as
  `(1−cos θ)/(1+cos θ)`. Not a perfectly matched layer, which would be right to four digits and is
  more machinery than this crate has earned — so the figure is **measured** rather than claimed. A
  line source in an open 24³ box is down to **0.149%** of its energy after two and a half crossings,
  against a conducting box's 101%.

  It unblocks the field formulation of optics: diffraction cannot be computed in a box whose walls
  send everything back.

  **The first source left 34% behind, and the boundary was not at fault.** It was a Gaussian blob
  of `Ey`, which has `∂Ey/∂y ≠ 0` and therefore `∇·D ≠ 0` — a *charge distribution*, whose near
  field does not propagate. No absorbing boundary can remove what was never going anywhere, and a
  measurement using such a source measures the electrostatics instead of the boundary. The source
  is a blob of `Hz` uniform along `z` now, which is divergence-free on the discrete grid exactly.

  Opening a face does not cost the divergence identity, and the test says so: `∂(∇·B)/∂t = −∇·(∇×E)`
  holds for *any* `E`, so changing `E` at a boundary cannot touch it. A measurement that said
  otherwise would mean the boundary was writing to `H`.

  `Cavity::mur_coefficient` is `(cΔt − Δ)/(cΔt + Δ)` — a function of the Courant number alone, which
  a test checks by getting the same number out of vacuum at 2 mm and glass at 6 mm.

- **A tenth domain: `dualis-fluid`, and the list `ARCHITECTURE.md` opened with is closed.**
  Incompressible Navier–Stokes by projection on a staggered grid — velocities on cell faces,
  pressure at centres, which is the arrangement Yee uses for electromagnetism and for the same
  reason.

  It was the hardest of the three to make checkable, and the crate is built around admitting why:
  fluids has few exact solutions, its schemes trade stability against numerical diffusion, and
  **"it looks like a fluid" is the easiest wrong answer in computational physics to accept.** So
  each test is chosen to be blind to a different mistake, and the docs say which.

  ```text
    Poiseuille      exact to 6e-15 — against the *discrete* parabola
    Couette         exact to 1.2e-14 — and blind to advection entirely
    Taylor–Green    0.272% → 0.065%, ratio 4.17 — the only one that sees the nonlinear term
    uniform flow    unchanged to 0.000e0 after 500 steps
    momentum        drifts 1.1e-15 in a periodic box, because the advection is in flux form
  ```

  Two limits on the step — viscous `dx²/6ν` and Courant `dx/|u|` — and one on the **mesh**: the
  cell Reynolds number `|u|dx/ν ≤ 2`, which no amount of shortening the step fixes. Refused rather
  than run, because past it central differences produce a sawtooth that reads as turbulence.

- **The discrete Poiseuille profile is not the continuum one, and the difference is a closed form.**

  The first test compared to `gh²/12ν` and measured `19/18` of it at six cells. That is
  `1 + 2/36` to the digit — and the scheme was right.

  Any quadratic satisfies the interior equation exactly, so the discretisation contributes nothing
  there. The **wall** does: a no-slip condition imposed by reflecting the first cell makes the
  *linear interpolation* between them vanish at the wall, and a parabola is not its own linear
  interpolation. The discrete answer is `u_j = (g/2ν)[(h²+Δ²)/4 − (y_j−h/2)²]`, whose mean is
  `(gh²/12ν)(1 + 2/n²)`.

  The test now makes two statements instead of one: the discrete form holds to machine precision at
  any mesh, and the gap to the continuum shrinks by exactly `(16/6)² = 7.111×` from six cells to
  sixteen. Comparing to the continuum alone and calling the difference a tolerance would have hidden
  both.

- **A periodic seam that was averaged instead of copied cost 4.7% of a decay rate.**

  With `Walls::None` the `y = 0` and `y = h` velocity faces **are** the same face. `apply_walls`
  set both to their mean — but the update writes the low one and leaves the high one stale, so the
  mean moved it only half as far as the physics did. A half-step lag dressed as a boundary
  condition.

  It showed up as a Taylor–Green decay 4.68% below the discrete Laplacian's own eigenvalue.
  Diagnosed by running the same vortex at an amplitude eight orders smaller: the rate was identical
  to six digits, which cleared the advection and left diffusion holding the bag. After the fix the
  measurement sits **1.02% and 0.26% above** that eigenvalue at the two meshes — which is forward
  Euler's `σΔt/2` to the digit — and 0.27% and 0.07% below the continuum, because the spatial error
  is larger and has the other sign. Two errors, each with a formula, which is what makes the
  comparison mean anything.

- **A ninth domain: `dualis-em`.** Maxwell's equations in three dimensions, on the grid Yee built
  for them, and the reason that grid is the grid.

  `E` on cell edges and `H` on cell faces, each half a cell from the others, so every curl is a
  difference of quantities that already sit where the curl belongs. The consequence is the point:
  the **discrete** divergence of the **discrete** curl is identically zero — every term appears
  twice with opposite signs — so `∇·B = 0` is a property of the update rather than something the
  scheme approximately respects.

  Demonstrated directly, from both sides. A released mode has `8.6e-17`; two thousand steps later
  it is `1.1e-12`, which is `f64` rounding accumulating linearly and twelve orders below the field.
  And injecting `1.000e-3` A/m of divergence by hand, then running 500 steps of live
  electromagnetism, leaves `1.000e-3` A/m — **1.000000000×**. The update cannot touch it in either
  direction.

  Beside that: cavity resonances against `(c/2)√((m/a)²+(n/b)²+(p/d)²)`, converging at 0.536% →
  0.134% → 0.033% for refinement ratios of **4.002 and 4.001**; the Courant limit `dx/(c√3)`,
  refused above and stable at; and equipartition to 0.998.

- **Four things the electromagnetic tests found, three of them in the tests.**

  **`H` has to be the discrete curl of `E`, not a sampling of the continuum one.** The identity is
  a property of the operators, and a continuum field that happens to be divergence-free does not
  inherit it: setting `H` analytically measured `3.3e-4` where the identity says zero. Applying the
  scheme's own half-step update to a zero field gives both the right initialisation *and* exact
  divergence-freedom, and deletes the hand-derived curl.

  **Joule heating is `∫σE²dt`, not the change in field energy.** The latter is the loss plus
  whatever flowed in from the magnetic field on the same step, and it reported dissipating
  thirty-six times the energy the cavity ever held.

  **Eight half-cycles cannot see second order.** A peak located to the nearest step carries an
  error of `dt`, and eight of them divides it by eight — leaving 0.6% of measurement noise on a
  0.55% signal. The first version reported the *same* frequency at 8, 16 and 32 cells, to seven
  digits, and it took a probe to establish that the scheme was fine and the ruler was not. Two
  hundred half-cycles with the peak interpolated through a parabola resolves it.

  **A 12.8% energy swing is not a defect.** `½εE² + ½μH²` is not what leapfrog conserves — `E` and
  `H` are half a step apart — and the naive sum oscillates by exactly `2 sin(ωΔt/2)`. The first
  version asserted "under 5%", measured 12.779%, and called it a failure where the formula says
  12.8167%. The test checks the closed form now, and that halving the step halves the swing:
  ratio 2.000. An arbitrary bound cannot tell a correct oscillation from an incorrect one.

  A fifth, smaller: a *normalised* divergence divides by `max |H|`, which oscillates through the
  cycle, so comparing it at two instants compares two denominators. The injection test read 1.099
  and then 0.5387 for a quantity that had not moved.

- **An eighth domain: `dualis-elastic`.** What a shape does under load — `∇·σ = 0` with
  `σ = λ tr(ε)I + 2μ ε` — solved rather than stated, so a stiffness is a property of a geometry the
  way `Conductor` made a resistance one.

  The third elliptic domain here and the first whose unknown is a **vector** at every node. That is
  the only thing that is new, and it is why the operator is assembled from the strain energy rather
  than differenced from the Navier–Cauchy equation: a Hessian is symmetric by construction, which
  conjugate gradients needs, and `∇(∇·u)` differenced on a collocated grid admits a checkerboard
  displacement that costs no energy.

  **Four moduli come out exactly**, at any mesh, because trilinear elements reproduce a linear
  displacement field exactly:

  ```text
    uniaxial stress    E = 68.90 GPa      and nu = 0.330000000 with it
    confined           M = 102.09 GPa     = lambda + 2mu, 1.48x E
    hydrostatic        K = 67.55 GPa
    simple shear       G = 25.90 GPa      where 2G would be the classic factor-of-two error
  ```

  Four combinations of two constants. A solver with `λ` and `μ` transposed reproduces none; one
  whose shear rows carried `2μ` passes the first three and fails the fourth. Checking one modulus
  is checking that a stiffness exists.

  Beside them: the **patch test** on a general linear field with a rotation mixed in — the interior
  lands on it to `1.9e-20 m` — **Clapeyron's `2U = Σf·u`** to `4e-16`, and a null space that is
  exactly the six rigid motions and nothing else, each costing `1e-32 J` against a reference
  strain's `9e-6`.

  **Bending only converges**, from the stiff side, and that is stated rather than tuned away: a
  fully integrated element develops shear where it should flex, so a cantilever runs 0.866, 0.933,
  0.960, 0.981 of `PL³/3EI` as the mesh refines. Reduced integration would cure it and buy hourglass
  modes, which is trading a stiffness error for a *singularity*.

- **Three things the elasticity tests found about their own design.**

  **A shear rig is not simple shear.** Clamp the base, drag the top, read `τ/γ` — and a cube gives
  **0.40 of `G`**, because the sides are free and the block bends. `prescribe_boundary` puts the
  analytic field on the whole surface instead, which is the patch test and removes the rig from the
  answer. Displacement control needed the solver to start its iterate at the prescribed values
  rather than at zero, and the relative residual to be measured against the initial one when there
  are no loads at all.

  **A reaction has to be read by component.** A node on an edge belongs to two faces, so
  `reaction(YLow)` sums the `y`-low share of the `x`-low rollers' 30 N and reports 2.5 N of a force
  it is not carrying. `normal_reaction` is the unambiguous one; the magnitude is not, and the test
  now asserts both — that the normal is zero *and* that the magnitude is not, so the distinction
  cannot quietly stop mattering.

  **Two domains wanted the name `Body`.** `dualis-mechanics` had it for a rigid one. The prelude is
  where that shows, and only because both are exported there. The elastic one is `Block`.

- **The native viewer opens the espresso run, and two things `runtime/viewer`'s README already
  claimed turned out not to be true of its renderer.**

  `portafilter_flow flow.json` writes the wire format, and `runtime/viewer` reads it — that viewer
  does not link `dualis` at all, so a run that opens there is the format demonstrating it carries
  enough to draw a run. It did, unchanged.

  What it drew was wrong in two ways, both found by pointing it at a run shaped unlike the
  fixtures: tall, thin, and spanning two orders of magnitude in value.

  **The colour scale was per frame.** `Run::scale_of` computes the range over the whole run, is
  documented as the thing that stops a decay looking like a steady state, is tested — and
  **nothing called it**. `segments` measured the values it had in hand. Water leaving a screen
  clean and arriving at a spout at 83 kg/m³ therefore rendered mid-ramp the whole way down. The
  new test checks the *consumption* rather than the accessor, which is the difference that let
  this survive.

  **The camera was set up for a cube.** `Framing` normalises by the longest side, so a tall thin
  run fills one axis and a fraction of the others: the portafilter came out at 15% of the frame
  height and 0.29% of its pixels. Backing in does not fix it — at this field of view a unit
  subject fills the frame at a distance of 0.35, inside its own bounding box. Distance sets how
  strong the perspective is and that was never the problem. `Camera::fit` sets the **focal
  length**, the projection is linear in it, and one pass is exact. 0.29% to 2.52%, with a cube, a
  plate and a portafilter all landing their furthest corner at 0.850 of the half-frame.

  `--snapshot` takes `--frame N` now, because a run that fills up over its length has nothing in
  it at `t = 0` and a snapshot of that still counts as "the renderer works".

- **The stream out of the spout is drawn as a stream.** A parcel falls the 75 mm from basket to cup
  in 0.12 s against a 0.7 s frame, so sampling parcels there caught roughly none: the first version
  drew the basket beautifully and had nothing between the spout and the cup. What is there is a
  continuous jet, coloured by the concentration leaving the basket — the domain's number, not the
  parcels'.

- **`portafilter_flow`**, which is the picture the domain was built to be able to draw: a shower
  screen, a basket, a body and a spout, with parcels of water leaving the screen, working down
  through the grounds and coming out the bottom darker than they went in.

  Nothing places a streamline. A parcel is advected by the **pore** velocity of the solved Darcy
  field — `u/ε`, not `u`, a factor of 2.2 — and its colour is what it has picked up from the cells
  it crossed. Two checks make that a readout rather than a drawing: the mean transit comes out at
  **0.981×** the closed form `εL/u`, and the load a parcel arrives with agrees with the domain's
  own outlet TDS to **0.9996×**. A Lagrangian tracer and an Eulerian reading, agreeing.

  The hardware and the water share one panel, so they occlude each other properly. What is solved
  and what is only drawn is stated in a table rather than left to look the same.

  With the wall gap the ring's water is through in 4.8 s against the core's 14.4 — **3.01× faster**
  — and arrives carrying 18.0 kg/m³ against the even bed's 82.9. The same water, a fifth of the
  coffee.

- **Two defects the picture found that no test had.**

  **`repack` did not re-solve.** The constructor does, because a quasi-static field read before its
  solve is a field of zeros wearing the shape of an answer. Leaving the solve merely *stale* is
  worse: the field left behind is the previous answer — smooth, bounded, the right order of
  magnitude, and wrong by 42%. Widening the ring to 0.60 left `flow_rate` reporting the even bed's
  flow to the last digit, 1.0000 where the closed form says 1.4188. Every existing test stepped the
  puck before reading it, and a step re-solves. Every mutator ends solved now, and
  `a_mutator_leaves_the_flow_solved` checks all three against their closed forms.

  **An out-of-range probe returned zero, which made the outlet an attractor.** RK2's midpoint lands
  further along than the parcel is, so a parcel one cell from the outlet probes past it. Zero there
  parks it a hair inside the last cell — measured, 19.979 mm of a 20 mm bed — for the rest of the
  run. Nothing about that looks wrong: the stream is steady, the colours are right, the flow rate is
  right, and the water never reaches the cup. Visible only because the example counts arrivals and
  the count was zero.

- **`Substance::stainless_304`, and a basket made of what a basket is made of.**

  `Basket::espresso` used aluminium because the catalogue had it. Steel conducts a tenth as well —
  16.2 W/m/K against 167 — while holding **more** heat per unit volume, 4.0 MJ/m³/K against 2.4:
  a better reservoir and a worse spreader, which is why a group head is brass and a basket is not.

  The wall sets the explicit solver's step, and the limit goes as the diffusivity, so steel's is
  seventeen times larger. `dualis-porous`'s suite went from **152 s to 5.8 s** in release and 52 s
  to 6.7 s in debug — the mode CI actually runs — and the model got more accurate in the same
  edit. It had made the debug matrix jobs ten times slower than the whole workspace had been.

  Reaching for the metal already in the catalogue cost an order of magnitude in run time *and*
  understated the thermal mass by 65%.

- **A seventh domain: `dualis-porous`.** Flow through a packed bed, the heat it carries, and the
  dissolution that rides on both. An espresso puck, and also a filter, a catalyst bed, a leaching
  heap and an aquifer.

  Darcy's law is elliptic — `∇·((k/μ)∇p) = 0` — which is the same operator `Conductor` solves for
  electric potential with mobility in place of conductivity. Reproduced exactly: a uniform bed
  gives `Q = kAΔp/(μL)` to 3×10⁻¹⁵, the way a uniform block gives `ρL/A`.

  On top of it, two things that are marched rather than solved: heat advected by the liquid and
  conducted through the bed, and solute dissolving out of the particles and carried away.

  **Grind is two lengths, not one.** Extraction is diffusion out of a particle, so it uses the
  sieve diameter: `K = 4π²D/d²`. Flow is not — a real grind is a coarse mode plus a tail of fines
  that lodge in the gaps and carry the pressure drop, and Kozeny–Carman from the sieve diameter
  says an espresso puck passes twenty litres a second. So the permeability uses a *hydraulic*
  diameter, 1/160 of the sieve one, and that ratio is one of exactly two fitted numbers.

  Those two `d²` point opposite ways, which is why espresso is hard: finer is less permeable
  *and* extracts faster. Neither is coded — they fall out of the two being derived separately.

  Seventeen tests, every one against a closed form, an exact limit or a conservation law.

- **The two things the first version of that crate got wrong**, because each was invisible until
  something was drawn or weighed.

  **Extraction had no equilibrium term**, so `dm/dt = −K·m` depended on nothing outside the
  particle: every cell at the same temperature extracted at the same rate however much liquid
  passed it. A bed with a channel through it then extracted *perfectly evenly*, which is the
  opposite of what a channel does. The driving force is a difference now — `dm/dt = −K(m −
  m₀C/C_sat)` — and a cell the flow rushes past is kept dilute and keeps going while a cell the
  flow avoids fills up its own pore liquid and stalls at `1/(1+β) = 0.577`. That ceiling is the
  entire reason a channelled shot under-extracts.

  **The dose was 52 g in an 18 g basket**, because `solid_density` was coffee's *skeletal* density
  of 1400 kg/m³. A ground particle is itself porous; its apparent density is 600, and with an
  inter-particle porosity of 0.45 that is a 330 kg/m³ puck, which is what 17.6 g in 20 mm of a
  58 mm basket weighs. Three densities, and the middle one is the one this model wants.

- **The statistic that looked like a channel detector and was not.** The obvious measure is the
  spread of the per-cell extraction. An evenly packed bed already sits at 0.105 — water that
  entered clean is loaded by the time it leaves — and a wall gap that halves the yield takes it
  only to 0.128. The signal is a fifth of the noise it rides on.

  `Puck::radial_contrast` divides the axial gradient out: the ring and the core span the same
  depths. It is 1.0000 on an even puck and 1.20 with the gap.

  Nor does the peak extraction rise, which is the story everybody tells. It **falls**, 0.936 to
  0.836, because the channelled basket reached the same weight in 15 s instead of 25. "The channel
  over-extracts" is about the ring relative to the core, and at equal weight the absolute numbers
  go the other way. Reported rather than asserted, so the test does not encode a plausible story.

- **"Finer grind, higher yield" is false at equal time.** Measured: at 25 s a 175 µm bed reaches
  12.9% against a 350 µm bed's 20.2%, its liquid sitting at 10.7% TDS and going nowhere, because
  four times less water crossed it. Pulled to the same *weight* the statement holds and both
  mechanisms show at once — 4.000× the time, exactly the permeability ratio, and 24.1% against
  11.0%.

- **`espresso_shot`**, which is what all of that is for. A basket from the pump to the cup:
  17.6 g in, 38.0 g out, 24.6 s, 19.6% yield, 8.3% TDS. Grind, temperature and pressure each swept
  on their own against the exponent each carries — time to a fixed weight goes as `1/d²` to 3×10⁻⁴
  and as `1/Δp` to 10⁻¹³. Then a gap at the wall, and a portafilter left on the counter.

  `espresso_shot shot.html` draws a vertical cut through three baskets in nine fields as the shot
  runs. A cut and not a volume: nine volumes on fifty frames is 45 MB to animate three baskets
  that are not moving, and the question was about a cross-section.

- **`dualis_scene::sample_field` is public.** One field per domain is not enough and the limit is
  the trait's rather than the physics'. A bed under flow has a temperature, a pressure, a speed,
  an extraction state and a concentration on the same grid, all true at once, and `as_field` can
  nominate one. A caller that wants the other four builds them here.

  Which one `Puck` nominates is a decision with the same shape. It is **extraction**, not
  temperature: a bed is isothermal unless somebody deliberately cooled the basket, so a
  temperature panel is a flat rectangle on every ordinary run — a picture that renders, looks
  fine, and carries nothing.

- **`Basket` replaces ten positional arguments.** Two of them were `Length` and two were
  temperatures, so transposing a pair compiled and ran; `clippy::too_many_arguments` had been
  silenced to allow it. `Basket::espresso()` is a conventional double basket and `..` covers the
  nine things you are not asking about.

  It also made the geometry statable. The basket radius was "the largest circle that fits", which
  puts the metal in the *corners* of the grid — the right heat capacity in the wrong shape, and a
  cut through the axis crosses none of it. That was invisible until the first cross-section came
  out flat.

- **`scenes/18-an-espresso-shot`**, so all seven domains have scenes, and `DynamicViscosity`,
  `MassFlow` and `Concentration` in `dualis-units`.

- **The scene format carries a version**, and `--check` validates a file without running it.

  `Scene::format`, where **absence means 1** — what every scene written before the field existed
  is, so old files are readable by construction rather than as a special case. A version this
  build cannot read is refused rather than half-run.

  `deny_unknown_fields` already catches a key that was *added*. It cannot catch a key whose
  **meaning changed** — same name, same type, different semantics — and that is the gap a version
  number closes. The promise attached is narrow and stated: within one version, a file that loads
  today loads tomorrow.

  `dualis-world --check scene.json` parses and builds without running, and reports a parse failure
  as `file:line:column` with the keys that were expected. CI runs it over every shipped scene, so
  it is not the one entry point nothing exercises.

  This was the editor's blocker: an editor writes files, and the format had changed under its own
  users once already, silently, when `mode` became `release`.

- **`runtime/gpu`**: `Solid3D`'s seven-point stencil as a WGSL compute shader, implementing
  `Domain` so it drops into a `Simulation`. Its own workspace, and the rule is that **the CPU
  domain is the reference and this is a cache of it**.

  Measured, 400 steps: 3.5× at 16³, 24× at 32³, 85× at 48³, **191× at 64³**. The GPU column is
  flat at ~0.055 s — it is bound by dispatch overhead, not the stencil, so what the table really
  shows is the CPU's `n³` growing away from a constant.

  WGSL has no `f64`, so this is a *different computation*, not a faster one. It conserves to
  `5.0e-11` where the CPU holds `9.1e-15`, which is below `Simulation`'s default `1e-9` audit — a
  scene using it must loosen `conservation_tolerance_for(ENERGY, ..)`, and `GpuSolid` declines
  `books_balance` for the same reason.

  Reductions stay on the CPU. A mean summed with atomics depends on which workgroup finished
  first, and addition is not associative, so `ledger` reads the grid back and sums in index order.

### Fixed

- **Single precision was never the problem; spending it on an offset was.** The GPU buffer first
  held absolute kelvin and diverged from the reference by `1.4e-3` after two hundred steps — a
  thousand times what accumulation predicts.

  The update is `centre + F·(sum − 6·centre)`. Near 293 K that `sum` is about 1759, where `f32`'s
  resolution is `1.2e-4`, and the difference being extracted is of order `1e-3` K. Subtracting two
  numbers that agree to five digits keeps **less than one digit** of the answer, every step.

  The buffer holds `T − T₀` now, where the same numbers sit near 1 K and the subtraction keeps
  about four digits. Divergence `1.449e-3` → `8.7e-7`, conservation drift `7.4e-7` → `1.2e-10`.
  The stencil is linear, so subtracting a constant commutes with it exactly and the fix cost
  nothing.

- **`dualis_view::gltf`**: a frame as glTF 2.0, so Blender, three.js, Omniverse, a USD pipeline
  or macOS Quick Look can open a result. **No new dependency** — glTF is JSON with the binary
  base64'd into a `data:` URI, and this crate already writes JSON by hand, so `dualis-view` still
  depends on exactly one thing and it is `dualis-scene`.

  Paths become `LINES`, points and 3D fields become `POINTS`, and a 1D or 2D field becomes
  **nothing, reported**: a line of samples is a graph, not something to put in a 3D viewer, and
  `Exported::skipped` says which panel and why. Eight of the seventeen shipped scenes export
  geometry; the other nine are refused with a reason rather than written as an empty file.

  One frame, not the run. glTF animates node transforms and morph targets — a thing moving or a
  mesh deforming between fixed vertex counts — and a field whose values change is neither.
  Encoding a run as an animation would mean choosing a lie about what is moving.

  Checked against what a loader enforces rather than against whether it looks plausible: buffer
  length against the decoded blob, accessor counts against the geometry, indices read back out of
  the bytes and range-checked, four-byte alignment on every view, `min`/`max` on every `POSITION`,
  and the hand-written base64 against the RFC 4648 vectors.

- **`runtime/viewer`**: a native wgpu window for a run — rotate, zoom, scrub — in **its own cargo
  workspace**, excluded from the library's. Measured: the library resolves 12 external crates, the
  python bindings 15, and this 86.

  It does **not depend on `dualis`**. It reads the JSON a run wrote and nothing else, so "the wire
  format carries enough to draw a run" is demonstrated rather than claimed.

  `viewer-core` holds everything a renderer gets wrong the same way twice — one colour scale
  across the run, one framing across the run, and a projection that clamps a point behind the eye
  instead of turning it into a streak — with seven tests against real run files.
  `--snapshot out.ppm` renders one frame headlessly, because a window nobody can photograph only
  proves the program did not panic.

- **`PanelData::Paths`** in `dualis-scene`, and a **`layout`** view in `dualis-view`: runs of
  connected points in space — a ray through a lens train, a trajectory, a field line — drawn
  rotatably and depth-sorted beside the field and body views.

  The third shape, and it took an optical bench to need it. A field is defined everywhere and a
  body is somewhere; a **path** is a thing that went from one place to another, and drawing a
  traced ray as a scatter of its vertices loses the one property that makes it a ray.

- **`optical_bench`**, an example that draws the *instrument* rather than a graph. A doublet, a
  fold mirror turning the axis through 90°, three field angles, and an image plane — prescribed,
  traced, refocused, and then **bent** until the spot lands inside the Airy disc:

  ```text
    effective focal length   97.8 mm   measured from a traced ray, not assumed
    RMS spot, as prescribed  209.2 um  47x the diffraction limit
    RMS spot, bent           0.57 um   0.13x — diffraction-limited
  ```

  `cargo run --release --example optical_bench bench.html` gives a layout you rotate in a browser.

- **`busbar_rating`**, an example shaped like an engineer's working day rather than a
  demonstration. A bolted busbar joint, geometry to production yield:

  ```text
    contact resistance   3.440 uohm, of which 37% is the joint itself
    thermal path         0.0140 W/K, from the network's solved balance
    continuous rating    445.8 A  (2.23 A/mm2, 105 C limit, 40 C ambient)
    thermal runaway      1018.8 A — 2.29x margin
    yield at nominal     36.6% of 20 000 units
    derate for 99.9%     384.6 A, 86% of nominal
  ```

  Every step has a closed form behind it: `rho L/A` for the bar, Maxwell's `rho/2a` for the
  constriction as a limit the solve is shown converging on, `dT = I²R20/(g − I²R20·alpha)` for the
  electro-thermal fixed point, and that expression's pole for the runaway current.

  The finding it is built to make is the last two lines. **A rating computed from nominal values
  is a coin toss in production** — 36.6% here — and the derating that fixes it is what the Monte
  Carlo is for.

- **Two 3D examples**, `heat_in_three_dimensions` and `room_in_three_dimensions`, run by CI like
  the rest.

  The first is built on the closed form for an instantaneous point source: the peak falls as
  `t^(-d/2)`, and **that exponent is the dimensionality**. A bar gives `-1/2`, a plate `-1`, a
  block `-3/2`. Fitted at `-1.514` inside the window where the source is still a point and the
  block still looks infinite, against `-1.362` before it and `-0.245` after — so the window is
  demonstrated rather than asserted.

  The second is `room_modes` with a ceiling: the floor-to-ceiling mode at 71 Hz that a floor plan
  does not have *at all*, and a mode count checked against Weyl's three-term estimate.

- **`Solid3D`, `Reading`, `Bodies` and `Tolerances` in the prelude.** Each was reachable only
  through `dualis::core` or its own crate, and each is a type a consumer meets while building a
  frame, setting an audit or writing a domain. Found by examples reaching for them.

- **A `volume` view** in `dualis-view`: a 3D field is raycast — trilinear sampling, front-to-back
  compositing, rotatable with the same camera the bodies view uses — **beside** the slice montage
  rather than instead of it. A render shows shape and a reader cannot get a number back out of it;
  a montage is quantitative and unreadable as a shape. `ARCHITECTURE.md` gap 6, and the answer
  turned out not to be a depth buffer.

  The opacity transfer function is chosen from the run's own range: transparent in the middle for
  a signed field, or a standing wave renders as a solid block; transparent at the low end for a
  one-sided one, or a block at ambient does the same for the opposite reason.

  When a feature occupies less than 3% of the frame the caption **says so**, with the figure. A
  single hot cell in a block of 729 is a small bright dot and everything else is transparent,
  which is correct and reads exactly like a broken renderer.

## [0.10.0] — 2026-08-10

Three domains gained a third dimension and electricity gained a field, so
`ARCHITECTURE.md`'s six gaps are down to one. **Breaking**, in two places, both in
`dualis-scene`: `Extent::new` takes an `nz`, and `Panel::grid` returns a triple. Everything
else is additive.

### Added

- **`Solid3D`** in `dualis-thermal`: conduction through a block in three dimensions. A seven-point
  stencil on cubic cells, insulated faces by mirroring, and the explicit limit `dx²/6α` — a third
  of `Bar1D`'s, because the limit tightens with every axis.

  The first 3D field domain, and the thing a bar cannot do: heat spreading *sideways* out of a
  hot spot is the whole job of a spreader plate, and a one-dimensional model has nowhere for it
  to go but along.

  Checked against the **exact eigenvalue of its own discrete operator**, not against a second
  implementation. A separable cosine mode on a cell-centred grid with mirrored faces decays by
  precisely the same factor every step, so the test is an equality at machine precision. Then the
  continuum: that discrete rate approaches `α·π²/L²` at second order, checked as a *rate*, since
  a first-order scheme also converges.

- **`Domain::books_balance`** in `dualis-core`: an opt-in claim that a domain's ledger changes by
  exactly what it took from the bus minus what it published. A domain that makes it is checked
  **on its own scale** every step, rather than inside the sum of every ledger.

  The failure that closes, demonstrated in `per_domain_books.rs`: a domain holding a microjoule
  beside one holding a kilojoule loses a fifth of itself, and the total moves by `2e-10`. No
  tolerance catches that — tightening to `1e-12` refuses the run for floating-point noise long
  before it can see a leak of that shape, because the problem is the scale and not the number.

  Opt-in because not every honest ledger is an exact one. `LumpedMass` loses heat to an
  environment that is not on the bus; that is a boundary being modelled, and a check that accused
  it would be the wrong check. Every other domain in the workspace takes the claim and passes.

  `Exchange::traffic` and `Exchange::total_published` are what make it attributable: the scheduler
  visits domains one at a time, so the bus traffic between the snapshot before and the snapshot
  after belongs to exactly one domain.

- **`Tolerances`** in `dualis-core`, and `Simulation::conservation_tolerance_for`: a relative
  tolerance **per conserved quantity** rather than one for the whole simulation.

  The failure it closes, demonstrated in `per_quantity_tolerances.rs` in both directions: a
  Barnes-Hut tree gives up exact momentum by construction, so at `1e-9` a correct run is refused;
  loosen to `1e-6` and a real energy leak passes. A quantity's achievable accuracy is a property
  of the scheme carrying it, and different quantities in one simulation are carried by different
  schemes.

  `audit_with` is the per-quantity form; `audit` keeps its signature and delegates, so no existing
  caller changes. `Violation` now carries the tolerance that *actually applied* rather than the
  default. A `BTreeMap` inside, because a violation's message must not depend on the order a
  builder was called in.

  A scene can set them too, as `tolerance_for`, and a channel name the kernel does not have is
  **refused** rather than ignored — the same failure `aluminum` for `aluminium` produced in this
  format once already, where one character turned off the check the library exists for.

- **`Conductor`** in `dualis-electrical`: current as a field. `∇·(σ∇φ) = 0` solved by conjugate
  gradients on a block with two electrodes, `J = −σ∇φ` read off it, and the dissipation as
  `∫σ|∇φ|²dV`. **Nobody states a resistance** — it comes out of the shape.

  For a uniform block it comes out as `ρL/A` to machine precision, which is what makes it
  checkable; for a notched one it comes out as whatever the notch gives, which is the point.
  Series adds resistances and parallel adds conductances, and neither is coded — both fall out of
  the same solve, on materials four orders of magnitude apart so a face conductivity that used an
  arithmetic mean instead of a harmonic one would show.

  `V·I` equals `∫σ|∇φ|²dV` to machine precision — Tellegen's theorem, and the sharpest single
  statement that the discretisation is self-consistent, since the two are different sums over
  different things.

  The first **elliptic** domain here, and the first whose failure mode is a solver rather than a
  stability limit: a solve stopped at its cap returns a field that is smooth, bounded and shaped
  exactly like an answer. `step` refuses one that did not converge, `residual` is a *reading*, and
  `with_solver` exists so the refusal path can be provoked from a test.

- **`Resistivity`, `Conductivity`, `ElectricField`, `CurrentDensity`** in `dualis-units`, with
  `product!` declarations for `J = σE` and `E = ρJ` — those lines compiling is the check that
  (S/m)·(V/m) is A/m².

- **`DomainSpec::Conductor`** and scene 17 — a copper busbar with a notch. Seventeen scenes.

- **`Hall`** in `dualis-acoustic`: the wave equation in three dimensions. A staggered grid with
  pressure on nodes and velocity on the faces between them, rigid surfaces, and `dx/(c√3)`.

  It is not a more accurate `Room`. A floor plan **does not have** the vertical and oblique modes
  — not less accurately, at all — and a 2.4 m ceiling puts the first one at 71 Hz. The mode count
  also grows as `f³` rather than `f²`, which is why a real room's resonances merge into a hiss
  where a two-dimensional model keeps them separable much further up.

  Checked against the rigid-wall mode frequencies, which are exact, and against a second-order
  convergence rate measured **across three doublings** — 13 to 97 nodes falls 54.6×, against 64×
  for second order and 8× for first. One doubling was tried first and proved nothing: the
  per-doubling ratios bounce between 2.9 and 5.6, because "worst departure over a run" is a
  maximum and therefore noisy. `Room`'s own convergence test reaches the same conclusion by the
  same route.

  It carries the leapfrog startup fix from its first line rather than inheriting the `O(h)` defect
  `Room` and `Tube` shipped with, and a mutation confirms three separate tests would catch it.

- **`DomainSpec::Hall`** and scene 16 — the same 4.4 × 3.1 m room with a 2.4 m ceiling, released
  in its oblique (1,1,1) mode. Sixteen scenes.

- **`Extent` and `PanelData::Field` gained a third axis** in `dualis-scene`. Breaking:
  `Extent::new` takes `nz`, `Panel::grid` returns a triple, `PanelData::Field` carries `nz`.
  `Extent::volume`, `Extent::count`, `Extent::dimensions` and `Panel::slice` are new.

- **`DomainSpec::Block`** in `dualis-world`, and scene 15 — a hot spot in a 9×9×9 aluminium block.
  Fifteen scenes now, all run by CI through the real binary.

- **A `slices` view** in `dualis-view`: a 3D field is drawn as every z-slice at once rather than
  one plane behind a slider, because a viewer who never touches the slider would see a picture of
  a solid that was really a picture of one plane. The filmstrip has no room for a montage, so it
  draws the middle slice and labels it `z-slice 5/9`.

### Removed

- **`Solid3D::as_bodies`.** It existed as cover for the capture gap below — a way to get a block's
  cells out as a point cloud when a field would have come back as one slice. With `Extent` now
  three-dimensional the cover is unnecessary, and it was never free: a domain that is two shapes
  at once makes the picture depend on whether somebody remembered to set an extent, which is a
  mode nothing announces. It is a field, and only a field.

### Fixed

- **`Solid3D::max_stable_dt` was documented as a limit and read as a recommendation.** At exactly
  `dx²/6α` the sharpest mode the grid can hold has an amplification factor of `-1`: marginally
  stable, so it flips sign every step and never decays. A point source excites it as hard as
  anything can, and the peak comes out **1.96×** the closed form — from a scheme that never
  diverges and whose conservation audit is exact to the last bit. At half the limit it is 1.005×.

  Nothing is wrong with the limit; it is a stability limit and stability is all it claims. The
  documentation now says so, and `heat_in_three_dimensions` runs at half of it with the numbers
  for all three cases in its header.

- **A quasi-static domain reported an answer before it had one.** `Conductor::new` left the
  potential at zeros, so the first captured frame reported a resistance 24× below the floor
  `ρL/A` puts under it — beside a residual of `inf` that nothing was reading. A quasi-static
  domain has no state before its solve, so it solves at construction now.

- **A three-dimensional field was captured as its `z = 0` face**, silently. `Extent::samples` was
  a pair and the sampler built its position as `(u, v, 0)`. For six domains that was exactly
  right; for the seventh it produced a 9×9 plane of a 9×9×9 block — a perfectly plausible picture
  of a block, two thirds of the samples missing, nothing anywhere to say so.

  `FRICTION.md` 23, and the lesson is not about `Extent`. A layer's assumptions are only visible
  from below: `dualis-scene` names no domain and succeeded at that, and could not have discovered
  that it assumed flatness, because everything it had ever been handed was flat.

- **A test that was green and measured nothing.** Replacing `Solid3D`'s stability constant `1/6`
  with the one-dimensional `1/2` left nine of ten closed-form tests passing — none of them ever
  excited a mode sharp enough to care. `the_limit_is_where_the_sharpest_mode_stops_growing`
  releases the fastest-alternating mode the grid can hold, steps it sixty times at exactly the
  reported limit, and measures. It also asserts the limit is *marginal*, since a scheme that
  damped that mode comfortably at its own limit would be leaving stability unused.

## [0.9.0] — 2026-08-10

### Added

- **`dualis-view`**, the eleventh crate and the top of the three layers: a filmstrip as SVG, a
  self-contained HTML report, a CSV of every domain's scalars, and the frames as JSON. No
  dependencies — SVG and HTML are text, so a `format!` and a file write is the whole renderer.

  **The view is chosen by the shape of the data**: scalars over time become a chart, a 1D field a
  profile, a 2D field a heatmap, points in space a rotatable scene. Its tests are driven by frames
  written out by hand rather than by a simulation, which is the only way to tell "a heatmap
  because the data is a 2D grid" apart from "a heatmap because that domain was a room".

  Every view holds one scale for the whole run. A picture that renormalises per frame makes a
  decay look like a steady state, and it is what you get if you do not think about it.

- **`dualis-scene`**, the tenth crate and the middle of the three layers `ARCHITECTURE.md`
  describes. `Placement`, `Extent`, `Frame`, `Panel`, `PanelData`, `capture`, `settle_framing`
  — where a domain sits, and what one instant of a run looks like.

  It **names no domain**, and `knows_no_physics.rs` demonstrates that rather than asserting it:
  the test defines a physics inside the test file — a field, two bodies and a reading, in a crate
  `dualis-scene` cannot possibly know about — places it and captures all three shapes. If a
  physics invented in a test comes back whole, a real one costs one crate and nothing else moves.

  Both crates were modules in `dualis-world`, which is `publish = false`. A consumer who could
  state a simulation and run it could reach neither the shape of the answer nor any view of it —
  the largest gap between what was built and what was usable.

- **`Domain::readings`** and **`Reading`**: the named scalars a domain has when it has no picture.
  Eight of the fourteen shipped scenes contain a domain that draws nothing — a heater, a lamp, a
  winding, a thermal network — and for several the scalar *is* the result.

- **`Domain::as_bodies`** and the **`Bodies`** trait: count, position, a value to colour by, and a
  *real* wall or `None`. The counterpart to `as_field`, which covered only the domains that are
  continua. `FRICTION.md` finding 11, recorded and unfixed for months, and paid the moment the
  layers were separated: a scene layer that must name three physics to find out where anything
  *is* needs editing every time a fourth arrives.

  The trait draws a line the old code could not. A periodic cell is a boundary condition and the
  domain reports it; an orbit's box is a property of the picture, and nothing physical sits at
  its edge, so a view measures that one over the whole run instead of being told.

- **`Simulation::domains`**: enumerate what is in a simulation. There was no way at all — a
  caller could ask for a domain *by name*, which is no use to a layer that must visit every one.

- **`ScalarField::unit`**: two characters for a legend, and the fifth place a layer had been
  matching on domain types to get something a domain already knew.

- **`Pose`** in `dualis-core`: a rigid motion, rotation and translation, no scale and no shear.
  An isometry preserves every distance and angle exactly, which is the only class of placement a
  physics can be moved by without its physics changing. The first test is
  `placing_something_cannot_change_a_distance`.

### Fixed

- **A second domain could empty a channel another had already emptied, and the audit could not
  see it.** `Exchange` counted takes but nothing compared the count across a turn, so two
  consumers of one channel each reported a consistent ledger while the amount was delivered
  twice. `Simulation::sweep` now compares takes per channel across the turn and raises a
  `Violation` naming the channel.

- **`Bar1D`'s field was labelled `"C"` and returns kelvin.** The application converted before
  drawing, so the offset and the label were applied in the same expression and nothing could
  disagree with anything. The field now says `"K"`, which is what the cells hold; the celsius a
  picture wants is a view's conversion, and `FRICTION.md` 22 records that the library gives it
  nowhere to live yet.

## [0.8.0] — 2026-08-09

### Added

- **`Ensemble`** in `dualis-core`: many independent samples, run in parallel, with an answer that
  does not depend on how many threads produced it. The other axis of parallelism — `TreeNBody`
  splits one evaluation across cores, this splits many evaluations, which is the shape a Monte
  Carlo study and a parameter sweep both have. Measured 8.74× on sixteen threads.

  It is bit-for-bit across thread counts because two existing decisions meet: `Rng::for_index` is
  stateless and index-addressed, so sample `i` draws the same numbers wherever it runs, and
  results land in a slot chosen by index rather than being appended. The failure this avoids is
  the usual one and is nearly undetectable — a Monte Carlo drawing from a shared generator gives
  a different answer per core, and the difference looks exactly like statistical noise.

  `estimate` folds in **fixed-size blocks**, so a study is bounded by its block count and not its
  sample count: ten million samples in kilobytes, tested. The block size is fixed rather than
  derived from the thread count on purpose — a per-thread split combines a different number of
  partial sums on four cores than on sixteen, and floating-point addition is not associative.
  `Ensemble::blocks` is public so a caller building a histogram or a quantile has the same
  discipline available.

  Welford within a block, Chan's merge between them, rather than `sum(x²) − n·mean²` — which
  subtracts two large nearly-equal numbers and loses every digit exactly when a Monte Carlo has
  converged. Pinned against a case with an exact answer.

- **`Estimate`**, with `mean`, `standard_error`, `samples`, and `standard_deviation()` kept
  distinct from the error on the mean, because confusing those is the usual way to misreport a
  Monte Carlo result.

### Changed

- **`Fluid` is about 1.95× faster**, in changes that move no result: the Lennard-Jones potential's
  loop invariants hoisted out of the pair loop, the force quotient taken once instead of twice,
  and the periodic wrap skipped where it provably does not apply. Verified on four platforms
  against a pinned digest. Two techniques that should have helped did not — a cell-ordered copy
  of the positions is *slower*, because the counting sort is stable in index order and the reads
  were never a random gather.

- **`detector_snr`** runs on `Ensemble`, which also fixes a variance it was computing as
  `sum(k²)/N − mean²` — at a mean of 900 that subtracts 1.6e11 from itself to reach 900.

- **`where_the_time_goes`**, a new dependency-free example, because this workspace had never
  measured itself and every claim about which loop mattered was a guess. It takes the best of
  five trials: consecutive runs vary by 8%, which is wider than several differences that were
  nearly reported as wins.

## [0.7.0] — 2026-08-09

### Added

- **`ThermalNetwork::path_conductance(node, at)`** — the conductance of the whole heat path from
  a node to ambient, as the slope of its own solved balance. Exact and operating-point
  independent when nothing radiates; the local slope when something does, which is the right
  answer because everything asking for this is asking a derivative question.

- **`Volume::cm3`/`mm3`/`m3`/`litres` and `Area::cm2`/`mm2`/`m2`.** Building a three-node network
  was six lines of `Volume::from_si(x * 1e-6)`, because the constructors take dimensioned types
  and the numbers a person has are cubic centimetres.

### Fixed

- **A number this workspace was quoting was wrong, and the fix is an API rather than an edit.**
  `runaway_current`'s documentation said a motor's threshold falls from 4.95 A to 4.11 A once the
  joints are counted. The 4.11 is a *convection-only* path: the real one is 0.220 W/K rather than
  0.203, because the housing also radiates at its operating temperature, and the true threshold
  is **4.28 A**. The hand-assembled formula understated the margin by 4%.

  Found by a sizing tool written against the published 0.6.0 — a consumer deliberately unlike
  `dualis-world`: no scenes, no rendering, no fields, asking for settled answers rather than
  stepping. It had to assemble that conductance out of numbers the network already held, which
  is `FRICTION.md` 20 and is what `path_conductance` now answers.

  It also cross-checked something worth keeping: the tool's fixed-point iteration on
  `steady_state` lands the winding at 99.0 °C, and scene 13's marching with a between-frames
  feedback lands at 99.02 °C. Two unrelated routes to the same coupled answer.

## [0.6.0] — 2026-08-09

### Added

- **`Domain::as_any_mut` and `Simulation::domain_as_mut`.** A caller could read a domain and not
  write one, which made a whole class of coupling closable from *nowhere*: not inside the step
  loop by design, and not outside it by omission.

  The case is a copper winding whose resistance rises 0.393%/K. Its temperature lives in
  `dualis-thermal`, its resistance in `dualis-electrical`, and neither can see the other's state
  — correctly, since domains meeting only on `Exchange` is the property the crate split defends.
  The caller between frames can see both.

  This does **not** weaken that rule: it is about what happens inside `step`, and this runs in
  code holding `&mut Simulation` that could drop the domain and rebuild it. It is also
  deliberately *not* a state channel on the bus, which a true in-loop coupling would need and
  which stays undecided.

  `as_any_mut` defaults to `None`, so a domain that forgets it is silently unwritable — the
  opt-in hazard of `FRICTION.md` findings 7 and 12, handled in the same change rather than
  rediscovered: all twelve implementors got the counterpart beside `as_any`, including the
  `Box<dyn Domain>` forwarding impl. `FRICTION.md` 18.

- **`Winding::dissipation_at(T)` and `resistance_at(T)`** — pure functions rather than `Domain`
  methods, so `P(T)` composes for whoever holds both sides. `dissipation()` is now
  `dissipation_at(its own temperature)` rather than a second copy of the arithmetic, checked on
  `to_bits()` at four temperatures.

- **`Winding::runaway_current(g)`**: `√(g/(R₂₀α))`, where `dP/dT` overtakes `dQ_out/dT`. The test
  measures the slope from two dissipations a kelvin apart and asserts the inequality *flips*
  across it, rather than reproducing the formula from itself. `None` for a voltage drive, which
  cannot run away because `V²/R` falls as it warms.

  `g` is the **whole** path to ambient. A winding reaching air through 0.9 and 2.4 W/K of joints
  and then 0.294 W/K of convection has a series conductance of 0.203, and the threshold falls
  from 4.95 A to 4.11 A — 17% of margin a lumped model reports as present when it is not.

- **Scene 13**, the feedback closed between frames and measured. The amplification is `1/(1−g)`,
  checked as a ratio against the same scene without it: 1.281 measured. Convection alone predicts
  1.310; including the housing's linearised radiative conductance at its operating point gives
  1.280, so the 2.2% is radiation stiffening the heat path rather than error.

## [0.5.0] — 2026-08-09

### Added

- **`dualis-electrical`**, the sixth domain and the tenth crate. `Winding` computes
  `R = ρ(T)·L/A` and publishes `I²R` onto the channel `dualis-thermal` already takes from, with
  neither crate naming the other.

  It closes a real gap rather than adding a sixth for its own sake. Every other producer of heat
  here answers a question about *something else* that happens to warm a thing — light landing on
  a mirror, a dashpot damping a bounce. A winding is the case where getting hot is the entire
  subject, and until now the workspace's own examples stood a stated number of watts in its
  place. A stated number cannot be wrong, which is another way of saying it is not a model.

  Two mistakes it made, both caught by something other than its author. It first declared its own
  `HEAT = "heat"` channel, which reads correctly and is a *different* channel from
  `quantity::ENERGY` — so it published joules nothing consumed, and the audit named it on the
  first step with the amount. And an infinite reserve turned out not to *fail* the audit but to
  **disable** it: `inf` before, `inf` after, `inf` equals itself, and a winding pouring joules
  into a plate runs green at any tolerance. `Winding::step` refuses that itself, because the
  audit structurally cannot.

  The electro-thermal feedback that causes thermal runaway is deliberately **not** expressible:
  a domain would have to read another's temperature inside the step loop, and `Exchange` carries
  amounts, not state. Resistance is evaluated at a temperature the caller states.

- **`Resistance`** in `dualis-units`, with `product!(Resistance, Current => Voltage)` — that line
  compiling is the check that ohms times amperes are volts. Plus `Current::a`, `Voltage::v` and
  `Resistance::ohm`/`milliohm` constructors.

- **Scene 12**, `12-winding-heats-a-motor`: the same motor as scene 11 with the watts computed
  from a length of wire rather than stated. The two settle within a fifth of a kelvin of each
  other, so the guess was good — and this is the scene that would have caught it if it had not
  been. Twelve scenes now cover all six domains.

- **Two tests that keep prose from aging.** `documented_version.rs` reads every `dualis = "x.y"`
  in the documentation and every "the tree is X.Y.Z" in `.claude/agents/`, and compares them
  against `CARGO_PKG_VERSION`. Both fail if they find *nothing*, because a check that stopped
  matching would pass forever. The `invariant-guard` line they now cover had been a release
  behind twice running, in the one file whose subject is checking things.

## [0.4.0] — 2026-08-09

### Added

- **`ThermalNetwork::steady_state(power)`** — where a network settles, solved rather than marched
  to. The question a designer actually asks is *will the winding survive*, and stepping to it is
  slow and approximate: reaching a part in a thousand takes about seven time constants of
  explicit Euler, each accumulating its own error. This solves the same balance the step loop
  converges to, exactly, and returns a `SteadyState` read with the same `Node` handles.

  Not the implicit stepping this workspace declines to have: `Domain::step` is unchanged, the
  kernel is untouched, no schedule learns anything, and the network is not modified by being
  asked. Newton, because radiation makes the balance `T⁴` and a single solve would answer the
  linearised problem — the mistake `LumpedMass::equilibrium_rise` exists to correct on one body.
  The Jacobian's radiative part is the `linearised_loss_conductance` the step limit already uses.

  Refuses a network where no node loses heat to an environment: it warms without limit, so there
  is no steady state, and a plausible finite number would be the worst possible answer.

  Exposed in Python as `steady_state(name, watts)`, where it is worth more than in Rust — a
  marching loop crosses the binding once per step.

### Fixed

- **The Newton bound was set from the wrong measurement and refused a kilowatt.** It was eight,
  twice the worst case the first tests exercised — but none of them loads a node hard enough for
  the `T⁴` term to dominate. At ambient the radiative slope `4εσAT³` is tiny against what the
  balance needs, so the first solve overshoots enormously and Newton walks down at the `3/4`
  ratio a quartic gives: twelve iterations at a kilowatt, sixty-six at a terawatt. Now 100, from
  counting them, with the table in the method's documentation.

  Found by instrumenting the iteration count, which also showed the loop had no way to *report*
  exhausting itself — it returned the last iterate, a plausible temperature for a balance that
  was never struck. It returns a `Violation` now, and the test that would have caught the
  original bound checks a radiation-dominated solve against a root found by **bisection**, which
  shares no arithmetic with Newton.

## [0.3.0] — 2026-08-09

### Added

- **`ThermalNetwork`**, the third domain in `dualis-thermal`: n lumped bodies joined by
  conductances, as *one* domain. Winding, stator, housing — and the **drop across each joint**,
  which is the number that decides whether a motor survives and the one a `LumpedMass` cannot
  give, because it reports the whole assembly as a single temperature. It also expresses a
  *contact* resistance between different materials, which `Bar1D`'s uniform grid cannot.

  One domain rather than a `conducting_to(peer)` on `LumpedMass`, because a conductance carries
  `UA(T₁ − T₂)` and needs **both** temperatures — and domains here meet on an `Exchange` that
  carries amounts rather than state, so neither side can compute the flux alone. Adding that
  method would have broken the property the crate split exists to hold. A network is a single
  coupled system of ODEs, which is what a thermal network physically is.

  Nodes are `Node` handles rather than names, and that is the load-bearing decision. A link
  contributes `+q` to one node and `−q` to another **in the same sum**, so they cancel
  identically and the conservation audit is blind to links *by construction*: a sign error, a
  transposed index or a link dropped altogether passes at machine precision, and the winding
  simply runs at a plausible wrong temperature forever. A handle can only come from a
  constructor, so a dangling link is not representable. `node_named` and `handles()` are the
  bridge for callers building from a file — the JSON scene format and the Python binding both
  resolve names once, at construction, and raise before any stepping happens.

  Seven tests, six of them closed forms, every one per-node or against a formula computed in the
  test file rather than on a total. `n = 1` reduces to `LumpedMass` bit for bit over a 4000-step
  trajectory including `max_stable_dt`, so the new domain inherits every check the old one
  already passes; `linearised_loss_conductance` is shared between them rather than written twice,
  because two copies is the obvious way for that to stop being true.

  Closes #2.

- **`Conductance`** and **`HeatCapacity::j_per_k`** in `dualis-units`. `Conductance × Temperature
  = Power` and `Conductance × Time = HeatCapacity` are declared, and the declarations compiling
  is itself the check that `UA·ΔT` is watts and `C/UA` is a time.

- **`ThermalNetwork` in the Python bindings**: `add_network(name, nodes=[…], links=[…],
  absorbing=…)`, with `node_temperatures`, `node_temperature` and `heat_flow_w` to read it back.
  A node given `ambient_k` without `area_m2` — or the reverse — is refused rather than quietly
  becoming an interior node that looks like it is cooling and is not. `temperature()` refuses a
  network rather than averaging it, and names the calls that answer.

- **The Python bindings are on PyPI**, as `dualis`. `pip install dualis` gets an abi3 wheel for
  Linux x86_64/aarch64, macOS x86_64/aarch64 or Windows x64, plus an sdist to fall back on.
  Built by `.github/workflows/release-python.yml` on a tag, because a wheel built on one
  machine is a wheel for one platform — uploading the Windows one alone would have made
  `pip install` fail on Linux and macOS in a shape that reads as an unsupported platform rather
  than a botched release. Trusted publishing, so no token lives in the repository.

- **A `network` domain in the scene format**, and scene 11: 12 W into a copper winding, out
  through electrical steel and an aluminium housing. The first scene with **nothing to draw** —
  `as_field` declines, because nodes have capacities rather than positions and a conductance is
  not a distance — so the scene test's "produced a panel" guard now takes an explicit list, and
  being on it costs a named check rather than buying a pass.

- **Python bindings**, in `bindings/python`, as their own cargo workspace. `pip install` the
  wheel and `import dualis`: a `Simulation`, the library's heater, bar and lumped-mass domains,
  and the conservation audit as a `dualis.Violation` carrying `quantity`, `site`, `before`,
  `after`, `scale` and `tolerance` — addressable rather than a sentence to parse. A refused step
  does not move the clock.

  SI floats at the boundary with the unit in the parameter name, because the dimensional types
  are a compile-time thing Python cannot have, and a runtime wrapper would cost per operation to
  catch an error a Python caller does not make. What crosses instead is the audit.

  A domain cannot be written in Python yet, and the reasons are in its README. Enough to *run and
  audit* coupled physics, not enough to *extend* it.

  Separate workspace because pyo3 brings about fifteen crates and links libpython, and the
  library's twelve external dependencies, its `deny.toml` allow-list and its WebAssembly jobs are
  promises that should not have to accommodate a Python extension. Verified rather than assumed:
  the library workspace still resolves to exactly twelve external crates. An abi3 wheel, so one
  build serves 3.10 upward; CI builds it, installs it and runs its ten tests, each against a
  number computed in the test file rather than read off the simulation.

### Fixed

- **An `O(h)` bias in `ThermalNetwork`'s steady state that the conservation audit could not
  see.** Heat arriving on the bus was added to the absorbing node's temperature *before* the flux
  snapshot was taken, so that node drove its link from an already-raised value. Explicit Euler
  otherwise reaches a steady state exactly — the fixed point is where the right-hand side
  vanishes, with no step-size dependence — so the joint next to the source sat `K·h/C` low:
  predicted 0.0031006, measured 0.0031005, while the far joint and the environment drop were
  exact to six figures. Every total stayed right throughout, because the excess simply landed in
  the neighbour. Found by the series-resistance closed form, not by the audit. The arriving heat
  is a term of the same right-hand side as the fluxes and now joins the same sum.

- **A NaN check in the Python bindings written as `!(x > 0.0)`**, which rejects NaN by the
  negation being true rather than by saying so. The nested binding workspace is excluded from the
  root one, so the `lint` CI job had never reached it and it had gone unlinted since it was
  written. `cargo fmt --check` and `clippy -D warnings` now run in the bindings CI job.

### Changed

- `FRICTION.md`'s header and footer disagreed on how many findings were fixed, and both
  disagreed with the file. Counted: twelve of seventeen. `AGENTS.md` gained a heat-model
  selection table and quotes a CI-run example function rather than a hand-written snippet, with
  the number it prints pinned by an assertion — prose stating a figure that nothing checks is
  how a document goes stale.

## [0.2.0] — 2026-08-08

Breaking, and almost nothing broke: `&str: Into<String>` meant not one existing call site
changed. Everything in it came from the workspace acquiring its first consumer, and then from
two subagents built out of what that consumer taught.

The headline is not the ergonomics. It is that a first-order accuracy defect in the kernel's own
scheduler — in the schedule chosen *for* accuracy — was found by an application comparing a
coupled run against the closed form of its own recursion, having survived every test the library
had while the conservation audit reported clean to 1e-12.

### Fixed — the kernel scheduler

- **`Schedule::Multirate` does not refine a coupled quantity.** `sweep` steps one domain to
  completion before the next, so a quasi-static publisher puts a whole outer step's joules on the
  bus at once and a subcycling consumer takes all of them on its *first* substep. Refining the
  substep therefore does not move the answer: the error is first order in the **outer** step and
  independent of the substep. Measured on a lumped plate under a lamp — 26.2% low at a 300 s
  outer step, 13.8% at 150 s, 7.1% at 75 s, whatever the substep count — and at 300 s the
  schedule chosen *for* accuracy is worse than `Staggered`, with the errors on opposite sides.

  Every one of those runs passes the conservation audit at ~1e-12. The total that crossed is
  exactly right; only its distribution in time is wrong, and a `Ledger` has no representation for
  *when*. The time-domain twin of the reason `audit_transfers` became a per-face check.

  **Fixed** by `Exchange::take_share(channel, dt)`, new in the kernel: `advance` tells the bus
  what interval the sweep covers and a subcycling consumer asks for its substep's share instead
  of the lot. `Bar1D` and `LumpedMass` use it. Error at a 300 s outer step went from 1.89 K to
  0.304 K — from the worse of the two schedules to fourteen times better than the alternative.

  The share is apportioned against the time *remaining* rather than the whole interval, which is
  what leaves the channel exactly empty: `A·dt/T` with both reduced keeps `A/T`, so the last
  substep takes the remainder. Against the whole interval, `n` shares strand `O(n·ε·A)` where
  `audit_transfers` uses an absolute tolerance.

  Built rather than deferred because the recommendation to wait for a second consumer did not
  survive checking: `04-heater-and-bar` already pairs a quasi-static heater with a bar that
  subcycles hard, so a shipped scene had the defect.

### Fixed — what the two new subagents found

- **`Violation::at`'s cases printed an ungrammatical sentence.** They carry a *message* in
  `quantity` and `Display` had no branch for them, so the first error a consumer ever saw read
  "substance has no heat capacity is not conserved at plate: inf". A third branch.
- **`Report` could not be named without a module path**, though it is what `advance` returns.
  Added to `dualis-core`'s root and the prelude.
- **`Substance` was in the prelude and unbuildable from it**: `bulk` leaves `thermal: None`,
  which `LumpedMass` refuses to step, and the three types needed to supply one were not exported.
  `ThermalProps`, `MechanicalProps`, `AcousticProps`, `ThermalConductivity` and
  `ThermalExpansion` now are.
- **The scene format discarded unknown keys.** `serde` does that by default, which is right for
  a wire protocol and wrong for a saved document: `main.rs`'s own built-in scene kept the
  pre-`release` spelling for two commits and nothing failed, because the keys were dropped and
  the field fell back to its `Default`. Editing them was a no-op that reported success.
  `deny_unknown_fields` on `Scene`, `DomainSpec`, `Release`, `ScheduleSpec` and `Boundary`.
- **The round-trip test could not have caught that.** Both sides of its byte assertion were
  serialiser output, so the hand-written spelling never entered any comparison. It now parses the
  text a person would type and requires each stated value to survive.
- **An unrecognised `finish` produced a silent zero-watt lamp.** The early return also skipped
  `with_reserve`, so the reserve stayed infinite, so `Light::ledger` reported nothing, so the
  audit had nothing to compare — and the scene ran green at `conservation_tolerance(0.0)`, the
  strictest setting expressible, with the lamp doing nothing. One character, `aluminium` against
  `aluminum`. `DomainSpec::build` is fallible now and names the finishes it knows.
- **Two domains could share a name.** `Simulation::domain` takes the first match, so the second
  was never sampled and the first was drawn twice under the second's label and geometry — a
  500 °C bar reported as 20 °C, twice. `World::build` refuses it.
- **A scene whose every domain lacked a field wrote a zero-byte SVG and exited 0**, and "0 KiB"
  could not distinguish that from a legitimate 937-byte strip. The report is a row per *domain*
  now, naming the ones with nothing to draw, with the run-wide extremum beside the final value —
  because a ball that bounced half a metre and one that never moved both end at zero. An empty
  picture is refused rather than written.
- **One colour scale spanned panels of different units**, so a 1 Pa room beside a 7546 m/s orbit
  rendered as an empty bordered square while the numbers beside it looked fine. One extent per
  panel now, still shared across frames.

### Changed — breaking

The first consumer went looking for the API's shape and found the same decision three times:
`&'static str` and `impl Domain` are free when every name is a literal in a test and every
domain type is known when you compile. They are not free for an application. All three are
reversed, and the cost was one order of magnitude smaller than the argument for keeping them —
**no existing call site changed**, because `&str: Into<String>`.

- **`Domain::name` returns `&str`.** Domains store a `String`; constructors take
  `impl Into<String>`, so `Bar1D::new("bar", ..)` and a name read out of a scene file both
  work. `Interface` followed, and `Exchange`'s spatial channel map is keyed by an owned
  interface name.
- **`Report::substeps` is `Vec<(String, u32)>`.** The only breakage across 349 tests: five
  comparisons against string literals. A report outlives the borrow it would otherwise hold.
- **`Simulation::with_boxed(Box<dyn Domain>)`**, plus `impl Domain for Box<dyn Domain>`. The
  simulation always stored boxes; now a caller who built one at run time can hand it over.
- **`Domain::as_field() -> Option<&dyn ScalarField>`** and `Simulation::field(name)`, both
  opt-in and `None` by default, in the style of `as_any`. `ScalarField` was written as the
  interface a visualiser reads a simulation through and was unreachable from `&dyn Domain`, so
  the workspace's own visualiser downcast to concrete types instead — precisely what the
  interface existed to avoid. It no longer names a single domain type.
- **`Room` is in the prelude**, which it should always have been; `Tube` already was.
- **`NBody`, `TreeNBody`, `RigidBody` and `Rolling` answer `as_any`.** All four returned the
  default `None`, so `Simulation::domain_as` could not reach any of them and a renderer got
  nothing back. The failure was a picture with no bodies in it — not an error, not a violation,
  just an empty frame, which is the least debuggable outcome there is. Optics, thermal,
  acoustic and molecular had all opted in because tests inside the workspace had reached for
  them; mechanics had not, because none had.
- **`Bar1D::exposing` takes `impl Into<String>`.** The sweep behind `Domain::name` had matched
  on the parameter being called `name` and missed this one, so a boundary name — the thing two
  domains have to agree on, and therefore exactly the kind of thing that comes from a file —
  was still compile-time only. The kernel's own worked examples were also still teaching
  `fn name(&self) -> &'static str`.

### Added

- **`dualis-world`** — the first consumer, and not published. Scenes described as JSON, built
  into a `Simulation`, run, and drawn as an SVG filmstrip with no dependency. It exists to use
  the SDK from outside rather than to be a good application, and it reports what that was like
  in `crates/dualis-world/FRICTION.md`: twelve findings, seven fixed. Five are recorded rather
  than actioned, and the reasons differ: sharing the examples' SVG plotting would mean
  committing to a public drawing API in a workspace whose scope excludes rendering; validating
  a scene's schedule against its domains at build time is something an application can already
  do with `Domain::max_stable_dt`; and a duplicated face count needs no format change because
  `Exchange::publish_on` already refuses the mismatch by name. The report also says what it
  does *not* cover.

  **Ten scenes ship**, in `crates/dualis-world/scenes/`, covering all five domains,
  and CI runs every one through the real binary as well as the test harness. Three acoustic —
  two room modes and a clap that reflects off all four walls; two thermal, which are the same
  heat told over a plain channel and over a shared boundary and are the argument for
  `Interface` in one picture; two mechanical — four satellites and a bouncing ball whose
  dashpot heat a thermal lump takes; and two molecular, the same 108 atoms at `T* = 0.15` and
  `T* = 1.4`, which side by side are melting. Each has one number asserted, chosen to be a
  property of the physics rather than of the file, and a scene that ships without a claim
  fails the test rather than passing quietly.

  A room can be released as a Gaussian pulse now and not only as a mode, and a panel can hold
  bodies as well as a sampled field — an orbit is a countable number of things at places, and
  rasterising one would invent a continuum it does not have.

  **Bodies are drawn in three dimensions.** The physics always had them: `NBody`,
  `ContactSystem` and `Fluid` all carry `DVec3`, and flattening to a plane was the renderer's
  simplification. They are projected axonometrically now, sorted back to front, with radius
  growing toward the viewer and colour mixed toward the plate for distance — all three are
  needed or the picture is flat however true the coordinates are. Periodic cells get a
  wireframe, because they are a real boundary. The orbit scene tilts its satellites out of one
  plane, which is what makes the third axis carry anything.

  **Optics has a scene**, and it is the fifth `Domain` written outside the library: a
  blackbody lamp on an aluminium mirror whose reflectance falls off in the blue, so the
  colour temperature decides how much of a hundred watts becomes heat. Asserted as the
  difference between 2800 K and 6500 K rather than as one number, because a flat reflectance
  would make `Spectrum` and `SurfaceOptics::absorptance` an expensive way to multiply by a
  constant.

  The scene format couples both ways. On a plain channel, a heater defined **in the
  application** publishes joules and a bar takes them. Over a shared boundary, a beam
  publishes a Gaussian `Flux` onto an `Interface` the bar exposes, and the bar ends up hotter
  in the middle than at the ends by a ratio the scene never states. Both are audited at 1e-9,
  and the totals are checked against `Q / ρVc_p` computed outside the library.

  That closes the gap the report itself had flagged — until then no `publish`, `take`,
  `Exchange`, `Interface` or `Flux` appeared anywhere in the consumer, so the part of the API
  this workspace exists for had only ever been driven from inside. Writing both domains from
  outside needed nothing beyond `dualis::prelude`.
  Excluded from the wasm, determinism and 1.78 jobs, which are promises the *library* makes to
  the people who depend on it.

### Found and fixed

- **A first-order startup error in `Room`, and in `Tube`.** A mode released from rest follows
  `|cos(2πft)|`, but the gap converged at first order against grid resolution where the
  scheme's interior is second. `released_from` left the velocity at `t = 0`; a staggered
  leapfrog carries it at `t = −h/2`, so the first velocity update travelled a whole step where
  it was owed half. `O(h)`, permanent, and `h` follows `dx` through the CFL condition.

  Fixed: the first velocity update takes half a step. Second order now, and the worst
  departure over 20 ms fell from 0.0528 to 0.00238 at 31 cells. Found by the workspace's own
  application checking itself against a closed form, because nothing inside the library was
  checking a rate.

  Two things came out of the fix. **A test had turned the bug into the specification** — one
  step from rest was asserted to move the pressure by `h²c²∇²p`, where Taylor gives `½h²c²∇²p`
  since `ṗ(0) = 0`; the test was missing the half because the scheme was, and `Tube` had the
  same pair. And **the old startup conserved energy exactly while the correct one does not**:
  with `v = 0` read as the half-step value, `Σ∇·(p∇p) = 0` at a rigid wall makes the first
  step's energy change cancel to the last bit. Starting correctly breaks that by `O(h²)` at
  the first step only — 0.42% at 31 cells, quartering on refinement, and 1e-15 thereafter. The
  old code had bought exact bookkeeping by making the scheme first order, which is this
  workspace's own documented trap appearing a second time in the same crate. Energy is now
  reported against the released state as its datum, with the difference available from
  `Room::startup_adjustment` and bounded at 25% so a real first-step bug cannot hide there.

- **The API is comfortable only when the set of domains is known at compile time.**
  `Simulation::with` takes `impl Domain` and there is no `impl Domain for Box<dyn Domain>`;
  domain names are `&'static str`, so a name from a file has to be leaked; and a renderer
  cannot get a `&dyn ScalarField` from a `&dyn Domain`, so it downcasts and knows every domain
  by name — which is what `ScalarField` existed to avoid. Three symptoms of one position. It
  may be the right position, but nobody chose it.

## [0.1.0] — 2026-08-07

First release. All eight crates published to crates.io together.

Everything below was in this release. Two notes on the publish itself, since they cost time
and are not obvious from the outside:

- crates.io requires a **verified** email address, not merely a registered one, and reports
  its absence as a `400` at the first upload rather than at login.
- New crates are rate limited to a burst of five, then roughly one every ten minutes. A
  workspace of eight publishes five, stops, and has to be resumed — so `cargo publish
  --workspace` is not atomic and a partial publish is the normal outcome, not a fault.

### Added

- **`dualis-units`** — dimensional analysis with the SI exponents in the type, so `Length +
  Time` does not compile. Const-generic `Qty` and `QVec3`, macro-generated products, and
  unit-bearing constructors as the only place a factor of a thousand may appear.
- **`dualis-core`** — the kernel. Conservation as an audit (`Ledger`, `audit`, `Violation`),
  fixed-step integrators including `velocity_verlet`, multi-domain scheduling with
  quasi-static, multirate and iterative coupling, deterministic sampling through
  `Rng::for_index`, an accurate discrete Fourier transform, closed-form rigid motion, scalar
  and vector fields, and shared boundaries (`Interface`, `Flux`).
- **`dualis-optics`** — spectral radiometry, surface optics with Fresnel and coatings,
  Sellmeier and Cauchy dispersion, ray geometry, Airy diffraction and the ideal MTF, Zernike
  wavefronts and aberrated PSFs, angular-spectrum propagation, partial coherence, and a
  detector with the four noises that come with counting photons.
- **`dualis-thermal`** — lumped masses and explicit one-dimensional conduction, with radiative
  and convective loss.
- **`dualis-mechanics`** — exact N-body, Barnes-Hut with a quadrupole, penalty contact with
  Coulomb friction, rolling, and rigid-body rotation with Euler's equations.
- **`dualis-acoustic`** — the linear wave equation on a staggered grid, in a tube and in a
  room, with impedance boundaries.
- **`dualis-molecular`** — Lennard-Jones fluids in periodic boxes, cell lists, a Langevin
  thermostat, virial pressure, and radial distribution functions.
- **Five examples**, each of which asserts its numbers and is run by CI. Give one a path and
  it writes an SVG; the plotting has no dependency. Two further examples are checks rather
  than showcases: `agents_quickstart` and `readme_check`.
- **`AGENTS.md` and `CLAUDE.md`** — the API on one page for a consumer, and the gate and
  conventions for a contributor. Written after an AI agent looked for dualis on `PATH`, as a
  Python package, and in a consuming repository, found it in none of the three, and used
  MuJoCo instead. That is a distribution failure and not a documentation one, but the
  quickstart it now lands on is `examples/agents_quickstart.rs`, which CI runs, so it cannot
  drift from the library the way a hand-written snippet does.
- **CI** across Linux, macOS, Windows, two WebAssembly targets and Rust 1.78, with formatting,
  clippy at `-D warnings`, rustdoc at `-D warnings`, licence and advisory checks, and the
  examples.

### Fixed

- **The acoustic wall weighting.** A pressure sample sits *on* a wall, so it owns half a cell,
  and both `Tube` and `Room` divided its divergence by the whole `dx`. Every mode read low —
  1.4% on an 89-cell room — and the scheme converged at first order despite a second-order
  interior. Found by the *rate* of convergence rather than by its size.
- **`Tube`'s absorbing ends** consequently needed their own step limit, `Z·dx/2ρc²`. At the
  full CFL limit the corrected boundary inverts a wave instead of absorbing it: stable, silent
  and wrong.
- **`Bar1D`'s enthalpy reference.** Measured from absolute zero, the bar in these tests holds
  1.42 kJ, so a millijoule arriving is a change in the seventh significant figure and the
  audit's relative check asked for precision the arithmetic had thrown away — a floor that
  grew from 1.6e-12 J at 41 cells to 7.3e-12 J at 161. Measured from the initial temperature,
  the number being summed *is* the change.
- **A statistical test that passed on one seed.** The ideal-gas check asserted that doubling
  the density doubles the departure from `PV = Nk_BT`; across four seeds the ratio came out
  1.35, 1.61, 2.34 and 2.92, averaging to 2.06. It now averages over seeds, because the fix for
  a noisy statistical test is more samples and not a wider tolerance.
- **Licence texts** now ship inside each crate. The packages declared `MIT OR Apache-2.0` and
  contained neither, which `cargo package` does not warn about.
- **Five tolerances that were loose rather than wrong.** A tail-correction ratio checked
  against the `rc⁻³` power law inside 0.02, where the exact ratio differs from that law by
  0.0188 — a known discrepancy filling 94% of the budget. It is now checked against the closed
  form to 1e-12, with the power law asserted separately as the limit it actually is. The
  Langevin settling test averaged one seed and allowed 5%; the seed-to-seed spread is 0.96%,
  so it now averages four seeds and allows 2% — more samples, and a *tighter* tolerance. Two
  gradients expected to be zero were checked with a relative comparison against zero. A
  crystal's negative pressure was asserted as `pressure.max(0.0)`, which an identically zero
  virial satisfies; it is now a band that excludes zero. And the bouncing ball's energy
  handover compared the heat that crossed against `start_energy - 0.0f64.max(0.0)` — the final
  energy was never computed, so the equality its comment promised was a one-sided bound.
  `ContactSystem` now answers `as_any` so the equality can be checked; publishing 5% extra
  heat fails it.

### Documented

- Every public item, with `#![deny(missing_docs)]` in all eight crates so it stays that way.
  244 items had no doc comment; the concentration was in `dualis-units`, which is the API
  nobody can avoid.
- A `compile_fail` doctest proving `Length + Time` does not build — the workspace's reason for
  existing, previously asserted only in prose.

[Unreleased]: https://github.com/YounghyeonPark/dualis/compare/v0.11.0...HEAD
[0.11.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.11.0
[0.10.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.10.0
[0.9.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.9.0
[0.8.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.8.0
[0.7.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.7.0
[0.6.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.6.0
[0.5.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.5.0
[0.4.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.4.0
[0.3.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.3.0
[0.2.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.2.0
[0.1.0]: https://github.com/YounghyeonPark/dualis/releases/tag/v0.1.0
