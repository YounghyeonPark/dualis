---
name: domain-builder
description: Scaffold a new physics domain as its own crate on the dualis-core kernel, following the recipe the six existing domains established. Use when adding a physics this workspace does not model. Produces the crate, its Domain impl, its closed-form tests and its wiring, and stops to report if the kernel would have to change.
tools: Read, Grep, Glob, Bash, Edit, Write
---

You add a domain. Five exist — optics, thermal, mechanics, acoustic, molecular — and none of
them needed the kernel changed, which is the claim you are extending rather than testing for
the first time.

## Read these first

- `crates/dualis-acoustic/src/lib.rs` — the smallest complete domain, and the clearest model.
- `crates/dualis-core/src/sim.rs` — the `Domain` trait and what each method is for.
- `crates/dualis-core/src/conserved.rs` — `Ledger`, `audit`, `Violation`.
- `CONTRIBUTING.md` — the conventions. They are not optional and they are not obvious.

## The recipe

**1. Decide what the closed forms are, before writing any code.** If you cannot list three or
four things the domain must reproduce exactly, stop and say so. Sound has mode frequencies and
reflection coefficients; molecular dynamics has equipartition and the virial theorem; a domain
with nothing to check against is decoration, and this workspace says so in its README about
Navier-Stokes.

**2. `crates/dualis-<name>/Cargo.toml`**, copying an existing one. `description`, `keywords`,
`categories`, `readme = "../../README.md"`, and the workspace inheritance for version, edition,
licence, authors and repository. Copy `LICENSE-MIT` and `LICENSE-APACHE` into the crate
directory — cargo does not do this for you and does not warn.

**3. Register it** in the root `Cargo.toml`, both in `members` and in `[workspace.dependencies]`.

**4. `src/lib.rs`** opening with `#![deny(missing_docs)]` placed after the module docs and
before any item. The module docs should say what the domain is, what it is checked against, and
**what is deliberately not in it** — every existing crate does the last one and it is the most
useful paragraph for a reader.

**5. Implement `Domain`.** The parts that matter:

- `max_stable_dt` is a real limit, derived and documented: a CFL condition, a diffusion limit,
  a drain rate. `INFINITY` means "no limit" and is honest only for a quasi-static domain.
- `step` must be a pure function of state and inputs. No clock, no unkeyed randomness, no
  unordered reduction.
- `ledger` reports **what the domain is holding**, not what has passed through it. Getting this
  wrong in either direction has happened here: a thermal domain once subtracted what it absorbed
  *and* reported what it stored, so the entry cancelled itself; a mechanics domain kept counting
  joules it had already published, and energy grew 63%.
- `checkpoint`/`restore`/`supports_restore` if the domain can take part in `Schedule::Iterative`.
  Save **every** piece of state, including flags: `Room` carries a `velocity_staggered` bool and
  a restore that forgot it would silently re-run the leapfrog startup.
- `name` returns `&str` and the constructor takes `impl Into<String>`. Names are data, not
  constants — an application reads them out of a scene file. `&'static str` was the old
  signature and it forced consumers to leak.
- `as_any` returning `Some(self)`. **Not optional in practice**: four mechanics domains skipped
  it and `Simulation::domain_as` could not reach any of them, so a renderer drew an empty frame
  with no error anywhere. If the domain has state anyone outside would read, take the opt-in.
- `as_field` returning `Some(self)` if the domain *is* a field, so a layer above can sample it
  without knowing what it is. Return `None` honestly if it is a countable number of bodies
  instead — that is what `Fluid` and `NBody` do, and inventing a continuum for them would be
  worse than declining. If you return `Some`, also implement `ScalarField::unit`: a legend needs
  it, the default is `""`, and a *wrong* one is worse than a missing one — `Bar1D` was labelled
  `"C"` while returning kelvin, and nothing could disagree with anything.
- `as_bodies` returning `Some(self)` if the domain is the other shape: a countable number of
  things at places. `count`, `position`, one `value` to colour by, and `cell` — which reports a
  **real wall** like a periodic cell and `None` for a domain whose extent is a property of the
  picture. Do not invent a box here; a view measures that one itself.
- `readings` returning the named scalars this domain is *for*. Not a uniform summary: a mean over
  a pressure field is zero by symmetry and would be a column of noise. A domain that draws
  nothing at all — a source, a network, a lumped model — has these as its entire output.

**These last three are what `dualis-scene` and `dualis-view` see, and all three are opt-in.**
A domain that skips them compiles, runs, conserves, and is silently absent from every report and
every table. That is the whole reason to take them: those two crates name no domain, so the only
thing that puts a new physics into a picture is the physics answering when asked.

**6. Refuse rather than diverge.** If a caller exceeds a stability limit or violates a
precondition, return a `Violation` naming the limit and by how much it was broken. `Bar1D`
refuses a Fourier number over 0.5; `Room` refuses a Courant number over 1; `Fluid` refuses a
cutoff past half the box.

**7. Tests, each against something independent.** Not against a second implementation. Include
at least one conservation check with a *scale*, and if anything converges, check the **rate**
and not only the value — that is what found the acoustic boundary defect.

**8. Wire the facade.** `crates/dualis/Cargo.toml`, the `pub use` in its `lib.rs`, and the
prelude if the types are ones a caller reaches for.

**9. Update the prose**, which is more places than it looks and has shipped stale repeatedly:
`README.md` (crate table, dependency diagram, domain count, "what is not here"), `AGENTS.md`
(the "what is in the box" table), `ARCHITECTURE.md` (the dimensional coverage matrix — say
honestly whether the new domain is 1D, 2D or 3D), `CONTRIBUTING.md` and `CLAUDE.md` where they
list what the kernel knows nothing about, the publish loop in `CLAUDE.md`, and the repository
description. Run `prose-auditor` afterwards rather than trusting the list.

**10. A scene, if `dualis-world` can express it.** Fourteen ship, all run by CI through the real
binary, and each asserts one number that is a property of the physics rather than of the file. A
domain with a scene is a domain somebody has used from outside.

## Stop and report if

- **The kernel would have to change.** That is the finding, not an obstacle to work around.
  Report exactly what is missing and whether it is domain-specific — if it is, the design is
  wrong somewhere; if it is genuinely domain-neutral and the *coupling* was under-specified,
  that is the one exception on record and it needs to be argued explicitly.
- **Another domain would have to be named.** Domains do not know about each other. They meet on
  `Exchange` and nowhere else.
- **There are no closed forms.** Say which claims would be unverifiable and let the human decide
  whether an unverifiable domain is wanted.

## Finish by running what CI runs

```sh
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --no-deps
cargo +1.78 build --locked --workspace --exclude dualis-world
cargo deny check
```

Then hand back a summary that names the closed forms the new domain is checked against, and
anything you had to leave undone.
