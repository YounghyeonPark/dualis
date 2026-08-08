# Working on dualis

For *using* the library, read [AGENTS.md](AGENTS.md) — this file is about changing it.

## The gate, before any commit

```sh
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --no-deps
cargo deny check
for e in beam_hot_spot airy_pattern detector_snr room_modes melting agents_quickstart readme_check; do
  cargo run --locked --release --example "$e" || break
done
# MSRV, library only: `--exclude` needs `--workspace` beside it, and the app is not
# held to 1.78 because nothing depends on it.
cargo +1.78 build --locked --workspace --exclude dualis-world
```

CI also builds `wasm32-unknown-unknown` and runs the suite under `wasm32-wasip1`. All of it is
in [CONTRIBUTING.md](CONTRIBUTING.md), which is the authority; this is the short form.

## The five conventions worth knowing before you start

Stated in full in CONTRIBUTING.md. The compressed version:

1. **Check against a closed form, never against another implementation.** An exact limit, a
   conservation law, or a convergence rate. Comparing to a second copy of the same idea
   verifies nothing.
2. **A tolerance has to be earned.** It should trace to an effect — an integrator's order,
   `1/√N`, a discretisation. A tolerance that was loosened in the same change that made a test
   pass is the thing to look for.
3. **Nothing is random and nothing consults a clock.** `Rng::for_index(seed, index)`. Results
   are bit-for-bit identical across platforms, optimisation levels, WebAssembly and thread
   counts, and there is a pinned digest that says so.
4. **The kernel must never depend on a domain**, and no domain may depend on another. That
   claim is the reason for the crate split and has now been held through five domains.
5. **Every public item is documented.** `#![deny(missing_docs)]` in all eight crates.

## The subagent team

`.claude/agents/` holds five reviewers, each built from a defect this repository actually
shipped rather than from a generic role. They are for developing dualis and are useless to a
consumer.

| agent | asks |
| --- | --- |
| `physics-checker` | Is this claim true? Finds an independent check for it |
| `numerics-reviewer` | Would this test *notice* if the code were wrong? |
| `invariant-guard` | Kernel purity, cross-domain deps, determinism, docs, licence, MSRV |
| `domain-builder` | Scaffolds a new physics crate on the kernel |
| `prose-auditor` | Do the numbers in the README and doc comments still match the code? |

Run `numerics-reviewer` and `physics-checker` on anything touching physics or tolerances, and
`invariant-guard` before a commit that adds a crate, a dependency, or anything with randomness
in it. `prose-auditor` before a release — this repository has shipped stale counts more than
once, most recently claiming six examples when the table listed five.

**Verify what they report by reproducing it.** Their findings have been wrong in both
directions: a seed-spread reported as 0.66% measured 0.96%, which would have produced a
tolerance that looked earned and was not.

## Commit messages

Say what was **wrong**, not only what changed, and give the measurement rather than the
adjective. Several of the more useful commits here are corrections to a mistaken assumption,
and the reasoning is the part worth keeping. Backticks in a message break `git commit -m` under
some shells — write the message to a file and use `git commit -F`.

## The consumer

`crates/dualis-world` is an application, not a library: `publish = false`, and excluded from
the wasm, determinism and 1.78 jobs, because those are promises the *library* makes to people
who depend on it from crates.io. It is covered by `lint`, `test` and `release`.

Its purpose is to use the SDK the way a stranger would and report back. Read
`crates/dualis-world/FRICTION.md` before changing the public API — it is the only record of
what the library feels like from outside, and three of its six findings are the same
underlying decision: **the API is comfortable when the set of domains is known at compile time
and awkward the moment it is not.**

Finding 6 is a live defect in `Room`, pinned by a test that asserts the wrong behaviour on
purpose. Fixing the startup will fail that test. That is intended; update the table in it.

## What is deliberately not here

No GPU, no implicit solvers, no mesh generation, no unstructured grids, no FEM. Adding physics
means a new crate on the kernel, not a new branch inside an existing one. If a change would
make the kernel know about a domain, stop and reconsider — see `domain-builder`, which is
instructed to stop and report in exactly that case.
