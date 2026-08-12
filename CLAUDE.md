# Working on dualis

This file is loaded every session, so it holds only what is needed *every* session. Everything read
once — a release, a change to the bindings, a look at the viewer — lives beside its subject, and the
table at the bottom says when to go there.

For *using* the library, read [AGENTS.md](AGENTS.md). For where it is all going, read
[ARCHITECTURE.md](ARCHITECTURE.md): three layers, the state of each, and the rules that make "add a
physics" cost one crate.

## The gate, before any commit

**Save it and run the file, or chain every step with `&&`.** `set -euo pipefail` in a pasted block
does **not** protect here, measured: `( set -e; false; echo reached )` prints `reached` and exits 0. It
works in a fresh `bash -c`, so the option is right and the paste is what defeats it.

The surest thing is one check per command with its exit code read. That is what caught the sixth time
this gate reported a pass it had not earned — three shell guards had not.
[CONTRIBUTING.md](CONTRIBUTING.md#run-what-ci-runs) is the authority and has all six with what each
cost.

```sh
# Correct in a file, inert when pasted -- see above.
set -euo pipefail

cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --no-deps
cargo deny check
for e in beam_hot_spot airy_pattern detector_snr room_modes melting lens_spots \
         heat_in_three_dimensions room_in_three_dimensions busbar_rating optical_bench \
         espresso_shot portafilter_flow agents_quickstart readme_check; do
  cargo run --locked --release --example "$e"
done
# MSRV, library only: `--exclude` needs `--workspace` beside it, and the app is not
# held to 1.78 because nothing depends on it.
cargo +1.78 build --locked --workspace --exclude dualis-world

echo "the gate passed"     # and if this line does not appear, it did not
```

CI also builds `wasm32-unknown-unknown` and runs the suite under `wasm32-wasip1`. It does **not**
cover `bindings/python` from this gate — that has its own job and its own procedure.

**Read a result from the thing that produced it, and read whether the check *ran* rather than what it
printed.** Six times this gate has said `ok` while failing, and once that reached `main`. A CI run's
roll-up has said `success` with a job still `queued`; ask each job for its own `conclusion`. A script
that edits several files can write the first and raise on the next, so check every anchor before
writing any of them.

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
   claim is the reason for the crate split and has now been held through ten domains.
5. **Every public item is documented.** `#![deny(missing_docs)]` in all sixteen crates.

## Commit messages

Say what was **wrong**, not only what changed, and give the measurement rather than the adjective.
Several of the more useful commits here are corrections to a mistaken assumption, and the reasoning
is the part worth keeping. Backticks in a message break `git commit -m` under some shells — write the
message to a file and use `git commit -F`.

## The subagent team

`.claude/agents/` holds seven reviewers, each built from a defect this repository actually shipped
rather than from a generic role. They are for developing dualis and are useless to a consumer.

| agent | asks |
| --- | --- |
| `physics-checker` | Is this claim true? Finds an independent check for it |
| `numerics-reviewer` | Would this test *notice* if the code were wrong? |
| `silent-failure-hunter` | What here can come out *empty* and look fine? |
| `consumer-advocate` | What is this API like from outside? |
| `invariant-guard` | Kernel purity, cross-domain deps, determinism, docs, licence, MSRV |
| `domain-builder` | Scaffolds a new physics crate on the kernel |
| `prose-auditor` | Do the numbers in the README and doc comments still match the code? |

Run `numerics-reviewer` and `physics-checker` on anything touching physics or tolerances, and
`invariant-guard` before a commit that adds a crate, a dependency, or anything with randomness in it.
`prose-auditor` before a release.

**Verify what they report by reproducing it.** Their findings have been wrong in both directions: a
seed-spread reported as 0.66% measured 0.96%, which would have produced a tolerance that looked
earned and was not.

## Where the rest of it is, and when to go there

| read | before |
| --- | --- |
| [RELEASING.md](RELEASING.md) | any release. Cadence, the seven places a version lives, the crate order, the wheel, and what the pipeline has actually been run through |
| [CONTRIBUTING.md](CONTRIBUTING.md) | changing a test or a tolerance. The authority on the gate and on the five conventions in full |
| [crates/dualis-world/FRICTION.md](crates/dualis-world/FRICTION.md) | changing the public API. Twenty-three findings from using the SDK as a stranger, five of them the same underlying decision |
| [bindings/python/README.md](bindings/python/README.md) | touching the bindings. Its own cargo workspace, its own gate, and the two boundary decisions not to relitigate |
| [runtime/viewer/README.md](runtime/viewer/README.md) | touching the viewer. Why it is a separate workspace and why it does not link `dualis` |
| [.claude/agents/README.md](.claude/agents/README.md) | adding a reviewer |

Three of those exist because a dependency tree does not belong in the library's lockfile. Measured:
the library resolves **12** external crates, `bindings/python` **15**, and the viewer's wgpu stack
**86**. `deny.toml` gates every one of the library's twelve, CI builds with `--locked`, and the same
crates go to `wasm32` and Rust 1.78 — none of which can carry a GPU stack or a libpython link.

`crates/dualis-world` is the fourth: an application, `publish = false`, whose purpose is to use the
SDK the way a stranger would and report back. Anything you find yourself adding to it that a
*consumer* would want is in the wrong crate.

## What is deliberately not here

No GPU in the library, no implicit solvers, no mesh generation, no unstructured grids, no FEM beyond
the trilinear element `dualis-elastic` uses. Adding physics means a new crate on the kernel, not a new
branch inside an existing one. If a change would make the kernel know about a domain, stop and
reconsider — see `domain-builder`, which is instructed to stop and report in exactly that case.
