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
   claim is the reason for the crate split and has now been held through six domains.
5. **Every public item is documented.** `#![deny(missing_docs)]` in all ten crates.

## The subagent team

`.claude/agents/` holds seven reviewers, each built from a defect this repository actually
shipped rather than from a generic role. They are for developing dualis and are useless to a
consumer.

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
who depend on it from crates.io. It is covered by `lint`, `test`, `release`, `examples` — which is where every scene is run
through the real binary — and `deny`.

Its purpose is to use the SDK the way a stranger would and report back. Read
`crates/dualis-world/FRICTION.md` before changing the public API — it is the only record of
what the library feels like from outside, and five of its twenty-one findings are the same
underlying decision: **the API is comfortable when the set of parts is known at compile time
and awkward the moment it is not.**

Finding 6 was a live defect in `Room` and in `Tube`, and it is fixed: the first velocity update
of the staggered leapfrog now travels half a step. Two tests in `dualis-acoustic` and one in
`dualis-world` pin the second-order rate. The startup no longer conserves energy exactly at the
first step, by design — `Room::startup_adjustment` reports the `O(h²)` difference, measured
0.42% at 31 cells and quartering on refinement.

## The Python bindings

`bindings/python` is a **separate cargo workspace**, excluded from the root one. That is the
escalation this repository wrote down for itself, and pyo3 is what triggered it: fifteen crates
and a libpython link, against a library that has kept to twelve external dependencies with every
one gated by `deny.toml`. Folding it in would put pyo3's tree in the library's lockfile and make
`cargo build --workspace` need a Python development environment.

The gate above therefore does not touch it. It has its own CI job, which builds a wheel,
**installs it**, and runs `tests/test_dualis.py` — because "the cdylib compiles" proves neither
half of "pip install then import". Its own gate:

```sh
cd bindings/python
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings
python -m maturin build --release
python -m pip install --force-reinstall target/wheels/dualis-0.7.0-*.whl
python tests/test_dualis.py
```

Two decisions worth not relitigating. The boundary is **SI floats with the unit in the parameter
name**, because the dimensional types are compile-time and a runtime wrapper would cost per
operation to catch an error a Python caller does not make. And **a `Domain` cannot be written in
Python**: that needs callbacks into the interpreter from inside the step loop, the GIL held
across it, and an answer for what an exception raised mid-sweep does to a half-advanced
simulation.

### How often to release, and the order

Nine crates are published together and share one version. A published version is permanent — it
can be yanked, never replaced — so the cost of a release is nine permanent version numbers on
crates.io, one on PyPI, and a prose sweep.

**Release on new public API that somebody outside would reach for.** A new crate, a new type, a
new method on an existing one. Not on a docs fix, not on a CI change, not to exercise the release
pipeline — batch those and let them ride along with the next real one. `main` being ahead of the
registries is the normal state, and the changelog's `[Unreleased]` section is where the batch
accumulates.

The order matters: each crate must be live on the index before the next one resolves it.

```sh
for c in dualis-units dualis-core dualis-acoustic dualis-mechanics dualis-molecular          dualis-optics dualis-thermal dualis-electrical dualis; do
  cargo publish -p "$c" --locked || break
done
git tag -a vX.Y.Z -F message.txt && git push origin vX.Y.Z   # the tag publishes the wheel
```

A *new* crate hits crates.io's new-crate rate limit — a burst of five, then roughly one per ten
minutes. Existing crates do not, so a release that adds no crate goes through in one pass.

### Releasing the wheel

**Never `maturin publish` from a workstation.** A local build produces a wheel for *one*
platform, and uploading only that makes `pip install dualis` fail everywhere else — a failure
shaped like the project not supporting Linux rather than like a release mistake.

`.github/workflows/release-python.yml` builds Linux x86_64 and aarch64, macOS x86_64 and
aarch64, Windows x64 and an sdist, installs and runs the tests on every wheel it can execute,
and warns on the cross-compiled ones rather than skipping them silently. It fires on a `v*` tag,
or on a manual dispatch with the `publish` box ticked.

It uses **PyPI trusted publishing**, so there is no token in the repository. Configured on PyPI
as owner `YounghyeonPark`, repository `dualis`, workflow `release-python.yml`, environment
`pypi`.

What has actually been exercised, since a release pipeline you have not run is a guess:

| path | run | result |
| --- | --- | --- |
| dispatch, `publish=false` | build only | six artefacts, both gates skipped |
| dispatch, `publish=true` | 0.3.0 to PyPI | five wheels and an sdist, installed from PyPI and tested |
| tag, version mismatch | `v9.9.9` | refused at `check-version`; nothing built, nothing uploaded |
| tag, version match | `v0.4.0` to PyPI | all four paths now run. `check-version` passed for the first time; five wheels, an sdist, and `pip install dualis==0.4.0` verified from a clean venv |

The publish job's `if` needs `always()` and the two results named. A skipped job propagates
**transitively**, and `wheels`/`sdist` opting out with their own `always()` does not opt out for
anything downstream of them — which cost two runs that built all six artefacts and then skipped
the upload with a condition that was correct.

The sdist is why `bindings/python/Cargo.toml` pins `dualis` with **both** a path and a version.
An sdist is a tarball rooted at that directory, so `../../crates/dualis` points outside it;
maturin vendors the whole crate tree in, and the version is what makes the manifest resolvable.
Verified by building the sdist, installing it into a clean venv with `--no-binary :all:`, and
running the test file against what came out.

## What is deliberately not here

No GPU, no implicit solvers, no mesh generation, no unstructured grids, no FEM. Adding physics
means a new crate on the kernel, not a new branch inside an existing one. If a change would
make the kernel know about a domain, stop and reconsider — see `domain-builder`, which is
instructed to stop and report in exactly that case.
