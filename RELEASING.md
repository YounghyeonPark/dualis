# Releasing dualis

Read this before a release and not otherwise. It was inside `CLAUDE.md`, which is loaded every
session, and a procedure you follow once per release does not need to be in front of you for the
hundred commits in between.

Fifteen crates are published together and share one version. **A published version is permanent** —
it can be yanked, never replaced — so the cost of a release is fifteen permanent version numbers on
crates.io, one on PyPI, and a prose sweep.

## When

**On new public API that somebody outside would reach for.** A new crate, a new type, a new method
on an existing one. Not on a docs fix, not on a CI change, not to exercise the release pipeline —
batch those and let them ride along with the next real one.

`main` being ahead of the registries is the normal state, and the changelog's `[Unreleased]` section
is where the batch accumulates.

## The seven places a version lives

All of them, or the release is broken in a way only one CI job can see:

| | |
| --- | --- |
| `Cargo.toml` | the workspace version and all fifteen path pins — sixteen occurrences |
| `bindings/python/Cargo.toml` | the crate's own version **and** the exact `dualis` pin |
| `bindings/python/pyproject.toml` | the wheel's version |
| `CLAUDE.md` | the `pip install ... .whl` line in the python gate |
| `AGENTS.md` | `dualis = "0.x"`, which `documented_version.rs` checks |
| `crates/dualis/src/lib.rs` | the same series, in the facade's own docs |
| `.claude/agents/invariant-guard.md` | which version is published against which is in the tree |

Then `cargo update --workspace --offline` in the root **and** in `bindings/python`, because
`--locked` refuses a stale lockfile and the bindings' lock has gone stale unnoticed before — nothing
in that job passes `--locked`, so cargo rewrites it silently on every build and the committed copy
drifts.

**The exact `dualis` pin in `bindings/python/Cargo.toml` is the trap.** Bumping the root workspace
and not that leaves it resolving a version that no longer exists. That failed the 0.9.0 release, and
only the `python bindings` job could have caught it — nothing in the main gate reads that directory.

## The prose sweep, which is the part that gets skipped

A release moves the test count, the crate count, the scene count, the FRICTION totals and the install
line across four documents no compiler reads. Three of those *are* under test —
`documented_version.rs` and `friction_counts.rs` — and the rest have shipped stale more than once.

Count them rather than remembering them:

```sh
ls crates | wc -l                                    # crates
ls crates/dualis/examples/*.rs | grep -vc common     # examples
ls crates/dualis-world/scenes/*.json | wc -l         # scenes
cargo test --locked --workspace --release 2>&1 | grep -E "test result:" \
  | awk -F'[; ]' '{p+=$4} END {print p}'             # tests
```

## Publishing, in order

Each crate must be live on the index before the next one resolves it.

```sh
set -euo pipefail
for c in dualis-units dualis-core dualis-acoustic dualis-mechanics dualis-molecular \
         dualis-optics dualis-thermal dualis-electrical dualis-elastic dualis-em \
         dualis-fluid dualis-porous dualis-scene dualis-view dualis; do
  cargo publish -p "$c" --locked      # once per crate. Twice publishes the first and stops on it
done
git tag -a vX.Y.Z -F message.txt && git push origin vX.Y.Z   # the tag publishes the wheel
```

A **new** crate hits crates.io's new-crate rate limit — a burst of five, then roughly one per ten
minutes. Existing crates do not, so a release that adds no crate goes through in one pass.

Verify by resolving from outside rather than by reading the output: `cargo new` a throwaway,
`cargo add dualis@X.Y.Z`, and call something the release added.

## The wheel

**Never `maturin publish` from a workstation.** A local build produces a wheel for *one* platform,
and uploading only that makes `pip install dualis` fail everywhere else — a failure shaped like the
project not supporting Linux rather than like a release mistake.

`.github/workflows/release-python.yml` builds Linux x86_64 and aarch64, macOS x86_64 and aarch64,
Windows x64 and an sdist, installs and runs the tests on every wheel it can execute, and warns on
the cross-compiled ones rather than skipping them silently. It fires on a `v*` tag, or on a manual
dispatch with the `publish` box ticked.

It uses **PyPI trusted publishing**, so there is no token in the repository. Configured on PyPI as
owner `YounghyeonPark`, repository `dualis`, workflow `release-python.yml`, environment `pypi`.

The publish job's `if` needs `always()` and the two results named. A skipped job propagates
**transitively**, and `wheels`/`sdist` opting out with their own `always()` does not opt out for
anything downstream of them — which cost two runs that built all six artefacts and then skipped the
upload with a condition that was correct.

The sdist is why `bindings/python/Cargo.toml` pins `dualis` with **both** a path and a version. An
sdist is a tarball rooted at that directory, so `../../crates/dualis` points outside it; maturin
vendors the whole crate tree in, and the version is what makes the manifest resolvable. Verified by
building the sdist, installing it into a clean venv with `--no-binary :all:`, and running the test
file against what came out.

### What has actually been exercised

A release pipeline you have not run is a guess.

| path | run | result |
| --- | --- | --- |
| dispatch, `publish=false` | build only | six artefacts, both gates skipped |
| dispatch, `publish=true` | 0.3.0 to PyPI | five wheels and an sdist, installed from PyPI and tested |
| tag, version mismatch | `v9.9.9` | refused at `check-version`; nothing built, nothing uploaded |
| tag, version match | `v0.4.0` to PyPI | all four paths now run. `check-version` passed for the first time; five wheels, an sdist, and `pip install dualis==0.4.0` verified from a clean venv |

## Reading the result

Ask each job for its own `conclusion`, not the run for its roll-up:

```sh
gh api "repos/YounghyeonPark/dualis/actions/runs/<id>/jobs?per_page=50" \
  -q '.jobs[] | "\(.name)\t\(.status)\t\(.conclusion // "-")"'
```

`gh run watch --exit-status` has returned zero with a job still `queued`, and the run reported
`success` while `examples` had not started. Re-run that job and wait for its own `status` to reach
`completed`.
