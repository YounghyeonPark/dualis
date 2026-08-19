# Releasing dualis

Read this before a release and not otherwise. It was inside `CLAUDE.md`, which is loaded every
session, and a procedure you follow once per release does not need to be in front of you for the
hundred commits in between.

Seventeen crates are published together and share one version. **A published version is permanent** —
it can be yanked, never replaced — so the cost of a release is seventeen permanent version numbers on
crates.io, one on PyPI, and a prose sweep.

## When

**On new public API that somebody outside would reach for.** A new crate, a new type, a new method
on an existing one. Not on a docs fix, not on a CI change, not to exercise the release pipeline —
batch those and let them ride along with the next real one.

`main` being ahead of the registries is the normal state, and the changelog's `[Unreleased]` section
is where the batch accumulates.

## The eight places a version lives

All of them, or the release is broken in a way only one CI job can see:

| | occurrences |
| --- | --- |
| `Cargo.toml` | 18 — the workspace version and all seventeen path pins |
| `bindings/python/Cargo.toml` | 2 — the crate's own version **and** the exact `dualis` pin |
| `bindings/python/pyproject.toml` | 1 — the wheel's version |
| `AGENTS.md` | 1 — `dualis = "0.x"`, which `documented_version.rs` checks |
| `crates/dualis/src/lib.rs` | 1 — the same series, in the facade's own docs |
| `.claude/agents/invariant-guard.md` | 1 — which version is published against which is in the tree |
| `CITATION.cff` | 1 — `version`. Also update `date-released`, which is not a version string and so is not caught by the grep below. The grep returns **2**: the other hit is a comment recording which version's Zenodo deposition failed, and bumping that would erase the history it is there for |
| `.zenodo.json` | 1 — `version`. **The row this table was missing**, and it gained it the way the last one did: the 0.15.0 release bumped the seven above and `citation_is_valid` refused, because it asserts the deposition's version *is* the crate's. A table that has now been wrong three times is a table to count against rather than to read |

Count them rather than trusting this table, because it has already been wrong in both directions. It
lost a row when the docs were split — `CLAUDE.md` carried the `pip install ... .whl` line and the
python gate moved to `bindings/python/README.md`, which installs `dualis-*.whl` by glob and needs no
bump at all — and gained one the same day when `CITATION.cff` arrived. Check each row against the file:

```sh
grep -c '0\.14' Cargo.toml bindings/python/Cargo.toml bindings/python/pyproject.toml \
    AGENTS.md crates/dualis/src/lib.rs .claude/agents/invariant-guard.md CITATION.cff .zenodo.json
```

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
find crates -path '*examples*' -name '*.rs' -not -path '*common*' | wc -l   # examples: 15
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
         dualis-fluid dualis-porous dualis-quantum dualis-shape dualis-scene dualis-view dualis; do
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

## The DOI, which is not turned on yet

`CITATION.cff` makes the repository citable by name and version. A **DOI** makes it citable by a
permanent identifier that resolves after the repository moves or disappears, which is what a reference
list actually wants. It is one switch and nobody has thrown it:

1. Sign in to [zenodo.org](https://zenodo.org) with the GitHub account.
2. Under *GitHub*, find `YounghyeonPark/dualis` and turn the toggle **on**.
3. Then publish a release **through GitHub's Releases page**, not by pushing a bare tag. Zenodo
   listens for the release webhook and a pushed tag alone does not fire it — the twelve tags already
   on the repository will therefore get no DOI, and the first release published after the toggle is
   the first one that does.

Zenodo mints two: a **version DOI** for that release and a **concept DOI** that always resolves to the
newest. Cite the concept DOI in prose and the version DOI when the result depends on which version ran
— which for this library it does, because the numbers in the changelog move.

Once the first one exists, add it to `CITATION.cff` as `doi:` and to the BibTeX block in `README.md`.
Neither can be written before there is a DOI to write, which is why they are not there now.

### Zenodo says one thing and it names no field

Every failed deposition reports exactly this, and nothing else:

```json
{"error_id": "b7cac9bfa9784775aa43f702472fcf73", "errors": "Citation metadata load failed"}
```

**It cannot distinguish a licence Zenodo will not store from a byte its YAML reader mishandles.** So the
statement below that v0.13.0 "died on the licence line" was an *inference* from this same message, not
something Zenodo reported — and it may have been wrong, which would make both of the first two fixes
corrections to something that was not the cause.

With a remote that only says "failed" there is nothing to bisect against, so the rule is: **remove every
plausible load hazard at once rather than one per attempt.** Three attempts have been spent one guess at
a time. `citation_is_valid.rs` now holds all of them — no BOM, no CRLF, no non-ASCII, one licence
identifier, and the same for `.zenodo.json`.

The non-ASCII one is worth naming because it was invisible: the file carried a single em dash in
`message`. Valid YAML, valid UTF-8, and handled correctly by every loader anybody would test with —
which is not the same as knowing Zenodo's handles it.

### What was fixed, and the test that had locked in the wrong shape

**v0.13.0**: `license: MIT OR Apache-2.0`. A valid SPDX *expression*, which CFF's schema does not take.
**v0.14.0**: the two-element *list*. Valid CFF — `cffconvert --validate` says so — and Zenodo rejected it
too, because **Zenodo stores exactly one licence per record** and a list is not one.

Valid CFF and valid Zenodo are different things, and `citation_is_valid.rs` had been asserting the first
while the release needed the second. Worse, it asserted the *list* specifically: it was written to lock
in the fix for 0.13.0 and it locked in the shape of the next failure. A test that encodes what the last
change did rather than what the consumer needs is the pattern this workspace has named twice before.

The fix, and the two things worth copying from it:

- **`.zenodo.json`**, which Zenodo prefers over `CITATION.cff`, written in Zenodo's own vocabulary. Its
  licence identifier is **lowercase** — `apache-2.0` — and that was read off
  `zenodo.org/api/vocabularies/licenses?q=apache` rather than guessed, because guessing it would have
  been a third failure on one field.
- **Both files kept correct**, rather than relying on the precedence. "Prefers" is documented behaviour
  and not a guarantee, and keeping both right costs one test.

`CITATION.cff` now names one identifier and says in `message` that it is naming one of two. The real
licensing is unchanged: `Cargo.toml`'s expression and the two LICENSE files are the authoritative
statement, and a DOI record that can hold one licence should say so rather than imply the project has
one.

### The first deposition failed, and the failure is silent from inside the repository

v0.13.0 went to crates.io and PyPI and got **no DOI**. Zenodo read `CITATION.cff`, rejected it, and
said so only as a red *Failed* on its own web page — nothing in the release, the tag, or CI knew.

One line: `license: MIT OR Apache-2.0`. That is a valid SPDX **expression** and `Cargo.toml` is right
to use it; CFF's schema takes an identifier or a **list** of them and an expression matches neither.

`crates/dualis-world/tests/citation_is_valid.rs` now checks that and the other fields a deposition is
built from, so the next one fails in the gate instead of on a web page. Before a release, also check
the page itself: **zenodo.org → GitHub → the repository** lists every release Zenodo has seen and what
it did with each.

**A failed deposition does not retry.** Delete the GitHub release and create it again on the same tag;
that sends a fresh `release: published` event. The tag stays where it is and nothing on crates.io or
PyPI is touched.

### What can and cannot be checked without a Zenodo login

Tried at 0.14.0, because "check the page yourself" is a poor instruction if something cheaper works.
Nothing cheaper does:

| route | result |
| --- | --- |
| `zenodo.org/api/records?q=…` | **useless.** A *failed* deposition is not a record, so zero hits means "failed" or "still queued" or "the index lags" and there is no way to tell which. Fifteen minutes of polling after 0.14.0 returned nothing |
| `zenodo.org/badge/latestdoi/<repo id>` | **useless.** 404s for a repository that has a DOI as readily as for one that does not |
| `curl` against crates.io | 403 — it rejects the default user agent. `cargo search` works |

So RELEASING.md's original instruction stands and is the only one that does: **look at the page.**

One signal *is* visible from the GitHub side and is worth reading, though it is not proof:

```sh
HOOK=$(gh api repos/YounghyeonPark/dualis/hooks -q '.[] | select(.config.url | contains("zenodo")) | .id')
gh api "repos/YounghyeonPark/dualis/hooks/$HOOK/deliveries?per_page=5"   -q '.[] | "\(.delivered_at)  \(.action)  \(.status) \(.status_code)"'
```

Three events fire per release — `created`, `published`, `released` — and Zenodo dedupes, so 409s are
expected on some of them. At 0.14.0 the `released` delivery returned **202 Accepted**; the two
deliveries recorded for 0.13.0 were both **409**. That is a difference and it is where to look first,
but a 202 means "queued", not "deposited", so it does not close the question either.

**v0.13.0 was left without one, deliberately.** Re-depositing means deleting and recreating a release
that is already public, and the fix was worth more than the DOI for a version already out. So the tag
list has a gap at 0.13.0 and **0.14.0 is the first version with a DOI** — recorded here because
otherwise it looks like a mistake later rather than a decision now.

### Where 0.14.0 was left, and what the next release should do

**0.14.0 has no DOI and the tag was not moved.** Three depositions failed, all reporting the same
fieldless message, and the fourth attempt would have needed force-pushing a published tag on a diagnosis
that is not confirmed. The trade was wrong: the DOI is worth something and it is not worth rewriting a
published ref to chase a guess.

What that costs is nothing, because **the fixes are already on `main`** — one licence identifier,
`.zenodo.json`, plain ASCII, and `citation_is_valid.rs` holding all four hazards. So the next release
carries them without a single extra step, and **0.15.0 is the first version that could get a DOI**. Two
versions now have a gap where one was expected; that is recorded here so it reads as a decision.

**If it fails a fourth time, stop guessing and get a real error.** Zenodo's GitHub integration reports
nothing usable, but its REST API validates a deposition field by field:

```sh
# zenodo.org -> Applications -> Personal access tokens, scope deposit:write
curl -s -X POST "https://zenodo.org/api/deposit/depositions?access_token=$ZENODO_TOKEN"   -H "Content-Type: application/json"   -d "{\"metadata\": $(python -c 'import json,sys; print(json.dumps(json.load(open(".zenodo.json"))))')}"
```

That returns the field and the reason, which is the thing three failed releases never produced. It needs
a token, which is why it has not been run — but it is the next step rather than a fifth guess.

**The order matters for the next release.** Throw the switch *before* tagging, or 0.13.0 is another tag
with no DOI and the first citable version waits for 0.14.0.

## Reading the result

Ask each job for its own `conclusion`, not the run for its roll-up:

```sh
gh api "repos/YounghyeonPark/dualis/actions/runs/<id>/jobs?per_page=50" \
  -q '.jobs[] | "\(.name)\t\(.status)\t\(.conclusion // "-")"'
```

`gh run watch --exit-status` has returned zero with a job still `queued`, and the run reported
`success` while `examples` had not started. Re-run that job and wait for its own `status` to reach
`completed`.

Select the run by `headSha` and never by recency. `gh run list --limit 1` straight after a push returns
the previous commit's run, because the new one does not exist yet — and with
`cancel-in-progress: true` on this workflow the previous one has probably just been **cancelled**,
which is neither a pass nor a failure:

```sh
SHA=$(git rev-parse HEAD)
RUN=$(gh run list --limit 10 --json databaseId,headSha \
  -q ".[] | select(.headSha==\"$SHA\") | .databaseId" | head -1)
```

This matters most at a release, where the tag push and the branch push are seconds apart and each
starts a run.
