# Releasing pantometry

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
| `bindings/python/Cargo.toml` | 2 — the crate's own version **and** the exact `pantometry` pin |
| `bindings/python/pyproject.toml` | 1 — the wheel's version |
| `AGENTS.md` | 1 — `pantometry = "0.x"`, which `documented_version.rs` checks |
| `crates/pantometry/src/lib.rs` | 1 — the same series, in the facade's own docs |
| `.claude/agents/invariant-guard.md` | 1 — which version is published against which is in the tree |
| `CITATION.cff` | 1 — `version`. Also update `date-released`, which is not a version string and so is not caught by the grep below. The grep returns **2**: the other hit is a comment recording which version's Zenodo deposition failed, and bumping that would erase the history it is there for |
| `.zenodo.json` | 1 — `version`. **The row this table was missing**, and it gained it the way the last one did: the 0.15.0 release bumped the seven above and `citation_is_valid` refused, because it asserts the deposition's version *is* the crate's. A table that has now been wrong three times is a table to count against rather than to read |

Count them rather than trusting this table, because it has already been wrong in both directions. It
lost a row when the docs were split — `CLAUDE.md` carried the `pip install ... .whl` line and the
python gate moved to `bindings/python/README.md`, which installs `pantometry-*.whl` by glob and needs no
bump at all — and gained one the same day when `CITATION.cff` arrived. Check each row against the file:

```sh
grep -c '0\.14' Cargo.toml bindings/python/Cargo.toml bindings/python/pyproject.toml \
    AGENTS.md crates/pantometry/src/lib.rs .claude/agents/invariant-guard.md CITATION.cff .zenodo.json
```

Then `cargo update --workspace --offline` in the root **and** in `bindings/python`, because
`--locked` refuses a stale lockfile and the bindings' lock has gone stale unnoticed before — nothing
in that job passes `--locked`, so cargo rewrites it silently on every build and the committed copy
drifts.

**The exact `pantometry` pin in `bindings/python/Cargo.toml` is the trap.** Bumping the root workspace
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
ls crates/pantometry-world/scenes/*.json | wc -l         # scenes
cargo test --locked --workspace --release 2>&1 | grep -E "test result:" \
  | awk -F'[; ]' '{p+=$4} END {print p}'             # tests
```

## Publishing, in order

Each crate must be live on the index before the next one resolves it.

```sh
set -euo pipefail
for c in pantometry-units pantometry-core pantometry-acoustic pantometry-mechanics pantometry-molecular \
         pantometry-optics pantometry-thermal pantometry-electrical pantometry-elastic pantometry-em \
         pantometry-fluid pantometry-porous pantometry-quantum pantometry-shape pantometry-scene pantometry-view pantometry; do
  cargo publish -p "$c" --locked      # once per crate. Twice publishes the first and stops on it
done
git tag -a vX.Y.Z -F message.txt && git push origin vX.Y.Z   # the tag publishes the wheel
```

A **new** crate hits crates.io's new-crate rate limit — a burst of five, then roughly one per ten
minutes. Existing crates do not, so a release that adds no crate goes through in one pass. A release
where *every* crate is new takes two hours and needs the retry loop under
[renaming](#renaming-the-project-which-is-not-a-rename); the loop above stops on the sixth.

Verify by resolving from outside rather than by reading the output: `cargo new` a throwaway,
`cargo add pantometry@X.Y.Z`, and call something the release added.

## The wheel

**Never `maturin publish` from a workstation.** A local build produces a wheel for *one* platform,
and uploading only that makes `pip install pantometry` fail everywhere else — a failure shaped like the
project not supporting Linux rather than like a release mistake.

`.github/workflows/release-python.yml` builds Linux x86_64 and aarch64, macOS x86_64 and aarch64,
Windows x64 and an sdist, installs and runs the tests on every wheel it can execute, and warns on
the cross-compiled ones rather than skipping them silently. It fires on a `v*` tag, or on a manual
dispatch with the `publish` box ticked.

It uses **PyPI trusted publishing**, so there is no token in the repository. Configured on PyPI as
owner `YounghyeonPark`, repository `pantometry`, workflow `release-python.yml`, environment `pypi`.

The publish job's `if` needs `always()` and the two results named. A skipped job propagates
**transitively**, and `wheels`/`sdist` opting out with their own `always()` does not opt out for
anything downstream of them — which cost two runs that built all six artefacts and then skipped the
upload with a condition that was correct.

The sdist is why `bindings/python/Cargo.toml` pins `pantometry` with **both** a path and a version. An
sdist is a tarball rooted at that directory, so `../../crates/pantometry` points outside it; maturin
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
| tag, version match | `v0.4.0` to PyPI | all four paths now run. `check-version` passed for the first time; five wheels, an sdist, and `pip install pantometry==0.4.0` verified from a clean venv |

## Renaming the project, which is not a rename

Done once, at 0.16.0, when `dualis` became `pantometry`. Everything below was measured that day and
none of it is in the sections above, because those describe a new *version* of the same name and
every assumption they make about credentials and configuration is an assumption about the name.

**A published name is permanent.** crates.io and PyPI let you yank a version and never release a
name. So this is not a rename — it is publishing a new project and retiring the old one, and the
seventeen `dualis-*` crates will sit on crates.io at 0.15.0 forever. Decide accordingly: the cost is
paid in full at the first publish under a name, and it doubles rather than transfers.

### Choosing the name, which is the part with no undo

Free on both registries is **necessary and nowhere near sufficient**. Seven candidates were free on
crates.io and PyPI and six of them failed the next check:

| candidate | what it collided with |
| --- | --- |
| `clapeyron` | **Clapeyron.jl**, an established fluid-thermodynamics toolkit with an ACS paper — the same field |
| `equipoise` | a live USPTO mark (Hi-Tech Pharmaceuticals), a steroid brand, and an engineering consultancy |
| `conserva` | Conserva Resources, Inc., an IT consultancy since 1986 |
| `virial` | Virial Ltd., a Russian materials company |
| `adiabat` | Adiabat, LLC and Adiabat Technologies SL |
| `holonomy` | an AI startup, Holonomy Systems (fusion), Holonomy Health, Holonomy Consulting |
| `speculum` | a medical-device company, and the common meaning is a surgical instrument |

So the check is three, in this order, and the third is the one that eliminates: **crates.io → PyPI →
a web search for the name as a company or a brand.** A term specific enough that no company would
use it is what survives; `pantometry` exists only as the title of two books, from 1571 and 1830.

None of that is a trademark clearance. A search finds obvious collisions and proves nothing about
the absence of one, and it reaches Korean marks poorly. Check KIPRIS or USPTO before publishing.

### Rename the GitHub repository *first*

`repository` in `Cargo.toml` becomes permanent registry metadata. Rename on GitHub before
publishing or seventeen crates point at a 404 forever. GitHub redirects the old URL to the new one,
so nothing that already exists breaks — the old clone URL, the old links and the old API path all
resolve after the rename.

### Two credentials stop working, and neither says so until it does

**The crates.io token is scoped to crate names.** The token that published seventeen `dualis-*`
crates an hour earlier answered `403 Forbidden: this token does not have the required permissions`
on `pantometry-units`. A new name needs a token with **`publish-new`** and either no crate-name
restriction or one covering the new pattern — and `yank` too, if the same token is to retire the old
crates. It failed on the first crate and published nothing, which is the good version of this
failure; a token that covered *some* of the new names would have stopped halfway.

**PyPI trusted publishing names the repository.** The claim is `owner / repository / workflow /
environment`, so renaming the repository invalidates it: `invalid-publisher: valid token, but no
corresponding publisher`. And a new project name has no publisher at all, so it needs a **pending
publisher** registered on PyPI *before* the first upload — Your account → Publishing → Add a new
pending publisher, with the PyPI project name, the owner, the new repository name,
`release-python.yml`, and environment `pypi`.

Both failed *after* everything else had gone right. All six wheels and the sdist built, the tag
matched the file, and only the upload job failed — so read the jobs and not the roll-up, and
re-run the failed job once the publisher exists rather than re-tagging.

### Seventeen new crates take two hours

crates.io rate-limits **new** crate names: a burst of five, then roughly one per ten minutes.
Existing crates do not, which is why the ordinary release loop above goes through in one pass and
this one does not. A plain loop stops on the sixth crate.

Retry on a rate limit and stop on everything else — a rate limit is a wait and every other failure
is a fault:

```sh
out=$(cargo publish -p "$c" --locked 2>&1)
if [ $? -ne 0 ]; then
  echo "$out" | grep -qi "too many\|rate limit\|try again" || exit 1   # a fault: stop
  sleep 620 && continue                                                # a wait: retry
fi
```

### The tree, and the things outside it

`git ls-files | xargs grep -l` finds every tracked file; there were 263 of them and 1,971
occurrences. Move the crate directories with `git mv` first so the rename shows as renames in the
history rather than as eighteen deletions.

Then regenerate **three** lockfiles — the root, `bindings/python`, and `runtime/editor` — because
`--locked` refuses a lockfile naming crates that no longer exist.

What the grep does not reach is anything outside the tree. The gate script kept
`--exclude dualis-world`, which fails **only at the last of twenty steps**: twenty-five minutes in,
after everything else has passed.

### What the changelog should say

Say once, at the top, that entries below the renaming version name crates as they are now and they
were published under the old name. The alternative is a changelog in two vocabularies, and the
reader needs the sentence and the fact that the registry holds both.

### Then yank

`cargo yank --version X.Y.Z <crate>` for each of the old ones. A yank stops new dependants and
leaves existing lockfiles resolving, so nobody using the old name breaks today and nobody starts
using it tomorrow. Verify from outside — `cargo add <old-name>@X.Y.Z` in a throwaway crate should be
refused — the same way the new name is verified by resolving and calling it.

### Zenodo

The GitHub integration is per repository, so a rename may need the toggle thrown again under the new
name. As ever, nothing inside the repository can see whether a deposition succeeded.

## The DOI, which works now, and how it was verified

`CITATION.cff` makes the repository citable by name and version. A **DOI** makes it citable by a
permanent identifier that resolves after the repository moves or disappears, which is what a reference
list actually wants.

**Minted for the first time at 0.16.0**, after failing silently at 0.13.0 and 0.14.0:

| | |
| --- | --- |
| concept DOI | `10.5281/zenodo.22024817` — always the newest version. Cite this in prose |
| 0.16.0 | `10.5281/zenodo.22024818` — cite this when the result depends on which version ran, which for this library it does |

Both are in `CITATION.cff` and `README.md`'s BibTeX block. The concept DOI is the `doi:` field,
because that is the one a reader following a reference wants; the version DOI lives on each Zenodo
record and in the BibTeX `doi`.

### What the third attempt did differently, and what it proves

Nothing about the *process*. The licence field had been fixed twice — an SPDX expression at 0.13.0,
then a two-element list at 0.14.0, and Zenodo stores exactly one licence — and
`citation_is_valid.rs` had been extended to refuse every plausible load hazard at once rather than
one guess per attempt. This was the first release published after all of that.

The steps, and they are the steps:

1. Sign in to [zenodo.org](https://zenodo.org) with the GitHub account.
2. Under *GitHub*, find `YounghyeonPark/pantometry` and turn the toggle **on**. **Per repository**, so
   a rename may need it thrown again.
3. Then publish a release **through GitHub's Releases page**, not by pushing a bare tag. Zenodo
   listens for the release webhook and a pushed tag alone does not fire it — the sixteen tags already
   on the repository will therefore get no DOI, and the first release published after the toggle is
   the first one that does. (Twelve when this was written, and nothing counts them: a number in
   prose beside a `git tag --list` that anybody can run is the shape this repository keeps finding
   stale.)

Zenodo mints two: a **version DOI** for that release and a **concept DOI** that always resolves to the
newest. Cite the concept DOI in prose and the version DOI when the result depends on which version ran
— which for this library it does, because the numbers in the changelog move.

Once the first one exists, add it to `CITATION.cff` as `doi:` and to the BibTeX block in `README.md`.
**Both are written now.** They could not be before 0.16.0, because there was no DOI to write — which
is the paragraph above surviving three failed depositions as an explanation for an absence, and is
worth keeping beside the fact that the absence is over.

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

`crates/pantometry-world/tests/citation_is_valid.rs` now checks that and the other fields a deposition is
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
| `zenodo.org/api/records?q=…` | **one-sided, and that is enough.** A *failed* deposition is not a record, so zero hits means "failed" or "still queued" or "the index lags" and there is no way to tell which — fifteen minutes of polling after 0.14.0 returned nothing. But a **hit is proof**: at 0.16.0 the first query returned the record, its version DOI, its concept DOI and the stored licence, seconds after the webhook. So poll it, and read a negative as *unknown* rather than as failed |
| the webhook's own deliveries | **worth reading, and not proof.** `gh api repos/OWNER/REPO/hooks/<id>/deliveries` shows what Zenodo answered. One GitHub release fires three events — `published`, `created`, `released` — and only the first is acted on: at 0.16.0 it returned **202 Accepted** and the other two returned 409, which is Zenodo refusing duplicates and not an error. A 202 is acceptance, not success; the deposition can still fail after it, which is exactly what 0.13.0 and 0.14.0 did |
| `zenodo.org/badge/latestdoi/<repo id>` | **useless.** 404s for a repository that has a DOI as readily as for one that does not |
| `curl` against crates.io | 403 — it rejects the default user agent. `cargo search` works |

So RELEASING.md's original instruction stands and is the only one that does: **look at the page.**

One signal *is* visible from the GitHub side and is worth reading, though it is not proof:

```sh
HOOK=$(gh api repos/YounghyeonPark/pantometry/hooks -q '.[] | select(.config.url | contains("zenodo")) | .id')
gh api "repos/YounghyeonPark/pantometry/hooks/$HOOK/deliveries?per_page=5"   -q '.[] | "\(.delivered_at)  \(.action)  \(.status) \(.status_code)"'
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

*Written before 0.13.0, and it took until **0.16.0**. The switch was not the only thing wrong: the
licence field failed twice more after it, each time reporting only "Citation metadata load failed".
Kept because the advice is right and the estimate of what it would cost was three releases short.*

## Reading the result

Ask each job for its own `conclusion`, not the run for its roll-up:

```sh
gh api "repos/YounghyeonPark/pantometry/actions/runs/<id>/jobs?per_page=50" \
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
