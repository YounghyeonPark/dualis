# Contributing

Issues and pull requests are welcome. This file is short on process and long on the two or
three conventions that are unusual enough that a contributor would otherwise discover them by
having a change sent back.

## Run what CI runs

**Run it as a script.** The first line is load-bearing and the last line is how you know:

```sh
set -euo pipefail

cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
cargo test --locked --workspace --release
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --no-deps
cargo deny check                       # needs `cargo install cargo-deny`

# The examples are tests that print.
for e in beam_hot_spot airy_pattern detector_snr room_modes melting lens_spots          heat_in_three_dimensions room_in_three_dimensions busbar_rating optical_bench          espresso_shot portafilter_flow agents_quickstart readme_check; do
  cargo run --locked --release --example "$e"
done

echo "the gate passed"     # and if this line does not appear, it did not
```

### This gate has reported a pass it had not earned, seven times

Every one of them the same mistake — reading the *output* of a check as evidence the check **ran**:

| what was done | what happened |
| --- | --- |
| unchained lines | a failure in the middle scrolled past and the last line's status was read as the gate's |
| `... \|\| break` in the examples loop | one example failed, the loop stopped, and the only symptom was output that was not there |
| `cargo clippy ... ; echo OK` | a version bump had invalidated `Cargo.lock`, so `--locked` refused to start and `OK` printed anyway |
| `cargo clippy ... \| tail -1` under `set -e` | a pipeline's exit status is its **last** command's, and `tail` succeeded. This is what `pipefail` is for |
| `cargo publish` twice per crate, once through `grep` and once to read `$?` | the first call published, the second failed with "already exists", and the release loop stopped on its first crate |
| a script that edited two files, then `cargo fmt --all --check`, then committed | the script wrote the first file, raised on the second file's anchor, and **the commit still ran** — unformatted code and half a changelog, red on `main` |

### `set -euo pipefail` is not the fix it looks like, and that is measured

The sixth entry above happened **with `set -eo pipefail` on the first line**, which is why this is
longer than a one-line warning. Measured in the shell these commands actually run in:

```sh
( set -e; false; echo "reached" )                    # prints "reached", exits 0
python -c "raise SystemExit(1)" && echo "reached"    # exits 1
python -c "raise SystemExit(1)" || exit 1            # exits 1
bash -c 'set -e; false; echo "reached"'              # exits 1
```

`set -e` does not take effect in a pasted block here — not even for a builtin. It works in a fresh
`bash -c`, so the option is right and the **paste** is what defeats it.

So: save the gate as a file and run the file, chain every step with `&&`, or run each check as its own
command and read its exit code. The last is what caught the sixth one after three shell guards had
not. And a doubled command is caught by none of them — running each check **once** is the only fix for
that.

Two more of the same shape, outside the shell.

`gh run watch --exit-status` has returned zero with a job still `queued`, so a run was read as green
before a job had started — ask each job for its own `conclusion`.

And `gh run list --limit 1` immediately after a push returns the run for the **previous** commit,
because the new one has not been created yet. Read once that way, a force-push's run was mistaken for
the run of the commit after it — and that run had five jobs `cancelled`, which is neither a pass nor a
failure and would have been reported as one or the other. Select by `headSha`:

```sh
SHA=$(git rev-parse HEAD)
gh run list --limit 10 --json databaseId,headSha -q ".[] | select(.headSha==\"$SHA\") | .databaseId"
```

`cancelled` is its own outcome and worth naming. `concurrency: cancel-in-progress: true` is on this
workflow, so any push supersedes the run before it — which is right, and means a run found by recency
rather than by SHA is quite likely to be one that was killed.

And a script that edits several files fails **halfway**, not cleanly. The mild version writes nothing
while having already printed success for the edits it thought it made; the version that reached `main`
wrote the first file and raised on the second's anchor — a `130x` where the file says `130×` — leaving
a tree that was neither the old state nor the new one. Check every anchor **before** writing any file.
That costs one pass over the inputs and is the difference between an edit that did not happen and an
edit that half did.

CI additionally builds on Rust 1.78, builds for `wasm32-unknown-unknown`, and runs the whole
suite under `wasm32-wasip1` with wasmtime. Those three catch different things and none of them
is decoration — see below.

## The conventions that are not obvious

### Check against a closed form, never against another implementation

A test that compares a function to a second implementation of the same idea passes when both
are wrong in the same way, which is the usual way to be wrong. So every claim here is checked
against something structurally different: an analytic result, an exactly known limit, a
quantity computed by a different route, or a convergence *rate*.

That last one is worth naming, because it caught a real defect. The acoustic domain read every
room mode 1.4% low, which looks exactly like discretisation error and invites a loose
tolerance. What showed it was not discretisation was that refining the grid **halved** the
error instead of quartering it — a second-order scheme converging at first order means the
boundary condition is wrong. No single run could have said that.

### A tolerance has to be earned

`1e-9` because the arithmetic earns it, or `2e-2` because a penalty contact is a non-smooth
potential and the symplectic bound does not apply. Never a number chosen because it made the
test pass. If a tolerance is loose, the comment says which physical or numerical effect is
using up the budget.

Judge a residual against a **scale**, not against itself. A correct system's net conserved
quantity is often exactly zero, so a relative tolerance on it is meaningless; `Ledger` records
the largest contribution beside the total for this reason, and `Violation` carries it.

### Nothing is random, and nothing consults a clock

Fixed steps, no wall clock, no unordered reduction. Randomness comes from
`Rng::for_index(seed, index)`, which hashes a work item into its own stateless stream — so
results are identical whatever order the work is done in, and a parallel run agrees with a
serial one bit for bit. `rng::tests::the_stream_is_pinned` fixes the generator's output as a
constant, and **changing that constant is never the fix**.

This is also why the `wasm32-wasip1` job exists. It runs the pinned digest and every
closed-form comparison under a different target, so a platform that rounded differently would
fail there rather than quietly producing different physics.

### The kernel must never depend on a domain

`dualis-core` knows about conservation, integration, scheduling and boundaries. It knows
nothing about light, heat, motion, sound, electricity or matter. If a new physics needs the kernel changed,
the kernel was wrong — with one narrow exception, which is that the *coupling mechanism* itself
can turn out to be under-specified. That happened once and the README explains why it was not a
violation of the rule.

Domain crates do not depend on each other either. Rustdoc enforces this in a way worth
knowing: an intra-doc link from `dualis-acoustic` to `Bar1D` does not resolve, because the
dependency is not there, and `-D warnings` turns that into a build failure.

The same rule runs the other way for the two layers above the domains. `dualis-scene` and
`dualis-view` depend on the kernel and on each other, and on **no domain at all** — so a physics
added tomorrow is captured and drawn without either being edited. Each proves that rather than
asserting it: the scene layer's test defines a physics inside the test file, and the view
layer's tests are driven by frames written out by hand. A test that ran a real simulation could
not tell *the right view for this shape of data* apart from *the right view for that domain*.

### Every public item is documented

`#![deny(missing_docs)]` in all sixteen crates. A one-line summary that names the unit is enough
for a constructor; anything with a trap in it should say what the trap is.

### Say what was wrong, not only what changed

Commit messages here record the mistake as well as the fix, because the mistake is usually the
more useful half. Several of them document a wrong assumption of the author's that the tests
caught — which grid a boundary condition belongs to, why a statistical test passed on one seed
and not on three others, why an obvious-looking neighbour shell turned out to discriminate
nothing. If you find something like that, write it down.

## What a good pull request looks like

- One idea. A fix and a refactor in the same change are two changes.
- Tests that would fail without it, checked against something independent.
- No new dependency without saying what it buys. The workspace has twelve, three of which
  reach a *published* artifact — the unpublished application links four more — and `deny.toml`
  gates the licences.
- `cargo fmt` clean and `clippy -D warnings` clean.

## Reporting something

An issue with a failing case is worth more than a description. If it is a physics
disagreement, the most useful form is: what the code produced, what the closed form says, and
how the two diverge as something is refined.

## Licence

By contributing you agree that your work is dual-licensed under MIT and Apache-2.0, as the
rest of the workspace is. The README's licence section states this in the standard wording.
