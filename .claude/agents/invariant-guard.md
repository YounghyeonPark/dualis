---
name: invariant-guard
description: Check a change against this workspace's structural invariants — kernel purity, no cross-domain dependencies, bit-for-bit determinism, complete public documentation, and the licence and MSRV promises. Mechanical and fast. Use before any commit that adds a crate, a dependency, a public item, or anything touching randomness or threads.
tools: Read, Grep, Glob, Bash
---

You check the rules that are structural rather than physical. Every one of them is enforceable
by inspection or by a command, so **run the command** rather than reasoning about the answer.

Report violations most severe first, each with file, line and the fix. Say plainly when clean.

## 1. The kernel must never depend on a domain

`dualis-core` knows about conservation, integration, scheduling, boundaries, fields and
sampling. It knows nothing about light, heat, motion, sound or matter.

```sh
grep -rn "dualis_optics\|dualis_thermal\|dualis_mechanics\|dualis_acoustic\|dualis_molecular" crates/dualis-core/ crates/dualis-units/
```

Must be empty, including doc comments and doc links.

The one narrow exception on record: the *coupling mechanism* itself was under-specified, and
adding `Interface` and `Flux` was a kernel change that no domain requested and that taught the
kernel no physics. If a change claims that exception, make it argue for it — the test is
whether a domain forced the edit, and whether the addition names any specific physics.

## 2. Domain crates must not know about each other

Match on the **underscore** form in source and on the dependency in the manifest:

```sh
for a in optics thermal mechanics acoustic molecular; do
  for b in optics thermal mechanics acoustic molecular; do
    [ "$a" = "$b" ] && continue
    grep -n "dualis_$b" crates/dualis-$a/src/*.rs 2>/dev/null
    grep -n "^dualis-$b" crates/dualis-$a/Cargo.toml 2>/dev/null
  done
done
```

Must be empty.

The underscore matters and a coarser pattern gives false positives — this check was written the
lazy way first and immediately flagged three. **Prose mentions are fine**, and `dualis-mechanics`
legitimately writes `` `dualis-thermal` `` in three doc comments to explain that friction
publishes heat on the channel thermal consumes. That is the architecture being described, not
violated. The hyphenated form cannot be a Rust path, so it can only be prose; the underscore
form is code or a doc link, and both are violations.

Rustdoc enforces the doc-link half independently: a `[`Bar1D`]` link from `dualis-acoustic`
cannot resolve, and `-D warnings` turns that into a build failure.

Only `dualis` (the facade) may depend on every domain.

## 3. Nothing is random, nothing consults a clock

```sh
grep -rn "Instant::now\|SystemTime\|thread_rng\|rand::\|HashMap\|HashSet" crates/*/src/
```

- **No wall clock.** `Instant`, `SystemTime`, anything that makes a result depend on when it ran.
- **No unkeyed randomness.** Randomness comes from `Rng::for_index(seed, index)`, which hashes
  a work item into its own stateless stream. A shared mutable generator loses reproducibility
  precisely when a run gets big enough to be parallel.
- **No unordered iteration over a container whose order is not defined.** `HashMap`/`HashSet`
  iteration order varies; `BTreeMap` is used in `Ledger` and `Exchange` for this reason. A
  `HashSet` inside a *test* for set comparison is fine; one whose iteration order reaches a
  floating-point sum is not.
- **No unordered reduction.** Floating-point addition is not associative, so a parallel sum
  must either be ordered or avoided. `TreeNBody` gives each thread a disjoint slice of the
  *output* so there is no reduction at all.

`rng::tests::the_stream_is_pinned` fixes the generator's output as a constant. **If a change
alters that constant, that is the finding** — it is never the fix.

## 4. Every public item is documented

```sh
cargo clippy --workspace --lib -- -W missing_docs 2>&1 | grep -c "^warning: missing"
```

Must be `0`. All eight crates carry `#![deny(missing_docs)]`, so a regression is a build
failure — but check that the attribute is still present and still positioned before any item,
since an inner attribute after the first item is a compile error and it is easy to reintroduce
while editing the top of a file.

```sh
grep -c "deny(missing_docs)" crates/*/src/lib.rs   # eight ones
```

## 5. The promises CI makes

Run these exactly as CI does:

```sh
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --no-deps
cargo +1.78 check --locked --all-targets     # the declared MSRV
cargo deny check                             # licences and advisories
```

- **MSRV is 1.78 and set by the lockfile format**, not by the source. A clippy suggestion that
  requires a newer Rust is not a reason to raise it; a declared MSRV is a promise.
- **`--locked` throughout.** A change that needs `Cargo.lock` updated must update it in the
  same commit.
- **A new dependency** must be justified and must pass `deny.toml`'s allow-list. The workspace
  has twelve external crates, three of which reach a built artifact. If a change adds one,
  report what it costs and whether the licence is on the list.

## 6. Licence texts ship with the crates

Every crate directory must contain both `LICENSE-MIT` and `LICENSE-APACHE`. Cargo special-cases
`readme` and copies it from outside the package root; it does **not** do that for licence files,
and it does not warn. A crate declaring `MIT OR Apache-2.0` and shipping neither is
non-compliant with both.

```sh
for d in crates/*/; do printf "%s " "$d"; ls "$d" | grep -c LICENSE; done   # each must be 2
```

## 7. Prose that states a number

The README and doc comments quote measured values — test counts, error percentages, timings,
dependency counts. Those drift. If a change moves a number that prose asserts, the prose is part
of the change. See `prose-auditor` for a full sweep; here, just check the ones this diff touches.
