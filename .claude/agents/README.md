# Agents

Seven, each built around work that recurs in this workspace and around defects that have
actually occurred in it. They are deliberately not generic: a "code reviewer" would have found
none of the bugs listed below, because every one of them was a check that passed while being
blind to the thing it was supposed to catch.

| Agent | Asks | Reach for it when |
| --- | --- | --- |
| `physics-checker` | Is this claim true? | Adding physics, or a number looks wrong |
| `numerics-reviewer` | Would the test notice if it were false? | Any diff touching tests or tolerances |
| `silent-failure-hunter` | What here can come out *empty* and look fine? | Opt-in trait methods, serialised formats, renderers |
| `consumer-advocate` | What is this API like from outside? | Before changing the public API or releasing |
| `invariant-guard` | Does this break a structural rule? | Before committing anything |
| `domain-builder` | — | Adding a whole new physics as a crate |
| `prose-auditor` | Do the documents still match the code? | Before a release or a publish |

## The review agents ask four different questions

They overlap in subject and not in method, and the distinction is the whole point.

`physics-checker` verifies a claim against an independent route — a closed form, an exact limit,
a conservation law, a convergence rate. It answers *is this right*.

`numerics-reviewer` looks at the check itself and asks what class of error it is structurally
unable to see. It answers *would we find out*. The most expensive bugs here passed every test
they had:

- The acoustic wall weighting made every mode read 1.4% low. The conservation audit passed at
  1e-9, because the energy functional and the update were consistent *with each other* and both
  wrong. Only the order of convergence gave it away.
- A spatial flux moved to the wrong part of a boundary keeps the total exactly right, which is
  why `Exchange::audit_transfers` had to become a per-face check.
- An ideal-gas linearity test passed on one seed and would have failed on half the others.

`silent-failure-hunter` asks a third thing: *what comes out absent*. Not a wrong value — no
value, and no complaint either. Four mechanics domains had never taken the `as_any` opt-in, so
the orbit scene parsed, ran, conserved to 1e-6, printed a clean summary and drew an empty frame,
with no error anywhere. A scene helper kept emitting config keys that `serde` silently ignored
and every test passed while the file said something the code no longer read.

`consumer-advocate` asks the question none of the others can, because it requires not knowing
how the library came to be: *what is this like to use*. It is not a review — it builds something
against the public API and reports what building was like. It is also the method with the best
record here by a distance: thirty-four findings, twenty-eight fixed, including the only real physics
defect found in the period. `Room` and `Tube` were starting a staggered leapfrog with the
velocity at the wrong time level, `O(h)` and permanent, and it survived 345 passing tests — two
of which had turned the bug into the specification.

Six came from sources the agent's own instructions did not anticipate. Five from **splitting the
application into layers**: building against an API finds what is awkward, and pulling a layer out
of one finds what was never stated. The sixth from **building the next domain**, which is the
cheapest of the four to run — a layer that names no domain still carries assumptions about every
domain it has met, and `dualis-scene` assumed fields were flat until one with a volume arrived.

Run the first three together on anything substantial; they fan out well in parallel. Run
`consumer-advocate` on its own, because it needs to write code and the others do not.

## Where the standard is written down

`CONTRIBUTING.md` states the conventions for humans; these agents encode the same ones with the
specific failure modes attached. If the two ever disagree, `CONTRIBUTING.md` is the source and
the agent needs updating. `crates/dualis-world/FRICTION.md` is the standing record of what the
API is like from outside, and `consumer-advocate` reads it first so it does not report what is
already known.

**The agents are prose, and prose goes stale.** `domain-builder` was still teaching
`fn name(&self) -> &'static str` after that signature was removed from every crate, and
`invariant-guard` had nothing about the version obligation that publishing to crates.io created.
Run `prose-auditor` over this directory when the public API moves — it is documentation like any
other, and it is documentation a subagent will act on.

## Adding to this directory

Only if a task recurs *and* has project-specific judgement in it. A generic agent adds a hop and
no knowledge. If you find yourself explaining the same convention to a subagent twice, that is
the signal — write it here rather than in the prompt.

And only on evidence. Every agent above points at a defect that happened; a release warden that
diffs the public API against the published version would be useful and is *not* here, because
the version has never yet been got wrong. That check went into `invariant-guard` as one more
structural item instead, which is proportionate to a risk that has not yet materialised.
