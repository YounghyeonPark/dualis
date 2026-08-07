# Agents

Five, each built around work that recurs in this workspace and around defects that have
actually occurred in it. They are deliberately not generic: a "code reviewer" would have found
none of the bugs listed below, because every one of them was a check that passed while being
blind to the thing it was supposed to catch.

| Agent | Asks | Reach for it when |
| --- | --- | --- |
| `physics-checker` | Is this claim true? | Adding physics, or a number looks wrong |
| `numerics-reviewer` | Would the test notice if it were false? | Any diff touching tests or tolerances |
| `invariant-guard` | Does this break a structural rule? | Before committing anything |
| `domain-builder` | — | Adding a whole new physics as a crate |
| `prose-auditor` | Do the documents still match the code? | Before a release or a publish |

## The two review agents are not the same question

`physics-checker` verifies a claim against an independent route — a closed form, an exact
limit, a conservation law, a convergence rate. It answers *is this right*.

`numerics-reviewer` looks at the check itself and asks what class of error it is structurally
unable to see. It answers *would we find out*. The distinction matters here because the most
expensive bugs in this repository passed every test they had:

- The acoustic wall weighting made every mode read 1.4% low. The conservation audit passed at
  1e-9, because the energy functional and the update were consistent *with each other* and both
  wrong. Only the order of convergence gave it away.
- A spatial flux moved to the wrong part of a boundary keeps the total exactly right, which is
  why `Exchange::audit_transfers` had to become a per-face check.
- An ideal-gas linearity test passed on one seed and would have failed on half the others.

Run them together on anything substantial. They fan out well in parallel.

## Where the standard is written down

`CONTRIBUTING.md` states the conventions for humans; these agents encode the same ones with the
specific failure modes attached. If the two ever disagree, `CONTRIBUTING.md` is the source and
the agent needs updating.

## Adding to this directory

Only if a task recurs *and* has project-specific judgement in it. A generic agent adds a hop
and no knowledge. If you find yourself explaining the same convention to a subagent twice, that
is the signal — write it here rather than in the prompt.
