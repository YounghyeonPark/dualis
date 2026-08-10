# A GPU accelerator, with the CPU domain as the reference

```sh
cd runtime/gpu
cargo test -- --nocapture
```

`GpuSolid` runs `Solid3D`'s seven-point stencil as a WGSL compute shader and implements `Domain`,
so it drops into a `Simulation` like anything else.

## The one rule

**`Solid3D` is the answer and this is a cache of it.** Where the two disagree the CPU is right.
That is what makes `tests/against_the_cpu.rs` the point of the crate rather than a nicety, and it
is the only arrangement under which the library's promises survive a GPU at all.

## What it costs and what it buys

Measured on this machine, 400 steps, the same deposit in the same cell:

| grid | CPU | GPU | |
| --- | --- | --- | --- |
| 16³ | 0.185 s | 0.053 s | 3.5× |
| 32³ | 1.308 s | 0.054 s | 24× |
| 48³ | 4.692 s | 0.055 s | 85× |
| 64³ | 11.418 s | 0.060 s | **191×** |

The GPU column is flat. At these sizes it is bound by dispatch overhead rather than by the
stencil, so the speedup is really the CPU's `n³` growing away from a constant — which is the shape
that makes a GPU worth having and also the reason 16³ barely pays.

Accuracy, after 60 steps with the field still structured:

| | |
| --- | --- |
| worst relative difference from the reference | `2.7e-7` |
| what single precision predicts over 60 steps | `7.7e-7` |
| conservation drift, CPU (`f64`) | `9.1e-15` |
| conservation drift, GPU (`f32`) | `5.0e-11` |

## Why this is a different computation, not a faster one

WGSL has no `f64`. This is single precision against the domain's double, so the two cannot agree —
the useful question is by how much, and the tests measure it rather than asserting it away.

The consequence is not cosmetic: `Simulation`'s conservation audit defaults to a relative `1e-9`,
and `f32` cannot hold that on a long run. A scene using `GpuSolid` needs
`conservation_tolerance_for(quantity::ENERGY, ..)` loosened, and choosing that number is choosing
what the run is allowed to lose. `GpuSolid` also declines `books_balance` for the same reason.

## The finding: single precision was not the problem

The first version stored absolute kelvin and diverged from the reference by `1.4e-3` after two
hundred steps — a thousand times what accumulation predicts. The cause was not accumulation.

The update is `centre + F·(sum − 6·centre)`. On absolute temperatures near 293 K that `sum` is
about 1759, where `f32`'s resolution is `1.2e-4`, and the difference being extracted from it is of
order `1e-3` K. Subtracting two numbers that agree to five digits **keeps less than one digit of
the answer**, every step, forever.

The buffer holds `T − T₀` now. The same numbers are near 1 K, resolution `1.2e-7`, and the
subtraction keeps about four digits:

| | absolute | deviation | |
| --- | --- | --- | --- |
| divergence, 200 steps | `1.449e-3` | `8.7e-7` | 1660× better |
| conservation drift | `7.4e-7` | `1.2e-10` | 6300× better |

The stencil is linear, so subtracting a constant commutes with it exactly and the fix cost
nothing. Single precision was adequate; spending it on an offset nobody needed was not.

## Determinism, and what is deliberately still on the CPU

Rule 5 of `ARCHITECTURE.md` is that results are bit-for-bit across platforms and thread counts.
A GPU cannot promise that in general, and the parts of this that could break it are handled
separately:

- **The stencil is safe.** Each cell reads six neighbours and writes itself, in no particular
  order. There is no reduction and therefore no ordering to depend on.
- **Reductions are not on the GPU.** A mean summed with atomics depends on the order workgroups
  finish, and floating-point addition is not associative, so the answer would change between runs
  on one machine. `ledger` reads the grid back and sums it on the CPU in index order.

`Ensemble` solved the same problem on threads with fixed-size blocks and the same discipline would
work here. A readback is simpler and correct; a faster deterministic reduction is worth writing
when a grid is large enough to need it.

## Two buffers, not one

The stencil ping-pongs between two buffers. Reading and writing one array means some neighbours
are already updated and some are not — a Gauss-Seidel sweep pretending to be Jacobi, which is a
different scheme with a different stability limit and an update order nobody chose.

## No GPU, no test — and it says so

Every test skips with a printed reason when there is no adapter, which is the usual case on a CI
runner. A software rasteriser would be checking a different implementation than anyone runs, so
the CI job lints and builds and lets the tests skip loudly rather than passing for the wrong
reason.
