# dualis for Python

Coupled physics whose conservation audit is an exception you can catch.

```sh
pip install dualis        # wheels for linux, macos and windows; abi3, so 3.10 upward
```

```python
import dualis

sim = dualis.Simulation(schedule="multirate", conservation_tolerance=1e-9)
sim.add_heater("element", watts=2.0, reserve_j=6.0)
sim.add_bar("bar", length_m=0.020, cells=41, area_m2=1e-4)

for _ in range(8):
    sim.advance(0.5)          # raises dualis.Violation if the books do not close

print(sim.temperature("bar"))     # 294.39 K
print(sim.profile("bar")[:3])     # warmer at the fed end
```

## Several bodies and the drop between them

`add_lump` gives you one temperature for a whole assembly. The number that decides whether a
motor survives is the *winding*, which is hotter than the case by however much the joint
between them resists — so for that, a network:

```python
sim = dualis.Simulation(schedule="staggered")
sim.add_heater("losses", watts=6.0, reserve_j=5_400.0)
sim.add_network("motor",
    nodes=[{"name": "winding", "material": "copper", "volume_m3": 18e-6,
            "thickness_m": 2e-3, "initial_k": 298.15},
           {"name": "case", "material": "aluminium", "volume_m3": 220e-6,
            "thickness_m": 4e-3, "initial_k": 298.15,
            "ambient_k": 298.15, "area_m2": 0.042}],   # both keys, or neither
    links=[{"from": "winding", "to": "case", "w_per_k": 0.9}],
    absorbing="winding")

for _ in range(900):
    sim.advance(1.0)

dict(sim.node_temperatures("motor"))     # {'winding': 311.2, 'case': 304.9}
sim.heat_flow_w("motor", "winding", "case")   # 5.585 W crossing the joint
```

Materials: `copper`, `aluminium`, `electrical_steel`, `fr4`, `pla`. A node with `ambient_k` and
`area_m2` loses heat to still air; a node with neither is interior. Giving one without the other
is refused rather than quietly producing a node that looks like it is cooling and is not.

Names are resolved once, when the network is built, so a link to a node that does not exist
raises before any stepping happens and lists the nodes that do. That matters more here than it
looks: a link adds `+q` to one node and `−q` to another in the same sum, so a missing or
misdirected link passes the conservation audit at machine precision and simply reports the wrong
temperature forever.

`temperature()` refuses a network rather than averaging it — the reason to build one is that its
nodes differ, so a mean would describe no part of it.

## Why this rather than numpy

Write physics in numpy and a wrong model runs happily, producing plausible output. There is no
signal; you find out when somebody notices the answer is silly.

Every `advance` here either moves the clock or raises, and the exception is addressable rather
than a sentence to parse:

```python
try:
    sim.advance(0.5)
except dualis.Violation as v:
    v.quantity     # 'energy'
    v.site         # 'bus (published but not consumed)'
    v.before, v.after, v.scale, v.tolerance
    sim.time       # unchanged — a refused step does not advance the clock
```

**Be clear about what it catches.** Quantities appearing or vanishing, amounts left unclaimed on
the bus, and a domain refusing a step it cannot take. It does *not* catch a model that is
internally consistent and physically wrong. For that, check against something the library did not
compute — which is what `aluminium_heat_capacity_j_per_k` is exported for, and what
`tests/test_dualis.py` does throughout.

## What crosses the boundary, and what does not

Values cross as **SI floats with the unit in the parameter name**: `length_m`, `watts`,
`initial_k`. dualis's dimensional types are a compile-time thing — `Length + Time` does not
compile in Rust — and Python cannot have that. A runtime `Quantity` wrapper would cost an
allocation and a check per operation to catch, at run time, a class of error a Python caller
mostly does not make, because it passes values in and reads them out rather than doing algebra on
them. So the types stay on the Rust side and the units stay in the names.

What crosses instead is the audit, which is the part worth having.

## What this cannot do

**You cannot run an ensemble from Python, and it is not an oversight.** `Ensemble` exists to run
many independent samples across threads, and a Python sample function cannot: the GIL serialises
it. A binding for it would look parallel, measure slower than a plain Python loop because of the
crossing cost, and have nothing left of the one thing it is for. Write the sampler in Rust, or
loop in Python and accept one core.

**You cannot write a domain in Python.** That means calling back into the interpreter from inside
the step loop, holding the GIL across it, and deciding what an exception raised mid-sweep does to
a half-advanced simulation. All three are answerable and none cheaply, so this binds the physics
the library has and says so rather than shipping a subtly wrong version of that.

So: enough to **run and audit** coupled physics from Python, not enough to **extend** it. The
extending lives in Rust, and [the workspace](https://github.com/YounghyeonPark/dualis) is where.

A `Simulation` is also bound to the thread that created it. `Domain` is not `Send`, and the
binding says so rather than asking the kernel for a bound it does not offer.

**There is no scene or view binding.** The Rust side has `dualis-scene` for capturing a run and
`dualis-view` for drawing one; neither crosses. Nothing forbids it in principle — the reason is
that a Python caller already has matplotlib, pandas and every plotting stack there is, and what
they want from this side is the *numbers*: `profile`, `node_temperatures`, `temperature`,
`ledger`. Handing them an SVG instead would be the worse half of both worlds. If you want the
HTML report, run the scene from Rust.

## Building

Its own cargo workspace, deliberately. This one resolves fifteen external crates where the
library resolves twelve, seven of them pyo3's and appearing nowhere in the library, plus a
libpython link — and the library's lockfile, licence allow-list and WebAssembly jobs are promises
that should not have to accommodate a Python extension.

```sh
cd bindings/python
pip install "maturin>=1.7,<2.0"
maturin build --release          # an abi3 wheel: one build serves 3.10 upward
pip install target/wheels/dualis-*.whl
python tests/test_dualis.py
```

`maturin develop` works too, in a virtualenv. CI builds the wheel, installs it and runs those
tests on every commit — because "it compiles" proves neither half of "pip install then import".
