# dualis, for an agent

Everything needed to write a working dualis program, on one page. If you are here to
*modify* dualis rather than use it, read [CLAUDE.md](CLAUDE.md) instead.

**dualis is a Rust library**, with Python bindings. It is not a CLI and there is nothing on
`PATH`, so `which dualis` will fail and that is not a broken installation.

From Python, `import dualis` works — see [`bindings/python`](bindings/python), built from this
repository and published as `dualis` on PyPI. It can *run and audit* the library's physics; it cannot *extend*
it, because writing a `Domain` in Python is unsupported and the reasons are written down there.
Everything below describes the Rust API.

```toml
[dependencies]
dualis = "0.9"
```

API docs: <https://docs.rs/dualis>. Source: <https://github.com/YounghyeonPark/dualis>.

---

## Why you might want this specifically

When you write physics in numpy or a general-purpose engine and get it wrong, the code runs
happily and produces plausible output. There is no signal. You find out when a human notices
the answer is silly.

dualis audits conservation on every step and returns a **`Violation` that names what went
missing and where**:

```
energy destroyed at simulation: 5.000000e2 became 4.995000e2,
a relative change of 1.000e-3 against a tolerance of 1.000e-9
```

`advance` returns `Result<Report, Violation>`, and the fields are machine-readable —
`quantity`, `site`, `before`, `after`, `scale`, `tolerance`. That is a correctness signal you
can loop on without asking anyone.

**Be clear about what it does not catch.** It catches quantities appearing or vanishing,
amounts left unclaimed on the bus, and fluxes that disagree face by face across a shared
boundary. It does not catch a model that is internally consistent and physically wrong:
publish a power where a joule was wanted, forgetting the factor of `dt`, and both sides agree
perfectly about a number that is off by `1/dt`. For that class, check against something the
code did not compute — a closed form, an exact limit, or a convergence rate.

---

## The whole thing, in three ideas

### 1. Units are types

Dimensions live in the type, so `Length + Time` does not compile. There is exactly one place
a factor of a thousand may appear — a unit-bearing constructor — and `to_si()` is the only way
back to a bare `f64`.

```rust
use dualis::prelude::*;

let area: Area = Length::mm(10.0) * Length::mm(10.0);
let absorbed: Power = Irradiance::mw_per_cm2(50.0) * area * 0.02;
let capacity: HeatCapacity = Mass::g(2.0) * SpecificHeat::j_per_kg_k(858.0);
let rise: Temperature = (absorbed * Time::s(1.0)) / capacity;
```

Do not reach for `f64` and a comment. If the algebra does not typecheck, the physics is
usually wrong.

### 2. A domain is anything that steps

Two methods are required. Everything else has a default.

```rust
impl Domain for MyThing {
    fn name(&self) -> &str { &self.name }

    fn step(&mut self, t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        // advance your state by dt, publish or take on the bus
        Ok(())
    }
}
```

Worth overriding: `ledger()` (your books — without it the audit has nothing to check),
`max_stable_dt()` (so a scheduler can subcycle you), `kind()` (`Kind::QuasiStatic` if you have
no state to march), `as_field()` (so a renderer can sample you without knowing what you are),
`as_any()` (so callers can downcast back to your concrete type).

Names are data, not constants. Constructors take `impl Into<String>`, so `"my-thing"` and a
`String` read out of a scene file both work, and `Simulation::with_boxed` takes a
`Box<dyn Domain>` for when the type is chosen at run time.

### 3. Domains never call each other

They meet on `Exchange`, a bus of named channels carrying **SI amounts, not rates** — joules,
not watts. A domain steps over an interval, so what crossed is an amount. That is what makes
the audit an equality rather than an approximation.

```rust
bus.publish(HEAT, joules);   // publisher, having multiplied by its own dt
let arrived = bus.take(HEAT); // consumer
```

**A ledger says what you are holding, not what has passed through you.** Joules you published
are gone from your books and are being reported by whoever took them. Adding them back is the
most common way to write a domain that audits green and is wrong.

Anything left unclaimed when the step ends is a `Violation`.

---

## A complete program

[`examples/agents_quickstart.rs`](crates/dualis/examples/agents_quickstart.rs) is a runnable
version of everything above — a publisher, a consumer, the books closing, and then the same
pair with a deliberate 10% leak so you can see the `Violation` it produces. CI runs it, so it
cannot drift.

```sh
cargo run --example agents_quickstart
```

Read that file before writing your own domain. It is about two hundred and fifty lines, a third
of them commentary.

---

## What is in the box

`dualis` is a facade; every name below is re-exported through `dualis::prelude::*`.

| crate | what it holds |
| --- | --- |
| `dualis-units` | `Length`, `Time`, `Mass`, `Energy`, `Power`, … and the vector forms. Dimensions in the type |
| `dualis-core` | The kernel: `Domain`, `Exchange`, `Simulation`, `Schedule`, `Ledger`, `Violation`, `Interface`, `Flux`, `Rng`, `Ensemble`, `Pose`, `Bodies` and `Reading` for what a domain offers a viewer, integrators, fields |
| `dualis-optics` | Radiometry, Fresnel and coatings, dispersion, rays, Airy diffraction, MTF, Zernike, PSFs, detector noise |
| `dualis-thermal` | `LumpedMass`, `Bar1D` conduction, `Solid3D` conduction in three dimensions on a cubic grid, `ThermalNetwork` of n bodies joined by conductances, radiative and convective loss |
| `dualis-mechanics` | `NBody`, `TreeNBody` (Barnes-Hut), `ContactSystem` with friction, `RigidBody` |
| `dualis-acoustic` | The wave equation on a staggered grid: `Tube` (1D), `Room` (2D), `Hall` (3D, with the vertical and oblique modes a floor plan cannot have), impedance boundaries |
| `dualis-molecular` | `Fluid` with Lennard-Jones, `PeriodicBox`, cell lists, Langevin thermostat, `RadialDistribution` |
| `dualis-electrical` | `Winding`: `I²R` onto the heat channel, copper rising 0.393%/K, and `runaway_current` — the exact threshold `√(g/(R₂₀α))` where the feedback overtakes the heat path. `Conductor`: `∇·(σ∇φ)=0` solved on a grid, so a resistance is a property of a *shape* — `ρL/A` exactly for a bar, and whatever a notch gives for a notch |
| `dualis-scene` | One layer up. `Placement`, `Extent`, `capture` — where a domain sits and what one instant of a run looks like, as `Frame`, `Panel`, `PanelData`. Names no domain |
| `dualis-view` | Two layers up. `svg` filmstrip, `html` report that opens in a browser with nothing installed, `readings_csv`, `to_json`. The view is chosen by the shape of the data |

`Schedule` picks how they interact: `OneWay`, `Staggered` (declaration order is execution
order), `Iterative { max_iter, tol }` for strong coupling, `Multirate` for domains with very
different stability limits.

### Picking a heat model

Three shapes, and choosing wrong is the most common way to get a plausible number that answers
the wrong question:

| you want | use | cost |
| --- | --- | --- |
| one body, one temperature | `LumpedMass` | one step per frame; check `biot_number() < 0.1` first |
| a gradient along one thing | `Bar1D` | pays `dt < dx²/2α`; milliseconds in aluminium |
| several bodies and the **drop between them** | `ThermalNetwork` | one step per frame per node |

The third is the one people reach for too late. A motor, a laser diode on a submount, a die in a
package: what fails is the winding or the junction, and what you can measure is the case. A
`LumpedMass` gives you one number for both, and it is the case's.

```rust
let mut motor = ThermalNetwork::new("motor");
let winding = motor.node("winding", Substance::copper(), Volume::from_si(18e-6),
                         Length::mm(2.0), Temperature::celsius(25.0));
let case = motor.node_losing_to("case", Substance::aluminium_6061(), Volume::from_si(220e-6),
                                Length::mm(4.0), Temperature::celsius(25.0),
                                Environment::still_air(Temperature::celsius(25.0),
                                                       Area::from_si(0.042)));
motor.link(winding, case, Conductance::w_per_k(0.9))?;
motor.absorbing(winding)?;                  // where heat off the bus lands

// after stepping — 6 W in, and the joint carries a 6.2 K drop at fifteen minutes:
motor.temperature(winding);                 // the number that decides survival
motor.heat_flow(winding, case);
```

That is `junction_to_case` in `examples/agents_quickstart.rs`, which CI runs — so it is a
snippet that compiled and produced that number, not one written into a document by hand.

If what you want is the *settled* answer rather than the trajectory, do not march to it:

```rust
let settled = motor.steady_state(Power::w(6.0))?;   // solves the balance directly
settled.temperature(winding);
```

Exact rather than "close enough after enough steps", and it refuses a network where heat has
nowhere to go — that has no steady state, and a finite number would be the wrong answer to a
question with none.

Nodes are `Node` handles rather than names on purpose — see below. `node_named` is the way in if
you are building from a file, and `handles()` walks a network you did not build.

---

## Rules you will otherwise break

These are enforced by CI and by the audit, so breaking them fails loudly rather than quietly.
They are listed here so it does not have to be loudly.

- **No wall clock, no unseeded randomness.** Results are bit-for-bit reproducible across
  platforms, optimisation levels, WebAssembly and thread counts. Use `Rng::for_index(seed, i)`,
  which gives the same value for the same index no matter what order you ask in. `Date::now`,
  `rand::thread_rng`, and reductions over unordered collections all break this.
- **Domains do not depend on each other.** If your new physics needs to `use dualis_thermal`,
  the design is wrong — publish on a channel instead. The kernel depends on no domain either.
- **Every public item is documented.** `#![deny(missing_docs)]` is set in all eleven crates.
- **MSRV is 1.78**, checked by CI.
- **Tolerances are earned.** A number in an `assert!` should trace to an effect — an
  integrator's order, `1/√N` for a sample count, a discretisation. If you cannot say which,
  it was chosen to make the test pass. See CONTRIBUTING.md.

---

## Two things that will trip you specifically

**A relative tolerance against a quantity that should be zero is meaningless.** A correct
system's net conserved quantity is usually exactly zero, so `(value - 0.0).abs() / 0.0` is a
100% error for anything nonzero. Supply a scale from outside — the energy that actually
crossed, the largest single contribution. `Ledger` records a `scale` beside each total for
exactly this reason, and `Violation` carries it.

**A conservation check that passes is necessary and nowhere near sufficient.** A flux
redistributed to the wrong part of a boundary keeps the total exactly right. Ask what class of
error your check is blind to.

The clearest case in this library is `ThermalNetwork`. A link adds `+q` to one node and `−q` to
another *in the same sum*, so they cancel identically: a sign error, a transposed index, or a
link you forgot to add passes the conservation audit at machine precision, and the winding just
runs at a plausible wrong temperature forever. That is why nodes are handles you cannot forge
and why every test for it is per node or against a closed form. It is also not hypothetical —
building it surfaced an `O(h)` bias on the joint next to the heat source that the audit never
saw, and a series-resistance formula found in one run.

The habit that generalises: after your audit passes, write down one wrong program that would
also pass it, then go and check *that*.

---

## Where to look next

- **`examples/`** — six worked problems that print their numbers and assert every one of
  them. `cargo run --example melting`. Give any of them a path and it writes an SVG.
- **[ARCHITECTURE.md](ARCHITECTURE.md)** — the three layers, what is built and what is missing,
  and the rules that make "add a physics" cost one crate.
- **[README.md](README.md)** — the long version, including what is deliberately *not* here.
- **[CONTRIBUTING.md](CONTRIBUTING.md)** — the conventions, and the gate CI runs.
- **[CLAUDE.md](CLAUDE.md)** — working on dualis rather than with it.
