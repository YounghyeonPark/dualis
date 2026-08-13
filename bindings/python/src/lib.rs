//! Python bindings for dualis.
//!
//! # What crosses and what does not
//!
//! dualis's central idea is that dimensions live in the type system, so `Length + Time` does not
//! compile. **Python cannot have that**, and pretending otherwise would be the wrong thing to
//! build: a runtime `Quantity` wrapper would cost an allocation and a check per operation to
//! catch, at run time, a class of error the Rust side catches at compile time and the Python
//! side would mostly not make — because a Python caller does not do algebra on the values, it
//! passes them in and reads them out.
//!
//! So the boundary is **SI floats with the unit in the parameter name**: `length_m`, `watts`,
//! `initial_k`. One place holds the conversion, which is the same rule the Rust API follows.
//!
//! What *does* cross, and is the reason to use this rather than numpy, is the **conservation
//! audit**. Every [`Simulation.advance`] either advances the clock or raises
//! [`Violation`](struct@Violation), and the exception carries the quantity, the site, the before
//! and after, the scale and the tolerance. Write physics in numpy and a wrong model runs happily
//! and produces plausible output; here it stops and says what went missing and where. That is a
//! correctness signal a caller — a person or an agent — can act on.
//!
//! # What this first version cannot do
//!
//! **You cannot write a `Domain` in Python.** Doing so means calling back into the interpreter
//! from inside `Simulation::step`, holding the GIL across a loop that is meant to be tight, and
//! deciding what a Python exception raised mid-sweep does to a half-advanced simulation. All
//! three are answerable and none is answerable cheaply, so this version binds the domains the
//! library already has and says so rather than shipping a version of that which is subtly wrong.
//!
//! The consequence is worth stating plainly: this is enough to *run and audit* coupled physics
//! from Python, and not enough to *extend* it. If the extending is what you want, the Rust side
//! is where it lives.

use pyo3::create_exception;
use pyo3::exceptions::{PyException, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use dualis::prelude::{quantity, Energy};
use dualis::prelude::{
    Area, Bar1D, Conductance, Domain, Environment, Exchange, Kind, Ledger, Length, LumpedMass,
    Power, Schedule, Simulation as RustSimulation, Substance, Temperature, ThermalNetwork, Time,
    Violation as RustViolation, Volume, HEAT,
};

create_exception!(
    dualis,
    Violation,
    PyException,
    "A conserved quantity did not balance, or a domain refused a step it could not take.\n\n\
     Carries `quantity`, `site`, `before`, `after`, `scale` and `tolerance`. The clock has not \
     advanced: the simulation is exactly where it was before the call."
);

/// Turn a kernel violation into the exception, keeping every field addressable.
fn raise(v: RustViolation) -> PyErr {
    let message = v.to_string();
    Python::attach(|py| {
        let err = Violation::new_err(message);
        // Every field addressable, not only the message. The point of the audit is that a
        // caller can *act* on it, and parsing a sentence is not acting on it.
        let obj = err.value(py);
        let _ = obj.setattr("quantity", v.quantity.clone());
        let _ = obj.setattr("site", v.site.clone());
        let _ = obj.setattr("before", v.before);
        let _ = obj.setattr("after", v.after);
        let _ = obj.setattr("scale", v.scale);
        let _ = obj.setattr("tolerance", v.tolerance);
        err
    })
}

/// A steady source paying joules onto the heat channel out of a finite tank.
///
/// The tank is what closes the books. A source with an unlimited supply creates energy from
/// nothing every step, and the audit is right to refuse that — so the choice is between saying
/// where the energy comes from and turning the check off.
struct Heater {
    name: String,
    watts: f64,
    reserve: f64,
}

impl Domain for Heater {
    fn name(&self) -> &str {
        &self.name
    }
    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }
    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), RustViolation> {
        let joules = (self.watts * dt.to_si()).min(self.reserve).max(0.0);
        self.reserve -= joules;
        bus.publish(HEAT, joules);
        Ok(())
    }
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.reserve)
    }
    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }
}

/// A set of domains sharing a clock, with the conservation audit live on every step.
///
/// **Bound to the thread that created it.** `Domain` is not `Send`, so a `Simulation` holding
/// `Vec<Box<dyn Domain>>` cannot cross threads, and `unsendable` says so rather than pyo3
/// demanding a bound the kernel does not offer. That is the right way round: a simulation is a
/// stateful thing driven from one place, and the library's parallelism claim is about
/// `Rng::for_index` making *work items* order-independent, not about moving a running
/// simulation between threads. Touching one from another thread raises rather than racing.
#[pyclass(name = "Simulation", unsendable)]
struct PySimulation {
    sim: Option<RustSimulation>,
    /// Names in the order they were added, so `advance` can report and `mean_temperature` can
    /// say which names it knows.
    names: Vec<String>,
}

#[pymethods]
impl PySimulation {
    /// `schedule` is one of `"one-way"`, `"staggered"`, `"multirate"`, and the choice matters.
    ///
    /// `"multirate"` subcycles each domain to its own stability limit, so a frame interval need
    /// not shrink to the stiffest domain's. It buys **stability**; accuracy is set by the
    /// interval you pass to `advance`.
    #[new]
    #[pyo3(signature = (schedule = "multirate", conservation_tolerance = 1e-9))]
    fn new(schedule: &str, conservation_tolerance: f64) -> PyResult<Self> {
        let s = match schedule {
            "one-way" => Schedule::OneWay,
            "staggered" => Schedule::Staggered,
            "multirate" => Schedule::Multirate,
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown schedule {other:?}; known are \"one-way\", \"staggered\", \
                     \"multirate\""
                )))
            }
        };
        Ok(PySimulation {
            sim: Some(RustSimulation::new(s).conservation_tolerance(conservation_tolerance)),
            names: Vec::new(),
        })
    }

    /// A source of `watts` with `reserve_j` joules to spend before it goes quiet.
    ///
    /// Order matters under the staggered schedules: a publisher added before its consumer is
    /// taken from on the same step rather than the next one.
    fn add_heater(&mut self, name: &str, watts: f64, reserve_j: f64) -> PyResult<()> {
        self.push(name, |sim| {
            sim.with(Heater {
                name: name.to_string(),
                watts: watts.max(0.0),
                reserve: reserve_j.max(0.0),
            })
        })
    }

    /// A one-dimensional conducting bar, insulated at both ends.
    ///
    /// `material` is any name in the kernel's catalogue and defaults to aluminium, which is what this
    /// took before it could be said. A default is right here — most callers reaching for a bar want a
    /// metal one and should not have to look up a spelling to get it — but it has to be *sayable*, and
    /// for `add_lump` below it was not even documented.
    #[pyo3(signature = (name, length_m, cells, area_m2, initial_k = 293.15, material = "aluminium"))]
    fn add_bar(
        &mut self,
        name: &str,
        length_m: f64,
        cells: usize,
        area_m2: f64,
        initial_k: f64,
        material: &str,
    ) -> PyResult<()> {
        if cells < 2 {
            return Err(PyValueError::new_err("a bar needs at least two cells"));
        }
        // `is_finite` first, and then the comparison. The original was `!(x > 0.0)`, which
        // rejects NaN by the negation being true rather than by saying so — correct, and
        // exactly the shape clippy asks you not to write, because the NaN case is invisible.
        if !length_m.is_finite() || length_m <= 0.0 || !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(PyValueError::new_err(
                "length_m and area_m2 must be positive",
            ));
        }
        let substance = substance_named(material, name)?;
        self.push(name, |sim| {
            sim.with(Bar1D::new(
                name.to_string(),
                substance,
                cells,
                Length::from_si(length_m / cells as f64),
                Area::from_si(area_m2),
                Temperature::from_si(initial_k),
            ))
        })
    }

    /// A body at one temperature, losing heat to still air.
    ///
    /// `material` is any name in the kernel's catalogue and defaults to aluminium. It **was** aluminium
    /// and the docstring did not say so, which is the worse half of a hardcoded material: a caller
    /// modelling a copper busbar got an answer for a different metal, and the only way to find out was to
    /// read the Rust.
    #[pyo3(signature = (name, volume_m3, thickness_m, area_m2, initial_k = 293.15, ambient_k = 293.15, material = "aluminium"))]
    #[allow(clippy::too_many_arguments)]
    fn add_lump(
        &mut self,
        name: &str,
        volume_m3: f64,
        thickness_m: f64,
        area_m2: f64,
        initial_k: f64,
        ambient_k: f64,
        material: &str,
    ) -> PyResult<()> {
        let substance = substance_named(material, name)?;
        self.push(name, |sim| {
            sim.with(LumpedMass::new(
                name.to_string(),
                substance,
                Volume::from_si(volume_m3),
                Length::from_si(thickness_m),
                Temperature::from_si(initial_k),
                Environment::still_air(Temperature::from_si(ambient_k), Area::from_si(area_m2)),
            ))
        })
    }

    /// Several bodies joined by conductances: winding, stator, housing.
    ///
    /// The one shape `add_lump` cannot express. A lump reports the whole thing as one
    /// temperature, and the number that decides whether a motor survives is the *winding*,
    /// which is hotter than the case by however much the joint between them resists.
    ///
    /// `nodes` is a list of dicts: `name`, `material`, `volume_m3`, `thickness_m`, `initial_k`,
    /// and optionally `ambient_k` **and** `area_m2` together to give the node somewhere to lose
    /// heat to.
    ///
    /// `material` is any name in the kernel's catalogue — the error lists them, because a caller who
    /// guessed has no other way to find out. It used to be a five-arm match written out here, which meant
    /// four of the nine were unnameable from Python and nothing said so. `links` is a list of `{"from": ..., "to": ..., "w_per_k": ...}`. `absorbing`
    /// names the node heat off the bus lands in.
    ///
    /// ```python
    /// sim.add_network("motor",
    ///     nodes=[{"name": "winding", "material": "copper", "volume_m3": 18e-6,
    ///             "thickness_m": 2e-3, "initial_k": 298.15},
    ///            {"name": "case", "material": "aluminium", "volume_m3": 220e-6,
    ///             "thickness_m": 4e-3, "initial_k": 298.15,
    ///             "ambient_k": 298.15, "area_m2": 0.042}],
    ///     links=[{"from": "winding", "to": "case", "w_per_k": 0.9}],
    ///     absorbing="winding")
    /// ```
    ///
    /// **Dicts rather than handles.** The Rust API addresses nodes by a `Node` value that can
    /// only come from a constructor, so a link to a node that does not exist cannot be written
    /// — which matters, because a link adds `+q` to one node and `−q` to another in the same
    /// sum, so a missing or misdirected link passes the conservation audit at machine precision.
    /// That guarantee cannot cross into Python, where a name is just a string. So the names are
    /// resolved here, once, at construction, and a name that is not a node raises **before** any
    /// stepping happens and lists the ones that are.
    #[pyo3(signature = (name, nodes, links, absorbing))]
    fn add_network(
        &mut self,
        name: &str,
        nodes: Vec<Bound<'_, PyDict>>,
        links: Vec<Bound<'_, PyDict>>,
        absorbing: &str,
    ) -> PyResult<()> {
        if nodes.is_empty() {
            return Err(PyValueError::new_err(format!(
                "{name:?}: a network needs at least one node"
            )));
        }

        let mut net = ThermalNetwork::new(name.to_string());
        for (i, n) in nodes.iter().enumerate() {
            let at = format!("{name:?} node {i}");
            let label = need_str(n, "name", &at)?;
            let material = need_str(n, "material", &at)?;
            // `substance_named`, not a match written out here. This was a five-arm match against a
            // catalogue of nine, so `borosilicate`, `ice`, `stainless_304` and `water` could not be named
            // from Python at all -- the third copy of the catalogue's spelling in this workspace, and the
            // third to go stale. `dualis-world` had the same defect and it cost eleven releases.
            let substance = substance_named(&material, &at)?;
            let volume = Volume::from_si(need_f64(n, "volume_m3", &at)?);
            let thickness = Length::from_si(need_f64(n, "thickness_m", &at)?);
            let initial = Temperature::from_si(need_f64(n, "initial_k", &at)?);

            // Both or neither. A node given an ambient with no area, or an area with no
            // ambient, is a caller who meant it to lose heat and will get a node that does not
            // — silently, and looking exactly like a hot interior node should look.
            let ambient = optional_f64(n, "ambient_k", &at)?;
            let area = optional_f64(n, "area_m2", &at)?;
            match (ambient, area) {
                (Some(a_k), Some(a_m2)) => {
                    net.node_losing_to(
                        label,
                        substance,
                        volume,
                        thickness,
                        initial,
                        Environment::still_air(Temperature::from_si(a_k), Area::from_si(a_m2)),
                    );
                }
                (None, None) => {
                    net.node(label, substance, volume, thickness, initial);
                }
                (a, _) => {
                    let (given, missing) = if a.is_some() {
                        ("ambient_k", "area_m2")
                    } else {
                        ("area_m2", "ambient_k")
                    };
                    return Err(PyValueError::new_err(format!(
                        "{at}: has {given} but not {missing}; a node loses heat only when it \
                         has both, and one alone would give you an interior node that looks \
                         like it is cooling and is not"
                    )));
                }
            }
        }

        // Resolved before anything is linked, so a bad name raises rather than half-building.
        let known = || {
            net.handles()
                .map(|(_, l)| l.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        };
        let mut wired = Vec::with_capacity(links.len());
        for (i, l) in links.iter().enumerate() {
            let at = format!("{name:?} link {i}");
            let from = need_str(l, "from", &at)?;
            let to = need_str(l, "to", &at)?;
            let w = need_f64(l, "w_per_k", &at)?;
            let a = net.node_named(&from).ok_or_else(|| {
                PyValueError::new_err(format!(
                    "{at}: no node named {from:?}; this network has {}",
                    known()
                ))
            })?;
            let b = net.node_named(&to).ok_or_else(|| {
                PyValueError::new_err(format!(
                    "{at}: no node named {to:?}; this network has {}",
                    known()
                ))
            })?;
            wired.push((a, b, w));
        }
        let sink = net.node_named(absorbing).ok_or_else(|| {
            PyValueError::new_err(format!(
                "{name:?}: absorbing names {absorbing:?}, which is not a node; this network \
                 has {}",
                known()
            ))
        })?;
        for (a, b, w) in wired {
            net.link(a, b, Conductance::w_per_k(w))
                .map_err(|e| PyValueError::new_err(format!("{name:?}: {e}")))?;
        }
        net.absorbing(sink)
            .map_err(|e| PyValueError::new_err(format!("{name:?}: {e}")))?;

        self.push(name, |sim| sim.with(net))
    }

    /// One node's temperature, in kelvin.
    fn node_temperature(&self, name: &str, node: &str) -> PyResult<f64> {
        let net = self.network(name)?;
        let handle = net.node_named(node).ok_or_else(|| {
            PyValueError::new_err(format!(
                "{name:?} has no node called {node:?}; it has {}",
                net.handles().map(|(_, l)| l).collect::<Vec<_>>().join(", ")
            ))
        })?;
        Ok(net.temperature(handle).to_si())
    }

    /// Every node's temperature, in kelvin, in the order the nodes were declared.
    ///
    /// Returned as pairs rather than a dict so the declaration order survives — which is the
    /// order heat flows along a chain, and therefore the order you want to read them in.
    fn node_temperatures(&self, name: &str) -> PyResult<Vec<(String, f64)>> {
        let net = self.network(name)?;
        Ok(net
            .handles()
            .map(|(n, l)| (l.to_string(), net.temperature(n).to_si()))
            .collect())
    }

    /// Heat flowing along a link right now, in watts, positive from `a` to `b`.
    ///
    /// Zero for a pair with no link between them, which is the same answer as a link of zero
    /// conductance and is not a mistake being hidden: ask `node_temperatures` what exists.
    fn heat_flow_w(&self, name: &str, a: &str, b: &str) -> PyResult<f64> {
        let net = self.network(name)?;
        let pick = |who: &str| {
            net.node_named(who).ok_or_else(|| {
                PyValueError::new_err(format!(
                    "{name:?} has no node called {who:?}; it has {}",
                    net.handles().map(|(_, l)| l).collect::<Vec<_>>().join(", ")
                ))
            })
        };
        Ok(net.heat_flow(pick(a)?, pick(b)?).to_si())
    }

    /// Where a network settles for a steady `watts`, without marching there.
    ///
    /// Returns `(name, kelvin)` pairs in declaration order. The network is not modified — this
    /// is a question asked of it, so asking mid-run does not disturb the run.
    ///
    /// Worth more from Python than from Rust: the loop that marches to a steady state crosses
    /// this boundary once per step, so `for _ in range(900): sim.advance(1.0)` pays nine hundred
    /// times for an answer one solve gives — and gives exactly, rather than to whatever the
    /// transient has decayed to.
    ///
    /// Raises if heat has nowhere to go, because then there is no steady state and the network
    /// warms without limit. A finite number would be the wrong answer to a question with none.
    fn steady_state(&self, name: &str, watts: f64) -> PyResult<Vec<(String, f64)>> {
        let net = self.network(name)?;
        let settled = net.steady_state(Power::from_si(watts)).map_err(raise)?;
        Ok(net
            .handles()
            .map(|(n, label)| (label.to_string(), settled.temperature(n).to_si()))
            .collect())
    }

    /// Advance every domain by `dt` seconds, auditing the crossing and the books.
    ///
    /// Raises `Violation` without moving the clock if a domain refuses the step, if anything
    /// published is left unconsumed, or if the totalled ledgers moved by more than the
    /// tolerance. Returns the substeps each domain took, which is how you find out whether the
    /// interval you chose was near anything's limit.
    fn advance(&mut self, dt: f64) -> PyResult<Vec<(String, u32)>> {
        let sim = self.borrow_mut()?;
        let report = sim.advance(Time::from_si(dt)).map_err(raise)?;
        Ok(report.substeps)
    }

    /// Simulation time, in seconds.
    #[getter]
    fn time(&self) -> PyResult<f64> {
        Ok(self.borrow()?.time().to_si())
    }

    /// The names this simulation knows, in the order they were added.
    #[getter]
    fn domains(&self) -> Vec<String> {
        self.names.clone()
    }

    /// A bar's mean temperature, or a lump's temperature, in kelvin.
    ///
    /// The *mean over cells*, for a bar, and not a sampled average — those differ by about
    /// `1/2n` on a cell-centred grid, which is a trap the Rust side documents and this avoids
    /// by not offering the sampled one.
    fn temperature(&self, name: &str) -> PyResult<f64> {
        let sim = self.borrow()?;
        if let Some(bar) = sim.domain_as::<Bar1D>(name) {
            return Ok(bar.mean_temperature().to_si());
        }
        if let Some(lump) = sim.domain_as::<LumpedMass>(name) {
            return Ok(lump.temperature().to_si());
        }
        // A network has no single temperature and will not invent one by averaging: the whole
        // reason to build one is that its nodes differ, and a mean would be a number that
        // describes no part of it. Refused by name, pointing at the two calls that answer.
        if sim.domain_as::<ThermalNetwork>(name).is_some() {
            return Err(PyValueError::new_err(format!(
                "{name:?} is a network, which has one temperature per node rather than one \
                 overall; use node_temperatures({name:?}) or node_temperature({name:?}, node)"
            )));
        }
        Err(self.unknown(name))
    }

    /// A bar's temperature profile, cell by cell, in kelvin.
    fn profile(&self, name: &str) -> PyResult<Vec<f64>> {
        let sim = self.borrow()?;
        let bar = sim
            .domain_as::<Bar1D>(name)
            .ok_or_else(|| self.unknown(name))?;
        Ok((0..bar.cell_count())
            .map(|i| bar.temperature_at(i).to_si())
            .collect())
    }

    /// Joules a domain has taken off the bus over the run.
    fn absorbed_j(&self, name: &str) -> PyResult<f64> {
        let sim = self.borrow()?;
        if let Some(bar) = sim.domain_as::<Bar1D>(name) {
            return Ok(bar.absorbed_energy().to_si());
        }
        if let Some(lump) = sim.domain_as::<LumpedMass>(name) {
            return Ok(lump.absorbed_energy().to_si());
        }
        if let Some(net) = sim.domain_as::<ThermalNetwork>(name) {
            return Ok(net.absorbed_energy().to_si());
        }
        Err(self.unknown(name))
    }

    /// Joules a heater has left to spend.
    fn reserve_j(&self, name: &str) -> PyResult<f64> {
        let sim = self.borrow()?;
        let h = sim
            .domain_as::<Heater>(name)
            .ok_or_else(|| self.unknown(name))?;
        Ok(h.reserve)
    }

    /// Every conserved quantity the domains report, summed. The audit's own view.
    fn ledger(&self) -> PyResult<Vec<(String, f64)>> {
        Ok(self
            .borrow()?
            .ledger()
            .quantities()
            .map(|(k, v)| (k.to_string(), v))
            .collect())
    }

    fn __repr__(&self) -> String {
        match &self.sim {
            Some(s) => format!(
                "<dualis.Simulation t={:.6}s domains={:?}>",
                s.time().to_si(),
                self.names
            ),
            None => "<dualis.Simulation poisoned>".to_string(),
        }
    }
}

impl PySimulation {
    /// The builder methods take `Simulation` by value, so the field is an `Option` and every
    /// path puts it back. A panic between the take and the put would poison the object rather
    /// than corrupt it, and `borrow` says so instead of unwrapping.
    fn push(
        &mut self,
        name: &str,
        add: impl FnOnce(RustSimulation) -> RustSimulation,
    ) -> PyResult<()> {
        if self.names.iter().any(|n| n == name) {
            return Err(PyValueError::new_err(format!(
                "there is already a domain called {name:?}; names are how they are looked up, \
                 so the second would be invisible"
            )));
        }
        let sim = self.sim.take().ok_or_else(Self::poisoned)?;
        self.sim = Some(add(sim));
        self.names.push(name.to_string());
        Ok(())
    }

    fn borrow(&self) -> PyResult<&RustSimulation> {
        self.sim.as_ref().ok_or_else(Self::poisoned)
    }

    fn borrow_mut(&mut self) -> PyResult<&mut RustSimulation> {
        self.sim.as_mut().ok_or_else(Self::poisoned)
    }

    fn poisoned() -> PyErr {
        PyValueError::new_err("this Simulation was left incomplete by an earlier failure")
    }

    fn unknown(&self, name: &str) -> PyErr {
        PyValueError::new_err(format!(
            "no domain called {name:?} of a type that answers this; known: {:?}",
            self.names
        ))
    }

    fn network(&self, name: &str) -> PyResult<&ThermalNetwork> {
        self.borrow()?
            .domain_as::<ThermalNetwork>(name)
            .ok_or_else(|| self.unknown(name))
    }
}

/// A material name from the kernel's catalogue, or an error listing every name there is.
///
/// One function, because a name resolved in three places is a catalogue copied three times — which is
/// what this binding had, and it held five of the nine.
fn substance_named(material: &str, at: &str) -> PyResult<Substance> {
    Substance::from_name(material).ok_or_else(|| {
        PyValueError::new_err(format!(
            "{at}: unknown material {material:?}; known are {}",
            Substance::CATALOGUE
                .iter()
                .map(|m| format!("{m:?}"))
                .collect::<Vec<_>>()
                .join(", ")
        ))
    })
}

/// A required key, with the error naming the dict it was missing from.
///
/// `PyDict::get_item` returning `None` and a key holding `None` are different things and both
/// are wrong here, so both take the same message rather than one of them extracting to a
/// default. A silent default is how a node ends up with a volume of zero and a capacity of
/// zero, which the network then refuses with a message about heat capacity rather than about
/// the key the caller forgot.
/// Concrete rather than generic over `FromPyObject`: the two types wanted are `str` and `f64`,
/// and the generic form needs a `where` clause on pyo3's associated error type that buys
/// nothing here.
fn need_str(d: &Bound<'_, PyDict>, key: &str, at: &str) -> PyResult<String> {
    match d.get_item(key)? {
        Some(v) if !v.is_none() => v
            .extract()
            .map_err(|_| PyValueError::new_err(format!("{at}: {key} should be a string"))),
        _ => Err(PyValueError::new_err(format!("{at}: missing {key:?}"))),
    }
}

fn need_f64(d: &Bound<'_, PyDict>, key: &str, at: &str) -> PyResult<f64> {
    match d.get_item(key)? {
        Some(v) if !v.is_none() => v
            .extract()
            .map_err(|_| PyValueError::new_err(format!("{at}: {key} should be a number"))),
        _ => Err(PyValueError::new_err(format!("{at}: missing {key:?}"))),
    }
}

/// An optional key. Absent and explicitly `None` are the same answer here, deliberately:
/// `{"area_m2": None}` reads as "no area" to anyone writing it.
fn optional_f64(d: &Bound<'_, PyDict>, key: &str, at: &str) -> PyResult<Option<f64>> {
    match d.get_item(key)? {
        Some(v) if !v.is_none() => v
            .extract()
            .map(Some)
            .map_err(|_| PyValueError::new_err(format!("{at}: {key} should be a number"))),
        _ => Ok(None),
    }
}

/// What a heat capacity is, for computing a closed form on the Python side.
///
/// Exposed because a caller checking its own results needs the same constants the simulation used, and
/// reading them off the simulation would be checking it against itself.
///
/// It was `aluminium_heat_capacity_j_per_k` and took no material, which was fine while a bar and a lump
/// were both aluminium and stopped being fine the moment they were not: a test of a copper bar would have
/// had to compare it against aluminium's capacity or hardcode copper's, and one of those is wrong and the
/// other is the constant this function exists to stop people copying.
#[pyfunction]
fn heat_capacity_j_per_k(material: &str, volume_m3: f64) -> PyResult<f64> {
    substance_named(material, "heat_capacity_j_per_k")?
        .heat_capacity(Volume::from_si(volume_m3))
        .map(|c| c.to_si())
        .ok_or_else(|| {
            PyValueError::new_err(format!(
                "{material:?} does not state a specific heat, so it has no heat capacity"
            ))
        })
}

/// One joule, in joules. A smoke test that the units layer crossed intact.
#[pyfunction]
fn one_joule() -> f64 {
    Energy::from_si(1.0).to_si()
}

#[pymodule]
#[pyo3(name = "dualis")]
fn dualis_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__doc__", include_str!("../README.md"))?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("Violation", m.py().get_type::<Violation>())?;
    m.add_class::<PySimulation>()?;
    m.add_function(wrap_pyfunction!(heat_capacity_j_per_k, m)?)?;
    m.add_function(wrap_pyfunction!(one_joule, m)?)?;
    Ok(())
}
