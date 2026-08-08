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

use dualis::prelude::{
    Area, Bar1D, Domain, Environment, Exchange, Kind, Ledger, Length, LumpedMass, Schedule,
    Simulation as RustSimulation, Substance, Temperature, Time, Violation as RustViolation, Volume,
    HEAT,
};
use dualis::prelude::{quantity, Energy};

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

    /// A one-dimensional conducting bar of aluminium, insulated at both ends.
    #[pyo3(signature = (name, length_m, cells, area_m2, initial_k = 293.15))]
    fn add_bar(
        &mut self,
        name: &str,
        length_m: f64,
        cells: usize,
        area_m2: f64,
        initial_k: f64,
    ) -> PyResult<()> {
        if cells < 2 {
            return Err(PyValueError::new_err("a bar needs at least two cells"));
        }
        if !(length_m > 0.0) || !(area_m2 > 0.0) {
            return Err(PyValueError::new_err(
                "length_m and area_m2 must be positive",
            ));
        }
        self.push(name, |sim| {
            sim.with(Bar1D::new(
                name.to_string(),
                Substance::aluminium_6061(),
                cells,
                Length::from_si(length_m / cells as f64),
                Area::from_si(area_m2),
                Temperature::from_si(initial_k),
            ))
        })
    }

    /// A body at one temperature, losing heat to still air.
    #[pyo3(signature = (name, volume_m3, thickness_m, area_m2, initial_k = 293.15, ambient_k = 293.15))]
    fn add_lump(
        &mut self,
        name: &str,
        volume_m3: f64,
        thickness_m: f64,
        area_m2: f64,
        initial_k: f64,
        ambient_k: f64,
    ) -> PyResult<()> {
        self.push(name, |sim| {
            sim.with(LumpedMass::new(
                name.to_string(),
                Substance::aluminium_6061(),
                Volume::from_si(volume_m3),
                Length::from_si(thickness_m),
                Temperature::from_si(initial_k),
                Environment::still_air(
                    Temperature::from_si(ambient_k),
                    Area::from_si(area_m2),
                ),
            ))
        })
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
}

/// What a heat capacity is, for computing a closed form on the Python side.
///
/// Exposed because a caller checking its own results needs the same constants the simulation
/// used, and reading them off the simulation would be checking it against itself.
#[pyfunction]
fn aluminium_heat_capacity_j_per_k(volume_m3: f64) -> PyResult<f64> {
    Substance::aluminium_6061()
        .heat_capacity(Volume::from_si(volume_m3))
        .map(|c| c.to_si())
        .ok_or_else(|| PyValueError::new_err("aluminium has a specific heat; this cannot happen"))
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
    m.add_function(wrap_pyfunction!(aluminium_heat_capacity_j_per_k, m)?)?;
    m.add_function(wrap_pyfunction!(one_joule, m)?)?;
    Ok(())
}
