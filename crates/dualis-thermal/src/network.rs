//! A lumped thermal network: n nodes, conductances between them, one domain.
//!
//! # Why this and not two coupled `LumpedMass` domains
//!
//! The standard electromechanical model is **junction → case → ambient**: a winding that gets
//! hot, a housing it conducts into, and a room the housing loses to. That is the model that
//! answers "will the winding survive", which is usually the only thermal question that matters.
//!
//! It cannot be built out of two [`LumpedMass`](crate::LumpedMass) domains, and the reason is
//! structural rather than an omission. A conductance carries `UA·(T₁ − T₂)`: it needs **both**
//! temperatures. Domains in this workspace never read each other — they meet on an
//! [`Exchange`], which carries *amounts* and not state, so neither side
//! can publish a temperature and neither can compute the flux alone. Any `conducting_to(peer)`
//! would have to break the property the rest of the design rests on.
//!
//! So the network is **one domain holding many nodes**, which is also what a thermal network
//! physically is: a single coupled system of ODEs, not independent bodies posting parcels to
//! each other. One `ledger()`, one stability limit, and the conservation audit unchanged.
//!
//! It buys something [`Bar1D`](crate::Bar1D) cannot express either — a **contact** resistance
//! between different materials. A bolted joint or a winding pressed into a stator is not bulk
//! conduction through one substance, and modelling it on a uniform grid needs a fictitious
//! conductivity standing in for the real interface.
//!
//! # What the audit cannot see here, and what covers it instead
//!
//! A link contributes `+q` to one node and `−q` to another **in the same sum**. They cancel
//! identically, so the ledger is blind to links by construction: a sign error, a transposed
//! index or a link dropped altogether passes the conservation audit at machine precision.
//!
//! That is not a gap in the kernel — [`audit_transfers`](dualis_core::Exchange::audit_transfers)
//! covers transfers *between* domains and this one is inside a domain — but it decides how the
//! tests are written. Every link check is per node or on a closed form, never on the total. It
//! is the same lesson as the per-face audit in space and the substep share in time, arriving a
//! third time.
//!
//! It also decided the API. A node is addressed by a [`Node`] handle rather than by name,
//! because a link naming a node that does not exist would be exactly the invisible case above:
//! the books balance, the winding runs hot forever, and the number is plausible. A handle can
//! only come from [`ThermalNetwork::node`], so a dangling reference is not representable.
//! Names still exist — every node carries a label, and [`ThermalNetwork::node_named`] is the
//! bridge for a caller building from a file, where the resolver is a fifteen-line loop that can
//! name the file's own vocabulary in its error.

use dualis_core::conserved::quantity;
use dualis_core::{Domain, Exchange, Kind, Ledger, Substance, Violation};
use dualis_units::{Conductance, Energy, Length, Power, Temperature, Time, Volume};

use crate::{Environment, HEAT};

/// A node in one particular network.
///
/// Carries the network's identity as well as the index, so a handle from one network used on
/// another is refused rather than silently addressing whatever sits at that index. The identity
/// is a hash of the network's name, which is deterministic: no counter, no clock, no global
/// state, and [`Simulation`](dualis_core::Simulation) already refuses two domains with one name.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Node {
    index: u32,
    network: u64,
}

/// One body in a network: a capacity, a temperature, and possibly somewhere to lose heat to.
struct NodeState {
    label: String,
    /// J/K. Held rather than recomputed, so a node is a capacity and not a substance — which is
    /// what lets a network mix copper, steel and air without the caller assembling a fiction.
    capacity: f64,
    /// For the Biot number, which is the only thing that still needs a length.
    thickness: f64,
    substance: Substance,
    temperature: f64,
    /// What `ledger` measures against. The *initial* temperature and not ambient, following
    /// `Bar1D`: differencing absolute enthalpies leaves a rounding floor that gets worse on
    /// refinement, and an interior node has no ambient to measure from at all.
    reference: f64,
    environment: Option<Environment>,
}

/// A conductance between two nodes, in W/K.
struct Link {
    a: usize,
    b: usize,
    ua: f64,
}

/// A network of lumped bodies joined by conductances.
pub struct ThermalNetwork {
    name: String,
    id: u64,
    nodes: Vec<NodeState>,
    links: Vec<Link>,
    /// Which node heat off the bus arrives at. `None` until set, and a network that is never
    /// told will leave anything published unclaimed, which the bus refuses by itself.
    absorbing: Option<usize>,
    absorbed: f64,
    lost: f64,
    saved: Option<Saved>,
}

/// Everything `ledger` reads. All of it, because saving the state and not the running totals is
/// how a rewound sweep comes to report heat it never shed — see `LumpedMass::checkpoint`, which
/// did exactly that until an iterative coupling was finally built that could reach the branch.
type Saved = (Vec<f64>, f64, f64);

/// A deterministic identity for a network, from its name.
fn identity(name: &str) -> u64 {
    // FNV-1a. Not for security; for telling two networks apart without a counter or a clock,
    // which the determinism rule forbids.
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in name.as_bytes() {
        h ^= *b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

impl ThermalNetwork {
    /// An empty network. Add nodes, then link them.
    pub fn new(name: impl Into<String>) -> ThermalNetwork {
        let name = name.into();
        let id = identity(&name);
        ThermalNetwork {
            name,
            id,
            nodes: Vec::new(),
            links: Vec::new(),
            absorbing: None,
            absorbed: 0.0,
            lost: 0.0,
            saved: None,
        }
    }

    /// An interior node: it conducts to its neighbours and loses to nothing.
    ///
    /// A winding inside a housing has no room to convect into — only metal to conduct through —
    /// and giving it an `Environment` with the numbers zeroed would need *two* knobs right
    /// (`convection` **and** `area`, or the radiative term still runs). Absence of a loss path
    /// is the absence of a thing, and this workspace spells that with `Option` everywhere else:
    /// `Substance::thermal`, `Bar1D::boundary`.
    pub fn node(
        &mut self,
        label: impl Into<String>,
        substance: Substance,
        volume: Volume,
        thickness: Length,
        initial: Temperature,
    ) -> Node {
        self.push_node(label, substance, volume, thickness, initial, None)
    }

    /// A node that also loses heat to its surroundings, like a [`LumpedMass`](crate::LumpedMass).
    pub fn node_losing_to(
        &mut self,
        label: impl Into<String>,
        substance: Substance,
        volume: Volume,
        thickness: Length,
        initial: Temperature,
        environment: Environment,
    ) -> Node {
        self.push_node(
            label,
            substance,
            volume,
            thickness,
            initial,
            Some(environment),
        )
    }

    fn push_node(
        &mut self,
        label: impl Into<String>,
        substance: Substance,
        volume: Volume,
        thickness: Length,
        initial: Temperature,
        environment: Option<Environment>,
    ) -> Node {
        let capacity = substance
            .heat_capacity(volume)
            .map(|c| c.to_si())
            .unwrap_or(f64::NAN);
        let index = self.nodes.len() as u32;
        self.nodes.push(NodeState {
            label: label.into(),
            capacity,
            thickness: thickness.to_si(),
            substance,
            temperature: initial.to_si(),
            reference: initial.to_si(),
            environment,
        });
        Node {
            index,
            network: self.id,
        }
    }

    /// Join two nodes by a conductance, in W/K.
    ///
    /// Refuses a self-link, a negative conductance and a handle from a different network. Two
    /// links between the same pair **accumulate**, because parallel conductances add — the same
    /// convention [`Exchange::publish`](dualis_core::Exchange::publish) follows for repeated
    /// offers on one channel.
    pub fn link(&mut self, a: Node, b: Node, ua: Conductance) -> Result<(), Violation> {
        let (i, j) = (self.resolve(a)?, self.resolve(b)?);
        if i == j {
            return Err(Violation::at(
                format!("{}/{}", self.name, self.nodes[i].label),
                "a node cannot conduct to itself",
                0.0,
            ));
        }
        let w = ua.to_si();
        // `!is_finite` first, so NaN is rejected by the branch that reads as rejecting it rather
        // than by a negated comparison that happens to be false.
        if !w.is_finite() || w < 0.0 {
            return Err(Violation::at(
                format!(
                    "{}/{}-{}",
                    self.name, self.nodes[i].label, self.nodes[j].label
                ),
                "conductance must be finite and not negative",
                w,
            ));
        }
        if let Some(existing) = self
            .links
            .iter_mut()
            .find(|l| (l.a == i && l.b == j) || (l.a == j && l.b == i))
        {
            existing.ua += w;
        } else {
            self.links.push(Link { a: i, b: j, ua: w });
        }
        Ok(())
    }

    /// Where heat taken off the bus arrives.
    ///
    /// A network takes a share of the [`HEAT`] channel like any other consumer, but the joules
    /// have no *place* — so, exactly as [`Bar1D`](crate::Bar1D) puts lumped heat in its first
    /// cell, the network needs telling which node it lands on. Naming it is better than
    /// defaulting to node zero, which would be a silent choice about the physics.
    pub fn absorbing(&mut self, node: Node) -> Result<(), Violation> {
        self.absorbing = Some(self.resolve(node)?);
        Ok(())
    }

    fn resolve(&self, node: Node) -> Result<usize, Violation> {
        if node.network != self.id {
            return Err(Violation::at(
                self.name.clone(),
                "a node handle from a different network",
                node.index as f64,
            ));
        }
        let i = node.index as usize;
        if i >= self.nodes.len() {
            return Err(Violation::at(
                self.name.clone(),
                "a node handle from a later state of this network",
                node.index as f64,
            ));
        }
        Ok(i)
    }

    /// A node's temperature. Total, because a handle cannot name a node that is not here.
    pub fn temperature(&self, node: Node) -> Temperature {
        Temperature::from_si(self.nodes[node.index as usize].temperature)
    }

    /// How far a node has moved from where it started.
    pub fn rise(&self, node: Node) -> Temperature {
        let n = &self.nodes[node.index as usize];
        Temperature::from_si(n.temperature - n.reference)
    }

    /// The label a node was given, for a violation or a legend.
    pub fn label(&self, node: Node) -> &str {
        &self.nodes[node.index as usize].label
    }

    /// A node by label, for a caller that built the network from a file and has names rather
    /// than handles. The seam between a name-shaped format and a handle-shaped API.
    pub fn node_named(&self, label: &str) -> Option<Node> {
        self.nodes
            .iter()
            .position(|n| n.label == label)
            .map(|i| Node {
                index: i as u32,
                network: self.id,
            })
    }

    /// Every node, in the order they were added, with its label.
    ///
    /// A [`Node`] can only come from [`node`](ThermalNetwork::node) or
    /// [`node_losing_to`](ThermalNetwork::node_losing_to), which is what makes a dangling link
    /// unrepresentable — and it also meant a caller holding a network it did not build could not
    /// walk it at all. [`nodes`](ThermalNetwork::nodes) returned a count and there was no way to
    /// turn a count into anything. Found the first time the consumer app tried to print a
    /// network's temperatures, which is a smaller version of `FRICTION.md`'s recurring finding:
    /// the API is comfortable when the parts are known at compile time and awkward the moment
    /// they are not.
    pub fn handles(&self) -> impl Iterator<Item = (Node, &str)> + '_ {
        let id = self.id;
        self.nodes.iter().enumerate().map(move |(i, n)| {
            (
                Node {
                    index: i as u32,
                    network: id,
                },
                n.label.as_str(),
            )
        })
    }

    /// How many nodes.
    pub fn nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Watts crossing the joint between two nodes right now, positive from `a` to `b`.
    ///
    /// Zero if they are not linked, which is a real answer rather than a missing one.
    pub fn heat_flow(&self, a: Node, b: Node) -> Power {
        let (Ok(i), Ok(j)) = (self.resolve(a), self.resolve(b)) else {
            return Power::from_si(0.0);
        };
        let w = self
            .links
            .iter()
            .find(|l| (l.a == i && l.b == j) || (l.a == j && l.b == i))
            .map(|l| l.ua)
            .unwrap_or(0.0);
        Power::from_si(w * (self.nodes[i].temperature - self.nodes[j].temperature))
    }

    /// Heat taken off the bus over the run.
    pub fn absorbed_energy(&self) -> Energy {
        Energy::from_si(self.absorbed)
    }

    /// Heat shed to the environments over the run.
    pub fn lost_energy(&self) -> Energy {
        Energy::from_si(self.lost)
    }

    /// Whether a node's lumped approximation applies, from its own conductivity and its own
    /// environment. `None` for an interior node, which has no film coefficient to compare
    /// against — and that is worth knowing rather than papering over.
    pub fn biot_number(&self, node: Node) -> Option<f64> {
        let n = &self.nodes[node.index as usize];
        let e = n.environment.as_ref()?;
        let k = n.substance.thermal.map(|t| t.conductivity.to_si())?;
        if k <= 0.0 {
            return None;
        }
        Some(e.convection_w_per_m2_k * n.thickness / k)
    }

    /// A node's local time constant: its capacity over everything that carries heat away from
    /// it, links included.
    pub fn time_constant(&self, node: Node) -> Time {
        let i = node.index as usize;
        let g = self.node_conductance(i);
        if g <= 0.0 || !self.nodes[i].capacity.is_finite() {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(self.nodes[i].capacity / g)
    }

    /// Everything conducting heat out of node `i`: its environment linearised at its current
    /// temperature, plus every link on it.
    fn node_conductance(&self, i: usize) -> f64 {
        let n = &self.nodes[i];
        let env = n
            .environment
            .as_ref()
            .map(|e| {
                crate::linearised_loss_conductance(
                    e,
                    Temperature::from_si(n.temperature),
                    n.substance.thermal.map(|t| t.emissivity).unwrap_or(0.0),
                )
            })
            .unwrap_or(0.0);
        env + self
            .links
            .iter()
            .filter(|l| l.a == i || l.b == i)
            .map(|l| l.ua)
            .sum::<f64>()
    }
}

impl Domain for ThermalNetwork {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// A tenth of the fastest node's time constant.
    ///
    /// The explicit limit is set by whichever node empties quickest, and in a network that is
    /// often **not** a node touching the outside: a small winding on a stiff link to a large
    /// housing is far faster than the housing is, and its limit comes entirely from the link.
    /// That is why [`ThermalNetwork`]'s violation names the node rather than only the network.
    ///
    /// A tenth for the same reason [`LumpedMass`](crate::LumpedMass) reports one: explicit Euler
    /// is stable to `2τ` and accurate nowhere near it, so the number a scheduler should honour
    /// is the accuracy limit. With one node and no links this is `LumpedMass::max_stable_dt`
    /// exactly.
    ///
    /// State-dependent, because the environment term is linearised at each node's current
    /// temperature — so the limit tightens as the network heats.
    fn max_stable_dt(&self, _now: Time) -> Time {
        let mut limit = f64::INFINITY;
        for i in 0..self.nodes.len() {
            let g = self.node_conductance(i);
            let c = self.nodes[i].capacity;
            if g > 0.0 && c.is_finite() && c > 0.0 {
                limit = limit.min(c / g);
            }
        }
        Time::from_si(limit / 10.0)
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let h = dt.to_si();
        if h <= 0.0 {
            return Ok(());
        }
        if self.nodes.is_empty() {
            return Err(Violation::at(
                self.name.clone(),
                "a network with no nodes",
                0.0,
            ));
        }

        // Refuse rather than diverge, and name the node. `F_i` is the fraction of a node's
        // capacity its own update moves in one step; past one the coefficient `1 - F_i` goes
        // negative and the temperature oscillates. Same criterion `Bar1D` applies to its
        // Fourier number, and the useful difference is that in a network the fast node is not
        // obvious, so the site says which it is.
        for i in 0..self.nodes.len() {
            let c = self.nodes[i].capacity;
            if !c.is_finite() || c <= 0.0 {
                return Err(Violation::at(
                    format!("{}/{}", self.name, self.nodes[i].label),
                    "substance has no heat capacity",
                    c,
                ));
            }
            let f = h * self.node_conductance(i) / c;
            if f > 1.0 + 1e-12 {
                return Err(Violation {
                    quantity: "network Fourier number".to_string(),
                    site: format!(
                        "{}/{} (explicit RC network)",
                        self.name, self.nodes[i].label
                    ),
                    before: 1.0,
                    after: f,
                    scale: 1.0,
                    tolerance: 1e-12,
                });
            }
        }

        // This substep's share of the channel, not the whole outer step's. See
        // `Exchange::take_share`.
        let gained = bus.take_share(HEAT, dt);
        self.absorbed += gained;
        if self.absorbing.is_none() && gained != 0.0 {
            return Err(Violation::at(
                self.name.clone(),
                "heat arrived but no node was named to absorb it",
                gained,
            ));
        }

        // Jacobi: every flux from the same snapshot, so the answer does not depend on the order
        // links were declared in. Gauss-Seidel would be a different scheme, not a rounding
        // difference.
        let before: Vec<f64> = self.nodes.iter().map(|n| n.temperature).collect();
        let mut delta = vec![0.0; self.nodes.len()];

        // The arriving heat is a term of the same right-hand side as the fluxes, so it joins
        // `delta` rather than being added to the temperature first.
        //
        // Applying it first is the obvious reading and it is wrong in a way conservation cannot
        // see: the absorbing node's outgoing flux is then computed from the already-raised
        // temperature, and the steady state — which explicit Euler otherwise reaches *exactly* —
        // acquires a bias of `K·h/C` on the joint next to the source. The excess simply lands in
        // the neighbour, so every total stays right. Measured on a three-node ladder before the
        // fix: the first joint sat 0.31% low against `P/K`, against a predicted `Kh/C = 0.0031006`
        // and an observed 0.0031005, while the far joint and the environment drop were exact to
        // six figures. That agreement is what identified it; the audit was silent throughout.
        if let Some(i) = self.absorbing {
            delta[i] += gained;
        }

        // Each link's flux computed **once** and applied twice with opposite signs. Computing it
        // from each side separately gives two values differing in the last bit, and the network
        // then leaks about 1e-16 per link per step — a drift that is invisible per step and is
        // not invisible over a long run.
        for l in &self.links {
            let q = l.ua * (before[l.a] - before[l.b]) * h;
            delta[l.a] -= q;
            delta[l.b] += q;
        }

        for (i, n) in self.nodes.iter().enumerate() {
            if let Some(e) = &n.environment {
                let lost = e
                    .loss_from(
                        Temperature::from_si(before[i]),
                        n.substance.thermal.map(|t| t.emissivity).unwrap_or(0.0),
                    )
                    .to_si()
                    * h;
                delta[i] -= lost;
                self.lost += lost;
            }
        }

        for (i, n) in self.nodes.iter_mut().enumerate() {
            n.temperature += delta[i] / n.capacity;
        }
        Ok(())
    }

    /// What every node is holding, plus what the environments have taken.
    ///
    /// One `add` per node rather than one for the sum, so `Ledger`'s scale is the largest single
    /// node's holding instead of a near-zero net — which is what the scale exists for, and the
    /// mistake `NBody::ledger` made until a test proved its momentum audit was inert.
    ///
    /// Measured from each node's **initial** temperature. An interior node has no ambient to
    /// measure from, and `Bar1D`'s reasoning applies anyway: differencing absolute enthalpies
    /// leaves a rounding floor that gets worse on refinement.
    fn ledger(&self) -> Ledger {
        let mut ledger = Ledger::new();
        for n in &self.nodes {
            ledger.add(quantity::ENERGY, n.capacity * (n.temperature - n.reference));
        }
        ledger.add(quantity::ENERGY, self.lost);
        ledger
    }

    /// Every node's temperature and both running totals.
    ///
    /// All of it, because `ledger` reads all of it. `LumpedMass` saved its temperature and not
    /// its `lost`, and a rewound iterative sweep therefore reported heat it had already shed —
    /// which went unnoticed for as long as it did because nothing in the workspace had a
    /// residual, so the restore branch never ran.
    fn checkpoint(&mut self) {
        self.saved = Some((
            self.nodes.iter().map(|n| n.temperature).collect(),
            self.absorbed,
            self.lost,
        ));
    }

    fn restore(&mut self) {
        if let Some((temps, absorbed, lost)) = self.saved.clone() {
            for (n, t) in self.nodes.iter_mut().zip(temps) {
                n.temperature = t;
            }
            self.absorbed = absorbed;
            self.lost = lost;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    /// **`None`, and not as an oversight.**
    ///
    /// A network is a graph with no embedding. Its nodes have capacities, not positions, and a
    /// conductance is not a distance — two nodes joined by 0.8 W/K are not "0.8 apart" in any
    /// space a field could be sampled over. Interpolating between them would invent a continuum
    /// with less justification than a box of atoms or a set of orbits would have, and those
    /// decline too.
    fn as_field(&self) -> Option<&dyn dualis_core::ScalarField> {
        None
    }
}
