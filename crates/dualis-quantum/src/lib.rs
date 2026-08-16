//! dualis-quantum: the Schrödinger equation, as a domain on the `dualis-core` kernel.
//!
//! One particle in one dimension, inside a well whose walls are infinitely hard:
//!
//! ```text
//! iℏ ∂ψ/∂t = −(ℏ²/2m) ∂²ψ/∂x² + V(x)ψ,      ψ = 0 at both walls
//! ```
//!
//! marched with Visscher's staggered real/imaginary scheme (P. B. Visscher, *A fast explicit
//! algorithm for the time-dependent Schrödinger equation*, Computers in Physics **5**, 596,
//! 1991): the real part of ψ lives at whole steps and the imaginary part half a step later,
//! each updated explicitly from the other through the discrete Hamiltonian. It is the same
//! staggered-leapfrog family as the acoustic domains, and it earns its place here the same
//! way the Yee grid earns its place in the electromagnetic one — by conserving the quantity
//! that matters **as an identity of the update** rather than as an accuracy claim. The
//! discrete probability
//!
//! ```text
//! P = Σⱼ [ Rⱼ(t)² + Iⱼ(t+dt/2)·Iⱼ(t−dt/2) ] · dx
//! ```
//!
//! telescopes to a boundary term under one step, and the boundary is a wall where ψ is zero,
//! so the update cannot change it at all — for a fixed step. `P` pairs `I` across two
//! consecutive steps, so changing `dt` re-pairs the levels and moves it once, by `O(dt²)`,
//! the same way starting the leapfrog does; the tests measure both. That is what this
//! domain's [`Ledger`] carries, under the name [`PROBABILITY`], normalised to 1.
//!
//! # What it is checked against
//!
//! Closed forms, never a second implementation:
//!
//! - **The discrete eigenvalues, exactly.** On interior nodes with ψ = 0 walls, the
//!   three-point Hamiltonian's eigenvectors *are* the sampled `sin(nπx/L)`, with eigenvalue
//!   `E = (2ℏ²/m·dx²)·sin²(nπ·dx/2L)` — so a sampled eigenstate's density is stationary
//!   under the march and its phase rotates at a rate the closed form predicts to rounding.
//! - **Second-order convergence** of those eigenvalues to `n²π²ℏ²/2mL²`, measured across
//!   grid doublings rather than asserted at one resolution.
//! - **Free wavepacket spreading**, `σ(t)² = σ₀²(1 + (ℏt/2mσ₀²)²)`.
//! - **Ehrenfest's theorem on a harmonic potential**: for quadratic `V`, `⟨x⟩` follows the
//!   classical cosine exactly in the continuum, and the discretisation error must fall at
//!   second order.
//!
//! # What is deliberately not in it
//!
//! One particle, which is the honest boundary of a grid method: a second particle squares
//! the state space rather than doubling the array, and pretending otherwise is how a
//! "quantum" module ends up simulating classical waves. No spin, no second or third
//! dimension, no open or absorbing boundaries (the walls are infinite and reflect — a
//! complex absorbing potential would make the probability books *supposed* to fall, which is
//! a different domain), no time-dependent potentials, and no measurement. Momentum is not on
//! the ledger, and neither is ⟨H⟩: the ledger carries exactly the one quantity the scheme
//! conserves identically, and the energy expectation is a [`Domain::readings`] entry checked
//! for drift in tests instead.
//!
//! # The eleventh domain, and what it had to add to the kernel
//!
//! Nothing. Optics, heat, mechanics, sound, molecules, electricity, percolation, elasticity,
//! electromagnetism and fluid flow each arrived without the kernel changing, and this crate
//! extends that claim rather than testing it for the first time: `probability` is a ledger
//! name the kernel has never heard of, carried by the same [`Ledger`] that carries joules.
//!
//! # Units and scales
//!
//! SI throughout, like everywhere else in the workspace: metres, kilograms, seconds, joules.
//! An electron in a nanometre well sits comfortably inside `f64` — energies around 10⁻¹⁹ J,
//! stable steps of a few attoseconds, periods of femtoseconds — so there is no internal unit
//! system to convert in and out of.

// Every public item carries a doc comment. Denied rather than warned: a public physics API
// whose constructor shows a blank summary in rustdoc is documented in the sense that a
// paragraph exists somewhere, and not in the sense a reader needs.
#![deny(missing_docs)]

use dualis_core::{Domain, Exchange, Kind, Ledger, Reading, ScalarField, Violation};
use dualis_units::{Energy, Frequency, Length, LengthVec, Mass, Time};

/// The reduced Planck constant, `1.054 571 817 × 10⁻³⁴ J·s` (CODATA, exact by the 2019 SI
/// definition of the kilogram: `h` is fixed and `ℏ = h/2π`).
const HBAR: f64 = 1.054_571_817e-34;

/// The ledger channel this domain conserves: total probability, dimensionless, normalised
/// to 1.
///
/// Not one of the kernel's named quantities, and that is the point: the kernel's
/// `conserved::quantity` module says a domain crate can add its own name without the kernel
/// being edited, and this is the first crate to hold it to that. Nothing else publishes or
/// consumes probability, so the channel never crosses the bus — the ledger entry exists for
/// the audit, which checks that the march creates and destroys nothing.
pub const PROBABILITY: &str = "probability";

/// A spatial angular frequency `k`, in radians per metre. `ℏk` is a momentum.
///
/// The wavenumber a Gaussian packet is launched with — see [`Well::with_gaussian`]. Defined
/// here the way `dualis-acoustic` defines its `Impedance`: a `Qty` alias for a dimension the
/// units crate has no everyday name for. Build one with `Wavenumber::from_si`.
pub type Wavenumber = dualis_units::Qty<-1, 0, 0, 0, 0, 0, 0>;

/// Everything [`Well::checkpoint`] has to put back: both time levels of both parts of ψ,
/// whether the leapfrog has been staggered yet, and the startup offset.
///
/// The flag is state. `Room` in the acoustic crate carries the same bool for the same
/// reason, and a restore that forgot it would silently re-run the leapfrog startup — a
/// second half-kick on a state that already had one, which no audit would see.
#[derive(Clone)]
struct Snapshot {
    re: Vec<f64>,
    re_prev: Vec<f64>,
    im: Vec<f64>,
    im_prev: Vec<f64>,
    staggered: bool,
    norm_offset: f64,
}

/// A single particle in a one-dimensional infinite square well, on a uniform grid.
///
/// # The grid, stated exactly
///
/// **Interior nodes only, node-centred**: `cells` samples at `x_j = j·dx` for
/// `j = 1 … cells`, with `dx = width/(cells+1)`, so the walls sit *on* the grid lines
/// `x = 0` and `x = width` just outside the stored array, where ψ is identically zero.
/// Node-centred like the acoustic `Room` rather than cell-centred like the thermal `Bar1D`,
/// and the reason is the boundary condition: a Dirichlet wall is a statement about a value
/// *at a point*, and putting that point on a grid line makes the sampled `sin(nπx/L)` an
/// **exact eigenvector of the discrete Hamiltonian** — the property
/// [`Well::in_eigenstate`] promises. A cell-centred grid would put the wall between
/// samples and the sampled sine would only be an approximate eigenvector.
///
/// # The scheme
///
/// Visscher's staggered leapfrog. `R = Re ψ` at whole steps, `I = Im ψ` half a step later:
///
/// ```text
/// R^{n+1}   = R^n     + (dt/ℏ)·H I^{n+1/2}
/// I^{n+3/2} = I^{n+1/2} − (dt/ℏ)·H R^{n+1}
/// ```
///
/// with `(Hψ)_j = −(ℏ²/2m)(ψ_{j+1} − 2ψ_j + ψ_{j−1})/dx² + V_j ψ_j` and ψ = 0 beyond the
/// walls. Explicit, second order in `dt` and `dx`, and at a fixed step it conserves the
/// staggered probability `Σ [R² + I⁺I⁻]·dx` identically — see the crate docs.
///
/// An initial condition supplies `I` at `t = 0`, not at `t = dt/2` where the scheme carries
/// it, so the first step half-kicks `I` before marching — the same startup the acoustic
/// domains perform on velocity, with the same property: for an eigenstate the half-kick is
/// not merely second-order accurate but *exact* for the discrete rotation, which is what
/// makes the stationarity claim above hold to rounding rather than to `O(dt)`.
pub struct Well {
    name: String,
    /// `Re ψ` at the interior nodes, at whole steps. Units m⁻¹ᐟ², like ψ itself in 1D.
    re: Vec<f64>,
    /// `Re ψ` one whole step earlier, kept so the field's reported rate can be exactly the
    /// step just taken. See [`ScalarField::rate`] on this type.
    re_prev: Vec<f64>,
    /// `Im ψ`, half a step *ahead* of `re` once staggered.
    im: Vec<f64>,
    /// `Im ψ` one whole step behind `im` — half a step *behind* `re` — which is the other
    /// factor in the conserved density `R² + I⁺I⁻`.
    im_prev: Vec<f64>,
    /// `V` at the interior nodes, in joules. Zero unless [`Well::with_harmonic`] set it,
    /// and never negative — the stability bound leans on that.
    potential: Vec<f64>,
    dx: f64,
    mass: f64,
    /// Whether `im` has been offset to the half step the scheme carries it at.
    staggered: bool,
    /// What the released state's probability was: the datum the ledger reports against.
    norm_datum: f64,
    /// The datum minus the staggered invariant, fixed at the first step. `O(dt²)`; see
    /// [`Well::startup_adjustment`].
    norm_offset: f64,
    saved: Option<Snapshot>,
}

impl Well {
    /// An empty well: a particle of the given mass confined to `width`, with `cells`
    /// interior samples and no wavefunction yet.
    ///
    /// Fewer than three cells is not a grid; the count is rounded up to three. Give the
    /// well a state with [`Well::in_eigenstate`] or [`Well::with_gaussian`] — an empty well
    /// steps, conserves its zero, and reports a norm of zero.
    pub fn new(name: impl Into<String>, mass: Mass, width: Length, cells: usize) -> Well {
        let cells = cells.max(3);
        Well {
            name: name.into(),
            re: vec![0.0; cells],
            re_prev: vec![0.0; cells],
            im: vec![0.0; cells],
            im_prev: vec![0.0; cells],
            potential: vec![0.0; cells],
            dx: width.to_si() / (cells + 1) as f64,
            mass: mass.to_si(),
            staggered: false,
            norm_datum: 0.0,
            norm_offset: 0.0,
            saved: None,
        }
    }

    /// Start in the `n`th eigenstate of the bare well: `ψ ∝ sin(nπx/L)`, real, normalised
    /// so `Σ|ψ|²·dx = 1`.
    ///
    /// `n` runs from 1 (the ground state, one hump) and is clamped to `1..=cells`, because
    /// `n` interior nodes carry exactly `n` distinct eigenvectors — a higher `n` would
    /// alias onto a lower one with sign flips rather than mean anything new.
    ///
    /// This is an **exact eigenvector of the discrete Hamiltonian**, not just a sample of
    /// the continuum one, so under the march its density does not move at all and its phase
    /// rotates at exactly [`Well::eigenstate_energy`]`/ℏ` (through the scheme's own
    /// dispersion — see [`Domain::max_stable_dt`] on this type). With a potential from
    /// [`Well::with_harmonic`] in place it is no longer an eigenstate of anything, just an
    /// initial condition.
    pub fn in_eigenstate(mut self, n: usize) -> Well {
        let count = self.re.len();
        let n = n.clamp(1, count);
        for j in 1..=count {
            let phase = n as f64 * std::f64::consts::PI * j as f64 / (count + 1) as f64;
            self.re[j - 1] = phase.sin();
            self.im[j - 1] = 0.0;
        }
        self.normalise()
    }

    /// Start with a normalised Gaussian wavepacket,
    /// `ψ ∝ exp(−(x−x₀)²/4σ²)·exp(ik₀x)`.
    ///
    /// `σ` is the *position* standard deviation of `|ψ|²`, and it is clamped below at one
    /// grid spacing — a packet narrower than a cell cannot be represented, and silently
    /// storing an alias of one would be worse than widening it. `k₀` gives the packet a
    /// mean momentum `ℏk₀`; zero is a packet at rest, which only spreads.
    ///
    /// Normalisation is done on the grid, over the well as it is — so a packet placed close
    /// enough to a wall for its tail to be clipped still integrates to one.
    pub fn with_gaussian(mut self, centre: Length, sigma: Length, k0: Wavenumber) -> Well {
        let count = self.re.len();
        let sigma = sigma.to_si().max(self.dx);
        let x0 = centre.to_si();
        let k = k0.to_si();
        for j in 1..=count {
            let x = j as f64 * self.dx;
            let envelope = (-(x - x0) * (x - x0) / (4.0 * sigma * sigma)).exp();
            self.re[j - 1] = envelope * (k * x).cos();
            self.im[j - 1] = envelope * (k * x).sin();
        }
        self.normalise()
    }

    /// Add a harmonic potential `V(x) = ½·m·ω²·(x − L/2)²` on the interior, centred on the
    /// well.
    ///
    /// `omega` is the **angular** frequency, in rad·s⁻¹ — the ω of `½mω²x²`, whose quantum
    /// energy spacing is `ℏω` — not a cycles-per-second frequency; the units crate carries
    /// both on the same `s⁻¹` dimension, so the distinction has to live here in words.
    /// Negative values are treated as zero. The potential is non-negative by construction,
    /// which the stability bound in [`Domain::max_stable_dt`] relies on, and it tightens
    /// that bound: a stiffer spring means faster phase rotation at the edges.
    pub fn with_harmonic(mut self, omega: Frequency) -> Well {
        let count = self.re.len();
        let w = omega.to_si().max(0.0);
        let centre = 0.5 * (count + 1) as f64 * self.dx;
        for j in 1..=count {
            let x = j as f64 * self.dx - centre;
            self.potential[j - 1] = 0.5 * self.mass * w * w * x * x;
        }
        self
    }

    /// Normalise to `Σ|ψ|²·dx = 1` and reset the leapfrog to an unstaggered start.
    fn normalise(mut self) -> Well {
        let mut total = 0.0;
        for (r, i) in self.re.iter().zip(self.im.iter()) {
            total += r * r + i * i;
        }
        total *= self.dx;
        if total > 0.0 {
            let inv = 1.0 / total.sqrt();
            for v in self.re.iter_mut().chain(self.im.iter_mut()) {
                *v *= inv;
            }
        }
        self.re_prev.copy_from_slice(&self.re);
        self.im_prev.copy_from_slice(&self.im);
        self.staggered = false;
        self.norm_offset = 0.0;
        self.norm_datum = self.invariant();
        self
    }

    /// How many interior samples the well carries.
    pub fn cells(&self) -> usize {
        self.re.len()
    }

    /// Width of the well, wall to wall: `(cells + 1)·dx`, because the walls sit one spacing
    /// outside the first and last stored nodes.
    pub fn width(&self) -> Length {
        Length::from_si((self.re.len() + 1) as f64 * self.dx)
    }

    /// The grid spacing.
    pub fn spacing(&self) -> Length {
        Length::from_si(self.dx)
    }

    /// Total probability, in the staggered form the scheme conserves identically, measured
    /// against the released state as its datum. One, for a normalised state; zero for an
    /// empty well.
    pub fn norm(&self) -> f64 {
        self.invariant() + self.norm_offset
    }

    /// How much the staggered invariant differs from the released state's `Σ|ψ|²·dx`.
    ///
    /// Zero until the first step, because the invariant pairs `I` across a step and there
    /// is no step before then. `O(dt²)` and converging — for an eigenstate it is exactly
    /// `sin²(E·dt/2ℏ)` — and it is taken as a datum rather than handed to the audit,
    /// the same trade the acoustic `Room` makes at its leapfrog startup and for the same
    /// reason: it is discretisation, not a leak, and it is bounded (at 25% of the total)
    /// so a real first-step bug cannot hide in it.
    ///
    /// The same pairing has a consequence mid-run: **changing `dt` moves the invariant
    /// once**, by the same `O(dt²)`, because the two `I` levels it multiplies were then
    /// made with different steps. The datum is not re-taken — the shift is real,
    /// discretisation-sized, and a caller who varies the outer step will see the norm
    /// carry it. Under a steady cadence, which is the scheduler's normal use, it never
    /// appears.
    pub fn startup_adjustment(&self) -> f64 {
        self.norm_offset
    }

    /// The `n`th eigenvalue of the **discrete** Hamiltonian of the bare well (potential
    /// excluded), exactly:
    ///
    /// ```text
    /// E_n = (2ℏ²/m·dx²)·sin²(nπ·dx/2L)
    /// ```
    ///
    /// This is the energy [`Well::in_eigenstate`] actually rotates at, and it converges to
    /// the continuum `n²π²ℏ²/2mL²` from below at second order in `dx` — the tests measure
    /// that rate rather than trusting this formula. `n` is clamped to `1..=cells` as in
    /// `in_eigenstate`.
    pub fn eigenstate_energy(&self, n: usize) -> Energy {
        let count = self.re.len();
        let n = n.clamp(1, count);
        if self.mass <= 0.0 {
            return Energy::ZERO;
        }
        let kinetic_sup = 2.0 * HBAR * HBAR / (self.mass * self.dx * self.dx);
        let phase = n as f64 * std::f64::consts::PI / (2.0 * (count + 1) as f64);
        Energy::from_si(kinetic_sup * phase.sin().powi(2))
    }

    /// The position expectation `⟨x⟩`, from the conserved density. An empty well reports
    /// its centre, which is where a symmetric nothing would be.
    pub fn position_expectation(&self) -> Length {
        let count = self.re.len();
        let (mut weight, mut moment) = (0.0, 0.0);
        for j in 1..=count {
            let rho = self.node_density(j as isize);
            weight += rho;
            moment += rho * j as f64 * self.dx;
        }
        if weight.abs() < 1e-300 {
            return Length::from_si(0.5 * (count + 1) as f64 * self.dx);
        }
        Length::from_si(moment / weight)
    }

    /// The energy expectation `⟨H⟩`, in the time-staggered form the scheme conserves:
    /// `Σ [R·HR + I⁺·H I⁻]·dx`, normalised by the probability.
    ///
    /// For a sampled eigenstate this is exactly [`Well::eigenstate_energy`], constant in
    /// time; for a general state it is constant to rounding, which the tests pin. Reported
    /// through [`Domain::readings`] and deliberately **not** on the ledger: the ledger
    /// carries the one quantity the update conserves as an identity, and mixing an `O(dt²)`
    /// estimator into the same books would dilute exactly the property that makes them
    /// worth auditing.
    pub fn energy_expectation(&self) -> Energy {
        let weight = self.invariant();
        if weight.abs() < 1e-300 {
            return Energy::ZERO;
        }
        let h_re = self.h_of(&self.re);
        let h_im_behind = self.h_of(&self.im_prev);
        let mut energy = 0.0;
        for j in 0..self.re.len() {
            energy += self.re[j] * h_re[j] + self.im[j] * h_im_behind[j];
        }
        Energy::from_si(energy * self.dx / weight)
    }

    /// The step against the stability limit, as a ratio that must not exceed 1 — this
    /// domain's Courant number.
    pub fn stability_ratio(&self, dt: Time) -> f64 {
        let limit = self.max_stable_dt(Time::ZERO).to_si();
        if limit <= 0.0 {
            return f64::INFINITY;
        }
        dt.to_si() / limit
    }

    /// The scheme's invariant: `Σ [R² + I⁺I⁻]·dx`, before the release datum is applied.
    fn invariant(&self) -> f64 {
        let mut total = 0.0;
        for j in 0..self.re.len() {
            total += self.re[j] * self.re[j] + self.im[j] * self.im_prev[j];
        }
        total * self.dx
    }

    /// The discrete Hamiltonian applied to a field: three-point kinetic term with ψ = 0
    /// walls, plus the potential.
    fn h_of(&self, f: &[f64]) -> Vec<f64> {
        let count = f.len();
        let k = HBAR * HBAR / (2.0 * self.mass * self.dx * self.dx);
        let mut out = vec![0.0; count];
        for j in 0..count {
            let left = if j == 0 { 0.0 } else { f[j - 1] };
            let right = if j + 1 == count { 0.0 } else { f[j + 1] };
            out[j] = -k * (left - 2.0 * f[j] + right) + self.potential[j] * f[j];
        }
        out
    }

    /// One raw Visscher step, staggering first if this is the first.
    ///
    /// No stability check here — that is [`Domain::step`]'s job, and the split is what lets
    /// the unit tests drive the scheme past the limit to show the divergence the public API
    /// refuses to produce.
    fn advance(&mut self, h: f64) {
        // Half a kick the first time: the initial condition put `I` at t = 0 and the
        // scheme wants it at t = dt/2. For an eigenstate this half-kick lands exactly on
        // the discrete rotation, not merely near it: the stationary value is
        // −sin(θ/2)·v with sin(θ/2) = E·dt/2ℏ, and the half-kick computes −(dt/2ℏ)·E·v,
        // which is the same number by the definition of θ.
        if !self.staggered {
            let h_re = self.h_of(&self.re);
            for (i, hr) in self.im.iter_mut().zip(h_re.iter()) {
                *i -= h / (2.0 * HBAR) * hr;
            }
            self.staggered = true;
        }

        let h_im = self.h_of(&self.im);
        self.re_prev.copy_from_slice(&self.re);
        for (r, hi) in self.re.iter_mut().zip(h_im.iter()) {
            *r += h / HBAR * hi;
        }

        let h_re = self.h_of(&self.re);
        self.im_prev.copy_from_slice(&self.im);
        for (i, hr) in self.im.iter_mut().zip(h_re.iter()) {
            *i -= h / HBAR * hr;
        }
    }

    /// The conserved probability density at interior node `j` (1-based), in m⁻¹; zero on
    /// and beyond the walls. `isize` so a stencil can ask about `j − 1` at the first node
    /// without a special case.
    fn node_density(&self, j: isize) -> f64 {
        if j < 1 || j > self.re.len() as isize {
            return 0.0;
        }
        let k = (j - 1) as usize;
        self.re[k] * self.re[k] + self.im[k] * self.im_prev[k]
    }

    /// Nearest interior node (1-based) to a position, clamped into the well.
    fn node_at(&self, p: LengthVec) -> usize {
        let u = p.to_si().x / self.dx;
        if u.is_nan() || u < 1.0 {
            return 1;
        }
        (u.round() as usize).min(self.re.len())
    }

    /// The exact change per second of the conserved density at node `j` over the step just
    /// taken — the discrete continuity equation. See [`ScalarField::rate`] on this type.
    ///
    /// Exact for a step taken from a staggered state. Across the very first step the
    /// observed density change additionally carries the one-time staggering of `I`, which
    /// is the same `O(dt²)` startup [`Well::startup_adjustment`] accounts for globally.
    fn node_rate(&self, j: usize) -> f64 {
        if self.mass <= 0.0 {
            return 0.0;
        }
        let count = self.re.len() as isize;
        let j = j as isize;
        // S = R + R_prev spans the step; I⁻ is the half-step level between them.
        let s = |i: isize| -> f64 {
            if i < 1 || i > count {
                0.0
            } else {
                self.re[(i - 1) as usize] + self.re_prev[(i - 1) as usize]
            }
        };
        let ib = |i: isize| -> f64 {
            if i < 1 || i > count {
                0.0
            } else {
                self.im_prev[(i - 1) as usize]
            }
        };
        let k = HBAR / (2.0 * self.mass * self.dx * self.dx);
        -k * (ib(j + 1) * s(j) + ib(j - 1) * s(j) - ib(j) * s(j + 1) - ib(j) * s(j - 1))
    }
}

impl Domain for Well {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// `2ℏ / (2ℏ²/m·dx² + V_max)`, derived and not decorated.
    ///
    /// # The derivation
    ///
    /// In the eigenbasis of the discrete Hamiltonian each mode with eigenvalue `E` obeys
    /// the scalar leapfrog `r ← r + ε·i`, `i ← i − ε·r` with `ε = E·dt/ℏ` — a rotation
    /// integrated explicitly, whose per-step phase `θ` satisfies `sin(θ/2) = ε/2`. That has
    /// a real solution, and therefore bounded orbits, only while `|ε| ≤ 2`: the step must
    /// keep `dt ≤ 2ℏ/|E|` for **every** eigenvalue. The three-point kinetic spectrum on
    /// interior nodes is `(2ℏ²/m·dx²)·sin²(kπ·dx/2L)`, strictly inside `(0, 2ℏ²/m·dx²)`,
    /// and the potential this crate builds is non-negative, shifting eigenvalues up by at
    /// most `V_max` — so `2ℏ²/m·dx² + V_max` bounds the spectrum and the limit follows.
    /// The bound is sharp to `O((π·dx/2L)²)`: the top eigenvalue sits at
    /// `cos²(π·dx/2L)` of the kinetic supremum, so almost nothing is given away.
    ///
    /// For the bare well this is `m·dx²/ℏ` — a **diffusion-shaped** limit, not a CFL one.
    /// Schrödinger is first order in time and second in space, like heat and unlike sound,
    /// so halving the grid costs four times the steps; that is why the stable step for an
    /// electron on a 20 pm grid is attoseconds against femtosecond periods.
    ///
    /// # What the limit guarantees, precisely: a *fixed* step
    ///
    /// Like every leapfrog, the scheme is stable under the bound for a constant `dt`, and
    /// that is the guarantee. A `dt` that **cycles periodically** can parametrically pump
    /// the top of the band even with every individual step under the limit: two alternating
    /// steps compose to `M(ε₂)M(ε₁)` with trace `2 − (ε₁+ε₂)² + ε₁²ε₂²`, which leaves
    /// `[−2, 2]` once `ε₁ + ε₂` nears 2 — and the growth seeds itself from rounding, so it
    /// does not need the state to hold any top-band amplitude. Measured: a 3-cycle of 0.9,
    /// 0.37 and 0.71 of the limit overflows a wavepacket run within a few thousand steps,
    /// while any one of the three held fixed is stable indefinitely. Changing `dt` *once*
    /// is harmless for stability (it does re-pair the invariant — see
    /// [`Well::startup_adjustment`]); equal substeps — which is what
    /// [`Schedule::Multirate`] takes — are what the limit is a promise about.
    ///
    /// [`Schedule::Multirate`]: dualis_core::Schedule::Multirate
    ///
    /// # Why the audit cannot stand in for this limit
    ///
    /// Past the limit the top mode's rotation turns hyperbolic and the field grows without
    /// bound — but the ledger's probability is conserved by the update *unconditionally*,
    /// growth and all: the exploding `R²` is cancelled by an exploding negative `I⁺I⁻` to
    /// within rounding on the entries' own scale. A diverging run balances its books, so
    /// stability has to be enforced by refusal in [`Domain::step`]; no conservation check
    /// would ever catch it. The tests demonstrate both halves.
    fn max_stable_dt(&self, _now: Time) -> Time {
        if self.mass <= 0.0 {
            // A massless particle's kinetic scale is unbounded: no positive step is
            // stable, and zero is the honest report.
            return Time::from_si(0.0);
        }
        let kinetic_sup = 2.0 * HBAR * HBAR / (self.mass * self.dx * self.dx);
        let v_max = self.potential.iter().fold(0.0f64, |a, v| a.max(*v));
        Time::from_si(2.0 * HBAR / (kinetic_sup + v_max))
    }

    fn step(&mut self, _t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let ratio = self.stability_ratio(dt);
        if ratio > 1.0 + 1e-12 {
            return Err(Violation {
                quantity: "stability ratio".to_string(),
                site: format!("{} (explicit Schrödinger leapfrog)", self.name),
                before: 1.0,
                after: ratio,
                scale: 1.0,
                tolerance: 1e-12,
            });
        }
        if dt.to_si() <= 0.0 {
            return Ok(());
        }
        let starting = !self.staggered;
        self.advance(dt.to_si());

        // Startup, once: staggering `I` moves the invariant off the released `Σ|ψ|²·dx`
        // by O(dt²) — for an eigenstate, by exactly sin²(E·dt/2ℏ). That is discretisation
        // and not a leak, so it is taken as a datum, bounded so a real first-step bug — a
        // sign, a factor of two, a wall — cannot shelter under it.
        if starting {
            let invariant = self.invariant();
            let offset = self.norm_datum - invariant;
            let scale = self.norm_datum.abs().max(invariant.abs());
            if scale > 0.0 && offset.abs() / scale > 0.25 {
                return Err(Violation {
                    quantity: PROBABILITY.to_string(),
                    site: format!("{} (starting the leapfrog)", self.name),
                    before: self.norm_datum,
                    after: invariant,
                    scale,
                    tolerance: 0.25,
                });
            }
            self.norm_offset = offset;
        }
        Ok(())
    }

    fn ledger(&self) -> Ledger {
        Ledger::new().with(PROBABILITY, self.norm())
    }

    /// True, and exactly: this domain touches no bus channel, and its ledger is the
    /// scheme-exact invariant, so its books balance to rounding every step.
    fn books_balance(&self) -> bool {
        true
    }

    fn checkpoint(&mut self) {
        self.saved = Some(Snapshot {
            re: self.re.clone(),
            re_prev: self.re_prev.clone(),
            im: self.im.clone(),
            im_prev: self.im_prev.clone(),
            staggered: self.staggered,
            norm_offset: self.norm_offset,
        });
    }

    fn restore(&mut self) {
        if let Some(s) = self.saved.clone() {
            self.re = s.re;
            self.re_prev = s.re_prev;
            self.im = s.im;
            self.im_prev = s.im_prev;
            self.staggered = s.staggered;
            self.norm_offset = s.norm_offset;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    /// The norm, the position expectation and the energy expectation — the three numbers
    /// this domain is *for* when nobody is drawing it.
    ///
    /// Not a mean of the field: a probability density's mean over the well is `1/L`
    /// whatever the state is doing, which would be a column of the same number.
    fn readings(&self) -> Vec<Reading> {
        vec![
            Reading::new(&self.name, "norm", self.norm(), ""),
            Reading::new(&self.name, "<x>", self.position_expectation().to_si(), "m"),
            Reading::new(&self.name, "<E>", self.energy_expectation().to_si(), "J"),
        ]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }

    /// The well reads as a probability-density field, so a renderer never has to know it
    /// is quantum.
    fn as_field(&self) -> Option<&dyn ScalarField> {
        Some(self)
    }
}

/// The well as a probability density `|ψ|²`, in m⁻¹.
///
/// The density reported is the **conserved** one, `R² + I⁺I⁻` — the staggered product whose
/// sum the update preserves identically — rather than `R² + I²` read off mismatched time
/// levels; the two agree to `O(dt²)` and only one of them is what the ledger holds. A
/// consequence a renderer should know: a single node's value can dip fractionally negative
/// where the phase swept past between the two half-steps. The total cannot.
///
/// # What is honest here and what is declined
///
/// The well lies along `x`; `y` and `z` are ignored, which is the claim a one-dimensional
/// model makes rather than an approximation it hides. `at` interpolates linearly between
/// nodes and through the exact zeros on both walls, and reports **zero beyond them** — an
/// infinite wall is not a place where a value is clamped, it is a place where the particle
/// cannot be.
///
/// `laplacian` uses the grid's own three-point stencil on the nodal densities, wall zeros
/// included, snapped to the nearest node — the same reason `Bar1D` gives: the interpolant
/// is piecewise linear, so *its* second derivative is zero between nodes and a spike at
/// them, and differencing it would report nonsense.
///
/// `rate` is exact, not estimated: the update's own discrete continuity equation, so it
/// reports precisely the change per second the step just taken produced — the test holds it
/// to rounding. Like the acoustic `Room`'s rate it is therefore half a statement about the
/// past rather than the future, and for the same staggered-leapfrog reason. Overriding it
/// matters here more than usual: `at` ignores the time argument (a marched domain holds
/// *now* and nothing else), so the default finite difference in time would return zero for
/// a field that is visibly moving.
///
/// `gradient` is the one operator left on the trait's default central difference over `at`.
/// Overriding it would need the kernel's vector type from `glam`, and this crate
/// deliberately depends on the kernel and the units crate alone; the default sampled at
/// `h ≥ dx` reports the interpolant's own slope, which is the grid's, and sampling finer
/// than the grid cannot add information it does not have.
impl ScalarField for Well {
    /// Per metre — a probability density in one dimension, whose integral over the well
    /// is the dimensionless norm.
    fn unit(&self) -> &'static str {
        "1/m"
    }

    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let u = p.to_si().x / self.dx;
        // NaN and both outsides read as zero: there is no wavefunction beyond an
        // infinite wall.
        if u.is_nan() || u <= 0.0 || u >= (self.re.len() + 1) as f64 {
            return 0.0;
        }
        let i = u.floor() as isize;
        let f = u - i as f64;
        self.node_density(i) * (1.0 - f) + self.node_density(i + 1) * f
    }

    /// The three-point stencil on the nodal densities, with the exact wall zeros, at the
    /// nearest node.
    fn laplacian(&self, p: LengthVec, _t: Time, _h: Length) -> f64 {
        let j = self.node_at(p) as isize;
        (self.node_density(j - 1) - 2.0 * self.node_density(j) + self.node_density(j + 1))
            / (self.dx * self.dx)
    }

    /// `∂|ψ|²/∂t`, exactly as the step just taken changed it: the discrete continuity
    /// equation
    ///
    /// ```text
    /// Δρ_j/dt = −(ℏ/2m·dx²)·[ I⁻ⱼ₊₁Sⱼ + I⁻ⱼ₋₁Sⱼ − I⁻ⱼSⱼ₊₁ − I⁻ⱼSⱼ₋₁ ],   S = R + R_prev
    /// ```
    ///
    /// whose telescoping sum is the reason the ledger's probability cannot move. The
    /// potential cancels out of it identically, which is the discrete image of `V` driving
    /// phase and not transport. Snapped to the nearest node, like the Laplacian, and `dt`
    /// drops out of the algebra so the argument is ignored.
    fn rate(&self, p: LengthVec, _t: Time, _dt: Time) -> f64 {
        self.node_rate(self.node_at(p))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// The electron mass, kg (CODATA).
    const M_E: f64 = 9.109_383_701_5e-31;

    fn electron_in(width_nm: f64, cells: usize) -> Well {
        Well::new("well", Mass::kg(M_E), Length::nm(width_nm), cells)
    }

    /// Total probability, mean position and variance from the conserved density.
    fn moments(w: &Well) -> (f64, f64, f64) {
        let n = w.cells();
        let (mut p, mut m1) = (0.0, 0.0);
        for j in 1..=n {
            let rho = w.node_density(j as isize);
            p += rho * w.dx;
            m1 += rho * w.dx * (j as f64 * w.dx);
        }
        let mean = m1 / p;
        let mut m2 = 0.0;
        for j in 1..=n {
            let x = j as f64 * w.dx - mean;
            m2 += w.node_density(j as isize) * w.dx * x * x;
        }
        (p, mean, m2 / p)
    }

    /// **The identity, and its exact edge.** At a fixed step the staggered probability is
    /// conserved by the update itself — like `∇·B` on a Yee grid — so it holds to
    /// accumulated rounding, not to a discretisation tolerance. The 1e-12 is arithmetic:
    /// each step recomputes a ~100-term sum whose entries are O(1/N), so the per-step
    /// wobble is a few 1e-16 and a random walk over 2000 steps stays far under it. This is
    /// not an accuracy claim about ψ, which carries the ordinary O(dx², dt²) error the
    /// other tests measure.
    ///
    /// The edge is measured rather than hidden: `P` multiplies `I` from two consecutive
    /// steps, so **changing dt re-pairs those levels and moves it once**, by O(dt²) —
    /// 1.3e-4 here for 0.9 → 0.55 of the limit — after which the identity resumes at the
    /// new step to the same 1e-12. (The first version of this test *cycled* three step
    /// sizes to show dt-independence and found parametric divergence instead; that
    /// hazard is documented on [`Domain::max_stable_dt`].)
    #[test]
    fn the_probability_invariant_is_an_identity_of_the_update() {
        let mut w = electron_in(2.0, 101)
            .with_harmonic(Frequency::from_si(1.5e15))
            .with_gaussian(
                Length::nm(1.2),
                Length::nm(0.15),
                Wavenumber::from_si(4.0e9),
            );
        let limit = w.max_stable_dt(Time::ZERO);
        let mut bus = Exchange::new();

        w.step(Time::ZERO, limit * 0.9, &mut bus).unwrap();
        let datum = w.norm();
        assert!(
            (datum - 1.0).abs() < 1e-9,
            "the datum is the released norm: {datum}"
        );

        for _ in 0..2000u32 {
            w.step(Time::ZERO, limit * 0.9, &mut bus).unwrap();
        }
        assert!(
            (w.norm() - datum).abs() < 1e-12,
            "P drifted from {datum} to {} at a fixed step",
            w.norm()
        );

        // A different step, once: a single bounded re-pairing, O(dt²) and not a leak...
        w.step(Time::ZERO, limit * 0.55, &mut bus).unwrap();
        let shift = (w.norm() - datum).abs();
        assert!(
            shift > 1e-8 && shift < 1e-3,
            "the re-pairing is real and discretisation-sized: {shift:e}"
        );

        // ...after which the identity holds again, as tightly as before.
        let datum = w.norm();
        for _ in 0..2000u32 {
            w.step(Time::ZERO, limit * 0.55, &mut bus).unwrap();
        }
        assert!(
            (w.norm() - datum).abs() < 1e-12,
            "P drifted from {datum} to {} after the switch",
            w.norm()
        );
    }

    /// **The discrete eigenvalue, exactly.** A sampled `sin` is an eigenvector of the
    /// discrete H, so its density must not move at all, and its phase must rotate at
    /// `E_disc/ℏ` *through the scheme's own dispersion*: the per-step phase θ satisfies
    /// `sin(θ/2) = E·dt/2ℏ`, so inverting that recovers `E_disc` with the time
    /// discretisation accounted for analytically rather than tolerated.
    ///
    /// The measurement is the three-term recurrence `Rⁿ⁺¹ + Rⁿ⁻¹ = 2cosθ·Rⁿ`, exact for a
    /// rotation. Taken at the first step so rounding has no time to accumulate: the probe
    /// values are exact to ~1e-16 relative, `1 − cosθ ≈ 2e-5` is therefore known to ~1e-11
    /// relative, and E to half that — the 1e-9 asserted is that cancellation budget with a
    /// factor ~20 of margin, not a physics tolerance.
    #[test]
    fn an_eigenstate_is_stationary_and_rotates_at_the_discrete_eigenvalue() {
        let mut w = electron_in(1.0, 49).in_eigenstate(2);
        let dt = w.max_stable_dt(Time::ZERO) * 0.8;
        let e_disc = w.eigenstate_energy(2).to_si();
        let probe = 12usize; // sin(2π·13/50) ≈ 0.998, nowhere near a node of the mode
        let mut bus = Exchange::new();

        let r0 = w.re[probe];
        w.step(Time::ZERO, dt, &mut bus).unwrap();
        let r1 = w.re[probe];
        let density_after_start: Vec<f64> = (1..=w.cells())
            .map(|j| w.node_density(j as isize))
            .collect();
        w.step(Time::ZERO, dt, &mut bus).unwrap();
        let r2 = w.re[probe];

        let cos_theta = (r2 + r0) / (2.0 * r1);
        let e_measured = 2.0 * HBAR / dt.to_si() * ((1.0 - cos_theta) / 2.0).sqrt();
        assert!(
            (e_measured / e_disc - 1.0).abs() < 1e-9,
            "marched {e_measured:e} J against the discrete closed form {e_disc:e} J"
        );

        // And the discrete eigenvalue is *not* the continuum one at this resolution —
        // otherwise the assertion above could not tell them apart.
        let e_cont = 4.0 * PI * PI * HBAR * HBAR / (2.0 * M_E * 1e-18);
        assert!(
            (e_disc / e_cont - 1.0).abs() > 1e-4,
            "50 samples should visibly separate discrete from continuum"
        );

        // Stationarity of |ψ|² under 2000 more steps, to rounding: the eigenstate claim
        // is exactness, so the tolerance is accumulated arithmetic and nothing else.
        for _ in 0..2000 {
            w.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        let peak = density_after_start.iter().fold(0.0f64, |a, r| a.max(*r));
        for (j, before) in density_after_start.iter().enumerate() {
            let now = w.node_density((j + 1) as isize);
            assert!(
                (now - before).abs() < peak * 1e-12,
                "node {}: density moved from {before:e} to {now:e}",
                j + 1
            );
        }
    }

    /// **The convergence rate, which is the check that catches boundary defects.** The
    /// marched eigenvalue must approach `n²π²ℏ²/2mL²` at second order: the closed-form
    /// error is `(π·dx/2L)²/3·E + O(dx⁴)`, so each grid doubling divides it by
    /// `4·(1 − O(dx²))`. Measured through the same recurrence as above, whose noise floor
    /// (~1e-10 relative) is six orders below the smallest error measured, so the ratio
    /// bounds (3.9, 4.1) are limited by the O(dx⁴) term — about 0.06% here — not by
    /// measurement.
    ///
    /// A wall handled wrong is exactly what this notices: the acoustic domain's boundary
    /// defect read as second-order code converging at first order, and no single-grid
    /// tolerance saw it.
    #[test]
    fn the_discrete_eigenvalue_converges_to_the_continuum_at_second_order() {
        let e_continuum = PI * PI * HBAR * HBAR / (2.0 * M_E * 1e-18);
        let measure = |cells: usize| -> f64 {
            let mut w = electron_in(1.0, cells).in_eigenstate(1);
            let dt = w.max_stable_dt(Time::ZERO) * 0.8;
            let probe = cells / 3;
            let mut bus = Exchange::new();
            let r0 = w.re[probe];
            w.step(Time::ZERO, dt, &mut bus).unwrap();
            let r1 = w.re[probe];
            w.step(Time::ZERO, dt, &mut bus).unwrap();
            let r2 = w.re[probe];
            let cos_theta = (r2 + r0) / (2.0 * r1);
            2.0 * HBAR / dt.to_si() * ((1.0 - cos_theta) / 2.0).sqrt()
        };

        // 20, 40 and 80 spacings across the well: dx halves twice.
        let errors: Vec<f64> = [19, 39, 79]
            .iter()
            .map(|&n| 1.0 - measure(n) / e_continuum)
            .collect();
        // The discrete eigenvalue approaches from below: sin²x < x².
        for (i, e) in errors.iter().enumerate() {
            assert!(
                *e > 0.0,
                "grid {i}: the discrete value must sit low, got {e:e}"
            );
        }
        // (π/40)²/3 = 2.06e-3 at the coarsest grid.
        assert!(
            (errors[0] / 2.06e-3 - 1.0).abs() < 0.01,
            "coarse error {:e} against the closed-form prediction 2.06e-3",
            errors[0]
        );
        for pair in errors.windows(2) {
            let fall = pair[0] / pair[1];
            assert!(
                (3.9..4.1).contains(&fall),
                "each doubling must quarter the error: {:e} then {:e}, a factor of {fall:.3}",
                pair[0],
                pair[1]
            );
        }
    }

    /// **Free spreading against the closed form**, `σ(t)² = σ₀²(1 + (ℏt/2mσ₀²)²)`, run to
    /// the moment the variance has exactly doubled, in a well wide enough (walls 10σ out
    /// at the end, truncation ~e⁻⁵⁰) that the boundary is irrelevant.
    ///
    /// The error budget has a closed form of its own. The discrete Laplacian's group
    /// velocity is `(ℏ/m)(k − k³dx²/6)`, so the variance's growth rate is low by
    /// `(σ_k·dx)²` relative and the doubled variance at `t*` is low by half that:
    /// `(σ_k·dx)²/2 = 1.953e-3` at 159 cells, with `σ_k = 1/2σ₀`. Measured 1.946e-3 —
    /// the coefficient, not merely the order — and 4.878e-4 at 319 cells, a fall of 3.99.
    /// The time term enters at `(E·dt/2ℏ)²/6 ≈ 1e-7` and, tied to dx² through the limit,
    /// falls sixteenfold per doubling, so it cannot blur the ratio; the moment sums over
    /// a smooth Gaussian are spectrally accurate and contribute nothing visible.
    #[test]
    fn a_free_gaussian_spreads_at_the_closed_form_rate() {
        let sigma0 = 0.4e-9;
        let t_star = 2.0 * M_E * sigma0 * sigma0 / HBAR; // σ² doubles here: 2.76 fs
        let spread_error = |cells: usize| -> f64 {
            let mut w = electron_in(8.0, cells).with_gaussian(
                Length::nm(4.0),
                Length::from_si(sigma0),
                Wavenumber::from_si(0.0),
            );
            let limit = w.max_stable_dt(Time::ZERO).to_si();
            let steps = (t_star / (0.9 * limit)).ceil() as u32;
            let dt = Time::from_si(t_star / steps as f64);
            let mut bus = Exchange::new();
            for _ in 0..steps {
                w.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            let (_, mean, var) = moments(&w);
            // The packet must not have drifted: no k₀, no asymmetry.
            assert!((mean - 4.0e-9).abs() < 1e-12, "⟨x⟩ moved to {mean:e}");
            (var / (2.0 * sigma0 * sigma0) - 1.0).abs()
        };

        let coarse = spread_error(159);
        let fine = spread_error(319);
        // Two-sided: the error must *be* the dispersion budget, not merely under a lid a
        // wrong coefficient could also fit beneath. The 5% window is the (σ_k·dx)⁴ tail
        // and the packet's k⁶ moments, both measured well inside it.
        assert!(
            (coarse / 1.953e-3 - 1.0).abs() < 0.05,
            "the coarse error should be (σ_k·dx)²/2 = 1.953e-3, got {coarse:e}"
        );
        let fall = coarse / fine;
        assert!(
            (3.8..4.2).contains(&fall),
            "halving dx must quarter the error: {coarse:e} then {fine:e}, a factor of {fall:.2}"
        );
    }

    /// **Ehrenfest on the harmonic potential.** For quadratic V the theorem is exact in
    /// the continuum: `⟨x⟩` of a displaced ground state — a coherent state — follows
    /// `x₀·cos(ωt)` with no quantum correction at all, so every deviation is
    /// discretisation and must fall at second order.
    ///
    /// Measured at the quarter period, where the cosine crosses zero and a frequency error
    /// shows up at first order in the phase — at the period's end the crossing error is
    /// quadratic in the phase error and a wrong frequency could hide. The budget is the
    /// Laplacian dispersion at the packet's peak momentum `k = mωd/ℏ`: `(k·dx)²/12`
    /// ≈ 1.5e-3 at 199 cells, times an O(1) factor for the packet's spread about that
    /// peak. Measured: 2.840e-3 of the amplitude at 199 cells, 7.160e-4 at 399 — a fall
    /// of 3.97, which is the second order the budget predicts.
    #[test]
    fn ehrenfest_holds_in_a_harmonic_potential() {
        let omega = 1.570_8e15; // period T = 4 fs
        let period = 2.0 * PI / omega;
        // The harmonic ground-state width, so the displaced packet is a coherent state
        // and its shape carries no breathing of its own.
        let sigma = (HBAR / (2.0 * M_E * omega)).sqrt(); // 0.192 nm
        let displacement = 0.5e-9;
        let centre = 2.0e-9;

        let crossing_error = |cells: usize| -> f64 {
            let mut w = electron_in(4.0, cells)
                .with_harmonic(Frequency::from_si(omega))
                .with_gaussian(
                    Length::from_si(centre + displacement),
                    Length::from_si(sigma),
                    Wavenumber::from_si(0.0),
                );
            let limit = w.max_stable_dt(Time::ZERO).to_si();
            // A step count divisible by four, so the quarter and half periods land on
            // steps exactly.
            let steps = 4 * (period / (4.0 * 0.9 * limit)).ceil() as u32;
            let dt = Time::from_si(period / steps as f64);
            let mut bus = Exchange::new();

            for _ in 0..steps / 4 {
                w.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            let quarter = (w.position_expectation().to_si() - centre).abs() / displacement;

            for _ in 0..steps / 4 {
                w.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            let half = w.position_expectation().to_si();
            assert!(
                (half - (centre - displacement)).abs() < 0.01 * displacement,
                "half a period should invert the displacement: ⟨x⟩ = {half:e}"
            );
            quarter
        };

        let coarse = crossing_error(199);
        let fine = crossing_error(399);
        assert!(
            coarse < 5e-3,
            "199 cells should hold a few 1e-3 of the amplitude, got {coarse:e}"
        );
        let fall = coarse / fine;
        assert!(
            (3.7..4.3).contains(&fall),
            "halving dx must quarter the crossing error: {coarse:e} then {fine:e}, \
             a factor of {fall:.2}"
        );
    }

    /// **The limit is derived, enforced, and real.** Its formula is checked against the
    /// bare-well value `m·dx²/ℏ`, a potential tightens it by exactly its `V_max`, a step
    /// past it is refused by name — and the raw scheme, driven past the limit through the
    /// private update the public API refuses to expose, actually diverges, while just
    /// under it the most dangerous state in the well holds its books to rounding.
    #[test]
    fn the_stability_limit_is_derived_and_real() {
        let cells = 60usize;
        let free = electron_in(1.0, cells);
        let dx = free.spacing().to_si();
        let limit = free.max_stable_dt(Time::ZERO).to_si();
        assert!(
            (limit / (M_E * dx * dx / HBAR) - 1.0).abs() < 1e-12,
            "the bare-well limit is m·dx²/ℏ, got {limit:e}"
        );

        // A potential tightens it by exactly its maximum over the nodes.
        let omega = 3.0e15;
        let stiff = electron_in(1.0, cells).with_harmonic(Frequency::from_si(omega));
        let v_max = stiff.potential.iter().fold(0.0f64, |a, v| a.max(*v));
        assert!(v_max > 0.0);
        let expected = 2.0 * HBAR / (2.0 * HBAR * HBAR / (M_E * dx * dx) + v_max);
        let tightened = stiff.max_stable_dt(Time::ZERO).to_si();
        assert!((tightened / expected - 1.0).abs() < 1e-12);
        assert!(tightened < limit);

        // Past the limit the step is refused, by name and by how much.
        let mut w = electron_in(1.0, cells).in_eigenstate(cells);
        let mut bus = Exchange::new();
        let err = w
            .step(Time::ZERO, Time::from_si(limit * 1.01), &mut bus)
            .expect_err("past the limit must be refused");
        assert_eq!(err.quantity, "stability ratio");
        assert!(
            (err.after - 1.01).abs() < 1e-9,
            "by how much: {}",
            err.after
        );

        // The top eigenstate at near-limit dt rotates almost π per step, so its startup
        // adjustment is sin²(θ/2) ≈ 96% of the total — and the public path *refuses* to
        // absorb that into the datum, because a step that cannot resolve a state's phase
        // should be reported, not accommodated.
        let mut w = electron_in(1.0, cells).in_eigenstate(cells);
        let err = w
            .step(Time::ZERO, Time::from_si(limit * 0.98), &mut bus)
            .expect_err("an unresolved startup must be refused");
        assert_eq!(err.quantity, PROBABILITY);
        assert!(err.site.contains("starting"), "{err}");

        // So the sharpness of the limit is shown on the raw update, same path as the
        // divergence below: the top of the discrete band — the mode the limit is *about* —
        // marches 3000 steps just under it with the invariant held to accumulated
        // rounding and the field bounded.
        let mut w = electron_in(1.0, cells).in_eigenstate(cells);
        let h = limit * 0.98;
        w.advance(h);
        let datum = w.invariant();
        for _ in 0..3000 {
            w.advance(h);
        }
        assert!(
            (w.invariant() - datum).abs() < 1e-11,
            "just under the limit the invariant holds: {datum} became {}",
            w.invariant()
        );
        let bounded = w.re.iter().fold(0.0f64, |a, r| a.max(r.abs()));

        // Just over it — reachable only through the private update, because `step`
        // refuses — the same state diverges. 2% over makes the top mode's half-angle
        // 1.02·cos²(π/122) ≈ 1.019 > 1: hyperbolic, ~1.5× per step.
        let mut w = electron_in(1.0, cells).in_eigenstate(cells);
        let start_peak = w.re.iter().fold(0.0f64, |a, r| a.max(r.abs()));
        let h = limit * 1.02;
        w.advance(h);
        let datum = w.invariant();
        for _ in 0..150 {
            w.advance(h);
        }
        let peak = w.re.iter().fold(0.0f64, |a, r| a.max(r.abs()));
        assert!(
            peak > 1e6 * start_peak,
            "2% past the limit must diverge: peak went from {start_peak:e} to {peak:e} \
             while the stable run stayed at {bounded:e}"
        );

        // And the invariant *still* balances — relative to the entries' scale, which has
        // exploded with the field. This is why the audit cannot stand in for the limit:
        // a diverging run keeps its books, on books nobody should accept, and it is the
        // reason `Ledger` records a scale beside every total.
        let gross: f64 = (1..=cells)
            .map(|j| {
                let k = j - 1;
                (w.re[k] * w.re[k]).abs() + (w.im[k] * w.im_prev[k]).abs()
            })
            .sum::<f64>()
            * dx;
        assert!(
            (w.invariant() - datum).abs() < 1e-9 * gross,
            "the identity holds even while diverging: {datum:e} against {:e} on a scale \
             of {gross:e}",
            w.invariant()
        );
    }

    /// **⟨H⟩ is exact on an eigenstate and does not drift on anything** — the check the
    /// task of keeping ⟨H⟩ *off* the ledger hands to the tests instead.
    ///
    /// The staggered estimator `Σ[R·HR + I⁺·H I⁻]·dx / P` is diagonal in the eigenbasis
    /// and constant mode by mode, so on a general state it is conserved to rounding; on an
    /// eigenstate it equals the discrete eigenvalue exactly, before the march and after
    /// it. Both tolerances are accumulated arithmetic, not physics.
    #[test]
    fn the_energy_expectation_is_exact_and_does_not_drift() {
        // Exactness on an eigenstate, before any step and after many.
        let mut w = electron_in(1.0, 49).in_eigenstate(3);
        let e_disc = w.eigenstate_energy(3).to_si();
        assert!(
            (w.energy_expectation().to_si() / e_disc - 1.0).abs() < 1e-12,
            "⟨H⟩ at rest: {:e} against {e_disc:e}",
            w.energy_expectation().to_si()
        );
        let dt = w.max_stable_dt(Time::ZERO) * 0.8;
        let mut bus = Exchange::new();
        for _ in 0..500 {
            w.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        assert!(
            (w.energy_expectation().to_si() / e_disc - 1.0).abs() < 1e-12,
            "⟨H⟩ after 500 steps: {:e} against {e_disc:e}",
            w.energy_expectation().to_si()
        );

        // No drift on a packet in a potential, where nothing is stationary.
        let mut w = electron_in(2.0, 101)
            .with_harmonic(Frequency::from_si(1.5e15))
            .with_gaussian(Length::nm(1.3), Length::nm(0.2), Wavenumber::from_si(2.0e9));
        let dt = w.max_stable_dt(Time::ZERO) * 0.9;
        w.step(Time::ZERO, dt, &mut bus).unwrap();
        let datum = w.energy_expectation().to_si();
        for _ in 0..3000 {
            w.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        let drift = (w.energy_expectation().to_si() / datum - 1.0).abs();
        assert!(drift < 1e-11, "⟨H⟩ drifted by {drift:e} over 3000 steps");
    }

    /// The rate the field reports is exactly what the step just taken did to the density —
    /// the same honesty `Room` pins for pressure, and the same half-step-stale trade.
    #[test]
    fn the_reported_rate_is_exactly_the_step_just_taken() {
        let mut w = electron_in(2.0, 61).with_gaussian(
            Length::nm(0.8),
            Length::nm(0.2),
            Wavenumber::from_si(6.0e9),
        );
        let dt = w.max_stable_dt(Time::ZERO) * 0.8;
        let mut bus = Exchange::new();

        // One priming step: the continuity form is exact for a step taken from a
        // *staggered* state, and across the startup step the observed change also carries
        // the one-time half-kick of I — the same O(dt²) the startup adjustment declares.
        // The first version of this test measured across that step and read a 6% mismatch
        // at the first node, which is that startup and not a defect in the rate.
        w.step(Time::ZERO, dt, &mut bus).unwrap();

        // Twice, so the check cannot be an accident of one particular state.
        for round in 0..2 {
            let before: Vec<f64> = (1..=w.cells())
                .map(|j| w.node_density(j as isize))
                .collect();
            w.step(Time::ZERO, dt, &mut bus).unwrap();
            let reported: Vec<f64> = (1..=w.cells())
                .map(|j| {
                    let x = j as f64 * w.dx;
                    w.rate(LengthVec::m(x, 0.0, 0.0), Time::ZERO, dt)
                })
                .collect();
            let scale = reported.iter().fold(0.0f64, |a, v| a.max(v.abs()));
            assert!(scale > 0.0, "round {round}: nothing moved, nothing tested");
            for j in 1..=w.cells() {
                let observed = (w.node_density(j as isize) - before[j - 1]) / dt.to_si();
                assert!(
                    (observed - reported[j - 1]).abs() < scale * 1e-12,
                    "round {round}, node {j}: reported {:e} but the step did {observed:e}",
                    reported[j - 1]
                );
            }
        }
    }

    /// The three readings and the field are wired, and honest about the walls.
    ///
    /// The opt-ins are the whole of how the scene and view layers can see a domain, and a
    /// domain that forgets one is silently absent rather than broken — which is why this
    /// is a test and not an assumption.
    #[test]
    fn readings_and_the_field_are_wired() {
        let w = electron_in(1.0, 49).in_eigenstate(1);
        let l = w.width().to_si();

        let readings = w.readings();
        let labels: Vec<&str> = readings.iter().map(|r| r.label.as_str()).collect();
        assert_eq!(labels, vec!["norm", "<x>", "<E>"]);
        assert!((readings[0].value - 1.0).abs() < 1e-12, "norm reads 1");
        assert_eq!(readings[1].unit, "m");
        assert_eq!(readings[2].unit, "J");
        // An eigenstate is symmetric, so ⟨x⟩ is the centre to rounding.
        assert!((readings[1].value / (l / 2.0) - 1.0).abs() < 1e-12);
        assert!(
            (readings[2].value / w.eigenstate_energy(1).to_si() - 1.0).abs() < 1e-12,
            "⟨E⟩ reads the discrete eigenvalue"
        );

        // The field: reachable through the trait, labelled, zero on and beyond the walls.
        let field = w.as_field().expect("the well is a field");
        assert_eq!(field.unit(), "1/m");
        assert_eq!(field.at(LengthVec::m(0.0, 0.0, 0.0), Time::ZERO), 0.0);
        assert_eq!(field.at(LengthVec::m(l, 0.0, 0.0), Time::ZERO), 0.0);
        assert_eq!(field.at(LengthVec::m(-1.0, 0.0, 0.0), Time::ZERO), 0.0);
        assert_eq!(field.at(LengthVec::m(2.0 * l, 0.0, 0.0), Time::ZERO), 0.0);
        assert!(field.at(LengthVec::m(f64::NAN, 0.0, 0.0), Time::ZERO) == 0.0);
        // The ground state peaks at the middle, ~2/L for sin² normalised over L.
        let mid = field.at(LengthVec::m(l / 2.0, 0.0, 0.0), Time::ZERO);
        assert!(
            (mid / (2.0 / l) - 1.0).abs() < 1e-3,
            "midwell density {mid:e} against 2/L = {:e}",
            2.0 / l
        );
        // Interpolation is linear between nodes, so a midpoint is the average of its ends.
        let dx = w.spacing().to_si();
        let a = field.at(LengthVec::m(3.0 * dx, 0.0, 0.0), Time::ZERO);
        let b = field.at(LengthVec::m(4.0 * dx, 0.0, 0.0), Time::ZERO);
        let between = field.at(LengthVec::m(3.5 * dx, 0.0, 0.0), Time::ZERO);
        assert!(((a + b) / 2.0 - between).abs() < 1e-12 * (a + b));
        // Uniform off-axis: the one-dimensional claim.
        assert_eq!(
            field.at(LengthVec::m(l / 3.0, 0.0, 0.0), Time::ZERO),
            field.at(LengthVec::m(l / 3.0, 55.0, -9.0), Time::ZERO)
        );
        // The Laplacian of sin² has a closed form: ρ = (2/L)sin²(πx/L) gives
        // ρ'' = (2/L)·2(π/L)²cos(2πx/L), which at the centre is −4π²/L³·... : just check
        // the sign structure — negative curvature at the peak, positive at the walls.
        assert!(field.laplacian(LengthVec::m(l / 2.0, 0.0, 0.0), Time::ZERO, w.spacing()) < 0.0);
        assert!(field.laplacian(LengthVec::m(dx, 0.0, 0.0), Time::ZERO, w.spacing()) > 0.0);
    }

    /// The eleventh domain under the scheduler: subcycled to its own limit, audited by the
    /// kernel — on a quantity the kernel has never heard of — with its books balancing
    /// per-domain because it claims they do.
    #[test]
    fn the_domain_runs_under_the_scheduler() {
        use dualis_core::{Schedule, Simulation};

        let w = electron_in(1.0, 49).in_eigenstate(1);
        let limit = w.max_stable_dt(Time::ZERO).to_si();
        // 49 interior nodes across 1 nm: 20 pm cells, and m·dx²/ℏ = 3.46 as.
        assert!((limit / 3.455e-18 - 1.0).abs() < 1e-3, "limit {limit:e} s");

        let outer = Time::from_si(2.0e-17);
        let expected = (outer.to_si() / limit).ceil() as u32;
        let mut sim = Simulation::new(Schedule::Multirate)
            // Tight, because the ledger carries the scheme-exact invariant: the identity
            // leaves nothing for a tolerance to forgive.
            .conservation_tolerance(1e-9)
            .conservation_tolerance_for(PROBABILITY, 1e-12)
            .with(w);

        let report = sim.advance(outer).unwrap();
        assert_eq!(report.substeps[0].0, "well");
        assert_eq!(report.substeps[0].1, expected);
        for _ in 0..20 {
            sim.advance(outer).expect("probability should hold");
        }
        assert!((sim.ledger().get(PROBABILITY).unwrap() - 1.0).abs() < 1e-9);

        // The domain is reachable back out through every opt-in.
        assert!(sim.domain_as::<Well>("well").is_some());
        assert!(sim.field("well").is_some());
        assert_eq!(sim.domain("well").unwrap().readings().len(), 3);
    }

    /// Restore puts back **every** piece of state, the staggering flag included, and the
    /// whole march is bit-for-bit deterministic.
    ///
    /// The flag is the trap: a restore that forgot it would re-run the leapfrog startup —
    /// a second half-kick on a state that already had one — and nothing else would notice.
    #[test]
    fn restore_puts_back_every_piece_of_state_and_the_march_is_deterministic() {
        let build = || {
            electron_in(2.0, 61)
                .with_harmonic(Frequency::from_si(1.2e15))
                .with_gaussian(Length::nm(1.1), Length::nm(0.2), Wavenumber::from_si(3.0e9))
        };
        let bits = |w: &Well| -> Vec<u64> {
            w.re.iter()
                .chain(w.im.iter())
                .chain(w.im_prev.iter())
                .map(|v| v.to_bits())
                .collect()
        };
        let mut bus = Exchange::new();

        // Checkpoint before the very first step, so the flag matters.
        let mut w = build();
        let dt = w.max_stable_dt(Time::ZERO) * 0.9;
        w.checkpoint();
        w.step(Time::ZERO, dt, &mut bus).unwrap();
        let first = bits(&w);
        w.restore();
        assert!(!w.staggered, "restore must rewind the staggering flag");
        w.step(Time::ZERO, dt, &mut bus).unwrap();
        assert_eq!(first, bits(&w), "a restored start must replay exactly");

        // Checkpoint mid-flight, run on, rewind, run again: identical.
        for _ in 0..100 {
            w.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        w.checkpoint();
        for _ in 0..50 {
            w.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        let onward = bits(&w);
        w.restore();
        for _ in 0..50 {
            w.step(Time::ZERO, dt, &mut bus).unwrap();
        }
        assert_eq!(onward, bits(&w));

        // And two independent runs agree to the bit — no clock, no randomness.
        let run = || {
            let mut w = build();
            let dt = w.max_stable_dt(Time::ZERO) * 0.9;
            let mut bus = Exchange::new();
            for _ in 0..500 {
                w.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            bits(&w)
        };
        assert_eq!(run(), run(), "not bit-identical across runs");
    }

    /// Degenerate wells are handled rather than being cases a caller must avoid.
    #[test]
    fn degenerate_wells_are_handled() {
        // Fewer than three cells is not a grid; the constructor rounds up.
        assert_eq!(electron_in(1.0, 1).cells(), 3);

        // An empty well steps, holds a norm of zero, and reports its centre.
        let mut empty = electron_in(1.0, 21);
        let mut bus = Exchange::new();
        let dt = empty.max_stable_dt(Time::ZERO) * 0.5;
        empty.step(Time::ZERO, dt, &mut bus).unwrap();
        assert_eq!(empty.norm(), 0.0);
        assert_eq!(empty.ledger().get(PROBABILITY), Some(0.0));
        assert!((empty.position_expectation().to_si() / 0.5e-9 - 1.0).abs() < 1e-12);
        assert_eq!(empty.energy_expectation().to_si(), 0.0);

        // A zero step changes nothing.
        let mut w = electron_in(1.0, 21).in_eigenstate(1);
        let before: Vec<u64> = w.re.iter().map(|v| v.to_bits()).collect();
        w.step(Time::ZERO, Time::ZERO, &mut bus).unwrap();
        let after: Vec<u64> = w.re.iter().map(|v| v.to_bits()).collect();
        assert_eq!(before, after);

        // A massless particle has no stable step at all, and is refused rather than
        // integrated: infinitely fast phase rotation is not a state this scheme can march.
        let mut massless = Well::new("m0", Mass::kg(0.0), Length::nm(1.0), 21);
        assert_eq!(massless.max_stable_dt(Time::ZERO).to_si(), 0.0);
        assert!(massless
            .step(Time::ZERO, Time::from_si(1e-20), &mut bus)
            .is_err());

        // Eigenstate indices clamp to what the grid can carry.
        let ground = electron_in(1.0, 21).in_eigenstate(0);
        assert!(
            (ground.energy_expectation().to_si() / ground.eigenstate_energy(1).to_si() - 1.0).abs()
                < 1e-12,
            "n = 0 clamps to the ground state"
        );
        let top = electron_in(1.0, 21).in_eigenstate(99);
        assert!(
            (top.energy_expectation().to_si() / top.eigenstate_energy(21).to_si() - 1.0).abs()
                < 1e-12,
            "n past the band clamps to the top"
        );

        // A packet narrower than a cell widens to one rather than aliasing, and a clipped
        // packet still normalises over the well as it is.
        let narrow = electron_in(1.0, 21).with_gaussian(
            Length::nm(0.5),
            Length::from_si(0.0),
            Wavenumber::from_si(0.0),
        );
        assert!((narrow.norm() - 1.0).abs() < 1e-12);
        let clipped = electron_in(1.0, 21).with_gaussian(
            Length::nm(0.05),
            Length::nm(0.2),
            Wavenumber::from_si(0.0),
        );
        assert!((clipped.norm() - 1.0).abs() < 1e-12);
    }
}
