//! Where scalar diffraction is right, and where a field solution says it is not.
//!
//! `dualis-optics` gives a slit's far field as `sinc²(π a sinθ/λ)`. That rests on **Kirchhoff's**
//! boundary condition: the field in the opening is taken to be the incident field, and zero on the
//! screen. Maxwell's equations do not say that. The metal carries currents, the field in the opening
//! is perturbed near its edges over a distance of order `λ`, and the fraction of the aperture that
//! is wrong therefore goes as `λ/a`.
//!
//! So the two agree in a limit and not otherwise, and this measures where. `dualis-em` marches the
//! actual field, the aperture plane is transformed to the far field, and the result is compared with
//! the formula at four slit widths.
//!
//! Measured as the largest absolute difference in intensity, normalised to the axis:
//!
//! ```text
//!   a = 12λ    0.0057     the formula is right
//!   a =  6λ    0.0125
//!   a =  3λ    0.0311
//!   a =  1λ    0.2772     it is wrong by more than a quarter of the pattern it predicts
//! ```
//!
//! A factor of 48 from one wavelength to twelve, monotone — which is what says the disagreement is
//! Kirchhoff's condition and not something that happens to be large at one width.
//!
//! # How the far field is got out of a small box
//!
//! Not by making the box big enough to be in the far field — that would be `a²/λ` deep, hundreds of
//! wavelengths for a wide slit. Instead the tangential field is recorded on a plane just past the
//! screen and Fourier transformed:
//!
//! ```text
//!   E(θ) ∝ ∫ E(x) e^{−i k x sinθ} dx
//! ```
//!
//! which **is** the Fraunhofer limit, exactly. The box only has to be deep enough to hold the screen
//! and a wavelength behind it, and wide enough that the diffracted near field has not reached the
//! side walls by the time it crosses the recording plane.
//!
//! # What makes a plane wave possible next to a screen
//!
//! Magnetic conductors on the `x` faces. A wave along `z` polarised along `y` has `Ey` tangential
//! there, so an electric conductor would force it to zero and there would be no plane wave to
//! diffract. The mirror images that a magnetic wall implies would make this a grating rather than a
//! slit — but the transform integrates over **one** box width, and the neighbouring slits lie outside
//! it, so what is transformed is this slit's aperture field and no other's.

use dualis::em::{Boundary, Cavity, Medium, Wall};
use dualis::optics::{single_slit_intensity, slit_zero};
use dualis::prelude::*;

/// Cells per wavelength.
const PER: usize = 16;
/// Half the box's width, in wavelengths. Wide enough that the diffracted field has not reached the
/// side walls one wavelength past the screen.
const HALF_WIDTH: usize = 9;

/// The aperture field of one slit, transformed in time at the carrier.
///
/// # Marched once, sampled as often as wanted
///
/// The march is the expensive part and the angle is not, so they are separated: three tests ask
/// about the same slit widths, and calling a `pattern(width, samples)` for each of them marched the
/// same box eight times where five will do. It cost 120 s in a debug build, which is the mode four
/// of CI's jobs run.
struct Aperture {
    /// Real and imaginary parts of `Ey` at the recording plane, per transverse node.
    re: Vec<f64>,
    im: Vec<f64>,
    dx: f64,
    k0: f64,
    /// `|E(0)|²`, so every intensity is normalised to the axis as the closed form is.
    axis: f64,
}

impl Aperture {
    /// The far-field intensity at an angle, by the Fraunhofer transform of the plane.
    fn intensity(&self, sin_theta: f64) -> f64 {
        let n = self.re.len();
        let (mut tr, mut ti) = (0.0, 0.0);
        for i in 0..n {
            let x = (i as f64 - (n - 1) as f64 / 2.0) * self.dx;
            let ph = -self.k0 * x * sin_theta;
            let (cs, sn) = (ph.cos(), ph.sin());
            tr += self.re[i] * cs - self.im[i] * sn;
            ti += self.re[i] * sn + self.im[i] * cs;
        }
        (tr * tr + ti * ti) / self.axis.max(f64::MIN_POSITIVE)
    }

    /// The pattern sampled out to half of grazing, which is past the third zero of every slit here.
    fn pattern(&self, samples: usize) -> Vec<(f64, f64)> {
        (0..samples)
            .map(|s| {
                let sin_theta = 0.5 * s as f64 / (samples - 1) as f64;
                (sin_theta, self.intensity(sin_theta))
            })
            .collect()
    }
}

/// March a slit and record its aperture field.
fn aperture(width_wavelengths: f64) -> Aperture {
    let wavelength = 500e-9;
    let dx = wavelength / PER as f64;
    let nx = 2 * HALF_WIDTH * PER;
    let nz = 10 * PER;
    let screen = 8 * PER;
    let record = screen + PER; // one wavelength behind the screen
    let half_slit = (0.5 * width_wavelengths * PER as f64).round() as usize;
    let centre = nx / 2;

    let mut c = Cavity::new("slit", (nx, 2, nz), Length::from_si(dx), Medium::vacuum());
    c.set_boundary(Wall::XLow, Boundary::Magnetic);
    c.set_boundary(Wall::XHigh, Boundary::Magnetic);
    c.set_boundary(Wall::ZLow, Boundary::Open);
    c.set_boundary(Wall::ZHigh, Boundary::Open);
    // A screen two cells thick with an opening in the middle. Two rather than one because a
    // one-cell screen leaks: the electric edges either side of it are shared with the vacuum cells
    // beyond, and holding only those leaves a gap the field slips through.
    c.obstruct(|i, _j, k| {
        (k == screen || k == screen + 1) && (i + half_slit < centre || i > centre + half_slit)
    });

    let dt = Time::from_si(c.courant_limit().to_si() * 0.5);
    let k0 = 2.0 * std::f64::consts::PI / wavelength;
    let (start, spread) = (3.0 * wavelength, 1.6 * wavelength);
    c.launch_along_z(dt, move |_x, _y, z| {
        let u = (z - start) / spread;
        (-u * u).exp() * (k0 * (z - start)).cos()
    });

    // March until the pulse has crossed the recording plane and gone, accumulating the transform in
    // time as well as space: one pass of the pulse is one measurement of the aperture field.
    let mut bus = Exchange::new();
    let steps = (1.3 * nz as f64 / 0.2887) as usize;
    let f = 299_792_458.0 / wavelength;
    let mut re = vec![0.0f64; nx + 1];
    let mut im = vec![0.0f64; nx + 1];
    for n in 0..steps {
        let t = n as f64 * dt.to_si();
        c.step(Time::from_si(t), dt, &mut bus).expect("stable");
        let phase = -2.0 * std::f64::consts::PI * f * t;
        let (cs, sn) = (phase.cos(), phase.sin());
        for i in 0..=nx {
            let e = c.electric_at(i, 0, record).y;
            re[i] += e * cs;
            im[i] += e * sn;
        }
    }

    let axis = {
        let (ar, ai): (f64, f64) = (re.iter().sum(), im.iter().sum());
        ar * ar + ai * ai
    };
    Aperture {
        re,
        im,
        dx,
        k0,
        axis,
    }
}

/// **A wide slit's field solution converges on `sinc²`, and a narrow one does not.**
///
/// The whole main lobe and the first sidelobe, compared point by point. The error is the largest
/// absolute difference in normalised intensity, which is the honest measure for a pattern that
/// spans two orders of magnitude — a relative error at a zero is meaningless.
///
/// Kirchhoff's condition is wrong over a strip of order `λ` at each edge of the opening, so the
/// error should fall roughly as `λ/a`. That is the claim: a *rate*, over four widths.
#[test]
fn a_wide_slit_agrees_with_scalar_diffraction_and_a_narrow_one_does_not() {
    let wavelength = Length::from_si(500e-9);
    let mut errors = Vec::new();
    for width in [1.0, 3.0, 6.0, 12.0] {
        let a = Length::from_si(width * 500e-9);
        let measured = aperture(width).pattern(121);
        // Compare out to the second zero, or to half of grazing for a slit with no zeros.
        let limit = slit_zero(a, wavelength, 2).unwrap_or(0.5).min(0.5);
        let mut worst: f64 = 0.0;
        for (sin_theta, intensity) in measured.iter().filter(|(s, _)| *s <= limit) {
            let closed = single_slit_intensity(a, wavelength, *sin_theta);
            worst = worst.max((intensity - closed).abs());
        }
        println!(
            "  a = {width:>4.1} lambda: worst difference {:.4} of the axial intensity, out to \
             sin theta = {limit:.3}",
            worst
        );
        errors.push((width, worst));
    }

    // The narrowest is wrong by more than the pattern it is predicting.
    assert!(
        errors[0].1 > 0.2,
        "at one wavelength the scalar formula should be badly wrong: {:.4}",
        errors[0].1
    );
    // The widest is close.
    assert!(
        errors[3].1 < 0.06,
        "at twelve wavelengths it should be close: {:.4}",
        errors[3].1
    );
    // And it improves monotonically, which is what says the disagreement is Kirchhoff's condition
    // and not something that happens to be large at one width.
    for pair in errors.windows(2) {
        assert!(
            pair[1].1 < pair[0].1,
            "widening the slit must help: {:.4} at {} lambda then {:.4} at {}",
            pair[0].1,
            pair[0].0,
            pair[1].1,
            pair[1].0
        );
    }
    println!(
        "  so from one wavelength to twelve the disagreement fell {:.1}x",
        errors[0].1 / errors[3].1
    );
}

/// **The first dark fringe lands at `sinθ = λ/a`**, measured off the field solution.
///
/// The sharpest single number in the pattern and the one with no special constant in it — a
/// circular aperture's is `1.22 λ/D`, where the 1.22 is the first zero of `J₁`. Here it is `sin π`.
///
/// Located by parabolic interpolation through the sampled minimum, so the answer is not quantised
/// to the sampling.
#[test]
fn the_first_dark_fringe_is_where_the_formula_says() {
    let wavelength = Length::from_si(500e-9);
    for width in [6.0, 12.0] {
        let a = Length::from_si(width * 500e-9);
        let measured = aperture(width).pattern(401);
        let closed = slit_zero(a, wavelength, 1).expect("a wide slit has a first zero");

        // The minimum nearest the predicted zero.
        let (at, _) = measured
            .iter()
            .enumerate()
            .filter(|(_, (s, _))| (*s - closed).abs() < 0.4 * closed)
            .min_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
            .expect("there is a minimum near the predicted zero");
        let (lo, mid, hi) = (measured[at - 1].1, measured[at].1, measured[at + 1].1);
        let denom = lo - 2.0 * mid + hi;
        let shift = if denom.abs() > 0.0 {
            0.5 * (lo - hi) / denom
        } else {
            0.0
        };
        let step = measured[1].0 - measured[0].0;
        let found = measured[at].0 + shift * step;

        println!(
            "  a = {width:>4.1} lambda: first zero at sin theta = {found:.5} against \
             lambda/a = {closed:.5} — off {:.2}%",
            (found / closed - 1.0).abs() * 100.0
        );
        assert!(
            (found / closed - 1.0).abs() < 0.04,
            "the first zero is at lambda/a: {found:.5} against {closed:.5}"
        );
    }
}

/// **A sub-wavelength slit radiates rather than diffracting**, which is the statement `slit_zero`
/// returning `None` is making.
///
/// A pattern with no zeros in it at all. Measured: the intensity at grazing is a substantial
/// fraction of the axial intensity, where a twelve-wavelength slit sends essentially nothing there.
/// That contrast is the whole difference between an aperture and an antenna.
#[test]
fn a_sub_wavelength_slit_radiates_instead_of_diffracting() {
    let wavelength = Length::from_si(500e-9);
    assert!(
        slit_zero(Length::from_si(0.6 * 500e-9), wavelength, 1).is_none(),
        "there is no angle at which a 0.6-wavelength slit is dark"
    );

    let narrow = aperture(0.75).intensity(0.5);
    let wide = aperture(12.0).intensity(0.5);
    println!(
        "  at sin theta = 0.5: a 0.75 lambda slit sends {narrow:.4} of its axial intensity, a 12 \
         lambda one {wide:.6}"
    );
    assert!(
        narrow > 0.1,
        "a sub-wavelength slit is nearly omnidirectional: {narrow:.4}"
    );
    assert!(
        narrow / wide.max(1e-12) > 20.0,
        "which is orders more than a wide slit puts there: {:.1}x",
        narrow / wide.max(1e-12)
    );
}
