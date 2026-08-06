//! Where two domains meet, and how a quantity crosses with its place intact.
//!
//! [`Exchange`](crate::sim::Exchange) carries one number per channel per step. Four
//! domains were built on that and none of them minded, but the reason is not that the
//! design was sufficient — it is that none of them ever had to share a *place*. A
//! dichroic absorbing 96 mW hands the whole lot to one lumped mass, because there is no
//! way to say which part of it got hot.
//!
//! That is the physics an instrument simulator actually needs. A coating heats where the
//! beam lands, the heat spreads through the glass, and the temperature field it leaves
//! changes the refractive index. Every step of that needs optics and heat to agree about
//! a surface rather than about a total.
//!
//! # One discretisation, shared, rather than two and an interpolation
//!
//! The temptation is to let each domain keep its own mesh and interpolate between them.
//! That is where energy goes missing: a resampling that is not conservative loses or
//! invents some, and it does so quietly, in a step that looks like bookkeeping.
//!
//! So an [`Interface`] is a boundary cut into faces that **both sides address**. A
//! publisher offers a [`Flux`] over those faces and a consumer takes it over the same
//! ones. A face-count mismatch is refused rather than papered over, and if a caller
//! genuinely needs to cross discretisations they say so with [`Flux::resample`], which
//! conserves the total by construction rather than by hoping.
//!
//! # What this makes auditable that was not
//!
//! The first design pass said the interface between two discretisations is exactly where
//! conservation breaks, and then built a coupling whose interface was a single number —
//! so the check it argued for could not be written. Now it can:
//! [`Exchange::audit_transfers`](crate::sim::Exchange::audit_transfers) names the *face*
//! that was left holding something, not just the channel.

use dualis_units::Area;

use crate::conserved::Violation;

/// A boundary two domains share, cut into faces they both address.
///
/// Faces carry their own areas, because a real boundary is not evenly divided — a
/// spherical cap cut into rings has a smaller innermost one, and a flux per unit area
/// means nothing without them.
///
/// # What it is not
///
/// An ordered sequence of faces with areas, and nothing more. No coordinates, no normals,
/// no connectivity — so a domain cannot ask where face 12 *is*, only how big it is and what
/// comes before it. A beam profile is therefore handed over in the boundary's own
/// coordinate (see [`Flux::profiled`]) rather than computed from geometry.
///
/// That order is what [`Flux::resample`] walks, which makes remapping an interval
/// intersection and keeps it conservative in a few lines. It also means a triangulated
/// surface does not fit: its faces have no sequence, and the overlaps between two
/// triangulations are not intervals. The conservation argument generalises to that case;
/// this implementation does not.
#[derive(Clone, Debug, PartialEq)]
pub struct Interface {
    name: &'static str,
    areas: Vec<f64>,
}

impl Interface {
    /// A boundary cut into equal faces.
    pub fn uniform(name: &'static str, faces: usize, face_area: Area) -> Interface {
        Interface {
            name,
            areas: vec![face_area.to_si().max(0.0); faces.max(1)],
        }
    }

    /// A boundary whose faces have their own areas.
    pub fn from_areas(name: &'static str, areas: Vec<Area>) -> Interface {
        let areas: Vec<f64> = areas.into_iter().map(|a| a.to_si().max(0.0)).collect();
        Interface {
            name,
            areas: if areas.is_empty() { vec![0.0] } else { areas },
        }
    }

    /// What this boundary is called. Both sides of a coupling must agree on it, and a
    /// mismatch is how they discover they meant different surfaces.
    pub fn name(&self) -> &'static str {
        self.name
    }

    /// How many faces it is cut into. This is the number both sides have to agree on.
    pub fn faces(&self) -> usize {
        self.areas.len()
    }

    /// Area of one face. Zero past the end, so a consumer walking the boundary need not
    /// bounds-check the kernel.
    pub fn area_of(&self, face: usize) -> Area {
        Area::from_si(self.areas.get(face).copied().unwrap_or(0.0))
    }

    /// Total area of the boundary.
    pub fn total_area(&self) -> Area {
        Area::from_si(self.areas.iter().sum())
    }

    /// Cumulative area up to and including each face, which is the coordinate a
    /// conservative resampling works in.
    fn cumulative(&self) -> Vec<f64> {
        let mut running = 0.0;
        let mut out = Vec::with_capacity(self.areas.len() + 1);
        out.push(0.0);
        for a in &self.areas {
            running += a;
            out.push(running);
        }
        out
    }
}

/// A quantity spread over an interface's faces, in SI units.
///
/// An amount per face rather than a density, so that summing is meaningful and the total
/// is a total. A density would need the areas to be carried alongside every arithmetic
/// operation, and the one thing this type must make easy is adding up.
#[derive(Clone, Debug, PartialEq)]
pub struct Flux {
    per_face: Vec<f64>,
}

impl Flux {
    /// Nothing, spread over the given number of faces. At least one face, so an empty flux
    /// is not a special case every consumer has to handle.
    pub fn zeros(faces: usize) -> Flux {
        Flux {
            per_face: vec![0.0; faces.max(1)],
        }
    }

    /// Amounts per face, in SI base units and in the interface's own order.
    pub fn from_faces(per_face: Vec<f64>) -> Flux {
        Flux {
            per_face: if per_face.is_empty() {
                vec![0.0]
            } else {
                per_face
            },
        }
    }

    /// One number spread over an interface in proportion to face area.
    ///
    /// The honest way to turn a lumped quantity into a distributed one: it says "evenly,
    /// because I do not know better" rather than silently putting everything on the first
    /// face. Which is what the lumped coupling was doing.
    pub fn spread_over(total: f64, interface: &Interface) -> Flux {
        let area = interface.total_area().to_si();
        if area <= 0.0 {
            let faces = interface.faces();
            return Flux::from_faces(vec![total / faces as f64; faces]);
        }
        Flux::from_faces(interface.areas.iter().map(|a| total * a / area).collect())
    }

    /// Distribute a total over an interface following a shape.
    ///
    /// `profile` is called with each face's centre as a fraction of the way along the
    /// boundary, from 0 to 1 in cumulative area, and returns an unnormalised weight. So a
    /// beam of waist `w` centred on a boundary of length `l` is
    /// `|u| (-2.0 * (((u - 0.5) * l / w).powi(2))).exp()`, written the way the physics is
    /// written rather than as a table of numbers.
    ///
    /// The weights are scaled so the faces sum to `total` exactly. Which splits the two
    /// claims deliberately: **the total is exact**, because a coupling that loses energy is
    /// a bug the audit must be able to trust, and **the shape is midpoint-accurate**,
    /// because a face gets its centre's weight rather than the profile's integral over it.
    /// The shape error falls as the boundary is refined; the total's does not exist.
    ///
    /// A profile summing to zero — all weights zero, or positives and negatives
    /// cancelling — has no scale to normalise against, so it falls back to spreading by
    /// area. That is a defined answer rather than infinities, and it is the same answer
    /// [`Flux::spread_over`] gives.
    pub fn profiled<F>(total: f64, interface: &Interface, mut profile: F) -> Flux
    where
        F: FnMut(f64) -> f64,
    {
        let span = interface.total_area().to_si();
        let mut weights = Vec::with_capacity(interface.faces());
        let mut running = 0.0;
        for area in &interface.areas {
            // The centre of this face in cumulative-area coordinates. Cumulative rather
            // than index, so an unevenly cut boundary is sampled where its faces actually
            // are — the point of letting faces carry their own areas.
            let centre = if span > 0.0 {
                (running + 0.5 * area) / span
            } else {
                0.5
            };
            running += area;
            weights.push(profile(centre));
        }
        let sum: f64 = weights.iter().sum();
        if !sum.is_finite() || sum == 0.0 {
            return Flux::spread_over(total, interface);
        }
        Flux::from_faces(weights.into_iter().map(|w| total * w / sum).collect())
    }

    /// How many faces this flux covers. Must match the interface it is published on.
    pub fn faces(&self) -> usize {
        self.per_face.len()
    }

    /// The amount on one face. Zero past the end.
    pub fn at(&self, face: usize) -> f64 {
        self.per_face.get(face).copied().unwrap_or(0.0)
    }

    /// Amounts per face, in the order the interface defines.
    pub fn per_face(&self) -> &[f64] {
        &self.per_face
    }

    /// Summed in index order, so the total is a function of the data and not of how it was
    /// visited.
    pub fn total(&self) -> f64 {
        self.per_face.iter().sum()
    }

    /// Largest single face's amount, which is the scale a rounding tolerance should be
    /// judged against for the same reason [`Ledger`](crate::Ledger) records one.
    pub fn largest(&self) -> f64 {
        self.per_face.iter().fold(0.0f64, |a, v| a.max(v.abs()))
    }

    /// Add another flux face by face. Refuses a mismatched face count rather than
    /// truncating or padding.
    pub fn add(&mut self, other: &Flux) -> Result<(), Violation> {
        if other.faces() != self.faces() {
            return Err(mismatch("flux addition", self.faces(), other.faces()));
        }
        for (mine, theirs) in self.per_face.iter_mut().zip(other.per_face.iter()) {
            *mine += theirs;
        }
        Ok(())
    }

    /// Scale every face.
    pub fn scaled(&self, by: f64) -> Flux {
        Flux::from_faces(self.per_face.iter().map(|v| v * by).collect())
    }

    /// Redistribute onto a different interface, conserving the total.
    ///
    /// Faces are treated as consecutive intervals in cumulative area, and each source
    /// face's amount is divided among the target faces it overlaps in proportion to how
    /// much of it each covers. The overlap fractions of any source face sum to one, so the
    /// total survives by construction rather than by a corrective scaling afterwards —
    /// which is the difference between a resampling that conserves and one that is
    /// checked and then adjusted.
    ///
    /// It conserves to summation rounding, not to the last bit: the pieces are added in a
    /// different order than they were split. That is a part in `10¹⁵`, and it is why the
    /// audit's tolerance is relative.
    ///
    /// Note what this cannot do. Redistributing by area assumes the two interfaces cover
    /// the same boundary in the same order, which is a statement about the geometry that
    /// this type has no way to check. It is a remap, not a projection between arbitrary
    /// meshes.
    pub fn resample(&self, from: &Interface, to: &Interface) -> Result<Flux, Violation> {
        if self.faces() != from.faces() {
            return Err(mismatch("resampling source", from.faces(), self.faces()));
        }
        let (source, target) = (from.cumulative(), to.cumulative());
        let source_span = source[source.len() - 1];
        let target_span = target[target.len() - 1];
        if source_span <= 0.0 || target_span <= 0.0 {
            // No area to redistribute over; spread evenly and say so by doing something
            // defined rather than dividing by zero.
            return Ok(Flux::spread_over(self.total(), to));
        }
        // Work in a normalised coordinate so the two boundaries need not have equal area —
        // a coating and the glass behind it are the same surface described twice.
        let mut out = vec![0.0; to.faces()];
        for i in 0..from.faces() {
            let (a0, a1) = (source[i] / source_span, source[i + 1] / source_span);
            let width = a1 - a0;
            if width <= 0.0 {
                continue;
            }
            for (j, slot) in out.iter_mut().enumerate() {
                let (b0, b1) = (target[j] / target_span, target[j + 1] / target_span);
                let overlap = a1.min(b1) - a0.max(b0);
                if overlap > 0.0 {
                    *slot += self.per_face[i] * (overlap / width);
                }
            }
        }
        Ok(Flux::from_faces(out))
    }
}

/// The violation a face-count disagreement raises.
///
/// Its own function because the message is the useful part: two numbers that should have
/// been the same, and where they were compared.
pub(crate) fn mismatch(site: &str, expected: usize, found: usize) -> Violation {
    Violation {
        quantity: format!("face count (expected {expected}, found {found})"),
        site: site.to_string(),
        before: expected as f64,
        after: found as f64,
        scale: expected.max(found) as f64,
        tolerance: 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cm2(v: f64) -> Area {
        Area::from_si(v * 1e-4)
    }

    #[test]
    fn an_interface_carries_its_faces_and_their_areas() {
        let uniform = Interface::uniform("plate", 8, cm2(0.5));
        assert_eq!(uniform.faces(), 8);
        assert!((uniform.total_area().to_si() - 8.0 * 0.5e-4).abs() < 1e-18);
        assert_eq!(uniform.area_of(0), uniform.area_of(7));
        assert_eq!(uniform.name(), "plate");

        // Unequal faces, as a cap cut into rings really is.
        let rings = Interface::from_areas("cap", vec![cm2(0.1), cm2(0.3), cm2(0.5)]);
        assert_eq!(rings.faces(), 3);
        assert!(rings.area_of(2) > rings.area_of(0));
        assert!((rings.total_area().to_si() - 0.9e-4).abs() < 1e-18);

        // Degenerate definitions give something usable rather than an empty vector that
        // every consumer would then have to guard against.
        assert_eq!(Interface::uniform("none", 0, cm2(1.0)).faces(), 1);
        assert_eq!(Interface::from_areas("none", vec![]).faces(), 1);
        assert!(Interface::from_areas("odd", vec![Area::from_si(-1.0)]).area_of(0) >= Area::ZERO);
    }

    /// Spreading a lumped total is area-weighted, which is what makes it the honest
    /// translation of "I do not know where it went".
    #[test]
    fn a_lumped_total_spreads_by_area() {
        let rings = Interface::from_areas("cap", vec![cm2(1.0), cm2(2.0), cm2(1.0)]);
        let flux = Flux::spread_over(8.0, &rings);
        assert_eq!(flux.faces(), 3);
        assert!((flux.total() - 8.0).abs() < 1e-12);
        // Twice the area, twice the share.
        assert!((flux.at(1) / flux.at(0) - 2.0).abs() < 1e-12);
        assert!((flux.at(0) - flux.at(2)).abs() < 1e-15);

        // An interface with no area at all divides evenly rather than dividing by zero.
        let empty = Interface::uniform("flat", 4, Area::ZERO);
        let flux = Flux::spread_over(8.0, &empty);
        assert!((flux.total() - 8.0).abs() < 1e-12);
        assert!((flux.at(0) - 2.0).abs() < 1e-12);
    }

    /// A face count that does not match is refused, with both numbers in the message. This
    /// is the refusal the whole design turns on: two domains either share a discretisation
    /// or say explicitly that they are crossing one.
    #[test]
    fn a_mismatched_face_count_is_refused_by_name() {
        let mut coarse = Flux::zeros(4);
        let fine = Flux::zeros(16);
        let err = coarse
            .add(&fine)
            .expect_err("4 and 16 must not silently combine");
        assert!(err.quantity.contains("expected 4"), "{err}");
        assert!(err.quantity.contains("found 16"), "{err}");

        // And resampling checks its source, since a flux carries no interface of its own.
        let from = Interface::uniform("a", 4, cm2(1.0));
        let err = Flux::zeros(7)
            .resample(&from, &Interface::uniform("b", 4, cm2(1.0)))
            .expect_err("a flux of 7 is not a flux over 4 faces");
        assert!(err.site.contains("resampling"), "{err}");
    }

    /// **The property resampling exists for.** Crossing discretisations conserves the
    /// total, because each source face's amount is divided among the targets it overlaps
    /// and those fractions sum to one.
    ///
    /// Checked both ways and at an awkward ratio, because a remap that happens to work for
    /// a factor of two can still be wrong for three into seven.
    #[test]
    fn resampling_conserves_the_total() {
        let cases = [
            (4usize, 16usize),
            (16, 4),
            (3, 7),
            (7, 3),
            (5, 5),
            (1, 9),
            (9, 1),
        ];
        for (n, m) in cases {
            let from = Interface::uniform("from", n, cm2(1.0));
            let to = Interface::uniform("to", m, cm2(1.0));
            // A ramp, so the answer is not uniform and a bug cannot hide in symmetry.
            let flux = Flux::from_faces((0..n).map(|i| 1.0 + i as f64).collect());
            let before = flux.total();

            let moved = flux.resample(&from, &to).unwrap();
            assert_eq!(moved.faces(), m, "{n} -> {m}: wrong face count");
            assert!(
                (moved.total() / before - 1.0).abs() < 1e-12,
                "{n} -> {m}: {} became {}",
                before,
                moved.total()
            );
        }
    }

    /// Resampling to the same interface changes nothing, which is the identity a remap has
    /// to satisfy before any of its other properties matter.
    #[test]
    fn resampling_onto_the_same_interface_is_the_identity() {
        let interface = Interface::from_areas("cap", vec![cm2(1.0), cm2(3.0), cm2(2.0)]);
        let flux = Flux::from_faces(vec![2.0, 7.0, -1.5]);
        let same = flux.resample(&interface, &interface).unwrap();
        for face in 0..interface.faces() {
            assert!(
                (same.at(face) - flux.at(face)).abs() < 1e-12,
                "face {face}: {} against {}",
                same.at(face),
                flux.at(face)
            );
        }
    }

    /// Refining and then coarsening back returns the original, when the coarse faces are
    /// unions of fine ones. A stronger statement than conservation: it says the
    /// distribution survived, not just its sum.
    #[test]
    fn a_round_trip_through_a_finer_grid_returns_the_distribution() {
        let coarse = Interface::uniform("coarse", 4, cm2(1.0));
        let fine = Interface::uniform("fine", 12, cm2(1.0) / 3.0);
        let flux = Flux::from_faces(vec![1.0, 5.0, 2.0, 9.0]);

        let refined = flux.resample(&coarse, &fine).unwrap();
        // Each coarse face became three fine ones holding a third each.
        assert!((refined.at(0) - 1.0 / 3.0).abs() < 1e-12);
        assert!(
            (refined.at(1) - refined.at(2)).abs() < 1e-15,
            "flat inside a coarse face"
        );
        // And the step is across the coarse boundary, between fine 2 and 3, not inside one.
        assert!(
            (refined.at(3) / refined.at(2) - 5.0).abs() < 1e-12,
            "the ramp survived"
        );

        let back = refined.resample(&fine, &coarse).unwrap();
        for face in 0..4 {
            assert!(
                (back.at(face) - flux.at(face)).abs() < 1e-12,
                "face {face}: {} against {}",
                back.at(face),
                flux.at(face)
            );
        }
    }

    /// Coarsening loses detail and cannot get it back, which is worth pinning so nobody
    /// treats a remap as lossless in both directions.
    #[test]
    fn coarsening_loses_what_it_averages() {
        let fine = Interface::uniform("fine", 8, cm2(1.0));
        let coarse = Interface::uniform("coarse", 2, cm2(4.0));
        // All the energy on one face.
        let spike = Flux::from_faces(vec![0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0]);

        let coarsened = spike.resample(&fine, &coarse).unwrap();
        assert!(
            (coarsened.total() - 10.0).abs() < 1e-12,
            "the total survives"
        );
        assert!(
            (coarsened.at(0) - 10.0).abs() < 1e-12,
            "and it stayed on its side"
        );

        // Going back spreads it over the whole coarse face: the spike is gone for good.
        let restored = coarsened.resample(&coarse, &fine).unwrap();
        assert!((restored.total() - 10.0).abs() < 1e-12);
        assert!(
            (restored.at(2) - 2.5).abs() < 1e-12,
            "a spike came back as an average, got {}",
            restored.at(2)
        );
        assert!(restored.at(0) > 0.0, "and it leaked onto its neighbours");
    }

    /// Adding fluxes face by face, and the scale a tolerance should use.
    #[test]
    fn fluxes_add_and_report_their_scale() {
        let mut a = Flux::from_faces(vec![1.0, -2.0, 0.5]);
        let b = Flux::from_faces(vec![0.5, 2.0, 0.5]);
        a.add(&b).unwrap();
        assert_eq!(a.per_face(), &[1.5, 0.0, 1.0]);
        assert!((a.total() - 2.5).abs() < 1e-15);
        // The largest single face, not the total, which is what a cancelling distribution
        // needs its rounding judged against.
        assert!((a.largest() - 1.5).abs() < 1e-15);

        let cancelling = Flux::from_faces(vec![1e6, -1e6]);
        assert!(cancelling.total().abs() < 1e-9);
        assert!((cancelling.largest() - 1e6).abs() < 1e-9);

        assert_eq!(a.scaled(2.0).per_face(), &[3.0, 0.0, 2.0]);
    }

    /// A profiled flux keeps the total exactly and puts the shape where it belongs, checked
    /// against the closed form the profile came from.
    #[test]
    fn a_profile_conserves_its_total_and_lands_where_the_beam_is() {
        // 20 mm of boundary in 41 faces, so face 20 is centred exactly at the middle.
        let faces = 41;
        let length = 20e-3;
        let plate = Interface::uniform("plate", faces, cm2(1.0));
        let waist = 3e-3;
        let profile = |u: f64| (-2.0 * (((u - 0.5) * length / waist).powi(2))).exp();

        let flux = Flux::profiled(0.096, &plate, profile);
        assert!((flux.total() - 0.096).abs() < 1e-15, "the total is exact");
        assert_eq!(flux.faces(), faces);

        // Peaked in the middle and symmetric about it, which is what a centred beam is.
        let peak = flux.largest();
        assert!(
            (flux.at(20) - peak).abs() < 1e-18,
            "the peak is at the centre face"
        );
        // Relatively, not to the last bit: the face centres come from a running sum of
        // areas, which rounds differently walking left and right from the middle. The
        // physics is symmetric; the arithmetic reaching it is not quite.
        for offset in 1..=20 {
            let (left, right) = (flux.at(20 - offset), flux.at(20 + offset));
            assert!(
                (left / right - 1.0).abs() < 1e-12,
                "asymmetric at offset {offset}: {left} against {right}"
            );
        }

        // And the ratio between two faces is the profile's own ratio, not something the
        // normalisation distorted: only the scale was changed, not the shape.
        let u = |i: usize| (i as f64 + 0.5) / faces as f64;
        assert!(
            (flux.at(14) / flux.at(20) - profile(u(14)) / profile(u(20))).abs() < 1e-12,
            "the normalisation must not bend the profile"
        );
        // A 3 mm waist on a 20 mm plate reaches the edges as e^(-2(10/3)^2) = 2e-10, so the
        // ends of the plate are dark. Which is the whole reason a lumped coupling is wrong
        // here: it would have warmed them equally.
        assert!(flux.at(0) / peak < 1e-9, "edge ratio {}", flux.at(0) / peak);
    }

    /// A profile with nothing to normalise against falls back to area rather than to
    /// infinity, and says so by giving the same answer as spreading.
    #[test]
    fn a_profile_with_no_scale_falls_back_to_area() {
        let rings = Interface::from_areas("cap", vec![cm2(1.0), cm2(3.0)]);
        for degenerate in [0.0, f64::NAN, f64::INFINITY] {
            let flux = Flux::profiled(4.0, &rings, |_| degenerate);
            let spread = Flux::spread_over(4.0, &rings);
            assert_eq!(
                flux.per_face(),
                spread.per_face(),
                "for weight {degenerate}"
            );
        }
        // Cancelling weights have a scale per face but none in sum, so they fall back too.
        let flux = Flux::profiled(4.0, &rings, |u| if u < 0.5 { 1.0 } else { -1.0 / 3.0 });
        assert!((flux.total() - 4.0).abs() < 1e-12);
    }

    /// An unevenly cut boundary is sampled where its faces actually are, in cumulative area
    /// rather than by index — which is the reason faces carry their own areas at all.
    #[test]
    fn an_uneven_boundary_is_sampled_at_its_faces() {
        // Three faces: a tenth, then eight tenths, then a tenth. Their centres in
        // cumulative area are 0.05, 0.5 and 0.95, nothing like 1/6, 1/2, 5/6.
        let uneven = Interface::from_areas("cap", vec![cm2(1.0), cm2(8.0), cm2(1.0)]);
        let mut seen = Vec::new();
        let _ = Flux::profiled(1.0, &uneven, |u| {
            seen.push(u);
            1.0
        });
        assert!((seen[0] - 0.05).abs() < 1e-12, "{:?}", seen);
        assert!((seen[1] - 0.50).abs() < 1e-12, "{:?}", seen);
        assert!((seen[2] - 0.95).abs() < 1e-12, "{:?}", seen);
    }

    /// Out-of-range faces read as nothing rather than panicking, since a consumer walking a
    /// boundary should not have to bounds-check the kernel's own data.
    #[test]
    fn reading_past_the_end_gives_nothing() {
        let flux = Flux::from_faces(vec![1.0, 2.0]);
        assert_eq!(flux.at(0), 1.0);
        assert_eq!(flux.at(5), 0.0);
        let interface = Interface::uniform("i", 2, cm2(1.0));
        assert_eq!(interface.area_of(9), Area::ZERO);
    }
}
