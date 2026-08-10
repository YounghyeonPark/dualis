//! The scene layer captures a domain it has never heard of.
//!
//! That is the property the layer split exists for, and it is checkable in one way that matters:
//! define a physics *in this test file* — a crate `dualis-scene` cannot possibly know — place it,
//! and capture it. If that works, adding a physics costs one crate and nothing here moves.

use dualis_core::{Bodies, Domain, Exchange, Pose, Reading, ScalarField, Schedule, Simulation};
use dualis_scene::{capture, settle_framing, Extent, Panel, PanelData, Placement};
use dualis_units::{Length, LengthVec, Time};
use std::collections::BTreeMap;

/// A physics nobody has written before: a scalar that decays, spread over a line, with two
/// motes drifting through it.
struct Invented {
    t: f64,
}

impl Domain for Invented {
    fn name(&self) -> &str {
        "invented"
    }
    fn step(
        &mut self,
        _t: Time,
        dt: Time,
        _bus: &mut Exchange,
    ) -> Result<(), dualis_core::Violation> {
        self.t += dt.to_si();
        Ok(())
    }
    fn as_field(&self) -> Option<&dyn ScalarField> {
        Some(self)
    }
    fn as_bodies(&self) -> Option<&dyn Bodies> {
        Some(self)
    }
    fn readings(&self) -> Vec<Reading> {
        vec![Reading::new("invented", "elapsed", self.t, "s")]
    }
}

impl ScalarField for Invented {
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        (-p.to_si().x).exp() * (1.0 + self.t)
    }
    fn unit(&self) -> &'static str {
        "widgets"
    }
}

impl Bodies for Invented {
    fn count(&self) -> usize {
        2
    }
    fn position(&self, i: usize) -> LengthVec {
        LengthVec::m(self.t * (i as f64 + 1.0), 0.0, 0.0)
    }
    fn value(&self, i: usize) -> f64 {
        i as f64
    }
    fn value_unit(&self) -> &'static str {
        "index"
    }
}

/// **A domain this crate cannot know is captured in all three shapes.**
#[test]
fn a_physics_invented_in_a_test_file_is_captured_whole() {
    let mut sim = Simulation::new(Schedule::Staggered).with(Invented { t: 0.0 });
    let mut placed = BTreeMap::new();
    placed.insert(
        "invented".to_string(),
        Placement::field(Extent::line(Length::m(2.0), 9)),
    );

    let mut frames = Vec::new();
    for _ in 0..5 {
        sim.advance(Time::s(0.5))
            .expect("it conserves nothing and claims nothing");
        frames.push(capture(&sim, &placed));
    }
    settle_framing(&mut frames);

    let last = frames.last().unwrap();

    // The field, sampled over the extent it was *placed* with — the crate was told how big it
    // is, because a `ScalarField` does not know where it stops.
    let field = last
        .panels
        .iter()
        .find(|p| matches!(p.data, PanelData::Field { .. }))
        .expect("a field panel");
    assert_eq!(field.grid(), Some((9, 1)));
    assert_eq!(field.unit, "widgets", "the field named its own unit");
    // Sampled at x = 0 and x = 2, so the first value is e^0 and the last e^-2, times (1 + t).
    let v = field.values();
    assert!((v[0] / (1.0 + last.time_s) - 1.0).abs() < 1e-12);
    assert!((v[8] / ((1.0 + last.time_s) * (-2.0f64).exp()) - 1.0).abs() < 1e-12);

    // The scalars, which every domain has whether or not it draws.
    assert_eq!(last.readings.len(), 1);
    assert_eq!(last.readings[0].label, "elapsed");
    assert!((last.readings[0].value - 2.5).abs() < 1e-12);
}

/// **A placed domain is captured where it was placed.**
///
/// `Pose` is physical, so bodies come back in world coordinates. Checked against a quarter turn
/// worked by hand: a mote on the local +x axis lands on world +y.
#[test]
fn bodies_come_back_in_world_coordinates() {
    let mut sim = Simulation::new(Schedule::Staggered).with(Invented { t: 0.0 });
    sim.advance(Time::s(1.0)).unwrap();

    let turned = Pose::new(
        LengthVec::m(0.0, 0.0, 5.0),
        glam::DQuat::from_rotation_z(std::f64::consts::FRAC_PI_2),
    );
    let mut placed = BTreeMap::new();
    placed.insert("invented".to_string(), Placement::at(turned));

    let frame = capture(&sim, &placed);
    let Some(Panel {
        data: PanelData::Points { positions, .. },
        ..
    }) = frame.panels.iter().find(|p| p.grid().is_none())
    else {
        panic!("expected bodies");
    };

    // Mote 0 sits at local (1, 0, 0) after one second; turned a quarter turn and lifted 5 m.
    let p = positions[0];
    assert!(p[0].abs() < 1e-12, "x should be 0, got {}", p[0]);
    assert!((p[1] - 1.0).abs() < 1e-12, "y should be 1, got {}", p[1]);
    assert!((p[2] - 5.0).abs() < 1e-12, "z should be 5, got {}", p[2]);
}

/// **A field with no extent is not drawn**, rather than drawn over a region nobody chose.
#[test]
fn a_field_nobody_sized_is_left_out() {
    let mut sim = Simulation::new(Schedule::Staggered).with(Invented { t: 0.0 });
    sim.advance(Time::s(1.0)).unwrap();

    // Placed, but without an extent.
    let mut placed = BTreeMap::new();
    placed.insert("invented".to_string(), Placement::default());
    let frame = capture(&sim, &placed);

    // It falls through to bodies, which need no extent — and no field panel appears.
    assert!(
        !frame.panels.iter().any(|p| p.grid().is_some()),
        "a field was drawn over a region nobody specified"
    );
    // The readings are still there: not drawing is not the same as not reporting.
    assert_eq!(frame.readings.len(), 1);
}
