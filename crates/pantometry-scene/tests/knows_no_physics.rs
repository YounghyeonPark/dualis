//! The scene layer captures a domain it has never heard of.
//!
//! That is the property the layer split exists for, and it is checkable in one way that matters:
//! define a physics *in this test file* — a crate `pantometry-scene` cannot possibly know — place it,
//! and capture it. If that works, adding a physics costs one crate and nothing here moves.

use pantometry_core::{Bodies, Domain, Exchange, Pose, Reading, ScalarField, Schedule, Simulation};
use pantometry_scene::{capture, settle_framing, Extent, Panel, PanelData, Placement};
use pantometry_units::{Length, LengthVec, Time};
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
    ) -> Result<(), pantometry_core::Violation> {
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
    assert_eq!(field.grid(), Some((9, 1, 1)));
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

/// A field that genuinely varies in all three directions, and knows what it should sample to.
struct Volume;

impl Domain for Volume {
    fn name(&self) -> &str {
        "volume"
    }
    fn step(
        &mut self,
        _t: Time,
        _dt: Time,
        _bus: &mut Exchange,
    ) -> Result<(), pantometry_core::Violation> {
        Ok(())
    }
    fn as_field(&self) -> Option<&dyn ScalarField> {
        Some(self)
    }
}

impl ScalarField for Volume {
    /// `100x + 10y + z`, in metres. Every sample is a different number, and the number *says*
    /// where it came from — so a capture that muddled two axes or collapsed one is legible in the
    /// values rather than only in a count.
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let v = p.to_si();
        100.0 * v.x + 10.0 * v.y + v.z
    }
    fn unit(&self) -> &'static str {
        "u"
    }
}

/// **A three-dimensional field is captured whole, not as a slice through it.**
///
/// This is the test that would have failed before `Extent` had a third axis, and it would have
/// failed *silently*: `samples` was a pair, the sampler built its position as `(u, v, 0)`, and a
/// solid came back as its `z = 0` face. A perfectly plausible picture of a block, one third of a
/// dimension short, with nothing anywhere to say so.
///
/// It took a domain with a real 3D field — `pantometry-thermal`'s `Solid3D` — to make the gap
/// visible. Nothing in this crate could have noticed on its own, because every field it had ever
/// been handed was flat.
#[test]
fn a_volume_is_captured_in_three_dimensions() {
    let sim = Simulation::new(Schedule::OneWay).with(Volume);
    let placed = BTreeMap::from([(
        "volume".to_string(),
        Placement::field(Extent::volume(LengthVec::m(1.0, 1.0, 1.0), 2, 3, 4)),
    )]);

    let frame = capture(&sim, &placed);
    let panel = &frame.panels[0];
    assert_eq!(panel.grid(), Some((2, 3, 4)));
    assert_eq!(panel.values().len(), 24, "every sample, not one face");

    // Each value names its own position, so the ordering is checkable rather than assumed:
    // x fastest, then y, then z.
    for k in 0..4 {
        for j in 0..3 {
            for i in 0..2 {
                let want = 100.0 * (i as f64 / 1.0) + 10.0 * (j as f64 / 2.0) + k as f64 / 3.0;
                let got = panel.values()[i + 2 * (j + 3 * k)];
                assert!(
                    (got - want).abs() < 1e-12,
                    "({i},{j},{k}): got {got}, expected {want} — the axes are muddled"
                );
            }
        }
    }

    // The z axis really varied, which is the whole point: before the fix every slice was slice 0.
    let first = panel.slice(0).expect("slice 0");
    let last = panel.slice(3).expect("slice 3");
    assert!(
        first != last,
        "every slice is identical, so z was never sampled"
    );
    assert!(panel.slice(4).is_none(), "there is no fifth slice");

    // And the shape is reported as three-dimensional rather than as a tall plane.
    let extent = placed["volume"].extent.unwrap();
    assert_eq!(extent.dimensions(), 3);
    assert_eq!(extent.count(), 24);
}

/// **An axis asked for at one sample is read in the middle of it, not at its corner.**
///
/// A flat extent samples the same point either way, so this only bites for an extent with real
/// thickness that a caller chose to summarise with one plane — asking for the low face there
/// would report a boundary as if it were the body.
#[test]
fn a_single_sample_across_a_thick_axis_is_taken_in_the_middle() {
    let sim = Simulation::new(Schedule::OneWay).with(Volume);
    let placed = BTreeMap::from([(
        "volume".to_string(),
        Placement::field(Extent::volume(LengthVec::m(1.0, 1.0, 1.0), 2, 1, 1)),
    )]);
    let frame = capture(&sim, &placed);
    let v = frame.panels[0].values();

    // y and z are one sample across a 1 m span, so both are read at 0.5 m: 10*0.5 + 1*0.5 = 5.5.
    assert!((v[0] - 5.5).abs() < 1e-12, "x=0 gave {}, wanted 5.5", v[0]);
    assert!(
        (v[1] - 105.5).abs() < 1e-12,
        "x=1 gave {}, wanted 105.5",
        v[1]
    );
}
