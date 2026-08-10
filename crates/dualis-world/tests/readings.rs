//! Every domain reports scalars, including the ones that have no picture.
//!
//! Eight of the fourteen shipped scenes contain a domain the filmstrip cannot draw — a heater, a
//! lamp, a winding, a thermal network — and for several of them the scalar *is* the result. The
//! winding temperature is what decides whether a motor survives, and until `Frame::readings` it
//! was reachable only by reading the terminal.

use dualis_world::{Scene, World};

fn run(path: &str) -> Vec<dualis_world::Frame> {
    let text = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("scenes")
            .join(path),
    )
    .expect("the scene is there");
    let scene: Scene = serde_json::from_str(&text).expect("it parses");
    let mut world = World::build(scene).expect("it builds");
    world.run().expect("it conserves")
}

/// **A scene with nothing to draw still reports everything it does.**
///
/// Scene 13 has two domains and neither has a field: a winding and a thermal network. The
/// filmstrip is empty for it — correctly, since a conductance is not a distance — and the run
/// is still entirely readable.
#[test]
fn the_undrawable_scene_is_fully_reported() {
    let frames = run("13-winding-that-heats-itself.json");
    let last = frames.last().expect("frames");
    let named: Vec<String> = last
        .readings
        .iter()
        .map(|r| format!("{}.{}", r.domain, r.label))
        .collect();

    for want in [
        "coil.dissipating",
        "coil.resistance",
        "coil.spent",
        "motor.winding",
        "motor.stator",
        "motor.housing",
    ] {
        assert!(named.iter().any(|n| n == want), "missing {want}: {named:?}");
    }

    let get = |q: &str| {
        last.readings
            .iter()
            .find(|r| format!("{}.{}", r.domain, r.label) == q)
            .map(|r| r.value)
            .unwrap_or_else(|| panic!("no {q}"))
    };

    // The feedback, visible in the table: a winding that heats up dissipates more. Both ends
    // computed here — 62 m of 0.35 mm² copper at 1.75 A is 12.46 W cold, and the ratio to the
    // settled value must be the resistance ratio, because `I` is fixed.
    let first = frames.first().expect("frames");
    let cold = first
        .readings
        .iter()
        .find(|r| r.label == "dissipating")
        .expect("cold reading")
        .value;
    assert!(
        (cold - 12.4558).abs() < 1e-3,
        "cold dissipation {cold:.4} W against I²R₂₀ at 25 °C"
    );
    let hot = get("coil.dissipating");
    let r_ratio = get("coil.resistance")
        / first
            .readings
            .iter()
            .find(|r| r.label == "resistance")
            .expect("cold resistance")
            .value;
    assert!(
        (hot / cold - r_ratio).abs() < 1e-9,
        "at fixed current the power ratio is the resistance ratio: {:.6} against {r_ratio:.6}",
        hot / cold
    );

    // And the network's own ordering, which is why it is three columns and not one.
    assert!(get("motor.winding") > get("motor.stator"));
    assert!(get("motor.stator") > get("motor.housing"));

    // Units travel with the values, so a header row needs no separate legend.
    let watts = last
        .readings
        .iter()
        .find(|r| r.label == "dissipating")
        .unwrap();
    assert_eq!(watts.unit, "W");
    assert_eq!(
        last.readings
            .iter()
            .find(|r| r.label == "winding")
            .unwrap()
            .unit,
        "C"
    );
}

/// **A drawable scene reports scalars too, and they agree with the picture.**
///
/// The two assets are of the same run, so the bar's peak in the readings and the bar's peak in
/// the panel must be the same number — otherwise one of them is describing a different frame.
#[test]
fn the_readings_and_the_panel_describe_the_same_frame() {
    let frames = run("14-a-world.json");
    for frame in &frames {
        let panel = frame
            .panels
            .iter()
            .find(|p| p.name == "bar")
            .expect("the bar is drawable");
        // The panel is **kelvin** and the reading is celsius, and neither is a mistake: the field
        // returns what the cells hold, `readings` returns the unit a column of temperatures is
        // read in, and each says which it is. So the conversion is done here, from the units the
        // two sides declare, rather than assumed — which is the only version of this test that
        // would notice if one of them changed.
        //
        // It changed once already. While the application sampled fields itself it applied a
        // -273.15 offset and labelled the result "C"; moving that to `dualis-scene` moved the
        // conversion to whoever wants celsius, and this assertion is what caught it.
        assert_eq!(panel.unit, "K", "the bar's field is what the cells hold");
        let from_panel = panel.values().iter().cloned().fold(f64::MIN, f64::max) - 273.15;
        let from_readings = frame
            .readings
            .iter()
            .find(|r| r.domain == "bar" && r.label == "peak")
            .expect("the bar reports a peak");
        assert_eq!(
            from_readings.unit, "C",
            "a column of temperatures is read in celsius"
        );
        let from_readings = from_readings.value;
        assert!(
            (from_panel - from_readings).abs() < 1e-9,
            "at t = {:.4}: panel says {from_panel:.6}, readings say {from_readings:.6}",
            frame.time_s
        );
    }

    // The sources report what is left, which the picture never showed at all.
    let last = frames.last().unwrap();
    for source in ["beam", "lamp"] {
        let left = last
            .readings
            .iter()
            .find(|r| r.domain == source && r.label == "reserve")
            .unwrap_or_else(|| panic!("{source} reports no reserve"));
        assert!(left.value >= 0.0, "{source} owes energy: {}", left.value);
        assert_eq!(left.unit, "J");
    }
}

/// **A body panel's framing does not move between frames.**
///
/// `PanelData::Points` promises it — "fixed for the whole run so a body moving is a body moving
/// and not the frame rescaling underneath it" — and the promise became breakable when framing
/// moved out of the scene format and into a measurement. `capture` sees one frame and cannot see
/// the future, so without a pass over the finished run each frame would be framed to itself.
///
/// Caught by writing this rather than by anything failing, which is the wrong order.
#[test]
fn a_moving_body_moves_and_the_picture_does_not() {
    let frames = run("06-orbits.json");
    let first = frames
        .first()
        .and_then(|f| f.panels.iter().find(|p| p.name == "orbits"))
        .or_else(|| frames.first().and_then(|f| f.panels.first()))
        .expect("the orbit draws bodies");
    let dualis_world::PanelData::Points { bounds: want, .. } = first.data else {
        panic!("expected bodies");
    };

    for (i, frame) in frames.iter().enumerate() {
        for panel in &frame.panels {
            if let dualis_world::PanelData::Points { bounds, boxed, .. } = panel.data {
                assert!(!boxed, "an orbit has no wall to report");
                assert_eq!(bounds, want, "frame {i} reframed {}", panel.name);
            }
        }
    }

    // And the bodies really did move inside that fixed frame, or the check is vacuous.
    let moved = frames
        .iter()
        .filter_map(|f| f.panels.first())
        .filter_map(|p| match &p.data {
            dualis_world::PanelData::Points { positions, .. } => positions.first().copied(),
            _ => None,
        })
        .collect::<Vec<_>>();
    let spread = moved
        .iter()
        .flat_map(|p| p.iter())
        .fold(f64::MIN, |m, v| m.max(*v))
        - moved
            .iter()
            .flat_map(|p| p.iter())
            .fold(f64::MAX, |m, v| m.min(*v));
    assert!(spread > 0.0, "nothing moved, so nothing was being framed");
}

/// **A field domain with no extent is not drawn, and nothing says so.**
///
/// Moving the capture into `dualis-scene` moved the extent out of the sampler and into
/// `DomainSpec::placement`, where it is now a `match` a future variant could join without being
/// noticed. The failure mode is not a panic: `capture` skips a field that nobody sized, and the
/// result is a report that is one panel short. Nobody counts panels in a report.
///
/// So it is counted here, across every scene that ships, from the two sides that must agree —
/// the simulation's own answer to "do you have a field", and the scene's answer to "how big is
/// it". A new field domain that forgets its extent fails this and nothing else.
#[test]
fn every_domain_with_a_field_was_given_an_extent() {
    let (mut checked, mut fields) = (0, 0);
    for entry in std::fs::read_dir(std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes"))
        .expect("the scenes are there")
    {
        let path = entry.expect("readable").path();
        if path.extension().and_then(|e| e.to_str()) != Some("json") {
            continue;
        }
        let scene: Scene = serde_json::from_str(&std::fs::read_to_string(&path).expect("readable"))
            .expect("it parses");
        let world = World::build(scene.clone()).expect("it builds");
        let placed = world.placements();

        for spec in &scene.domains {
            let has_field = world
                .simulation()
                .domain(spec.name())
                .and_then(|d| d.as_field())
                .is_some();
            let sized = placed[spec.name()].extent.is_some();
            assert_eq!(
                has_field,
                sized,
                "{}: {} has a field: {has_field}, was given an extent: {sized}",
                path.display(),
                spec.name()
            );
            checked += 1;
            fields += usize::from(has_field);
        }
    }
    // Both counters, because the assertion above is an equality and passes perfectly on a run
    // where nothing has a field at all. That is the vacuous version of this test, and it is the
    // one a rewire that broke `as_field` everywhere would leave behind.
    assert!(checked >= 15, "only {checked} domains were checked");
    assert!(fields >= 9, "only {fields} of them had a field to size");
}
