//! The report picks a view from the data's shape, for a researcher who should not have to.
//!
//! The claim is not that the pictures are pretty. It is that a run of *any* domain produces a
//! file that opens in a browser and shows something correct — without the caller knowing that a
//! bar wants a profile, a room wants a heatmap and a fluid wants a rotatable scene.

use dualis_world::{report, Scene, World};

fn frames_of(scene_file: &str) -> (String, Vec<dualis_world::Frame>) {
    let text = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("scenes")
            .join(scene_file),
    )
    .expect("the scene is there");
    let scene: Scene = serde_json::from_str(&text).expect("it parses");
    let title = scene.title.clone();
    let mut world = World::build(scene).expect("it builds");
    (title, world.run().expect("it conserves"))
}

/// **Each shape of data gets its own view, and the choice is made from the shape.**
///
/// Scene 14 has all four at once: a 1D bar, a 2D room, bodies in space, and scalars from every
/// domain including the two with no picture at all.
#[test]
fn every_shape_of_data_gets_the_view_that_suits_it() {
    let (title, frames) = frames_of("14-a-world.json");
    let html = report::html(&title, &frames);

    for kind in ["profile", "heatmap", "scene", "series"] {
        assert!(
            html.contains(&format!("data-kind=\"{kind}\"")),
            "no {kind} view in the report"
        );
    }
    // One card per drawable domain, plus one for the readings.
    assert_eq!(
        html.matches("class=\"card\"").count(),
        frames[0].panels.len() + 1
    );

    // Self-contained: no network, no library, nothing to install.
    for forbidden in ["http://", "https://", "src=", "@import"] {
        assert!(
            !html.contains(forbidden),
            "the report reaches outside itself with {forbidden:?}"
        );
    }
    assert!(html.starts_with("<!doctype html>"), "not a whole document");
    assert!(html.trim_end().ends_with("</html>"));
}

/// **A scene with nothing to draw still produces a report worth opening.**
///
/// Scene 13 is a winding and a thermal network: no field, no bodies, and the filmstrip is empty
/// for it. The report is not — every domain has readings, so it gets the series view, which for
/// this scene is the entire result.
#[test]
fn the_undrawable_scene_still_reports() {
    let (title, frames) = frames_of("13-winding-that-heats-itself.json");
    let html = report::html(&title, &frames);

    assert!(frames[0].panels.is_empty(), "scene 13 should draw nothing");
    assert!(html.contains("data-kind=\"series\""), "no series view");
    assert_eq!(html.matches("class=\"card\"").count(), 1);

    // The numbers themselves are in the file, so it is a record and not a viewer needing a
    // second file beside it.
    assert!(
        html.contains("\"winding\""),
        "the node names should travel with it"
    );
    assert!(html.contains("\"dissipating\""));
}
