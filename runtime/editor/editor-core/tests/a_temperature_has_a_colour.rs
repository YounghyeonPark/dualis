//! The colours a run's temperatures are drawn in, checked without a window.
//!
//! The viewer learned this once already: a window nobody can photograph proves the program did
//! not panic, which is a much weaker claim than it looks. The editor's viewport paints a
//! temperature field in the colour Planck says that temperature is, and the claim worth pinning
//! is not that it painted *something* — it is that the colours are the ones physics computed and
//! that the fallback fires exactly when physics has no colour to give.
//!
//! These live here rather than in the shell because the shell is a window and this is
//! arithmetic. What the shell adds — a depth sort and a splat radius — is geometry that
//! `viewer-core`'s own tests already cover in the form that matters.

/// **A room-temperature block gets no computed colour, and that is a physical fact rather than
/// a threshold.** Nothing at 300 K emits visible light; a renderer that gave it one would be
/// inventing an appearance the physics does not have, so the viewport falls back to false
/// colour and says so on the canvas.
#[test]
fn a_cool_field_has_no_colour_of_its_own() {
    // The fraction of radiated power that lands in the visible band, which is what the shell
    // tests to choose its colouring.
    assert!(dualis::view::glow_fraction(293.15) < 1e-20);
    assert!(dualis::view::glow_fraction(500.0) < 1e-6);
    // And the shell's threshold sits above those: a warm block is drawn in false colour.
    assert!(dualis::view::glow_fraction(700.0) < 1e-6);
}

/// **A melting or glowing field gets the colour it actually has.** Above about 900 K the visible
/// share is past the shell's threshold and the colours are Planck's: dull red first, then orange,
/// then white — the sequence a blacksmith names, arrived at by integrating a spectrum rather
/// than by anybody choosing three colours.
#[test]
fn a_glowing_field_is_drawn_the_colour_it_is() {
    assert!(dualis::view::glow_fraction(1000.0) > 1e-6);

    let dull_red = dualis::view::blackbody_srgb(1000.0);
    let orange = dualis::view::blackbody_srgb(1800.0);
    let white = dualis::view::blackbody_srgb(6000.0);

    // Red at the bottom: the red channel saturated, no blue at all — the colour is outside
    // sRGB's gamut on that side, which is why a photograph of hot iron looks like this.
    assert_eq!(dull_red[0], 255);
    assert_eq!(dull_red[2], 0);
    // Orange in the middle: green has climbed. **Not** "blue has opened" — this test asserted
    // that and failed, because at 1800 K the chromaticity is still outside sRGB's gamut on the
    // blue side and the channel is clamped to zero exactly as at 1000 K. The gamut is a
    // property of the monitor, not of the fire; the physics is in the chromaticity, and that
    // is where the claim belongs.
    assert!(orange[1] > dull_red[1]);
    let (x_dull, _) = dualis::view::planckian_chromaticity(1000.0);
    let (x_orange, _) = dualis::view::planckian_chromaticity(1800.0);
    assert!(
        x_orange < x_dull,
        "1800 K should sit bluer along the locus than 1000 K"
    );
    // White at the top: all three channels within a quarter of each other.
    let lo = white.iter().copied().min().unwrap() as f64;
    let hi = white.iter().copied().max().unwrap() as f64;
    assert!(hi / lo < 1.25, "6000 K should be near white, got {white:?}");
}

/// **Steel's melting point comes out the colour steel looks at its melting point.**
///
/// 1811 K for iron, and the visual is not a matter of opinion: mill scale at that temperature is
/// the bright yellow-orange every foundry photograph shows. Checked as a hue relation — red the
/// largest channel, green well up, blue present but least — rather than as three exact bytes,
/// because the byte values carry the colour-matching fit's percent-level error and the hue
/// ordering does not.
#[test]
fn iron_at_its_melting_point_is_the_colour_iron_is() {
    let c = dualis::view::blackbody_srgb(1811.0);
    assert!(
        c[0] > c[1] && c[1] > c[2],
        "molten iron should run red > green > blue, got {c:?}"
    );
    // Bright rather than dim: the red channel is saturated and green is past half.
    assert_eq!(c[0], 255);
    assert!(
        (128..=220).contains(&c[1]),
        "the yellow-orange of a melt has green well up but not equal, got {c:?}"
    );
}

/// **The colour a scene's own run produces is the colour of its own temperatures**, end to end:
/// build a scene with a block hot enough to glow, run it, read the field's values back out of
/// the wire format the viewport draws from, and check the colours those values map to.
///
/// This is the check that would fail if the panel's unit ever stopped being kelvin, if the
/// values arrived scaled, or if the field stopped reaching the frame at all — none of which the
/// unit tests above can see, because they never touch a simulation.
#[test]
fn a_run_carries_temperatures_the_viewport_can_colour() {
    let scene = r#"{
  "title": "a block hot enough to glow",
  "duration_s": 0.5,
  "frames": 2,
  "domains": [
    { "kind": "block", "name": "block", "cells": [4, 4, 4], "cell_mm": 5.0,
      "initial_c": 1200.0 }
  ]
}"#;
    let json = editor_core::run(scene).expect("the scene runs");
    let run = viewer_core::Run::from_json(&json).expect("the viewer's reader accepts it");
    let panel = run.frames[0]
        .panels
        .iter()
        .find(|p| p.name() == "block")
        .expect("the block has a panel");

    // The unit is what the viewport dispatches on — data, not a domain name.
    assert_eq!(panel.unit(), "K", "the viewport colours on this");
    let hottest = panel.values().iter().copied().fold(f64::MIN, f64::max);
    assert!(
        (1470.0..1480.0).contains(&hottest),
        "1200 C should reach the frame as about 1473 K, got {hottest}"
    );

    // And at that temperature the field glows, so the viewport draws Planck's colour.
    assert!(dualis::view::glow_fraction(hottest) > 1e-6);
    let c = dualis::view::blackbody_srgb(hottest);
    assert!(
        c[0] > c[1] && c[1] > c[2],
        "a 1473 K block is orange: {c:?}"
    );
}
