//! Run a scene and draw it.
//!
//! The asset is chosen by the output file's extension, because a run has several shapes and
//! only one of them is a picture.
//!
//! ```sh
//! dualis-world                        # the built-in scene, checked, nothing written
//! dualis-world scene.json out.html    # a report: a view per domain, chosen from its shape
//! dualis-world scene.json out.svg     # a filmstrip: every frame, one page, still
//! dualis-world scene.json out.csv     # every domain's scalars over time, one row per frame
//! dualis-world scene.json out.json    # the frames themselves — fields, bodies and readings
//! dualis-world --emit-default s.json  # write the built-in scene out to start from
//! ```
//!
//! `.html` is the one to reach for if you do not already know how you want this drawn. It picks
//! a view from each domain's *shape* — a profile for a 1D field, a heatmap for a 2D one, a
//! rotatable scene for bodies, a line chart for scalars — and opens in a browser with nothing
//! installed. `.csv` is the one that changed what this crate can say. Eight of the fourteen shipped scenes
//! have a domain the filmstrip cannot draw, and for several of them the scalar *is* the result —
//! scene 13's whole subject is a winding whose resistance follows its own temperature, and it
//! drew nothing at all.

use dualis::prelude::ThermalNetwork;
// The views are a library now — `dualis-view`, above `dualis-scene`, above the kernel. They
// were modules in this binary, which is `publish = false`, so a consumer who could state a
// simulation and run it could not draw one.
use dualis::view::{html, readings_csv, svg as filmstrip, to_json};
use dualis_world::{Scene, World};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();

    if args.first().map(String::as_str) == Some("--emit-default") {
        let path = args.get(1).ok_or("--emit-default needs a path")?;
        std::fs::write(path, serde_json::to_string_pretty(&default_scene())? + "\n")?;
        println!("wrote {path}");
        return Ok(());
    }

    let scene = match args.first() {
        Some(path) => serde_json::from_str::<Scene>(&std::fs::read_to_string(path)?)
            .map_err(|e| format!("{path}: {e}"))?,
        None => default_scene(),
    };
    let out = args.get(1);

    println!("{}", scene.title);
    println!(
        "  {} domain(s), {:.3} s in {} frames, drift budget {:.0e}",
        scene.domains.len(),
        scene.duration_s,
        scene.frames,
        scene.conservation_tolerance
    );

    let mut world = World::build(scene)?;
    let frames = match world.run() {
        Ok(frames) => frames,
        // The whole reason to use this library: a run that stopped conserving says which
        // quantity and by how much, instead of quietly drawing something plausible.
        Err(v) => {
            eprintln!(
                "\nthe audit stopped the run at t = {:.4} s",
                world.time().to_si()
            );
            eprintln!("  {v}");
            std::process::exit(1);
        }
    };

    // A row per *domain*, not per panel. A domain with nothing to draw used to print no line
    // at all, so a scene of two coupled domains reported one row and looked complete — and a
    // scene where every domain was undrawable printed a header, a zero-byte SVG and exit 0.
    // A gap a reader can see beats an absence they cannot.
    let last = &frames[frames.len() - 1];
    for spec in &world.scene().domains {
        let Some(panel) = last.panels.iter().find(|p| p.name == spec.name()) else {
            // A network has no field on purpose — a conductance is not a distance — but the
            // number it exists to produce is the *drop across a joint*, and "not drawn" reports
            // none of it. The picture is not the only output.
            if let Some(net) = world.simulation().domain_as::<ThermalNetwork>(spec.name()) {
                let nodes: Vec<_> = net.handles().collect();
                println!("  {:<14} {:<12} node temperatures", spec.name(), "network");
                for (i, (node, label)) in nodes.iter().enumerate() {
                    // The drop against the previous node, which is the number a network exists
                    // to give and the one a single lumped mass cannot: it reports the housing
                    // and the winding as one temperature.
                    let drop = match i.checked_sub(1).map(|j| nodes[j]) {
                        Some((up, up_label)) => format!(
                            ",  {:.2} K below {up_label}",
                            net.temperature(up).to_si() - net.temperature(*node).to_si()
                        ),
                        None => String::new(),
                    };
                    println!(
                        "    {:<12} {:>8.2} C{}",
                        label,
                        net.temperature(*node).to_si() - 273.15,
                        drop
                    );
                }
                continue;
            }
            println!(
                "  {:<14} {:<12} no field and no bodies — not drawn",
                spec.name(),
                "—"
            );
            continue;
        };
        let shape = match panel.grid() {
            Some((nx, ny, 1)) => format!("{nx} x {ny}"),
            Some((nx, ny, nz)) => format!("{nx} x {ny} x {nz}"),
            None => format!("{} bodies", panel.values().len()),
        };
        // The run-wide extremum beside the final value. The final value alone cannot tell a
        // ball that bounced half a metre from one that never moved: both end at zero.
        let now = panel.values().iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let over_run = frames
            .iter()
            .flat_map(|f| f.panels.iter().filter(|p| p.name == spec.name()))
            .flat_map(|p| p.values().iter())
            .fold(0.0f64, |m, v| m.max(v.abs()));
        println!(
            "  {:<14} {:<12} |{}| {:.4} now, {:.4} peak over the run",
            panel.name, shape, panel.unit, now, over_run
        );
    }

    // The asset to write is chosen by the extension, not by a flag.
    //
    // A run has several shapes and only one of them is a picture. A field is a raster, a body is
    // a point in space, and a source is a number over time — and for eight of the fourteen
    // shipped scenes at least one domain is that last thing, which the SVG could not draw at
    // all. `.csv` is the asset those domains always had and nothing collected.
    //
    // Extension rather than `--format` because the caller is already naming the file, and two
    // ways to say the same thing is one way to disagree with yourself.
    match out {
        Some(path) => {
            let ext = std::path::Path::new(path)
                .extension()
                .and_then(|e| e.to_str())
                .unwrap_or("")
                .to_ascii_lowercase();
            match ext.as_str() {
                "svg" => {
                    let svg = filmstrip(&world.scene().title, &frames, 6);
                    // A zero-byte file used to be written and reported as "0 KiB" — which a
                    // legitimate 937-byte strip also reports, so the one number on the line
                    // could not tell them apart. Bytes now, and an empty picture is refused
                    // rather than saved.
                    if svg.is_empty() {
                        eprintln!(
                            "\nnothing to draw: none of the {} domain(s) has a field or bodies. \
                             Try a .csv, which every domain can fill.",
                            world.scene().domains.len()
                        );
                        std::process::exit(1);
                    }
                    std::fs::write(path, &svg)?;
                    println!("  wrote {path} ({} bytes, filmstrip)", svg.len());
                }
                "csv" => {
                    let csv = readings_csv(&frames);
                    std::fs::write(path, &csv)?;
                    println!(
                        "  wrote {path} ({} bytes, {} columns over {} rows)",
                        csv.len(),
                        frames.first().map_or(0, |f| f.readings.len()),
                        frames.len()
                    );
                }
                "html" => {
                    let page = html(world.scene().title.as_str(), &frames);
                    std::fs::write(path, &page)?;
                    println!(
                        "  wrote {path} ({} bytes, a report that opens in a browser)",
                        page.len()
                    );
                }
                "gltf" => {
                    // The last frame, because glTF is geometry and not an animation: it moves
                    // node transforms and morph targets, and a field whose values change is
                    // neither. The frames themselves are what `.json` is for.
                    let Some(last) = frames.last() else {
                        eprintln!(
                            "
nothing to export: the run produced no frames"
                        );
                        std::process::exit(1);
                    };
                    let out = dualis::view::gltf(world.scene().title.as_str(), last);
                    for note in &out.skipped {
                        println!("  not exported: {note}");
                    }
                    if out.document.contains("\"nodes\":[]") {
                        eprintln!(
                            "\nnothing in this scene is geometry, so the file would open \
                             empty. Try .html or .csv."
                        );
                        std::process::exit(1);
                    }
                    std::fs::write(path, &out.document)?;
                    println!(
                        "  wrote {path} ({} bytes, geometry for Blender, three.js or a USD tool)",
                        out.document.len()
                    );
                }
                "json" => {
                    let json = to_json(world.scene().title.as_str(), &frames);
                    std::fs::write(path, &json)?;
                    println!(
                        "  wrote {path} ({} bytes, {} frames)",
                        json.len(),
                        frames.len()
                    );
                }
                other => {
                    eprintln!(
                        "\ndon't know how to write {other:?}. Known: .svg a filmstrip, \
                         .csv every domain's scalars over time, .json the frames themselves,                          .gltf the geometry."
                    );
                    std::process::exit(1);
                }
            }
        }
        None => println!("  name an output file: .html, .svg, .csv, .json or .gltf"),
    }
    Ok(())
}

/// A room ringing in its (1,1) mode, which is the cheapest scene that is worth looking at.
fn default_scene() -> Scene {
    serde_json::from_str(
        r#"{
  "title": "a small room ringing in its (1,1) mode",
  "schedule": "multirate",
  "duration_s": 0.02,
  "frames": 11,
  "conservation_tolerance": 1e-6,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.4, "height_m": 3.1,
      "cells_across": 61,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } }
  ]
}"#,
    )
    .expect("the built-in scene parses")
}
