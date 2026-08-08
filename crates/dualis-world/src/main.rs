//! Run a scene and draw it.
//!
//! ```sh
//! dualis-world                        # the built-in scene, checked, no output written
//! dualis-world scene.json out.svg     # a scene from a file, drawn
//! dualis-world --emit-default s.json  # write the built-in scene out to start from
//! ```

use dualis::prelude::ThermalNetwork;
use dualis_world::{render, Scene, World};

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
            Some((nx, ny)) => format!("{nx} x {ny}"),
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

    match out {
        Some(path) => {
            let svg = render::filmstrip(&world.scene().title, &frames, 6);
            // A zero-byte file used to be written and reported as "0 KiB" — which a legitimate
            // 937-byte strip also reports, so the one number on the line could not tell them
            // apart. Bytes now, and an empty picture is refused rather than saved.
            if svg.is_empty() {
                eprintln!(
                    "
nothing to draw: none of the {} domain(s) has a field or bodies,                      so {path} was not written",
                    world.scene().domains.len()
                );
                std::process::exit(1);
            }
            std::fs::write(path, &svg)?;
            println!("  wrote {path} ({} bytes)", svg.len());
        }
        None => println!("  give a second argument to write an SVG"),
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
