//! Run a scene and draw it.
//!
//! ```sh
//! dualis-world                        # the built-in scene, checked, no output written
//! dualis-world scene.json out.svg     # a scene from a file, drawn
//! dualis-world --emit-default s.json  # write the built-in scene out to start from
//! ```

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

    for panel in &frames[frames.len() - 1].panels {
        let peak = panel.values.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        println!(
            "  {:<14} {:>3} x {:<3} peak |{}| = {:.4}",
            panel.name, panel.nx, panel.ny, panel.unit, peak
        );
    }

    match out {
        Some(path) => {
            let svg = render::filmstrip(&world.scene().title, &frames, 6);
            std::fs::write(path, &svg)?;
            println!("  wrote {path} ({} KiB)", svg.len() / 1024);
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
      "cells_across": 61, "mode": [1, 1], "amplitude_pa": 1.0 }
  ]
}"#,
    )
    .expect("the built-in scene parses")
}
