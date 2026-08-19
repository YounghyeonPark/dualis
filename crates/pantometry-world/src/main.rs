//! Run a scene and draw it.
//!
//! The asset is chosen by the output file's extension, because a run has several shapes and
//! only one of them is a picture.
//!
//! ```sh
//! pantometry-world                        # the built-in scene, checked, nothing written
//! pantometry-world scene.json out.html    # a report: a view per domain, chosen from its shape
//! pantometry-world scene.json out.svg     # a filmstrip: every frame, one page, still
//! pantometry-world scene.json out.csv     # every domain's scalars over time, one row per frame
//! pantometry-world scene.json out.json    # the frames themselves — fields, bodies and readings
//! pantometry-world --check s.json         # does it parse and build, without running it
//! pantometry-world --emit-default s.json  # write the built-in scene out to start from
//! pantometry-world verify s.json          # the battery: margins, determinism, both sweeps
//! pantometry-world verify s.json --deep   # and a third run per sweep, measuring the order
//! ```
//!
//! `.html` is the one to reach for if you do not already know how you want this drawn. It picks
//! a view from each domain's *shape* — a profile for a 1D field, a heatmap for a 2D one, a
//! rotatable scene for bodies, a line chart for scalars — and opens in a browser with nothing
//! installed. `.csv` is the one that changed what this crate can say. Eight of the fourteen shipped scenes
//! have a domain the filmstrip cannot draw, and for several of them the scalar *is* the result —
//! scene 13's whole subject is a winding whose resistance follows its own temperature, and it
//! drew nothing at all.

use pantometry::prelude::ThermalNetwork;
// The views are a library now — `pantometry-view`, above `pantometry-scene`, above the kernel. They
// were modules in this binary, which is `publish = false`, so a consumer who could state a
// simulation and run it could not draw one.
use pantometry::view::{html, readings_csv, svg as filmstrip, to_json};
use pantometry_world::{Scene, World};

/// Read a file, saying **which** file when it cannot be read.
///
/// `std::fs::read_to_string(path)?` propagates an `io::Error` that names the reason and not the
/// file, so a mistyped or missing path exits with
/// `Os { code: 2, kind: NotFound, message: "The system cannot find the file specified." }` and
/// nothing in it to act on. That is what a first run of this binary actually printed, and the
/// parse path one line further on had been reporting `file:line:column` the whole time — the two
/// halves of reading a file disagreed about whether the reader deserved to know its name.
///
/// A missing file also gets the sentence that fixes it, because "not found" is a complete
/// diagnosis and an incomplete answer: nothing in this repository is called `scene.json` until
/// somebody writes one.
fn read(path: &str) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| {
        if e.kind() == std::io::ErrorKind::NotFound {
            format!(
                "{path}: {e}\n  write one to start from with `--emit-default {path}`, or run a \
                 shipped scene from crates/pantometry-world/scenes/"
            )
        } else {
            format!("{path}: {e}")
        }
    })
}

/// Write a file, saying which file when it cannot be written. The counterpart to [`read`], and
/// here for the same reason: a full disk or a read-only directory should not report itself as a
/// bare `Os { code: 13 }`.
fn write(path: &str, bytes: &str) -> Result<(), String> {
    std::fs::write(path, bytes).map_err(|e| format!("{path}: {e}"))
}

/// Print what went wrong and exit non-zero, rather than letting the runtime `Debug`-print it.
///
/// Returning `Result` from `main` is idiomatic and prints the **Debug** form: a `String` error
/// comes out wrapped in quotes with its newlines escaped, so the hint this binary attaches to a
/// missing file arrived as a literal backslash-n. Every other error path here already printed
/// and exited; this makes the last one agree with them.
fn main() {
    if let Err(e) = run() {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();

    // Validate without running. What an editor needs while somebody is typing: the same checks
    // `World::build` makes — the format version, the domain names, a `tracks` that points at a
    // node the scene defines — and none of the seconds a run costs.
    //
    // A separate path rather than a flag on the run, because "did this file parse" and "what did
    // this file do" are different questions and an editor asks the first one constantly.
    if args.first().map(String::as_str) == Some("--check") {
        let path = args.get(1).ok_or("--check needs a path")?;
        let text = read(path)?;
        let scene: Scene = match serde_json::from_str(&text) {
            Ok(s) => s,
            Err(e) => {
                // Line and column, because that is what an editor puts a squiggle under.
                eprintln!("{path}:{}:{}: {e}", e.line(), e.column());
                std::process::exit(1);
            }
        };
        match World::build(scene.clone()) {
            Ok(world) => {
                println!(
                    "{path}: format {}, {} domain(s), {:.3} s in {} frames",
                    scene.format,
                    scene.domains.len(),
                    scene.duration_s,
                    scene.frames
                );
                for spec in &world.scene().domains {
                    let placement = spec.placement();
                    let shape = match placement.extent {
                        Some(e) => {
                            let (nx, ny, nz) = e.samples;
                            format!("a field sampled {nx} x {ny} x {nz}")
                        }
                        None => "no field".to_string(),
                    };
                    println!("  {:<14} {shape}", spec.name());
                }
                // The dismissals: a stated condition a domain correctly ignored, with the
                // measurement that earns it. Printed because a dismissal nobody can see is
                // the silence it exists to replace.
                for note in world.notes() {
                    println!("  note: {note}");
                }
                return Ok(());
            }
            Err(why) => {
                eprintln!("{path}: {why}");
                std::process::exit(1);
            }
        }
    }

    // Choosing a grid, which `ARCHITECTURE.md` left open as the third assembly gap because it
    // "would be the first place in the workspace that guesses". It does not guess: it rasterises
    // at every candidate and prints what each one cost. See `pantometry_world::fit`.
    if args.first().map(String::as_str) == Some("fit") {
        let mut paths = Vec::new();
        let mut budget = 2_000_000usize;
        let mut material = "aluminium".to_string();
        let mut rest = args[1..].iter();
        while let Some(a) = rest.next() {
            match a.as_str() {
                "--cells" => {
                    budget = rest
                        .next()
                        .ok_or("--cells needs a number")?
                        .parse()
                        .map_err(|e| format!("--cells: {e}"))?;
                }
                "--material" => material = rest.next().ok_or("--material needs a name")?.clone(),
                other if other.starts_with("--") => {
                    return Err(format!(
                        "fit takes STL paths, --cells N and --material NAME, not {other}"
                    )
                    .into());
                }
                other => paths.push(other.to_string()),
            }
        }
        if paths.is_empty() {
            return Err("fit needs at least one STL".into());
        }
        let mut parts = Vec::new();
        for path in &paths {
            let bytes = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?;
            let mesh =
                pantometry::shape::Mesh::from_stl(&bytes).map_err(|e| format!("{path}: {e}"))?;
            // The name a scene will use is the file name, not the path it was read from, so a
            // fragment printed here is one a browser can use as well as a terminal.
            let name = std::path::Path::new(path)
                .file_name()
                .map_or_else(|| path.clone(), |n| n.to_string_lossy().into_owned());
            parts.push((name, mesh));
        }
        let fit = pantometry_world::fit::propose(&parts, budget)?;
        print!("{}", fit.render());
        match fit.recommended(0.5) {
            Some(c) => {
                println!(
                    "  the coarsest grid that holds every part with under half its volume in \
                     boundary cells:\n"
                );
                println!("{}", fit.scene_fragment(c, &material));
            }
            None => println!(
                "  no row holds every part under half boundary — raise --cells, or the parts have \
                 detail this grid cannot carry"
            ),
        }
        return Ok(());
    }

    // The measurements a passing audit does not make: margins, determinism, and what moves
    // when the coupling window halves or every grid refines. See `pantometry_world::verify` for
    // what each is and which shipped error it exists because of.
    if args.first().map(String::as_str) == Some("verify") {
        let path = args.get(1).ok_or("verify needs a path")?;
        // The whole tail is matched, not just the third argument — `verify s.json --deep
        // --whatever` silently ignoring the rest would run something other than what was
        // typed and report it as verified.
        let deep = match &args[2..] {
            [] => false,
            [d] if d == "--deep" => true,
            rest => {
                eprintln!("verify takes a scene and optionally --deep, not {rest:?}");
                std::process::exit(1);
            }
        };
        let text = read(path)?;
        let scene: Scene = serde_json::from_str(&text).map_err(|e| format!("{path}: {e}"))?;
        println!("{}", scene.title);
        let battery = match pantometry_world::verify::verify(&scene, deep) {
            Ok(b) => b,
            Err(why) => {
                eprintln!("{path}: {why}");
                std::process::exit(1);
            }
        };
        print!("{}", battery.render());
        if !battery.findings.is_empty() {
            std::process::exit(1);
        }
        return Ok(());
    }

    if args.first().map(String::as_str) == Some("--emit-default") {
        let path = args.get(1).ok_or("--emit-default needs a path")?;
        write(
            path,
            &(serde_json::to_string_pretty(&default_scene())? + "\n"),
        )?;
        println!("wrote {path}");
        return Ok(());
    }

    let scene = match args.first() {
        Some(path) => {
            serde_json::from_str::<Scene>(&read(path)?).map_err(|e| format!("{path}: {e}"))?
        }
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
                    write(path, &svg)?;
                    println!("  wrote {path} ({} bytes, filmstrip)", svg.len());
                }
                "csv" => {
                    let csv = readings_csv(&frames);
                    write(path, &csv)?;
                    println!(
                        "  wrote {path} ({} bytes, {} columns over {} rows)",
                        csv.len(),
                        frames.first().map_or(0, |f| f.readings.len()),
                        frames.len()
                    );
                }
                "html" => {
                    let page = html(world.scene().title.as_str(), &frames);
                    write(path, &page)?;
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
                    let out = pantometry::view::gltf(world.scene().title.as_str(), last);
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
                    write(path, &out.document)?;
                    println!(
                        "  wrote {path} ({} bytes, geometry for Blender, three.js or a USD tool)",
                        out.document.len()
                    );
                }
                "json" => {
                    let json = to_json(world.scene().title.as_str(), &frames);
                    write(path, &json)?;
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
