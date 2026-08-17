//! A frame as glTF 2.0, so somebody else's renderer can draw it.
//!
//! # Why export rather than render
//!
//! Blender, three.js, Omniverse, Isaac Sim, macOS Quick Look and every USD pipeline read glTF.
//! Competing with those on rendering is a losing use of effort; what this workspace has that they
//! do not is physics that is audited, deterministic and checked against closed forms. So the move
//! is to make a result *reachable from* them.
//!
//! # It costs nothing
//!
//! glTF is JSON with the binary data base64'd into a `data:` URI, and this crate already writes
//! JSON by hand on purpose. So there is no new dependency, no encoder, and no build step — the
//! same reason SVG was chosen over a raster format. `dualis-view` still has exactly one
//! dependency, and it is `dualis-scene`.
//!
//! # One frame, and why not the run
//!
//! [`gltf`] takes a single [`Frame`]. glTF animates **node transforms and morph targets** — a
//! thing moving or a mesh deforming between fixed vertex counts. A field whose values change, a
//! ray bundle that retraces, a body count that is constant only by luck: none of those is either.
//! Encoding a run as an animation would mean choosing a lie about what is moving.
//!
//! Exporting the frame a caller asks for is the honest shape. For a whole run there is
//! [`to_json`](crate::to_json), which carries every frame and is what the native viewer reads.
//!
//! # What each shape becomes
//!
//! | panel | glTF |
//! | --- | --- |
//! | paths | `LINES`, one segment per adjacent pair, coloured per path |
//! | points | `POINTS`, coloured per body |
//! | a 3D field | `POINTS` at the cell centres, coloured by value |
//! | a 1D or 2D field | nothing — a line of samples is a graph, not geometry |
//!
//! The last row is a refusal rather than an omission: a profile drawn as a row of dots in a 3D
//! scene is a picture of the sampling, not of the physics, and [`gltf`] reports what it left out.

use dualis_scene::{Frame, Panel, PanelData};

/// One frame as a glTF 2.0 document, and what was left out of it.
pub struct Exported {
    /// The document. Write it to a `.gltf` file; the binary data is embedded, so there is no
    /// sidecar to lose.
    pub document: String,
    /// Panels that produced no geometry, and why — so a caller is never left wondering whether an
    /// empty scene is the physics or the exporter.
    pub skipped: Vec<String>,
}

/// Colour ramp shared with every other view here, so one run looks like one run wherever it is
/// opened.
fn ramp(t: f64) -> [f32; 4] {
    let t = t.clamp(0.0, 1.0);
    let stops = [
        [24.0, 52.0, 110.0],
        [61.0, 139.0, 255.0],
        [255.0, 194.0, 71.0],
        [255.0, 82.0, 33.0],
    ];
    let (a, b, u) = if t < 0.34 {
        (stops[0], stops[1], t / 0.34)
    } else if t < 0.67 {
        (stops[1], stops[2], (t - 0.34) / 0.33)
    } else {
        (stops[2], stops[3], (t - 0.67) / 0.33)
    };
    [
        ((a[0] + (b[0] - a[0]) * u) / 255.0) as f32,
        ((a[1] + (b[1] - a[1]) * u) / 255.0) as f32,
        ((a[2] + (b[2] - a[2]) * u) / 255.0) as f32,
        1.0,
    ]
}

/// One mesh under construction: positions, per-vertex colours, and line indices if it has any.
struct Mesh {
    name: String,
    positions: Vec<[f32; 3]>,
    colours: Vec<[f32; 4]>,
    /// Empty for a point cloud.
    lines: Vec<u32>,
}

/// Export one frame.
///
/// Every panel that is geometry becomes a node under one scene, named after its domain, so a
/// reader opening the file sees the same names the run reported.
pub fn gltf(title: &str, frame: &Frame) -> Exported {
    let mut meshes = Vec::new();
    let mut skipped = Vec::new();

    for panel in &frame.panels {
        match &panel.data {
            PanelData::Paths {
                vertices,
                starts,
                values,
                ..
            } => {
                let (lo, hi) = span(values);
                let mut mesh = Mesh {
                    name: panel.name.clone(),
                    positions: Vec::with_capacity(vertices.len()),
                    colours: Vec::with_capacity(vertices.len()),
                    lines: Vec::new(),
                };
                for (k, start) in starts.iter().enumerate() {
                    let from = *start;
                    let to = starts.get(k + 1).copied().unwrap_or(vertices.len());
                    let colour = ramp((values.get(k).copied().unwrap_or(lo) - lo) / (hi - lo));
                    let base = mesh.positions.len() as u32;
                    for v in &vertices[from..to] {
                        mesh.positions.push([v[0] as f32, v[1] as f32, v[2] as f32]);
                        mesh.colours.push(colour);
                    }
                    // One segment per adjacent pair. Indices rather than duplicated vertices, so a
                    // ray of n points costs n positions and not 2(n-1).
                    for i in 0..(to - from).saturating_sub(1) as u32 {
                        mesh.lines.push(base + i);
                        mesh.lines.push(base + i + 1);
                    }
                }
                if !mesh.lines.is_empty() {
                    meshes.push(mesh);
                }
            }
            PanelData::Points {
                positions, values, ..
            } => {
                let (lo, hi) = span(values);
                let mesh = Mesh {
                    name: panel.name.clone(),
                    positions: positions
                        .iter()
                        .map(|p| [p[0] as f32, p[1] as f32, p[2] as f32])
                        .collect(),
                    colours: values.iter().map(|v| ramp((v - lo) / (hi - lo))).collect(),
                    lines: Vec::new(),
                };
                if !mesh.positions.is_empty() {
                    meshes.push(mesh);
                }
            }
            PanelData::Field { nx, ny, nz, values } if *nz > 1 => {
                // Cell centres in grid units. The extent the field was sampled over is not in the
                // frame — `Placement` held it and did not travel — so this is the grid, and a
                // reader scaling the node is doing what the writer could not.
                let (lo, hi) = span(values);
                let mut mesh = Mesh {
                    name: panel.name.clone(),
                    positions: Vec::with_capacity(values.len()),
                    colours: Vec::with_capacity(values.len()),
                    lines: Vec::new(),
                };
                for k in 0..*nz {
                    for j in 0..*ny {
                        for i in 0..*nx {
                            let v = values[i + nx * (j + ny * k)];
                            // **A sample that is not a number is not a point.** This is the one
                            // output that leaves the workspace, and a grid with a clearance in it
                            // was exported as a solid brick: every empty cell arrived in Blender
                            // as a vertex with a colour. Skipped rather than drawn dark, because
                            // dark is a temperature and absence is not.
                            if !v.is_finite() {
                                continue;
                            }
                            mesh.positions
                                .push([i as f32 + 0.5, j as f32 + 0.5, k as f32 + 0.5]);
                            mesh.colours.push(ramp((v - lo) / (hi - lo)));
                        }
                    }
                }
                meshes.push(mesh);
            }
            PanelData::Field { nx, ny, nz, .. } => {
                skipped.push(format!(
                    "{} is a {nx}x{ny}x{nz} field: a line or a plane of samples is a graph, not \
                     geometry, so it is not in the scene",
                    panel.name
                ));
            }
        }
    }

    Exported {
        document: document(title, &meshes),
        skipped,
    }
}

fn span(values: &[f64]) -> (f64, f64) {
    let lo = values.iter().copied().fold(f64::MAX, f64::min);
    let hi = values.iter().copied().fold(f64::MIN, f64::max);
    if hi > lo {
        (lo, hi)
    } else {
        (lo, lo + 1.0)
    }
}

/// Assemble the document.
///
/// The binary layout is positions, then colours, then indices, per mesh, each padded to a
/// four-byte boundary — glTF requires an accessor's offset to be a multiple of its component
/// size, and every component here is four bytes wide.
fn document(title: &str, meshes: &[Mesh]) -> String {
    let mut blob: Vec<u8> = Vec::new();
    let mut views = Vec::new();
    let mut accessors = Vec::new();
    let mut nodes = Vec::new();
    let mut mesh_json = Vec::new();

    for mesh in meshes {
        let pos_view = push_view(&mut blob, &mut views, floats3(&mesh.positions), 34962);
        let (min, max) = extent(&mesh.positions);
        let pos_accessor = accessors.len();
        accessors.push(format!(
            "{{\"bufferView\":{pos_view},\"componentType\":5126,\"count\":{},\"type\":\"VEC3\",\
             \"min\":[{:.9},{:.9},{:.9}],\"max\":[{:.9},{:.9},{:.9}]}}",
            mesh.positions.len(),
            min[0],
            min[1],
            min[2],
            max[0],
            max[1],
            max[2]
        ));

        let col_view = push_view(&mut blob, &mut views, floats4(&mesh.colours), 34962);
        let col_accessor = accessors.len();
        accessors.push(format!(
            "{{\"bufferView\":{col_view},\"componentType\":5126,\"count\":{},\"type\":\"VEC4\"}}",
            mesh.colours.len()
        ));

        let (mode, index_accessor) = if mesh.lines.is_empty() {
            (0, None)
        } else {
            let view = push_view(&mut blob, &mut views, u32s(&mesh.lines), 34963);
            let k = accessors.len();
            accessors.push(format!(
                "{{\"bufferView\":{view},\"componentType\":5125,\"count\":{},\"type\":\"SCALAR\"}}",
                mesh.lines.len()
            ));
            (1, Some(k))
        };

        let indices = index_accessor
            .map(|k| format!(",\"indices\":{k}"))
            .unwrap_or_default();
        mesh_json.push(format!(
            "{{\"name\":{},\"primitives\":[{{\"attributes\":{{\"POSITION\":{pos_accessor},\
             \"COLOR_0\":{col_accessor}}}{indices},\"mode\":{mode},\"material\":0}}]}}",
            quote(&mesh.name)
        ));
        nodes.push(format!(
            "{{\"name\":{},\"mesh\":{}}}",
            quote(&mesh.name),
            mesh_json.len() - 1
        ));
    }

    let node_list: Vec<String> = (0..nodes.len()).map(|i| i.to_string()).collect();
    format!(
        "{{\n\
         \"asset\":{{\"version\":\"2.0\",\"generator\":{}}},\n\
         \"scene\":0,\n\
         \"scenes\":[{{\"name\":{},\"nodes\":[{}]}}],\n\
         \"nodes\":[{}],\n\
         \"meshes\":[{}],\n\
         \"materials\":[{{\"name\":\"vertex colour\",\"pbrMetallicRoughness\":\
         {{\"baseColorFactor\":[1,1,1,1],\"metallicFactor\":0,\"roughnessFactor\":1}},\
         \"doubleSided\":true}}],\n\
         \"buffers\":[{{\"byteLength\":{},\"uri\":\"data:application/octet-stream;base64,{}\"}}],\n\
         \"bufferViews\":[{}],\n\
         \"accessors\":[{}]\n\
         }}\n",
        quote(&format!("dualis-view {}", env!("CARGO_PKG_VERSION"))),
        quote(title),
        node_list.join(","),
        nodes.join(","),
        mesh_json.join(","),
        blob.len(),
        base64(&blob),
        views.join(","),
        accessors.join(",")
    )
}

/// Append bytes to the blob, pad to four, and record a view over them.
fn push_view(blob: &mut Vec<u8>, views: &mut Vec<String>, bytes: Vec<u8>, target: u32) -> usize {
    while blob.len() % 4 != 0 {
        blob.push(0);
    }
    let offset = blob.len();
    let length = bytes.len();
    blob.extend(bytes);
    views.push(format!(
        "{{\"buffer\":0,\"byteOffset\":{offset},\"byteLength\":{length},\"target\":{target}}}"
    ));
    views.len() - 1
}

fn floats3(v: &[[f32; 3]]) -> Vec<u8> {
    v.iter()
        .flat_map(|p| p.iter().flat_map(|f| f.to_le_bytes()))
        .collect()
}

fn floats4(v: &[[f32; 4]]) -> Vec<u8> {
    v.iter()
        .flat_map(|p| p.iter().flat_map(|f| f.to_le_bytes()))
        .collect()
}

fn u32s(v: &[u32]) -> Vec<u8> {
    v.iter().flat_map(|i| i.to_le_bytes()).collect()
}

/// The tightest box round a set of positions. **Required** by the spec on a `POSITION` accessor,
/// and a viewer that cannot frame a scene is the usual symptom of leaving it out.
fn extent(positions: &[[f32; 3]]) -> ([f32; 3], [f32; 3]) {
    let mut min = [f32::MAX; 3];
    let mut max = [f32::MIN; 3];
    for p in positions {
        for a in 0..3 {
            min[a] = min[a].min(p[a]);
            max[a] = max[a].max(p[a]);
        }
    }
    if positions.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    (min, max)
}

/// Standard base64, written out rather than depended on. Twenty lines against a crate.
fn base64(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let b = [
            chunk[0],
            chunk.get(1).copied().unwrap_or(0),
            chunk.get(2).copied().unwrap_or(0),
        ];
        let n = ((b[0] as u32) << 16) | ((b[1] as u32) << 8) | b[2] as u32;
        out.push(ALPHABET[(n >> 18) as usize & 63] as char);
        out.push(ALPHABET[(n >> 12) as usize & 63] as char);
        out.push(if chunk.len() > 1 {
            ALPHABET[(n >> 6) as usize & 63] as char
        } else {
            '='
        });
        out.push(if chunk.len() > 2 {
            ALPHABET[n as usize & 63] as char
        } else {
            '='
        });
    }
    out
}

fn quote(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if c.is_control() => out.push(' '),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// One panel's geometry, for a caller that wants a single domain rather than the frame.
///
/// The same document, with one node in it.
pub fn gltf_panel(title: &str, panel: &Panel) -> Exported {
    gltf(
        title,
        &Frame {
            time_s: 0.0,
            panels: vec![panel.clone()],
            readings: Vec::new(),
        },
    )
}
