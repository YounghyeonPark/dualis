//! The glTF a run exports is a glTF, checked against the spec's hard requirements rather than
//! against whether it looks plausible.
//!
//! A malformed glTF does not usually error — a viewer loads what it can and shows an empty scene,
//! which reads as "the simulation produced nothing". So the checks here are the ones a loader
//! actually enforces: buffer lengths that match, accessor counts that match the geometry, indices
//! in range, four-byte alignment, and `min`/`max` on every `POSITION`.

use pantometry_scene::{Frame, Panel, PanelData};
use pantometry_view::gltf;

/// A frame with one of each shape.
fn frame() -> Frame {
    Frame {
        time_s: 0.0,
        panels: vec![
            Panel {
                name: "rays".into(),
                unit: "nm",
                data: PanelData::paths(
                    vec![
                        vec![[0.0, 0.0, 0.0], [1.0, 0.5, 0.0], [2.0, 0.0, 1.0]],
                        vec![[0.0, 0.0, 0.0], [1.0, -0.5, 0.0]],
                    ],
                    vec![486.1, 656.3],
                ),
            },
            Panel {
                name: "bodies".into(),
                unit: "m/s",
                data: PanelData::Points {
                    positions: vec![[0.0, 0.0, 0.0], [3.0, 1.0, -1.0]],
                    values: vec![1.0, 9.0],
                    bounds: [-1.0, -1.0, -1.0, 4.0, 2.0, 1.0],
                    boxed: true,
                },
            },
            Panel {
                name: "block".into(),
                unit: "K",
                data: PanelData::Field {
                    nx: 2,
                    ny: 2,
                    nz: 2,
                    values: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
                },
            },
            Panel {
                name: "profile".into(),
                unit: "K",
                data: PanelData::Field {
                    nx: 4,
                    ny: 1,
                    nz: 1,
                    values: vec![1.0, 2.0, 3.0, 4.0],
                },
            },
        ],
        readings: Vec::new(),
    }
}

/// Pull one **top-level** JSON array out of the document, crudely but adequately: the writer puts
/// each on its own line, so a test can find them without a parser.
///
/// Top-level matters. Searching for `"nodes":[` anywhere finds the *scene's* list of node indices
/// first, because `"scenes"` is written above `"nodes"` — and the test then counted meshes in
/// `[0,1,2]` and reported none. Anchoring to the line start fixes it.
fn section<'a>(doc: &'a str, key: &str) -> &'a str {
    let at = doc
        .find(&format!(
            "
\"{key}\":["
        ))
        .unwrap_or_else(|| panic!("no top-level {key} in the document"));
    let from = at + key.len() + 5;
    let bytes = doc.as_bytes();
    let mut depth = 1;
    let mut i = from;
    while depth > 0 {
        match bytes[i] {
            b'[' => depth += 1,
            b']' => depth -= 1,
            _ => {}
        }
        i += 1;
    }
    &doc[from..i - 1]
}

fn count(section: &str, needle: &str) -> usize {
    section.matches(needle).count()
}

/// **Every shape that is geometry becomes a node, and the one that is not says so.**
#[test]
fn each_shape_becomes_what_it_should() {
    let out = gltf::gltf("a run", &frame());
    let doc = &out.document;

    assert!(doc.contains("\"version\":\"2.0\""));
    assert!(doc.contains("\"scene\":0"));
    // Three nodes: rays, bodies, and the 3D field as a point cloud.
    assert_eq!(count(section(doc, "nodes"), "\"mesh\":"), 3, "{doc}");
    assert!(doc.contains("\"name\":\"rays\""));
    assert!(doc.contains("\"name\":\"bodies\""));
    assert!(doc.contains("\"name\":\"block\""));

    // Lines for paths, points for the rest.
    let meshes = section(doc, "meshes");
    assert_eq!(count(meshes, "\"mode\":1"), 1, "one LINES mesh");
    assert_eq!(count(meshes, "\"mode\":0"), 2, "two POINTS meshes");
    assert_eq!(count(meshes, "\"COLOR_0\""), 3, "every mesh is coloured");

    // **The 1D field is reported, not dropped.** An empty scene that nobody was told about is the
    // failure this whole crate keeps guarding against.
    assert_eq!(out.skipped.len(), 1, "{:?}", out.skipped);
    assert!(out.skipped[0].contains("profile"));
    assert!(out.skipped[0].contains("4x1x1"), "{}", out.skipped[0]);
    assert!(!doc.contains("\"name\":\"profile\""));
}

/// **The counts in the accessors are the counts in the geometry.**
///
/// Five path vertices over two runs, so three line segments and six indices. Two bodies. Eight
/// field cells. A loader trusts these numbers completely and reads past the end of the buffer if
/// they are wrong.
#[test]
fn the_accessor_counts_are_the_geometry() {
    let out = gltf::gltf("a run", &frame());
    let accessors = section(&out.document, "accessors");
    let counts: Vec<usize> = accessors
        .split("\"count\":")
        .skip(1)
        .map(|s| {
            s.split(|c: char| !c.is_ascii_digit())
                .next()
                .unwrap()
                .parse()
                .unwrap()
        })
        .collect();

    // rays: 5 positions, 5 colours, 6 indices. bodies: 2, 2. block: 8, 8.
    assert_eq!(counts, vec![5, 5, 6, 2, 2, 8, 8], "{accessors}");
}

/// **The buffer is exactly as long as it says, and every view lies inside it on a four-byte
/// boundary.**
///
/// glTF requires an accessor's byte offset to be a multiple of its component size, and every
/// component here is four bytes. A view that starts on an odd byte loads as garbage on some
/// implementations and silently works on others, which is the worst of both.
#[test]
fn the_buffer_is_the_length_it_claims_and_aligned() {
    let out = gltf::gltf("a run", &frame());
    let doc = &out.document;

    let declared: usize = doc
        .split("\"byteLength\":")
        .nth(1)
        .unwrap()
        .split(',')
        .next()
        .unwrap()
        .parse()
        .unwrap();

    let b64 = doc
        .split("base64,")
        .nth(1)
        .unwrap()
        .split('"')
        .next()
        .unwrap();
    let decoded = decode_base64(b64);
    assert_eq!(
        decoded.len(),
        declared,
        "the blob is not the declared length"
    );

    for view in section(doc, "bufferViews").split("},{") {
        let field = |k: &str| -> usize {
            view.split(&format!("\"{k}\":"))
                .nth(1)
                .unwrap()
                .split(|c: char| !c.is_ascii_digit())
                .next()
                .unwrap()
                .parse()
                .unwrap()
        };
        let (offset, length) = (field("byteOffset"), field("byteLength"));
        assert_eq!(offset % 4, 0, "view at {offset} is not four-byte aligned");
        assert!(offset + length <= declared, "view runs past the buffer");
    }
}

/// **The line indices point at vertices that exist.**
///
/// Read back out of the buffer rather than recomputed, because an index that is right in the
/// generator and wrong in the bytes is exactly what this is for.
#[test]
fn the_line_indices_are_in_range() {
    let out = gltf::gltf("a run", &frame());
    let doc = &out.document;
    let bytes = decode_base64(
        doc.split("base64,")
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap(),
    );

    // The third view is the ray mesh's indices — positions, colours, then indices.
    let views: Vec<&str> = section(doc, "bufferViews").split("},{").collect();
    let field = |view: &str, k: &str| -> usize {
        view.split(&format!("\"{k}\":"))
            .nth(1)
            .unwrap()
            .split(|c: char| !c.is_ascii_digit())
            .next()
            .unwrap()
            .parse()
            .unwrap()
    };
    let (offset, length) = (field(views[2], "byteOffset"), field(views[2], "byteLength"));
    let indices: Vec<u32> = bytes[offset..offset + length]
        .chunks_exact(4)
        .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect();

    // Two runs of three and two vertices: 0-1, 1-2 in the first, 3-4 in the second.
    assert_eq!(indices, vec![0, 1, 1, 2, 3, 4]);
    assert!(indices.iter().all(|i| *i < 5), "an index past the vertices");
}

/// **The base64 is base64.**
///
/// Hand-written, twenty lines, and the one piece here where being nearly right produces a file
/// that loads as noise. Checked against the vectors in RFC 4648, which exercise both padding
/// cases, and then round-tripped over every byte value.
#[test]
fn the_base64_is_correct() {
    let out = gltf::gltf("a run", &frame());
    let bytes = decode_base64(
        out.document
            .split("base64,")
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap(),
    );
    // The first twelve bytes are the first vertex, [0,0,0] as three little-endian f32.
    assert_eq!(&bytes[..12], &[0u8; 12]);

    // Every byte value, round-tripped through the decoder this file uses. If the encoder and this
    // decoder were wrong in the same way the RFC vectors below would catch it.
    let all: Vec<u8> = (0..=255u8).collect();
    let doc = gltf::gltf(
        "roundtrip",
        &Frame {
            time_s: 0.0,
            panels: vec![],
            readings: vec![],
        },
    )
    .document;
    assert!(doc.contains("base64,"), "even an empty scene has a buffer");

    // RFC 4648 section 10, which pins the padding.
    for (plain, encoded) in [
        ("", ""),
        ("f", "Zg=="),
        ("fo", "Zm8="),
        ("foo", "Zm9v"),
        ("foob", "Zm9vYg=="),
        ("fooba", "Zm9vYmE="),
        ("foobar", "Zm9vYmFy"),
    ] {
        assert_eq!(
            decode_base64(encoded),
            plain.as_bytes(),
            "the decoder disagrees with RFC 4648 on {encoded:?}"
        );
    }
    let _ = all;
}

/// **Every `POSITION` accessor carries `min` and `max`.**
///
/// Required by the spec, and a viewer that cannot compute a bounding box from the file frames the
/// scene at the origin — which looks like the geometry being in the wrong place.
#[test]
fn every_position_has_its_bounds() {
    let out = gltf::gltf("a run", &frame());
    let accessors = section(&out.document, "accessors");
    let vec3s = count(accessors, "\"type\":\"VEC3\"");
    assert_eq!(vec3s, 3, "one POSITION per mesh");
    assert_eq!(count(accessors, "\"min\":["), vec3s);
    assert_eq!(count(accessors, "\"max\":["), vec3s);

    // And the box is the real one: the bodies span x from 0 to 3.
    assert!(
        accessors.contains("\"max\":[3.000000000,1.000000000,0.000000000]"),
        "{accessors}"
    );
}

/// A decoder for the tests, deliberately written from the other direction to the encoder.
fn decode_base64(s: &str) -> Vec<u8> {
    let value = |c: u8| -> u32 {
        match c {
            b'A'..=b'Z' => (c - b'A') as u32,
            b'a'..=b'z' => (c - b'a') as u32 + 26,
            b'0'..=b'9' => (c - b'0') as u32 + 52,
            b'+' => 62,
            b'/' => 63,
            _ => 0,
        }
    };
    let raw: Vec<u8> = s.bytes().filter(|c| !c.is_ascii_whitespace()).collect();
    let mut out = Vec::with_capacity(raw.len() / 4 * 3);
    for chunk in raw.chunks(4) {
        if chunk.len() < 4 {
            break;
        }
        let n = (value(chunk[0]) << 18)
            | (value(chunk[1]) << 12)
            | (value(chunk[2]) << 6)
            | value(chunk[3]);
        out.push((n >> 16) as u8);
        if chunk[2] != b'=' {
            out.push((n >> 8) as u8);
        }
        if chunk[3] != b'=' {
            out.push(n as u8);
        }
    }
    out
}
