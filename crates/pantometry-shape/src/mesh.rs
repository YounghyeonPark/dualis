//! A bag of triangles, read from what a CAD tool exported, and what can be measured from one.

use glam::DVec3;
use pantometry_units::{Area, Length, LengthVec, Volume};

/// A vertex as its three coordinates' bit patterns, which is how [`Mesh::is_closed`] decides that two
/// triangles are talking about the same point. See that method for why it is exact.
type Vertex = (u64, u64, u64);

/// An edge as its two vertices, ordered so the same edge keys the same however it was wound.
type Edge = (Vertex, Vertex);

/// One triangle, in metres.
///
/// Wound counter-clockwise seen from outside, which is what makes [`Mesh::volume`] positive. STL states
/// a normal per facet as well, and it is **ignored**: exporters disagree about it often enough that the
/// winding is the more reliable of the two, and carrying a normal that might contradict the vertices
/// would mean choosing between them at every use.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Triangle {
    /// First vertex.
    pub a: DVec3,
    /// Second vertex.
    pub b: DVec3,
    /// Third vertex.
    pub c: DVec3,
}

impl Triangle {
    /// Twice the area as a vector along the normal — the cross product of two edges.
    pub fn normal_area(&self) -> DVec3 {
        (self.b - self.a).cross(self.c - self.a)
    }

    /// The triangle's area.
    pub fn area(&self) -> Area {
        Area::from_si(0.5 * self.normal_area().length())
    }
}

/// A closed surface as triangles, in metres.
///
/// # STL carries no topology, and that shapes what can be checked
///
/// The format is a flat list of facets, each with its three vertices written out in full. Two triangles
/// sharing an edge repeat those two vertices, and nothing in the file says they are the same points. So
/// [`Mesh::is_closed`] has to *infer* the topology by matching coordinates, and it matches them exactly —
/// bit for bit.
///
/// That is the right strictness for the thing being asked. A mesh whose shared vertices differ in the
/// last bit is not closed for any purpose that matters here: a ray can pass between the two triangles,
/// and the rasteriser below will see it. Reporting it as open is the true answer, and a tolerance would
/// turn a real defect into a silent one.
#[derive(Clone, Debug, Default)]
pub struct Mesh {
    triangles: Vec<Triangle>,
}

impl Mesh {
    /// A mesh from triangles, in metres.
    pub fn new(triangles: Vec<Triangle>) -> Mesh {
        Mesh { triangles }
    }

    /// Read an STL, binary or ASCII.
    ///
    /// # Which one it is, decided by arithmetic rather than by the first word
    ///
    /// The usual test is whether the file starts with `solid`, and it is wrong: a binary STL's header is
    /// eighty arbitrary bytes and plenty of exporters write `solid` into it. The reliable test is the
    /// length — a binary file is exactly `84 + 50n` bytes for its own declared `n`, and no ASCII file of
    /// that content is. This uses that, and falls back to ASCII.
    ///
    /// Lengths are taken as **millimetres**, because STL has no units and every mechanical CAD tool
    /// writes millimetres. That is a convention rather than a fact about the format, so it is stated
    /// here and nowhere else has to guess.
    pub fn from_stl(bytes: &[u8]) -> Result<Mesh, String> {
        if bytes.len() >= 84 {
            let count = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
            if let Some(expected) = count.checked_mul(50).and_then(|n| n.checked_add(84)) {
                if expected == bytes.len() {
                    return Mesh::from_binary_stl(bytes, count);
                }
            }
        }
        Mesh::from_ascii_stl(bytes)
    }

    fn from_binary_stl(bytes: &[u8], count: usize) -> Result<Mesh, String> {
        let mut triangles = Vec::with_capacity(count);
        let f32_at = |at: usize| -> f64 {
            f32::from_le_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]]) as f64
        };
        for n in 0..count {
            // 84-byte header and count, then 50 per facet: a normal this ignores, three vertices, and
            // two attribute bytes nothing standard uses.
            let base = 84 + 50 * n + 12;
            let v = |k: usize| {
                DVec3::new(
                    f32_at(base + 12 * k),
                    f32_at(base + 12 * k + 4),
                    f32_at(base + 12 * k + 8),
                ) * 1e-3
            };
            triangles.push(Triangle {
                a: v(0),
                b: v(1),
                c: v(2),
            });
        }
        Ok(Mesh { triangles })
    }

    fn from_ascii_stl(bytes: &[u8]) -> Result<Mesh, String> {
        let text = std::str::from_utf8(bytes)
            .map_err(|e| format!("not a binary STL by length, and not UTF-8 text either: {e}"))?;
        let mut vertices: Vec<DVec3> = Vec::new();
        let mut triangles = Vec::new();
        for (n, line) in text.lines().enumerate() {
            let mut word = line.split_whitespace();
            if word.next() != Some("vertex") {
                continue;
            }
            let mut coordinate = || -> Result<f64, String> {
                word.next()
                    .ok_or_else(|| format!("line {}: a vertex needs three numbers", n + 1))?
                    .parse::<f64>()
                    .map_err(|e| format!("line {}: {e}", n + 1))
            };
            let (x, y, z) = (coordinate()?, coordinate()?, coordinate()?);
            vertices.push(DVec3::new(x, y, z) * 1e-3);
            if vertices.len() == 3 {
                triangles.push(Triangle {
                    a: vertices[0],
                    b: vertices[1],
                    c: vertices[2],
                });
                vertices.clear();
            }
        }
        if !vertices.is_empty() {
            return Err(format!(
                "the file ends with {} vertices left over, so a facet is incomplete",
                vertices.len()
            ));
        }
        if triangles.is_empty() {
            return Err("no facets found; is this an STL?".to_string());
        }
        Ok(Mesh { triangles })
    }

    /// The triangles, in the order read.
    pub fn triangles(&self) -> &[Triangle] {
        &self.triangles
    }

    /// The axis-aligned bounds, as `(low, high)`.
    ///
    /// `None` for a mesh with no triangles, because the bounds of nothing are not a box at the origin.
    pub fn bounds(&self) -> Option<(LengthVec, LengthVec)> {
        let first = self.triangles.first()?;
        let mut low = first.a;
        let mut high = first.a;
        for t in &self.triangles {
            for v in [t.a, t.b, t.c] {
                low = low.min(v);
                high = high.max(v);
            }
        }
        Some((LengthVec::from_si(low), LengthVec::from_si(high)))
    }

    /// The enclosed volume, by the divergence theorem.
    ///
    /// `Σ a · (b × c) / 6` — the signed volume of the tetrahedron each triangle makes with the origin,
    /// summed. Everything outside the surface cancels exactly, so for a **closed** mesh this is the
    /// enclosed volume and it is exact to floating point, with no tolerance and no sampling.
    ///
    /// That exactness is what makes it the reference [`Loss::volume_error`](crate::Loss::volume_error) measures
    /// against: comparing a rasterisation to an analytic sphere would conflate two errors, the
    /// tessellation's and the grid's. Comparing it to the mesh's own volume isolates the one being
    /// measured.
    ///
    /// For an **open** mesh the number is meaningless rather than approximate — check [`Mesh::is_closed`].
    /// A negative volume means the winding is inside out, which is a real and common export defect.
    pub fn volume(&self) -> Volume {
        Volume::from_si(
            self.triangles
                .iter()
                .map(|t| t.a.dot(t.b.cross(t.c)) / 6.0)
                .sum::<f64>(),
        )
    }

    /// The total surface area.
    pub fn area(&self) -> Area {
        Area::from_si(self.triangles.iter().map(|t| t.area().to_si()).sum())
    }

    /// Whether every edge is shared by exactly two triangles.
    ///
    /// Matched on the vertices' **bit patterns**, for the reason in this type's documentation: STL stores
    /// no topology, so shared vertices are shared only if they were written identically, and a ray passes
    /// through a gap of one bit as readily as through a gap of one millimetre.
    ///
    /// With one exception, and it is not a tolerance. **Negative zero is folded onto zero**, because
    /// `-0.0` and `0.0` are the *same point* — the distance between them is nothing, and no ray passes
    /// between them. Their bit patterns differ, so a raw comparison reports a watertight mesh as open,
    /// and that is a false alarm rather than a strict answer. It arises constantly on anything symmetric
    /// about an axis, where one side's coordinate is a product that happened to carry a minus sign.
    ///
    /// A mesh that is not closed has no enclosed volume and cannot be rasterised by parity, so
    /// [`Voxels::of`](crate::Voxels::of) refuses one rather than producing a shape with holes in it.
    pub fn is_closed(&self) -> bool {
        // `+ 0.0` would do it in one operation, but this says what is meant and does not read as a
        // no-op that a later reader deletes.
        let zeroed = |c: f64| if c == 0.0 { 0.0 } else { c };
        let key = |v: DVec3| {
            (
                zeroed(v.x).to_bits(),
                zeroed(v.y).to_bits(),
                zeroed(v.z).to_bits(),
            )
        };
        let mut edges: std::collections::HashMap<Edge, i32> = std::collections::HashMap::new();
        for t in &self.triangles {
            for (p, q) in [(t.a, t.b), (t.b, t.c), (t.c, t.a)] {
                let (p, q) = (key(p), key(q));
                // Undirected, so the two triangles sharing an edge meet on the same key however they
                // wound it.
                let edge = if p <= q { (p, q) } else { (q, p) };
                *edges.entry(edge).or_insert(0) += 1;
            }
        }
        !edges.is_empty() && edges.values().all(|n| *n == 2)
    }

    /// How many triangles are smaller than one face of a cell of side `cell`.
    ///
    /// A feature the grid cannot hold, counted before anything is rasterised. It is not a proof that
    /// something is lost — a large flat face can be tessellated into small triangles and lose nothing —
    /// but a mesh where many facets are below the cell's own area is a mesh whose detail is finer than
    /// the grid, and that is worth being told before the run rather than after.
    pub fn triangles_below(&self, cell: Length) -> usize {
        let face = cell.to_si() * cell.to_si();
        self.triangles
            .iter()
            .filter(|t| t.area().to_si() < face)
            .count()
    }
}
