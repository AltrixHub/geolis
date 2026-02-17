use crate::geometry::pline::PlineVertex;
use crate::math::arc_2d::{arc_from_bulge, arc_point_at};

use super::self_intersect::Intersection;

/// A slice of a polyline between two intersection points.
#[derive(Debug, Clone)]
pub struct PlineSlice {
    pub vertices: Vec<PlineVertex>,
    /// Index of the starting intersection in the original intersection list.
    pub start_idx: usize,
    /// Index of the ending intersection in the original intersection list.
    pub end_idx: usize,
}

/// Slices a closed polyline at all intersection points, producing sub-paths.
///
/// Uses a BTreeMap-style approach: for each segment, collect the intersection
/// parameters in order, then split the polyline at those points.
///
/// Handles arc segments by computing sub-arc bulge values and interpolating
/// points on arcs.
#[must_use]
pub fn build(
    vertices: &[PlineVertex],
    n_segs: usize,
    intersections: &[Intersection],
) -> Vec<PlineSlice> {
    if intersections.is_empty() || vertices.is_empty() {
        return Vec::new();
    }

    let n = vertices.len();

    // Build a list of split points per segment: (segment_index, t, intersection_index).
    let mut splits: Vec<(usize, f64, usize)> = Vec::new();
    for (idx, ix) in intersections.iter().enumerate() {
        splits.push((ix.seg_i, ix.t_i, idx));
        splits.push((ix.seg_j, ix.t_j, idx));
    }
    splits.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });

    // Walk around the polyline, emitting slices between consecutive split points.
    let mut slices = Vec::new();
    let total_splits = splits.len();

    for sp_idx in 0..total_splits {
        let (seg_start, t_start, ix_start) = splits[sp_idx];
        let (seg_end, t_end, ix_end) = splits[(sp_idx + 1) % total_splits];

        let verts = build_slice_verts(vertices, n, n_segs, seg_start, t_start, seg_end, t_end);

        if verts.len() >= 2 {
            slices.push(PlineSlice {
                vertices: verts,
                start_idx: ix_start,
                end_idx: ix_end,
            });
        }
    }

    slices
}

/// Builds the vertices for a single slice from `(seg_start, t_start)` to `(seg_end, t_end)`.
///
/// Correctly handles arc segments by:
/// - Interpolating points on arcs (not linear)
/// - Computing sub-arc bulge for partial segments
/// - Preserving original bulge for full interior segments
#[allow(clippy::too_many_arguments)]
fn build_slice_verts(
    vertices: &[PlineVertex],
    n: usize,
    n_segs: usize,
    seg_start: usize,
    t_start: f64,
    seg_end: usize,
    t_end: f64,
) -> Vec<PlineVertex> {
    let mut verts = Vec::new();

    if seg_start == seg_end {
        // Both split points on the same segment: single sub-segment.
        let start_pos = point_on_segment(vertices, n, seg_start, t_start);
        let end_pos = point_on_segment(vertices, n, seg_end, t_end);
        let bulge = sub_bulge(vertices[seg_start].bulge, t_start, t_end);
        verts.push(PlineVertex::new(start_pos.0, start_pos.1, bulge));
        verts.push(PlineVertex::line(end_pos.0, end_pos.1));
        return verts;
    }

    // Start vertex: connects to the end of seg_start (or next full vertex).
    let start_pos = point_on_segment(vertices, n, seg_start, t_start);
    let start_bulge = sub_bulge(vertices[seg_start].bulge, t_start, 1.0);
    verts.push(PlineVertex::new(start_pos.0, start_pos.1, start_bulge));

    // Walk full interior segments from seg_start+1 to seg_end-1.
    let mut seg = (seg_start + 1) % n_segs;
    while seg != seg_end {
        let vi = seg;
        let v = &vertices[vi];
        // Full segment: preserve original bulge.
        verts.push(PlineVertex::new(v.x, v.y, v.bulge));
        seg = (seg + 1) % n_segs;
    }

    // Last original vertex before the end split point.
    let v_end_start = &vertices[seg_end];
    let end_bulge = sub_bulge(v_end_start.bulge, 0.0, t_end);
    verts.push(PlineVertex::new(v_end_start.x, v_end_start.y, end_bulge));

    // End point.
    let end_pos = point_on_segment(vertices, n, seg_end, t_end);
    verts.push(PlineVertex::line(end_pos.0, end_pos.1));

    verts
}

/// Computes the position of a point at parameter `t` on a segment.
///
/// For line segments (bulge=0): linear interpolation.
/// For arc segments (bulge≠0): point on the arc.
fn point_on_segment(
    vertices: &[PlineVertex],
    n: usize,
    seg_idx: usize,
    t: f64,
) -> (f64, f64) {
    let v0 = &vertices[seg_idx];
    let v1 = &vertices[(seg_idx + 1) % n];

    if v0.bulge.abs() < 1e-12 {
        // Line segment.
        (v0.x + t * (v1.x - v0.x), v0.y + t * (v1.y - v0.y))
    } else {
        // Arc segment.
        let (cx, cy, r, sa, sw) = arc_from_bulge(v0.x, v0.y, v1.x, v1.y, v0.bulge);
        arc_point_at(cx, cy, r, sa, sw, t)
    }
}

/// Computes the bulge for a sub-arc spanning parameter range `[t_start, t_end]`.
///
/// For line segments (bulge ≈ 0), returns 0.
/// For arc segments, computes `tan(sub_sweep / 4)` where `sub_sweep = sweep * (t_end - t_start)`.
fn sub_bulge(original_bulge: f64, t_start: f64, t_end: f64) -> f64 {
    if original_bulge.abs() < 1e-12 {
        return 0.0;
    }
    let sweep = 4.0 * original_bulge.atan();
    let sub_sweep = sweep * (t_end - t_start);
    (sub_sweep / 4.0).tan()
}
