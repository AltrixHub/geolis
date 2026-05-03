//! 2D polygon boolean-union via directed half-edge arrangement.
//!
//! [`union_all_with_holes`] is the public entry point. It returns the
//! boolean-union outline of one or more `PolygonWithHoles` inputs.
//! Output edges are exactly the boundary that separates **filled** from
//! **empty** points, where "filled" follows the OR-of-PWH-filled rule
//! (a point is filled iff there exists at least one input PWH whose
//! outer contains the point AND none of whose holes contains it).
//!
//! ## Algorithm
//! 1. **Collect raw segments** from every ring of every input PWH.
//! 2. **Arrangement split** every segment at every transverse crossing,
//!    T-junction (an endpoint of one segment falling on the interior of
//!    another within `WALL_EPS`), and collinear-overlap endpoint with
//!    every other segment.
//! 3. **Vertex snap** by grid quantization (cell size = `WALL_EPS`) plus
//!    a bounded cross-cell merge that preserves cluster diameter
//!    `≤ 2·WALL_EPS`. Cluster representatives are sorted-mean lex; class
//!    ids are assigned in the lex order of representatives so the result
//!    is order-independent.
//! 4. **Canonicalize** sub-edges (sorted-endpoint-tuple dedup, lex sort).
//! 5. **Classify directed half-edges** by bilateral perpendicular-ε
//!    sampling against the OR-of-PWH-filled oracle. Keep half-edges
//!    where the LEFT side is `Filled` AND the RIGHT side is `Empty`.
//!    Adaptive ε retry (3 halvings) on `AmbiguousOnBoundary`; if
//!    classification remains ambiguous after retries, the function
//!    returns `Err(OperationError::Failed)`.
//! 6. **Face-walk trace** with the polar-angle-Δ successor rule:
//!    incoming half-edge `u→v`, outgoing successor `v→w` chosen as the
//!    one whose clockwise Δ from the reverse direction `θ_in` is
//!    smallest (the self-reverse half-edge at Δ=0 is treated as 2π).
//!
//! Determinism: outputs are topologically identical (and float-equivalent
//! within `WALL_EPS` precision) regardless of input order.

use crate::error::{OperationError, Result};
use crate::math::distance_2d::point_to_segment_dist;
use std::collections::HashMap;

/// Single epsilon for all geometric decisions in the wall outline pipeline.
pub const WALL_EPS: f64 = 1e-6;
pub const WALL_EPS_SQ: f64 = WALL_EPS * WALL_EPS;

/// A simple 2D polygon (closed, no holes).
pub type Polygon = Vec<(f64, f64)>;

/// A polygon with optional holes.
pub struct PolygonWithHoles {
    pub outer: Polygon,
    pub holes: Vec<Polygon>,
}

/// Result of a polygon union: closed boundary loops.
pub struct UnionResult {
    pub boundaries: Vec<Polygon>,
}

/// Three-valued classification of a point relative to a single ring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PointClass {
    Inside,
    Outside,
    Boundary,
}

/// Identifies which ring of an input PWH a `Boundary` classification touched.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum BoundaryRef {
    Outer,
    Hole(usize),
}

/// Three-valued OR-of-PWH-filled classification.
///
/// `AmbiguousOnBoundary` is NOT folded into Filled or Empty — it forces
/// the caller to shrink ε and re-sample. The `touched` field is purely
/// diagnostic: when retries exhaust, the panic / `Err` message includes
/// it so the caller can pinpoint which input PWH boundary the sample hit.
#[derive(Debug, Clone)]
pub(crate) enum FilledClass {
    Filled,
    Empty,
    AmbiguousOnBoundary {
        #[allow(dead_code, reason = "diagnostic via Debug formatting")]
        touched: Vec<(usize, BoundaryRef)>,
    },
}

/// OR-of-PWH-filled rule.
///
/// Returns `Filled` if there exists an input PWH whose outer contains
/// `p` strictly Inside AND every hole has `p` strictly Outside. Returns
/// `AmbiguousOnBoundary` if any per-ring classification call returned
/// `Boundary`. Otherwise `Empty`.
pub(crate) fn classify_filled(p: (f64, f64), inputs: &[PolygonWithHoles]) -> FilledClass {
    let mut touched: Vec<(usize, BoundaryRef)> = Vec::new();
    let mut any_filled = false;
    for (idx, pwh) in inputs.iter().enumerate() {
        let outer_class = point_in_polygon_class(p, &pwh.outer);
        if outer_class == PointClass::Boundary {
            touched.push((idx, BoundaryRef::Outer));
        }
        let hole_classes: Vec<PointClass> = pwh
            .holes
            .iter()
            .map(|h| point_in_polygon_class(p, h))
            .collect();
        for (hi, hc) in hole_classes.iter().enumerate() {
            if *hc == PointClass::Boundary {
                touched.push((idx, BoundaryRef::Hole(hi)));
            }
        }
        if outer_class == PointClass::Inside
            && hole_classes.iter().all(|c| *c == PointClass::Outside)
        {
            any_filled = true;
        }
    }
    if !touched.is_empty() {
        FilledClass::AmbiguousOnBoundary { touched }
    } else if any_filled {
        FilledClass::Filled
    } else {
        FilledClass::Empty
    }
}

/// Compute the boolean-union outline of `inputs`. See module docs for
/// the algorithm. Output boundary loops are closed implicitly (vertex
/// list `[v0, v1, ..., vn-1]` represents the closed loop
/// `v0 → v1 → ... → vn-1 → v0`).
///
/// # Errors
/// `OperationError::Failed` if a half-edge's bilateral classification
/// remains `AmbiguousOnBoundary` after 3 ε-shrink retries (typically
/// indicates degenerate input where multiple inputs share a tangent
/// boundary at the sampled edge midpoint).
pub fn union_all_with_holes(inputs: &[PolygonWithHoles]) -> Result<UnionResult> {
    let raw_segments = collect_raw_segments(inputs);
    if raw_segments.is_empty() {
        return Ok(UnionResult {
            boundaries: Vec::new(),
        });
    }
    let split_segments = arrangement_split(&raw_segments);
    let (vertex_table, snapped_segments) = vertex_snap(&split_segments);
    let undirected = canonicalize_undirected(snapped_segments);
    let kept_half_edges = classify_and_filter(&undirected, &vertex_table, inputs)?;
    let mut boundaries = face_walk(&kept_half_edges, &vertex_table);
    boundaries.retain(|p| signed_area(p).abs() > WALL_EPS_SQ);
    debug_assert_cdt_safe(&boundaries);
    Ok(UnionResult { boundaries })
}

/// Defense-in-depth post-condition: verifies that the output boundary set
/// is safe for ingestion by `spade::ConstrainedDelaunayTriangulation`'s
/// `try_add_constraint`. Active only in debug builds; release builds
/// skip the check entirely. Panics on failure rather than returning Err
/// so the bug is surfaced loudly during development.
#[cfg(debug_assertions)]
fn debug_assert_cdt_safe(boundaries: &[Polygon]) {
    use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
    let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> =
        ConstrainedDelaunayTriangulation::new();
    for (bi, boundary) in boundaries.iter().enumerate() {
        let n = boundary.len();
        if n < 3 {
            continue;
        }
        let mut handles = Vec::with_capacity(n);
        for (vi, &(x, y)) in boundary.iter().enumerate() {
            match cdt.insert(Point2::new(x, y)) {
                Ok(h) => handles.push(h),
                Err(e) => panic!(
                    "polygon_union post-condition: CDT vertex insert rejected \
                     (b={bi}, v={vi}, p=({x:.6},{y:.6})): {e:?}"
                ),
            }
        }
        for k in 0..n {
            let from = handles[k];
            let to = handles[(k + 1) % n];
            if from == to {
                continue;
            }
            assert!(
                !cdt.try_add_constraint(from, to).is_empty(),
                "polygon_union post-condition: CDT constraint rejected \
                 (b={bi}, edge={k}): would intersect an existing constraint",
            );
        }
    }
}

#[cfg(not(debug_assertions))]
#[inline]
fn debug_assert_cdt_safe(_boundaries: &[Polygon]) {}

// === Internals ===

type RawSegment = ((f64, f64), (f64, f64));

fn collect_raw_segments(inputs: &[PolygonWithHoles]) -> Vec<RawSegment> {
    let mut out = Vec::new();
    for pwh in inputs {
        for ring in std::iter::once(&pwh.outer).chain(pwh.holes.iter()) {
            let n = ring.len();
            if n < 3 {
                continue;
            }
            for i in 0..n {
                out.push((ring[i], ring[(i + 1) % n]));
            }
        }
    }
    out
}

/// Split every segment at every transverse crossing, T-junction
/// (endpoint-on-edge), and collinear-overlap endpoint with every other
/// segment. Returns the resulting list of sub-segments.
fn arrangement_split(segs: &[RawSegment]) -> Vec<RawSegment> {
    let mut out = Vec::new();
    for (si, &(a0, a1)) in segs.iter().enumerate() {
        let mut params: Vec<f64> = vec![0.0, 1.0];
        for (sj, &(b0, b1)) in segs.iter().enumerate() {
            if si == sj {
                continue;
            }
            if let Some((t, _u)) = seg_seg_intersect(a0, a1, b0, b1) {
                if t > WALL_EPS && t < 1.0 - WALL_EPS {
                    params.push(t);
                }
            } else {
                for t in collinear_overlap_params(a0, a1, b0, b1) {
                    params.push(t);
                }
            }
            for ep in [b0, b1] {
                if let Some(t) = project_endpoint_on_interior(a0, a1, ep) {
                    params.push(t);
                }
            }
        }
        params.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        params.dedup_by(|x, y| (*x - *y).abs() < WALL_EPS);
        for w in params.windows(2) {
            if (w[1] - w[0]).abs() < WALL_EPS {
                continue;
            }
            let p0 = lerp(a0, a1, w[0]);
            let p1 = lerp(a0, a1, w[1]);
            out.push((p0, p1));
        }
    }
    out
}

/// Returns `Some(t)` if `p` lies on the interior of segment `(a0, a1)`
/// within perpendicular distance `< WALL_EPS` and parameter `t` in
/// `(WALL_EPS, 1 - WALL_EPS)`. Returns `None` otherwise.
fn project_endpoint_on_interior(a0: (f64, f64), a1: (f64, f64), p: (f64, f64)) -> Option<f64> {
    let dx = a1.0 - a0.0;
    let dy = a1.1 - a0.1;
    let len_sq = dx * dx + dy * dy;
    if len_sq < WALL_EPS_SQ {
        return None;
    }
    let t = ((p.0 - a0.0) * dx + (p.1 - a0.1) * dy) / len_sq;
    if t <= WALL_EPS || t >= 1.0 - WALL_EPS {
        return None;
    }
    let proj = (a0.0 + t * dx, a0.1 + t * dy);
    let perp_sq = (p.0 - proj.0).powi(2) + (p.1 - proj.1).powi(2);
    if perp_sq < WALL_EPS_SQ {
        Some(t)
    } else {
        None
    }
}

/// `(vertex_table, snapped_segments)` returned by [`vertex_snap`].
type SnapResult = (Vec<(f64, f64)>, Vec<(usize, usize)>);

/// Grid-quantized vertex snap with bounded cross-cell merging.
///
/// Output `vertex_table[class_id]` is the representative point (sorted
/// arithmetic mean of cluster members). `snapped_segments[i]` is the
/// `(start_class_id, end_class_id)` for sub-edge `i`; zero-length
/// sub-edges are dropped. Cluster diameter is bounded by `2·WALL_EPS`.
#[allow(
    clippy::too_many_lines,
    reason = "single coherent pass: build vertices → grid quantize → bounded merge → \
              representative-mean → class id assignment → re-snap. Splitting just for \
              line count would obscure data flow."
)]
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    reason = "f64→i64 cast for grid cell key is intentional and bounded by typical \
              CAD coordinate range (no overflow); usize→f64 for cluster size mean \
              loses no relevant precision (cluster sizes are O(10) at most)."
)]
fn vertex_snap(segs: &[RawSegment]) -> SnapResult {
    if segs.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let mut vertices: Vec<(f64, f64)> = Vec::with_capacity(segs.len() * 2);
    let mut seg_endpoints: Vec<(usize, usize)> = Vec::with_capacity(segs.len());
    for &(s, e) in segs {
        let si = vertices.len();
        vertices.push(s);
        let ei = vertices.len();
        vertices.push(e);
        seg_endpoints.push((si, ei));
    }
    let n = vertices.len();
    let cells: Vec<(i64, i64)> = vertices
        .iter()
        .map(|&v| {
            (
                (v.0 / WALL_EPS).round() as i64,
                (v.1 / WALL_EPS).round() as i64,
            )
        })
        .collect();
    let mut parent: Vec<usize> = (0..n).collect();

    let mut by_cell: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (i, &c) in cells.iter().enumerate() {
        by_cell.entry(c).or_default().push(i);
    }

    // (a) Unconditional intra-cell merge.
    for members in by_cell.values() {
        if members.len() <= 1 {
            continue;
        }
        let anchor = members[0];
        for &m in &members[1..] {
            let ra = find_root(&mut parent, anchor);
            let rb = find_root(&mut parent, m);
            if ra != rb {
                parent[rb] = ra;
            }
        }
    }

    // (b) Bounded cross-cell merge. 4 forward neighbors per cell to avoid
    // double-checking; the remaining 4 of the 8-neighborhood are covered
    // by symmetry.
    let neighbor_offsets: [(i64, i64); 4] = [(0, 1), (1, -1), (1, 0), (1, 1)];
    let cell_keys: Vec<(i64, i64)> = by_cell.keys().copied().collect();
    for &(cx, cy) in &cell_keys {
        for &(dx, dy) in &neighbor_offsets {
            let neighbor_key = (cx + dx, cy + dy);
            let Some(my_members) = by_cell.get(&(cx, cy)).cloned() else {
                continue;
            };
            let Some(other_members) = by_cell.get(&neighbor_key).cloned() else {
                continue;
            };
            for &i in &my_members {
                for &j in &other_members {
                    let ddx = vertices[i].0 - vertices[j].0;
                    let ddy = vertices[i].1 - vertices[j].1;
                    if ddx * ddx + ddy * ddy >= WALL_EPS_SQ {
                        continue;
                    }
                    let ri = find_root(&mut parent, i);
                    let rj = find_root(&mut parent, j);
                    if ri == rj {
                        continue;
                    }
                    let combined: Vec<usize> = (0..n)
                        .filter(|&k| {
                            let r = find_root(&mut parent, k);
                            r == ri || r == rj
                        })
                        .collect();
                    let mut max_d_sq = 0.0_f64;
                    for ai in 0..combined.len() {
                        for bi in (ai + 1)..combined.len() {
                            let ix = combined[ai];
                            let jx = combined[bi];
                            let dxs = vertices[ix].0 - vertices[jx].0;
                            let dys = vertices[ix].1 - vertices[jx].1;
                            let d_sq = dxs * dxs + dys * dys;
                            if d_sq > max_d_sq {
                                max_d_sq = d_sq;
                            }
                        }
                    }
                    let limit_sq = (2.0 * WALL_EPS) * (2.0 * WALL_EPS);
                    if max_d_sq <= limit_sq {
                        parent[rj] = ri;
                    }
                }
            }
        }
    }

    // Collect cluster members.
    let mut clusters: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let r = find_root(&mut parent, i);
        clusters.entry(r).or_default().push(i);
    }

    // Compute representative points (sorted-mean for determinism).
    let mut cluster_reps: Vec<(usize, (f64, f64))> = Vec::with_capacity(clusters.len());
    for (root, members) in clusters {
        let mut sorted_members = members;
        sorted_members.sort_by(|&a, &b| {
            vertices[a]
                .partial_cmp(&vertices[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let count = sorted_members.len() as f64;
        let mut sx = 0.0_f64;
        let mut sy = 0.0_f64;
        for &m in &sorted_members {
            sx += vertices[m].0;
            sy += vertices[m].1;
        }
        cluster_reps.push((root, (sx / count, sy / count)));
    }
    cluster_reps.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut class_of_root: HashMap<usize, usize> = HashMap::new();
    let mut points: Vec<(f64, f64)> = Vec::with_capacity(cluster_reps.len());
    for (cls_id, (root, rep)) in cluster_reps.into_iter().enumerate() {
        class_of_root.insert(root, cls_id);
        points.push(rep);
    }

    let mut snapped: Vec<(usize, usize)> = Vec::with_capacity(seg_endpoints.len());
    for (si, ei) in seg_endpoints {
        let cs = class_of_root[&find_root(&mut parent, si)];
        let ce = class_of_root[&find_root(&mut parent, ei)];
        if cs != ce {
            snapped.push((cs, ce));
        }
    }
    (points, snapped)
}

fn find_root(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

/// Dedup directed sub-edges into a deterministic undirected list. Each
/// sub-edge becomes a sorted endpoint-index tuple `(min, max)`; the
/// returned vector is sorted lexicographically. Combined with
/// `vertex_snap`'s order-independent class id assignment, this guarantees
/// that the output of `union_all_with_holes` is topologically identical
/// regardless of input order.
fn canonicalize_undirected(snapped: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    let mut canon: Vec<(usize, usize)> = snapped
        .into_iter()
        .map(|(a, b)| if a <= b { (a, b) } else { (b, a) })
        .collect();
    canon.sort_unstable();
    canon.dedup();
    canon
}

/// Generate two directed half-edges per undirected sub-edge, classify
/// each by bilateral perpendicular-ε sampling, and return only those
/// half-edges with `Filled` on the LEFT and `Empty` on the RIGHT.
fn classify_and_filter(
    undirected: &[(usize, usize)],
    vertex_table: &[(f64, f64)],
    inputs: &[PolygonWithHoles],
) -> Result<Vec<(usize, usize)>> {
    let mut out: Vec<(usize, usize)> = Vec::with_capacity(undirected.len() * 2);
    for &(u, v) in undirected {
        for (a, b) in [(u, v), (v, u)] {
            let pa = vertex_table[a];
            let pb = vertex_table[b];
            let dx = pb.0 - pa.0;
            let dy = pb.1 - pa.1;
            let len = (dx * dx + dy * dy).sqrt();
            if len < WALL_EPS {
                continue;
            }
            let mid = ((pa.0 + pb.0) * 0.5, (pa.1 + pb.1) * 0.5);
            let nx = -dy / len;
            let ny = dx / len;

            let initial_eps = (WALL_EPS * 10.0).min(len * 0.1);
            let mut eps = initial_eps;
            let mut tries = 0;
            let (left, right) = loop {
                let l = classify_filled((mid.0 + eps * nx, mid.1 + eps * ny), inputs);
                let r = classify_filled((mid.0 - eps * nx, mid.1 - eps * ny), inputs);
                let l_amb = matches!(l, FilledClass::AmbiguousOnBoundary { .. });
                let r_amb = matches!(r, FilledClass::AmbiguousOnBoundary { .. });
                if !l_amb && !r_amb {
                    break (l, r);
                }
                tries += 1;
                if tries >= 3 {
                    return Err(OperationError::Failed(format!(
                        "polygon_union: ambiguous half-edge classification at \
                         edge ({pa:?} → {pb:?}) (mid=({:.6}, {:.6})); ε exhausted; \
                         left={l:?} right={r:?}",
                        mid.0, mid.1
                    ))
                    .into());
                }
                eps *= 0.5;
            };
            if matches!(left, FilledClass::Filled) && matches!(right, FilledClass::Empty) {
                out.push((a, b));
            }
        }
    }
    Ok(out)
}

/// Trace closed boundary loops by walking kept half-edges.
///
/// At each vertex, the successor is picked by the polar-angle Δ rule:
/// for incoming half-edge `u→v`, the outgoing successor is the kept
/// half-edge `v→w` minimizing the clockwise Δ from the reverse
/// direction `θ_in = atan2(u.y-v.y, u.x-v.x)`. The self-reverse
/// half-edge (Δ ≈ 0) is treated as 2π so it is picked last; in practice
/// the self-reverse rarely exists since each undirected sub-edge
/// contributes only one kept direction.
fn face_walk(kept: &[(usize, usize)], vertex_table: &[(f64, f64)]) -> Vec<Polygon> {
    if kept.is_empty() {
        return Vec::new();
    }
    let n_classes = vertex_table.len();
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); n_classes];
    for (idx, &(a, _)) in kept.iter().enumerate() {
        adjacency[a].push(idx);
    }
    let mut used: Vec<bool> = vec![false; kept.len()];
    let mut boundaries: Vec<Polygon> = Vec::new();

    for start in 0..kept.len() {
        if used[start] {
            continue;
        }
        let mut path: Vec<usize> = Vec::new();
        let mut current = start;
        let max_steps = kept.len() + 1;
        for _ in 0..max_steps {
            if used[current] {
                break;
            }
            used[current] = true;
            let (a, b) = kept[current];
            path.push(a);
            let pa = vertex_table[a];
            let pb = vertex_table[b];
            let theta_in = (pa.1 - pb.1).atan2(pa.0 - pb.0);
            let two_pi = 2.0 * std::f64::consts::PI;
            let mut best: Option<(usize, f64)> = None;
            for &idx2 in &adjacency[b] {
                if used[idx2] {
                    continue;
                }
                let target = kept[idx2].1;
                let pt = vertex_table[target];
                let theta_k = (pt.1 - pb.1).atan2(pt.0 - pb.0);
                let mut delta = (theta_in - theta_k).rem_euclid(two_pi);
                if delta < WALL_EPS {
                    delta = two_pi;
                }
                if best.is_none_or(|(_, d)| delta < d) {
                    best = Some((idx2, delta));
                }
            }
            match best {
                Some((next, _)) => {
                    if next == start {
                        used[next] = true;
                        break;
                    }
                    current = next;
                }
                None => break,
            }
        }
        if path.len() >= 3 {
            let polygon: Polygon = path.iter().map(|&i| vertex_table[i]).collect();
            boundaries.push(polygon);
        }
    }
    boundaries
}

// === Geometry primitives ===

/// When two segments are parallel and collinear, return the parameter
/// values on segment A where segment B's endpoints project (clamped to
/// the interior of A). Returns an empty vec when the segments are not
/// collinear or have no interior overlap.
fn collinear_overlap_params(
    a0: (f64, f64),
    a1: (f64, f64),
    b0: (f64, f64),
    b1: (f64, f64),
) -> Vec<f64> {
    let d1x = a1.0 - a0.0;
    let d1y = a1.1 - a0.1;
    let d2x = b1.0 - b0.0;
    let d2y = b1.1 - b0.1;
    let dir_cross = d1x * d2y - d1y * d2x;
    if dir_cross.abs() >= WALL_EPS {
        return Vec::new();
    }
    let bax = b0.0 - a0.0;
    let bay = b0.1 - a0.1;
    let colinear_cross = d1x * bay - d1y * bax;
    if colinear_cross.abs() >= WALL_EPS {
        return Vec::new();
    }
    let len_sq = d1x * d1x + d1y * d1y;
    if len_sq < WALL_EPS_SQ {
        return Vec::new();
    }
    let t0 = (bax * d1x + bay * d1y) / len_sq;
    let t1 = ((b1.0 - a0.0) * d1x + (b1.1 - a0.1) * d1y) / len_sq;
    let mut out = Vec::new();
    if t0 > WALL_EPS && t0 < 1.0 - WALL_EPS {
        out.push(t0);
    }
    if t1 > WALL_EPS && t1 < 1.0 - WALL_EPS {
        out.push(t1);
    }
    out
}

fn seg_seg_intersect(
    a0: (f64, f64),
    a1: (f64, f64),
    b0: (f64, f64),
    b1: (f64, f64),
) -> Option<(f64, f64)> {
    let d1x = a1.0 - a0.0;
    let d1y = a1.1 - a0.1;
    let d2x = b1.0 - b0.0;
    let d2y = b1.1 - b0.1;
    let cross = d1x * d2y - d1y * d2x;
    if cross.abs() < WALL_EPS {
        return None;
    }
    let d3x = b0.0 - a0.0;
    let d3y = b0.1 - a0.1;
    let t = (d3x * d2y - d3y * d2x) / cross;
    let u = (d3x * d1y - d3y * d1x) / cross;
    if t > -WALL_EPS && t < 1.0 + WALL_EPS && u > -WALL_EPS && u < 1.0 + WALL_EPS {
        Some((t.clamp(0.0, 1.0), u.clamp(0.0, 1.0)))
    } else {
        None
    }
}

pub(crate) fn point_in_polygon_class(p: (f64, f64), poly: &Polygon) -> PointClass {
    let n = poly.len();
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        let dist = point_to_segment_dist(p.0, p.1, a.0, a.1, b.0, b.1);
        if dist < WALL_EPS {
            return PointClass::Boundary;
        }
    }
    let mut winding = 0i32;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        if a.1 <= p.1 {
            if b.1 > p.1 && cross_2d(a, b, p) > 0.0 {
                winding += 1;
            }
        } else if b.1 <= p.1 && cross_2d(a, b, p) < 0.0 {
            winding -= 1;
        }
    }
    if winding != 0 {
        PointClass::Inside
    } else {
        PointClass::Outside
    }
}

fn cross_2d(a: (f64, f64), b: (f64, f64), p: (f64, f64)) -> f64 {
    (b.0 - a.0) * (p.1 - a.1) - (b.1 - a.1) * (p.0 - a.0)
}

pub fn signed_area(poly: &Polygon) -> f64 {
    let n = poly.len();
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += poly[i].0 * poly[j].1;
        area -= poly[j].0 * poly[i].1;
    }
    area * 0.5
}

fn lerp(a: (f64, f64), b: (f64, f64), t: f64) -> (f64, f64) {
    (a.0 + t * (b.0 - a.0), a.1 + t * (b.1 - a.1))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;

    fn rect(x: f64, y: f64, w: f64, h: f64) -> Polygon {
        vec![(x, y), (x + w, y), (x + w, y + h), (x, y + h)]
    }

    fn no_hole_inputs(polys: Vec<Polygon>) -> Vec<PolygonWithHoles> {
        polys
            .into_iter()
            .map(|outer| PolygonWithHoles {
                outer,
                holes: Vec::new(),
            })
            .collect()
    }

    fn segment_to_rect(a: (f64, f64), b: (f64, f64), lw: f64, rw: f64) -> Polygon {
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let len = (dx * dx + dy * dy).sqrt();
        let (nx, ny) = (-dy / len, dx / len);
        vec![
            (a.0 + lw * nx, a.1 + lw * ny),
            (b.0 + lw * nx, b.1 + lw * ny),
            (b.0 - rw * nx, b.1 - rw * ny),
            (a.0 - rw * nx, a.1 - rw * ny),
        ]
    }

    #[test]
    fn union_non_overlapping() {
        let result = union_all_with_holes(&no_hole_inputs(vec![
            rect(0.0, 0.0, 2.0, 2.0),
            rect(5.0, 0.0, 2.0, 2.0),
        ]))
        .expect("union must succeed");
        assert_eq!(result.boundaries.len(), 2);
    }

    #[test]
    fn union_overlapping() {
        let result = union_all_with_holes(&no_hole_inputs(vec![
            rect(0.0, 0.0, 3.0, 2.0),
            rect(2.0, 0.0, 3.0, 2.0),
        ]))
        .expect("union must succeed");
        assert_eq!(result.boundaries.len(), 1);
        let area = signed_area(&result.boundaries[0]);
        assert!((area - 10.0).abs() < 0.1, "area={area}");
    }

    #[test]
    fn union_shared_edge() {
        let result = union_all_with_holes(&no_hole_inputs(vec![
            rect(0.0, 0.0, 4.0, 3.0),
            rect(4.0, 0.0, 4.0, 3.0),
        ]))
        .expect("union must succeed");
        assert_eq!(result.boundaries.len(), 1);
        let area = signed_area(&result.boundaries[0]);
        assert!((area - 24.0).abs() < 0.1, "area={area}");
    }

    #[test]
    fn union_contained() {
        let result = union_all_with_holes(&no_hole_inputs(vec![
            rect(0.0, 0.0, 6.0, 6.0),
            rect(1.0, 1.0, 2.0, 2.0),
        ]))
        .expect("union must succeed");
        assert_eq!(result.boundaries.len(), 1);
        let area = signed_area(&result.boundaries[0]);
        assert!((area - 36.0).abs() < 0.1, "area={area}");
    }

    #[test]
    fn union_t_shape() {
        let result = union_all_with_holes(&no_hole_inputs(vec![
            rect(0.0, -1.0, 8.0, 2.0),
            rect(3.0, -1.0, 2.0, 5.0),
        ]))
        .expect("union must succeed");
        assert_eq!(result.boundaries.len(), 1);
        let expected_area = 8.0 * 2.0 + 2.0 * 5.0 - 2.0 * 2.0;
        let area = signed_area(&result.boundaries[0]);
        assert!((area - expected_area).abs() < 0.1, "area={area}");
    }

    #[test]
    fn union_cross_shape() {
        let result = union_all_with_holes(&no_hole_inputs(vec![
            rect(0.0, 1.0, 6.0, 2.0),
            rect(2.0, 0.0, 2.0, 4.0),
        ]))
        .expect("union must succeed");
        assert_eq!(result.boundaries.len(), 1);
        let expected_area = 6.0 * 2.0 + 2.0 * 4.0 - 2.0 * 2.0;
        let area = signed_area(&result.boundaries[0]);
        assert!((area - expected_area).abs() < 0.1, "area={area}");
    }

    #[test]
    fn union_donut_from_four_rects() {
        // Four rectangles forming a closed square wall → outer (CCW) + hole (CW).
        let d = 0.3;
        let inputs = no_hole_inputs(vec![
            segment_to_rect((0.0, 0.0), (10.0, 0.0), d, d),
            segment_to_rect((10.0, 0.0), (10.0, 10.0), d, d),
            segment_to_rect((10.0, 10.0), (0.0, 10.0), d, d),
            segment_to_rect((0.0, 10.0), (0.0, 0.0), d, d),
        ]);
        let result = union_all_with_holes(&inputs).expect("union must succeed");
        assert_eq!(result.boundaries.len(), 2, "expected outer + hole");
        let areas: Vec<f64> = result.boundaries.iter().map(signed_area).collect();
        assert!(areas.iter().any(|a| *a > 0.0), "needs CCW outer");
        assert!(areas.iter().any(|a| *a < 0.0), "needs CW hole");
    }

    #[test]
    fn union_wall_segments_t_junction() {
        let d = 0.15;
        let result = union_all_with_holes(&no_hole_inputs(vec![
            segment_to_rect((0.0, 0.0), (4.0, 0.0), d, d),
            segment_to_rect((4.0, 0.0), (4.0, 3.0), d, d),
            segment_to_rect((4.0, 0.0), (8.0, 0.0), d, d),
        ]))
        .expect("union must succeed");
        assert!(!result.boundaries.is_empty());
        for b in &result.boundaries {
            for &(x, y) in b {
                assert!((-0.5..=8.5).contains(&x), "x={x} out of range");
                assert!((-0.5..=3.5).contains(&y), "y={y} out of range");
            }
        }
    }

    #[test]
    fn union_angled_wall_segments() {
        let d = 0.15;
        let inputs = no_hole_inputs(vec![
            segment_to_rect((-3.217, -4.144), (-2.635, 2.085), d, d),
            segment_to_rect((-3.217, -4.144), (2.002, -4.631), d, d),
            segment_to_rect((-2.635, 2.085), (2.578, 1.534), d, d),
            segment_to_rect((2.002, -4.631), (2.578, 1.534), d, d),
            segment_to_rect((2.002, -4.631), (6.473, -5.049), d, d),
            segment_to_rect((2.578, 1.534), (6.861, -0.896), d, d),
            segment_to_rect((6.473, -5.049), (6.861, -0.896), d, d),
        ]);
        let result = union_all_with_holes(&inputs).expect("union must succeed");
        assert!(!result.boundaries.is_empty());
        for b in &result.boundaries {
            for &(x, y) in b {
                assert!(
                    (-4.0..=8.0).contains(&x) && (-6.0..=3.0).contains(&y),
                    "vertex ({x:.3}, {y:.3}) out of expected range"
                );
            }
        }
    }

    // --- Tests for union_all_with_holes (production path) ---

    #[test]
    fn with_holes_single_ring() {
        let pwh = PolygonWithHoles {
            outer: rect(0.0, 0.0, 5.0, 3.0),
            holes: Vec::new(),
        };
        let result = union_all_with_holes(&[pwh]).expect("union must succeed");
        assert_eq!(result.boundaries.len(), 1);
        let area = signed_area(&result.boundaries[0]);
        assert!(area > 0.0, "outer should be CCW, area={area}");
    }

    #[test]
    fn with_holes_donut() {
        let outer = rect(0.0, 0.0, 10.0, 10.0);
        let hole = vec![(2.0, 2.0), (2.0, 8.0), (8.0, 8.0), (8.0, 2.0)]; // CW
        let pwh = PolygonWithHoles {
            outer,
            holes: vec![hole],
        };
        let result = union_all_with_holes(&[pwh]).expect("union must succeed");
        assert_eq!(result.boundaries.len(), 2, "outer + hole");
        let areas: Vec<f64> = result.boundaries.iter().map(signed_area).collect();
        assert!(areas.iter().any(|a| *a > 0.0), "needs CCW outer");
        assert!(areas.iter().any(|a| *a < 0.0), "needs CW hole");
    }

    #[test]
    fn with_holes_two_rings_union() {
        let a = PolygonWithHoles {
            outer: rect(0.0, 0.0, 4.0, 3.0),
            holes: Vec::new(),
        };
        let b = PolygonWithHoles {
            outer: rect(2.0, 0.0, 4.0, 3.0),
            holes: Vec::new(),
        };
        let result = union_all_with_holes(&[a, b]).expect("union must succeed");
        assert_eq!(result.boundaries.len(), 1);
        let area = signed_area(&result.boundaries[0]);
        assert!((area.abs() - 18.0).abs() < 1.0, "area={area}");
    }

    #[test]
    fn with_holes_two_donuts_union() {
        let a = PolygonWithHoles {
            outer: rect(0.0, 0.0, 6.0, 6.0),
            holes: vec![rect(1.0, 1.0, 4.0, 4.0)],
        };
        let b = PolygonWithHoles {
            outer: rect(3.0, 0.0, 6.0, 6.0),
            holes: vec![rect(4.0, 1.0, 4.0, 4.0)],
        };
        let result = union_all_with_holes(&[a, b]).expect("union must succeed");
        assert!(!result.boundaries.is_empty());
        let ccw_count = result
            .boundaries
            .iter()
            .filter(|b| signed_area(b) > 0.0)
            .count();
        assert!(ccw_count >= 1, "needs at least one outer");
    }

    #[test]
    fn with_holes_rings_sharing_a_colinear_face() {
        let a_outer = rect(0.0, 0.0, 5.0, 3.0);
        let a_hole = vec![(0.5, 0.5), (0.5, 2.5), (4.5, 2.5), (4.5, 0.5)];
        let b_outer = rect(4.0, 0.0, 5.0, 3.0);
        let b_hole = vec![(4.5, 0.5), (4.5, 2.5), (8.5, 2.5), (8.5, 0.5)];
        let a = PolygonWithHoles {
            outer: a_outer,
            holes: vec![a_hole],
        };
        let b = PolygonWithHoles {
            outer: b_outer,
            holes: vec![b_hole],
        };
        let result = union_all_with_holes(&[a, b]).expect("union must succeed");
        let outers: Vec<&Polygon> = result
            .boundaries
            .iter()
            .filter(|b| signed_area(b) > 0.0)
            .collect();
        assert_eq!(outers.len(), 1, "one combined outer");
        for &(x, y) in outers[0] {
            let on_south = (y - 0.0).abs() < WALL_EPS;
            let on_north = (y - 3.0).abs() < WALL_EPS;
            let on_west = (x - 0.0).abs() < WALL_EPS;
            let on_east = (x - 9.0).abs() < WALL_EPS;
            assert!(
                on_south || on_north || on_west || on_east,
                "vertex ({x:.3}, {y:.3}) lies off the combined rectangle boundary"
            );
        }
    }

    #[test]
    fn with_holes_two_open_wall_strokes_t_junction() {
        let horiz = PolygonWithHoles {
            outer: rect(0.0, -0.15, 4.0, 0.30),
            holes: Vec::new(),
        };
        let vert = PolygonWithHoles {
            outer: rect(1.85, 0.0, 0.30, 3.0),
            holes: Vec::new(),
        };
        let result = union_all_with_holes(&[horiz, vert]).expect("union must succeed");
        assert_eq!(
            result.boundaries.len(),
            1,
            "two overlapping wall strokes must union into one T boundary",
        );
    }

    #[test]
    fn with_holes_two_adjacent_zones_produce_one_outer_two_holes() {
        let d = 0.15;
        let a = PolygonWithHoles {
            outer: vec![(-d, -d), (5.0 + d, -d), (5.0 + d, 3.0 + d), (-d, 3.0 + d)],
            holes: vec![vec![(d, d), (d, 3.0 - d), (5.0 - d, 3.0 - d), (5.0 - d, d)]],
        };
        let b = PolygonWithHoles {
            outer: vec![
                (5.0 - d, -d),
                (8.0 + d, -d),
                (8.0 + d, 3.0 + d),
                (5.0 - d, 3.0 + d),
            ],
            holes: vec![vec![
                (5.0 + d, d),
                (5.0 + d, 3.0 - d),
                (8.0 - d, 3.0 - d),
                (8.0 - d, d),
            ]],
        };
        let result = union_all_with_holes(&[a, b]).expect("union must succeed");
        let outers: Vec<&Polygon> = result
            .boundaries
            .iter()
            .filter(|b| signed_area(b) > 0.0)
            .collect();
        let holes: Vec<&Polygon> = result
            .boundaries
            .iter()
            .filter(|b| signed_area(b) < 0.0)
            .collect();
        assert_eq!(outers.len(), 1);
        assert_eq!(holes.len(), 2);
    }

    #[test]
    fn with_holes_non_overlapping() {
        let a = PolygonWithHoles {
            outer: rect(0.0, 0.0, 2.0, 2.0),
            holes: Vec::new(),
        };
        let b = PolygonWithHoles {
            outer: rect(5.0, 0.0, 2.0, 2.0),
            holes: Vec::new(),
        };
        let result = union_all_with_holes(&[a, b]).expect("union must succeed");
        assert_eq!(result.boundaries.len(), 2, "two separate outers");
        let ccw_count = result
            .boundaries
            .iter()
            .filter(|b| signed_area(b) > 0.0)
            .count();
        assert_eq!(ccw_count, 2, "both should be CCW outers");
    }

    // --- P3.2 algorithm verification fixtures ---

    /// Concentric outer + hole as a single PWH input. The polar-angle Δ
    /// face-walk rule must walk the outer in CCW order (signed area > 0)
    /// and the hole in CW order (signed area < 0).
    #[test]
    fn arrangement_concentric_square_outer_ccw_hole_cw() {
        let pwh = PolygonWithHoles {
            outer: rect(0.0, 0.0, 10.0, 10.0),
            holes: vec![vec![(2.0, 2.0), (2.0, 8.0), (8.0, 8.0), (8.0, 2.0)]],
        };
        let r = union_all_with_holes(&[pwh]).expect("union must succeed");
        assert_eq!(r.boundaries.len(), 2, "outer + hole");
        let outer_count = r.boundaries.iter().filter(|b| signed_area(b) > 0.0).count();
        let hole_count = r.boundaries.iter().filter(|b| signed_area(b) < 0.0).count();
        assert_eq!(outer_count, 1, "outer must be CCW (signed area > 0)");
        assert_eq!(hole_count, 1, "hole must be CW (signed area < 0)");
    }

    /// Two squares touching at a single vertex (degree-4 in the
    /// arrangement). The face-walk must NOT emit the same undirected
    /// edge in two different output loops.
    #[test]
    #[allow(
        clippy::cast_possible_truncation,
        reason = "input coordinates are bounded ints (0..8); quantizing by 1/WALL_EPS \
                  yields values well within i64 range"
    )]
    fn arrangement_degree_4_two_squares_share_one_vertex() {
        use std::collections::HashSet;
        const Q: f64 = 1.0 / WALL_EPS;

        let a = PolygonWithHoles {
            outer: rect(0.0, 0.0, 4.0, 4.0),
            holes: Vec::new(),
        };
        let b = PolygonWithHoles {
            outer: rect(4.0, 4.0, 4.0, 4.0),
            holes: Vec::new(),
        };
        let r = union_all_with_holes(&[a, b]).expect("union must succeed");
        let mut seen: HashSet<((i64, i64), (i64, i64))> = HashSet::new();
        for boundary in &r.boundaries {
            let n = boundary.len();
            for i in 0..n {
                let p0 = boundary[i];
                let p1 = boundary[(i + 1) % n];
                let qa = ((p0.0 * Q).round() as i64, (p0.1 * Q).round() as i64);
                let qb = ((p1.0 * Q).round() as i64, (p1.1 * Q).round() as i64);
                let key = if qa <= qb { (qa, qb) } else { (qb, qa) };
                assert!(
                    seen.insert(key),
                    "edge {qa:?}-{qb:?} appears in two output loops",
                );
            }
        }
    }
}
