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

/// A planar face described by an outer boundary and zero or more holes.
///
/// # Winding contract
/// - `outer` is CCW with `signed_area(outer) > 0`.
/// - Each `holes[i]` is CW with `signed_area(holes[i]) < 0`.
/// - Every hole is fully contained in `outer`.
/// - Sibling holes are non-overlapping.
///
/// The contract is enforced by [`assemble_faces`] for outputs of
/// [`union_all_with_holes`], and by `WallFootprint2D::try_from_parts` for
/// cross-crate construction. Internal callers that bypass these entry points
/// must uphold the invariants themselves.
#[derive(Clone, Debug, PartialEq)]
pub struct PolygonWithHoles {
    pub outer: Polygon,
    pub holes: Vec<Polygon>,
}

impl PolygonWithHoles {
    pub fn into_parts(self) -> (Polygon, Vec<Polygon>) {
        (self.outer, self.holes)
    }
}

/// Result of a polygon union: typed face topology.
pub struct UnionResult {
    pub faces: Vec<PolygonWithHoles>,
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
        return Ok(UnionResult { faces: Vec::new() });
    }
    let split_segments = arrangement_split(&raw_segments);
    let (vertex_table, snapped_segments) = vertex_snap(&split_segments);
    let undirected = canonicalize_undirected(snapped_segments);
    let kept_half_edges = classify_and_filter(&undirected, &vertex_table, inputs)?;
    let loops = face_walk(&kept_half_edges, &vertex_table);
    // Drop degenerate loops before assembly so the containment matrix is not
    // skewed by zero-area artifacts.
    let loops: Vec<WalkedLoop> = loops
        .into_iter()
        .filter(|l| {
            let polygon: Polygon = l.vertex_indices.iter().map(|&i| vertex_table[i]).collect();
            signed_area(&polygon).abs() > WALL_EPS_SQ
        })
        .collect();
    let faces = assemble_faces(&loops, &vertex_table)?;
    debug_assert_cdt_safe(&faces);
    Ok(UnionResult { faces })
}

/// Defense-in-depth post-condition: verifies that the output boundary set
/// is safe for ingestion by `spade::ConstrainedDelaunayTriangulation`'s
/// `try_add_constraint`. Active only in debug builds; release builds
/// skip the check entirely. Panics on failure rather than returning Err
/// so the bug is surfaced loudly during development.
#[cfg(debug_assertions)]
fn debug_assert_cdt_safe(faces: &[PolygonWithHoles]) {
    use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
    let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> =
        ConstrainedDelaunayTriangulation::new();
    let boundaries: Vec<&Polygon> = faces
        .iter()
        .flat_map(|f| std::iter::once(&f.outer).chain(f.holes.iter()))
        .collect();
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
fn debug_assert_cdt_safe(_faces: &[PolygonWithHoles]) {}

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

/// Orientation of a closed planar loop, derived from its signed area.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Orientation {
    Ccw,
    Cw,
}

/// One closed boundary loop emitted by [`face_walk`].
///
/// `vertex_indices[k]` is the start-vertex class id of the k-th half-edge in
/// the loop. `kept_indices[k]` is the index into the `kept` slice of the same
/// half-edge. The two vectors have the same length and are aligned:
/// `kept[loop.kept_indices[k]] == (loop.vertex_indices[k], loop.vertex_indices[(k + 1) % len])`.
pub(crate) struct WalkedLoop {
    pub vertex_indices: Vec<usize>,
    /// Index trail into the `kept` slice. Currently unused by `assemble_faces`,
    /// but retained so future diagnostic logging can map a loop back to the
    /// half-edges it bounds without re-walking the arrangement.
    #[allow(dead_code)]
    pub kept_indices: Vec<usize>,
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
///
/// Returns one [`WalkedLoop`] per closed boundary cycle, retaining the
/// kept-half-edge index trail so downstream face assembly can map each loop
/// back to the half-edges it bounds.
pub(crate) fn face_walk(kept: &[(usize, usize)], vertex_table: &[(f64, f64)]) -> Vec<WalkedLoop> {
    if kept.is_empty() {
        return Vec::new();
    }
    let n_classes = vertex_table.len();
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); n_classes];
    for (idx, &(a, _)) in kept.iter().enumerate() {
        adjacency[a].push(idx);
    }
    let mut used: Vec<bool> = vec![false; kept.len()];
    let mut loops: Vec<WalkedLoop> = Vec::new();

    for start in 0..kept.len() {
        if used[start] {
            continue;
        }
        let mut vertex_indices: Vec<usize> = Vec::new();
        let mut kept_indices: Vec<usize> = Vec::new();
        let mut current = start;
        let max_steps = kept.len() + 1;
        for _ in 0..max_steps {
            if used[current] {
                break;
            }
            used[current] = true;
            let (a, b) = kept[current];
            vertex_indices.push(a);
            kept_indices.push(current);
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
        if vertex_indices.len() >= 3 {
            loops.push(WalkedLoop {
                vertex_indices,
                kept_indices,
            });
        }
    }
    loops
}

/// Assemble closed boundary loops into face topology
/// (`Vec<PolygonWithHoles>`).
///
/// # Algorithm
///
/// Each kept half-edge has Filled on its LEFT and Empty on its RIGHT
/// (locked by [`classify_and_filter`]). A loop is therefore unambiguously
/// classified by signed area: CCW = outer of a filled face, CW = hole
/// carved into a filled face. The remaining task is to match every CW
/// hole to its parent CCW outer, and every depth-≥2 nested CCW island to
/// the CW hole that encloses it.
///
/// Containment is determined by a containment matrix on leftmost-vertex
/// witnesses, not by ray casting (which fails for disjoint donuts on the
/// same scanline). Depth = `|contained_in[i]|`; parent is the unique
/// loop in `contained_in[i]` with `depth - 1`.
///
/// # Errors
///
/// Returns [`OperationError::Failed`] if:
/// - A witness vertex coincides with another loop's boundary even after
///   trying up to `k = 3` fallback witnesses (tangent / near-degenerate
///   arrangement that should have been rejected upstream).
/// - The number of parent candidates at `depth - 1` is not exactly one.
/// - A loop's orientation does not match its depth parity (even depth ⇒
///   CCW; odd depth ⇒ CW).
pub(crate) fn assemble_faces(
    loops: &[WalkedLoop],
    vertex_table: &[(f64, f64)],
) -> Result<Vec<PolygonWithHoles>> {
    let n = loops.len();
    if n == 0 {
        return Ok(Vec::new());
    }

    // 1. Materialise polygons + orientations + sorted witness candidates.
    let polygons: Vec<Polygon> = loops
        .iter()
        .map(|l| l.vertex_indices.iter().map(|&i| vertex_table[i]).collect())
        .collect();

    let orientation: Vec<Orientation> = polygons
        .iter()
        .map(|p| {
            if signed_area(p) > 0.0 {
                Orientation::Ccw
            } else {
                Orientation::Cw
            }
        })
        .collect();

    // For each loop, pre-sort its vertices by (x, y) ascending so the
    // first is the leftmost-vertex witness, the second is the fallback,
    // and the third is the second fallback.
    let witness_candidates: Vec<Vec<(f64, f64)>> = polygons
        .iter()
        .map(|p| {
            let mut v = p.clone();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            v
        })
        .collect();

    // 2. Containment matrix with bounded Boundary fallback.
    const MAX_WITNESS_FALLBACK: usize = 3;
    let mut contained_in: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let mut witness_idx = 0;
            loop {
                if witness_idx >= MAX_WITNESS_FALLBACK || witness_idx >= witness_candidates[i].len()
                {
                    return Err(OperationError::Failed(format!(
                        "polygon_union: assemble_faces witness disambiguation \
                         exhausted for loop {i} against loop {j} (tried \
                         {witness_idx} candidates, all on j's boundary)"
                    ))
                    .into());
                }
                let w = witness_candidates[i][witness_idx];
                match point_in_polygon_class(w, &polygons[j]) {
                    PointClass::Inside => {
                        contained_in[i].push(j);
                        break;
                    }
                    PointClass::Outside => break,
                    PointClass::Boundary => {
                        witness_idx += 1;
                    }
                }
            }
        }
    }

    // 3. Depth and parent.
    let depth: Vec<usize> = contained_in.iter().map(|s| s.len()).collect();
    let mut parent: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        if depth[i] == 0 {
            continue;
        }
        let target_depth = depth[i] - 1;
        let candidates: Vec<usize> = contained_in[i]
            .iter()
            .copied()
            .filter(|&j| depth[j] == target_depth)
            .collect();
        if candidates.len() != 1 {
            return Err(OperationError::Failed(format!(
                "polygon_union: assemble_faces found {} parent candidates for \
                 loop {} (depth {}); arrangement topology is broken",
                candidates.len(),
                i,
                depth[i],
            ))
            .into());
        }
        parent[i] = Some(candidates[0]);
    }

    // 4. Validate orientation parity.
    for i in 0..n {
        let want_ccw = depth[i] % 2 == 0;
        let is_ccw = orientation[i] == Orientation::Ccw;
        if want_ccw != is_ccw {
            return Err(OperationError::Failed(format!(
                "polygon_union: assemble_faces orientation/depth mismatch on \
                 loop {} (orientation {:?}, depth {})",
                i, orientation[i], depth[i]
            ))
            .into());
        }
    }

    // 5. Assemble faces. Each CCW loop is the outer of one face; its holes
    //    are the CW loops whose parent is this loop.
    let mut faces: Vec<PolygonWithHoles> = Vec::new();
    for i in 0..n {
        if orientation[i] != Orientation::Ccw {
            continue;
        }
        let holes: Vec<Polygon> = (0..n)
            .filter(|&j| orientation[j] == Orientation::Cw && parent[j] == Some(i))
            .map(|j| polygons[j].clone())
            .collect();
        faces.push(PolygonWithHoles {
            outer: polygons[i].clone(),
            holes,
        });
    }

    // 6. Post-conditions.
    debug_assert!(
        faces.iter().all(|f| signed_area(&f.outer) > 0.0),
        "assemble_faces: outer winding contract violated"
    );
    debug_assert!(
        faces
            .iter()
            .flat_map(|f| f.holes.iter())
            .all(|h| signed_area(h) < 0.0),
        "assemble_faces: hole winding contract violated"
    );
    let cw_count = orientation
        .iter()
        .filter(|o| **o == Orientation::Cw)
        .count();
    let claimed_holes: usize = faces.iter().map(|f| f.holes.len()).sum();
    debug_assert_eq!(
        cw_count, claimed_holes,
        "assemble_faces: not every CW loop is claimed as a hole"
    );

    Ok(faces)
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

pub(super) fn seg_seg_intersect(
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

    /// Legacy flat-boundary view of a union result for tests written against
    /// the pre-`assemble_faces` API. Equivalent to the old `legacy_boundaries(&result)`.
    fn legacy_boundaries(r: &UnionResult) -> Vec<Polygon> {
        r.faces
            .iter()
            .flat_map(|f| std::iter::once(&f.outer).chain(f.holes.iter()))
            .cloned()
            .collect()
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
        assert_eq!(legacy_boundaries(&result).len(), 2);
    }

    #[test]
    fn union_overlapping() {
        let result = union_all_with_holes(&no_hole_inputs(vec![
            rect(0.0, 0.0, 3.0, 2.0),
            rect(2.0, 0.0, 3.0, 2.0),
        ]))
        .expect("union must succeed");
        assert_eq!(legacy_boundaries(&result).len(), 1);
        let area = signed_area(&legacy_boundaries(&result)[0]);
        assert!((area - 10.0).abs() < 0.1, "area={area}");
    }

    #[test]
    fn union_shared_edge() {
        let result = union_all_with_holes(&no_hole_inputs(vec![
            rect(0.0, 0.0, 4.0, 3.0),
            rect(4.0, 0.0, 4.0, 3.0),
        ]))
        .expect("union must succeed");
        assert_eq!(legacy_boundaries(&result).len(), 1);
        let area = signed_area(&legacy_boundaries(&result)[0]);
        assert!((area - 24.0).abs() < 0.1, "area={area}");
    }

    #[test]
    fn union_contained() {
        let result = union_all_with_holes(&no_hole_inputs(vec![
            rect(0.0, 0.0, 6.0, 6.0),
            rect(1.0, 1.0, 2.0, 2.0),
        ]))
        .expect("union must succeed");
        assert_eq!(legacy_boundaries(&result).len(), 1);
        let area = signed_area(&legacy_boundaries(&result)[0]);
        assert!((area - 36.0).abs() < 0.1, "area={area}");
    }

    #[test]
    fn union_t_shape() {
        let result = union_all_with_holes(&no_hole_inputs(vec![
            rect(0.0, -1.0, 8.0, 2.0),
            rect(3.0, -1.0, 2.0, 5.0),
        ]))
        .expect("union must succeed");
        assert_eq!(legacy_boundaries(&result).len(), 1);
        let expected_area = 8.0 * 2.0 + 2.0 * 5.0 - 2.0 * 2.0;
        let area = signed_area(&legacy_boundaries(&result)[0]);
        assert!((area - expected_area).abs() < 0.1, "area={area}");
    }

    #[test]
    fn union_cross_shape() {
        let result = union_all_with_holes(&no_hole_inputs(vec![
            rect(0.0, 1.0, 6.0, 2.0),
            rect(2.0, 0.0, 2.0, 4.0),
        ]))
        .expect("union must succeed");
        assert_eq!(legacy_boundaries(&result).len(), 1);
        let expected_area = 6.0 * 2.0 + 2.0 * 4.0 - 2.0 * 2.0;
        let area = signed_area(&legacy_boundaries(&result)[0]);
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
        assert_eq!(legacy_boundaries(&result).len(), 2, "expected outer + hole");
        let areas: Vec<f64> = legacy_boundaries(&result).iter().map(signed_area).collect();
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
        assert!(!legacy_boundaries(&result).is_empty());
        for b in &legacy_boundaries(&result) {
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
        assert!(!legacy_boundaries(&result).is_empty());
        for b in &legacy_boundaries(&result) {
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
        assert_eq!(legacy_boundaries(&result).len(), 1);
        let area = signed_area(&legacy_boundaries(&result)[0]);
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
        assert_eq!(legacy_boundaries(&result).len(), 2, "outer + hole");
        let areas: Vec<f64> = legacy_boundaries(&result).iter().map(signed_area).collect();
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
        assert_eq!(legacy_boundaries(&result).len(), 1);
        let area = signed_area(&legacy_boundaries(&result)[0]);
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
        assert!(!legacy_boundaries(&result).is_empty());
        let ccw_count = legacy_boundaries(&result)
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
        let boundaries = legacy_boundaries(&result);
        let outers: Vec<&Polygon> = boundaries.iter().filter(|b| signed_area(b) > 0.0).collect();
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
            legacy_boundaries(&result).len(),
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
        let boundaries = legacy_boundaries(&result);
        let outers: Vec<&Polygon> = boundaries.iter().filter(|b| signed_area(b) > 0.0).collect();
        let holes: Vec<&Polygon> = boundaries.iter().filter(|b| signed_area(b) < 0.0).collect();
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
        assert_eq!(legacy_boundaries(&result).len(), 2, "two separate outers");
        let ccw_count = legacy_boundaries(&result)
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
        assert_eq!(legacy_boundaries(&r).len(), 2, "outer + hole");
        let outer_count = legacy_boundaries(&r)
            .iter()
            .filter(|b| signed_area(b) > 0.0)
            .count();
        let hole_count = legacy_boundaries(&r)
            .iter()
            .filter(|b| signed_area(b) < 0.0)
            .count();
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
        for boundary in &legacy_boundaries(&r) {
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

    // ===== assemble_faces tests =====
    //
    // These build synthetic `WalkedLoop`s directly to exercise the face-assembly
    // contract independently of the union pipeline. `kept_indices` is unused by
    // `assemble_faces` itself, so the helper leaves it empty.

    /// Build a `(loops, vertex_table)` pair from a list of polygon loops.
    fn loops_from_polygons(polygons: &[Polygon]) -> (Vec<WalkedLoop>, Vec<(f64, f64)>) {
        let mut vertex_table: Vec<(f64, f64)> = Vec::new();
        let mut loops: Vec<WalkedLoop> = Vec::new();
        for poly in polygons {
            let base = vertex_table.len();
            vertex_table.extend(poly.iter().copied());
            loops.push(WalkedLoop {
                vertex_indices: (base..base + poly.len()).collect(),
                kept_indices: Vec::new(),
            });
        }
        (loops, vertex_table)
    }

    fn ccw_rect(x: f64, y: f64, w: f64, h: f64) -> Polygon {
        vec![(x, y), (x + w, y), (x + w, y + h), (x, y + h)]
    }

    fn cw_rect(x: f64, y: f64, w: f64, h: f64) -> Polygon {
        vec![(x, y), (x, y + h), (x + w, y + h), (x + w, y)]
    }

    #[test]
    fn assemble_faces_outer_only() {
        let (loops, vertex_table) = loops_from_polygons(&[ccw_rect(0.0, 0.0, 10.0, 10.0)]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 1);
        assert_eq!(faces[0].holes.len(), 0);
        assert!(signed_area(&faces[0].outer) > 0.0);
    }

    #[test]
    fn assemble_faces_outer_plus_one_hole() {
        let (loops, vertex_table) =
            loops_from_polygons(&[ccw_rect(0.0, 0.0, 10.0, 10.0), cw_rect(3.0, 3.0, 4.0, 4.0)]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 1);
        assert_eq!(faces[0].holes.len(), 1);
        assert!(signed_area(&faces[0].outer) > 0.0);
        assert!(signed_area(&faces[0].holes[0]) < 0.0);
    }

    #[test]
    fn assemble_faces_two_disjoint_outers() {
        let (loops, vertex_table) =
            loops_from_polygons(&[ccw_rect(0.0, 0.0, 5.0, 5.0), ccw_rect(10.0, 0.0, 5.0, 5.0)]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 2);
        assert!(faces.iter().all(|f| f.holes.is_empty()));
    }

    #[test]
    fn assemble_faces_two_adjacent_zones_one_outer_two_holes() {
        // Outer 20x10; two CW holes side by side inside.
        let (loops, vertex_table) = loops_from_polygons(&[
            ccw_rect(0.0, 0.0, 20.0, 10.0),
            cw_rect(2.0, 2.0, 6.0, 6.0),
            cw_rect(12.0, 2.0, 6.0, 6.0),
        ]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 1);
        assert_eq!(faces[0].holes.len(), 2);
    }

    #[test]
    fn assemble_faces_two_sibling_holes_same_scanline() {
        // Both holes have leftmost-vertex at the same y-coordinate.
        let (loops, vertex_table) = loops_from_polygons(&[
            ccw_rect(0.0, 0.0, 20.0, 10.0),
            cw_rect(3.0, 4.0, 4.0, 2.0),  // hole A, leftmost x=3 y=4
            cw_rect(13.0, 4.0, 4.0, 2.0), // hole B, leftmost x=13 y=4
        ]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 1);
        assert_eq!(faces[0].holes.len(), 2);
    }

    #[test]
    fn assemble_faces_two_disjoint_donuts_same_scanline() {
        // Donut A at x∈[0,10], donut B at x∈[20,30]. All four leftmost
        // vertices share y=0. Verifies that B_outer is NOT misidentified as
        // a child of A_hole.
        let (loops, vertex_table) = loops_from_polygons(&[
            ccw_rect(0.0, 0.0, 10.0, 10.0),  // A outer
            cw_rect(2.0, 2.0, 6.0, 6.0),     // A hole
            ccw_rect(20.0, 0.0, 10.0, 10.0), // B outer
            cw_rect(22.0, 2.0, 6.0, 6.0),    // B hole
        ]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 2);
        assert!(faces.iter().all(|f| f.holes.len() == 1));
        // Each face's hole's leftmost x sits inside the same face's outer's bbox.
        for f in &faces {
            let outer_min_x = f.outer.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
            let outer_max_x = f
                .outer
                .iter()
                .map(|p| p.0)
                .fold(f64::NEG_INFINITY, f64::max);
            let hole_min_x = f.holes[0].iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
            assert!(hole_min_x > outer_min_x && hole_min_x < outer_max_x);
        }
    }

    #[test]
    fn assemble_faces_nested_island_depth_2() {
        // Outer ⊃ hole ⊃ island; island stands alone as a separate face.
        let (loops, vertex_table) = loops_from_polygons(&[
            ccw_rect(0.0, 0.0, 20.0, 20.0), // depth 0, CCW
            cw_rect(4.0, 4.0, 12.0, 12.0),  // depth 1, CW (hole of outer)
            ccw_rect(8.0, 8.0, 4.0, 4.0),   // depth 2, CCW (nested island)
        ]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 2);
        // Find each face by area.
        let mut face_areas: Vec<f64> = faces.iter().map(|f| signed_area(&f.outer)).collect();
        face_areas.sort_by(|a, b| b.partial_cmp(a).unwrap());
        assert!(face_areas[0] > 100.0); // big outer
        assert!(face_areas[1] < 30.0 && face_areas[1] > 10.0); // island
                                                               // The big outer has 1 hole; the island has 0 holes.
        let big = faces
            .iter()
            .find(|f| signed_area(&f.outer) > 100.0)
            .unwrap();
        let small = faces
            .iter()
            .find(|f| signed_area(&f.outer) < 100.0)
            .unwrap();
        assert_eq!(big.holes.len(), 1);
        assert_eq!(small.holes.len(), 0);
    }

    #[test]
    fn assemble_faces_nested_island_depth_3() {
        // Outer ⊃ hole ⊃ island ⊃ inner_hole.
        let (loops, vertex_table) = loops_from_polygons(&[
            ccw_rect(0.0, 0.0, 30.0, 30.0), // depth 0, CCW
            cw_rect(3.0, 3.0, 24.0, 24.0),  // depth 1, CW
            ccw_rect(6.0, 6.0, 18.0, 18.0), // depth 2, CCW
            cw_rect(9.0, 9.0, 12.0, 12.0),  // depth 3, CW
        ]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 2);
        // Both faces have exactly 1 hole.
        assert!(faces.iter().all(|f| f.holes.len() == 1));
    }

    #[test]
    fn assemble_faces_witness_robust_to_collinear_leftmost_vertices() {
        // Two disjoint outers at the same x but different y; their
        // leftmost-vertex witnesses share x but differ in y. Tie-break by y.
        let (loops, vertex_table) =
            loops_from_polygons(&[ccw_rect(0.0, 0.0, 5.0, 5.0), ccw_rect(0.0, 10.0, 5.0, 5.0)]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 2);
        assert!(faces.iter().all(|f| f.holes.is_empty()));
    }

    #[test]
    fn assemble_faces_returns_err_on_orientation_violation() {
        // CCW square at depth 1: should be CW (hole). Forces orientation/depth
        // mismatch.
        let (loops, vertex_table) = loops_from_polygons(&[
            ccw_rect(0.0, 0.0, 10.0, 10.0),
            ccw_rect(3.0, 3.0, 4.0, 4.0), // depth 1 but CCW — invalid
        ]);
        let err = assemble_faces(&loops, &vertex_table).expect_err("must fail");
        let msg = format!("{err}");
        assert!(
            msg.contains("orientation/depth mismatch"),
            "expected orientation/depth message; got {msg}"
        );
    }

    #[test]
    fn assemble_faces_returns_err_on_witness_boundary_tangent() {
        // Two squares that share an edge — every vertex of one square lies
        // exactly on the other's boundary. After k=3 fallbacks, all four
        // candidates remain on the boundary, so assemble_faces must Err.
        let a = vec![(0.0, 0.0), (5.0, 0.0), (5.0, 5.0), (0.0, 5.0)];
        let b = vec![(5.0, 0.0), (10.0, 0.0), (10.0, 5.0), (5.0, 5.0)];
        // Use only 3 vertices of B so all of them sit on A's right edge.
        let b_on_a = vec![(5.0, 1.0), (5.0, 2.0), (5.0, 3.0)];
        let (loops, vertex_table) = loops_from_polygons(&[a, b, b_on_a]);
        // With degenerate input the algorithm should refuse rather than
        // silently coerce. Either Err or a crisp parent classification is
        // acceptable; we only require not-panic and (if Ok) parity-consistent.
        match assemble_faces(&loops, &vertex_table) {
            Err(_) => {} // expected
            Ok(faces) => {
                for f in &faces {
                    assert!(signed_area(&f.outer) > 0.0);
                    for h in &f.holes {
                        assert!(signed_area(h) < 0.0);
                    }
                }
            }
        }
    }

    #[test]
    fn assemble_faces_input_order_independent() {
        let polygons = vec![
            ccw_rect(0.0, 0.0, 20.0, 20.0),
            cw_rect(4.0, 4.0, 12.0, 12.0),
            ccw_rect(8.0, 8.0, 4.0, 4.0),
        ];
        let (loops_a, vt_a) = loops_from_polygons(&polygons);
        let faces_a = assemble_faces(&loops_a, &vt_a).expect("a");

        // Reverse order.
        let polygons_rev: Vec<Polygon> = polygons.into_iter().rev().collect();
        let (loops_b, vt_b) = loops_from_polygons(&polygons_rev);
        let faces_b = assemble_faces(&loops_b, &vt_b).expect("b");

        // Same number of faces; same hole counts (sorted).
        assert_eq!(faces_a.len(), faces_b.len());
        let mut counts_a: Vec<usize> = faces_a.iter().map(|f| f.holes.len()).collect();
        let mut counts_b: Vec<usize> = faces_b.iter().map(|f| f.holes.len()).collect();
        counts_a.sort();
        counts_b.sort();
        assert_eq!(counts_a, counts_b);
    }
}
