//! 2D arrangement engine — segment splitting, vertex snapping, half-edge
//! classification, face walking, and face assembly.
//!
//! The engine itself is **oracle-driven**. Operation modules (`union`,
//! `subtract`, ...) supply a [`FillOracle`] that decides whether a
//! probe point lies in the operation's result region. The engine then:
//!
//! 1. Collects raw segments from every input ring.
//! 2. Splits them into a planar arrangement (transverse crossings,
//!    T-junctions, collinear overlaps).
//! 3. Snaps endpoints to a `WALL_EPS` grid with bounded cross-cell
//!    merging (cluster diameter ≤ `2·WALL_EPS`).
//! 4. Canonicalises undirected sub-edges.
//! 5. Classifies each directed half-edge by bilateral perpendicular-ε
//!    sampling against the supplied oracle; keeps half-edges where the
//!    LEFT side is `Filled` AND the RIGHT side is `Empty`.
//! 6. Traces closed boundary loops by the polar-angle Δ successor rule.
//! 7. Assembles loops into face topology by containment matrix.
//!
//! The output is **always** the boundary that separates filled from
//! empty points for whatever oracle the caller supplied — that is the
//! whole point of the abstraction.

use std::collections::HashMap;

use crate::error::{OperationError, Result};

use super::types::{
    point_in_polygon_class, signed_area, PointClass, Polygon, PolygonWithHoles, WALL_EPS,
    WALL_EPS_SQ,
};

/// Identifies which ring of an input PWH a `Boundary` classification touched.
///
/// Diagnostic-only: surfaced in error messages when bilateral sampling
/// remains ambiguous after ε-shrink retries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BoundaryRef {
    Outer,
    Hole(usize),
}

/// Three-valued classification of a probe point against the operation's
/// result region.
///
/// `AmbiguousOnBoundary` is **never** silently folded into Filled or
/// Empty — it forces the engine to shrink ε and re-sample. `touched`
/// records every input ring the probe hit so the error message can
/// pinpoint the degeneracy.
#[derive(Debug, Clone)]
pub enum FilledClass {
    Filled,
    Empty,
    AmbiguousOnBoundary {
        #[allow(dead_code, reason = "diagnostic via Debug formatting")]
        touched: Vec<(usize, BoundaryRef)>,
    },
}

/// Trait for the fill oracle the engine consults during half-edge
/// classification. Implementors decide whether a probe point lies in
/// the operation's result region.
///
/// The engine guarantees `classify` is called only with probe points
/// produced by bilateral sampling around a half-edge midpoint — no
/// caller assumptions about probe distribution are required.
pub trait FillOracle {
    fn classify(&self, p: (f64, f64)) -> FilledClass;
}

/// OR-of-PWH-filled rule used by the union operation.
///
/// A point is `Filled` iff there exists at least one input PWH whose
/// outer contains it strictly Inside AND every hole has it strictly
/// Outside.
pub struct UnionOracle<'a> {
    pub inputs: &'a [PolygonWithHoles],
}

impl FillOracle for UnionOracle<'_> {
    fn classify(&self, p: (f64, f64)) -> FilledClass {
        let mut touched: Vec<(usize, BoundaryRef)> = Vec::new();
        let mut any_filled = false;
        for (idx, pwh) in self.inputs.iter().enumerate() {
            let (filled, ring_touches) = classify_pwh_filled(p, pwh, idx);
            touched.extend(ring_touches);
            if filled {
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
}

/// Subtract rule: `Filled` iff `base` is filled AND every entry in
/// `subtracts` is empty.
///
/// `base` is treated as `inputs[0]` and `subtracts` as `inputs[1..]`
/// for diagnostic indexing in `BoundaryRef::touched`, matching the
/// position those PWHs occupy in the segment-collection input list.
pub struct SubtractOracle<'a> {
    pub base: &'a PolygonWithHoles,
    pub subtracts: &'a [PolygonWithHoles],
}

impl FillOracle for SubtractOracle<'_> {
    fn classify(&self, p: (f64, f64)) -> FilledClass {
        let mut touched: Vec<(usize, BoundaryRef)> = Vec::new();
        let (base_filled, base_touches) = classify_pwh_filled(p, self.base, 0);
        touched.extend(base_touches);
        let mut any_subtract_filled = false;
        for (idx, pwh) in self.subtracts.iter().enumerate() {
            let (filled, ring_touches) = classify_pwh_filled(p, pwh, idx + 1);
            touched.extend(ring_touches);
            if filled {
                any_subtract_filled = true;
            }
        }
        if !touched.is_empty() {
            FilledClass::AmbiguousOnBoundary { touched }
        } else if base_filled && !any_subtract_filled {
            FilledClass::Filled
        } else {
            FilledClass::Empty
        }
    }
}

/// Two-valued PWH classification used internally by both oracles.
///
/// Returns `(filled, touched_rings)` where `filled` is true iff the
/// outer strictly contains the point AND every hole strictly excludes
/// it. `touched_rings` lists every ring whose classification returned
/// `Boundary`. If `touched_rings` is non-empty the boolean verdict is
/// unreliable and the engine must retry with a smaller ε.
fn classify_pwh_filled(
    p: (f64, f64),
    pwh: &PolygonWithHoles,
    input_idx: usize,
) -> (bool, Vec<(usize, BoundaryRef)>) {
    let mut touched: Vec<(usize, BoundaryRef)> = Vec::new();
    let outer_class = point_in_polygon_class(p, &pwh.outer);
    if outer_class == PointClass::Boundary {
        touched.push((input_idx, BoundaryRef::Outer));
    }
    let hole_classes: Vec<PointClass> = pwh
        .holes
        .iter()
        .map(|h| point_in_polygon_class(p, h))
        .collect();
    for (hi, hc) in hole_classes.iter().enumerate() {
        if *hc == PointClass::Boundary {
            touched.push((input_idx, BoundaryRef::Hole(hi)));
        }
    }
    let filled =
        outer_class == PointClass::Inside && hole_classes.iter().all(|c| *c == PointClass::Outside);
    (filled, touched)
}

/// Run the arrangement engine over `segment_inputs` and classify the
/// resulting half-edges with `oracle`.
///
/// `segment_inputs` is the set of PWHs whose rings are collected into
/// the raw segment list. For union this is identical to the oracle's
/// input set. For subtract it is `[base, subtracts...]` (base ∪
/// subtracts — every ring whose boundary may appear in the output).
///
/// # Errors
///
/// Returns [`OperationError::Failed`] if:
/// - Bilateral half-edge classification remains `AmbiguousOnBoundary`
///   after 3 ε-shrink retries.
/// - `assemble_faces` cannot pick a unique parent for a nested loop, or
///   detects an orientation/depth parity violation.
pub fn run_arrangement(
    segment_inputs: &[PolygonWithHoles],
    oracle: &impl FillOracle,
) -> Result<Vec<PolygonWithHoles>> {
    let raw_segments = collect_raw_segments(segment_inputs);
    if raw_segments.is_empty() {
        return Ok(Vec::new());
    }
    let split_segments = arrangement_split(&raw_segments);
    let (vertex_table, snapped_segments) = vertex_snap(&split_segments);
    let undirected = canonicalize_undirected(snapped_segments);
    let kept_half_edges = classify_and_filter(&undirected, &vertex_table, oracle)?;
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
    Ok(faces)
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
                    "boolean_2d post-condition: CDT vertex insert rejected \
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
                "boolean_2d post-condition: CDT constraint rejected \
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
/// that engine outputs are topologically identical regardless of input
/// order.
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
/// each by bilateral perpendicular-ε sampling against `oracle`, and
/// return only those half-edges with `Filled` on the LEFT and `Empty`
/// on the RIGHT.
fn classify_and_filter(
    undirected: &[(usize, usize)],
    vertex_table: &[(f64, f64)],
    oracle: &impl FillOracle,
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
                let l = oracle.classify((mid.0 + eps * nx, mid.1 + eps * ny));
                let r = oracle.classify((mid.0 - eps * nx, mid.1 - eps * ny));
                let l_amb = matches!(l, FilledClass::AmbiguousOnBoundary { .. });
                let r_amb = matches!(r, FilledClass::AmbiguousOnBoundary { .. });
                if !l_amb && !r_amb {
                    break (l, r);
                }
                tries += 1;
                if tries >= 3 {
                    return Err(OperationError::Failed(format!(
                        "boolean_2d: ambiguous half-edge classification at \
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

/// Builds the containment matrix on leftmost-vertex witnesses with bounded
/// `Boundary` fallback: `contained_in[i]` lists every loop that strictly
/// contains loop `i`.
///
/// # Errors
///
/// Returns [`OperationError::Failed`] when a loop's witness candidates are
/// exhausted (every tried witness lies on the other loop's boundary).
fn containment_matrix(
    polygons: &[Polygon],
    witness_candidates: &[Vec<(f64, f64)>],
) -> Result<Vec<Vec<usize>>> {
    const MAX_WITNESS_FALLBACK: usize = 3;
    let n = polygons.len();
    let mut contained_in: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for (j, polygon) in polygons.iter().enumerate() {
            if i == j {
                continue;
            }
            let mut witness_idx = 0;
            loop {
                if witness_idx >= MAX_WITNESS_FALLBACK || witness_idx >= witness_candidates[i].len()
                {
                    return Err(OperationError::Failed(format!(
                        "boolean_2d: assemble_faces witness disambiguation \
                         exhausted for loop {i} against loop {j} (tried \
                         {witness_idx} candidates, all on j's boundary)"
                    ))
                    .into());
                }
                let w = witness_candidates[i][witness_idx];
                match point_in_polygon_class(w, polygon) {
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
    Ok(contained_in)
}

/// Computes each loop's parent: the unique containing loop at `depth − 1`.
///
/// # Errors
///
/// Returns [`OperationError::Failed`] when the parent candidate at
/// `depth − 1` is not unique (broken arrangement topology).
fn parent_loops(contained_in: &[Vec<usize>], depth: &[usize]) -> Result<Vec<Option<usize>>> {
    let n = depth.len();
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
                "boolean_2d: assemble_faces found {} parent candidates for \
                 loop {} (depth {}); arrangement topology is broken",
                candidates.len(),
                i,
                depth[i],
            ))
            .into());
        }
        parent[i] = Some(candidates[0]);
    }
    Ok(parent)
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
    let contained_in = containment_matrix(&polygons, &witness_candidates)?;

    // 3. Depth and parent.
    let depth: Vec<usize> = contained_in.iter().map(Vec::len).collect();
    let parent = parent_loops(&contained_in, &depth)?;

    // 4. Validate orientation parity.
    for i in 0..n {
        let want_ccw = depth[i].is_multiple_of(2);
        let is_ccw = orientation[i] == Orientation::Ccw;
        if want_ccw != is_ccw {
            return Err(OperationError::Failed(format!(
                "boolean_2d: assemble_faces orientation/depth mismatch on \
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

/// Segment-segment intersection in parameter space. Returns
/// `Some((t, u))` with `t ∈ [0, 1]` on segment A and `u ∈ [0, 1]` on
/// segment B when the two segments cross transversely.
///
/// Returns `None` for parallel or non-intersecting pairs. Collinear
/// overlap is **not** reported here — use [`collinear_overlap_params`]
/// for that case.
pub(crate) fn seg_seg_intersect(
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

fn lerp(a: (f64, f64), b: (f64, f64), t: f64) -> (f64, f64) {
    (a.0 + t * (b.0 - a.0), a.1 + t * (b.1 - a.1))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    //! Engine-internal tests (`assemble_faces` synthetic loop fixtures).
    //!
    //! The end-to-end union / subtract behaviour is exercised by tests
    //! in the `union` / `subtract` submodules; this file only covers
    //! engine internals that take pre-built `WalkedLoop` inputs.

    use super::*;

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
        let (loops, vertex_table) = loops_from_polygons(&[
            ccw_rect(0.0, 0.0, 20.0, 10.0),
            cw_rect(3.0, 4.0, 4.0, 2.0),
            cw_rect(13.0, 4.0, 4.0, 2.0),
        ]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 1);
        assert_eq!(faces[0].holes.len(), 2);
    }

    #[test]
    fn assemble_faces_two_disjoint_donuts_same_scanline() {
        let (loops, vertex_table) = loops_from_polygons(&[
            ccw_rect(0.0, 0.0, 10.0, 10.0),
            cw_rect(2.0, 2.0, 6.0, 6.0),
            ccw_rect(20.0, 0.0, 10.0, 10.0),
            cw_rect(22.0, 2.0, 6.0, 6.0),
        ]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 2);
        assert!(faces.iter().all(|f| f.holes.len() == 1));
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
        let (loops, vertex_table) = loops_from_polygons(&[
            ccw_rect(0.0, 0.0, 20.0, 20.0),
            cw_rect(4.0, 4.0, 12.0, 12.0),
            ccw_rect(8.0, 8.0, 4.0, 4.0),
        ]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 2);
        let mut face_areas: Vec<f64> = faces.iter().map(|f| signed_area(&f.outer)).collect();
        face_areas.sort_by(|a, b| b.partial_cmp(a).unwrap());
        assert!(face_areas[0] > 100.0);
        assert!(face_areas[1] < 30.0 && face_areas[1] > 10.0);
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
        let (loops, vertex_table) = loops_from_polygons(&[
            ccw_rect(0.0, 0.0, 30.0, 30.0),
            cw_rect(3.0, 3.0, 24.0, 24.0),
            ccw_rect(6.0, 6.0, 18.0, 18.0),
            cw_rect(9.0, 9.0, 12.0, 12.0),
        ]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 2);
        assert!(faces.iter().all(|f| f.holes.len() == 1));
    }

    #[test]
    fn assemble_faces_witness_robust_to_collinear_leftmost_vertices() {
        let (loops, vertex_table) =
            loops_from_polygons(&[ccw_rect(0.0, 0.0, 5.0, 5.0), ccw_rect(0.0, 10.0, 5.0, 5.0)]);
        let faces = assemble_faces(&loops, &vertex_table).expect("assemble must succeed");
        assert_eq!(faces.len(), 2);
        assert!(faces.iter().all(|f| f.holes.is_empty()));
    }

    #[test]
    fn assemble_faces_returns_err_on_orientation_violation() {
        let (loops, vertex_table) =
            loops_from_polygons(&[ccw_rect(0.0, 0.0, 10.0, 10.0), ccw_rect(3.0, 3.0, 4.0, 4.0)]);
        let err = assemble_faces(&loops, &vertex_table).expect_err("must fail");
        let msg = format!("{err}");
        assert!(
            msg.contains("orientation/depth mismatch"),
            "expected orientation/depth message; got {msg}"
        );
    }

    #[test]
    fn assemble_faces_returns_err_on_witness_boundary_tangent() {
        let a = vec![(0.0, 0.0), (5.0, 0.0), (5.0, 5.0), (0.0, 5.0)];
        let b = vec![(5.0, 0.0), (10.0, 0.0), (10.0, 5.0), (5.0, 5.0)];
        let b_on_a = vec![(5.0, 1.0), (5.0, 2.0), (5.0, 3.0)];
        let (loops, vertex_table) = loops_from_polygons(&[a, b, b_on_a]);
        match assemble_faces(&loops, &vertex_table) {
            Err(_) => {}
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

        let polygons_rev: Vec<Polygon> = polygons.into_iter().rev().collect();
        let (loops_b, vt_b) = loops_from_polygons(&polygons_rev);
        let faces_b = assemble_faces(&loops_b, &vt_b).expect("b");

        assert_eq!(faces_a.len(), faces_b.len());
        let mut counts_a: Vec<usize> = faces_a.iter().map(|f| f.holes.len()).collect();
        let mut counts_b: Vec<usize> = faces_b.iter().map(|f| f.holes.len()).collect();
        counts_a.sort_unstable();
        counts_b.sort_unstable();
        assert_eq!(counts_a, counts_b);
    }
}
