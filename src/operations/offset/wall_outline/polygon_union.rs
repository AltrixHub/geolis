//! 2D polygon union via edge-graph construction.
//!
//! Each polygon is a CCW-wound simple polygon represented as `Vec<(f64, f64)>`.
//! The union of N polygons is computed by iterative pairwise union.
//! Result: CCW boundaries = outer shells, CW boundaries = holes.

use crate::math::distance_2d::point_to_segment_dist;

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

/// Result of a polygon union: outer boundaries + holes.
pub struct UnionResult {
    pub boundaries: Vec<Polygon>,
}

/// Computes the union of polygons-with-holes.
///
/// Flattens outers and holes into a single edge graph, classifies edges,
/// and traces boundaries. A sub-edge's midpoint is "inside" the union if
/// it is inside any outer AND outside all holes of at least one input polygon.
pub fn union_all_with_holes(inputs: &[PolygonWithHoles]) -> UnionResult {
    // Collect all rings: outers as filled, holes as voids
    let mut outers: Vec<(usize, &Polygon)> = Vec::new();
    let mut holes: Vec<(usize, &Polygon)> = Vec::new();
    for (i, pwh) in inputs.iter().enumerate() {
        outers.push((i, &pwh.outer));
        for hole in &pwh.holes {
            holes.push((i, hole));
        }
    }

    // Collect all rings as polygons for edge splitting
    let all_rings: Vec<&Polygon> = outers.iter().chain(holes.iter()).map(|(_, p)| *p).collect();

    // Build sub-edges from all rings
    let mut all_sub_edges: Vec<(usize, bool, (f64, f64), (f64, f64))> = Vec::new(); // (input_idx, is_outer, start, end)

    for (ring_idx, ring) in all_rings.iter().enumerate() {
        let (input_idx, is_outer) = if ring_idx < outers.len() {
            (outers[ring_idx].0, true)
        } else {
            (holes[ring_idx - outers.len()].0, false)
        };

        let n = ring.len();
        for ei in 0..n {
            let a0 = ring[ei];
            let a1 = ring[(ei + 1) % n];

            let mut params: Vec<f64> = vec![0.0, 1.0];
            for (other_idx, other_ring) in all_rings.iter().enumerate() {
                if other_idx == ring_idx {
                    continue;
                }
                let m = other_ring.len();
                for ej in 0..m {
                    let b0 = other_ring[ej];
                    let b1 = other_ring[(ej + 1) % m];
                    if let Some((t, _u)) = seg_seg_intersect(a0, a1, b0, b1) {
                        if t > WALL_EPS && t < 1.0 - WALL_EPS {
                            params.push(t);
                        }
                    } else {
                        // Parallel edges: `seg_seg_intersect` bails, but if
                        // the edges are also collinear and their projections
                        // overlap, we must still split at the other edge's
                        // endpoints so the shared portion deduplicates down
                        // to a single edge. Without this, adjacent rings
                        // that share a face leave both rings' edges intact
                        // and `trace_loops` cannot follow the combined
                        // boundary through the overlap.
                        for t in collinear_overlap_params(a0, a1, b0, b1) {
                            params.push(t);
                        }
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
                all_sub_edges.push((input_idx, is_outer, p0, p1));
            }
        }
    }

    if all_sub_edges.is_empty() {
        return UnionResult {
            boundaries: Vec::new(),
        };
    }

    // Classify: keep edges on the union boundary.
    // An edge is on the boundary if its midpoint is NOT strictly inside
    // any OTHER input's filled region (outer minus holes).
    all_sub_edges.retain(|&(_src, _is_outer, start, end)| {
        let mid = midpoint(start, end);
        // Remove edge if its midpoint is strictly inside ANY input's filled region.
        // Same-input is included: an outer edge's midpoint is Boundary on its own
        // outer (not Inside), so point_in_filled returns false. A hole edge's
        // midpoint is Boundary on its own hole, also handled correctly.
        for pwh in inputs {
            if point_in_filled(mid, pwh) {
                return false;
            }
        }
        true
    });

    let mut edges: Vec<((f64, f64), (f64, f64))> =
        all_sub_edges.iter().map(|&(_, _, s, e)| (s, e)).collect();
    dedup_edges(&mut edges);

    let mut boundaries = trace_loops(&edges);
    boundaries.retain(|p| signed_area(p).abs() > WALL_EPS_SQ);

    // Normalize winding: largest-area boundary is outer (CCW),
    // smaller boundaries contained in it are holes (CW).
    normalize_boundary_winding(&mut boundaries);

    UnionResult { boundaries }
}

fn point_in_filled(p: (f64, f64), pwh: &PolygonWithHoles) -> bool {
    if point_in_polygon_class(p, &pwh.outer) != PointClass::Inside {
        return false;
    }
    for hole in &pwh.holes {
        let c = point_in_polygon_class(p, hole);
        if c == PointClass::Inside || c == PointClass::Boundary {
            return false;
        }
    }
    true
}

/// Computes the union of all input polygons (simple, no holes).
/// Used by tests; production code uses `union_all_with_holes`.
#[cfg(test)]
fn union_all(polygons: &[Polygon]) -> UnionResult {
    if polygons.is_empty() {
        return UnionResult {
            boundaries: Vec::new(),
        };
    }
    if polygons.len() == 1 {
        return UnionResult {
            boundaries: vec![polygons[0].clone()],
        };
    }

    // Step 1: Split every polygon's edges at all intersection points
    // with every other polygon. Track source polygon index.
    let mut all_sub_edges: Vec<(usize, (f64, f64), (f64, f64))> = Vec::new();

    for (i, poly) in polygons.iter().enumerate() {
        let n = poly.len();
        for ei in 0..n {
            let a0 = poly[ei];
            let a1 = poly[(ei + 1) % n];

            let mut params: Vec<f64> = vec![0.0, 1.0];
            for (j, other) in polygons.iter().enumerate() {
                if i == j {
                    continue;
                }
                let m = other.len();
                for ej in 0..m {
                    let b0 = other[ej];
                    let b1 = other[(ej + 1) % m];
                    if let Some((t, _u)) = seg_seg_intersect(a0, a1, b0, b1) {
                        if t > WALL_EPS && t < 1.0 - WALL_EPS {
                            params.push(t);
                        }
                    }
                }
            }

            params.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
            params.dedup_by(|x, y| (*x - *y).abs() < WALL_EPS);

            for w in params.windows(2) {
                let t0 = w[0];
                let t1 = w[1];
                if (t1 - t0).abs() < WALL_EPS {
                    continue;
                }
                let p0 = lerp(a0, a1, t0);
                let p1 = lerp(a0, a1, t1);
                all_sub_edges.push((i, p0, p1));
            }
        }
    }

    if all_sub_edges.is_empty() {
        return UnionResult {
            boundaries: Vec::new(),
        };
    }

    // Step 2: Classify each sub-edge. Keep it if its midpoint is NOT
    // strictly inside any OTHER polygon (the edge is on the union boundary).
    all_sub_edges.retain(|&(src, start, end)| {
        let mid = midpoint(start, end);
        for (j, poly) in polygons.iter().enumerate() {
            if j == src {
                continue;
            }
            if point_in_polygon_class(mid, poly) == PointClass::Inside {
                return false;
            }
        }
        true
    });

    // Strip source index, deduplicate
    let mut edges: Vec<((f64, f64), (f64, f64))> =
        all_sub_edges.iter().map(|&(_, s, e)| (s, e)).collect();
    dedup_edges(&mut edges);

    // Step 3: Trace closed loops from remaining edges.
    let mut boundaries = trace_loops(&edges);

    // Remove zero-area results
    boundaries.retain(|p| signed_area(p).abs() > WALL_EPS_SQ);

    UnionResult { boundaries }
}

fn normalize_boundary_winding(boundaries: &mut [Polygon]) {
    if boundaries.is_empty() {
        return;
    }

    // Determine containment: for each boundary, find the smallest
    // enclosing boundary using point-in-polygon on the first vertex.
    let n = boundaries.len();
    let mut parent: Vec<Option<usize>> = vec![None; n];
    let areas: Vec<f64> = boundaries.iter().map(|b| signed_area(b).abs()).collect();

    for i in 0..n {
        let test_pt = boundaries[i][0];
        let mut best_parent: Option<usize> = None;
        let mut best_area = f64::MAX;
        for j in 0..n {
            if i == j {
                continue;
            }
            if areas[j] <= areas[i] {
                continue;
            }
            if point_in_polygon_class(test_pt, &boundaries[j]) == PointClass::Inside
                && areas[j] < best_area
            {
                best_parent = Some(j);
                best_area = areas[j];
            }
        }
        parent[i] = best_parent;
    }

    // Depth 0 (no parent) = outer → CCW
    // Depth 1 (has parent, parent has no parent) = hole → CW
    // Depth 2 = outer again, etc.
    for i in 0..n {
        let depth = containment_depth(i, &parent);
        if depth % 2 == 0 {
            if signed_area(&boundaries[i]) < 0.0 {
                boundaries[i].reverse();
            }
        } else if signed_area(&boundaries[i]) > 0.0 {
            boundaries[i].reverse();
        }
    }
}

fn containment_depth(idx: usize, parent: &[Option<usize>]) -> usize {
    let mut depth = 0;
    let mut current = parent[idx];
    while let Some(p) = current {
        depth += 1;
        current = parent[p];
    }
    depth
}

fn dedup_edges(edges: &mut Vec<((f64, f64), (f64, f64))>) {
    let mut i = 0;
    while i < edges.len() {
        let mut j = i + 1;
        let mut found_dup = false;
        while j < edges.len() {
            if (points_eq(edges[i].0, edges[j].0) && points_eq(edges[i].1, edges[j].1))
                || (points_eq(edges[i].0, edges[j].1) && points_eq(edges[i].1, edges[j].0))
            {
                edges.remove(j);
                found_dup = true;
                break;
            }
            j += 1;
        }
        if !found_dup {
            i += 1;
        }
    }
}

fn trace_loops(edges: &[((f64, f64), (f64, f64))]) -> Vec<Polygon> {
    if edges.is_empty() {
        return Vec::new();
    }

    // Build adjacency: point → list of (edge_idx, target_point)
    let mut points: Vec<(f64, f64)> = Vec::new();
    let mut edge_list: Vec<(usize, usize)> = Vec::new();

    for &(start, end) in edges {
        let si = ensure_point(&mut points, start);
        let ei = ensure_point(&mut points, end);
        if si != ei {
            edge_list.push((si, ei));
        }
    }

    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); points.len()];
    for (edge_idx, &(si, ei)) in edge_list.iter().enumerate() {
        adjacency[si].push((edge_idx, ei));
    }

    let mut used = vec![false; edge_list.len()];
    let mut results: Vec<Polygon> = Vec::new();

    loop {
        // Find an unused edge starting from the lowest-y point
        let start = find_unused_start(&edge_list, &used, &points);
        let Some(start_edge) = start else { break };

        let boundary = trace_one_loop(start_edge, &edge_list, &adjacency, &points, &mut used);
        if boundary.len() >= 3 {
            results.push(boundary);
        }
    }

    results
}

fn trace_one_loop(
    start_edge: usize,
    edge_list: &[(usize, usize)],
    adjacency: &[Vec<(usize, usize)>],
    points: &[(f64, f64)],
    used: &mut [bool],
) -> Polygon {
    let mut boundary: Polygon = Vec::new();
    let mut current_edge = start_edge;
    let start_point = edge_list[current_edge].0;

    for _ in 0..edge_list.len() + 1 {
        if used[current_edge] {
            break;
        }
        used[current_edge] = true;
        let (si, ei) = edge_list[current_edge];
        boundary.push(points[si]);

        // Pick next edge: minimum CCW angle
        let incoming_angle = (points[ei].1 - points[si].1).atan2(points[ei].0 - points[si].0);
        let reverse_angle = normalize_angle(incoming_angle + std::f64::consts::PI);

        let mut best: Option<(usize, f64)> = None;
        for &(next_edge, next_target) in &adjacency[ei] {
            if used[next_edge] {
                continue;
            }
            let next_angle =
                (points[next_target].1 - points[ei].1).atan2(points[next_target].0 - points[ei].0);
            let mut delta = normalize_angle(next_angle - reverse_angle);
            if delta.abs() < WALL_EPS {
                delta = 2.0 * std::f64::consts::PI;
            }
            if best.is_none() || delta < best.map_or(f64::MAX, |(_, d)| d) {
                best = Some((next_edge, delta));
            }
        }

        if let Some((next, _)) = best {
            current_edge = next;
            if edge_list[current_edge].0 == start_point {
                used[current_edge] = true;
                break;
            }
        } else {
            break;
        }
    }

    // Remove closing duplicate if present
    if boundary.len() > 1 && points_eq(boundary[0], *boundary.last().unwrap_or(&boundary[0])) {
        boundary.pop();
    }

    boundary
}

fn find_unused_start(
    edge_list: &[(usize, usize)],
    used: &[bool],
    points: &[(f64, f64)],
) -> Option<usize> {
    let mut best: Option<(usize, f64, f64)> = None;
    for (idx, &(si, _)) in edge_list.iter().enumerate() {
        if used[idx] {
            continue;
        }
        let p = points[si];
        if let Some((_, by, bx)) = best {
            if p.1 < by - WALL_EPS || ((p.1 - by).abs() < WALL_EPS && p.0 < bx) {
                best = Some((idx, p.1, p.0));
            }
        } else {
            best = Some((idx, p.1, p.0));
        }
    }
    best.map(|(idx, _, _)| idx)
}

// --- Geometry helpers ---

/// When two segments are parallel and collinear, return the parameter values
/// on segment A where segment B's endpoints project (clamped to the interior
/// of A). Returns an empty vec when the segments are not collinear or have no
/// interior overlap.
///
/// This complements `seg_seg_intersect` which bails on parallel segments —
/// adjacent rings that share a face produce parallel overlapping edges that
/// still need to be split so the overlapping portion can be deduplicated.
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

    // Parallel direction vectors.
    let dir_cross = d1x * d2y - d1y * d2x;
    if dir_cross.abs() >= WALL_EPS {
        return Vec::new();
    }
    // Collinear: b0 must lie on the infinite line through a0-a1.
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PointClass {
    Inside,
    Outside,
    Boundary,
}

fn point_in_polygon_class(p: (f64, f64), poly: &Polygon) -> PointClass {
    let n = poly.len();
    // First check if on any edge
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        let dist = point_to_segment_dist(p.0, p.1, a.0, a.1, b.0, b.1);
        if dist < WALL_EPS {
            return PointClass::Boundary;
        }
    }

    // Winding number
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

fn midpoint(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    ((a.0 + b.0) * 0.5, (a.1 + b.1) * 0.5)
}

fn points_eq(a: (f64, f64), b: (f64, f64)) -> bool {
    (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2) < WALL_EPS_SQ
}

fn normalize_angle(a: f64) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut r = a % two_pi;
    if r < 0.0 {
        r += two_pi;
    }
    r
}

fn ensure_point(points: &mut Vec<(f64, f64)>, p: (f64, f64)) -> usize {
    for (i, q) in points.iter().enumerate() {
        if points_eq(*q, p) {
            return i;
        }
    }
    points.push(p);
    points.len() - 1
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rect(x: f64, y: f64, w: f64, h: f64) -> Polygon {
        vec![(x, y), (x + w, y), (x + w, y + h), (x, y + h)]
    }

    #[test]
    fn union_non_overlapping() {
        let a = rect(0.0, 0.0, 2.0, 2.0);
        let b = rect(5.0, 0.0, 2.0, 2.0);
        let result = union_all(&[a, b]);
        assert_eq!(result.boundaries.len(), 2);
    }

    #[test]
    fn union_overlapping() {
        let a = rect(0.0, 0.0, 3.0, 2.0);
        let b = rect(2.0, 0.0, 3.0, 2.0);
        let result = union_all(&[a, b]);
        assert_eq!(result.boundaries.len(), 1);
        let area = signed_area(&result.boundaries[0]);
        assert!((area - 10.0).abs() < 0.1, "expected area ~10, got {area}");
    }

    #[test]
    fn union_shared_edge() {
        let a = rect(0.0, 0.0, 4.0, 3.0);
        let b = rect(4.0, 0.0, 4.0, 3.0);
        let result = union_all(&[a, b]);
        assert_eq!(result.boundaries.len(), 1);
        let area = signed_area(&result.boundaries[0]);
        assert!((area - 24.0).abs() < 0.1, "expected area ~24, got {area}");
    }

    #[test]
    fn union_contained() {
        let a = rect(0.0, 0.0, 6.0, 6.0);
        let b = rect(1.0, 1.0, 2.0, 2.0);
        let result = union_all(&[a, b]);
        assert_eq!(result.boundaries.len(), 1);
        let area = signed_area(&result.boundaries[0]);
        assert!((area - 36.0).abs() < 0.1, "expected area ~36, got {area}");
    }

    #[test]
    fn union_t_shape() {
        let horiz = rect(0.0, -1.0, 8.0, 2.0);
        let vert = rect(3.0, -1.0, 2.0, 5.0);
        let result = union_all(&[horiz, vert]);
        assert_eq!(result.boundaries.len(), 1);
        let expected_area = 8.0 * 2.0 + 2.0 * 5.0 - 2.0 * 2.0; // 16+10-4=22
        let area = signed_area(&result.boundaries[0]);
        assert!(
            (area - expected_area).abs() < 0.1,
            "expected area ~{expected_area}, got {area}"
        );
    }

    #[test]
    fn union_cross_shape() {
        let horiz = rect(0.0, 1.0, 6.0, 2.0);
        let vert = rect(2.0, 0.0, 2.0, 4.0);
        let result = union_all(&[horiz, vert]);
        assert_eq!(result.boundaries.len(), 1);
        let expected_area = 6.0 * 2.0 + 2.0 * 4.0 - 2.0 * 2.0;
        let area = signed_area(&result.boundaries[0]);
        assert!(
            (area - expected_area).abs() < 0.1,
            "expected area ~{expected_area}, got {area}"
        );
    }

    #[test]
    fn union_donut_from_four_rects() {
        // Four rectangles forming a closed square wall → should produce
        // an outer boundary (CCW) and a hole (CW). Exercises the production
        // `union_all_with_holes` path with hole-free inputs.
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
        let d = 0.3;
        let inputs: Vec<PolygonWithHoles> = vec![
            segment_to_rect((0.0, 0.0), (10.0, 0.0), d, d),
            segment_to_rect((10.0, 0.0), (10.0, 10.0), d, d),
            segment_to_rect((10.0, 10.0), (0.0, 10.0), d, d),
            segment_to_rect((0.0, 10.0), (0.0, 0.0), d, d),
        ]
        .into_iter()
        .map(|outer| PolygonWithHoles {
            outer,
            holes: Vec::new(),
        })
        .collect();
        let result = union_all_with_holes(&inputs);
        assert_eq!(
            result.boundaries.len(),
            2,
            "expected 2 boundaries (outer + hole)"
        );
        let areas: Vec<f64> = result.boundaries.iter().map(|b| signed_area(b)).collect();
        assert!(
            areas.iter().any(|a| *a > 0.0),
            "should have a CCW outer boundary"
        );
        assert!(areas.iter().any(|a| *a < 0.0), "should have a CW hole");
    }

    #[test]
    fn union_wall_segments_t_junction() {
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
        let d = 0.15;
        let rects: Vec<Polygon> = vec![
            segment_to_rect((0.0, 0.0), (4.0, 0.0), d, d),
            segment_to_rect((4.0, 0.0), (4.0, 3.0), d, d),
            segment_to_rect((4.0, 0.0), (8.0, 0.0), d, d),
        ];
        let result = union_all(&rects);
        assert!(!result.boundaries.is_empty());
        // All boundary vertices should be reasonable
        for b in &result.boundaries {
            for &(x, y) in b {
                assert!(x >= -0.5 && x <= 8.5, "x={x} out of range");
                assert!(y >= -0.5 && y <= 3.5, "y={y} out of range");
            }
        }
    }

    #[test]
    fn union_angled_wall_segments() {
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
        let d = 0.15;
        let rects: Vec<Polygon> = vec![
            segment_to_rect((-3.217, -4.144), (-2.635, 2.085), d, d),
            segment_to_rect((-3.217, -4.144), (2.002, -4.631), d, d),
            segment_to_rect((-2.635, 2.085), (2.578, 1.534), d, d),
            segment_to_rect((2.002, -4.631), (2.578, 1.534), d, d),
            segment_to_rect((2.002, -4.631), (6.473, -5.049), d, d),
            segment_to_rect((2.578, 1.534), (6.861, -0.896), d, d),
            segment_to_rect((6.473, -5.049), (6.861, -0.896), d, d),
        ];
        let result = union_all(&rects);
        assert!(!result.boundaries.is_empty());
        for b in &result.boundaries {
            for &(x, y) in b {
                assert!(
                    x >= -4.0 && x <= 8.0 && y >= -6.0 && y <= 3.0,
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
        let result = union_all_with_holes(&[pwh]);
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
        let result = union_all_with_holes(&[pwh]);
        assert_eq!(result.boundaries.len(), 2, "outer + hole");
        let areas: Vec<f64> = result.boundaries.iter().map(|b| signed_area(b)).collect();
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
        let result = union_all_with_holes(&[a, b]);
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
        let result = union_all_with_holes(&[a, b]);
        assert!(!result.boundaries.is_empty());
        // Should have outers and potentially holes
        let ccw_count = result
            .boundaries
            .iter()
            .filter(|b| signed_area(b) > 0.0)
            .count();
        assert!(ccw_count >= 1, "needs at least one outer");
    }

    #[test]
    fn with_holes_rings_sharing_a_colinear_face() {
        // Two zone-stroke rings whose outer faces lie on the same horizontal
        // line (y = 0): emulates two adjacent zones A and B stroked into
        // wall rings with a shared south wall segment. Before
        // `collinear_overlap_params` was introduced, the parallel overlap
        // left the trace_loops adjacency ambiguous and the output boundary
        // was fragmented. After the fix, the union yields a single outer
        // L/T-shape boundary that passes through the shared face cleanly.
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
        let result = union_all_with_holes(&[a, b]);
        // Expect one combined outer boundary + some interior holes for each
        // room. The critical property: every outer boundary vertex lies on
        // one of the four combined-shape sides (no stray vertices floating
        // off the rectangle).
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
        // Two stroke-expanded open walls forming a T-junction. The horizontal
        // wall's rectangle and the vertical stem rectangle overlap. The union
        // must collapse to a single T-shape boundary.
        let horiz = PolygonWithHoles {
            outer: rect(0.0, -0.15, 4.0, 0.30),
            holes: Vec::new(),
        };
        let vert = PolygonWithHoles {
            outer: rect(1.85, 0.0, 0.30, 3.0),
            holes: Vec::new(),
        };
        let result = union_all_with_holes(&[horiz, vert]);
        assert_eq!(
            result.boundaries.len(),
            1,
            "two overlapping wall strokes must union into one T boundary",
        );
    }

    #[test]
    fn with_holes_two_adjacent_zones_produce_one_outer_two_holes() {
        // Two adjacent zones with the SAME Y extent (pure T configuration).
        // Each zone's stroke is an annulus (outer rectangle + interior hole).
        // The combined geometry should be:
        //   - 1 outer rectangle (combined perimeter)
        //   - 2 holes (one per room)
        //
        // This mirrors `BoundarySolver` emitting one closed-ring `WallBaseline`
        // per zone into WallLayer's Rings slot.
        let d = 0.15;
        // Zone A: (0, 0) to (5, 3) → stroke outer (-d, -d) to (5+d, 3+d), hole CW.
        let a = PolygonWithHoles {
            outer: vec![(-d, -d), (5.0 + d, -d), (5.0 + d, 3.0 + d), (-d, 3.0 + d)],
            // Hole must be CW for polygon_union to treat it as a hole.
            holes: vec![vec![(d, d), (d, 3.0 - d), (5.0 - d, 3.0 - d), (5.0 - d, d)]],
        };
        // Zone B: (5, 0) to (8, 3) → stroke outer (5-d, -d) to (8+d, 3+d), hole CW.
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
        let result = union_all_with_holes(&[a, b]);
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
        assert_eq!(
            outers.len(),
            1,
            "two adjacent zones must produce exactly one outer boundary",
        );
        assert_eq!(
            holes.len(),
            2,
            "two adjacent zones must produce exactly two holes (one per room)",
        );
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
        let result = union_all_with_holes(&[a, b]);
        assert_eq!(result.boundaries.len(), 2, "two separate outers");
        let ccw_count = result
            .boundaries
            .iter()
            .filter(|b| signed_area(b) > 0.0)
            .count();
        assert_eq!(ccw_count, 2, "both should be CCW outers");
    }
}
