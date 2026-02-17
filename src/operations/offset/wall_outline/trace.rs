use crate::geometry::pline::{Pline, PlineVertex};
use crate::math::TOLERANCE;

use super::offset_edges::OffsetEdge;

/// Traces outer boundaries from a set of offset edges.
///
/// Uses an angle-based walk: starting from the bottom-left point,
/// always turns to the minimum CCW angle at each junction.
///
/// Returns one or more closed `Pline` boundaries.
#[allow(clippy::too_many_lines)]
pub fn trace_boundaries(edges: &[OffsetEdge]) -> Vec<Pline> {
    if edges.is_empty() {
        return Vec::new();
    }

    let tol_sq = TOLERANCE * 1e4;
    let mut points: Vec<(f64, f64)> = Vec::new();
    let mut edge_list: Vec<(usize, usize)> = Vec::new();

    for e in edges {
        let si = ensure_point(&mut points, e.start, tol_sq);
        let ei = ensure_point(&mut points, e.end, tol_sq);
        if si != ei {
            edge_list.push((si, ei));
        }
    }

    if edge_list.is_empty() {
        return Vec::new();
    }

    // Build adjacency list: point_idx -> list of (edge_idx, target_point_idx).
    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); points.len()];
    for (edge_idx, &(si, ei)) in edge_list.iter().enumerate() {
        adjacency[si].push((edge_idx, ei));
    }

    let mut used = vec![false; edge_list.len()];
    let mut results: Vec<Pline> = Vec::new();

    loop {
        let start_edge_idx = find_start_edge(&edge_list, &used, &points);
        let Some(start_edge_idx) = start_edge_idx else {
            break;
        };

        let boundary = trace_one_boundary(
            start_edge_idx,
            &edge_list,
            &adjacency,
            &points,
            &mut used,
            tol_sq,
        );

        if boundary.len() >= 3 {
            let vertices: Vec<PlineVertex> = boundary
                .iter()
                .map(|&(x, y)| PlineVertex::line(x, y))
                .collect();
            results.push(Pline {
                vertices,
                closed: true,
            });
        }
    }

    results
}

/// Finds the next unused edge starting from the point with lowest y (then x).
fn find_start_edge(
    edge_list: &[(usize, usize)],
    used: &[bool],
    points: &[(f64, f64)],
) -> Option<usize> {
    let mut best: Option<(usize, f64, f64)> = None;
    for (edge_idx, &(si, _)) in edge_list.iter().enumerate() {
        if used[edge_idx] {
            continue;
        }
        let p = points[si];
        if let Some((_, by, bx)) = best {
            if p.1 < by - TOLERANCE || ((p.1 - by).abs() < TOLERANCE && p.0 < bx) {
                best = Some((edge_idx, p.1, p.0));
            }
        } else {
            best = Some((edge_idx, p.1, p.0));
        }
    }
    best.map(|(idx, _, _)| idx)
}

/// Traces a single boundary loop starting from the given edge.
fn trace_one_boundary(
    start_edge_idx: usize,
    edge_list: &[(usize, usize)],
    adjacency: &[Vec<(usize, usize)>],
    points: &[(f64, f64)],
    used: &mut [bool],
    tol_sq: f64,
) -> Vec<(f64, f64)> {
    let mut boundary: Vec<(f64, f64)> = Vec::new();
    let mut current_edge = start_edge_idx;
    let start_point = edge_list[current_edge].0;

    loop {
        if used[current_edge] {
            break;
        }
        used[current_edge] = true;
        let (si, ei) = edge_list[current_edge];
        boundary.push(points[si]);

        let incoming_angle =
            (points[ei].1 - points[si].1).atan2(points[ei].0 - points[si].0);

        let best_next = pick_next_edge(ei, incoming_angle, adjacency, points, used);

        if let Some(next) = best_next {
            current_edge = next;
            if edge_list[current_edge].0 == start_point {
                used[current_edge] = true;
                boundary.push(points[edge_list[current_edge].0]);
                // Remove closing duplicate.
                if boundary.len() > 1 {
                    let first = boundary[0];
                    let last = boundary[boundary.len() - 1];
                    if (first.0 - last.0).powi(2) + (first.1 - last.1).powi(2) < tol_sq {
                        boundary.pop();
                    }
                }
                break;
            }
        } else {
            break;
        }
    }

    boundary
}

/// Picks the next edge from a node using minimum-CCW-angle selection.
fn pick_next_edge(
    node: usize,
    incoming_angle: f64,
    adjacency: &[Vec<(usize, usize)>],
    points: &[(f64, f64)],
    used: &[bool],
) -> Option<usize> {
    let reverse_angle = normalize_angle(incoming_angle + std::f64::consts::PI);
    let mut best: Option<(usize, f64)> = None;

    for &(next_edge_idx, next_target) in &adjacency[node] {
        if used[next_edge_idx] {
            continue;
        }
        let next_angle = (points[next_target].1 - points[node].1)
            .atan2(points[next_target].0 - points[node].0);

        let mut delta = normalize_angle(next_angle - reverse_angle);
        if delta.abs() < TOLERANCE {
            delta = 2.0 * std::f64::consts::PI;
        }

        if best.is_none() || delta < best.map_or(f64::MAX, |(_, bd)| bd) {
            best = Some((next_edge_idx, delta));
        }
    }

    best.map(|(idx, _)| idx)
}

/// Normalizes an angle to [0, 2pi).
fn normalize_angle(a: f64) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut r = a % two_pi;
    if r < 0.0 {
        r += two_pi;
    }
    r
}

/// Finds or inserts a point, returning its index.
fn ensure_point(points: &mut Vec<(f64, f64)>, p: (f64, f64), tol_sq: f64) -> usize {
    for (i, q) in points.iter().enumerate() {
        if (q.0 - p.0).powi(2) + (q.1 - p.1).powi(2) < tol_sq {
            return i;
        }
    }
    points.push(p);
    points.len() - 1
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn simple_rectangle_trace() {
        let edges = vec![
            OffsetEdge { start: (0.0, 0.0), end: (5.0, 0.0) },
            OffsetEdge { start: (5.0, 0.0), end: (5.0, 3.0) },
            OffsetEdge { start: (5.0, 3.0), end: (0.0, 3.0) },
            OffsetEdge { start: (0.0, 3.0), end: (0.0, 0.0) },
        ];
        let result = trace_boundaries(&edges);
        assert_eq!(result.len(), 1, "expected 1 boundary");
        assert_eq!(result[0].vertices.len(), 4, "expected 4 vertices");
        assert!(result[0].closed);
    }
}
