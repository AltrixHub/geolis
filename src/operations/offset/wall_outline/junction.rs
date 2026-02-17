use crate::math::intersect_2d::segment_segment_intersect_2d;
use crate::math::{Point3, TOLERANCE};

use super::decompose::UniqueSegment;

/// A node in the centerline network — either a junction (intersection)
/// or a dead end (segment endpoint with no neighbor).
#[derive(Debug, Clone)]
pub struct Node {
    #[allow(dead_code)]
    pub point: (f64, f64),
    pub kind: NodeKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeKind {
    Junction,
    Interior,
    DeadEnd,
}

/// A sub-segment produced by splitting unique segments at junctions.
#[derive(Debug, Clone)]
pub struct SubSegment {
    pub start: (f64, f64),
    pub end: (f64, f64),
    /// Index of start node.
    pub start_node: usize,
    /// Index of end node.
    pub end_node: usize,
}

/// The centerline network: nodes + sub-segments.
#[derive(Debug)]
pub struct Network {
    pub nodes: Vec<Node>,
    pub sub_segments: Vec<SubSegment>,
}

/// Builds a network from unique segments by detecting junctions and splitting.
pub fn build_network(segments: &[UniqueSegment]) -> Network {
    // Step 1: Find all junction points (intersections between segments).
    let mut junction_points: Vec<(f64, f64)> = Vec::new();

    for i in 0..segments.len() {
        for j in (i + 1)..segments.len() {
            let a0 = Point3::new(segments[i].start.0, segments[i].start.1, 0.0);
            let a1 = Point3::new(segments[i].end.0, segments[i].end.1, 0.0);
            let b0 = Point3::new(segments[j].start.0, segments[j].start.1, 0.0);
            let b1 = Point3::new(segments[j].end.0, segments[j].end.1, 0.0);

            if let Some((pt, _t, _u)) = segment_segment_intersect_2d(&a0, &a1, &b0, &b1) {
                add_unique_point(&mut junction_points, (pt.x, pt.y));
            }
        }
    }

    // Step 2: Split each segment at junction points.
    let mut all_nodes: Vec<(f64, f64)> = Vec::new();
    let mut sub_segments: Vec<SubSegment> = Vec::new();

    for seg in segments {
        // Collect junction points that lie on this segment.
        let mut splits: Vec<f64> = Vec::new();
        let dx = seg.end.0 - seg.start.0;
        let dy = seg.end.1 - seg.start.1;

        for &jp in &junction_points {
            let t = project_on_segment(seg, jp);
            if t > TOLERANCE * 10.0 && t < 1.0 - TOLERANCE * 10.0 {
                // Verify the junction point is actually on the segment
                // (not just a parameter match with large perpendicular distance).
                let foot = (seg.start.0 + dx * t, seg.start.1 + dy * t);
                let dist_sq = (foot.0 - jp.0).powi(2) + (foot.1 - jp.1).powi(2);
                if dist_sq < TOLERANCE * 100.0 {
                    splits.push(t);
                }
            }
        }

        // Sort split parameters.
        splits.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        splits.dedup_by(|a, b| (*a - *b).abs() < TOLERANCE * 10.0);

        // Build sub-segments.
        let prev_node = ensure_node(&mut all_nodes, seg.start);

        let mut current_start = seg.start;
        let mut current_start_node = prev_node;

        for &t in &splits {
            let mid = (
                seg.start.0 + dx * t,
                seg.start.1 + dy * t,
            );
            let mid_node = ensure_node(&mut all_nodes, mid);
            sub_segments.push(SubSegment {
                start: current_start,
                end: mid,
                start_node: current_start_node,
                end_node: mid_node,
            });
            current_start = mid;
            current_start_node = mid_node;
        }

        // Final sub-segment to the end.
        let end_node = ensure_node(&mut all_nodes, seg.end);
        sub_segments.push(SubSegment {
            start: current_start,
            end: seg.end,
            start_node: current_start_node,
            end_node,
        });
    }

    // Step 3: Compute valence (number of connected sub-segments) for each node.
    let mut valence = vec![0_usize; all_nodes.len()];
    for ss in &sub_segments {
        valence[ss.start_node] += 1;
        valence[ss.end_node] += 1;
    }

    // Step 4: Classify nodes by junction status and valence.
    let nodes: Vec<Node> = all_nodes
        .iter()
        .enumerate()
        .map(|(i, &p)| {
            let is_junction = point_is_junction(&junction_points, p);
            let kind = if is_junction || valence[i] >= 3 {
                NodeKind::Junction
            } else if valence[i] == 2 {
                NodeKind::Interior
            } else {
                NodeKind::DeadEnd
            };
            Node { point: p, kind }
        })
        .collect();

    Network {
        nodes,
        sub_segments,
    }
}

/// Projects a point onto a segment, returning parameter t in [0, 1].
fn project_on_segment(seg: &UniqueSegment, p: (f64, f64)) -> f64 {
    let dx = seg.end.0 - seg.start.0;
    let dy = seg.end.1 - seg.start.1;
    let len_sq = dx * dx + dy * dy;
    if len_sq < TOLERANCE * TOLERANCE {
        return 0.0;
    }
    let t = ((p.0 - seg.start.0) * dx + (p.1 - seg.start.1) * dy) / len_sq;
    t.clamp(0.0, 1.0)
}

/// Adds a point to the list if not already present (within tolerance).
fn add_unique_point(points: &mut Vec<(f64, f64)>, p: (f64, f64)) {
    let tol_sq = TOLERANCE * 100.0;
    let exists = points
        .iter()
        .any(|q| (q.0 - p.0).powi(2) + (q.1 - p.1).powi(2) < tol_sq);
    if !exists {
        points.push(p);
    }
}

/// Finds or inserts a node, returning its index.
fn ensure_node(nodes: &mut Vec<(f64, f64)>, p: (f64, f64)) -> usize {
    let tol_sq = TOLERANCE * 100.0;
    for (i, n) in nodes.iter().enumerate() {
        if (n.0 - p.0).powi(2) + (n.1 - p.1).powi(2) < tol_sq {
            return i;
        }
    }
    nodes.push(p);
    nodes.len() - 1
}

/// Checks if a point is in the junction list.
fn point_is_junction(junctions: &[(f64, f64)], p: (f64, f64)) -> bool {
    let tol_sq = TOLERANCE * 100.0;
    junctions
        .iter()
        .any(|j| (j.0 - p.0).powi(2) + (j.1 - p.1).powi(2) < tol_sq)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn double_cross_network() {
        // 4 unique segments for the 井 shape.
        let segments = vec![
            UniqueSegment { start: (3.0, 0.0), end: (3.0, 10.0) },
            UniqueSegment { start: (0.0, 7.0), end: (10.0, 7.0) },
            UniqueSegment { start: (7.0, 0.0), end: (7.0, 10.0) },
            UniqueSegment { start: (0.0, 3.0), end: (10.0, 3.0) },
        ];
        let net = build_network(&segments);

        // 4 junctions: (3,3), (3,7), (7,3), (7,7)
        let junction_count = net.nodes.iter().filter(|n| n.kind == NodeKind::Junction).count();
        assert_eq!(junction_count, 4, "expected 4 junctions, got {junction_count}");

        // 8 dead ends: endpoints of the 4 segments.
        let dead_end_count = net.nodes.iter().filter(|n| n.kind == NodeKind::DeadEnd).count();
        assert_eq!(dead_end_count, 8, "expected 8 dead ends, got {dead_end_count}");

        // 12 sub-segments: 4 segments × 3 sub-segments each.
        assert_eq!(net.sub_segments.len(), 12, "expected 12 sub-segments, got {}", net.sub_segments.len());
    }

    #[test]
    fn closed_square_no_dead_ends() {
        // 4-segment closed square: corners are junctions (shared endpoints),
        // no dead ends.
        let segments = vec![
            UniqueSegment { start: (0.0, 0.0), end: (10.0, 0.0) },
            UniqueSegment { start: (10.0, 0.0), end: (10.0, 10.0) },
            UniqueSegment { start: (10.0, 10.0), end: (0.0, 10.0) },
            UniqueSegment { start: (0.0, 10.0), end: (0.0, 0.0) },
        ];
        let net = build_network(&segments);

        assert_eq!(net.nodes.len(), 4, "expected 4 nodes, got {}", net.nodes.len());
        assert_eq!(net.sub_segments.len(), 4, "expected 4 sub-segments, got {}", net.sub_segments.len());

        let dead_end_count = net.nodes.iter().filter(|n| n.kind == NodeKind::DeadEnd).count();
        assert_eq!(dead_end_count, 0, "closed square should have no dead ends, got {dead_end_count}");
    }

    #[test]
    fn single_line_no_junctions() {
        let segments = vec![
            UniqueSegment { start: (0.0, 0.0), end: (5.0, 0.0) },
        ];
        let net = build_network(&segments);
        assert_eq!(net.nodes.len(), 2);
        assert_eq!(net.sub_segments.len(), 1);
        assert!(net.nodes.iter().all(|n| n.kind == NodeKind::DeadEnd));
    }

    #[test]
    fn l_shape_one_junction() {
        let segments = vec![
            UniqueSegment { start: (0.0, 0.0), end: (5.0, 0.0) },
            UniqueSegment { start: (5.0, 0.0), end: (5.0, 5.0) },
        ];
        let net = build_network(&segments);

        // (5,0) is at the endpoint of both segments — it's an endpoint intersection.
        // segment_segment_intersect_2d includes endpoints.
        // The junction at (5,0) is at t=1.0 for seg0 and t=0.0 for seg1,
        // so it won't be split (it's already an endpoint).
        assert_eq!(net.sub_segments.len(), 2, "expected 2 sub-segments, got {}", net.sub_segments.len());
    }
}
