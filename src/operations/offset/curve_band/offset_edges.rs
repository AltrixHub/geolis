use crate::math::intersect_2d::line_line_intersect_2d;
use crate::math::{Point3, Vector3, TOLERANCE};

use super::junction::{Network, NodeKind, SubSegment};

/// An offset edge segment — the fundamental unit for boundary tracing.
#[derive(Debug, Clone)]
pub struct OffsetEdge {
    pub start: (f64, f64),
    pub end: (f64, f64),
}

/// Builds offset edges from a centerline network at the given widths.
///
/// `left_width` is the offset distance to the left of the segment direction;
/// `right_width` is the offset distance to the right.  Pass equal values for a
/// centred wall; pass `(0, thickness)` or `(thickness, 0)` for a wall that
/// extends only to one side of the baseline.
///
/// Edge orientation convention (produces CW outer boundary):
/// - Right offset edges go FORWARD (segment start to end direction)
/// - Left offset edges go BACKWARD (segment end to start direction)
/// - Start caps: left start to right start
/// - End caps: right end to left end
pub fn build(network: &Network, left_width: f64, right_width: f64) -> Vec<OffsetEdge> {
    let sub_segs = &network.sub_segments;

    // Precompute offset line data for each sub-segment.
    let offset_data: Vec<OffsetLineData> = sub_segs
        .iter()
        .map(|ss| compute_offset_lines(ss, left_width, right_width))
        .collect();

    // Resolve endpoint positions at junctions.
    let max_width = left_width.max(right_width);
    let resolved = resolve_all_endpoints(network, &offset_data, max_width);

    let mut result: Vec<OffsetEdge> = Vec::new();

    // Left offset edge (backward) and right offset edge (forward) for each sub-segment.
    for (i, _ss) in sub_segs.iter().enumerate() {
        let r = &resolved[i];
        // Left edge goes backward: end → start.
        result.push(OffsetEdge {
            start: r.left_end,
            end: r.left_start,
        });
        // Right edge goes forward: start → end.
        result.push(OffsetEdge {
            start: r.right_start,
            end: r.right_end,
        });
    }

    // Dead-end cap edges.
    for (node_idx, node) in network.nodes.iter().enumerate() {
        if node.kind != NodeKind::DeadEnd {
            continue;
        }
        for (seg_idx, ss) in sub_segs.iter().enumerate() {
            let r = &resolved[seg_idx];
            if ss.start_node == node_idx {
                // Start cap: left_start → right_start.
                result.push(OffsetEdge {
                    start: r.left_start,
                    end: r.right_start,
                });
            }
            if ss.end_node == node_idx {
                // End cap: right_end → left_end.
                result.push(OffsetEdge {
                    start: r.right_end,
                    end: r.left_end,
                });
            }
        }
    }

    result
}

#[derive(Debug, Clone)]
struct OffsetLineData {
    dir: (f64, f64),
    left_start: (f64, f64),
    left_end: (f64, f64),
    right_start: (f64, f64),
    right_end: (f64, f64),
}

#[derive(Debug, Clone)]
struct ResolvedEndpoints {
    left_start: (f64, f64),
    left_end: (f64, f64),
    right_start: (f64, f64),
    right_end: (f64, f64),
}

fn compute_offset_lines(ss: &SubSegment, left_width: f64, right_width: f64) -> OffsetLineData {
    let dx = ss.end.0 - ss.start.0;
    let dy = ss.end.1 - ss.start.1;
    let len = (dx * dx + dy * dy).sqrt();
    let (nx, ny) = if len > TOLERANCE {
        (dx / len, dy / len)
    } else {
        (1.0, 0.0)
    };

    // Left normal: (-ny, nx)
    let ln = (-ny, nx);

    let left_start = (
        ss.start.0 + left_width * ln.0,
        ss.start.1 + left_width * ln.1,
    );
    let left_end = (ss.end.0 + left_width * ln.0, ss.end.1 + left_width * ln.1);
    let right_start = (
        ss.start.0 - right_width * ln.0,
        ss.start.1 - right_width * ln.1,
    );
    let right_end = (ss.end.0 - right_width * ln.0, ss.end.1 - right_width * ln.1);

    OffsetLineData {
        dir: (nx, ny),
        left_start,
        left_end,
        right_start,
        right_end,
    }
}

fn neg_dir(dir: (f64, f64)) -> (f64, f64) {
    (-dir.0, -dir.1)
}

/// Minimum absolute sine of the angle between two arm directions below
/// which they are considered parallel and merged into one group.
/// sin(10°) ≈ 0.174. Generous threshold to catch nearly-collinear segments
/// from different polylines sharing an endpoint.
const PARALLEL_SIN: f64 = 0.174;

/// Resolves all offset edge endpoints by computing junction corner
/// intersections.
///
/// Nearly-parallel arms are merged into groups before computing corners.
/// For each pair of adjacent groups in CCW order, we compute the
/// intersection of the representative LEFT and RIGHT offset lines.
/// All arms within a group receive the same resolved corner.
fn resolve_all_endpoints(
    network: &Network,
    offset_data: &[OffsetLineData],
    _max_width: f64,
) -> Vec<ResolvedEndpoints> {
    let sub_segs = &network.sub_segments;
    let mut resolved: Vec<ResolvedEndpoints> = offset_data
        .iter()
        .map(|od| ResolvedEndpoints {
            left_start: od.left_start,
            left_end: od.left_end,
            right_start: od.right_start,
            right_end: od.right_end,
        })
        .collect();

    for (node_idx, node) in network.nodes.iter().enumerate() {
        if node.kind == NodeKind::DeadEnd {
            continue;
        }

        // Collect connected arms with their outgoing angle from the junction.
        let mut arms: Vec<(f64, usize, bool)> = Vec::new();
        for (seg_idx, ss) in sub_segs.iter().enumerate() {
            if ss.start_node == node_idx {
                let (dx, dy) = (ss.end.0 - ss.start.0, ss.end.1 - ss.start.1);
                arms.push((dy.atan2(dx), seg_idx, true));
            }
            if ss.end_node == node_idx {
                let (dx, dy) = (ss.start.0 - ss.end.0, ss.start.1 - ss.end.1);
                arms.push((dy.atan2(dx), seg_idx, false));
            }
        }
        arms.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        let n = arms.len();
        if n < 2 {
            continue;
        }

        // Group nearly-parallel adjacent arms. Each group is represented
        // by its first arm for LEFT and last arm for RIGHT offset lines.
        let groups = merge_parallel_arms(&arms, offset_data);
        let gn = groups.len();

        if gn < 2 {
            continue;
        }

        for k in 0..gn {
            let group_k = &groups[k];
            let group_next = &groups[(k + 1) % gn];

            // LEFT of the last arm in group_k
            let (_, seg_left, out_left) = group_k.last;
            let od_left = &offset_data[seg_left];
            let (left_base, left_dir) = offset_line_at_node(od_left, out_left, true);

            // RIGHT of the first arm in group_next
            let (_, seg_right, out_right) = group_next.first;
            let od_right = &offset_data[seg_right];
            let (right_base, right_dir) = offset_line_at_node(od_right, out_right, false);

            let p1 = Point3::new(left_base.0, left_base.1, 0.0);
            let d1 = Vector3::new(left_dir.0, left_dir.1, 0.0);
            let p2 = Point3::new(right_base.0, right_base.1, 0.0);
            let d2 = Vector3::new(right_dir.0, right_dir.1, 0.0);

            let corner = if let Some((t, _u)) = line_line_intersect_2d(&p1, &d1, &p2, &d2) {
                (p1.x + d1.x * t, p1.y + d1.y * t)
            } else {
                ((left_base.0 + right_base.0) * 0.5, (left_base.1 + right_base.1) * 0.5)
            };

            // Assign LEFT corner to ALL arms in group_k
            for &(_, seg_i, out_i) in &group_k.members {
                if out_i {
                    resolved[seg_i].left_start = corner;
                } else {
                    resolved[seg_i].right_end = corner;
                }
            }

            // Assign RIGHT corner to ALL arms in group_next
            for &(_, seg_j, out_j) in &group_next.members {
                if out_j {
                    resolved[seg_j].right_start = corner;
                } else {
                    resolved[seg_j].left_end = corner;
                }
            }
        }
    }

    resolved
}

struct ArmGroup {
    first: (f64, usize, bool),
    last: (f64, usize, bool),
    members: Vec<(f64, usize, bool)>,
}

fn merge_parallel_arms(
    arms: &[(f64, usize, bool)],
    offset_data: &[OffsetLineData],
) -> Vec<ArmGroup> {
    let n = arms.len();
    let mut groups: Vec<ArmGroup> = Vec::new();
    let mut i = 0;

    while i < n {
        let mut members = vec![arms[i]];
        let mut j = i + 1;

        while j < n {
            let dir_a = arm_direction(arms[j - 1], offset_data);
            let dir_b = arm_direction(arms[j], offset_data);
            let cross = dir_a.0 * dir_b.1 - dir_a.1 * dir_b.0;
            let dot = dir_a.0 * dir_b.0 + dir_a.1 * dir_b.1;
            if cross.abs() < PARALLEL_SIN && dot > 0.0 {
                members.push(arms[j]);
                j += 1;
            } else {
                break;
            }
        }

        groups.push(ArmGroup {
            first: members[0],
            last: *members.last().unwrap_or(&members[0]),
            members,
        });
        i = j;
    }

    // Check wrap-around: if the last group and first group are parallel,
    // merge them.
    if groups.len() >= 2 {
        let last_dir = arm_direction(groups.last().unwrap().last, offset_data);
        let first_dir = arm_direction(groups[0].first, offset_data);
        let cross = last_dir.0 * first_dir.1 - last_dir.1 * first_dir.0;
        let dot = last_dir.0 * first_dir.0 + last_dir.1 * first_dir.1;
        if cross.abs() < PARALLEL_SIN && dot > 0.0 {
            let mut last_group = groups.pop().unwrap();
            last_group.members.extend_from_slice(&groups[0].members);
            let new_first = groups[0].first;
            groups[0] = ArmGroup {
                first: last_group.first,
                last: new_first,
                members: last_group.members,
            };
        }
    }

    groups
}

fn arm_direction(arm: (f64, usize, bool), offset_data: &[OffsetLineData]) -> (f64, f64) {
    let (_, seg_idx, outgoing) = arm;
    let od = &offset_data[seg_idx];
    if outgoing { od.dir } else { neg_dir(od.dir) }
}

/// Returns the base point and direction for an offset line at a junction node.
///
/// For incoming arms (outgoing=false), left and right are swapped relative to
/// the original segment direction, because the arm's outgoing direction from
/// the junction is the reverse of the segment direction.
fn offset_line_at_node(
    od: &OffsetLineData,
    outgoing: bool,
    is_left: bool,
) -> ((f64, f64), (f64, f64)) {
    let dir = if outgoing { od.dir } else { neg_dir(od.dir) };

    let base = match (outgoing, is_left) {
        (true, true) => od.left_start,   // left of outgoing at start
        (true, false) => od.right_start, // right of outgoing at start
        (false, true) => od.right_end,   // left of reversed = right of original, at end
        (false, false) => od.left_end,   // right of reversed = left of original, at end
    };

    (base, dir)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::super::decompose::UniqueSegment;
    use super::super::junction;
    use super::*;

    #[test]
    fn double_cross_junction_corners() {
        let segments = vec![
            UniqueSegment {
                start: (3.0, 0.0),
                end: (3.0, 10.0),
            },
            UniqueSegment {
                start: (0.0, 7.0),
                end: (10.0, 7.0),
            },
            UniqueSegment {
                start: (7.0, 0.0),
                end: (7.0, 10.0),
            },
            UniqueSegment {
                start: (0.0, 3.0),
                end: (10.0, 3.0),
            },
        ];
        let net = junction::build_network(&segments);
        let edges = build(&net, 0.3, 0.3);
        // Should produce offset edges for all sub-segments + caps.
        assert!(!edges.is_empty(), "should produce offset edges");
    }

    #[test]
    fn single_segment_offset_edges() {
        let segments = vec![UniqueSegment {
            start: (0.0, 0.0),
            end: (5.0, 0.0),
        }];
        let net = junction::build_network(&segments);
        let edges = build(&net, 0.3, 0.3);
        // 1 sub-seg → 2 side edges + 2 cap edges = 4.
        assert_eq!(edges.len(), 4, "expected 4 edges, got {}", edges.len());
    }

    #[test]
    fn closed_square_no_caps() {
        // 4-segment closed square: all corners are junctions, no dead ends.
        let segments = vec![
            UniqueSegment {
                start: (0.0, 0.0),
                end: (10.0, 0.0),
            },
            UniqueSegment {
                start: (10.0, 0.0),
                end: (10.0, 10.0),
            },
            UniqueSegment {
                start: (10.0, 10.0),
                end: (0.0, 10.0),
            },
            UniqueSegment {
                start: (0.0, 10.0),
                end: (0.0, 0.0),
            },
        ];
        let net = junction::build_network(&segments);
        let edges = build(&net, 0.3, 0.3);
        // 4 sub-segments × 2 side edges = 8. No cap edges (no dead ends).
        assert_eq!(
            edges.len(),
            8,
            "expected 8 edges (no caps), got {}",
            edges.len()
        );
    }

    #[test]
    fn single_segment_rectangle_boundary() {
        let segments = vec![UniqueSegment {
            start: (0.0, 0.0),
            end: (5.0, 0.0),
        }];
        let net = junction::build_network(&segments);
        let edges = build(&net, 0.3, 0.3);

        // Verify edge coordinates for a horizontal segment at y=0, half_width=0.3.
        // Left offset at y=0.3, right offset at y=-0.3.
        // Left (backward): (5, 0.3) → (0, 0.3)
        // Right (forward): (0, -0.3) → (5, -0.3)
        // Start cap: (0, 0.3) → (0, -0.3)
        // End cap: (5, -0.3) → (5, 0.3)

        let has_left = edges.iter().any(|e| {
            (e.start.0 - 5.0).abs() < 0.01
                && (e.start.1 - 0.3).abs() < 0.01
                && (e.end.0).abs() < 0.01
                && (e.end.1 - 0.3).abs() < 0.01
        });
        assert!(has_left, "missing left backward edge");

        let has_right = edges.iter().any(|e| {
            (e.start.0).abs() < 0.01
                && (e.start.1 + 0.3).abs() < 0.01
                && (e.end.0 - 5.0).abs() < 0.01
                && (e.end.1 + 0.3).abs() < 0.01
        });
        assert!(has_right, "missing right forward edge");
    }
}
