//! 2D polygon boolean-union via the shared arrangement engine.
//!
//! [`union_all_with_holes`] returns the boolean-union outline of one or
//! more [`PolygonWithHoles`] inputs as typed face topology. Output
//! edges are exactly the boundary that separates **filled** from
//! **empty** points, where "filled" follows the OR-of-PWH-filled rule
//! (a point is filled iff there exists at least one input PWH whose
//! outer contains the point AND none of whose holes contains it).
//!
//! All arrangement work (segment split, vertex snap, half-edge
//! classification, face walk, face assembly) lives in
//! [`super::engine`]; this module supplies the union-specific fill
//! oracle and the thin user-facing entry point.

use crate::error::Result;

use super::engine::{run_arrangement, run_arrangement_traced, TracedFace, UnionOracle};
use super::types::{PolygonWithHoles, UnionResult};

/// Compute the boolean-union outline of `inputs`. Output boundary loops
/// are closed implicitly (vertex list `[v0, v1, ..., vn-1]` represents
/// the closed loop `v0 → v1 → ... → vn-1 → v0`).
///
/// Determinism: outputs are topologically identical (and float-equivalent
/// within `WALL_EPS` precision) regardless of input order.
///
/// # Errors
///
/// [`crate::error::OperationError::Failed`] propagated from the
/// engine when:
/// - A half-edge's bilateral classification remains
///   `AmbiguousOnBoundary` after 3 ε-shrink retries (typically
///   indicates degenerate input where multiple inputs share a tangent
///   boundary at the sampled edge midpoint).
/// - The face-assembly stage cannot pick a unique parent for a nested
///   loop, or detects an orientation/depth parity violation.
pub fn union_all_with_holes(inputs: &[PolygonWithHoles]) -> Result<UnionResult> {
    let oracle = UnionOracle { inputs };
    let faces = run_arrangement(inputs, &oracle)?;
    Ok(UnionResult { faces })
}

/// [`union_all_with_holes`] variant that additionally reports, per output
/// edge, the [`super::engine::SegmentSite`] of the input edge it came
/// from. Sites are threaded through the arrangement (never recovered by
/// geometric matching), so they are exact and deterministic.
///
/// # Errors
///
/// Same failure modes as [`union_all_with_holes`].
pub(crate) fn union_all_with_holes_traced(inputs: &[PolygonWithHoles]) -> Result<Vec<TracedFace>> {
    let oracle = UnionOracle { inputs };
    run_arrangement_traced(inputs, &oracle)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::super::types::{signed_area, Polygon, WALL_EPS};
    use super::*;

    fn rect(x: f64, y: f64, w: f64, h: f64) -> Polygon {
        vec![(x, y), (x + w, y), (x + w, y + h), (x, y + h)]
    }

    /// Flat-boundary view of a union result, used by tests written
    /// against the pre-`UnionResult` API. Equivalent to concatenating
    /// every face's outer + holes.
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
        let hole = vec![(2.0, 2.0), (2.0, 8.0), (8.0, 8.0), (8.0, 2.0)];
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

    // ===== union_all_with_holes_traced — per-edge source sites =====

    /// Every traced output edge must lie on the supporting line of the
    /// input edge its site names (within the snap tolerance), and the
    /// site arrays must align 1:1 with the ring vertex counts.
    fn assert_sites_exact(traced: &[super::super::TracedFace], inputs: &[PolygonWithHoles]) {
        use super::super::engine::RingRef;
        for tf in traced {
            let rings: Vec<(&Polygon, &Vec<super::super::SegmentSite>)> =
                std::iter::once((&tf.face.outer, &tf.outer_sites))
                    .chain(tf.face.holes.iter().zip(tf.hole_sites.iter()))
                    .collect();
            for (pts, sites) in rings {
                assert_eq!(pts.len(), sites.len(), "site array must align with ring");
                for (e, site) in sites.iter().enumerate() {
                    let src_ring: &Polygon = match site.ring {
                        RingRef::Outer => &inputs[site.input].outer,
                        RingRef::Hole(h) => &inputs[site.input].holes[h],
                    };
                    let a = src_ring[site.edge];
                    let b = src_ring[(site.edge + 1) % src_ring.len()];
                    let dx = b.0 - a.0;
                    let dy = b.1 - a.1;
                    let len = (dx * dx + dy * dy).sqrt();
                    for p in [pts[e], pts[(e + 1) % pts.len()]] {
                        let perp = ((p.0 - a.0) * dy - (p.1 - a.1) * dx) / len;
                        assert!(
                            perp.abs() < WALL_EPS * 5.0,
                            "output edge {e} not on supporting line of its \
                             site {site:?}: perp={perp}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn traced_t_shape_sites_are_exact() {
        let inputs = no_hole_inputs(vec![rect(0.0, -1.0, 8.0, 2.0), rect(3.0, -1.0, 2.0, 5.0)]);
        let traced = union_all_with_holes_traced(&inputs).expect("union must succeed");
        assert_eq!(traced.len(), 1);
        assert_sites_exact(&traced, &inputs);
        // Both inputs must be represented on the merged boundary.
        let mut input_ids: Vec<usize> = traced[0].outer_sites.iter().map(|s| s.input).collect();
        input_ids.sort_unstable();
        input_ids.dedup();
        assert_eq!(input_ids, vec![0, 1]);
    }

    #[test]
    fn traced_identical_inputs_attribute_to_smallest_site() {
        // Two bit-identical squares: the dedup tie-break must attribute
        // every surviving edge to input 0, deterministically.
        let inputs = no_hole_inputs(vec![rect(0.0, 0.0, 4.0, 3.0), rect(0.0, 0.0, 4.0, 3.0)]);
        let traced = union_all_with_holes_traced(&inputs).expect("union must succeed");
        assert_eq!(traced.len(), 1);
        assert_sites_exact(&traced, &inputs);
        for site in &traced[0].outer_sites {
            assert_eq!(site.input, 0, "coincident edges must resolve to input 0");
        }
    }

    #[test]
    fn traced_donut_hole_ring_carries_sites() {
        use super::super::engine::RingRef;
        let outer = rect(0.0, 0.0, 10.0, 10.0);
        let hole = vec![(2.0, 2.0), (2.0, 8.0), (8.0, 8.0), (8.0, 2.0)];
        let inputs = vec![PolygonWithHoles {
            outer,
            holes: vec![hole],
        }];
        let traced = union_all_with_holes_traced(&inputs).expect("union must succeed");
        assert_eq!(traced.len(), 1);
        assert_eq!(traced[0].face.holes.len(), 1);
        assert_sites_exact(&traced, &inputs);
        for site in &traced[0].hole_sites[0] {
            assert_eq!(site.ring, RingRef::Hole(0));
        }
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
}
