mod decompose;
mod junction;
mod offset_edges;
mod trace;

use crate::error::{OperationError, Result};
use crate::geometry::pline::Pline;

/// Generates wall outlines from one or more centerline polylines.
///
/// Given a collection of `Pline`s representing wall centerlines (potentially with
/// self-intersecting paths), produces closed outline polygons at the
/// specified half-width distance. When multiple polylines are provided,
/// their segments are merged into a single network so that intersections
/// between separate walls are properly trimmed.
#[derive(Debug)]
pub struct WallOutline2D {
    plines: Vec<Pline>,
    half_width: f64,
}

impl WallOutline2D {
    /// Creates a new wall outline operation from one or more centerline polylines.
    #[must_use]
    pub fn new(plines: Vec<Pline>, half_width: f64) -> Self {
        Self { plines, half_width }
    }

    /// Executes the wall outline generation.
    ///
    /// # Errors
    ///
    /// Returns `OperationError::InvalidInput` if no polyline has at least
    /// 2 vertices, or `OperationError::Failed` if no outline can be generated.
    pub fn execute(&self) -> Result<Vec<Pline>> {
        let valid: Vec<&Pline> = self.plines.iter()
            .filter(|p| p.vertices.len() >= 2)
            .collect();

        if valid.is_empty() {
            return Err(OperationError::InvalidInput(
                "at least 2 vertices required for wall outline".to_owned(),
            )
            .into());
        }

        if self.half_width.abs() < crate::math::TOLERANCE {
            return Ok(self.plines.clone());
        }

        // Step 1: Decompose all polylines into unique segments.
        let segments = decompose::decompose(&valid);
        if segments.is_empty() {
            return Err(
                OperationError::Failed("no valid segments in plines".to_owned()).into(),
            );
        }

        // Step 2: Detect junctions and split segments.
        let network = junction::build_network(&segments);

        // Step 3: Generate offset edges with junction resolution.
        let edges = offset_edges::build(&network, self.half_width);

        // Step 4: Trace outer boundaries.
        let outlines = trace::trace_boundaries(&edges);

        if outlines.is_empty() {
            return Err(
                OperationError::Failed("wall outline trace produced no results".to_owned())
                    .into(),
            );
        }

        Ok(outlines)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::pline::PlineVertex;

    fn double_cross_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(3.0, 0.0),
                PlineVertex::line(3.0, 10.0),
                PlineVertex::line(3.0, 7.0),
                PlineVertex::line(0.0, 7.0),
                PlineVertex::line(10.0, 7.0),
                PlineVertex::line(7.0, 7.0),
                PlineVertex::line(7.0, 10.0),
                PlineVertex::line(7.0, 0.0),
                PlineVertex::line(7.0, 3.0),
                PlineVertex::line(10.0, 3.0),
                PlineVertex::line(0.0, 3.0),
            ],
            closed: false,
        }
    }

    fn closed_square_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(10.0, 10.0),
                PlineVertex::line(0.0, 10.0),
            ],
            closed: true,
        }
    }

    fn closed_l_room_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(5.0, 3.0),
                PlineVertex::line(3.0, 3.0),
                PlineVertex::line(3.0, 5.0),
                PlineVertex::line(0.0, 5.0),
            ],
            closed: true,
        }
    }

    #[test]
    fn closed_square_wall_outline() {
        let wall = WallOutline2D::new(vec![closed_square_pline()], 0.3);
        let result = wall.execute().unwrap();
        // Expect 2 closed boundaries: outer and inner.
        assert_eq!(result.len(), 2, "expected 2 boundaries, got {}", result.len());
        assert!(result.iter().all(|p| p.closed), "all boundaries should be closed");
    }

    #[test]
    fn closed_l_room_wall_outline() {
        let wall = WallOutline2D::new(vec![closed_l_room_pline()], 0.3);
        let result = wall.execute().unwrap();
        // Expect 2 closed boundaries: outer and inner.
        assert_eq!(result.len(), 2, "expected 2 boundaries, got {}", result.len());
        assert!(result.iter().all(|p| p.closed), "all boundaries should be closed");
    }

    /// Closed square room + corridor extending outward from bottom.
    fn closed_room_with_corridor_pline() -> Pline {
        // Room: (0,0)-(10,0)-(10,10)-(0,10), corridor at (5,0) going down to (5,-5).
        Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(5.0, -5.0),
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(10.0, 10.0),
                PlineVertex::line(0.0, 10.0),
            ],
            closed: true,
        }
    }

    /// Closed square room divided by a horizontal partition at y=5.
    fn closed_room_with_partition_pline() -> Pline {
        // Room: (0,0)-(10,0)-(10,10)-(0,10), partition at y=5.
        Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(10.0, 5.0),
                PlineVertex::line(0.0, 5.0),
                PlineVertex::line(10.0, 5.0),
                PlineVertex::line(10.0, 10.0),
                PlineVertex::line(0.0, 10.0),
            ],
            closed: true,
        }
    }

    #[test]
    fn closed_room_with_corridor() {
        let wall = WallOutline2D::new(vec![closed_room_with_corridor_pline()], 0.3);
        let result = wall.execute().unwrap();
        // Expect 2 boundaries: outer (room+corridor) + inner (room).
        assert_eq!(result.len(), 2, "expected 2 boundaries, got {}", result.len());
        assert!(result.iter().all(|p| p.closed), "all boundaries should be closed");
    }

    #[test]
    fn closed_room_with_partition() {
        let wall = WallOutline2D::new(vec![closed_room_with_partition_pline()], 0.3);
        let result = wall.execute().unwrap();
        // Expect 3 boundaries: outer + 2 inner rooms.
        assert_eq!(result.len(), 3, "expected 3 boundaries, got {}", result.len());
        assert!(result.iter().all(|p| p.closed), "all boundaries should be closed");
    }

    /// Closed room with a wall penetrating through both sides at y=5.
    /// Wall extends from x=-3 to x=13, passing through the room (0,0)-(10,10).
    ///
    /// The path must explicitly traverse the interior partition (0,5)→(10,5)
    /// so the decompose step produces a continuous segment from (-3,5) to (13,5).
    fn closed_room_with_penetrating_wall_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(10.0, 5.0),
                PlineVertex::line(13.0, 5.0),
                PlineVertex::line(10.0, 5.0),
                PlineVertex::line(0.0, 5.0),  // interior partition: (10,5)→(0,5)
                PlineVertex::line(10.0, 5.0),  // backtrack: (0,5)→(10,5)
                PlineVertex::line(10.0, 10.0),
                PlineVertex::line(0.0, 10.0),
                PlineVertex::line(0.0, 5.0),
                PlineVertex::line(-3.0, 5.0),
                PlineVertex::line(0.0, 5.0),
            ],
            closed: true,
        }
    }

    #[test]
    fn closed_room_with_penetrating_wall() {
        let wall = WallOutline2D::new(vec![closed_room_with_penetrating_wall_pline()], 0.3);
        let result = wall.execute().unwrap();
        // Expect 3 boundaries: outer (room+extensions) + 2 inner (top/bottom rooms).
        assert_eq!(result.len(), 3, "expected 3 boundaries, got {}", result.len());
        assert!(result.iter().all(|p| p.closed), "all boundaries should be closed");
    }

    /// Closed room with a diagonal wall penetrating through both sides.
    /// Diagonal from (-5,0) to (15,10) passes through the room at (0,2.5) and (10,7.5).
    /// The diagonal is encoded as a single line (15,10)→(-5,0).
    fn closed_room_with_diagonal_wall_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(10.0, 7.5),
                PlineVertex::line(15.0, 10.0),  // right extension
                PlineVertex::line(-5.0, 0.0),   // full diagonal in one line
                PlineVertex::line(0.0, 2.5),
                PlineVertex::line(0.0, 10.0),
                PlineVertex::line(10.0, 10.0),
                PlineVertex::line(10.0, 7.5),
                PlineVertex::line(0.0, 2.5),
            ],
            closed: true,
        }
    }

    #[test]
    fn closed_room_with_diagonal_wall() {
        let wall = WallOutline2D::new(vec![closed_room_with_diagonal_wall_pline()], 0.3);
        let result = wall.execute().unwrap();
        // Expect 3 boundaries: outer (room+diagonal extensions) + 2 inner rooms.
        assert_eq!(result.len(), 3, "expected 3 boundaries, got {}", result.len());
        assert!(result.iter().all(|p| p.closed), "all boundaries should be closed");
    }

    #[test]
    fn debug_double_cross_wall_outline() {
        let wall = WallOutline2D::new(vec![double_cross_pline()], 0.3);
        let result = wall.execute().unwrap();
        assert!(!result.is_empty(), "expected at least 1 boundary");
    }

    // ── Spoke / bilateral buffer tests (migrated from PlineOffset2D) ──

    /// Checks that every expected point appears somewhere in the result vertices.
    fn assert_vertices_match(result: &[PlineVertex], expected: &[(f64, f64)], tol: f64) {
        assert_eq!(
            result.len(),
            expected.len(),
            "vertex count mismatch: got {}, expected {}",
            result.len(),
            expected.len()
        );
        for &(ex, ey) in expected {
            let found = result
                .iter()
                .any(|v| (v.x - ex).abs() < tol && (v.y - ey).abs() < tol);
            assert!(
                found,
                "expected vertex ({ex:.4}, {ey:.4}) not found in result"
            );
        }
    }

    /// Open cross: 4 arms from center with 180° reversals.
    fn open_cross_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(-1.5, 0.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(0.0, 1.5),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(1.5, 0.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(0.0, -1.5),
            ],
            closed: false,
        }
    }

    /// Expected cross outline at distance d (12 vertices).
    fn open_cross_expected(d: f64) -> Vec<(f64, f64)> {
        vec![
            (-1.5, -d),
            (-1.5, d),
            (-d, d),
            (-d, 1.5),
            (d, 1.5),
            (d, d),
            (1.5, d),
            (1.5, -d),
            (d, -d),
            (d, -1.5),
            (-d, -1.5),
            (-d, -d),
        ]
    }

    #[test]
    fn open_cross_d03() {
        let wall = WallOutline2D::new(vec![open_cross_pline()], 0.3);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 closed polygon");
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        assert_vertices_match(&poly.vertices, &open_cross_expected(0.3), 0.05);
    }

    #[test]
    fn open_cross_d05() {
        let wall = WallOutline2D::new(vec![open_cross_pline()], 0.5);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 closed polygon");
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        assert_vertices_match(&poly.vertices, &open_cross_expected(0.5), 0.05);
    }

    /// X-cross: 2 diagonal lines crossing at center with 180° reversals.
    fn x_cross_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(-3.0, -3.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(3.0, 3.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(-3.0, 3.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(3.0, -3.0),
            ],
            closed: false,
        }
    }

    /// Expected X-cross outline at distance d (12 vertices).
    fn x_cross_expected(a: f64, d: f64) -> Vec<(f64, f64)> {
        let s2 = std::f64::consts::SQRT_2;
        let h = d * s2 / 2.0;
        let d2 = d * s2;
        vec![
            (-a - h, -a + h),
            (-a + h, -a - h),
            (0.0, -d2),
            (a - h, -a - h),
            (a + h, -a + h),
            (d2, 0.0),
            (a + h, a - h),
            (a - h, a + h),
            (0.0, d2),
            (-a + h, a + h),
            (-a - h, a - h),
            (-d2, 0.0),
        ]
    }

    #[test]
    fn x_cross_d05() {
        let wall = WallOutline2D::new(vec![x_cross_pline()], 0.5);
        let result = wall.execute().unwrap();
        let expected = x_cross_expected(3.0, 0.5);
        assert_eq!(result.len(), 1, "expected 1 closed polygon");
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        assert_vertices_match(&poly.vertices, &expected, 0.05);
    }

    /// Fork (Y-shape): stem + 2 branches with reversal at junction.
    fn fork_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(5.0, 4.0),
                PlineVertex::line(0.0, 9.0),
                PlineVertex::line(5.0, 4.0),
                PlineVertex::line(10.0, 9.0),
            ],
            closed: false,
        }
    }

    /// Expected fork outline at distance d (9 vertices).
    fn fork_expected(d: f64) -> Vec<(f64, f64)> {
        let s2 = std::f64::consts::SQRT_2;
        let h = d * s2 / 2.0;
        let jy = 4.0 + d * (1.0 - s2);
        vec![
            (5.0 - d, 0.0),
            (5.0 + d, 0.0),
            (5.0 + d, jy),
            (10.0 + h, 9.0 - h),
            (10.0 - h, 9.0 + h),
            (5.0, 4.0 + d * s2),
            (h, 9.0 + h),
            (-h, 9.0 - h),
            (5.0 - d, jy),
        ]
    }

    #[test]
    fn fork_d05() {
        let wall = WallOutline2D::new(vec![fork_pline()], 0.5);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 closed polygon");
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        assert_vertices_match(&poly.vertices, &fork_expected(0.5), 0.05);
    }

    /// Expected double-cross outline at distance d (28 vertices).
    fn double_cross_expected(d: f64) -> Vec<(f64, f64)> {
        vec![
            (3.0 - d, 0.0),
            (3.0 + d, 0.0),
            (3.0 + d, 3.0 - d),
            (7.0 - d, 3.0 - d),
            (7.0 - d, 0.0),
            (7.0 + d, 0.0),
            (7.0 + d, 3.0 - d),
            (10.0, 3.0 - d),
            (10.0, 3.0 + d),
            (7.0 + d, 3.0 + d),
            (7.0 + d, 7.0 - d),
            (10.0, 7.0 - d),
            (10.0, 7.0 + d),
            (7.0 + d, 7.0 + d),
            (7.0 + d, 10.0),
            (7.0 - d, 10.0),
            (7.0 - d, 7.0 + d),
            (3.0 + d, 7.0 + d),
            (3.0 + d, 10.0),
            (3.0 - d, 10.0),
            (3.0 - d, 7.0 + d),
            (0.0, 7.0 + d),
            (0.0, 7.0 - d),
            (3.0 - d, 7.0 - d),
            (3.0 - d, 3.0 + d),
            (0.0, 3.0 + d),
            (0.0, 3.0 - d),
            (3.0 - d, 3.0 - d),
        ]
    }

    #[test]
    fn double_cross_d03() {
        let wall = WallOutline2D::new(vec![double_cross_pline()], 0.3);
        let result = wall.execute().unwrap();
        assert!(!result.is_empty(), "expected at least 1 polygon");
        let outer = result
            .iter()
            .max_by_key(|p| p.vertices.len())
            .unwrap();
        assert!(outer.closed, "outer boundary should be closed");
        assert_vertices_match(&outer.vertices, &double_cross_expected(0.3), 0.05);
    }

    #[test]
    fn double_cross_d08() {
        let wall = WallOutline2D::new(vec![double_cross_pline()], 0.8);
        let result = wall.execute().unwrap();
        assert!(!result.is_empty(), "expected at least 1 polygon");
        let outer = result
            .iter()
            .max_by_key(|p| p.vertices.len())
            .unwrap();
        assert!(outer.closed, "outer boundary should be closed");
        assert_vertices_match(&outer.vertices, &double_cross_expected(0.8), 0.05);
    }

    // ── Vertex-intersection stress tests (cases 18–22) ──

    /// Case 18: T-shape, arm shorter than d (concave notch).
    /// Arm length 0.5, d=1.0 — arm cap below spine wall, concave notch.
    #[test]
    fn t_very_short_arm_d10() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 3.0),
                PlineVertex::line(8.0, 3.0),
                PlineVertex::line(4.0, 3.0),
                PlineVertex::line(4.0, 3.5),
            ],
            closed: false,
        };
        let wall = WallOutline2D::new(vec![pline], 1.0);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 boundary, got {}", result.len());
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        let expected = vec![
            (0.0, 2.0), (4.0, 2.0), (8.0, 2.0), (8.0, 4.0),
            (5.0, 4.0), (5.0, 3.5), (3.0, 3.5), (3.0, 4.0), (0.0, 4.0),
        ];
        assert_vertices_match(&poly.vertices, &expected, 0.05);
    }

    /// Case 19: T-shape, arm = 2d (baseline comparison).
    /// Arm length 2, d=1.0 — clean T-shape outline.
    #[test]
    fn t_arm_2d_d10() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 3.0),
                PlineVertex::line(8.0, 3.0),
                PlineVertex::line(4.0, 3.0),
                PlineVertex::line(4.0, 5.0),
            ],
            closed: false,
        };
        let wall = WallOutline2D::new(vec![pline], 1.0);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 boundary, got {}", result.len());
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        let expected = vec![
            (0.0, 2.0), (4.0, 2.0), (8.0, 2.0), (8.0, 4.0),
            (5.0, 4.0), (5.0, 5.0), (3.0, 5.0), (3.0, 4.0), (0.0, 4.0),
        ];
        assert_vertices_match(&poly.vertices, &expected, 0.05);
    }

    /// Case 20: Cross, arm length = d — degenerates to square.
    /// Arm length 2, d=2.0 — all arm side edges degenerate.
    #[test]
    fn cross_short_d20() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(4.0, 2.0),
                PlineVertex::line(4.0, 6.0),
                PlineVertex::line(4.0, 4.0),
                PlineVertex::line(2.0, 4.0),
                PlineVertex::line(6.0, 4.0),
            ],
            closed: false,
        };
        let wall = WallOutline2D::new(vec![pline], 2.0);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 boundary, got {}", result.len());
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        let expected = vec![(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)];
        assert_vertices_match(&poly.vertices, &expected, 0.05);
    }

    /// Case 21: L-shape, d > horizontal leg → miter extends past original.
    /// horizontal leg=2, d=2.5 → miter pushes far past original vertex.
    #[test]
    fn l_large_d_d25() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(2.0, 0.0),
                PlineVertex::line(2.0, 4.0),
            ],
            closed: false,
        };
        let wall = WallOutline2D::new(vec![pline], 2.5);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 boundary, got {}", result.len());
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        let expected = vec![
            (0.0, -2.5), (4.5, -2.5), (4.5, 4.0),
            (-0.5, 4.0), (-0.5, 2.5), (0.0, 2.5),
        ];
        assert_vertices_match(&poly.vertices, &expected, 0.05);
    }

    /// Case 22: T-shape, arm = d (exact degeneration boundary).
    /// Arm length 1, d=1.0 — arm side edges degenerate to zero length.
    #[test]
    fn t_arm_eq_d_d10() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 3.0),
                PlineVertex::line(8.0, 3.0),
                PlineVertex::line(4.0, 3.0),
                PlineVertex::line(4.0, 4.0),
            ],
            closed: false,
        };
        let wall = WallOutline2D::new(vec![pline], 1.0);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 boundary, got {}", result.len());
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        let expected = vec![
            (0.0, 2.0), (4.0, 2.0), (8.0, 2.0), (8.0, 4.0),
            (5.0, 4.0), (3.0, 4.0), (0.0, 4.0),
        ];
        assert_vertices_match(&poly.vertices, &expected, 0.05);
    }

    // ── Diagonal / non-orthogonal junction tests (cases 23–27) ──

    /// Case 23: Open L-shape at 45°, d=0.3 — non-orthogonal miter.
    #[test]
    fn l_shape_45_d03() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(8.0, 3.0),
            ],
            closed: false,
        };
        let d = 0.3;
        let ds = d / std::f64::consts::SQRT_2;
        let dm = d * (std::f64::consts::SQRT_2 - 1.0);
        let wall = WallOutline2D::new(vec![pline], d);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 boundary, got {}", result.len());
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        let expected = vec![
            (0.0, -d), (5.0 + dm, -d), (8.0 + ds, 3.0 - ds),
            (8.0 - ds, 3.0 + ds), (5.0 - dm, d), (0.0, d),
        ];
        assert_vertices_match(&poly.vertices, &expected, 0.05);
    }

    /// Case 24: Open T with 45° upward branch, d=0.3.
    #[test]
    fn t_diagonal_branch_d03() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 3.0),
                PlineVertex::line(10.0, 3.0),
                PlineVertex::line(5.0, 3.0),
                PlineVertex::line(7.0, 5.0),
            ],
            closed: false,
        };
        let d = 0.3;
        let ds = d / std::f64::consts::SQRT_2;
        let s2 = std::f64::consts::SQRT_2;
        let wall = WallOutline2D::new(vec![pline], d);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 boundary, got {}", result.len());
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        let expected = vec![
            (0.0, 3.0 - d), (5.0, 3.0 - d), (10.0, 3.0 - d),
            (10.0, 3.0 + d), (5.0 + d * (1.0 + s2), 3.0 + d),
            (7.0 + ds, 5.0 - ds), (7.0 - ds, 5.0 + ds),
            (5.0 + d * (1.0 - s2), 3.0 + d), (0.0, 3.0 + d),
        ];
        assert_vertices_match(&poly.vertices, &expected, 0.05);
    }

    /// Case 25: Open Y-junction (L-shape + diagonal), d=0.3.
    #[test]
    fn y_mixed_junction_d03() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(5.0, 5.0),
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(8.0, 3.0),
            ],
            closed: false,
        };
        let d = 0.3;
        let ds = d / std::f64::consts::SQRT_2;
        let dm = d * (std::f64::consts::SQRT_2 - 1.0);
        let s2 = std::f64::consts::SQRT_2;
        let wall = WallOutline2D::new(vec![pline], d);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 boundary, got {}", result.len());
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        let expected = vec![
            (0.0, -d), (5.0 + dm, -d),
            (8.0 + ds, 3.0 - ds), (8.0 - ds, 3.0 + ds),
            (5.0 + d, d * (1.0 + s2)),
            (5.0 + d, 5.0), (5.0 - d, 5.0),
            (5.0 - d, d), (0.0, d),
        ];
        assert_vertices_match(&poly.vertices, &expected, 0.05);
    }

    /// Case 26: Closed room + diagonal stub from corner, d=0.3.
    #[test]
    fn room_corner_stub_d03() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(-3.0, -3.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(8.0, 0.0),
                PlineVertex::line(8.0, 8.0),
                PlineVertex::line(0.0, 8.0),
            ],
            closed: true,
        };
        let d = 0.3;
        let ds = d / std::f64::consts::SQRT_2;
        let dm = d * (std::f64::consts::SQRT_2 - 1.0);
        let wall = WallOutline2D::new(vec![pline], d);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 2, "expected 2 boundaries, got {}", result.len());

        let outer = result.iter().max_by_key(|p| p.vertices.len()).unwrap();
        let outer_expected = vec![
            (8.0 + d, -d), (8.0 + d, 8.0 + d), (-d, 8.0 + d),
            (-d, dm),
            (-3.0 - ds, -3.0 + ds), (-3.0 + ds, -3.0 - ds),
            (dm, -d),
        ];
        assert_vertices_match(&outer.vertices, &outer_expected, 0.05);

        let inner = result.iter().min_by_key(|p| p.vertices.len()).unwrap();
        let inner_expected = vec![
            (d, d), (8.0 - d, d), (8.0 - d, 8.0 - d), (d, 8.0 - d),
        ];
        assert_vertices_match(&inner.vertices, &inner_expected, 0.05);
    }

    /// Case 27: Closed room + diagonal partition through corner, d=0.3.
    #[test]
    fn room_corner_diagonal_d03() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(10.0, 8.0),
                PlineVertex::line(8.0, 8.0),
                PlineVertex::line(11.0, 11.0),
                PlineVertex::line(-3.0, -3.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(0.0, 8.0),
                PlineVertex::line(8.0, 8.0),
            ],
            closed: true,
        };
        let d = 0.3;
        let ds = d / std::f64::consts::SQRT_2;
        let dm = d * (std::f64::consts::SQRT_2 - 1.0);
        let s2 = std::f64::consts::SQRT_2;
        let wall = WallOutline2D::new(vec![pline], d);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 3, "expected 3 boundaries, got {}", result.len());

        let outer = result.iter().max_by_key(|p| p.vertices.len()).unwrap();
        let outer_expected = vec![
            (10.0 + d, -d), (10.0 + d, 8.0 + d),
            (8.0 + d * (1.0 + s2), 8.0 + d),
            (11.0 + ds, 11.0 - ds), (11.0 - ds, 11.0 + ds),
            (8.0 + d * (1.0 - s2), 8.0 + d),
            (-d, 8.0 + d), (-d, dm),
            (-3.0 - ds, -3.0 + ds), (-3.0 + ds, -3.0 - ds),
            (dm, -d),
        ];
        assert_vertices_match(&outer.vertices, &outer_expected, 0.05);
    }

    /// Case 28: Closed room + diagonal stub near corner, junction at (0, 0.5).
    #[test]
    fn room_near_corner_stub_d03() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.5),
                PlineVertex::line(-3.0, -2.5),
                PlineVertex::line(0.0, 0.5),
                PlineVertex::line(0.0, 8.0),
                PlineVertex::line(8.0, 8.0),
                PlineVertex::line(8.0, 0.0),
                PlineVertex::line(0.0, 0.0),
            ],
            closed: true,
        };
        let d = 0.3;
        let ds = d / std::f64::consts::SQRT_2;
        let dm = d * (std::f64::consts::SQRT_2 - 1.0);
        let s2 = std::f64::consts::SQRT_2;
        let wall = WallOutline2D::new(vec![pline], d);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 2, "expected 2 boundaries, got {}", result.len());

        let outer = result.iter().max_by_key(|p| p.vertices.len()).unwrap();
        let outer_expected = vec![
            (8.0 + d, -d), (8.0 + d, 8.0 + d), (-d, 8.0 + d),
            (-d, 0.5 + dm),
            (-3.0 - ds, -2.5 + ds), (-3.0 + ds, -2.5 - ds),
            (-d, 0.5 - d * (1.0 + s2)),
            (-d, -d),
        ];
        assert_vertices_match(&outer.vertices, &outer_expected, 0.05);

        // Inner boundary includes a collinear junction corner at (d, 0.5).
        let inner = result.iter().min_by_key(|p| p.vertices.len()).unwrap();
        let inner_expected = vec![
            (d, d), (d, 0.5), (d, 8.0 - d), (8.0 - d, 8.0 - d), (8.0 - d, d),
        ];
        assert_vertices_match(&inner.vertices, &inner_expected, 0.05);
    }

    /// Case 29: Closed room + diagonal partition near corner, junctions at (0, 0.5) and (7.5, 8).
    #[test]
    fn room_near_corner_diagonal_d03() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(10.0, 8.0),
                PlineVertex::line(7.5, 8.0),
                PlineVertex::line(10.5, 11.0),
                PlineVertex::line(-3.0, -2.5),
                PlineVertex::line(0.0, 0.5),
                PlineVertex::line(0.0, 8.0),
                PlineVertex::line(7.5, 8.0),
                PlineVertex::line(0.0, 0.5),
            ],
            closed: true,
        };
        let d = 0.3;
        let ds = d / std::f64::consts::SQRT_2;
        let dm = d * (std::f64::consts::SQRT_2 - 1.0);
        let s2 = std::f64::consts::SQRT_2;
        let wall = WallOutline2D::new(vec![pline], d);
        let result = wall.execute().unwrap();
        assert_eq!(result.len(), 3, "expected 3 boundaries, got {}", result.len());

        let outer = result.iter().max_by_key(|p| p.vertices.len()).unwrap();
        let outer_expected = vec![
            (10.0 + d, -d), (10.0 + d, 8.0 + d),
            (7.5 + d * (1.0 + s2), 8.0 + d),
            (10.5 + ds, 11.0 - ds), (10.5 - ds, 11.0 + ds),
            (7.5 + d * (1.0 - s2), 8.0 + d),
            (-d, 8.0 + d), (-d, 0.5 + dm),
            (-3.0 - ds, -2.5 + ds), (-3.0 + ds, -2.5 - ds),
            (-d, 0.5 - d * (1.0 + s2)),
            (-d, -d),
        ];
        assert_vertices_match(&outer.vertices, &outer_expected, 0.05);
    }
}
