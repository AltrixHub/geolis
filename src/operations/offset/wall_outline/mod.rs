mod decompose;
mod junction;
mod offset_edges;
mod trace;

use crate::error::{OperationError, Result};
use crate::geometry::pline::Pline;

/// Generates wall outlines from a centerline network.
///
/// Given a `Pline` representing wall centerlines (potentially with
/// self-intersecting paths), produces closed outline polygons at the
/// specified half-width distance.
#[derive(Debug)]
pub struct WallOutline2D {
    pline: Pline,
    half_width: f64,
}

impl WallOutline2D {
    /// Creates a new wall outline operation.
    #[must_use]
    pub fn new(pline: Pline, half_width: f64) -> Self {
        Self { pline, half_width }
    }

    /// Executes the wall outline generation.
    ///
    /// # Errors
    ///
    /// Returns `OperationError::InvalidInput` if the polyline has fewer than
    /// 2 vertices, or `OperationError::Failed` if no outline can be generated.
    pub fn execute(&self) -> Result<Vec<Pline>> {
        if self.pline.vertices.len() < 2 {
            return Err(OperationError::InvalidInput(
                "at least 2 vertices required for wall outline".to_owned(),
            )
            .into());
        }

        if self.half_width.abs() < crate::math::TOLERANCE {
            return Ok(vec![self.pline.clone()]);
        }

        // Step 1: Decompose into unique segments.
        let segments = decompose::decompose(&self.pline);
        if segments.is_empty() {
            return Err(
                OperationError::Failed("no valid segments in pline".to_owned()).into(),
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
        let wall = WallOutline2D::new(closed_square_pline(), 0.3);
        let result = wall.execute().unwrap();
        // Expect 2 closed boundaries: outer and inner.
        assert_eq!(result.len(), 2, "expected 2 boundaries, got {}", result.len());
        assert!(result.iter().all(|p| p.closed), "all boundaries should be closed");
    }

    #[test]
    fn closed_l_room_wall_outline() {
        let wall = WallOutline2D::new(closed_l_room_pline(), 0.3);
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
        let wall = WallOutline2D::new(closed_room_with_corridor_pline(), 0.3);
        let result = wall.execute().unwrap();
        // Expect 2 boundaries: outer (room+corridor) + inner (room).
        assert_eq!(result.len(), 2, "expected 2 boundaries, got {}", result.len());
        assert!(result.iter().all(|p| p.closed), "all boundaries should be closed");
    }

    #[test]
    fn closed_room_with_partition() {
        let wall = WallOutline2D::new(closed_room_with_partition_pline(), 0.3);
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
        let wall = WallOutline2D::new(closed_room_with_penetrating_wall_pline(), 0.3);
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
        let wall = WallOutline2D::new(closed_room_with_diagonal_wall_pline(), 0.3);
        let result = wall.execute().unwrap();
        // Expect 3 boundaries: outer (room+diagonal extensions) + 2 inner rooms.
        assert_eq!(result.len(), 3, "expected 3 boundaries, got {}", result.len());
        assert!(result.iter().all(|p| p.closed), "all boundaries should be closed");
    }

    #[test]
    fn debug_double_cross_wall_outline() {
        let wall = WallOutline2D::new(double_cross_pline(), 0.3);
        let result = wall.execute().unwrap();
        assert!(!result.is_empty(), "expected at least 1 boundary");
    }
}
