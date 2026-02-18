use crate::error::{OperationError, Result};
use crate::geometry::curve::Line;
use crate::math::{Point3, TOLERANCE};
use crate::topology::{EdgeCurve, EdgeData, OrientedEdge, TopologyStore, VertexData, WireData, WireId};

/// Creates a wire from a sequence of 3D points.
pub struct MakeWire {
    points: Vec<Point3>,
    close: bool,
}

impl MakeWire {
    /// Creates a new `MakeWire` operation.
    #[must_use]
    pub fn new(points: Vec<Point3>, close: bool) -> Self {
        Self { points, close }
    }

    /// Executes the operation, creating the wire in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if fewer than 2 points are provided, or if consecutive
    /// points are coincident (distance < `TOLERANCE`).
    pub fn execute(&self, store: &mut TopologyStore) -> Result<WireId> {
        let n = self.points.len();
        if n < 2 {
            return Err(OperationError::InvalidInput(
                "at least 2 points are required to create a wire".into(),
            )
            .into());
        }

        // Validate no consecutive duplicates
        for i in 0..n - 1 {
            let dist = (self.points[i + 1] - self.points[i]).norm();
            if dist < TOLERANCE {
                return Err(OperationError::InvalidInput(format!(
                    "consecutive points {i} and {} are coincident",
                    i + 1
                ))
                .into());
            }
        }

        // Validate closing edge if closed
        if self.close {
            let dist = (self.points[0] - self.points[n - 1]).norm();
            if dist < TOLERANCE {
                return Err(OperationError::InvalidInput(
                    "last and first points are coincident for a closed wire".into(),
                )
                .into());
            }
        }

        // Create vertices
        let vertex_ids: Vec<_> = self
            .points
            .iter()
            .map(|p| store.add_vertex(VertexData::new(*p)))
            .collect();

        // Create edges between consecutive vertices
        let edge_count = if self.close { n } else { n - 1 };
        let mut oriented_edges = Vec::with_capacity(edge_count);

        for i in 0..edge_count {
            let start_idx = i;
            let end_idx = (i + 1) % n;
            let start_v = vertex_ids[start_idx];
            let end_v = vertex_ids[end_idx];
            let a = self.points[start_idx];
            let b = self.points[end_idx];
            let direction = b - a;
            let t_end = direction.norm();
            let line = Line::new(a, direction)?;

            let edge_id = store.add_edge(EdgeData {
                start: start_v,
                end: end_v,
                curve: EdgeCurve::Line(line),
                t_start: 0.0,
                t_end,
            });

            oriented_edges.push(OrientedEdge::new(edge_id, true));
        }

        let wire_id = store.add_wire(WireData {
            edges: oriented_edges,
            is_closed: self.close,
        });

        Ok(wire_id)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64) -> Point3 {
        Point3::new(x, y, 0.0)
    }

    #[test]
    fn triangle_closed_creates_3_edges() {
        let mut store = TopologyStore::new();
        let pts = vec![p(0.0, 0.0), p(4.0, 0.0), p(2.0, 3.0)];
        let wire_id = MakeWire::new(pts, true).execute(&mut store).unwrap();
        let wire = store.wire(wire_id).unwrap();
        assert_eq!(wire.edges.len(), 3);
        assert!(wire.is_closed);
    }

    #[test]
    fn open_3_points_creates_2_edges() {
        let mut store = TopologyStore::new();
        let pts = vec![p(0.0, 0.0), p(1.0, 0.0), p(2.0, 0.0)];
        let wire_id = MakeWire::new(pts, false).execute(&mut store).unwrap();
        let wire = store.wire(wire_id).unwrap();
        assert_eq!(wire.edges.len(), 2);
        assert!(!wire.is_closed);
    }

    #[test]
    fn two_points_closed_creates_2_edges() {
        let mut store = TopologyStore::new();
        let pts = vec![p(0.0, 0.0), p(1.0, 1.0)];
        let wire_id = MakeWire::new(pts, true).execute(&mut store).unwrap();
        let wire = store.wire(wire_id).unwrap();
        assert_eq!(wire.edges.len(), 2);
        assert!(wire.is_closed);
    }

    #[test]
    fn edge_parameters_match_distance() {
        let mut store = TopologyStore::new();
        let pts = vec![p(0.0, 0.0), p(3.0, 4.0)];
        let wire_id = MakeWire::new(pts, false).execute(&mut store).unwrap();
        let wire = store.wire(wire_id).unwrap();
        let edge = store.edge(wire.edges[0].edge).unwrap();
        assert!((edge.t_start).abs() < f64::EPSILON);
        assert!((edge.t_end - 5.0).abs() < 1e-10);
    }

    #[test]
    fn single_point_fails() {
        let mut store = TopologyStore::new();
        let result = MakeWire::new(vec![p(0.0, 0.0)], false).execute(&mut store);
        assert!(result.is_err());
    }

    #[test]
    fn empty_points_fails() {
        let mut store = TopologyStore::new();
        let result = MakeWire::new(vec![], false).execute(&mut store);
        assert!(result.is_err());
    }

    #[test]
    fn duplicate_consecutive_fails() {
        let mut store = TopologyStore::new();
        let pts = vec![p(0.0, 0.0), p(0.0, 0.0), p(1.0, 0.0)];
        let result = MakeWire::new(pts, false).execute(&mut store);
        assert!(result.is_err());
    }

    #[test]
    fn closing_edge_degenerate_fails() {
        let mut store = TopologyStore::new();
        let pts = vec![p(1.0, 0.0), p(2.0, 0.0), p(1.0, 0.0)];
        let result = MakeWire::new(pts, true).execute(&mut store);
        assert!(result.is_err());
    }
}
