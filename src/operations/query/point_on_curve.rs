use crate::error::Result;
use crate::geometry::curve::Curve;
use crate::math::Point3;
use crate::topology::{EdgeCurve, EdgeId, TopologyStore};

/// Evaluates a point on a curve at a given parameter.
pub struct PointOnCurve {
    edge: EdgeId,
    t: f64,
}

impl PointOnCurve {
    /// Creates a new `PointOnCurve` query.
    #[must_use]
    pub fn new(edge: EdgeId, t: f64) -> Self {
        Self { edge, t }
    }

    /// Executes the query, returning the 3D point.
    ///
    /// # Errors
    ///
    /// Returns an error if the edge is not found or the parameter is invalid.
    pub fn execute(&self, store: &TopologyStore) -> Result<Point3> {
        let edge = store.edge(self.edge)?;
        match &edge.curve {
            EdgeCurve::Line(line) => line.evaluate(self.t),
            EdgeCurve::Arc(arc) => arc.evaluate(self.t),
            EdgeCurve::Circle(circle) => circle.evaluate(self.t),
            EdgeCurve::Ellipse(ellipse) => ellipse.evaluate(self.t),
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::MakeWire;
    use crate::topology::TopologyStore;

    #[test]
    fn point_at_start_of_line_edge() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(3.0, 4.0, 0.0),
            ],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        let pt = PointOnCurve::new(edge_id, 0.0)
            .execute(&store)
            .unwrap();
        assert!((pt.x).abs() < 1e-10);
        assert!((pt.y).abs() < 1e-10);
    }

    #[test]
    fn point_at_end_of_line_edge() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(3.0, 4.0, 0.0),
            ],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        // t_end = 5.0 (distance from (0,0) to (3,4))
        let pt = PointOnCurve::new(edge_id, 5.0)
            .execute(&store)
            .unwrap();
        assert!((pt.x - 3.0).abs() < 1e-10);
        assert!((pt.y - 4.0).abs() < 1e-10);
    }
}
