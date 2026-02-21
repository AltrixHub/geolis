use crate::error::{OperationError, Result};
use crate::geometry::curve::{Curve, Line};
use crate::math::{Point3, TOLERANCE};
use crate::topology::{EdgeCurve, EdgeData, EdgeId, TopologyStore, VertexData};

/// Trims an edge at the given parameter values.
///
/// Creates a new edge that spans only the `[t_start, t_end]` portion
/// of the original edge's curve, with new vertices at the trimmed endpoints.
pub struct Trim {
    edge: EdgeId,
    t_start: f64,
    t_end: f64,
}

impl Trim {
    /// Creates a new `Trim` operation.
    #[must_use]
    pub fn new(edge: EdgeId, t_start: f64, t_end: f64) -> Self {
        Self {
            edge,
            t_start,
            t_end,
        }
    }

    /// Executes the trim, creating the trimmed edge in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The edge is not found
    /// - The trim range is degenerate (start ≈ end)
    /// - The parameters are outside the edge's domain
    pub fn execute(&self, store: &mut TopologyStore) -> Result<EdgeId> {
        let edge = store.edge(self.edge)?;
        let curve = edge.curve.clone();
        let orig_start = edge.t_start;
        let orig_end = edge.t_end;

        // Validate range
        if (self.t_end - self.t_start).abs() < TOLERANCE {
            return Err(OperationError::InvalidInput(
                "trim range is degenerate (start ≈ end)".into(),
            )
            .into());
        }

        // Validate parameters are within the edge's domain (with tolerance)
        let t_min = orig_start.min(orig_end) - TOLERANCE;
        let t_max = orig_start.max(orig_end) + TOLERANCE;
        if self.t_start < t_min || self.t_start > t_max || self.t_end < t_min || self.t_end > t_max
        {
            return Err(OperationError::InvalidInput(format!(
                "trim parameters [{}, {}] outside edge domain [{}, {}]",
                self.t_start, self.t_end, orig_start, orig_end,
            ))
            .into());
        }

        // Evaluate new endpoint positions
        let start_point = evaluate_curve(&curve, self.t_start)?;
        let end_point = evaluate_curve(&curve, self.t_end)?;

        // Create new vertices
        let start_vertex = store.add_vertex(VertexData::new(start_point));
        let end_vertex = store.add_vertex(VertexData::new(end_point));

        // Build new curve for the trimmed edge
        let new_edge_data = match &curve {
            EdgeCurve::Line(_) => {
                let direction = end_point - start_point;
                let t_end = direction.norm();
                let line = Line::new(start_point, direction)?;
                EdgeData {
                    start: start_vertex,
                    end: end_vertex,
                    curve: EdgeCurve::Line(line),
                    t_start: 0.0,
                    t_end,
                }
            }
            EdgeCurve::Arc(arc) => {
                // Keep the same arc geometry, just update the parameter bounds
                EdgeData {
                    start: start_vertex,
                    end: end_vertex,
                    curve: EdgeCurve::Arc(arc.clone()),
                    t_start: self.t_start,
                    t_end: self.t_end,
                }
            }
        };

        Ok(store.add_edge(new_edge_data))
    }
}

/// Evaluates a point on an edge curve at parameter t.
fn evaluate_curve(curve: &EdgeCurve, t: f64) -> Result<Point3> {
    match curve {
        EdgeCurve::Line(line) => line.evaluate(t),
        EdgeCurve::Arc(arc) => arc.evaluate(t),
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
    fn trim_line_midpoint() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        // Trim to the middle third [3.333, 6.666]
        let trimmed = Trim::new(edge_id, 3.0, 7.0)
            .execute(&mut store)
            .unwrap();

        let edge = store.edge(trimmed).unwrap();
        let start = store.vertex(edge.start).unwrap().point;
        let end = store.vertex(edge.end).unwrap().point;
        assert!((start.x - 3.0).abs() < 1e-10);
        assert!((end.x - 7.0).abs() < 1e-10);
    }

    #[test]
    fn trim_degenerate_returns_error() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        let result = Trim::new(edge_id, 5.0, 5.0).execute(&mut store);
        assert!(result.is_err());
    }

    #[test]
    fn trim_out_of_range_returns_error() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        let result = Trim::new(edge_id, -5.0, 5.0).execute(&mut store);
        assert!(result.is_err());
    }
}
