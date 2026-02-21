use crate::error::Result;
use crate::geometry::curve::Curve;
use crate::topology::{EdgeCurve, EdgeId, TopologyStore};

/// Computes the length of a curve (edge).
pub struct Length {
    edge: EdgeId,
}

impl Length {
    /// Creates a new `Length` query.
    #[must_use]
    pub fn new(edge: EdgeId) -> Self {
        Self { edge }
    }

    /// Executes the query, returning the curve length.
    ///
    /// For a `Line`, this is `|t_end - t_start|`.
    /// For an `Arc`, this is `radius * |sweep_angle|`.
    ///
    /// # Errors
    ///
    /// Returns an error if the edge is not found.
    pub fn execute(&self, store: &TopologyStore) -> Result<f64> {
        let edge = store.edge(self.edge)?;
        match &edge.curve {
            EdgeCurve::Line(_) => {
                // Line is parameterized by arc length, so length = |t_end - t_start|
                Ok((edge.t_end - edge.t_start).abs())
            }
            EdgeCurve::Arc(arc) => {
                let domain = arc.domain();
                let sweep = (domain.t_max - domain.t_min).abs();
                Ok(arc.radius() * sweep)
            }
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
    fn line_length_3_4_5() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 4.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        let len = Length::new(edge_id).execute(&store).unwrap();
        assert!((len - 5.0).abs() < 1e-10);
    }

    #[test]
    fn unit_edge_length() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        let len = Length::new(edge_id).execute(&store).unwrap();
        assert!((len - 1.0).abs() < 1e-10);
    }
}
