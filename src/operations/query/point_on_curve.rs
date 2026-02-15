use crate::error::Result;
use crate::math::Point3;
use crate::topology::{EdgeId, TopologyStore};

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
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &TopologyStore) -> Result<Point3> {
        let _ = (self.edge, self.t);
        todo!()
    }
}
