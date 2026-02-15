use crate::error::Result;
use crate::math::Point3;
use crate::topology::{EdgeId, TopologyStore};

/// Result of a closest point query.
#[derive(Debug, Clone, Copy)]
pub struct ClosestPointResult {
    /// The closest point on the curve.
    pub point: Point3,
    /// The parameter value at the closest point.
    pub parameter: f64,
    /// The distance from the query point to the closest point.
    pub distance: f64,
}

/// Finds the closest point on a curve to a given point.
pub struct ClosestPointOnCurve {
    edge: EdgeId,
    point: Point3,
}

impl ClosestPointOnCurve {
    /// Creates a new `ClosestPointOnCurve` query.
    #[must_use]
    pub fn new(edge: EdgeId, point: Point3) -> Self {
        Self { edge, point }
    }

    /// Executes the query, returning the closest point result.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &TopologyStore) -> Result<ClosestPointResult> {
        let _ = (self.edge, self.point);
        todo!()
    }
}
