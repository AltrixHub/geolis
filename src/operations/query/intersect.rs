use crate::error::Result;
use crate::math::Point3;
use crate::topology::{EdgeId, TopologyStore};

/// Result of a curve-curve intersection.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionResult {
    /// The intersection point.
    pub point: Point3,
    /// Parameter on the first curve.
    pub t1: f64,
    /// Parameter on the second curve.
    pub t2: f64,
}

/// Computes intersections between two curves.
pub struct CurveCurveIntersect {
    edge_a: EdgeId,
    edge_b: EdgeId,
}

impl CurveCurveIntersect {
    /// Creates a new `CurveCurveIntersect` query.
    #[must_use]
    pub fn new(edge_a: EdgeId, edge_b: EdgeId) -> Self {
        Self { edge_a, edge_b }
    }

    /// Executes the query, returning all intersection points.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &TopologyStore) -> Result<Vec<IntersectionResult>> {
        let _ = (self.edge_a, self.edge_b);
        todo!()
    }
}
