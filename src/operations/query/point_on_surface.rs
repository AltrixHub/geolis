use crate::error::Result;
use crate::math::Point3;
use crate::topology::{FaceId, TopologyStore};

/// Evaluates a point on a surface at given parameters.
pub struct PointOnSurface {
    face: FaceId,
    u: f64,
    v: f64,
}

impl PointOnSurface {
    /// Creates a new `PointOnSurface` query.
    #[must_use]
    pub fn new(face: FaceId, u: f64, v: f64) -> Self {
        Self { face, u, v }
    }

    /// Executes the query, returning the 3D point.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &TopologyStore) -> Result<Point3> {
        let _ = (self.face, self.u, self.v);
        todo!()
    }
}
