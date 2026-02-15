use crate::error::Result;
use crate::math::Point3;
use crate::topology::{SolidId, TopologyStore};

/// An axis-aligned bounding box.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Minimum corner of the bounding box.
    pub min: Point3,
    /// Maximum corner of the bounding box.
    pub max: Point3,
}

/// Computes the axis-aligned bounding box of a solid.
pub struct BoundingBox {
    solid: SolidId,
}

impl BoundingBox {
    /// Creates a new `BoundingBox` query.
    #[must_use]
    pub fn new(solid: SolidId) -> Self {
        Self { solid }
    }

    /// Executes the query, returning the AABB.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &TopologyStore) -> Result<Aabb> {
        let _ = self.solid;
        todo!()
    }
}
