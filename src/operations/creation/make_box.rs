use crate::error::Result;
use crate::math::Point3;
use crate::topology::{SolidId, TopologyStore};

/// Creates a box solid from two corner points.
pub struct MakeBox {
    min_corner: Point3,
    max_corner: Point3,
}

impl MakeBox {
    /// Creates a new `MakeBox` operation.
    #[must_use]
    pub fn new(min_corner: Point3, max_corner: Point3) -> Self {
        Self {
            min_corner,
            max_corner,
        }
    }

    /// Executes the operation, creating the box in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<SolidId> {
        let _ = (self.min_corner, self.max_corner);
        todo!()
    }
}
