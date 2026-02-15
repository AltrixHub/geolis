use crate::error::Result;
use crate::math::Point3;
use crate::topology::{SolidId, TopologyStore};

/// Scales a solid uniformly from a center point.
pub struct Scale {
    solid: SolidId,
    center: Point3,
    factor: f64,
}

impl Scale {
    /// Creates a new `Scale` operation.
    #[must_use]
    pub fn new(solid: SolidId, center: Point3, factor: f64) -> Self {
        Self {
            solid,
            center,
            factor,
        }
    }

    /// Executes the scaling, modifying the solid in-place.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<()> {
        let _ = (self.solid, self.center, self.factor);
        todo!()
    }
}
