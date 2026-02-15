use crate::error::Result;
use crate::math::Matrix4;
use crate::topology::{SolidId, TopologyStore};

/// Applies an arbitrary 4x4 transformation matrix to a solid.
pub struct GeneralTransform {
    solid: SolidId,
    matrix: Matrix4,
}

impl GeneralTransform {
    /// Creates a new `GeneralTransform` operation.
    #[must_use]
    pub fn new(solid: SolidId, matrix: Matrix4) -> Self {
        Self { solid, matrix }
    }

    /// Executes the transformation, modifying the solid in-place.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<()> {
        let _ = (self.solid, self.matrix);
        todo!()
    }
}
