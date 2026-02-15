use crate::error::Result;
use crate::math::{Point3, Vector3};
use crate::topology::{SolidId, TopologyStore};

/// Rotates a solid around an axis.
pub struct Rotate {
    solid: SolidId,
    axis_origin: Point3,
    axis_direction: Vector3,
    angle: f64,
}

impl Rotate {
    /// Creates a new `Rotate` operation.
    ///
    /// * `angle` - Rotation angle in radians.
    #[must_use]
    pub fn new(solid: SolidId, axis_origin: Point3, axis_direction: Vector3, angle: f64) -> Self {
        Self {
            solid,
            axis_origin,
            axis_direction,
            angle,
        }
    }

    /// Executes the rotation, modifying the solid in-place.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<()> {
        let _ = (self.solid, self.axis_origin, self.axis_direction, self.angle);
        todo!()
    }
}
