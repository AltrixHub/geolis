use crate::error::Result;
use crate::math::Vector3;
use crate::topology::{SolidId, TopologyStore};

/// Translates a solid by a displacement vector.
pub struct Translate {
    solid: SolidId,
    displacement: Vector3,
}

impl Translate {
    /// Creates a new `Translate` operation.
    #[must_use]
    pub fn new(solid: SolidId, displacement: Vector3) -> Self {
        Self {
            solid,
            displacement,
        }
    }

    /// Executes the translation, modifying the solid in-place.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<()> {
        let _ = (self.solid, self.displacement);
        todo!()
    }
}
