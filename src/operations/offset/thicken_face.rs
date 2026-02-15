use crate::error::Result;
use crate::topology::{FaceId, SolidId, TopologyStore};

/// Thickens a face into a solid by offsetting in the normal direction.
pub struct ThickenFace {
    face: FaceId,
    thickness: f64,
}

impl ThickenFace {
    /// Creates a new `ThickenFace` operation.
    #[must_use]
    pub fn new(face: FaceId, thickness: f64) -> Self {
        Self { face, thickness }
    }

    /// Executes the thickening, creating the result solid in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<SolidId> {
        let _ = (self.face, self.thickness);
        todo!()
    }
}
