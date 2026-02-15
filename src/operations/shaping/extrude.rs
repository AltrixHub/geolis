use crate::error::Result;
use crate::math::Vector3;
use crate::topology::{FaceId, SolidId, TopologyStore};

/// Extrudes a face along a direction vector to create a solid.
pub struct Extrude {
    face: FaceId,
    direction: Vector3,
}

impl Extrude {
    /// Creates a new `Extrude` operation.
    #[must_use]
    pub fn new(face: FaceId, direction: Vector3) -> Self {
        Self { face, direction }
    }

    /// Executes the extrusion, creating the solid in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<SolidId> {
        let _ = (self.face, self.direction);
        todo!()
    }
}
