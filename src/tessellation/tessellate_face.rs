use crate::error::Result;
use crate::topology::{FaceId, TopologyStore};

use super::{TessellationParams, TriangleMesh};

/// Tessellates a face into a triangle mesh.
pub struct TessellateFace {
    face: FaceId,
    params: TessellationParams,
}

impl TessellateFace {
    /// Creates a new `TessellateFace` operation.
    #[must_use]
    pub fn new(face: FaceId, params: TessellationParams) -> Self {
        Self { face, params }
    }

    /// Executes the tessellation, returning a triangle mesh.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &TopologyStore) -> Result<TriangleMesh> {
        let _ = (self.face, self.params);
        todo!()
    }
}
