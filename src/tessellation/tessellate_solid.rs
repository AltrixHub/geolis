use crate::error::Result;
use crate::topology::{SolidId, TopologyStore};

use super::{TessellationParams, TriangleMesh};

/// Tessellates all faces of a solid into a combined triangle mesh.
pub struct TessellateSolid {
    solid: SolidId,
    params: TessellationParams,
}

impl TessellateSolid {
    /// Creates a new `TessellateSolid` operation.
    #[must_use]
    pub fn new(solid: SolidId, params: TessellationParams) -> Self {
        Self { solid, params }
    }

    /// Executes the tessellation, returning a combined triangle mesh.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &TopologyStore) -> Result<TriangleMesh> {
        let _ = (self.solid, self.params);
        todo!()
    }
}
