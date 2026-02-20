use crate::error::Result;
use crate::topology::{SolidId, TopologyStore};

use super::{TessellateFace, TessellationParams, TriangleMesh};

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
    /// Returns an error if the solid or any of its faces cannot be tessellated.
    pub fn execute(&self, store: &TopologyStore) -> Result<TriangleMesh> {
        let solid = store.solid(self.solid)?;
        let shell = store.shell(solid.outer_shell)?;

        let mut combined = TriangleMesh::default();
        for &face_id in &shell.faces {
            let face_mesh = TessellateFace::new(face_id, self.params).execute(store)?;
            combined.merge(&face_mesh);
        }

        Ok(combined)
    }
}
