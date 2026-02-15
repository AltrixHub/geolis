use crate::error::Result;
use crate::topology::{EdgeId, TopologyStore};

use super::{Polyline, TessellationParams};

/// Tessellates a curve (edge) into a polyline.
pub struct TessellateCurve {
    edge: EdgeId,
    params: TessellationParams,
}

impl TessellateCurve {
    /// Creates a new `TessellateCurve` operation.
    #[must_use]
    pub fn new(edge: EdgeId, params: TessellationParams) -> Self {
        Self { edge, params }
    }

    /// Executes the tessellation, returning a polyline.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &TopologyStore) -> Result<Polyline> {
        let _ = (self.edge, self.params);
        todo!()
    }
}
