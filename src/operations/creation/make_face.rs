use crate::error::Result;
use crate::topology::{FaceId, TopologyStore, WireId};

/// Creates a face from a wire boundary and a surface.
pub struct MakeFace {
    outer_wire: WireId,
    inner_wires: Vec<WireId>,
}

impl MakeFace {
    /// Creates a new `MakeFace` operation.
    #[must_use]
    pub fn new(outer_wire: WireId, inner_wires: Vec<WireId>) -> Self {
        Self {
            outer_wire,
            inner_wires,
        }
    }

    /// Executes the operation, creating the face in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<FaceId> {
        let _ = (self.outer_wire, &self.inner_wires);
        todo!()
    }
}
