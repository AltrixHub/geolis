use crate::error::Result;
use crate::topology::{TopologyStore, WireId};

/// Offsets a 2D wire by a given distance.
pub struct WireOffset2D {
    wire: WireId,
    distance: f64,
}

impl WireOffset2D {
    /// Creates a new `WireOffset2D` operation.
    #[must_use]
    pub fn new(wire: WireId, distance: f64) -> Self {
        Self { wire, distance }
    }

    /// Executes the offset, creating the result wire in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<WireId> {
        let _ = (self.wire, self.distance);
        todo!()
    }
}
