use crate::error::Result;
use crate::topology::{EdgeId, TopologyStore};

/// Offsets a 2D curve by a given distance.
pub struct CurveOffset2D {
    edge: EdgeId,
    distance: f64,
}

impl CurveOffset2D {
    /// Creates a new `CurveOffset2D` operation.
    #[must_use]
    pub fn new(edge: EdgeId, distance: f64) -> Self {
        Self { edge, distance }
    }

    /// Executes the offset, creating the result edge in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<EdgeId> {
        let _ = (self.edge, self.distance);
        todo!()
    }
}
