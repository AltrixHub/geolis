use crate::error::Result;
use crate::topology::{EdgeId, TopologyStore};

/// Computes the length of a curve (edge).
pub struct Length {
    edge: EdgeId,
}

impl Length {
    /// Creates a new `Length` query.
    #[must_use]
    pub fn new(edge: EdgeId) -> Self {
        Self { edge }
    }

    /// Executes the query, returning the curve length.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &TopologyStore) -> Result<f64> {
        let _ = self.edge;
        todo!()
    }
}
