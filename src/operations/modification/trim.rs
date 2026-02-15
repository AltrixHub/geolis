use crate::error::Result;
use crate::topology::{EdgeId, TopologyStore};

/// Trims an edge at the given parameter values.
pub struct Trim {
    edge: EdgeId,
    t_start: f64,
    t_end: f64,
}

impl Trim {
    /// Creates a new `Trim` operation.
    #[must_use]
    pub fn new(edge: EdgeId, t_start: f64, t_end: f64) -> Self {
        Self {
            edge,
            t_start,
            t_end,
        }
    }

    /// Executes the trim, creating the trimmed edge in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<EdgeId> {
        let _ = (self.edge, self.t_start, self.t_end);
        todo!()
    }
}
