use crate::error::Result;
use crate::topology::{SolidId, TopologyStore};

use super::engine::boolean_execute;
use super::select::BooleanOp;

/// Computes the boolean union of two solids.
pub struct Union {
    solid_a: SolidId,
    solid_b: SolidId,
}

impl Union {
    /// Creates a new `Union` operation.
    #[must_use]
    pub fn new(solid_a: SolidId, solid_b: SolidId) -> Self {
        Self { solid_a, solid_b }
    }

    /// Executes the union, creating the result solid in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        boolean_execute(store, self.solid_a, self.solid_b, BooleanOp::Union)
    }
}
