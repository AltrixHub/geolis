use crate::error::Result;
use crate::topology::{SolidId, TopologyStore};

/// Computes the boolean subtraction of one solid from another.
pub struct Subtract {
    solid_a: SolidId,
    solid_b: SolidId,
}

impl Subtract {
    /// Creates a new `Subtract` operation (A - B).
    #[must_use]
    pub fn new(solid_a: SolidId, solid_b: SolidId) -> Self {
        Self { solid_a, solid_b }
    }

    /// Executes the subtraction, creating the result solid in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<SolidId> {
        let _ = (self.solid_a, self.solid_b);
        todo!()
    }
}
