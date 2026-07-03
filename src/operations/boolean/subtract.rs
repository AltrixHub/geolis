use crate::error::Result;
use crate::topology::{SolidId, TopologyStore};

use super::select::BooleanOp;

/// Computes the boolean subtraction of one solid from another.
pub struct Subtract {
    solid_a: SolidId,
    solid_b: SolidId,
    op_id: Option<crate::topology::OpId>,
}

impl Subtract {
    /// Creates a new `Subtract` operation (A - B).
    #[must_use]
    pub fn new(solid_a: SolidId, solid_b: SolidId) -> Self {
        Self {
            solid_a,
            solid_b,
            op_id: None,
        }
    }

    /// Evolves persistent names through this boolean under the
    /// caller-supplied operation identity.
    #[must_use]
    pub fn with_op_id(mut self, op: crate::topology::OpId) -> Self {
        self.op_id = Some(op);
        self
    }

    /// Executes the subtraction, creating the result solid in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        crate::operations::boolean::engine::boolean_execute_named(
            store,
            self.solid_a,
            self.solid_b,
            BooleanOp::Subtract,
            self.op_id.as_ref(),
        )
    }
}
