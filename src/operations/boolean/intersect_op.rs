use crate::error::Result;
use crate::topology::{SolidId, TopologyStore};

use super::select::BooleanOp;

/// Computes the boolean intersection of two solids.
pub struct Intersect {
    solid_a: SolidId,
    solid_b: SolidId,
    op_id: Option<crate::topology::OpId>,
}

impl Intersect {
    /// Creates a new `Intersect` operation.
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

    /// Executes the intersection, creating the result solid in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the solids don't overlap or the operation fails.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        crate::operations::boolean::engine::boolean_execute_named(
            store,
            self.solid_a,
            self.solid_b,
            BooleanOp::Intersect,
            self.op_id.as_ref(),
        )
    }
}
