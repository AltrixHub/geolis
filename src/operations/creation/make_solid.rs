use crate::error::Result;
use crate::topology::{ShellId, SolidId, TopologyStore};

/// Creates a solid from shells.
pub struct MakeSolid {
    outer_shell: ShellId,
    inner_shells: Vec<ShellId>,
}

impl MakeSolid {
    /// Creates a new `MakeSolid` operation.
    #[must_use]
    pub fn new(outer_shell: ShellId, inner_shells: Vec<ShellId>) -> Self {
        Self {
            outer_shell,
            inner_shells,
        }
    }

    /// Executes the operation, creating the solid in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<SolidId> {
        let _ = (self.outer_shell, &self.inner_shells);
        todo!()
    }
}
