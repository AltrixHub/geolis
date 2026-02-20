use crate::error::Result;
use crate::topology::{ShellId, SolidData, SolidId, TopologyStore};

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
    /// Returns [`TopologyError::EntityNotFound`] if any shell ID is invalid.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        // Validate that all shells exist
        let _ = store.shell(self.outer_shell)?;
        for &shell_id in &self.inner_shells {
            let _ = store.shell(shell_id)?;
        }

        let solid_id = store.add_solid(SolidData {
            outer_shell: self.outer_shell,
            inner_shells: self.inner_shells.clone(),
        });

        Ok(solid_id)
    }
}
