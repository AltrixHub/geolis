use crate::topology::{SolidId, TopologyStore};

/// Validates the topological and geometric consistency of a solid.
pub struct IsValid {
    solid: SolidId,
}

impl IsValid {
    /// Creates a new `IsValid` query.
    #[must_use]
    pub fn new(solid: SolidId) -> Self {
        Self { solid }
    }

    /// Executes the validation, returning `true` if the solid is valid.
    #[must_use]
    pub fn execute(&self, _store: &TopologyStore) -> bool {
        let _ = self.solid;
        todo!()
    }
}
