use super::shell::ShellId;

slotmap::new_key_type! {
    /// Unique identifier for a solid in the topology store.
    pub struct SolidId;
}

/// Data associated with a topological solid.
///
/// A solid is a bounded volume enclosed by one or more shells.
/// The first shell is the outer shell; additional shells represent voids.
#[derive(Debug, Clone)]
pub struct SolidData {
    /// The outer shell of the solid.
    pub outer_shell: ShellId,
    /// Inner shells representing voids within the solid.
    pub inner_shells: Vec<ShellId>,
}
