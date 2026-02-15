use super::face::FaceId;

slotmap::new_key_type! {
    /// Unique identifier for a shell in the topology store.
    pub struct ShellId;
}

/// Data associated with a topological shell.
///
/// A shell is a connected set of faces forming a surface boundary.
/// It may be open or closed.
#[derive(Debug, Clone)]
pub struct ShellData {
    /// The faces that make up this shell.
    pub faces: Vec<FaceId>,
    /// Whether this shell is closed (watertight).
    pub is_closed: bool,
}
