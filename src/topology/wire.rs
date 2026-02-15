use super::edge::EdgeId;

slotmap::new_key_type! {
    /// Unique identifier for a wire in the topology store.
    pub struct WireId;
}

/// An edge with orientation information within a wire.
#[derive(Debug, Clone, Copy)]
pub struct OrientedEdge {
    /// The edge identifier.
    pub edge: EdgeId,
    /// If `true`, the edge is traversed in its natural direction (start → end).
    /// If `false`, the edge is traversed in reverse (end → start).
    pub forward: bool,
}

impl OrientedEdge {
    /// Creates a new oriented edge.
    #[must_use]
    pub fn new(edge: EdgeId, forward: bool) -> Self {
        Self { edge, forward }
    }
}

/// Data associated with a topological wire.
///
/// A wire is an ordered sequence of oriented edges forming a connected path.
/// It may be open or closed.
#[derive(Debug, Clone)]
pub struct WireData {
    /// The ordered sequence of oriented edges.
    pub edges: Vec<OrientedEdge>,
    /// Whether this wire forms a closed loop.
    pub is_closed: bool,
}
