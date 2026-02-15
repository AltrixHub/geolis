use crate::math::Point3;

slotmap::new_key_type! {
    /// Unique identifier for a vertex in the topology store.
    pub struct VertexId;
}

/// Data associated with a topological vertex.
#[derive(Debug, Clone)]
pub struct VertexData {
    /// The 3D position of the vertex.
    pub point: Point3,
}

impl VertexData {
    /// Creates a new vertex at the given point.
    #[must_use]
    pub fn new(point: Point3) -> Self {
        Self { point }
    }
}
