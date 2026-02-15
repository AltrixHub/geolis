use crate::geometry::curve::{Arc, Line};

use super::vertex::VertexId;

slotmap::new_key_type! {
    /// Unique identifier for an edge in the topology store.
    pub struct EdgeId;
}

/// The geometric curve associated with an edge.
#[derive(Debug, Clone)]
pub enum EdgeCurve {
    /// A line segment.
    Line(Line),
    /// A circular arc.
    Arc(Arc),
}

/// Data associated with a topological edge.
///
/// An edge connects two vertices and carries a geometric curve
/// that defines the shape of the edge between them.
#[derive(Debug, Clone)]
pub struct EdgeData {
    /// Start vertex of the edge.
    pub start: VertexId,
    /// End vertex of the edge.
    pub end: VertexId,
    /// The geometric curve defining this edge's shape.
    pub curve: EdgeCurve,
    /// Parameter on the curve corresponding to the start vertex.
    pub t_start: f64,
    /// Parameter on the curve corresponding to the end vertex.
    pub t_end: f64,
}
