use crate::geometry::surface::Plane;

use super::wire::WireId;

slotmap::new_key_type! {
    /// Unique identifier for a face in the topology store.
    pub struct FaceId;
}

/// The geometric surface associated with a face.
#[derive(Debug, Clone)]
pub enum FaceSurface {
    /// A planar surface.
    Plane(Plane),
}

/// Data associated with a topological face.
///
/// A face is a bounded region on a surface, defined by an outer wire
/// and optionally inner wires (holes).
#[derive(Debug, Clone)]
pub struct FaceData {
    /// The geometric surface on which this face lies.
    pub surface: FaceSurface,
    /// The outer boundary wire.
    pub outer_wire: WireId,
    /// Inner boundary wires (holes).
    pub inner_wires: Vec<WireId>,
    /// If `true`, the face normal agrees with the surface normal.
    pub same_sense: bool,
}
