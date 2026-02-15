pub mod edge;
pub mod face;
pub mod shell;
pub mod solid;
pub mod vertex;
pub mod wire;

pub use edge::{EdgeCurve, EdgeData, EdgeId};
pub use face::{FaceData, FaceId, FaceSurface};
pub use shell::{ShellData, ShellId};
pub use solid::{SolidData, SolidId};
pub use vertex::{VertexData, VertexId};
pub use wire::{OrientedEdge, WireData, WireId};

use crate::error::TopologyError;
use slotmap::SlotMap;

/// Central arena that owns all topological entities.
///
/// Entities reference each other via typed IDs (generational indices),
/// avoiding self-referential structures and enabling safe mutation.
#[derive(Debug, Default)]
pub struct TopologyStore {
    vertices: SlotMap<VertexId, VertexData>,
    edges: SlotMap<EdgeId, EdgeData>,
    wires: SlotMap<WireId, WireData>,
    faces: SlotMap<FaceId, FaceData>,
    shells: SlotMap<ShellId, ShellData>,
    solids: SlotMap<SolidId, SolidData>,
}

impl TopologyStore {
    /// Creates a new, empty topology store.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    // --- Vertex operations ---

    /// Inserts a vertex and returns its ID.
    pub fn add_vertex(&mut self, data: VertexData) -> VertexId {
        self.vertices.insert(data)
    }

    /// Returns a reference to the vertex data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn vertex(&self, id: VertexId) -> Result<&VertexData, TopologyError> {
        self.vertices
            .get(id)
            .ok_or_else(|| TopologyError::EntityNotFound("vertex".into()))
    }

    /// Returns a mutable reference to the vertex data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn vertex_mut(&mut self, id: VertexId) -> Result<&mut VertexData, TopologyError> {
        self.vertices
            .get_mut(id)
            .ok_or_else(|| TopologyError::EntityNotFound("vertex".into()))
    }

    // --- Edge operations ---

    /// Inserts an edge and returns its ID.
    pub fn add_edge(&mut self, data: EdgeData) -> EdgeId {
        self.edges.insert(data)
    }

    /// Returns a reference to the edge data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn edge(&self, id: EdgeId) -> Result<&EdgeData, TopologyError> {
        self.edges
            .get(id)
            .ok_or_else(|| TopologyError::EntityNotFound("edge".into()))
    }

    /// Returns a mutable reference to the edge data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn edge_mut(&mut self, id: EdgeId) -> Result<&mut EdgeData, TopologyError> {
        self.edges
            .get_mut(id)
            .ok_or_else(|| TopologyError::EntityNotFound("edge".into()))
    }

    // --- Wire operations ---

    /// Inserts a wire and returns its ID.
    pub fn add_wire(&mut self, data: WireData) -> WireId {
        self.wires.insert(data)
    }

    /// Returns a reference to the wire data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn wire(&self, id: WireId) -> Result<&WireData, TopologyError> {
        self.wires
            .get(id)
            .ok_or_else(|| TopologyError::EntityNotFound("wire".into()))
    }

    /// Returns a mutable reference to the wire data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn wire_mut(&mut self, id: WireId) -> Result<&mut WireData, TopologyError> {
        self.wires
            .get_mut(id)
            .ok_or_else(|| TopologyError::EntityNotFound("wire".into()))
    }

    // --- Face operations ---

    /// Inserts a face and returns its ID.
    pub fn add_face(&mut self, data: FaceData) -> FaceId {
        self.faces.insert(data)
    }

    /// Returns a reference to the face data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn face(&self, id: FaceId) -> Result<&FaceData, TopologyError> {
        self.faces
            .get(id)
            .ok_or_else(|| TopologyError::EntityNotFound("face".into()))
    }

    /// Returns a mutable reference to the face data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn face_mut(&mut self, id: FaceId) -> Result<&mut FaceData, TopologyError> {
        self.faces
            .get_mut(id)
            .ok_or_else(|| TopologyError::EntityNotFound("face".into()))
    }

    // --- Shell operations ---

    /// Inserts a shell and returns its ID.
    pub fn add_shell(&mut self, data: ShellData) -> ShellId {
        self.shells.insert(data)
    }

    /// Returns a reference to the shell data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn shell(&self, id: ShellId) -> Result<&ShellData, TopologyError> {
        self.shells
            .get(id)
            .ok_or_else(|| TopologyError::EntityNotFound("shell".into()))
    }

    /// Returns a mutable reference to the shell data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn shell_mut(&mut self, id: ShellId) -> Result<&mut ShellData, TopologyError> {
        self.shells
            .get_mut(id)
            .ok_or_else(|| TopologyError::EntityNotFound("shell".into()))
    }

    // --- Solid operations ---

    /// Inserts a solid and returns its ID.
    pub fn add_solid(&mut self, data: SolidData) -> SolidId {
        self.solids.insert(data)
    }

    /// Returns a reference to the solid data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn solid(&self, id: SolidId) -> Result<&SolidData, TopologyError> {
        self.solids
            .get(id)
            .ok_or_else(|| TopologyError::EntityNotFound("solid".into()))
    }

    /// Returns a mutable reference to the solid data, or an error if not found.
    ///
    /// # Errors
    ///
    /// Returns an error if the entity is not found in the store.
    pub fn solid_mut(&mut self, id: SolidId) -> Result<&mut SolidData, TopologyError> {
        self.solids
            .get_mut(id)
            .ok_or_else(|| TopologyError::EntityNotFound("solid".into()))
    }
}
