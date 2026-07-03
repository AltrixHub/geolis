//! Persistent topology names (topological naming).
//!
//! Slotmap ids are fresh on every rebuild, so cross-rebuild references (an
//! opening on a wall face, a per-face material) need rebuild-stable names.
//! Names are **derivational**: a pure function of the creating operation's
//! identity ([`OpId`], supplied by the caller — geolis never invents one),
//! the entity's role within that operation, and — for boolean products — the
//! parent names. Same inputs, same names, independent of allocation order.
//!
//! A boolean's result carries the target's names forward UNCHANGED (a punched
//! wall face is still "the wall's outer face"); new faces (band, pocket
//! floor) and new edges (hole rims) get names composed from their parents.
//! Resolution failure is `None` — no geometric best-match heuristics
//! (Kripac 1997 / Chen & Hoffmann 1995 / OCCT TNaming motivate the problem;
//! geolis's small deterministic op set lets derivational names replace their
//! heavyweight history machinery).

use std::collections::HashMap;
use std::fmt;
use std::sync::Arc;

use super::edge::EdgeId;
use super::face::FaceId;

/// Identity of the graph operation that created an entity, supplied by the
/// caller (revion passes its cognet node id). Rebuild-stable by construction.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub struct OpId(Arc<str>);

impl OpId {
    /// Creates an operation id from the caller's stable identifier.
    pub fn new(id: impl Into<Arc<str>>) -> Self {
        Self(id.into())
    }

    /// The raw identifier.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for OpId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// The role of a face within its creation operation.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum FaceRole {
    /// A side face, indexed by the op's deterministic side order (e.g. the
    /// curved wall: 0 = inner, 1 = outer, 2 = start end, 3 = end end).
    Side(u8),
    /// The cap at the extrusion start (`v0` end).
    CapStart,
    /// The cap at the extrusion end (`v1` end).
    CapEnd,
    /// The top face (slab / wall).
    Top,
    /// The bottom face (slab / wall).
    Bottom,
    /// The revolved wall surface.
    Wall,
}

/// The role of an edge within its creation operation.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum EdgeRole {
    /// The shared ring edge at the extrusion start (`v0` / first profile end).
    RingStart,
    /// The shared ring edge at the extrusion end (`v1` / last profile end).
    RingEnd,
}

/// A persistent, rebuild-stable name for a face.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum FaceName {
    /// A face born in a creation operation.
    Created {
        /// The creating operation.
        op: OpId,
        /// The face's role within that operation.
        role: FaceRole,
    },
    /// The band (hole / pocket wall) a boolean carved with a tool face.
    Band {
        /// The boolean operation.
        op: OpId,
        /// The tool side face the band lies on.
        tool_face: Box<FaceName>,
        /// Deterministic loop index (loops sorted by mean tool-`v`).
        loop_index: u32,
    },
    /// The pocket floor: the buried tool cap kept (sense-flipped) by a
    /// pocket subtract.
    Floor {
        /// The boolean operation.
        op: OpId,
        /// The buried cap's name.
        cap: Box<FaceName>,
    },
}

/// A persistent, rebuild-stable name for an edge.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum EdgeName {
    /// An edge born in a creation operation.
    Created {
        /// The creating operation.
        op: OpId,
        /// The edge's role within that operation.
        role: EdgeRole,
    },
    /// A hole-rim ring a boolean cut into a target face.
    CutRim {
        /// The boolean operation.
        op: OpId,
        /// The punched target face's name.
        target: Box<FaceName>,
        /// Deterministic loop index (loops sorted by mean tool-`v`).
        loop_index: u32,
    },
}

/// Bidirectional registry `PersistentName ↔ current slotmap id`.
///
/// Both directions stay bijective: registering a name that is already bound
/// rebinds it (the previous holder drops out), which is exactly the boolean
/// move semantics — the newest result owns the name.
#[derive(Debug, Default)]
pub struct NameRegistry {
    face_names: HashMap<FaceId, FaceName>,
    faces_by_name: HashMap<FaceName, FaceId>,
    edge_names: HashMap<EdgeId, EdgeName>,
    edges_by_name: HashMap<EdgeName, EdgeId>,
}

impl NameRegistry {
    /// Binds `name` to `face`, unbinding any previous holder of the name and
    /// any previous name of the face.
    pub fn bind_face(&mut self, face: FaceId, name: FaceName) {
        if let Some(old_face) = self.faces_by_name.remove(&name) {
            self.face_names.remove(&old_face);
        }
        if let Some(old_name) = self.face_names.remove(&face) {
            self.faces_by_name.remove(&old_name);
        }
        self.face_names.insert(face, name.clone());
        self.faces_by_name.insert(name, face);
    }

    /// Binds `name` to `edge`, unbinding any previous holders.
    pub fn bind_edge(&mut self, edge: EdgeId, name: EdgeName) {
        if let Some(old_edge) = self.edges_by_name.remove(&name) {
            self.edge_names.remove(&old_edge);
        }
        if let Some(old_name) = self.edge_names.remove(&edge) {
            self.edges_by_name.remove(&old_name);
        }
        self.edge_names.insert(edge, name.clone());
        self.edges_by_name.insert(name, edge);
    }

    /// Moves the name of `from` (if any) onto `to` — the boolean carry-over.
    pub fn transfer_face(&mut self, from: FaceId, to: FaceId) {
        if let Some(name) = self.face_names.remove(&from) {
            self.faces_by_name.remove(&name);
            self.bind_face(to, name);
        }
    }

    /// Resolves a face name to the current face id.
    #[must_use]
    pub fn face(&self, name: &FaceName) -> Option<FaceId> {
        self.faces_by_name.get(name).copied()
    }

    /// The current name of a face, if registered.
    #[must_use]
    pub fn name_of_face(&self, face: FaceId) -> Option<&FaceName> {
        self.face_names.get(&face)
    }

    /// Resolves an edge name to the current edge id.
    #[must_use]
    pub fn edge(&self, name: &EdgeName) -> Option<EdgeId> {
        self.edges_by_name.get(name).copied()
    }

    /// The current name of an edge, if registered.
    #[must_use]
    pub fn name_of_edge(&self, edge: EdgeId) -> Option<&EdgeName> {
        self.edge_names.get(&edge)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::surface::Plane;
    use crate::math::{Point3, Vector3};
    use crate::topology::{FaceData, FaceSurface, TopologyStore, WireId};

    /// Adds a minimal placeholder face (registry tests never dereference the
    /// wire, so a null wire id is fine).
    fn dummy_face(store: &mut TopologyStore) -> FaceId {
        let plane = Plane::new(Point3::origin(), Vector3::z(), Vector3::x()).unwrap();
        store.add_face(FaceData {
            surface: FaceSurface::Plane(plane),
            outer_wire: WireId::default(),
            inner_wires: vec![],
            same_sense: true,
            trim: None,
            pcurves: Vec::new(),
        })
    }

    fn outer(op: &str) -> FaceName {
        FaceName::Created {
            op: OpId::new(op),
            role: FaceRole::Side(1),
        }
    }

    #[test]
    fn bind_resolves_both_directions() {
        let mut store = TopologyStore::new();
        let face = dummy_face(&mut store);
        let mut reg = NameRegistry::default();
        reg.bind_face(face, outer("wall1"));
        assert_eq!(reg.face(&outer("wall1")), Some(face));
        assert_eq!(reg.name_of_face(face), Some(&outer("wall1")));
        assert_eq!(reg.face(&outer("wall2")), None);
    }

    #[test]
    fn rebinding_a_name_moves_it_off_the_old_face() {
        let mut store = TopologyStore::new();
        let old = dummy_face(&mut store);
        let new = dummy_face(&mut store);
        let mut reg = NameRegistry::default();
        reg.bind_face(old, outer("wall1"));
        reg.bind_face(new, outer("wall1"));
        assert_eq!(reg.face(&outer("wall1")), Some(new));
        assert_eq!(reg.name_of_face(old), None, "old holder must be unbound");
    }

    #[test]
    fn transfer_moves_the_name_to_the_copy() {
        let mut store = TopologyStore::new();
        let original = dummy_face(&mut store);
        let copy = dummy_face(&mut store);
        let mut reg = NameRegistry::default();
        reg.bind_face(original, outer("wall1"));
        reg.transfer_face(original, copy);
        assert_eq!(reg.face(&outer("wall1")), Some(copy));
        assert_eq!(reg.name_of_face(original), None);
        // Transferring from an unnamed face is a no-op.
        let unrelated = dummy_face(&mut store);
        reg.transfer_face(unrelated, original);
        assert_eq!(reg.face(&outer("wall1")), Some(copy));
    }
}
