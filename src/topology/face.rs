use crate::geometry::nurbs::{NurbsCurve2D, NurbsSurface};
use crate::geometry::surface::{Cone, Cylinder, Plane, Sphere, Torus};

use super::edge::EdgeId;
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
    /// A cylindrical surface.
    Cylinder(Cylinder),
    /// A conical surface.
    Cone(Cone),
    /// A spherical surface.
    Sphere(Sphere),
    /// A toroidal surface.
    Torus(Torus),
    /// A free-form NURBS surface.
    Nurbs(NurbsSurface),
}

/// The UV image of a boundary edge in one face's parameter space.
///
/// A shared edge is referenced by faces with different surfaces, so the
/// edge→UV mapping is stored per face, not on the edge or the (shared) wire.
///
/// Same-parameter convention: for every `t` in the edge's parameter domain,
/// `surface.point_at(pcurve.point_at(t)) == edge_curve.point_at(t)` (within
/// tolerance). Creation ops guarantee this exactly by construction — extrude /
/// revolve / ruled surfaces preserve the profile parameterization on their
/// boundary isocurves.
#[derive(Debug, Clone)]
pub struct FacePcurve {
    /// The shared boundary edge this pcurve maps.
    pub edge: EdgeId,
    /// The UV curve, parameterized identically to the edge curve.
    pub curve: NurbsCurve2D,
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
    /// UV-space trim geometry. `None` for analytic faces and untrimmed NURBS
    /// faces (full parameter domain); `Some` for trimmed NURBS faces.
    pub trim: Option<super::trim::FaceTrim>,
    /// Per-edge UV images for this face's boundary edges. Empty for faces
    /// whose builder predates shared-edge topology; consumers fall back to
    /// geometric boundary sampling in that case.
    pub pcurves: Vec<FacePcurve>,
}

impl FaceData {
    /// Returns this face's pcurve for `edge`, if recorded.
    #[must_use]
    pub fn pcurve_for(&self, edge: EdgeId) -> Option<&NurbsCurve2D> {
        self.pcurves
            .iter()
            .find(|p| p.edge == edge)
            .map(|p| &p.curve)
    }
}
