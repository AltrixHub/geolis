//! UV-space trim geometry for NURBS faces.
//!
//! A trimmed NURBS face restricts the surface to a region of its parameter
//! domain bounded by closed loops of 2D (UV-space) curves — pcurves. The outer
//! loop runs counter-clockwise and bounds the kept region; each hole loop runs
//! clockwise. Trim geometry lives on the face because the face owns its UV
//! space; the 3D wires/edges remain the source of truth for topological
//! connectivity, while [`FaceTrim`] is the UV-space companion required for
//! trimmed tessellation (and, in later phases, UV-space splitting).

use crate::geometry::nurbs::NurbsCurve2D;

/// A closed loop of UV-space trim curves bounding a region of a NURBS face.
///
/// The curves are traversed head-to-tail in UV space and together form a
/// closed boundary. Counter-clockwise winding denotes an outer boundary,
/// clockwise winding denotes a hole.
#[derive(Debug, Clone)]
pub struct TrimLoop {
    /// Closed 2D curves traversed head-to-tail in UV space.
    pub curves: Vec<NurbsCurve2D>,
}

impl TrimLoop {
    /// Creates a new trim loop from a sequence of UV-space curves.
    #[must_use]
    pub fn new(curves: Vec<NurbsCurve2D>) -> Self {
        Self { curves }
    }
}

/// Trim data for a NURBS face: one outer loop plus zero or more holes.
#[derive(Debug, Clone)]
pub struct FaceTrim {
    /// The outer boundary loop (counter-clockwise).
    pub outer: TrimLoop,
    /// Hole loops (clockwise).
    pub holes: Vec<TrimLoop>,
}

impl FaceTrim {
    /// Creates a new face trim from an outer loop and a set of holes.
    #[must_use]
    pub fn new(outer: TrimLoop, holes: Vec<TrimLoop>) -> Self {
        Self { outer, holes }
    }
}
