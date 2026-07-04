//! Through-cut Subtract boolean for NURBS-faced solids (P6 + F5 Phase B).
//!
//! Scope: the **through-cut Subtract** class. The tool passes fully through
//! the target (or ends inside it — pocket). A tube-like tool with one
//! periodic side face cuts each target face in closed SSI loops; a multi-face
//! (box-like) tool with planar side faces meeting at kink edges cuts OPEN
//! branches that are chained across adjacent tool faces into closed loops
//! ([`stitch`]). The result is the target with trim holes punched where the
//! tool passes through, plus the tool's band faces (one per periodic tool
//! face, or one fragment per crossed tool side face) forming the hole walls.
//!
//! Every precondition violation (open branch off the tool kinks, partial cut,
//! Union, NURBS↔planar splitting beyond caps, cap-face intersection) returns
//! an explicit [`crate::error::OperationError::Failed`] naming the
//! unsupported case. No silent wrong geometry, no panics.
//!
//! Pipeline:
//! 1. [`loops`] — SSI over face pairs → validated closed loops grouped per
//!    tool side face, plus open branches chained by [`stitch`] into closed
//!    multi-face loops grouped per face set.
//! 2. [`punch`] — convert each loop's target-UV trace into a trim hole on the
//!    target face plus a 3D hole wire (one edge per chain segment).
//! 3. [`band`] — build the tool's hole-wall face(s) between each loop pair.
//! 4. [`assemble`] — collect result faces into a new shell + solid.

pub(crate) mod assemble;
pub(crate) mod band;
pub(crate) mod intersect;
pub(crate) mod loops;
pub(crate) mod pocket;
pub(crate) mod punch;
pub(crate) mod stitch;

use crate::error::{OperationError, Result};
use crate::topology::{SolidId, TopologyStore};

use super::select::BooleanOp;

/// Routes a boolean operation on (at least one) NURBS-faced solid.
///
/// The through-cut [`BooleanOp::Subtract`] (keep-outside) and
/// [`BooleanOp::Intersect`] (keep-inside) are supported; Union returns an
/// explicit unsupported error, as does any operation that violates the
/// through-cut preconditions.
pub(crate) fn try_boolean(
    store: &mut TopologyStore,
    solid_a: SolidId,
    solid_b: SolidId,
    op: BooleanOp,
    op_id: Option<&crate::topology::OpId>,
) -> Result<SolidId> {
    match op {
        BooleanOp::Subtract => assemble::subtract_through_cut(store, solid_a, solid_b, op_id),
        BooleanOp::Intersect => intersect::intersect_through_cut(store, solid_a, solid_b, op_id),
        BooleanOp::Union => Err(OperationError::Failed(
            "union of NURBS-faced solids is not supported (through-cut subtract/intersect only)"
                .into(),
        )
        .into()),
    }
}
