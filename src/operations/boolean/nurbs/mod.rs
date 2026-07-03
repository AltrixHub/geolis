//! Through-cut Subtract boolean for NURBS-faced solids (P6).
//!
//! Scope: the **through-cut Subtract** class only. The tool (a tube-like
//! NURBS-faced solid) passes fully through the target; each tool side face cuts
//! the target in two closed SSI loops (entry + exit). The result is the target
//! with trim holes punched where the tool passes through, plus the tool's band
//! faces forming the hole walls.
//!
//! Every precondition violation (open branch, partial cut, Union/Intersect,
//! NURBS↔planar splitting beyond caps, cap-face intersection) returns an
//! explicit [`crate::error::OperationError::Failed`] naming the unsupported
//! case. No silent wrong geometry, no panics.
//!
//! Pipeline:
//! 1. [`loops`] — SSI over face pairs → validated closed loops grouped per tool
//!    side face (exactly 2 each).
//! 2. [`punch`] — convert each loop's target-UV trace into a trim hole on the
//!    target face plus a 3D hole wire.
//! 3. [`band`] — build the tool's hole-wall face between each loop pair.
//! 4. [`assemble`] — collect result faces into a new shell + solid.

pub(crate) mod assemble;
pub(crate) mod band;
pub(crate) mod intersect;
pub(crate) mod loops;
pub(crate) mod pocket;
pub(crate) mod punch;

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
