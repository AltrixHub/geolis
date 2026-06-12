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

pub(crate) mod loops;
pub(crate) mod punch;
