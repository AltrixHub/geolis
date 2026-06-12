//! NURBS intersection machinery.
//!
//! Provides transversal intersections needed by trimming (P5) and solid
//! booleans (P6): 2D curve×curve, 3D curve×surface, surface×plane, and
//! surface×surface (SSI). Every result carries both 3D geometry and UV
//! (parameter-space) data.
//!
//! ## Pipeline
//!
//! All solvers share the same two-stage shape: bounding-box subdivision
//! ([`bbox`]) produces Newton seeds from the control-hull convex-hull property,
//! then a Newton/Gauss-Newton refinement (or marching, for the surface curves)
//! converges to the true intersection.
//!
//! ## Unsupported degeneracies
//!
//! Per the P4 scope, only *transversal* intersections are a quality target.
//! Tangential contacts, overlapping (coincident) regions, and self-intersecting
//! input are required only to **terminate cleanly** (no hang, no panic) — they
//! may return empty, partial, or boundary-only output. The `max_points` guard
//! in [`IntersectionOptions`] bounds every marching loop.

mod bbox;
mod curve_curve;
mod curve_surface;
mod surface_plane;
mod surface_surface;
mod types;

pub use curve_curve::intersect_curves_2d;
pub use curve_surface::intersect_curve_surface;
pub use surface_plane::intersect_surface_plane;
pub use surface_surface::intersect_surfaces;
pub use types::{
    CurveCurveIntersection2D, CurveSurfaceIntersection, IntersectionOptions,
    SurfaceIntersectionCurve,
};
