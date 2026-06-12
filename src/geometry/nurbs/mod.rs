//! NURBS curves and surfaces.
//!
//! Algorithms follow Piegl & Tiller, *The NURBS Book* (2nd ed.); function
//! documentation cites algorithm numbers (e.g. A2.2, A4.3).

mod basis;
mod curve;
mod knot;
mod surface;

mod construct;
mod intersect;

pub use basis::{basis_function_derivatives, basis_functions, binomial};
pub use curve::{NurbsCurve, NurbsCurve2D, NurbsCurve3D};
pub use intersect::{
    intersect_curve_surface, intersect_curves_2d, CurveCurveIntersection2D,
    CurveSurfaceIntersection, IntersectionOptions, SurfaceIntersectionCurve,
};
pub use knot::KnotVector;
pub use surface::{InversionOptions, NurbsSurface, SurfaceInversion};
