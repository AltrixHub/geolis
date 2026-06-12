//! Result and option types shared by the NURBS intersection solvers.

use crate::math::{Point2, Point3};

/// A transversal intersection point between two 2D curves.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveCurveIntersection2D {
    /// The intersection point in 2D space.
    pub point: Point2,
    /// Parameter on curve `a`.
    pub t_a: f64,
    /// Parameter on curve `b`.
    pub t_b: f64,
}

/// An intersection point between a 3D curve and a surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveSurfaceIntersection {
    /// The intersection point in 3D space.
    pub point: Point3,
    /// Curve parameter.
    pub t: f64,
    /// Surface u parameter.
    pub u: f64,
    /// Surface v parameter.
    pub v: f64,
}

/// One branch of a surface-surface (or surface-plane) intersection: a polyline
/// in 3D with synchronized UV traces on both surfaces.
///
/// Curve fitting of these polylines into NURBS pcurves happens in P5. For a
/// surface-plane intersection `uv_b` holds the plane's own 2D coordinates in an
/// orthonormal in-plane basis.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceIntersectionCurve {
    /// 3D points along the branch (ordered).
    pub points: Vec<Point3>,
    /// UV trace on surface `a`, one entry per point.
    pub uv_a: Vec<Point2>,
    /// UV trace on surface `b` (or plane coordinates), one entry per point.
    pub uv_b: Vec<Point2>,
    /// Whether the branch forms a closed loop.
    pub closed: bool,
}

/// Numerical options shared by the intersection solvers.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionOptions {
    /// Geometric coincidence tolerance.
    pub tolerance: f64,
    /// Maximum Newton refinement iterations.
    pub max_iterations: usize,
    /// Marching step size as a fraction of the local feature scale.
    pub step_factor: f64,
    /// Marching runaway guard: maximum points emitted per branch.
    pub max_points: usize,
}

impl Default for IntersectionOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-9,
            max_iterations: 64,
            step_factor: 0.5,
            max_points: 10_000,
        }
    }
}
