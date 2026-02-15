mod tessellate_curve;
mod tessellate_face;
mod tessellate_solid;
mod tessellate_with_holes;

pub use tessellate_curve::TessellateCurve;
pub use tessellate_face::TessellateFace;
pub use tessellate_solid::TessellateSolid;
pub use tessellate_with_holes::TessellateWithHoles;

use crate::math::{Point2, Point3, Vector3};

/// Parameters controlling tessellation quality.
#[derive(Debug, Clone, Copy)]
pub struct TessellationParams {
    /// Maximum allowed deviation from the true geometry.
    pub tolerance: f64,
    /// Minimum number of segments for curves.
    pub min_segments: usize,
    /// Maximum number of segments for curves.
    pub max_segments: usize,
}

impl Default for TessellationParams {
    fn default() -> Self {
        Self {
            tolerance: 0.01,
            min_segments: 4,
            max_segments: 256,
        }
    }
}

/// A polyline approximation of a curve.
#[derive(Debug, Clone, Default)]
pub struct Polyline {
    /// The ordered vertices of the polyline.
    pub points: Vec<Point3>,
}

/// A triangle mesh approximation of a surface.
#[derive(Debug, Clone, Default)]
pub struct TriangleMesh {
    /// Vertex positions.
    pub vertices: Vec<Point3>,
    /// Vertex normals.
    pub normals: Vec<Vector3>,
    /// UV coordinates.
    pub uvs: Vec<Point2>,
    /// Triangle indices (each triple defines a triangle).
    pub indices: Vec<[u32; 3]>,
}
