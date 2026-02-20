mod stroke_style;
mod tessellate_curve;
mod tessellate_face;
mod tessellate_solid;
mod tessellate_stroke;
mod tessellate_with_holes;

pub use stroke_style::{LineJoin, StrokeStyle};
pub use tessellate_curve::TessellateCurve;
pub use tessellate_face::TessellateFace;
pub use tessellate_solid::TessellateSolid;
pub use tessellate_stroke::TessellateStroke;
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

impl TriangleMesh {
    /// Merges another mesh into this one, offsetting indices appropriately.
    #[allow(clippy::cast_possible_truncation)]
    pub fn merge(&mut self, other: &Self) {
        let offset = self.vertices.len() as u32;
        self.vertices.extend_from_slice(&other.vertices);
        self.normals.extend_from_slice(&other.normals);
        self.uvs.extend_from_slice(&other.uvs);
        for tri in &other.indices {
            self.indices
                .push([tri[0] + offset, tri[1] + offset, tri[2] + offset]);
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn make_triangle_mesh(base_vertex: f64, base_index: u32) -> TriangleMesh {
        TriangleMesh {
            vertices: vec![
                Point3::new(base_vertex, 0.0, 0.0),
                Point3::new(base_vertex + 1.0, 0.0, 0.0),
                Point3::new(base_vertex, 1.0, 0.0),
            ],
            normals: vec![
                Vector3::new(0.0, 0.0, 1.0),
                Vector3::new(0.0, 0.0, 1.0),
                Vector3::new(0.0, 0.0, 1.0),
            ],
            uvs: vec![
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 0.0),
                Point2::new(0.0, 1.0),
            ],
            indices: vec![[base_index, base_index + 1, base_index + 2]],
        }
    }

    #[test]
    fn merge_offsets_indices() {
        let mut a = make_triangle_mesh(0.0, 0);
        let b = make_triangle_mesh(2.0, 0);
        a.merge(&b);

        assert_eq!(a.vertices.len(), 6);
        assert_eq!(a.normals.len(), 6);
        assert_eq!(a.uvs.len(), 6);
        assert_eq!(a.indices.len(), 2);
        assert_eq!(a.indices[0], [0, 1, 2]);
        assert_eq!(a.indices[1], [3, 4, 5]); // offset by 3
    }

    #[test]
    fn merge_into_empty() {
        let mut a = TriangleMesh::default();
        let b = make_triangle_mesh(0.0, 0);
        a.merge(&b);

        assert_eq!(a.vertices.len(), 3);
        assert_eq!(a.indices.len(), 1);
        assert_eq!(a.indices[0], [0, 1, 2]); // offset 0
    }

    #[test]
    fn merge_empty_into_existing() {
        let mut a = make_triangle_mesh(0.0, 0);
        let b = TriangleMesh::default();
        a.merge(&b);

        assert_eq!(a.vertices.len(), 3);
        assert_eq!(a.indices.len(), 1);
    }
}
