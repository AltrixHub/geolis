use crate::error::{Result, TessellationError};
use crate::math::{Point2, Point3, Vector3};

use super::stroke_style::StrokeStyle;
use super::TriangleMesh;

/// Maximum miter scale factor to prevent spikes at sharp angles.
const MAX_MITER_SCALE: f64 = 2.0;

/// Up direction for the flat ribbon (Z+).
const UP: Vector3 = Vector3::new(0.0, 0.0, 1.0);

/// Generates a flat ribbon triangle mesh from a polyline and stroke style.
///
/// The ribbon lies in the XY plane with normals pointing in the Z+ direction.
#[derive(Debug)]
pub struct TessellateStroke {
    points: Vec<Point3>,
    style: StrokeStyle,
    closed: bool,
}

impl TessellateStroke {
    /// Creates a new stroke tessellation operation.
    #[must_use]
    pub fn new(points: Vec<Point3>, style: StrokeStyle, closed: bool) -> Self {
        Self {
            points,
            style,
            closed,
        }
    }

    /// Executes the tessellation, producing a ribbon mesh.
    ///
    /// # Errors
    ///
    /// Returns an error if fewer than 2 points are provided, or if consecutive
    /// points are coincident (zero-length segment).
    pub fn execute(&self) -> Result<TriangleMesh> {
        let n = self.points.len();
        if n < 2 {
            return Err(TessellationError::InvalidParameters(
                "at least 2 points are required for stroke tessellation".to_owned(),
            )
            .into());
        }

        let half_w = self.style.half_width();

        // Compute tangent directions at each vertex.
        let tangents = self.compute_tangents()?;

        // Compute offset directions and miter scales.
        let offsets = self.compute_offsets(&tangents)?;

        // Generate left/right vertices, normals, and UVs.
        let vertex_count = n * 2;
        let mut vertices = Vec::with_capacity(vertex_count);
        let mut normals = Vec::with_capacity(vertex_count);
        let mut uvs = Vec::with_capacity(vertex_count);

        // Compute cumulative arc lengths for UV V-coordinate.
        let arc_lengths = self.cumulative_arc_lengths();
        let total_length = arc_lengths[n - 1];
        // Avoid division by zero for degenerate polylines.
        let inv_total = if total_length > f64::EPSILON {
            1.0 / total_length
        } else {
            0.0
        };

        let normal = UP;

        for i in 0..n {
            let (offset_dir, miter_scale) = offsets[i];
            let offset = offset_dir * half_w * miter_scale;

            // Left vertex (U=0)
            let left = Point3::new(
                self.points[i].x + offset.x,
                self.points[i].y + offset.y,
                self.points[i].z + offset.z,
            );
            // Right vertex (U=1)
            let right = Point3::new(
                self.points[i].x - offset.x,
                self.points[i].y - offset.y,
                self.points[i].z - offset.z,
            );

            let v = arc_lengths[i] * inv_total;

            vertices.push(left);
            vertices.push(right);
            normals.push(normal);
            normals.push(normal);
            uvs.push(Point2::new(0.0, v));
            uvs.push(Point2::new(1.0, v));
        }

        // Generate triangle indices.
        let segment_count = if self.closed { n } else { n - 1 };
        let mut indices = Vec::with_capacity(segment_count * 2);

        for i in 0..segment_count {
            let j = (i + 1) % n;
            #[allow(clippy::cast_possible_truncation)]
            let i0 = (i * 2) as u32;
            #[allow(clippy::cast_possible_truncation)]
            let i1 = (i * 2 + 1) as u32;
            #[allow(clippy::cast_possible_truncation)]
            let j0 = (j * 2) as u32;
            #[allow(clippy::cast_possible_truncation)]
            let j1 = (j * 2 + 1) as u32;

            // Two triangles per quad segment.
            indices.push([i0, j0, i1]);
            indices.push([i1, j0, j1]);
        }

        Ok(TriangleMesh {
            vertices,
            normals,
            uvs,
            indices,
        })
    }

    /// Computes tangent direction at each vertex.
    fn compute_tangents(&self) -> Result<Vec<Vector3>> {
        let n = self.points.len();
        let mut tangents = Vec::with_capacity(n);

        for i in 0..n {
            let tangent = if self.closed {
                let prev = if i == 0 { n - 1 } else { i - 1 };
                let next = (i + 1) % n;
                let d_prev = self.segment_direction(prev, i)?;
                let d_next = self.segment_direction(i, next)?;
                average_direction(d_prev, d_next)
            } else if i == 0 {
                self.segment_direction(0, 1)?
            } else if i == n - 1 {
                self.segment_direction(n - 2, n - 1)?
            } else {
                let d_prev = self.segment_direction(i - 1, i)?;
                let d_next = self.segment_direction(i, i + 1)?;
                average_direction(d_prev, d_next)
            };
            tangents.push(tangent);
        }

        Ok(tangents)
    }

    /// Computes the normalized direction from point `a` to point `b`.
    fn segment_direction(&self, a: usize, b: usize) -> Result<Vector3> {
        let d = self.points[b] - self.points[a];
        let len = d.norm();
        if len < f64::EPSILON {
            return Err(TessellationError::InvalidParameters(format!(
                "zero-length segment between points {a} and {b}"
            ))
            .into());
        }
        Ok(d / len)
    }

    /// Computes offset direction and miter scale at each vertex.
    fn compute_offsets(&self, tangents: &[Vector3]) -> Result<Vec<(Vector3, f64)>> {
        let n = self.points.len();
        let mut offsets = Vec::with_capacity(n);

        for (i, tangent) in tangents.iter().enumerate() {
            let offset_dir = tangent.cross(&UP);
            let offset_len = offset_dir.norm();
            if offset_len < f64::EPSILON {
                return Err(TessellationError::InvalidParameters(
                    "tangent is parallel to the up direction (Z axis)".to_owned(),
                )
                .into());
            }
            let offset_dir = offset_dir / offset_len;

            // Compute miter scale: 1 / cos(half_angle) between adjacent segments.
            let miter_scale = if self.closed || (i > 0 && i < n - 1) {
                let prev = if i == 0 { n - 1 } else { i - 1 };
                let next = (i + 1) % n;
                let d_prev = (self.points[i] - self.points[prev]).normalize();
                let d_next = (self.points[next] - self.points[i]).normalize();
                let cos_angle = d_prev.dot(&d_next);
                // miter_scale = 1 / cos(half_angle)
                // cos(half_angle) = sqrt((1 + cos_angle) / 2)
                let cos_half = f64::midpoint(1.0, cos_angle).sqrt();
                if cos_half > f64::EPSILON {
                    (1.0 / cos_half).min(MAX_MITER_SCALE)
                } else {
                    MAX_MITER_SCALE
                }
            } else {
                1.0
            };

            offsets.push((offset_dir, miter_scale));
        }

        Ok(offsets)
    }

    /// Returns cumulative arc lengths at each vertex.
    fn cumulative_arc_lengths(&self) -> Vec<f64> {
        let n = self.points.len();
        let mut lengths = Vec::with_capacity(n);
        lengths.push(0.0);
        for i in 1..n {
            let seg_len = (self.points[i] - self.points[i - 1]).norm();
            lengths.push(lengths[i - 1] + seg_len);
        }
        lengths
    }
}

/// Returns the normalized average of two direction vectors.
fn average_direction(a: Vector3, b: Vector3) -> Vector3 {
    let avg = a + b;
    let len = avg.norm();
    if len < f64::EPSILON {
        // Opposite directions: fall back to the first direction.
        a
    } else {
        avg / len
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::float_cmp)]
mod tests {
    use super::*;

    fn style(width: f64) -> StrokeStyle {
        StrokeStyle::new(width).unwrap()
    }

    #[test]
    fn straight_line_two_points() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
        ];
        let op = TessellateStroke::new(points, style(2.0), false);
        let mesh = op.execute().unwrap();

        // 2 points -> 4 vertices, 2 triangles
        assert_eq!(mesh.vertices.len(), 4);
        assert_eq!(mesh.indices.len(), 2);

        // Check width: vertices should be offset by half_width in Y.
        // tangent=(1,0,0) x UP=(0,0,1) = (0,-1,0), so offset points -Y.
        let v0 = &mesh.vertices[0]; // point + offset => y = -1
        let v1 = &mesh.vertices[1]; // point - offset => y = +1
        let spread = (v1.y - v0.y).abs();
        assert!((spread - 2.0).abs() < 1e-10, "total width should be 2.0");
    }

    #[test]
    fn l_shape_three_points() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(5.0, 5.0, 0.0),
        ];
        let op = TessellateStroke::new(points, style(1.0), false);
        let mesh = op.execute().unwrap();

        // 3 points -> 6 vertices, 4 triangles
        assert_eq!(mesh.vertices.len(), 6);
        assert_eq!(mesh.indices.len(), 4);
    }

    #[test]
    fn closed_polyline() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
        ];
        let op = TessellateStroke::new(points, style(1.0), true);
        let mesh = op.execute().unwrap();

        // 3 points -> 6 vertices, 6 triangles (3 segments x 2)
        assert_eq!(mesh.vertices.len(), 6);
        assert_eq!(mesh.indices.len(), 6);
    }

    #[test]
    fn normals_point_up() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
        ];
        let op = TessellateStroke::new(points, style(1.0), false);
        let mesh = op.execute().unwrap();

        for normal in &mesh.normals {
            assert!((normal.z - 1.0).abs() < 1e-10);
            assert!(normal.x.abs() < 1e-10);
            assert!(normal.y.abs() < 1e-10);
        }
    }

    #[test]
    fn uvs_range() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
        ];
        let op = TessellateStroke::new(points, style(1.0), false);
        let mesh = op.execute().unwrap();

        // First pair: V=0
        assert!((mesh.uvs[0].y).abs() < 1e-10);
        assert!((mesh.uvs[1].y).abs() < 1e-10);
        // Last pair: V=1
        assert!((mesh.uvs[4].y - 1.0).abs() < 1e-10);
        assert!((mesh.uvs[5].y - 1.0).abs() < 1e-10);
        // U: left=0, right=1
        assert!((mesh.uvs[0].x).abs() < 1e-10);
        assert!((mesh.uvs[1].x - 1.0).abs() < 1e-10);
    }

    #[test]
    fn too_few_points_fails() {
        let points = vec![Point3::new(0.0, 0.0, 0.0)];
        let op = TessellateStroke::new(points, style(1.0), false);
        assert!(op.execute().is_err());
    }

    #[test]
    fn zero_length_segment_fails() {
        let points = vec![
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
        ];
        let op = TessellateStroke::new(points, style(1.0), false);
        assert!(op.execute().is_err());
    }
}
