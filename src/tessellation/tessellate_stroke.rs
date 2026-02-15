use crate::error::{Result, TessellationError};
use crate::math::{Point2, Point3, Vector3};

use super::stroke_style::{LineJoin, StrokeStyle};
use super::TriangleMesh;

/// When the miter scale exceeds this limit, switch to a bevel join.
const BEVEL_THRESHOLD: f64 = 2.0;

/// Up direction for the flat ribbon (Z+).
const UP: Vector3 = Vector3::new(0.0, 0.0, 1.0);

/// How a polyline vertex maps to mesh vertices at the join.
enum JoinKind {
    /// Miter join or endpoint — 2 mesh vertices (left, right).
    Miter { dir: Vector3, scale: f64 },
    /// Bevel join — 3 mesh vertices (1 shared inside + 2 split outside).
    Bevel {
        inside_dir: Vector3,
        inside_scale: f64,
        /// `true` when the inside of the bend is the right (−offset) side.
        inside_is_right: bool,
        outside_in_dir: Vector3,
        outside_out_dir: Vector3,
    },
}

/// Mesh vertex indices associated with a single polyline vertex.
struct VertexSlot {
    /// Left index for connecting to the *incoming* segment.
    in_left: u32,
    /// Right index for connecting to the *incoming* segment.
    in_right: u32,
    /// Left index for connecting to the *outgoing* segment.
    out_left: u32,
    /// Right index for connecting to the *outgoing* segment.
    out_right: u32,
}

/// Generates a flat ribbon triangle mesh from a polyline and stroke style.
///
/// The ribbon lies in the XY plane with normals pointing in the Z+ direction.
/// At sharp angles (miter scale > [`BEVEL_THRESHOLD`]) a bevel join is used
/// instead of a miter to prevent spikes.
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
    #[allow(clippy::cast_possible_truncation, clippy::too_many_lines)]
    pub fn execute(&self) -> Result<TriangleMesh> {
        let n = self.points.len();
        if n < 2 {
            return Err(TessellationError::InvalidParameters(
                "at least 2 points are required for stroke tessellation".to_owned(),
            )
            .into());
        }

        let half_w = self.style.half_width();
        let joins = self.compute_joins()?;

        let mut vertices = Vec::new();
        let mut normals = Vec::new();
        let mut uvs = Vec::new();
        let mut slots = Vec::with_capacity(n);
        let mut bevel_tris: Vec<[u32; 3]> = Vec::new();

        // Compute cumulative arc lengths for UV V-coordinate.
        let arc_lengths = self.cumulative_arc_lengths();
        let total_length = arc_lengths[n - 1];
        let inv_total = if total_length > f64::EPSILON {
            1.0 / total_length
        } else {
            0.0
        };

        let normal = UP;

        for i in 0..n {
            let p = &self.points[i];
            let v = arc_lengths[i] * inv_total;

            match &joins[i] {
                JoinKind::Miter { dir, scale } => {
                    let idx = vertices.len() as u32;
                    let off = *dir * half_w * *scale;
                    vertices.push(Point3::new(p.x + off.x, p.y + off.y, p.z + off.z));
                    vertices.push(Point3::new(p.x - off.x, p.y - off.y, p.z - off.z));
                    normals.extend_from_slice(&[normal, normal]);
                    uvs.extend_from_slice(&[
                        Point2::new(0.0, v),
                        Point2::new(1.0, v),
                    ]);
                    slots.push(VertexSlot {
                        in_left: idx,
                        in_right: idx + 1,
                        out_left: idx,
                        out_right: idx + 1,
                    });
                }
                JoinKind::Bevel {
                    inside_dir,
                    inside_scale,
                    inside_is_right,
                    outside_in_dir,
                    outside_out_dir,
                } => {
                    let idx = vertices.len() as u32;
                    let in_off = *inside_dir * half_w * *inside_scale;
                    let off_in = *outside_in_dir * half_w;
                    let off_out = *outside_out_dir * half_w;

                    if *inside_is_right {
                        // Left turn: inside = right (−offset), outside = left (+offset).
                        let v_out_in = Point3::new(p.x + off_in.x, p.y + off_in.y, p.z + off_in.z);
                        let v_in = Point3::new(p.x - in_off.x, p.y - in_off.y, p.z - in_off.z);
                        let v_out_out = Point3::new(p.x + off_out.x, p.y + off_out.y, p.z + off_out.z);

                        vertices.extend_from_slice(&[v_out_in, v_in, v_out_out]);
                        normals.extend_from_slice(&[normal, normal, normal]);
                        uvs.extend_from_slice(&[
                            Point2::new(0.0, v),
                            Point2::new(1.0, v),
                            Point2::new(0.0, v),
                        ]);

                        slots.push(VertexSlot {
                            in_left: idx,
                            in_right: idx + 1,
                            out_left: idx + 2,
                            out_right: idx + 1,
                        });
                        // Bevel triangle: inside → outside_in → outside_out (CCW).
                        bevel_tris.push([idx + 1, idx, idx + 2]);
                    } else {
                        // Right turn: inside = left (+offset), outside = right (−offset).
                        let v_in = Point3::new(p.x + in_off.x, p.y + in_off.y, p.z + in_off.z);
                        let v_out_in = Point3::new(p.x - off_in.x, p.y - off_in.y, p.z - off_in.z);
                        let v_out_out = Point3::new(p.x - off_out.x, p.y - off_out.y, p.z - off_out.z);

                        vertices.extend_from_slice(&[v_in, v_out_in, v_out_out]);
                        normals.extend_from_slice(&[normal, normal, normal]);
                        uvs.extend_from_slice(&[
                            Point2::new(0.0, v),
                            Point2::new(1.0, v),
                            Point2::new(1.0, v),
                        ]);

                        slots.push(VertexSlot {
                            in_left: idx,
                            in_right: idx + 1,
                            out_left: idx,
                            out_right: idx + 2,
                        });
                        // Bevel triangle: inside → outside_out → outside_in (CCW).
                        bevel_tris.push([idx, idx + 2, idx + 1]);
                    }
                }
            }
        }

        // Generate segment quads.
        let segment_count = if self.closed { n } else { n - 1 };
        let mut indices = Vec::with_capacity(segment_count * 2 + bevel_tris.len());

        for i in 0..segment_count {
            let j = (i + 1) % n;
            let si = &slots[i];
            let sj = &slots[j];
            indices.push([si.out_left, sj.in_left, si.out_right]);
            indices.push([si.out_right, sj.in_left, sj.in_right]);
        }

        // Append bevel triangles.
        indices.extend_from_slice(&bevel_tris);

        Ok(TriangleMesh {
            vertices,
            normals,
            uvs,
            indices,
        })
    }

    /// Determines the join kind (miter or bevel) at each polyline vertex.
    fn compute_joins(&self) -> Result<Vec<JoinKind>> {
        let n = self.points.len();
        let line_join = self.style.line_join();
        let mut joins = Vec::with_capacity(n);

        for i in 0..n {
            let is_interior = self.closed || (i > 0 && i < n - 1);

            if is_interior {
                let prev = if i == 0 { n - 1 } else { i - 1 };
                let next = (i + 1) % n;
                let d_prev = self.segment_direction(prev, i)?;
                let d_next = self.segment_direction(i, next)?;

                let cos_angle = d_prev.dot(&d_next);
                let cos_half = f64::midpoint(1.0, cos_angle).sqrt();
                let miter_scale = if cos_half > f64::EPSILON {
                    1.0 / cos_half
                } else {
                    f64::MAX
                };

                let use_bevel = match line_join {
                    LineJoin::Miter => false,
                    LineJoin::Bevel => true,
                    LineJoin::Auto => miter_scale > BEVEL_THRESHOLD,
                };

                if use_bevel {
                    let turn_z = d_prev.cross(&d_next).z;
                    let inside_is_right = turn_z > 0.0;

                    let tangent = average_direction(d_prev, d_next);
                    let inside_dir = normalize_perp(tangent)?;
                    let inside_scale = miter_scale.min(BEVEL_THRESHOLD);

                    let outside_in_dir = normalize_perp(d_prev)?;
                    let outside_out_dir = normalize_perp(d_next)?;

                    joins.push(JoinKind::Bevel {
                        inside_dir,
                        inside_scale,
                        inside_is_right,
                        outside_in_dir,
                        outside_out_dir,
                    });
                } else {
                    let tangent = average_direction(d_prev, d_next);
                    let perp = normalize_perp(tangent)?;
                    joins.push(JoinKind::Miter {
                        dir: perp,
                        scale: miter_scale.min(BEVEL_THRESHOLD),
                    });
                }
            } else {
                // Endpoint: perpendicular to the single adjacent segment.
                let seg_dir = if i == 0 {
                    self.segment_direction(0, 1)?
                } else {
                    self.segment_direction(n - 2, n - 1)?
                };
                let perp = normalize_perp(seg_dir)?;
                joins.push(JoinKind::Miter {
                    dir: perp,
                    scale: 1.0,
                });
            }
        }

        Ok(joins)
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

/// Returns the normalized perpendicular of `dir` in the XY plane (`dir × UP`).
fn normalize_perp(dir: Vector3) -> Result<Vector3> {
    let perp = dir.cross(&UP);
    let len = perp.norm();
    if len < f64::EPSILON {
        return Err(TessellationError::InvalidParameters(
            "tangent is parallel to the up direction (Z axis)".to_owned(),
        )
        .into());
    }
    Ok(perp / len)
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

        // 2 endpoints -> 4 vertices, 2 triangles
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
        // 90° turn: miter_scale = 1.414 < BEVEL_THRESHOLD → miter join
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(5.0, 5.0, 0.0),
        ];
        let op = TessellateStroke::new(points, style(1.0), false);
        let mesh = op.execute().unwrap();

        // 2 endpoints + 1 miter -> 6 vertices, 4 triangles
        assert_eq!(mesh.vertices.len(), 6);
        assert_eq!(mesh.indices.len(), 4);
    }

    #[test]
    fn closed_triangle_uses_bevel() {
        // Closed right triangle: two vertices exceed BEVEL_THRESHOLD.
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
        ];
        let op = TessellateStroke::new(points, style(1.0), true);
        let mesh = op.execute().unwrap();

        // Vertex 0 (0,0): bevel → 3, Vertex 1 (10,0): miter → 2,
        // Vertex 2 (10,10): bevel → 3. Total: 8 vertices.
        assert_eq!(mesh.vertices.len(), 8);
        // 3 segment quads (6) + 2 bevel triangles = 8.
        assert_eq!(mesh.indices.len(), 8);
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
    fn hairpin_uses_bevel() {
        // Near-reversal: miter_scale >> BEVEL_THRESHOLD → bevel
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 4.0, 0.0),
            Point3::new(0.2, 0.2, 0.0),
        ];
        let op = TessellateStroke::new(points, style(0.5), false);
        let mesh = op.execute().unwrap();

        // 2 endpoints (2 each) + 1 bevel (3) = 7 vertices
        assert_eq!(mesh.vertices.len(), 7);
        // 2 segment quads (4) + 1 bevel triangle = 5
        assert_eq!(mesh.indices.len(), 5);
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
