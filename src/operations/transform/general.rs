use std::collections::HashSet;

use crate::error::Result;
use crate::geometry::curve::Line;
use crate::math::{Matrix4, Point3, Vector3};
use crate::topology::{EdgeCurve, SolidId, TopologyStore, VertexId};

/// Applies an arbitrary 4x4 transformation matrix to a solid.
pub struct GeneralTransform {
    solid: SolidId,
    matrix: Matrix4,
}

impl GeneralTransform {
    /// Creates a new `GeneralTransform` operation.
    #[must_use]
    pub fn new(solid: SolidId, matrix: Matrix4) -> Self {
        Self { solid, matrix }
    }

    /// Executes the transformation, modifying the solid in-place.
    ///
    /// Transforms all vertex positions using the 4x4 matrix, then rebuilds
    /// each edge's curve from the updated vertex positions.
    ///
    /// # Errors
    ///
    /// Returns an error if any topology entity is missing or curve
    /// reconstruction fails (e.g. zero-length edge after transform).
    pub fn execute(&self, store: &mut TopologyStore) -> Result<()> {
        // Collect all unique vertex IDs from the solid
        let vertex_ids = collect_solid_vertices(store, self.solid)?;

        // Transform all vertices
        for &vid in &vertex_ids {
            let vertex = store.vertex_mut(vid)?;
            vertex.point = transform_point(&self.matrix, &vertex.point);
        }

        // Collect all unique edge IDs and rebuild curves
        let edge_ids = collect_solid_edges(store, self.solid)?;
        for edge_id in edge_ids {
            let edge = store.edge(edge_id)?;
            let start_point = store.vertex(edge.start)?.point;
            let end_point = store.vertex(edge.end)?.point;

            let new_curve = match &edge.curve {
                EdgeCurve::Line(_) => {
                    let direction = end_point - start_point;
                    let t_end = direction.norm();
                    let line = Line::new(start_point, direction)?;
                    (EdgeCurve::Line(line), 0.0, t_end)
                }
                EdgeCurve::Arc(arc) => {
                    // Transform the arc's geometric properties
                    let center = transform_point(&self.matrix, arc.center());
                    let normal = transform_direction(&self.matrix, arc.normal());
                    let normal_len = normal.norm();
                    let normal = normal / normal_len;

                    // Recompute ref_dir perpendicular to transformed normal
                    let to_start = start_point - center;
                    let to_start_len = to_start.norm();
                    let ref_dir = to_start / to_start_len;
                    let radius = to_start_len;

                    let domain = {
                        use crate::geometry::curve::Curve;
                        arc.domain()
                    };

                    let new_arc = crate::geometry::curve::Arc::new(
                        center,
                        radius,
                        normal,
                        ref_dir,
                        domain.t_min,
                        domain.t_max,
                    )?;
                    (EdgeCurve::Arc(new_arc), domain.t_min, domain.t_max)
                }
            };

            let edge = store.edge_mut(edge_id)?;
            edge.curve = new_curve.0;
            edge.t_start = new_curve.1;
            edge.t_end = new_curve.2;
        }

        // Rebuild face surfaces from updated vertex positions
        let face_ids = collect_solid_faces(store, self.solid)?;
        for face_id in face_ids {
            let face = store.face(face_id)?;
            let outer_wire_id = face.outer_wire;
            let wire = store.wire(outer_wire_id)?;

            // Collect points from the outer wire
            let mut points = Vec::with_capacity(wire.edges.len());
            for oe in &wire.edges.clone() {
                let edge = store.edge(oe.edge)?;
                let vid = if oe.forward { edge.start } else { edge.end };
                points.push(store.vertex(vid)?.point);
            }

            if points.len() >= 3 {
                let plane = compute_plane_from_points(&points)?;
                let face = store.face_mut(face_id)?;
                face.surface = crate::topology::FaceSurface::Plane(plane);
            }
        }

        Ok(())
    }
}

/// Transforms a point by a 4x4 matrix (homogeneous coordinates).
fn transform_point(matrix: &Matrix4, point: &Point3) -> Point3 {
    let v = matrix * nalgebra::Vector4::new(point.x, point.y, point.z, 1.0);
    Point3::new(v.x, v.y, v.z)
}

/// Transforms a direction vector by a 4x4 matrix (ignoring translation).
fn transform_direction(matrix: &Matrix4, dir: &Vector3) -> Vector3 {
    let v = matrix * nalgebra::Vector4::new(dir.x, dir.y, dir.z, 0.0);
    Vector3::new(v.x, v.y, v.z)
}

/// Collects all unique vertex IDs referenced by a solid.
fn collect_solid_vertices(
    store: &TopologyStore,
    solid_id: SolidId,
) -> Result<HashSet<VertexId>> {
    let mut vertices = HashSet::new();
    let edge_ids = collect_solid_edges(store, solid_id)?;
    for edge_id in edge_ids {
        let edge = store.edge(edge_id)?;
        vertices.insert(edge.start);
        vertices.insert(edge.end);
    }
    Ok(vertices)
}

/// Collects all unique edge IDs referenced by a solid.
fn collect_solid_edges(
    store: &TopologyStore,
    solid_id: SolidId,
) -> Result<Vec<crate::topology::EdgeId>> {
    let mut edges = HashSet::new();
    let face_ids = collect_solid_faces(store, solid_id)?;
    for face_id in face_ids {
        let face = store.face(face_id)?;
        let wire = store.wire(face.outer_wire)?;
        for oe in &wire.edges {
            edges.insert(oe.edge);
        }
        for &inner_wire_id in &face.inner_wires {
            let inner_wire = store.wire(inner_wire_id)?;
            for oe in &inner_wire.edges {
                edges.insert(oe.edge);
            }
        }
    }
    Ok(edges.into_iter().collect())
}

/// Collects all face IDs from a solid (outer + inner shells).
fn collect_solid_faces(
    store: &TopologyStore,
    solid_id: SolidId,
) -> Result<Vec<crate::topology::FaceId>> {
    let solid = store.solid(solid_id)?;
    let outer_shell = solid.outer_shell;
    let inner_shells = solid.inner_shells.clone();

    let mut faces = Vec::new();
    let shell = store.shell(outer_shell)?;
    faces.extend_from_slice(&shell.faces);
    for &inner_shell_id in &inner_shells {
        let shell = store.shell(inner_shell_id)?;
        faces.extend_from_slice(&shell.faces);
    }
    Ok(faces)
}

/// Computes a plane from a set of points using Newell's method.
fn compute_plane_from_points(points: &[Point3]) -> Result<crate::geometry::surface::Plane> {
    use crate::error::OperationError;
    use crate::math::TOLERANCE;

    let n = points.len();
    let mut normal = Vector3::new(0.0, 0.0, 0.0);
    for i in 0..n {
        let curr = &points[i];
        let next = &points[(i + 1) % n];
        normal.x += (curr.y - next.y) * (curr.z + next.z);
        normal.y += (curr.z - next.z) * (curr.x + next.x);
        normal.z += (curr.x - next.x) * (curr.y + next.y);
    }

    let len = normal.norm();
    if len < TOLERANCE {
        return Err(
            OperationError::Failed("degenerate polygon: cannot compute normal".into()).into(),
        );
    }

    #[allow(clippy::cast_precision_loss)]
    let inv_n = 1.0 / n as f64;
    let centroid = Point3::new(
        points.iter().map(|p| p.x).sum::<f64>() * inv_n,
        points.iter().map(|p| p.y).sum::<f64>() * inv_n,
        points.iter().map(|p| p.z).sum::<f64>() * inv_n,
    );

    crate::geometry::surface::Plane::from_normal(centroid, normal)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Vector3;
    use crate::operations::creation::{MakeFace, MakeWire};
    use crate::operations::shaping::Extrude;
    use crate::topology::TopologyStore;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn make_unit_cube(store: &mut TopologyStore) -> SolidId {
        let wire = MakeWire::new(
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)],
            true,
        )
        .execute(store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(store).unwrap();
        Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(store)
            .unwrap()
    }

    #[test]
    fn identity_transform_preserves_vertices() {
        let mut store = TopologyStore::new();
        let solid = make_unit_cube(&mut store);

        let matrix = Matrix4::identity();
        GeneralTransform::new(solid, matrix)
            .execute(&mut store)
            .unwrap();

        let verts = collect_solid_vertices(&store, solid).unwrap();
        for vid in verts {
            let pt = store.vertex(vid).unwrap().point;
            assert!(pt.x >= -1e-10 && pt.x <= 1.0 + 1e-10);
            assert!(pt.y >= -1e-10 && pt.y <= 1.0 + 1e-10);
            assert!(pt.z >= -1e-10 && pt.z <= 1.0 + 1e-10);
        }
    }

    #[test]
    fn translation_shifts_all_vertices() {
        let mut store = TopologyStore::new();
        let solid = make_unit_cube(&mut store);

        let mut matrix = Matrix4::identity();
        matrix[(0, 3)] = 5.0;
        matrix[(1, 3)] = 3.0;
        matrix[(2, 3)] = 2.0;
        GeneralTransform::new(solid, matrix)
            .execute(&mut store)
            .unwrap();

        let verts = collect_solid_vertices(&store, solid).unwrap();
        for vid in verts {
            let pt = store.vertex(vid).unwrap().point;
            assert!(pt.x >= 5.0 - 1e-10 && pt.x <= 6.0 + 1e-10);
            assert!(pt.y >= 3.0 - 1e-10 && pt.y <= 4.0 + 1e-10);
            assert!(pt.z >= 2.0 - 1e-10 && pt.z <= 3.0 + 1e-10);
        }
    }

    #[test]
    fn uniform_scale_doubles_size() {
        let mut store = TopologyStore::new();
        let solid = make_unit_cube(&mut store);

        let matrix = Matrix4::new_scaling(2.0);
        GeneralTransform::new(solid, matrix)
            .execute(&mut store)
            .unwrap();

        let verts = collect_solid_vertices(&store, solid).unwrap();
        for vid in verts {
            let pt = store.vertex(vid).unwrap().point;
            assert!(pt.x >= -1e-10 && pt.x <= 2.0 + 1e-10);
            assert!(pt.y >= -1e-10 && pt.y <= 2.0 + 1e-10);
            assert!(pt.z >= -1e-10 && pt.z <= 2.0 + 1e-10);
        }
    }
}
