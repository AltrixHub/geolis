use crate::error::{OperationError, Result};
use crate::geometry::curve::Line;
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::topology::{
    EdgeCurve, EdgeData, FaceData, FaceSurface, OrientedEdge, ShellData, SolidId, TopologyStore,
    VertexData, WireData,
};

/// Mirrors a solid across a plane defined by a point and normal.
pub struct Mirror {
    solid: SolidId,
    plane_origin: Point3,
    plane_normal: Vector3,
}

impl Mirror {
    /// Creates a new `Mirror` operation.
    #[must_use]
    pub fn new(solid: SolidId, plane_origin: Point3, plane_normal: Vector3) -> Self {
        Self {
            solid,
            plane_origin,
            plane_normal,
        }
    }

    /// Executes the mirror, creating a mirrored copy in the topology store.
    ///
    /// Mirror reflection reverses handedness, so face winding must be reversed
    /// to keep normals pointing outward.
    ///
    /// # Errors
    ///
    /// Returns an error if the plane normal is zero-length.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        let len = self.plane_normal.norm();
        if len < TOLERANCE {
            return Err(
                OperationError::InvalidInput("mirror plane normal must be non-zero".into()).into(),
            );
        }
        let normal = self.plane_normal / len;

        // Collect the source solid structure
        let solid = store.solid(self.solid)?;
        let outer_shell_id = solid.outer_shell;
        let inner_shell_ids = solid.inner_shells.clone();

        let new_outer_shell = mirror_shell(store, outer_shell_id, &self.plane_origin, &normal)?;
        let mut new_inner_shells = Vec::with_capacity(inner_shell_ids.len());
        for &shell_id in &inner_shell_ids {
            new_inner_shells.push(mirror_shell(store, shell_id, &self.plane_origin, &normal)?);
        }

        crate::operations::creation::MakeSolid::new(new_outer_shell, new_inner_shells)
            .execute(store)
    }
}

/// Reflects a point across a plane defined by origin and unit normal.
fn mirror_point(point: &Point3, plane_origin: &Point3, plane_normal: &Vector3) -> Point3 {
    let d = (point - plane_origin).dot(plane_normal);
    point - 2.0 * d * plane_normal
}

/// Mirrors an entire shell, creating new topology entities.
fn mirror_shell(
    store: &mut TopologyStore,
    shell_id: crate::topology::ShellId,
    plane_origin: &Point3,
    plane_normal: &Vector3,
) -> Result<crate::topology::ShellId> {
    let shell = store.shell(shell_id)?;
    let face_ids = shell.faces.clone();
    let is_closed = shell.is_closed;

    // Collect unique vertex and edge IDs from the shell (read-only pass)
    let (unique_verts, unique_edges) =
        collect_unique_ids(store, &face_ids)?;

    // Mirror all unique vertices
    let vertex_map = mirror_vertices(store, &unique_verts, plane_origin, plane_normal)?;

    // Mirror all unique edges
    let edge_map = mirror_edges(store, &unique_edges, &vertex_map, plane_origin, plane_normal)?;

    // Create faces with reversed winding
    build_mirrored_faces(store, &face_ids, &edge_map, is_closed)
}

/// Collects all unique vertex and edge IDs from a set of faces.
fn collect_unique_ids(
    store: &TopologyStore,
    face_ids: &[crate::topology::FaceId],
) -> Result<(
    Vec<crate::topology::VertexId>,
    Vec<crate::topology::EdgeId>,
)> {
    use std::collections::HashSet;

    let mut vert_set = HashSet::new();
    let mut edge_set = HashSet::new();

    for &face_id in face_ids {
        let face = store.face(face_id)?;
        let wire_ids: Vec<_> = std::iter::once(face.outer_wire)
            .chain(face.inner_wires.iter().copied())
            .collect();

        for wire_id in wire_ids {
            let wire = store.wire(wire_id)?;
            for oe in &wire.edges {
                edge_set.insert(oe.edge);
                let edge = store.edge(oe.edge)?;
                vert_set.insert(edge.start);
                vert_set.insert(edge.end);
            }
        }
    }

    Ok((
        vert_set.into_iter().collect(),
        edge_set.into_iter().collect(),
    ))
}

/// Creates mirrored vertex copies, returning old-to-new ID mapping.
fn mirror_vertices(
    store: &mut TopologyStore,
    vertex_ids: &[crate::topology::VertexId],
    plane_origin: &Point3,
    plane_normal: &Vector3,
) -> Result<std::collections::HashMap<crate::topology::VertexId, crate::topology::VertexId>> {
    let mut map = std::collections::HashMap::with_capacity(vertex_ids.len());
    for &vid in vertex_ids {
        let pt = store.vertex(vid)?.point;
        let mirrored = mirror_point(&pt, plane_origin, plane_normal);
        let new_id = store.add_vertex(VertexData::new(mirrored));
        map.insert(vid, new_id);
    }
    Ok(map)
}

/// Creates mirrored edge copies, returning old-to-new ID mapping.
#[allow(clippy::similar_names)]
fn mirror_edges(
    store: &mut TopologyStore,
    edge_ids: &[crate::topology::EdgeId],
    vertex_map: &std::collections::HashMap<crate::topology::VertexId, crate::topology::VertexId>,
    plane_origin: &Point3,
    plane_normal: &Vector3,
) -> Result<std::collections::HashMap<crate::topology::EdgeId, crate::topology::EdgeId>> {
    use crate::geometry::curve::Curve;

    let mut map = std::collections::HashMap::with_capacity(edge_ids.len());

    for &eid in edge_ids {
        let edge = store.edge(eid)?;
        let new_start = vertex_map[&edge.start];
        let new_end = vertex_map[&edge.end];
        let curve = edge.curve.clone();

        let start_pt = store.vertex(new_start)?.point;
        let end_pt = store.vertex(new_end)?.point;

        let new_edge_data = match &curve {
            EdgeCurve::Line(_) => {
                let dir = end_pt - start_pt;
                let t_end = dir.norm();
                let line = Line::new(start_pt, dir)?;
                EdgeData {
                    start: new_start,
                    end: new_end,
                    curve: EdgeCurve::Line(line),
                    t_start: 0.0,
                    t_end,
                }
            }
            EdgeCurve::Arc(arc) => {
                let center = mirror_point(arc.center(), plane_origin, plane_normal);
                let arc_normal = *arc.normal();
                let reflected = arc_normal - 2.0 * arc_normal.dot(plane_normal) * plane_normal;
                let mirrored_normal = -reflected;
                let to_start = start_pt - center;
                let radius = to_start.norm();
                let ref_dir = to_start / radius;
                let domain = arc.domain();

                let new_arc = crate::geometry::curve::Arc::new(
                    center,
                    radius,
                    mirrored_normal,
                    ref_dir,
                    domain.t_min,
                    domain.t_max,
                )?;
                EdgeData {
                    start: new_start,
                    end: new_end,
                    curve: EdgeCurve::Arc(new_arc),
                    t_start: domain.t_min,
                    t_end: domain.t_max,
                }
            }
        };
        let new_eid = store.add_edge(new_edge_data);
        map.insert(eid, new_eid);
    }
    Ok(map)
}

/// Builds mirrored faces with reversed winding, assembles into a shell.
fn build_mirrored_faces(
    store: &mut TopologyStore,
    face_ids: &[crate::topology::FaceId],
    edge_map: &std::collections::HashMap<crate::topology::EdgeId, crate::topology::EdgeId>,
    is_closed: bool,
) -> Result<crate::topology::ShellId> {
    let mut new_faces = Vec::with_capacity(face_ids.len());

    for &face_id in face_ids {
        let face = store.face(face_id)?;
        let outer_wire_id = face.outer_wire;
        let inner_wire_ids = face.inner_wires.clone();

        let new_outer_wire = mirror_wire_reversed(store, outer_wire_id, edge_map)?;

        let mut new_inner_wires = Vec::with_capacity(inner_wire_ids.len());
        for &iw in &inner_wire_ids {
            new_inner_wires.push(mirror_wire_reversed(store, iw, edge_map)?);
        }

        // Compute new face plane from the mirrored outer wire points
        let points = collect_wire_points(store, new_outer_wire)?;

        let plane = if points.len() >= 3 {
            compute_plane_from_points(&points)?
        } else {
            match &store.face(face_id)?.surface {
                FaceSurface::Plane(p) => p.clone(),
            }
        };

        let new_face_id = store.add_face(FaceData {
            surface: FaceSurface::Plane(plane),
            outer_wire: new_outer_wire,
            inner_wires: new_inner_wires,
            same_sense: true,
        });
        new_faces.push(new_face_id);
    }

    Ok(store.add_shell(ShellData {
        faces: new_faces,
        is_closed,
    }))
}

/// Collects vertex points from a wire in traversal order.
fn collect_wire_points(
    store: &TopologyStore,
    wire_id: crate::topology::WireId,
) -> Result<Vec<Point3>> {
    let wire = store.wire(wire_id)?;
    let edges = wire.edges.clone();
    let mut points = Vec::with_capacity(edges.len());
    for oe in &edges {
        let edge = store.edge(oe.edge)?;
        let vid = if oe.forward { edge.start } else { edge.end };
        points.push(store.vertex(vid)?.point);
    }
    Ok(points)
}

/// Creates a mirrored wire with reversed edge order and flipped orientations.
fn mirror_wire_reversed(
    store: &mut TopologyStore,
    wire_id: crate::topology::WireId,
    edge_map: &std::collections::HashMap<crate::topology::EdgeId, crate::topology::EdgeId>,
) -> Result<crate::topology::WireId> {
    let wire = store.wire(wire_id)?;
    let is_closed = wire.is_closed;
    let edges = wire.edges.clone();

    // Reverse the order and flip forward/backward to reverse winding
    let new_edges: Vec<OrientedEdge> = edges
        .iter()
        .rev()
        .map(|oe| OrientedEdge::new(edge_map[&oe.edge], !oe.forward))
        .collect();

    Ok(store.add_wire(WireData {
        edges: new_edges,
        is_closed,
    }))
}

/// Computes a plane from points using Newell's method.
fn compute_plane_from_points(points: &[Point3]) -> Result<crate::geometry::surface::Plane> {
    use crate::error::OperationError;

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

    #[test]
    fn mirror_across_yz_plane() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(1.0, 0.0, 0.0), p(2.0, 0.0, 0.0), p(2.0, 1.0, 0.0), p(1.0, 1.0, 0.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        // Mirror across YZ plane (x = 0)
        let mirrored = Mirror::new(solid, p(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0))
            .execute(&mut store)
            .unwrap();

        // Mirrored solid should have x in [-2, -1]
        let solid_data = store.solid(mirrored).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        assert_eq!(shell.faces.len(), 6);

        for &fid in &shell.faces {
            let face = store.face(fid).unwrap();
            let wire = store.wire(face.outer_wire).unwrap();
            for oe in &wire.edges {
                let edge = store.edge(oe.edge).unwrap();
                let pt = store.vertex(edge.start).unwrap().point;
                assert!(
                    pt.x >= -2.0 - 1e-6 && pt.x <= -1.0 + 1e-6,
                    "x={} out of range",
                    pt.x
                );
            }
        }
    }

    #[test]
    fn mirror_preserves_original() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(1.0, 0.0, 0.0), p(2.0, 0.0, 0.0), p(2.0, 1.0, 0.0), p(1.0, 1.0, 0.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        let _ = Mirror::new(solid, p(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0))
            .execute(&mut store)
            .unwrap();

        // Original should still be in [1, 2]
        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        for &fid in &shell.faces {
            let face = store.face(fid).unwrap();
            let wire = store.wire(face.outer_wire).unwrap();
            for oe in &wire.edges {
                let edge = store.edge(oe.edge).unwrap();
                let pt = store.vertex(edge.start).unwrap().point;
                assert!(
                    pt.x >= 1.0 - 1e-6 && pt.x <= 2.0 + 1e-6,
                    "original x={} out of range",
                    pt.x
                );
            }
        }
    }

    #[test]
    fn zero_normal_returns_error() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        let result = Mirror::new(solid, p(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0))
            .execute(&mut store);
        assert!(result.is_err());
    }
}
