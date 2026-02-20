use std::collections::{HashMap, HashSet, VecDeque};

use spade::handles::FixedFaceHandle;
use spade::{
    ConstrainedDelaunayTriangulation, InsertionError, Point2 as SpadePoint2, Triangulation,
};

use crate::error::{Result, TessellationError};
use crate::math::Point2;
use crate::topology::{FaceId, FaceSurface, TopologyStore};

use super::{TessellationParams, TriangleMesh};

/// Tessellates a face into a triangle mesh.
pub struct TessellateFace {
    face: FaceId,
    params: TessellationParams,
}

impl TessellateFace {
    /// Creates a new `TessellateFace` operation.
    #[must_use]
    pub fn new(face: FaceId, params: TessellationParams) -> Self {
        Self { face, params }
    }

    /// Executes the tessellation, returning a triangle mesh.
    ///
    /// # Errors
    ///
    /// Returns an error if the face cannot be tessellated.
    #[allow(clippy::cast_possible_truncation)]
    pub fn execute(&self, store: &TopologyStore) -> Result<TriangleMesh> {
        let _ = self.params;
        let face = store.face(self.face)?;
        let same_sense = face.same_sense;

        let plane = match &face.surface {
            FaceSurface::Plane(p) => p.clone(),
        };

        let outer_wire_id = face.outer_wire;
        let inner_wire_ids = face.inner_wires.clone();

        // Collect 3D points from wires
        let outer_3d = collect_wire_points_3d(store, outer_wire_id)?;
        let mut inner_3d_list = Vec::new();
        for &wire_id in &inner_wire_ids {
            inner_3d_list.push(collect_wire_points_3d(store, wire_id)?);
        }

        let origin = plane.origin();
        let u_dir = plane.u_dir();
        let v_dir = plane.v_dir();
        let normal = if same_sense {
            *plane.plane_normal()
        } else {
            -*plane.plane_normal()
        };

        // Project to 2D (plane UV space)
        let project = |p: &crate::math::Point3| -> SpadePoint2<f64> {
            let d = p - origin;
            SpadePoint2::new(d.dot(u_dir), d.dot(v_dir))
        };

        let outer_2d: Vec<_> = outer_3d.iter().map(&project).collect();
        let inner_2d_list: Vec<Vec<_>> = inner_3d_list
            .iter()
            .map(|pts| pts.iter().map(&project).collect())
            .collect();

        // Build CDT
        let mut cdt = ConstrainedDelaunayTriangulation::<SpadePoint2<f64>>::new();

        insert_constraint_loop(&mut cdt, &outer_2d)?;
        for inner_2d in &inner_2d_list {
            insert_constraint_loop(&mut cdt, inner_2d)?;
        }

        // Flood-fill to classify interior/exterior triangles
        let interior_faces = classify_interior_faces(&cdt);

        // Build TriangleMesh from interior triangles
        let mut mesh = TriangleMesh::default();
        let mut vertex_map: HashMap<usize, u32> = HashMap::new();

        for face_handle in cdt.inner_faces() {
            let fix = face_handle.fix();
            if !interior_faces.contains(&fix.index()) {
                continue;
            }

            let verts = face_handle.vertices();
            let mut tri_indices = [0u32; 3];

            for (i, vh) in verts.iter().enumerate() {
                let idx = vh.fix().index();
                let mesh_idx = if let Some(&existing) = vertex_map.get(&idx) {
                    existing
                } else {
                    let pos = vh.position();
                    let u = pos.x;
                    let v = pos.y;
                    let p3 = *origin + *u_dir * u + *v_dir * v;
                    let new_idx = mesh.vertices.len() as u32;
                    mesh.vertices.push(p3);
                    mesh.normals.push(normal);
                    mesh.uvs.push(Point2::new(u, v));
                    vertex_map.insert(idx, new_idx);
                    new_idx
                };
                tri_indices[i] = mesh_idx;
            }

            mesh.indices.push(tri_indices);
        }

        Ok(mesh)
    }
}

/// Collects 3D vertex positions from a wire in traversal order.
fn collect_wire_points_3d(
    store: &TopologyStore,
    wire_id: crate::topology::WireId,
) -> Result<Vec<crate::math::Point3>> {
    let edges = store.wire(wire_id)?.edges.clone();
    let mut points = Vec::with_capacity(edges.len());
    for oe in &edges {
        let edge = store.edge(oe.edge)?;
        let vertex_id = if oe.forward { edge.start } else { edge.end };
        points.push(store.vertex(vertex_id)?.point);
    }
    Ok(points)
}

/// Inserts a closed polygon as constraint edges into the CDT.
fn insert_constraint_loop(
    cdt: &mut ConstrainedDelaunayTriangulation<SpadePoint2<f64>>,
    points: &[SpadePoint2<f64>],
) -> Result<()> {
    if points.len() < 3 {
        return Err(
            TessellationError::Failed("constraint loop needs at least 3 points".into()).into(),
        );
    }

    let mut handles = Vec::with_capacity(points.len());
    for &pt in points {
        let h = cdt
            .insert(pt)
            .map_err(|e: InsertionError| TessellationError::Failed(format!("CDT insert: {e}")))?;
        handles.push(h);
    }

    for i in 0..handles.len() {
        let from = handles[i];
        let to = handles[(i + 1) % handles.len()];
        if from != to {
            cdt.add_constraint(from, to);
        }
    }

    Ok(())
}

/// Classifies which inner faces of the CDT are inside the polygon using flood-fill.
///
/// Starts from faces adjacent to the outer (infinite) face at depth 0. Each time
/// a constraint edge is crossed, depth increments. Odd depth = interior.
fn classify_interior_faces(
    cdt: &ConstrainedDelaunayTriangulation<SpadePoint2<f64>>,
) -> HashSet<usize> {
    let mut interior = HashSet::new();
    let mut depth_map: HashMap<usize, u32> = HashMap::new();
    let mut queue: VecDeque<(FixedFaceHandle<spade::handles::InnerTag>, u32)> = VecDeque::new();

    let outer_fix = cdt.outer_face().fix();

    // Seed: find inner faces adjacent to the outer face via directed edges
    for edge in cdt.directed_edges() {
        if edge.face().fix() == outer_fix {
            let rev_face = edge.rev().face();
            if let Some(inner) = rev_face.as_inner() {
                let idx = inner.fix().index();
                if depth_map.contains_key(&idx) {
                    continue;
                }
                let depth = u32::from(cdt.is_constraint_edge(edge.as_undirected().fix()));
                depth_map.insert(idx, depth);
                if depth % 2 == 1 {
                    interior.insert(idx);
                }
                queue.push_back((inner.fix(), depth));
            }
        }
    }

    // BFS flood-fill
    while let Some((face_fix, depth)) = queue.pop_front() {
        let face = cdt.face(face_fix);
        for edge in face.adjacent_edges() {
            let neighbor = edge.rev().face();
            if let Some(inner_neighbor) = neighbor.as_inner() {
                let n_idx = inner_neighbor.fix().index();
                if depth_map.contains_key(&n_idx) {
                    continue;
                }
                let new_depth = if cdt.is_constraint_edge(edge.as_undirected().fix()) {
                    depth + 1
                } else {
                    depth
                };
                depth_map.insert(n_idx, new_depth);
                if new_depth % 2 == 1 {
                    interior.insert(n_idx);
                }
                queue.push_back((inner_neighbor.fix(), new_depth));
            }
        }
    }

    interior
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::{MakeFace, MakeWire};

    fn p(x: f64, y: f64) -> Point3 {
        Point3::new(x, y, 0.0)
    }

    fn make_face_from_points(
        store: &mut crate::topology::TopologyStore,
        points: Vec<Point3>,
    ) -> FaceId {
        let wire = MakeWire::new(points, true).execute(store).unwrap();
        MakeFace::new(wire, vec![]).execute(store).unwrap()
    }

    #[test]
    fn triangle_produces_1_triangle() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_face_from_points(&mut store, vec![p(0.0, 0.0), p(4.0, 0.0), p(2.0, 3.0)]);
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert_eq!(mesh.indices.len(), 1);
        assert_eq!(mesh.vertices.len(), 3);
        assert_eq!(mesh.normals.len(), 3);
        assert_eq!(mesh.uvs.len(), 3);
    }

    #[test]
    fn square_produces_2_triangles() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_face_from_points(
            &mut store,
            vec![p(0.0, 0.0), p(4.0, 0.0), p(4.0, 4.0), p(0.0, 4.0)],
        );
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert_eq!(mesh.indices.len(), 2);
        assert_eq!(mesh.vertices.len(), 4);
    }

    #[test]
    fn l_shape_concave_tessellates() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_face_from_points(
            &mut store,
            vec![
                p(0.0, 0.0),
                p(4.0, 0.0),
                p(4.0, 2.0),
                p(2.0, 2.0),
                p(2.0, 4.0),
                p(0.0, 4.0),
            ],
        );
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        // L-shape (6 vertices, concave) → should produce 4 triangles
        assert_eq!(mesh.indices.len(), 4);
        assert_eq!(mesh.vertices.len(), 6);
    }

    #[test]
    fn face_with_hole_excludes_interior() {
        let mut store = crate::topology::TopologyStore::new();
        let outer = MakeWire::new(
            vec![p(0.0, 0.0), p(10.0, 0.0), p(10.0, 10.0), p(0.0, 10.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let inner = MakeWire::new(
            vec![p(3.0, 3.0), p(7.0, 3.0), p(7.0, 7.0), p(3.0, 7.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(outer, vec![inner])
            .execute(&mut store)
            .unwrap();
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();

        // No triangle center should be inside the hole (3..7, 3..7)
        for tri in &mesh.indices {
            let cx = (mesh.vertices[tri[0] as usize].x
                + mesh.vertices[tri[1] as usize].x
                + mesh.vertices[tri[2] as usize].x)
                / 3.0;
            let cy = (mesh.vertices[tri[0] as usize].y
                + mesh.vertices[tri[1] as usize].y
                + mesh.vertices[tri[2] as usize].y)
                / 3.0;
            let in_hole = cx > 3.0 && cx < 7.0 && cy > 3.0 && cy < 7.0;
            assert!(
                !in_hole,
                "triangle centroid ({cx}, {cy}) is inside the hole"
            );
        }
    }

    #[test]
    fn normals_match_plane_normal() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_face_from_points(
            &mut store,
            vec![p(0.0, 0.0), p(4.0, 0.0), p(4.0, 4.0), p(0.0, 4.0)],
        );
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        for n in &mesh.normals {
            assert!(n.z.abs() > 0.99, "normal z should be ~1.0, got {}", n.z);
        }
    }
}
