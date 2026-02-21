use crate::error::{OperationError, Result};
use crate::geometry::curve::Line;
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::operations::creation::{MakeFace, MakeSolid};
use crate::topology::{
    EdgeCurve, EdgeData, EdgeId, FaceId, OrientedEdge, ShellData, SolidId, TopologyStore,
    VertexData, VertexId, WireData, WireId,
};

/// Extrudes a face along a direction vector to create a solid.
pub struct Extrude {
    face: FaceId,
    direction: Vector3,
}

impl Extrude {
    /// Creates a new `Extrude` operation.
    #[must_use]
    pub fn new(face: FaceId, direction: Vector3) -> Self {
        Self { face, direction }
    }

    /// Executes the extrusion, creating the solid in the topology store.
    ///
    /// Builds a proper `BRep` solid where adjacent faces share edges via
    /// `OrientedEdge`. Supports faces with inner wires (holes): for each hole,
    /// additional side faces are created to form the hole walls.
    ///
    /// # Errors
    ///
    /// Returns [`OperationError::InvalidInput`] if the direction is zero-length.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        // Validate direction is non-zero
        if self.direction.norm() < TOLERANCE {
            return Err(
                OperationError::InvalidInput("extrude direction must be non-zero".into()).into(),
            );
        }

        let face = store.face(self.face)?;
        let outer_wire = face.outer_wire;
        let inner_wire_ids = face.inner_wires.clone();

        // Collect base points from the outer wire
        let base_points = collect_wire_points(store, outer_wire)?;

        // Compute Newell normal of the base polygon
        let normal = newell_normal(&base_points)?;

        // Ensure base_points are ordered so their Newell normal aligns with the
        // extrude direction. Then:
        //   - bottom face = reversed winding → normal ≈ -direction (outward below)
        //   - top face = same winding translated → normal ≈ +direction (outward above)
        //   - side quads naturally face outward
        let should_reverse = normal.dot(&self.direction) < 0.0;
        let base_points = if should_reverse {
            base_points.into_iter().rev().collect()
        } else {
            base_points
        };

        let n = base_points.len();

        // --- Create outer wire vertices and edges ---
        let bottom_verts: Vec<VertexId> = base_points
            .iter()
            .map(|p| store.add_vertex(VertexData::new(*p)))
            .collect();
        let top_points: Vec<Point3> = base_points.iter().map(|p| p + self.direction).collect();
        let top_verts: Vec<VertexId> = top_points
            .iter()
            .map(|p| store.add_vertex(VertexData::new(*p)))
            .collect();

        let bottom_edges = create_loop_edges(store, &bottom_verts, &base_points)?;
        let top_edges = create_loop_edges(store, &top_verts, &top_points)?;
        let vert_edges = create_vertical_edges(store, &bottom_verts, &top_verts, &base_points, &top_points)?;

        let mut all_faces = Vec::with_capacity(n + 2 + inner_wire_ids.len() * 4);

        // --- Process inner wires (holes) ---
        let (bottom_inner_wires, top_inner_wires, hole_faces) =
            extrude_inner_wires(store, &inner_wire_ids, should_reverse, &self.direction)?;
        all_faces.extend(hole_faces);

        // Bottom face: reversed winding with inner wires as holes
        let bottom_wire_edges: Vec<OrientedEdge> = (0..n)
            .rev()
            .map(|i| OrientedEdge::new(bottom_edges[i], false))
            .collect();
        let bottom_wire = create_closed_wire(store, bottom_wire_edges);
        let bottom_face = MakeFace::new(bottom_wire, bottom_inner_wires).execute(store)?;
        all_faces.push(bottom_face);

        // Top face: same winding with inner wires as holes
        let top_wire_edges: Vec<OrientedEdge> = (0..n)
            .map(|i| OrientedEdge::new(top_edges[i], true))
            .collect();
        let top_wire = create_closed_wire(store, top_wire_edges);
        let top_face = MakeFace::new(top_wire, top_inner_wires).execute(store)?;
        all_faces.push(top_face);

        // Outer side faces
        for i in 0..n {
            let j = (i + 1) % n;
            let side_wire_edges = vec![
                OrientedEdge::new(bottom_edges[i], true),
                OrientedEdge::new(vert_edges[j], true),
                OrientedEdge::new(top_edges[i], false),
                OrientedEdge::new(vert_edges[i], false),
            ];
            let side_wire = create_closed_wire(store, side_wire_edges);
            let side_face = MakeFace::new(side_wire, vec![]).execute(store)?;
            all_faces.push(side_face);
        }

        // Create shell (closed) and solid
        let shell_id = store.add_shell(ShellData {
            faces: all_faces,
            is_closed: true,
        });
        MakeSolid::new(shell_id, vec![]).execute(store)
    }
}

/// Processes inner wires (holes) for extrusion, creating hole side faces
/// and inner wires for the cap faces.
#[allow(clippy::similar_names)]
fn extrude_inner_wires(
    store: &mut TopologyStore,
    inner_wire_ids: &[WireId],
    should_reverse: bool,
    direction: &Vector3,
) -> Result<(Vec<WireId>, Vec<WireId>, Vec<FaceId>)> {
    let mut bottom_inner_wires = Vec::with_capacity(inner_wire_ids.len());
    let mut top_inner_wires = Vec::with_capacity(inner_wire_ids.len());
    let mut hole_faces = Vec::new();

    for &inner_wire_id in inner_wire_ids {
        let inner_points_raw = collect_wire_points(store, inner_wire_id)?;
        let inner_points = if should_reverse {
            inner_points_raw.into_iter().rev().collect::<Vec<_>>()
        } else {
            inner_points_raw
        };

        let m = inner_points.len();

        let ib_verts: Vec<VertexId> = inner_points
            .iter()
            .map(|p| store.add_vertex(VertexData::new(*p)))
            .collect();
        let it_points: Vec<Point3> = inner_points.iter().map(|p| p + direction).collect();
        let it_verts: Vec<VertexId> = it_points
            .iter()
            .map(|p| store.add_vertex(VertexData::new(*p)))
            .collect();

        let ib_edges = create_loop_edges(store, &ib_verts, &inner_points)?;
        let it_edges = create_loop_edges(store, &it_verts, &it_points)?;
        let iv_edges = create_vertical_edges(store, &ib_verts, &it_verts, &inner_points, &it_points)?;

        let bottom_inner_wire_edges: Vec<OrientedEdge> = (0..m)
            .map(|i| OrientedEdge::new(ib_edges[i], true))
            .collect();
        let biw = create_closed_wire(store, bottom_inner_wire_edges);
        bottom_inner_wires.push(biw);

        let top_inner_wire_edges: Vec<OrientedEdge> = (0..m)
            .rev()
            .map(|i| OrientedEdge::new(it_edges[i], false))
            .collect();
        let tiw = create_closed_wire(store, top_inner_wire_edges);
        top_inner_wires.push(tiw);

        for i in 0..m {
            let j = (i + 1) % m;
            let side_wire_edges = vec![
                OrientedEdge::new(ib_edges[i], false),
                OrientedEdge::new(iv_edges[i], true),
                OrientedEdge::new(it_edges[i], true),
                OrientedEdge::new(iv_edges[j], false),
            ];
            let side_wire = create_closed_wire(store, side_wire_edges);
            let side_face = MakeFace::new(side_wire, vec![]).execute(store)?;
            hole_faces.push(side_face);
        }
    }

    Ok((bottom_inner_wires, top_inner_wires, hole_faces))
}

/// Collects vertex positions from a wire in traversal order.
fn collect_wire_points(store: &TopologyStore, wire_id: crate::topology::WireId) -> Result<Vec<Point3>> {
    let edges = store.wire(wire_id)?.edges.clone();
    let mut points = Vec::with_capacity(edges.len());

    for oe in &edges {
        let edge = store.edge(oe.edge)?;
        let vertex_id = if oe.forward { edge.start } else { edge.end };
        let vertex = store.vertex(vertex_id)?;
        points.push(vertex.point);
    }

    Ok(points)
}

/// Computes the normal of a polygon using Newell's method.
fn newell_normal(points: &[Point3]) -> Result<Vector3> {
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
    Ok(normal / len)
}

/// Creates edges forming a closed loop from vertex/point arrays.
fn create_loop_edges(
    store: &mut TopologyStore,
    verts: &[VertexId],
    points: &[Point3],
) -> Result<Vec<EdgeId>> {
    let n = verts.len();
    let mut edges = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        edges.push(create_line_edge(store, verts[i], verts[j], points[i], points[j])?);
    }
    Ok(edges)
}

/// Creates vertical edges connecting bottom and top vertex arrays.
fn create_vertical_edges(
    store: &mut TopologyStore,
    bottom_verts: &[VertexId],
    top_verts: &[VertexId],
    bottom_points: &[Point3],
    top_points: &[Point3],
) -> Result<Vec<EdgeId>> {
    let n = bottom_verts.len();
    let mut edges = Vec::with_capacity(n);
    for i in 0..n {
        edges.push(create_line_edge(
            store,
            bottom_verts[i],
            top_verts[i],
            bottom_points[i],
            top_points[i],
        )?);
    }
    Ok(edges)
}

/// Creates a line edge between two existing vertices.
fn create_line_edge(
    store: &mut TopologyStore,
    start: VertexId,
    end: VertexId,
    start_point: Point3,
    end_point: Point3,
) -> Result<EdgeId> {
    let direction = end_point - start_point;
    let t_end = direction.norm();
    let line = Line::new(start_point, direction)?;
    Ok(store.add_edge(EdgeData {
        start,
        end,
        curve: EdgeCurve::Line(line),
        t_start: 0.0,
        t_end,
    }))
}

/// Creates a closed wire from a sequence of oriented edges.
fn create_closed_wire(store: &mut TopologyStore, edges: Vec<OrientedEdge>) -> WireId {
    store.add_wire(WireData {
        edges,
        is_closed: true,
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use std::collections::{HashMap, HashSet};

    use super::*;
    use crate::operations::creation::MakeWire;
    use crate::tessellation::{TessellateFace, TessellateSolid, TessellationParams};
    use crate::topology::FaceSurface;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn make_face(store: &mut TopologyStore, points: Vec<Point3>) -> FaceId {
        let wire = MakeWire::new(points, true).execute(store).unwrap();
        MakeFace::new(wire, vec![]).execute(store).unwrap()
    }

    // ── Unit cube ──────────────────────────────────────────────

    #[test]
    fn unit_cube_has_6_faces() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)],
        );
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        assert_eq!(shell.faces.len(), 6); // top + bottom + 4 sides
        assert!(shell.is_closed);
    }

    // ── Triangle prism ─────────────────────────────────────────

    #[test]
    fn triangle_prism_has_5_faces() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(0.0, 0.0, 0.0), p(3.0, 0.0, 0.0), p(1.5, 2.0, 0.0)],
        );
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        assert_eq!(shell.faces.len(), 5); // top + bottom + 3 sides
    }

    // ── L-shape ────────────────────────────────────────────────

    #[test]
    fn l_shape_has_8_faces() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![
                p(0.0, 0.0, 0.0), p(4.0, 0.0, 0.0), p(4.0, 2.0, 0.0),
                p(2.0, 2.0, 0.0), p(2.0, 4.0, 0.0), p(0.0, 4.0, 0.0),
            ],
        );
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        assert_eq!(shell.faces.len(), 8); // top + bottom + 6 sides
    }

    // ── Normals point outward ──────────────────────────────────

    #[test]
    fn all_face_normals_point_outward() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(0.0, 0.0, 0.0), p(2.0, 0.0, 0.0), p(2.0, 2.0, 0.0), p(0.0, 2.0, 0.0)],
        );
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();

        // Compute centroid of the solid
        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        let centroid = p(1.0, 1.0, 1.5); // center of 2x2x3 box

        for &face_id in &shell.faces {
            let face_data = store.face(face_id).unwrap();
            let FaceSurface::Plane(plane) = &face_data.surface;
            let face_normal = plane.plane_normal();
            let face_origin = plane.origin();

            // Vector from solid centroid to face origin
            let to_face = face_origin - centroid;
            // Normal should point away from centroid (same direction as to_face)
            assert!(
                face_normal.dot(&to_face) > 0.0,
                "face normal {face_normal:?} should point outward (dot with {to_face:?} was {})",
                face_normal.dot(&to_face)
            );
        }
    }

    // ── Error cases ────────────────────────────────────────────

    #[test]
    fn zero_direction_returns_error() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0)],
        );
        let result = Extrude::new(face, Vector3::new(0.0, 0.0, 0.0)).execute(&mut store);
        assert!(result.is_err());
    }

    #[test]
    fn face_with_hole_creates_solid() {
        let mut store = TopologyStore::new();
        let outer = vec![
            p(0.0, 0.0, 0.0), p(10.0, 0.0, 0.0),
            p(10.0, 10.0, 0.0), p(0.0, 10.0, 0.0),
        ];
        let inner = vec![
            p(2.0, 2.0, 0.0), p(8.0, 2.0, 0.0),
            p(8.0, 8.0, 0.0), p(2.0, 8.0, 0.0),
        ];
        let outer_wire = MakeWire::new(outer, true).execute(&mut store).unwrap();
        let inner_wire = MakeWire::new(inner, true).execute(&mut store).unwrap();
        let face = MakeFace::new(outer_wire, vec![inner_wire])
            .execute(&mut store)
            .unwrap();

        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 5.0))
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        // 2 cap faces (with holes) + 4 outer sides + 4 inner sides = 10
        assert_eq!(shell.faces.len(), 10);
        assert!(shell.is_closed);
    }

    #[test]
    fn face_with_hole_all_edges_shared() {
        let mut store = TopologyStore::new();
        let outer = vec![
            p(0.0, 0.0, 0.0), p(10.0, 0.0, 0.0),
            p(10.0, 10.0, 0.0), p(0.0, 10.0, 0.0),
        ];
        let inner = vec![
            p(3.0, 3.0, 0.0), p(7.0, 3.0, 0.0),
            p(7.0, 7.0, 0.0), p(3.0, 7.0, 0.0),
        ];
        let outer_wire = MakeWire::new(outer, true).execute(&mut store).unwrap();
        let inner_wire = MakeWire::new(inner, true).execute(&mut store).unwrap();
        let face = MakeFace::new(outer_wire, vec![inner_wire])
            .execute(&mut store)
            .unwrap();

        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        let counts = count_edge_usage(&store, shell);
        for (edge_id, count) in &counts {
            assert_eq!(
                *count, 2,
                "edge {edge_id:?} should be used exactly 2 times, got {count}"
            );
        }
    }

    // ── TessellateSolid ────────────────────────────────────────

    #[test]
    fn tessellate_cube_produces_12_triangles() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)],
        );
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();

        // 6 faces × 2 triangles each = 12 triangles
        assert_eq!(mesh.indices.len(), 12);
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
        assert_eq!(mesh.vertices.len(), mesh.uvs.len());
    }

    #[test]
    fn tessellate_prism_normals_are_unit() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(0.0, 0.0, 0.0), p(3.0, 0.0, 0.0), p(1.5, 2.0, 0.0)],
        );
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();

        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();

        for normal in &mesh.normals {
            let len = normal.norm();
            assert!(
                (len - 1.0).abs() < 1e-6,
                "normal {normal:?} has length {len}, expected 1.0"
            );
        }
    }

    #[test]
    fn tessellate_face_via_extrude_produces_correct_counts() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(0.0, 0.0, 0.0), p(3.0, 0.0, 0.0), p(1.5, 2.0, 0.0)],
        );
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();

        // Each face should tessellate independently
        for &face_id in &shell.faces {
            let mesh = TessellateFace::new(face_id, TessellationParams::default())
                .execute(&store)
                .unwrap();
            assert!(!mesh.indices.is_empty());
            assert!(!mesh.vertices.is_empty());
        }
    }

    // ── BRep shared topology helpers ─────────────────────────────

    /// Collects unique vertex and edge IDs from all faces in a shell.
    fn collect_shell_topology(
        store: &TopologyStore,
        shell: &crate::topology::ShellData,
    ) -> (HashSet<VertexId>, HashSet<EdgeId>) {
        let mut vertices = HashSet::new();
        let mut edges = HashSet::new();
        for &face_id in &shell.faces {
            let face = store.face(face_id).unwrap();
            let wire = store.wire(face.outer_wire).unwrap();
            for oe in &wire.edges {
                let edge = store.edge(oe.edge).unwrap();
                vertices.insert(edge.start);
                vertices.insert(edge.end);
                edges.insert(oe.edge);
            }
        }
        (vertices, edges)
    }

    /// Counts how many times each edge appears across all face wires.
    fn count_edge_usage(
        store: &TopologyStore,
        shell: &crate::topology::ShellData,
    ) -> HashMap<EdgeId, usize> {
        let mut counts = HashMap::new();
        for &face_id in &shell.faces {
            let face = store.face(face_id).unwrap();
            let wire = store.wire(face.outer_wire).unwrap();
            for oe in &wire.edges {
                *counts.entry(oe.edge).or_insert(0) += 1;
            }
            for &inner_wire_id in &face.inner_wires {
                let inner_wire = store.wire(inner_wire_id).unwrap();
                for oe in &inner_wire.edges {
                    *counts.entry(oe.edge).or_insert(0) += 1;
                }
            }
        }
        counts
    }

    // ── BRep shared topology verification ────────────────────────

    #[test]
    fn cube_shared_vertices() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)],
        );
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        let (vertices, _) = collect_shell_topology(&store, shell);
        assert_eq!(vertices.len(), 8, "cube should have 8 unique vertices");
    }

    #[test]
    fn cube_shared_edges() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)],
        );
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        let (_, edges) = collect_shell_topology(&store, shell);
        assert_eq!(edges.len(), 12, "cube should have 12 unique edges");
    }

    #[test]
    fn cube_each_edge_used_twice() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)],
        );
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        let counts = count_edge_usage(&store, shell);
        for (edge_id, count) in &counts {
            assert_eq!(
                *count, 2,
                "edge {edge_id:?} should be used exactly 2 times, got {count}"
            );
        }
    }

    #[test]
    fn prism_shared_topology() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(0.0, 0.0, 0.0), p(3.0, 0.0, 0.0), p(1.5, 2.0, 0.0)],
        );
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        let (vertices, edges) = collect_shell_topology(&store, shell);
        assert_eq!(vertices.len(), 6, "prism should have 6 unique vertices");
        assert_eq!(edges.len(), 9, "prism should have 9 unique edges");

        let counts = count_edge_usage(&store, shell);
        for (edge_id, count) in &counts {
            assert_eq!(
                *count, 2,
                "edge {edge_id:?} should be used exactly 2 times, got {count}"
            );
        }
    }
}
