use std::f64::consts::TAU;

use crate::error::{OperationError, Result};
use crate::geometry::curve::{Circle, Line};
use crate::geometry::surface::{Cone, Cylinder, Plane};
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::operations::creation::MakeSolid;
use crate::topology::{
    EdgeCurve, EdgeData, EdgeId, FaceData, FaceId, FaceSurface, OrientedEdge, ShellData, SolidId,
    TopologyStore, VertexData, VertexId, WireData, WireId,
};

/// Revolves a planar face 360 degrees around an axis to create a solid of revolution.
///
/// The profile face must be planar with line edges only (polygon profile).
/// Full revolution only (no partial angles).
pub struct Revolve {
    face: FaceId,
    axis_origin: Point3,
    axis_dir: Vector3,
}

impl Revolve {
    /// Creates a new `Revolve` operation.
    #[must_use]
    pub fn new(face: FaceId, axis_origin: Point3, axis_dir: Vector3) -> Self {
        Self {
            face,
            axis_origin,
            axis_dir,
        }
    }

    /// Executes the revolution, creating a solid in the topology store.
    ///
    /// For each profile edge, a side face is generated:
    /// - Edges parallel to the axis produce Cylinder faces
    /// - Edges at an angle to the axis produce Cone faces
    /// - Vertices on the axis are degenerate (zero-radius circle)
    ///
    /// # Errors
    ///
    /// Returns an error if the axis direction is zero-length, the face doesn't exist,
    /// or any profile edge is not a line.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        let axis_len = self.axis_dir.norm();
        if axis_len < TOLERANCE {
            return Err(
                OperationError::InvalidInput("revolve axis direction must be non-zero".into())
                    .into(),
            );
        }
        let axis = self.axis_dir / axis_len;

        let face = store.face(self.face)?;
        let outer_wire_id = face.outer_wire;

        // Collect profile vertices in order
        let profile_points = collect_wire_points(store, outer_wire_id)?;
        let n = profile_points.len();
        if n < 3 {
            return Err(OperationError::InvalidInput(
                "revolve profile must have at least 3 vertices".into(),
            )
            .into());
        }

        // Compute per-vertex distance from axis and axis-projected height
        let vert_info: Vec<VertexInfo> = profile_points
            .iter()
            .map(|p| compute_vertex_info(p, &self.axis_origin, &axis))
            .collect();

        // Compute a stable reference direction for circles/surfaces.
        // Use the first vertex that is NOT on the axis.
        let ref_dir = compute_ref_dir(&vert_info, &axis)?;

        // Create topology vertices (full revolution: start = end, so one vertex per profile point)
        let verts: Vec<VertexId> = profile_points
            .iter()
            .map(|p| store.add_vertex(VertexData::new(*p)))
            .collect();

        // Create circle edges for each vertex (full revolution, start == end)
        let circle_edges: Vec<Option<EdgeId>> = vert_info
            .iter()
            .zip(&verts)
            .map(|(info, &vid)| {
                if info.radius < TOLERANCE {
                    // Vertex is on the axis: degenerate, no circle edge
                    None
                } else {
                    let circle = make_circle_on_axis(
                        &info.axis_foot,
                        info.radius,
                        &axis,
                        &ref_dir,
                    );
                    match circle {
                        Ok(c) => Some(store.add_edge(EdgeData {
                            start: vid,
                            end: vid,
                            curve: EdgeCurve::Circle(c),
                            t_start: 0.0,
                            t_end: TAU,
                        })),
                        Err(_) => None,
                    }
                }
            })
            .collect();

        // Create seam line edges (connecting profile vertices along the axis seam).
        // Full revolution: each seam edge is shared by two adjacent side faces.
        let seam_edges: Vec<EdgeId> = (0..n)
            .map(|i| {
                let j = (i + 1) % n;
                create_line_edge(store, verts[i], verts[j], profile_points[i], profile_points[j])
            })
            .collect::<Result<Vec<_>>>()?;

        // Create side faces
        let mut all_faces = Vec::with_capacity(n);

        for i in 0..n {
            let j = (i + 1) % n;
            let face_id = create_side_face(
                store,
                &vert_info[i],
                &vert_info[j],
                verts[i],
                verts[j],
                circle_edges[i],
                circle_edges[j],
                seam_edges[i],
                &self.axis_origin,
                &axis,
                &ref_dir,
            )?;
            all_faces.push(face_id);
        }

        // Create shell and solid
        let shell_id = store.add_shell(ShellData {
            faces: all_faces,
            is_closed: true,
        });
        MakeSolid::new(shell_id, vec![]).execute(store)
    }
}

/// Per-vertex geometric information relative to the revolution axis.
struct VertexInfo {
    /// The original 3D point.
    point: Point3,
    /// Distance from the revolution axis.
    radius: f64,
    /// Height (signed projection onto the axis).
    height: f64,
    /// Foot of perpendicular on the axis (center for the circle).
    axis_foot: Point3,
}

fn compute_vertex_info(point: &Point3, axis_origin: &Point3, axis: &Vector3) -> VertexInfo {
    let dp = point - axis_origin;
    let height = dp.dot(axis);
    let axis_foot = axis_origin + axis * height;
    let radial = point - axis_foot;
    let radius = radial.norm();
    VertexInfo {
        point: *point,
        radius,
        height,
        axis_foot,
    }
}

/// Finds a reference direction perpendicular to the axis, pointing towards the first
/// off-axis profile vertex.
fn compute_ref_dir(vert_info: &[VertexInfo], _axis: &Vector3) -> Result<Vector3> {
    for vi in vert_info {
        if vi.radius > TOLERANCE {
            let radial = vi.point - vi.axis_foot;
            return Ok(radial / vi.radius);
        }
    }
    // All vertices on axis: degenerate profile
    Err(OperationError::InvalidInput("all profile vertices lie on the revolution axis".into())
        .into())
}

/// Creates a Circle curve on the axis plane at the given foot point.
fn make_circle_on_axis(
    center: &Point3,
    radius: f64,
    axis: &Vector3,
    ref_dir: &Vector3,
) -> Result<Circle> {
    // The ref_dir may not be exactly perpendicular to the axis at this center
    // (it's a global ref_dir), but it IS perpendicular to axis since it was
    // computed as a radial direction. We need to re-derive it for this center.
    Circle::new(*center, radius, *axis, *ref_dir)
}

/// Creates a line edge between two vertices.
fn create_line_edge(
    store: &mut TopologyStore,
    start: VertexId,
    end: VertexId,
    start_pt: Point3,
    end_pt: Point3,
) -> Result<EdgeId> {
    let direction = end_pt - start_pt;
    let t_end = direction.norm();
    let line = Line::new(start_pt, direction)?;
    Ok(store.add_edge(EdgeData {
        start,
        end,
        curve: EdgeCurve::Line(line),
        t_start: 0.0,
        t_end,
    }))
}

/// Creates a closed wire from oriented edges.
fn create_closed_wire(store: &mut TopologyStore, edges: Vec<OrientedEdge>) -> WireId {
    store.add_wire(WireData {
        edges,
        is_closed: true,
    })
}

/// Creates a side face for one profile edge revolved 360 degrees.
#[allow(clippy::too_many_arguments)]
fn create_side_face(
    store: &mut TopologyStore,
    vi: &VertexInfo,
    vj: &VertexInfo,
    _vid_i: VertexId,
    _vid_j: VertexId,
    circle_i: Option<EdgeId>,
    circle_j: Option<EdgeId>,
    seam_edge: EdgeId,
    axis_origin: &Point3,
    axis: &Vector3,
    ref_dir: &Vector3,
) -> Result<FaceId> {
    // Determine the surface type based on the radii
    let on_axis_i = vi.radius < TOLERANCE;
    let on_axis_j = vj.radius < TOLERANCE;

    // Build the wire for this side face
    let wire_edges = match (on_axis_i, on_axis_j, circle_i, circle_j) {
        // Both vertices off-axis: normal case
        (false, false, Some(ci), Some(cj)) => {
            vec![
                OrientedEdge::new(ci, true),    // circle at i (forward)
                OrientedEdge::new(seam_edge, false), // seam edge (reverse = return seam)
                OrientedEdge::new(cj, false),   // circle at j (reverse)
                OrientedEdge::new(seam_edge, true),  // seam edge (forward)
            ]
        }
        // i on axis, j off-axis: cone tip at i
        (true, false, None, Some(cj)) => {
            vec![
                OrientedEdge::new(seam_edge, true),  // seam from i to j
                OrientedEdge::new(cj, true),         // circle at j
                OrientedEdge::new(seam_edge, false), // seam from j back to i
            ]
        }
        // i off-axis, j on axis: cone tip at j
        (false, true, Some(ci), None) => {
            vec![
                OrientedEdge::new(ci, true),         // circle at i
                OrientedEdge::new(seam_edge, true),  // seam from i to j
                OrientedEdge::new(seam_edge, false), // seam from j back to i
            ]
        }
        // Both on axis: degenerate (zero-area face), skip
        _ => {
            return Err(OperationError::Failed(
                "both vertices on axis: degenerate face".into(),
            )
            .into());
        }
    };

    let wire = create_closed_wire(store, wire_edges);

    // Determine surface
    let surface = compute_side_surface(vi, vj, on_axis_i, on_axis_j, axis_origin, axis, ref_dir)?;

    // Determine same_sense: the surface normal should point outward.
    // For a CCW profile (looking from outside towards axis), the outward normal
    // of the surface should agree with the surface's natural normal.
    let same_sense = determine_same_sense(vi, vj, axis, &surface);

    Ok(store.add_face(FaceData {
        surface,
        outer_wire: wire,
        inner_wires: vec![],
        same_sense,
    }))
}

/// Computes the surface for a side face.
fn compute_side_surface(
    vi: &VertexInfo,
    vj: &VertexInfo,
    on_axis_i: bool,
    on_axis_j: bool,
    axis_origin: &Point3,
    axis: &Vector3,
    ref_dir: &Vector3,
) -> Result<FaceSurface> {
    let r_i = vi.radius;
    let r_j = vj.radius;
    let h_i = vi.height;
    let h_j = vj.height;

    // Check if the edge is parallel to the axis (same radius => cylinder)
    if !on_axis_i && !on_axis_j && (r_i - r_j).abs() < TOLERANCE {
        let cyl = Cylinder::new(*axis_origin, r_i, *axis, *ref_dir)?;
        Ok(FaceSurface::Cylinder(cyl))
    } else if (h_j - h_i).abs() < TOLERANCE {
        // Edge is perpendicular to axis (same height, different radii).
        // Revolution produces a flat annular disc → Plane surface.
        let origin = axis_origin + axis * h_i;
        let plane = Plane::from_normal(origin, *axis)?;
        Ok(FaceSurface::Plane(plane))
    } else if on_axis_i {
        // Apex at i, cone opens toward j
        let half_angle = (r_j / (h_j - h_i).abs()).atan();
        let cone_axis = if h_j > h_i { *axis } else { -*axis };
        let cone = Cone::new(vi.point, cone_axis, half_angle, *ref_dir)?;
        Ok(FaceSurface::Cone(cone))
    } else if on_axis_j {
        // Apex at j, cone opens toward i
        let half_angle = (r_i / (h_i - h_j).abs()).atan();
        let cone_axis = if h_i > h_j { *axis } else { -*axis };
        let cone = Cone::new(vj.point, cone_axis, half_angle, *ref_dir)?;
        Ok(FaceSurface::Cone(cone))
    } else {
        // General cone: find apex by intersecting the generator line with the axis
        let dh = h_j - h_i;
        let dr = r_j - r_i;
        let t_apex = -r_i / dr;
        let apex_height = h_i + t_apex * dh;
        let apex = axis_origin + axis * apex_height;
        let half_angle = (r_i / (h_i - apex_height).abs()).atan();
        let cone_axis = if h_i > apex_height { *axis } else { -*axis };
        let cone = Cone::new(apex, cone_axis, half_angle, *ref_dir)?;
        Ok(FaceSurface::Cone(cone))
    }
}

/// Determines whether the face normal should agree with the surface normal.
///
/// For a revolution around the axis, outward normals point away from the axis.
/// The surface's natural normal direction depends on its parametrization.
fn determine_same_sense(
    vi: &VertexInfo,
    vj: &VertexInfo,
    axis: &Vector3,
    surface: &FaceSurface,
) -> bool {
    // Use the midpoint of the profile edge to test
    let mid = Point3::new(
        f64::midpoint(vi.point.x, vj.point.x),
        f64::midpoint(vi.point.y, vj.point.y),
        f64::midpoint(vi.point.z, vj.point.z),
    );

    // Profile edge tangent (from i to j)
    let edge_dir = vj.point - vi.point;
    let edge_len = edge_dir.norm();
    if edge_len < TOLERANCE {
        return true;
    }
    let edge_tangent = edge_dir / edge_len;

    // Outward direction = edge_tangent × axis (for CCW profile viewed from outside)
    // Actually: for a profile in a plane containing the axis, the outward normal
    // in the revolution is: (axis × edge_tangent) or -(axis × edge_tangent).
    // We need to check which points away from the axis.
    let mid_foot_h = (mid - Point3::origin()).dot(axis);
    let mid_foot = Point3::origin() + axis * mid_foot_h;
    let radial = mid - mid_foot;
    let radial_len = radial.norm();

    // Cross product of axis and edge tangent gives a direction
    let cross = axis.cross(&edge_tangent);
    let outward_agrees = if radial_len > TOLERANCE {
        cross.dot(&(radial / radial_len)) > 0.0
    } else {
        true
    };

    // For Cylinder/Cone, the natural outward normal points radially out.
    // So same_sense = outward_agrees (the outward direction matches surface normal).
    match surface {
        FaceSurface::Cylinder(_) | FaceSurface::Cone(_) => outward_agrees,
        FaceSurface::Plane(_) => {
            // For an annular disc face: determine if the plane normal should
            // point along or against the axis. Use cross product of axis and
            // the profile edge direction to find the outward-pointing sense.
            let edge_dir = vj.point - vi.point;
            let cross = axis.cross(&edge_dir);
            let mid_foot = Point3::origin() + axis * f64::midpoint(vi.height, vj.height);
            let mid_pt = Point3::new(
                f64::midpoint(vi.point.x, vj.point.x),
                f64::midpoint(vi.point.y, vj.point.y),
                f64::midpoint(vi.point.z, vj.point.z),
            );
            let radial_dir = mid_pt - mid_foot;
            let rl = radial_dir.norm();
            if rl > TOLERANCE {
                cross.dot(&(radial_dir / rl)) > 0.0
            } else {
                true
            }
        }
        _ => true,
    }
}

/// Collects vertex positions from a wire in traversal order.
fn collect_wire_points(
    store: &TopologyStore,
    wire_id: WireId,
) -> Result<Vec<Point3>> {
    let edges = store.wire(wire_id)?.edges.clone();
    let mut points = Vec::with_capacity(edges.len());
    for oe in &edges {
        let edge = store.edge(oe.edge)?;
        let vertex_id = if oe.forward { edge.start } else { edge.end };
        points.push(store.vertex(vertex_id)?.point);
    }
    Ok(points)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::creation::{MakeFace, MakeWire};
    use crate::tessellation::{TessellateFace, TessellateSolid, TessellationParams};
    use std::collections::{HashMap, HashSet};

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn make_face(store: &mut TopologyStore, points: Vec<Point3>) -> FaceId {
        let wire = MakeWire::new(points, true).execute(store).unwrap();
        MakeFace::new(wire, vec![]).execute(store).unwrap()
    }

    // ── Square profile → cylinder solid ────────────────────────

    #[test]
    fn square_revolve_has_4_faces() {
        let mut store = TopologyStore::new();
        // Square profile in the XZ plane, offset from axis
        let face = make_face(
            &mut store,
            vec![
                p(2.0, 0.0, 0.0),
                p(4.0, 0.0, 0.0),
                p(4.0, 0.0, 3.0),
                p(2.0, 0.0, 3.0),
            ],
        );
        let solid = Revolve::new(face, Point3::origin(), Vector3::z())
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        // 4 profile edges → 4 side faces
        assert_eq!(shell.faces.len(), 4);
        assert!(shell.is_closed);
    }

    #[test]
    fn square_revolve_surfaces_are_correct() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![
                p(2.0, 0.0, 0.0),
                p(4.0, 0.0, 0.0),
                p(4.0, 0.0, 3.0),
                p(2.0, 0.0, 3.0),
            ],
        );
        let solid = Revolve::new(face, Point3::origin(), Vector3::z())
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();

        let mut cyl_count = 0;
        let mut plane_count = 0;
        for &fid in &shell.faces {
            let f = store.face(fid).unwrap();
            match &f.surface {
                FaceSurface::Cylinder(_) => cyl_count += 1,
                FaceSurface::Plane(_) => plane_count += 1,
                FaceSurface::Cone(_) => {}
                _ => panic!("unexpected surface type"),
            }
        }
        // 2 vertical edges (parallel to axis) → Cylinder
        // 2 horizontal edges (perpendicular to axis) → Plane (annular disc)
        assert_eq!(cyl_count, 2, "expected 2 cylinder faces, got {cyl_count}");
        assert_eq!(plane_count, 2, "expected 2 plane faces, got {plane_count}");
    }

    // ── Triangle with vertex on axis → cone faces ──────────────

    #[test]
    fn triangle_on_axis_has_3_faces() {
        let mut store = TopologyStore::new();
        // Triangle with one vertex on the Z axis
        let face = make_face(
            &mut store,
            vec![
                p(0.0, 0.0, 5.0),  // on axis
                p(3.0, 0.0, 0.0),
                p(3.0, 0.0, 5.0),
            ],
        );
        let solid = Revolve::new(face, Point3::origin(), Vector3::z())
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        assert_eq!(shell.faces.len(), 3);
    }

    #[test]
    fn triangle_on_axis_has_cone_face() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![
                p(0.0, 0.0, 5.0),
                p(3.0, 0.0, 0.0),
                p(3.0, 0.0, 5.0),
            ],
        );
        let solid = Revolve::new(face, Point3::origin(), Vector3::z())
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();

        let mut has_cone = false;
        for &fid in &shell.faces {
            let f = store.face(fid).unwrap();
            if matches!(&f.surface, FaceSurface::Cone(_)) {
                has_cone = true;
            }
        }
        assert!(has_cone, "expected at least one cone face");
    }

    // ── Edge sharing ───────────────────────────────────────────

    #[test]
    fn square_revolve_edges_shared() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![
                p(2.0, 0.0, 0.0),
                p(4.0, 0.0, 0.0),
                p(4.0, 0.0, 3.0),
                p(2.0, 0.0, 3.0),
            ],
        );
        let solid = Revolve::new(face, Point3::origin(), Vector3::z())
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

    // ── Tessellation ───────────────────────────────────────────

    #[test]
    fn revolve_tessellates() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![
                p(2.0, 0.0, 0.0),
                p(4.0, 0.0, 0.0),
                p(4.0, 0.0, 3.0),
                p(2.0, 0.0, 3.0),
            ],
        );
        let solid = Revolve::new(face, Point3::origin(), Vector3::z())
            .execute(&mut store)
            .unwrap();

        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
    }

    #[test]
    fn revolve_each_face_tessellates() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![
                p(2.0, 0.0, 0.0),
                p(4.0, 0.0, 0.0),
                p(4.0, 0.0, 3.0),
                p(2.0, 0.0, 3.0),
            ],
        );
        let solid = Revolve::new(face, Point3::origin(), Vector3::z())
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        for &fid in &shell.faces {
            let mesh = TessellateFace::new(fid, TessellationParams::default())
                .execute(&store)
                .unwrap();
            assert!(!mesh.indices.is_empty(), "face {fid:?} produced empty mesh");
        }
    }

    #[test]
    fn revolve_solid_has_many_triangles() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![
                p(2.0, 0.0, 0.0),
                p(4.0, 0.0, 0.0),
                p(4.0, 0.0, 3.0),
                p(2.0, 0.0, 3.0),
            ],
        );
        let solid = Revolve::new(face, Point3::origin(), Vector3::z())
            .execute(&mut store)
            .unwrap();

        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();

        // 4 faces (2 cylinder + 2 plane), each should produce triangles.
        // Cylinder faces with full TAU sweep produce many triangles from the UV grid.
        assert!(
            mesh.indices.len() > 50,
            "expected many triangles for revolved solid, got {}",
            mesh.indices.len()
        );

        // Verify each face individually produces a reasonable number of triangles
        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        for &fid in &shell.faces {
            let face_mesh = TessellateFace::new(fid, TessellationParams::default())
                .execute(&store)
                .unwrap();
            assert!(
                face_mesh.indices.len() >= 2,
                "face {fid:?} has only {} triangles",
                face_mesh.indices.len()
            );
        }
    }

    // ── Error cases ────────────────────────────────────────────

    #[test]
    fn zero_axis_returns_error() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![p(2.0, 0.0, 0.0), p(4.0, 0.0, 0.0), p(3.0, 0.0, 3.0)],
        );
        let result = Revolve::new(face, Point3::origin(), Vector3::new(0.0, 0.0, 0.0))
            .execute(&mut store);
        assert!(result.is_err());
    }

    // ── Helpers ────────────────────────────────────────────────

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
        }
        counts
    }

    fn _collect_shell_topology(
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
}
