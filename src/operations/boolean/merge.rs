use std::collections::{HashMap, HashSet, VecDeque};

use crate::error::{OperationError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::operations::creation::{MakeFace, MakeSolid, MakeWire};
use crate::topology::{FaceId, FaceSurface, ShellData, SolidId, TopologyStore};

/// Tolerance for grouping coplanar faces (plane distance comparison).
const COPLANAR_DISTANCE_TOL: f64 = 1e-6;

/// Tolerance for grouping coplanar faces (normal direction comparison).
const COPLANAR_NORMAL_TOL: f64 = 1e-6;

/// Info extracted from one face for the merge algorithm.
struct FaceInfo {
    face_id: FaceId,
    /// Directed edges as (`start_point`, `end_point`) pairs.
    edges: Vec<(Point3, Point3)>,
    /// Directed edges from inner wires, grouped by wire.
    inner_edges: Vec<Vec<(Point3, Point3)>>,
    /// Effective outward normal of the face.
    normal: Vector3,
    /// Signed distance from origin to the face plane along the normal.
    plane_dist: f64,
}

/// Merges coplanar adjacent faces in a solid after boolean operations.
///
/// Groups faces that share the same plane, finds connected components
/// (faces sharing edges), and merges each component into a single face
/// with potential inner wires (holes).
pub fn merge_coplanar_faces(store: &mut TopologyStore, solid_id: SolidId) -> Result<SolidId> {
    let face_infos = collect_all_face_info(store, solid_id)?;

    if face_infos.is_empty() {
        return Ok(solid_id);
    }

    let groups = group_coplanar(&face_infos);

    let mut merged_face_ids: Vec<FaceId> = Vec::new();
    let mut consumed: HashSet<usize> = HashSet::new();

    for group in &groups {
        if group.len() < 2 {
            continue;
        }

        let components = find_connected_components(group, &face_infos);

        for component in &components {
            if component.len() < 2 {
                continue;
            }

            // The planar-arrangement merge returns multiple `FaceId`s
            // when the union of the component is a disjoint set of
            // regions on the plane (e.g. two separated wall sections
            // sharing a coplanar group). Empty result = "skip this
            // component, keep originals" so cancellation / classification
            // failures don't strip legitimate geometry.
            let new_faces = merge_component(store, component, &face_infos)?;
            if !new_faces.is_empty() {
                merged_face_ids.extend(new_faces);
                for &idx in component {
                    consumed.insert(idx);
                }
            }
        }
    }

    // If nothing was merged, return the original solid
    if consumed.is_empty() {
        return Ok(solid_id);
    }

    // Collect non-merged faces
    for (i, info) in face_infos.iter().enumerate() {
        if !consumed.contains(&i) {
            merged_face_ids.push(info.face_id);
        }
    }

    // Build new shell and solid
    let shell_id = store.add_shell(ShellData {
        faces: merged_face_ids,
        is_closed: true,
    });
    MakeSolid::new(shell_id, vec![]).execute(store)
}

/// Collects face info for all faces of a solid.
fn collect_all_face_info(store: &TopologyStore, solid_id: SolidId) -> Result<Vec<FaceInfo>> {
    let solid = store.solid(solid_id)?;
    let shell = store.shell(solid.outer_shell)?;
    let face_ids = shell.faces.clone();

    let mut infos = Vec::with_capacity(face_ids.len());
    for face_id in face_ids {
        infos.push(collect_face_info(store, face_id)?);
    }
    Ok(infos)
}

/// Extracts directed edges and plane info from a single face.
fn collect_face_info(store: &TopologyStore, face_id: FaceId) -> Result<FaceInfo> {
    let face = store.face(face_id)?;
    let FaceSurface::Plane(ref plane) = face.surface else {
        if matches!(face.surface, FaceSurface::Nurbs(_)) {
            return Err(OperationError::Failed(
                "boolean operations on NURBS faces are not yet supported".into(),
            )
            .into());
        }
        todo!("Coplanar face merge for non-planar faces")
    };
    let wire = store.wire(face.outer_wire)?;

    let mut edges = Vec::with_capacity(wire.edges.len());
    for oe in &wire.edges {
        let edge = store.edge(oe.edge)?;
        let (start_vid, end_vid) = if oe.forward {
            (edge.start, edge.end)
        } else {
            (edge.end, edge.start)
        };
        let start_pt = store.vertex(start_vid)?.point;
        let end_pt = store.vertex(end_vid)?.point;
        edges.push((start_pt, end_pt));
    }

    // Collect inner wire edges
    let mut inner_edges = Vec::with_capacity(face.inner_wires.len());
    for &inner_wire_id in &face.inner_wires {
        let inner_wire = store.wire(inner_wire_id)?;
        let mut wire_edges = Vec::with_capacity(inner_wire.edges.len());
        for oe in &inner_wire.edges {
            let edge = store.edge(oe.edge)?;
            let (start_vid, end_vid) = if oe.forward {
                (edge.start, edge.end)
            } else {
                (edge.end, edge.start)
            };
            let start_pt = store.vertex(start_vid)?.point;
            let end_pt = store.vertex(end_vid)?.point;
            wire_edges.push((start_pt, end_pt));
        }
        inner_edges.push(wire_edges);
    }

    // Effective normal: if same_sense is false, flip the normal
    let surface_normal = *plane.plane_normal();
    let normal = if face.same_sense {
        surface_normal
    } else {
        -surface_normal
    };

    // Signed distance from origin to the plane along the normal
    let plane_dist = plane.origin().coords.dot(&normal);

    Ok(FaceInfo {
        face_id,
        edges,
        inner_edges,
        normal,
        plane_dist,
    })
}

/// Groups face indices by coplanarity (same normal direction AND same plane distance).
fn group_coplanar(face_infos: &[FaceInfo]) -> Vec<Vec<usize>> {
    let n = face_infos.len();
    let mut visited = vec![false; n];
    let mut groups: Vec<Vec<usize>> = Vec::new();

    for i in 0..n {
        if visited[i] {
            continue;
        }
        visited[i] = true;
        let mut group = vec![i];

        for j in (i + 1)..n {
            if visited[j] {
                continue;
            }
            if are_coplanar(&face_infos[i], &face_infos[j]) {
                visited[j] = true;
                group.push(j);
            }
        }

        groups.push(group);
    }

    groups
}

/// Checks if two faces are coplanar (same normal direction and same plane distance).
fn are_coplanar(a: &FaceInfo, b: &FaceInfo) -> bool {
    // Normals must point in the same direction
    let dot = a.normal.dot(&b.normal);
    if (dot - 1.0).abs() > COPLANAR_NORMAL_TOL {
        return false;
    }

    // Same signed distance from origin
    (a.plane_dist - b.plane_dist).abs() < COPLANAR_DISTANCE_TOL
}

/// Finds connected components within a group of coplanar faces.
///
/// Two faces are adjacent if one has directed edge (A→B) and the other has (B→A).
fn find_connected_components(group: &[usize], face_infos: &[FaceInfo]) -> Vec<Vec<usize>> {
    let n = group.len();
    if n <= 1 {
        return vec![group.to_vec()];
    }

    // Build adjacency: for each directed edge, record which group-local index owns it
    let mut edge_to_face: HashMap<EdgeKey, Vec<usize>> = HashMap::new();
    for (local_idx, &global_idx) in group.iter().enumerate() {
        for (start, end) in &face_infos[global_idx].edges {
            let key = EdgeKey::new(start, end);
            edge_to_face.entry(key).or_default().push(local_idx);
        }
        for inner_wire_edges in &face_infos[global_idx].inner_edges {
            for (start, end) in inner_wire_edges {
                let key = EdgeKey::new(start, end);
                edge_to_face.entry(key).or_default().push(local_idx);
            }
        }
    }

    // Build adjacency graph: face i is adjacent to face j if they share a reversed edge
    let mut adj: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for (local_idx, &global_idx) in group.iter().enumerate() {
        let info = &face_infos[global_idx];
        let outer_iter = info.edges.iter();
        let inner_iter = info.inner_edges.iter().flat_map(|w| w.iter());
        for &(start, end) in outer_iter.chain(inner_iter) {
            // Look for the reverse edge (end→start)
            let reverse_key = EdgeKey::new(&end, &start);
            if let Some(neighbors) = edge_to_face.get(&reverse_key) {
                for &neighbor in neighbors {
                    if neighbor != local_idx {
                        adj[local_idx].insert(neighbor);
                        adj[neighbor].insert(local_idx);
                    }
                }
            }
        }
    }

    // BFS to find connected components
    let mut visited = vec![false; n];
    let mut components: Vec<Vec<usize>> = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut component = vec![group[start]]; // store global indices
        let mut queue = VecDeque::new();
        queue.push_back(start);

        while let Some(curr) = queue.pop_front() {
            for &neighbor in &adj[curr] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    component.push(group[neighbor]);
                    queue.push_back(neighbor);
                }
            }
        }

        components.push(component);
    }

    components
}

/// Merges a connected component of coplanar faces into one or more
/// faces via the shared 2D planar arrangement engine.
///
/// The previous implementation cancelled exact reverse-matched edges
/// and walked the remainder into closed loops with a "pick max-area as
/// outer, rest as inner" heuristic. That produced PSLG-invalid faces
/// whenever the cancellation left T-junctions in place (e.g. a door
/// reveal whose vertical edges touch the wall's bottom outline) —
/// `chain_into_loops` chose an arbitrary outgoing edge at the junction
/// and split a single notched outer into two touching loops, which
/// later panicked `spade::cdt`.
///
/// The new implementation projects every input face fragment to the
/// component's plane UV, hands the polygon-with-holes set to
/// `boolean_2d::engine::run_arrangement` with `UnionOracle`, and
/// converts each resulting union polygon back to a BRep face. The 2D
/// engine handles segment splitting at T-junctions / proper crossings,
/// vertex snapping, bilateral half-edge classification, and face
/// walking — so the BRep faces it produces are guaranteed simple in
/// the plane and CDT-safe (a debug post-condition in the engine
/// asserts the latter).
///
/// May return multiple `FaceId`s when the union of the component
/// happens to be a disjoint set of regions on the plane. Returns
/// `Ok(Vec::new())` to signal "skip this component, keep originals" —
/// either because the component is empty / degenerate, or because the
/// arrangement engine declined the input.
fn merge_component(
    store: &mut TopologyStore,
    component: &[usize],
    face_infos: &[FaceInfo],
) -> Result<Vec<FaceId>> {
    use crate::operations::boolean_2d::{
        run_arrangement, Polygon as B2dPolygon, PolygonWithHoles as B2dPwh, UnionOracle,
    };

    if component.is_empty() {
        return Ok(Vec::new());
    }

    // Plane frame for projection.
    let normal_raw = face_infos[component[0]].normal;
    let n_unit = {
        let norm = normal_raw.norm();
        if norm <= TOLERANCE {
            return Ok(Vec::new());
        }
        normal_raw / norm
    };
    let seed = if n_unit.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let u_axis = n_unit.cross(&seed).normalize();
    let v_axis = n_unit.cross(&u_axis);
    let origin = match face_infos[component[0]].edges.first() {
        Some(&(start, _)) => start,
        None => return Ok(Vec::new()),
    };

    let project = |p: &Point3| -> (f64, f64) {
        let d = *p - origin;
        (d.dot(&u_axis), d.dot(&v_axis))
    };
    let unproject = |uv: (f64, f64)| -> Point3 { origin + u_axis * uv.0 + v_axis * uv.1 };

    // Build PolygonWithHoles input for the arrangement engine: one PWH
    // per input face, outer polygon from the outer wire's start points,
    // holes from each inner wire's start points.
    let mut input_pwhs: Vec<B2dPwh> = Vec::with_capacity(component.len());
    for &idx in component {
        let info = &face_infos[idx];
        let outer: B2dPolygon = info.edges.iter().map(|(s, _)| project(s)).collect();
        if outer.len() < 3 {
            continue;
        }
        let holes: Vec<B2dPolygon> = info
            .inner_edges
            .iter()
            .filter_map(|edges| {
                let ring: B2dPolygon = edges.iter().map(|(s, _)| project(s)).collect();
                (ring.len() >= 3).then_some(ring)
            })
            .collect();
        input_pwhs.push(B2dPwh { outer, holes });
    }
    if input_pwhs.is_empty() {
        return Ok(Vec::new());
    }

    // Run the planar arrangement union. If the engine rejects the
    // input (extreme degeneracy after 3 ε-shrink retries), leave the
    // component unmerged.
    let union_pwhs = {
        let oracle = UnionOracle {
            inputs: &input_pwhs,
        };
        match run_arrangement(&input_pwhs, &oracle) {
            Ok(r) => r,
            Err(_) => return Ok(Vec::new()),
        }
    };

    // Lift each union PWH back to 3D and materialise it as a BRep face.
    let mut face_ids = Vec::with_capacity(union_pwhs.len());
    for pwh in union_pwhs {
        let outer_3d: Vec<Point3> = pwh.outer.iter().copied().map(unproject).collect();
        if outer_3d.len() < 3 {
            continue;
        }
        let outer_wire = MakeWire::new(outer_3d, true).execute(store)?;
        let mut inner_wires = Vec::with_capacity(pwh.holes.len());
        for hole in pwh.holes {
            let hole_3d: Vec<Point3> = hole.iter().copied().map(unproject).collect();
            if hole_3d.len() < 3 {
                continue;
            }
            inner_wires.push(MakeWire::new(hole_3d, true).execute(store)?);
        }
        let face_id = MakeFace::new(outer_wire, inner_wires).execute(store)?;
        face_ids.push(face_id);
    }

    Ok(face_ids)
}

/// Quantised directed-edge key used by `find_connected_components` to
/// match adjacent face fragments. Uses a 1-micron grid — fine enough to
/// distinguish independently-placed vertices, coarse enough to merge
/// `f64` rounding noise from the boolean engine's output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct EdgeKey {
    start_key: (i64, i64, i64),
    end_key: (i64, i64, i64),
}

impl EdgeKey {
    fn new(start: &Point3, end: &Point3) -> Self {
        Self {
            start_key: quantize(start),
            end_key: quantize(end),
        }
    }
}

#[expect(
    clippy::cast_possible_truncation,
    reason = "values are millimetre-scale Point3 coordinates × 1e6 = single-micron grid; never exceed i64 range in practice"
)]
fn quantize(p: &Point3) -> (i64, i64, i64) {
    const INV_GRID: f64 = 1e6;
    (
        (p.x * INV_GRID).round() as i64,
        (p.y * INV_GRID).round() as i64,
        (p.z * INV_GRID).round() as i64,
    )
}
