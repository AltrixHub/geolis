use std::collections::{HashMap, HashSet, VecDeque};

use crate::error::{OperationError, Result};
use crate::math::polygon_3d::polygon_area_3d;
use crate::math::{Point3, Vector3};
use crate::operations::creation::{MakeFace, MakeSolid, MakeWire};
use crate::topology::{FaceId, FaceSurface, ShellData, SolidId, TopologyStore};

/// Tolerance for grouping coplanar faces (plane distance comparison).
const COPLANAR_DISTANCE_TOL: f64 = 1e-6;

/// Tolerance for grouping coplanar faces (normal direction comparison).
const COPLANAR_NORMAL_TOL: f64 = 1e-6;

/// Tolerance for collinearity simplification.
const COLLINEAR_TOL: f64 = 1e-8;

/// Info extracted from one face for the merge algorithm.
struct FaceInfo {
    face_id: FaceId,
    /// Directed edges as (`start_point`, `end_point`) pairs.
    edges: Vec<(Point3, Point3)>,
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

            let face_id = merge_component(store, component, &face_infos)?;
            merged_face_ids.push(face_id);
            for &idx in component {
                consumed.insert(idx);
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
    let FaceSurface::Plane(ref plane) = face.surface;
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
    }

    // Build adjacency graph: face i is adjacent to face j if they share a reversed edge
    let mut adj: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for (local_idx, &global_idx) in group.iter().enumerate() {
        for (start, end) in &face_infos[global_idx].edges {
            // Look for the reverse edge (end→start)
            let reverse_key = EdgeKey::new(end, start);
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

/// Merges a connected component of coplanar faces into a single face.
fn merge_component(
    store: &mut TopologyStore,
    component: &[usize],
    face_infos: &[FaceInfo],
) -> Result<FaceId> {
    // Collect all directed edges and classify as internal or boundary
    let mut edge_count: HashMap<EdgeKey, usize> = HashMap::new();
    let mut edge_points: HashMap<EdgeKey, (Point3, Point3)> = HashMap::new();

    for &idx in component {
        for &(start, end) in &face_infos[idx].edges {
            let key = EdgeKey::new(&start, &end);
            *edge_count.entry(key).or_insert(0) += 1;
            edge_points.entry(key).or_insert((start, end));
        }
    }

    // Boundary edges: those whose reverse does NOT appear
    let mut boundary_edges: Vec<(Point3, Point3)> = Vec::new();
    for (key, &count) in &edge_count {
        if count > 1 {
            continue; // Duplicate forward edge (shouldn't happen, but skip)
        }
        let reverse = EdgeKey::new_from_ints(key.end_key, key.start_key);
        if edge_count.contains_key(&reverse) {
            // Internal edge (both directions exist) — skip
            continue;
        }
        let &(start, end) = edge_points
            .get(key)
            .ok_or_else(|| OperationError::Failed("missing edge points".into()))?;
        boundary_edges.push((start, end));
    }

    if boundary_edges.is_empty() {
        return Err(OperationError::Failed("merge produced no boundary edges".into()).into());
    }

    // Chain boundary edges into closed loops
    let loops = chain_into_loops(&boundary_edges)?;

    if loops.is_empty() {
        return Err(OperationError::Failed("merge produced no loops".into()).into());
    }

    // Get the face normal from the first face in the component
    let normal = face_infos[component[0]].normal;

    // Classify loops: largest area = outer, rest = inner holes
    let (outer_loop, inner_loops) = classify_loops(&loops, &normal)?;

    // Simplify collinear vertices
    let outer_simplified = simplify_collinear(&outer_loop);
    let inners_simplified: Vec<Vec<Point3>> =
        inner_loops.iter().map(|l| simplify_collinear(l)).collect();

    // Create outer wire
    let outer_wire = MakeWire::new(outer_simplified, true).execute(store)?;

    // Create inner wires
    let mut inner_wires = Vec::with_capacity(inners_simplified.len());
    for inner in &inners_simplified {
        let wire = MakeWire::new(inner.clone(), true).execute(store)?;
        inner_wires.push(wire);
    }

    // Create the merged face
    MakeFace::new(outer_wire, inner_wires).execute(store)
}

/// Chains directed boundary edges into closed loops.
fn chain_into_loops(edges: &[(Point3, Point3)]) -> Result<Vec<Vec<Point3>>> {
    // Build adjacency: start_point → list of (end_point, used_index)
    let mut start_map: HashMap<PointKey, Vec<(usize, Point3)>> = HashMap::new();
    for (i, &(start, end)) in edges.iter().enumerate() {
        let key = PointKey::from_point(&start);
        start_map.entry(key).or_default().push((i, end));
    }

    let mut used = vec![false; edges.len()];
    let mut loops: Vec<Vec<Point3>> = Vec::new();

    for seed_idx in 0..edges.len() {
        if used[seed_idx] {
            continue;
        }

        used[seed_idx] = true;
        let mut chain = vec![edges[seed_idx].0];
        let mut current_end = edges[seed_idx].1;
        let start_key = PointKey::from_point(&edges[seed_idx].0);

        loop {
            let end_key = PointKey::from_point(&current_end);

            // Check if we've closed the loop
            if end_key == start_key {
                break;
            }

            chain.push(current_end);

            // Find the next edge starting from current_end
            let next = start_map.get(&end_key).and_then(|candidates| {
                candidates
                    .iter()
                    .find(|&&(idx, _)| !used[idx])
                    .copied()
            });

            match next {
                Some((idx, end_pt)) => {
                    used[idx] = true;
                    current_end = end_pt;
                }
                None => {
                    return Err(OperationError::Failed(
                        "boundary edges do not form a closed loop".into(),
                    )
                    .into());
                }
            }
        }

        if chain.len() >= 3 {
            loops.push(chain);
        }
    }

    Ok(loops)
}

/// Classifies loops into outer boundary (largest area) and inner holes.
fn classify_loops(
    loops: &[Vec<Point3>],
    normal: &Vector3,
) -> Result<(Vec<Point3>, Vec<Vec<Point3>>)> {
    if loops.is_empty() {
        return Err(OperationError::Failed("no loops to classify".into()).into());
    }

    if loops.len() == 1 {
        return Ok((loops[0].clone(), Vec::new()));
    }

    // Find the loop with the largest area — that's the outer boundary
    let mut max_area = f64::NEG_INFINITY;
    let mut max_idx = 0;
    for (i, loop_pts) in loops.iter().enumerate() {
        let area = polygon_area_3d(loop_pts, normal);
        if area > max_area {
            max_area = area;
            max_idx = i;
        }
    }

    let outer = loops[max_idx].clone();
    let inners: Vec<Vec<Point3>> = loops
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != max_idx)
        .map(|(_, l)| l.clone())
        .collect();

    Ok((outer, inners))
}

/// Removes collinear mid-vertices from a polygon loop.
pub(crate) fn simplify_collinear(points: &[Point3]) -> Vec<Point3> {
    let n = points.len();
    if n < 3 {
        return points.to_vec();
    }

    let mut result: Vec<Point3> = Vec::with_capacity(n);
    for i in 0..n {
        let prev = points[(i + n - 1) % n];
        let curr = points[i];
        let next = points[(i + 1) % n];

        if !is_collinear(&prev, &curr, &next) {
            result.push(curr);
        }
    }

    // Edge case: if all points are collinear (degenerate polygon), keep at least the originals
    if result.len() < 3 {
        return points.to_vec();
    }

    result
}

/// Checks if three points are collinear.
fn is_collinear(a: &Point3, b: &Point3, c: &Point3) -> bool {
    let ab = b - a;
    let ac = c - a;
    let cross = ab.cross(&ac);
    cross.norm_squared() < COLLINEAR_TOL * COLLINEAR_TOL
}

/// Key for hashing directed edges by quantized start/end coordinates.
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

    fn new_from_ints(start_key: (i64, i64, i64), end_key: (i64, i64, i64)) -> Self {
        Self { start_key, end_key }
    }
}

/// Key for hashing points by quantized coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct PointKey {
    x: i64,
    y: i64,
    z: i64,
}

impl PointKey {
    fn from_point(p: &Point3) -> Self {
        let (x, y, z) = quantize(p);
        Self { x, y, z }
    }
}

/// Quantizes a point to integer coordinates for hashing.
/// Uses a grid resolution that's fine enough to distinguish separate vertices
/// but coarse enough to merge coincident ones.
#[allow(clippy::cast_possible_truncation)]
fn quantize(p: &Point3) -> (i64, i64, i64) {
    const INV_GRID: f64 = 1e6; // 1 micron resolution
    (
        (p.x * INV_GRID).round() as i64,
        (p.y * INV_GRID).round() as i64,
        (p.z * INV_GRID).round() as i64,
    )
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::TOLERANCE;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn simplify_collinear_removes_mid_vertex() {
        // Square with an extra collinear point on one edge
        let points = vec![
            p(0.0, 0.0, 0.0),
            p(2.0, 0.0, 0.0), // collinear mid-point on bottom edge
            p(4.0, 0.0, 0.0),
            p(4.0, 4.0, 0.0),
            p(0.0, 4.0, 0.0),
        ];
        let simplified = simplify_collinear(&points);
        assert_eq!(simplified.len(), 4, "should remove the collinear mid-point");
        // Verify the mid-point (2,0,0) was removed
        assert!(!simplified.iter().any(|p| (p.x - 2.0).abs() < TOLERANCE && p.y.abs() < TOLERANCE));
    }

    #[test]
    fn simplify_collinear_keeps_non_collinear() {
        let points = vec![
            p(0.0, 0.0, 0.0),
            p(4.0, 0.0, 0.0),
            p(4.0, 4.0, 0.0),
            p(0.0, 4.0, 0.0),
        ];
        let simplified = simplify_collinear(&points);
        assert_eq!(simplified.len(), 4);
    }

    #[test]
    fn chain_two_triangles_into_rectangle() {
        // Two triangles forming a rectangle, boundary edges only
        // Rectangle: (0,0) (4,0) (4,3) (0,3)
        // Shared diagonal: (0,0)→(4,3) is internal
        // Boundary edges of the combined shape:
        let boundary = vec![
            (p(0.0, 0.0, 0.0), p(4.0, 0.0, 0.0)),
            (p(4.0, 0.0, 0.0), p(4.0, 3.0, 0.0)),
            (p(4.0, 3.0, 0.0), p(0.0, 3.0, 0.0)),
            (p(0.0, 3.0, 0.0), p(0.0, 0.0, 0.0)),
        ];
        let loops = chain_into_loops(&boundary).unwrap();
        assert_eq!(loops.len(), 1);
        assert_eq!(loops[0].len(), 4);
    }

    #[test]
    fn chain_with_hole_produces_two_loops() {
        // Outer boundary (CCW square)
        // Inner boundary (CW square hole)
        let boundary = vec![
            // Outer
            (p(0.0, 0.0, 0.0), p(10.0, 0.0, 0.0)),
            (p(10.0, 0.0, 0.0), p(10.0, 10.0, 0.0)),
            (p(10.0, 10.0, 0.0), p(0.0, 10.0, 0.0)),
            (p(0.0, 10.0, 0.0), p(0.0, 0.0, 0.0)),
            // Inner hole
            (p(3.0, 3.0, 0.0), p(3.0, 7.0, 0.0)),
            (p(3.0, 7.0, 0.0), p(7.0, 7.0, 0.0)),
            (p(7.0, 7.0, 0.0), p(7.0, 3.0, 0.0)),
            (p(7.0, 3.0, 0.0), p(3.0, 3.0, 0.0)),
        ];
        let loops = chain_into_loops(&boundary).unwrap();
        assert_eq!(loops.len(), 2);
    }

    #[test]
    fn classify_loops_largest_is_outer() {
        let outer = vec![
            p(0.0, 0.0, 0.0),
            p(10.0, 0.0, 0.0),
            p(10.0, 10.0, 0.0),
            p(0.0, 10.0, 0.0),
        ];
        let inner = vec![
            p(3.0, 3.0, 0.0),
            p(7.0, 3.0, 0.0),
            p(7.0, 7.0, 0.0),
            p(3.0, 7.0, 0.0),
        ];
        let normal = Vector3::new(0.0, 0.0, 1.0);
        let (classified_outer, classified_inners) =
            classify_loops(&[outer.clone(), inner.clone()], &normal).unwrap();

        let outer_area = polygon_area_3d(&classified_outer, &normal);
        let inner_area = polygon_area_3d(&classified_inners[0], &normal);
        assert!(outer_area > inner_area);
    }
}
