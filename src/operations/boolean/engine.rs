use std::collections::HashMap;

use crate::error::{OperationError, Result};
use crate::math::{Point3, TOLERANCE};
use crate::topology::{FaceId, FaceSurface, SolidId, TopologyStore};

use super::assemble::assemble_result;
use super::classify::{classify_point_in_solid, PointClassification};
use super::face_intersection::intersect_face_face;
use super::select::{should_keep_fragment, BooleanOp, KeepDecision};
use super::split::{split_face, FaceFragment, SolidSource};

/// Executes a boolean operation on two solids.
///
/// Orchestrates the full pipeline: face-face intersection, splitting,
/// classification, selection, and assembly.
pub fn boolean_execute(
    store: &mut TopologyStore,
    solid_a: SolidId,
    solid_b: SolidId,
    op: BooleanOp,
) -> Result<SolidId> {
    // Step 1: AABB early-out
    let aabb_a = compute_solid_aabb(store, solid_a)?;
    let aabb_b = compute_solid_aabb(store, solid_b)?;

    if !aabb_overlap(&aabb_a, &aabb_b) {
        return handle_disjoint(store, solid_a, solid_b, op);
    }

    // Step 2: Collect faces from both solids
    let faces_a = collect_solid_faces(store, solid_a)?;
    let faces_b = collect_solid_faces(store, solid_b)?;

    // Step 3: Compute all face-face intersections
    let mut cuts_by_face: HashMap<FaceId, Vec<(Point3, Point3)>> = HashMap::new();

    for &fa in &faces_a {
        for &fb in &faces_b {
            let intersections = intersect_face_face(store, fa, fb)?;
            for isect in intersections {
                cuts_by_face
                    .entry(fa)
                    .or_default()
                    .push((isect.start, isect.end));
                cuts_by_face
                    .entry(fb)
                    .or_default()
                    .push((isect.start, isect.end));
            }
        }
    }

    // If no intersections found, handle like disjoint or contained
    if cuts_by_face.is_empty() {
        return handle_no_intersection(store, solid_a, solid_b, op);
    }

    // Step 4: Split faces into fragments
    let mut all_fragments: Vec<(FaceFragment, KeepDecision)> = Vec::new();

    // Split faces from solid A
    for &face_id in &faces_a {
        let cuts = cuts_by_face.get(&face_id).map_or(&[][..], |v| v.as_slice());
        let fragments = split_face(store, face_id, cuts, SolidSource::A)?;
        for frag in fragments {
            let classification = classify_fragment_centroid(store, &frag, solid_b)?;
            let decision = should_keep_fragment(frag.source, classification, op);
            all_fragments.push((frag, decision));
        }
    }

    // Split faces from solid B
    for &face_id in &faces_b {
        let cuts = cuts_by_face.get(&face_id).map_or(&[][..], |v| v.as_slice());
        let fragments = split_face(store, face_id, cuts, SolidSource::B)?;
        for frag in fragments {
            let classification = classify_fragment_centroid(store, &frag, solid_a)?;
            let decision = should_keep_fragment(frag.source, classification, op);
            all_fragments.push((frag, decision));
        }
    }

    // Step 5: Check if we have any kept fragments
    let kept_count = all_fragments
        .iter()
        .filter(|(_, d)| *d != KeepDecision::Discard)
        .count();

    if kept_count == 0 {
        return Err(
            OperationError::Failed("boolean operation produced empty result".into()).into(),
        );
    }

    // Step 6: Assemble the result
    let assembled = assemble_result(store, &all_fragments)?;

    // Step 7: Merge coplanar adjacent faces
    super::merge::merge_coplanar_faces(store, assembled)
}

/// Classifies a fragment's centroid against the other solid.
fn classify_fragment_centroid(
    store: &TopologyStore,
    fragment: &FaceFragment,
    other_solid: SolidId,
) -> Result<PointClassification> {
    let centroid = polygon_centroid(&fragment.boundary);
    // Offset centroid slightly inward from the face plane to avoid boundary issues
    let normal = fragment.plane.plane_normal();
    // Use same_sense to determine inward direction
    let inward_dir = if fragment.same_sense {
        -normal
    } else {
        *normal
    };
    let test_point = centroid + inward_dir * (TOLERANCE * 100.0);
    classify_point_in_solid(&test_point, other_solid, store)
}

/// Computes the centroid of a polygon.
fn polygon_centroid(points: &[Point3]) -> Point3 {
    let n = points.len();
    if n == 0 {
        return Point3::new(0.0, 0.0, 0.0);
    }
    #[allow(clippy::cast_precision_loss)]
    let inv_n = 1.0 / n as f64;
    Point3::new(
        points.iter().map(|p| p.x).sum::<f64>() * inv_n,
        points.iter().map(|p| p.y).sum::<f64>() * inv_n,
        points.iter().map(|p| p.z).sum::<f64>() * inv_n,
    )
}

/// Collects all face IDs from a solid's outer shell.
fn collect_solid_faces(store: &TopologyStore, solid_id: SolidId) -> Result<Vec<FaceId>> {
    let solid = store.solid(solid_id)?;
    let shell = store.shell(solid.outer_shell)?;
    Ok(shell.faces.clone())
}

/// Axis-Aligned Bounding Box.
#[derive(Debug)]
struct Aabb {
    min: Point3,
    max: Point3,
}

/// Computes the AABB of a solid.
fn compute_solid_aabb(store: &TopologyStore, solid_id: SolidId) -> Result<Aabb> {
    let solid = store.solid(solid_id)?;
    let shell = store.shell(solid.outer_shell)?;

    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut min_z = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut max_z = f64::NEG_INFINITY;

    for &face_id in &shell.faces {
        let face = store.face(face_id)?;
        let wire = store.wire(face.outer_wire)?;
        for oe in &wire.edges {
            let edge = store.edge(oe.edge)?;
            for &vid in &[edge.start, edge.end] {
                let v = store.vertex(vid)?;
                min_x = min_x.min(v.point.x);
                min_y = min_y.min(v.point.y);
                min_z = min_z.min(v.point.z);
                max_x = max_x.max(v.point.x);
                max_y = max_y.max(v.point.y);
                max_z = max_z.max(v.point.z);
            }
        }
    }

    Ok(Aabb {
        min: Point3::new(min_x, min_y, min_z),
        max: Point3::new(max_x, max_y, max_z),
    })
}

/// Checks if two AABBs overlap.
fn aabb_overlap(a: &Aabb, b: &Aabb) -> bool {
    a.min.x <= b.max.x + TOLERANCE
        && a.max.x >= b.min.x - TOLERANCE
        && a.min.y <= b.max.y + TOLERANCE
        && a.max.y >= b.min.y - TOLERANCE
        && a.min.z <= b.max.z + TOLERANCE
        && a.max.z >= b.min.z - TOLERANCE
}

/// Handles the case where solids are disjoint (AABBs don't overlap).
fn handle_disjoint(
    store: &mut TopologyStore,
    solid_a: SolidId,
    _solid_b: SolidId,
    op: BooleanOp,
) -> Result<SolidId> {
    match op {
        BooleanOp::Union => {
            // For disjoint union, we'd need multi-shell solids.
            // For now, return error — not common in architecture.
            Err(OperationError::Failed(
                "union of disjoint solids is not yet supported".into(),
            )
            .into())
        }
        BooleanOp::Subtract => {
            // A - B where they don't overlap = A (copy the faces)
            copy_solid(store, solid_a)
        }
        BooleanOp::Intersect => Err(OperationError::Failed(
            "intersection of disjoint solids produces empty result".into(),
        )
        .into()),
    }
}

/// Handles the case where solids overlap in AABB but no face-face intersections.
/// This means one solid is entirely inside the other, or they share a face.
fn handle_no_intersection(
    store: &mut TopologyStore,
    solid_a: SolidId,
    solid_b: SolidId,
    op: BooleanOp,
) -> Result<SolidId> {
    // Sample a point from solid B to see if it's inside A
    let b_sample = get_solid_vertex(store, solid_b)?;
    let b_in_a = classify_point_in_solid(&b_sample, solid_a, store)?;

    // Sample a point from solid A to see if it's inside B
    let a_sample = get_solid_vertex(store, solid_a)?;
    let a_in_b = classify_point_in_solid(&a_sample, solid_b, store)?;

    match (a_in_b, b_in_a) {
        (PointClassification::Inside, _) => {
            // A is inside B
            match op {
                BooleanOp::Union => copy_solid(store, solid_b),
                BooleanOp::Subtract => Err(OperationError::Failed(
                    "subtraction where A is inside B produces empty result".into(),
                )
                .into()),
                BooleanOp::Intersect => copy_solid(store, solid_a),
            }
        }
        (_, PointClassification::Inside) => {
            // B is inside A
            match op {
                BooleanOp::Union => copy_solid(store, solid_a),
                BooleanOp::Subtract => {
                    // A - B where B is inside A => A with a void
                    // For now, return error — inner shells need more work
                    Err(OperationError::Failed(
                        "subtraction where B is fully inside A (void creation) is not yet supported".into(),
                    )
                    .into())
                }
                BooleanOp::Intersect => copy_solid(store, solid_b),
            }
        }
        _ => {
            // No clear containment — treat as disjoint
            handle_disjoint(store, solid_a, solid_b, op)
        }
    }
}

/// Gets a vertex position from a solid (for containment testing).
fn get_solid_vertex(store: &TopologyStore, solid_id: SolidId) -> Result<Point3> {
    let solid = store.solid(solid_id)?;
    let shell = store.shell(solid.outer_shell)?;
    let face = store.face(shell.faces[0])?;
    let wire = store.wire(face.outer_wire)?;
    let oe = &wire.edges[0];
    let edge = store.edge(oe.edge)?;
    let vertex = store.vertex(edge.start)?;
    Ok(vertex.point)
}

/// Creates a copy of a solid by duplicating all its topology.
fn copy_solid(store: &mut TopologyStore, solid_id: SolidId) -> Result<SolidId> {
    use super::face_intersection::collect_face_polygon;
    use super::split::SolidSource;

    let faces = collect_solid_faces(store, solid_id)?;
    let mut fragments = Vec::new();

    for face_id in &faces {
        let polygon = collect_face_polygon(store, *face_id)?;
        let face = store.face(*face_id)?;
        let FaceSurface::Plane(ref plane) = face.surface;

        fragments.push((
            FaceFragment {
                boundary: polygon,
                plane: plane.clone(),
                same_sense: face.same_sense,
                source_face: *face_id,
                source: SolidSource::A,
            },
            KeepDecision::Keep,
        ));
    }

    assemble_result(store, &fragments)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Vector3;
    use crate::operations::creation::{MakeFace, MakeWire};
    use crate::operations::shaping::Extrude;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn make_box(
        store: &mut TopologyStore,
        x: f64,
        y: f64,
        z: f64,
        dx: f64,
        dy: f64,
        dz: f64,
    ) -> SolidId {
        let pts = vec![
            p(x, y, z),
            p(x + dx, y, z),
            p(x + dx, y + dy, z),
            p(x, y + dy, z),
        ];
        let wire = MakeWire::new(pts, true).execute(store).unwrap();
        let face = MakeFace::new(wire, vec![]).execute(store).unwrap();
        Extrude::new(face, Vector3::new(0.0, 0.0, dz))
            .execute(store)
            .unwrap()
    }

    #[test]
    fn subtract_overlapping_boxes() {
        let mut store = TopologyStore::new();
        // Large box: 0..4 x 0..4 x 0..4
        let a = make_box(&mut store, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0);
        // Small box: 1..3 x 1..3 x -0.5..4.5 (extends beyond wall to avoid coplanar faces)
        let b = make_box(&mut store, 1.0, 1.0, -0.5, 2.0, 2.0, 5.0);

        let result = boolean_execute(&mut store, a, b, BooleanOp::Subtract);
        assert!(result.is_ok(), "subtract failed: {result:?}");

        let solid_id = result.unwrap();
        let solid = store.solid(solid_id).unwrap();
        let shell = store.shell(solid.outer_shell).unwrap();
        // After merge: 4 outer walls + 1 top (with hole) + 1 bottom (with hole) + 4 inner walls = 10
        assert_eq!(shell.faces.len(), 10, "expected 10 faces after merge");

        // Verify that exactly 2 faces have inner wires (the top and bottom with holes)
        let faces_with_holes = shell
            .faces
            .iter()
            .filter(|&&fid| {
                let face = store.face(fid).unwrap();
                !face.inner_wires.is_empty()
            })
            .count();
        assert_eq!(faces_with_holes, 2, "expected 2 faces with inner wires (holes)");
    }

    #[test]
    fn union_overlapping_boxes() {
        let mut store = TopologyStore::new();
        let a = make_box(&mut store, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let b = make_box(&mut store, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0);

        let result = boolean_execute(&mut store, a, b, BooleanOp::Union);
        assert!(result.is_ok(), "union failed: {result:?}");

        let solid_id = result.unwrap();
        let solid = store.solid(solid_id).unwrap();
        let shell = store.shell(solid.outer_shell).unwrap();
        // Union of two overlapping boxes produces more than 6 faces
        assert!(
            shell.faces.len() > 6,
            "expected >6 faces, got {}",
            shell.faces.len()
        );
    }

    #[test]
    fn intersect_overlapping_boxes() {
        let mut store = TopologyStore::new();
        let a = make_box(&mut store, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let b = make_box(&mut store, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0);

        let result = boolean_execute(&mut store, a, b, BooleanOp::Intersect);
        assert!(result.is_ok(), "intersect failed: {result:?}");

        let solid_id = result.unwrap();
        let solid = store.solid(solid_id).unwrap();
        let shell = store.shell(solid.outer_shell).unwrap();
        // Intersection should produce a box-like shape
        assert!(
            shell.faces.len() >= 6,
            "expected >=6 faces, got {}",
            shell.faces.len()
        );
    }

    #[test]
    fn subtract_disjoint_returns_copy() {
        let mut store = TopologyStore::new();
        let a = make_box(&mut store, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = make_box(&mut store, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0);

        let result = boolean_execute(&mut store, a, b, BooleanOp::Subtract);
        assert!(result.is_ok());

        let solid_id = result.unwrap();
        let solid = store.solid(solid_id).unwrap();
        let shell = store.shell(solid.outer_shell).unwrap();
        assert_eq!(shell.faces.len(), 6);
    }

    #[test]
    fn intersect_disjoint_returns_error() {
        let mut store = TopologyStore::new();
        let a = make_box(&mut store, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = make_box(&mut store, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0);

        let result = boolean_execute(&mut store, a, b, BooleanOp::Intersect);
        assert!(result.is_err());
    }
}
