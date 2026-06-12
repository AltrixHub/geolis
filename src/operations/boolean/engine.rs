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
    // NURBS routing: if either solid has a NURBS face, the planar pipeline does
    // not apply. The through-cut subtract handles it; everything else returns an
    // explicit unsupported error.
    if solid_has_nurbs_face(store, solid_a)? || solid_has_nurbs_face(store, solid_b)? {
        return super::nurbs::try_boolean(store, solid_a, solid_b, op);
    }

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

/// Whether any face of the solid's outer shell is a NURBS face.
fn solid_has_nurbs_face(store: &TopologyStore, solid_id: SolidId) -> Result<bool> {
    let shell = store.shell(store.solid(solid_id)?.outer_shell)?;
    for &fid in &shell.faces {
        if matches!(store.face(fid)?.surface, FaceSurface::Nurbs(_)) {
            return Ok(true);
        }
    }
    Ok(false)
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
            Err(
                OperationError::Failed("union of disjoint solids is not yet supported".into())
                    .into(),
            )
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
    use super::face_intersection::{collect_face_polygon, collect_inner_wire_polygons};
    use super::split::SolidSource;

    let faces = collect_solid_faces(store, solid_id)?;
    let mut fragments = Vec::new();

    for face_id in &faces {
        let polygon = collect_face_polygon(store, *face_id)?;
        let inner_polygons = collect_inner_wire_polygons(store, *face_id)?;
        let face = store.face(*face_id)?;
        let FaceSurface::Plane(ref plane) = face.surface else {
            if matches!(face.surface, FaceSurface::Nurbs(_)) {
                return Err(OperationError::Failed(
                    "boolean operations on NURBS faces are not yet supported".into(),
                )
                .into());
            }
            todo!("Boolean operations for non-planar faces")
        };

        fragments.push((
            FaceFragment {
                boundary: polygon,
                inner_boundaries: inner_polygons,
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
        assert_eq!(
            faces_with_holes, 2,
            "expected 2 faces with inner wires (holes)"
        );
    }

    /// Union/Intersect on a NURBS-faced solid stay unsupported (only the
    /// through-cut Subtract is handled). Replaces the previous pin that asserted
    /// every NURBS boolean errors.
    #[test]
    fn nurbs_faced_solid_union_intersect_unsupported() {
        use crate::geometry::nurbs::{KnotVector, NurbsSurface};

        let build = || {
            let mut store = TopologyStore::new();
            let a = make_box(&mut store, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0);
            let b = make_box(&mut store, 1.0, 1.0, -0.5, 2.0, 2.0, 5.0);
            let face_a = collect_solid_faces(&store, a).unwrap()[0];
            let patch = NurbsSurface::from_unweighted(
                vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(0.0, 4.0, 0.0),
                    Point3::new(4.0, 0.0, 0.0),
                    Point3::new(4.0, 4.0, 0.0),
                ],
                2,
                2,
                KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
                KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
                1,
                1,
            )
            .unwrap();
            store.face_mut(face_a).unwrap().surface = FaceSurface::Nurbs(patch);
            (store, a, b)
        };

        let (mut store, a, b) = build();
        assert!(
            boolean_execute(&mut store, a, b, BooleanOp::Union).is_err(),
            "Union on a NURBS-faced solid must stay unsupported"
        );
        let (mut store, a, b) = build();
        assert!(
            boolean_execute(&mut store, a, b, BooleanOp::Intersect).is_err(),
            "Intersect on a NURBS-faced solid must stay unsupported"
        );
        // A NURBS Subtract that is NOT a clean through-cut also errors (no panic).
        let (mut store, a, b) = build();
        assert!(
            boolean_execute(&mut store, a, b, BooleanOp::Subtract).is_err(),
            "non-through-cut NURBS Subtract must error, not panic"
        );
    }

    /// The supported case: a curved NURBS slab minus a NURBS tube punches a real
    /// through-hole and returns a manifold solid via `boolean_execute`.
    #[test]
    fn nurbs_through_cut_subtract_succeeds() {
        use crate::operations::creation::{MakeCurvedSlab, MakeNurbsTube};

        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -1.5), 0.7, 5.0)
            .execute(&mut store)
            .unwrap();
        let result = boolean_execute(&mut store, slab, tube, BooleanOp::Subtract);
        assert!(result.is_ok(), "through-cut subtract failed: {result:?}");
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

    /// Builds an axis-aligned-footprint box rotated about the Z axis by `angle`
    /// (radians) around its planar center. Mirrors how revion authors a wall
    /// solid: an oriented rectangular footprint extruded along +Z.
    #[allow(clippy::too_many_arguments)]
    fn make_rotated_box(
        store: &mut TopologyStore,
        cx: f64,
        cy: f64,
        z: f64,
        length: f64,
        thickness: f64,
        height: f64,
        angle: f64,
    ) -> SolidId {
        let (sin_a, cos_a) = angle.sin_cos();
        // Tangent (along length) and perpendicular (along thickness).
        let tx = cos_a;
        let ty = sin_a;
        let nx = -sin_a;
        let ny = cos_a;
        let half_l = length * 0.5;
        let half_t = thickness * 0.5;

        let corner = |sl: f64, st: f64| p(cx + tx * sl + nx * st, cy + ty * sl + ny * st, z);

        let pts = vec![
            corner(-half_l, -half_t),
            corner(half_l, -half_t),
            corner(half_l, half_t),
            corner(-half_l, half_t),
        ];
        let wire = MakeWire::new(pts, true).execute(store).unwrap();
        let face = MakeFace::new(wire, vec![]).execute(store).unwrap();
        Extrude::new(face, Vector3::new(0.0, 0.0, height))
            .execute(store)
            .unwrap()
    }

    /// Regression test for the "all points collinear, cannot define a plane"
    /// failure observed when subtracting a clean opening box from a rotated
    /// wall box. The split step previously emitted collinear sliver fragments
    /// (4 coincident-edge vertices) whose Newell normal vanished, breaking
    /// downstream face construction.
    ///
    /// Geometry mirrors the failing revion case: wall centered at (-3.0, -1.4)
    /// with tangent (-0.94, 0.33) (~atan2(0.33, -0.94)), thickness 0.18,
    /// height 2.7, length 5; opening centered at (-0.8, -1.2), depth 0.9
    /// (1.2x wall thickness), width 0.9, sill 0, height 2.1.
    #[test]
    fn subtract_clean_opening_from_rotated_wall_succeeds() {
        let mut store = TopologyStore::new();
        let angle = (0.33_f64).atan2(-0.94);
        let wall = make_rotated_box(&mut store, -3.0, -1.4, 0.0, 5.0, 0.18, 2.7, angle);
        // Opening box aligned along the same wall tangent: sill 0, height 2.1,
        // width 0.9 (length-direction). Depth is 1.2 x wall thickness — small
        // overshoot — so the LONG side faces of the opening sit nearly
        // parallel to (and slightly outside of) the wall faces. This is the
        // near-parallel configuration that previously generated collinear
        // sliver fragments inside `split_face` and made the downstream
        // `MakeFace` fail with "all points collinear".
        let opening = make_rotated_box(&mut store, -0.8, -1.2, 0.0, 0.9, 0.18 * 1.2, 2.1, angle);

        let result = boolean_execute(&mut store, wall, opening, BooleanOp::Subtract);
        assert!(
            result.is_ok(),
            "subtract should succeed, got error: {result:?}",
        );
    }

    /// Axis-aligned version of the wall-opening subtraction. The opening's
    /// LONG sides sit slightly outside the wall thickness (1.2x overshoot)
    /// and the opening shares the wall's bottom face (z=0). Both are the
    /// near-parallel face configurations that tripped the "all points
    /// collinear" failure before the Newell sliver-rejection fix.
    #[test]
    fn subtract_clean_opening_through_axis_aligned_wall() {
        let mut store = TopologyStore::new();
        // Wall: 5.0 long (x), 0.18 thick (y), 2.7 tall (z), bottom at z=0.
        let wall = make_box(&mut store, -2.5, -0.09, 0.0, 5.0, 0.18, 2.7);
        // Opening: 0.9 wide (x), 0.216 deep (y, 1.2x wall thickness so it
        // pokes both faces by only 18 mm), 2.1 tall (z), bottom flush at z=0.
        let depth = 0.18 * 1.2;
        let opening = make_box(&mut store, -0.45, -depth * 0.5, 0.0, 0.9, depth, 2.1);

        let result = boolean_execute(&mut store, wall, opening, BooleanOp::Subtract);
        assert!(
            result.is_ok(),
            "subtract with coplanar bottom should succeed, got error: {result:?}",
        );
    }

    /// Regression test pinning the EXACT runtime geometry that revion's window-
    /// placement flow hits: a wall solid and an opening solid captured straight
    /// from the live application logs.
    ///
    /// The wall is a thin (~0.18 m), gently tilted cuboid; the opening is a
    /// smaller cuboid roughly parallel to the wall axis and positioned across
    /// the wall's centerline. The depth of the opening pokes slightly out of
    /// each wall face by ~50 mm (the production overshoot constant). At runtime
    /// `Subtract::execute` consistently fails with
    /// "all points are collinear, cannot define a plane" — this test pins the
    /// coordinates so we can reproduce + fix the bug.
    ///
    /// If this test passes the fix is correct. If it fails (asserts the bug)
    /// we've successfully reproduced the runtime failure locally.
    #[test]
    fn subtract_runtime_wall_opening_does_not_collinear() {
        let mut store = TopologyStore::new();

        // Wall bottom (z=0): captured XY from runtime logs, CCW order.
        let wall_bottom = vec![
            p(-2.651, -0.963, 0.0),
            p(6.165, -2.830, 0.0),
            p(6.128, -3.006, 0.0),
            p(-2.689, -1.139, 0.0),
        ];
        let wall_wire = MakeWire::new(wall_bottom, true)
            .execute(&mut store)
            .unwrap();
        let wall_face = MakeFace::new(wall_wire, vec![])
            .execute(&mut store)
            .unwrap();
        // Wall height 2.4 m (from runtime logs).
        let wall = Extrude::new(wall_face, Vector3::new(0.0, 0.0, 2.4))
            .execute(&mut store)
            .unwrap();

        // Opening bottom (z=0): captured XY from runtime logs, CCW order.
        let opening_bottom = vec![
            p(-0.238, -1.423, 0.0),
            p(-0.296, -1.696, 0.0),
            p(-1.176, -1.510, 0.0),
            p(-1.118, -1.236, 0.0),
        ];
        let opening_wire = MakeWire::new(opening_bottom, true)
            .execute(&mut store)
            .unwrap();
        let opening_face = MakeFace::new(opening_wire, vec![])
            .execute(&mut store)
            .unwrap();
        // Opening height 2.1 m (from runtime logs).
        let opening = Extrude::new(opening_face, Vector3::new(0.0, 0.0, 2.1))
            .execute(&mut store)
            .unwrap();

        let result = boolean_execute(&mut store, wall, opening, BooleanOp::Subtract);
        assert!(
            result.is_ok(),
            "runtime wall/opening subtract should succeed, got: {result:?}",
        );
    }

    /// Same as [`subtract_runtime_wall_opening_does_not_collinear`] but
    /// reconstructs the opening corners via the exact tangent / perpendicular
    /// trigonometry `OpeningSolid3DNode::build_opening_solids` runs at
    /// runtime (rather than copying the truncated 3-decimal corners from the
    /// log). The wall centerline is the midline of the two short captured
    /// wall faces; the opening center sits at the same fractional position
    /// along that axis we observed in the failing session.
    ///
    /// Pinning the box construction this way removes the "log truncation"
    /// degree of freedom — if the failing runtime case lived in this
    /// trigonometric path it would surface here.
    #[test]
    fn subtract_runtime_wall_opening_full_precision() {
        let mut store = TopologyStore::new();

        let wall_bottom = vec![
            p(-2.651, -0.963, 0.0),
            p(6.165, -2.830, 0.0),
            p(6.128, -3.006, 0.0),
            p(-2.689, -1.139, 0.0),
        ];
        let wall_wire = MakeWire::new(wall_bottom, true)
            .execute(&mut store)
            .unwrap();
        let wall_face = MakeFace::new(wall_wire, vec![])
            .execute(&mut store)
            .unwrap();
        let wall = Extrude::new(wall_face, Vector3::new(0.0, 0.0, 2.4))
            .execute(&mut store)
            .unwrap();

        let axis_start_x: f64 = (-2.651_f64 + -2.689_f64) * 0.5;
        let axis_start_y: f64 = (-0.963_f64 + -1.139_f64) * 0.5;
        let axis_end_x: f64 = (6.165_f64 + 6.128_f64) * 0.5;
        let axis_end_y: f64 = (-2.830_f64 + -3.006_f64) * 0.5;
        let dx = axis_end_x - axis_start_x;
        let dy = axis_end_y - axis_start_y;
        let len = (dx * dx + dy * dy).sqrt();
        let tx = dx / len;
        let ty = dy / len;
        let px = -ty;
        let py = tx;
        let wall_thickness = 0.181_f64;
        let overshoot = 0.05_f64;
        let half_d = wall_thickness * 0.5 + overshoot;
        let half_w = 0.9_f64 * 0.5;
        let t = 0.222_f64;
        let cx = axis_start_x + dx * t;
        let cy = axis_start_y + dy * t;
        let sill_z = 0.0_f64;
        let opening_bottom = vec![
            p(
                cx - tx * half_w - px * half_d,
                cy - ty * half_w - py * half_d,
                sill_z,
            ),
            p(
                cx + tx * half_w - px * half_d,
                cy + ty * half_w - py * half_d,
                sill_z,
            ),
            p(
                cx + tx * half_w + px * half_d,
                cy + ty * half_w + py * half_d,
                sill_z,
            ),
            p(
                cx - tx * half_w + px * half_d,
                cy - ty * half_w + py * half_d,
                sill_z,
            ),
        ];
        let opening_wire = MakeWire::new(opening_bottom, true)
            .execute(&mut store)
            .unwrap();
        let opening_face = MakeFace::new(opening_wire, vec![])
            .execute(&mut store)
            .unwrap();
        let opening = Extrude::new(opening_face, Vector3::new(0.0, 0.0, 2.1))
            .execute(&mut store)
            .unwrap();

        let result = boolean_execute(&mut store, wall, opening, BooleanOp::Subtract);
        assert!(
            result.is_ok(),
            "full-precision runtime wall/opening subtract should succeed, got: {result:?}",
        );
    }

    /// Cascade regression: Subtract(wall, door) succeeds, then
    /// Subtract(that_result, window1) fails at runtime with
    /// "invalid input: consecutive points 7 and 8 are coincident".
    ///
    /// Coordinates captured verbatim (`{:.17e}`) from a live modeling
    /// session that placed a door followed by two windows on the same
    /// wall. Runtime imports **all three cutters** into the base store
    /// before the first subtract runs (revion's `BRepSubtract` batches
    /// the imports). This test mirrors that import order — wall + door
    /// + window1 + window2 all materialised first, then the subtract
    /// chain — so a state-pollution failure (e.g. the boolean engine
    /// snapping new fragments onto vertices of an *unrelated* solid
    /// that happens to sit in the same store) shows up here too.
    ///
    /// The failure originates in `merge_coplanar_faces::merge_component`
    /// where the merged outer loop carries two consecutive points whose
    /// distance falls below `MakeWire`'s tolerance.
    #[test]
    fn subtract_runtime_wall_door_then_window_cascade() {
        let mut store = TopologyStore::new();

        // Wall — gently tilted cuboid, 8.78 m long × 0.18 m thick × 2.4 m tall.
        let wall_bottom = vec![
            p(-3.269_483_891_652_068_32, -6.843_402_330_123_231_62e-1, 0.0),
            p(5.258_328_608_347_930_81, -3.538_324_608_012_322_51, 0.0),
            p(5.201_202_641_652_067_80, -3.709_019_141_987_676_79, 0.0),
            p(-3.326_609_858_347_931_77, -8.550_347_669_876_767_75e-1, 0.0),
        ];
        let wall = build_extruded(&mut store, wall_bottom, 2.4);

        // Door — perpendicular cuboid, ~0.9 m wide × ~0.22 m deep × ~2.1 m tall.
        let door_bottom = vec![
            p(-7.497_596_405_132_468_39e-1, -1.508_629_561_657_746_31, 0.0),
            p(-8.183_108_003_719_504_75e-1, -1.713_463_006_253_083_12, 0.0),
            p(-1.671_783_463_576_549_17, -1.427_833_181_075_091_27, 0.0),
            p(-1.603_232_303_717_845_75, -1.222_999_736_479_754_46, 0.0),
        ];
        let door = build_extruded(&mut store, door_bottom, 2.099_999_904_632_568_36);

        // Window 1 — sits at z ∈ [0.9, 2.1] (sill 0.9 m, height 1.2 m).
        // Built BEFORE the first subtract so the store carries it
        // alongside the wall/door — mirroring how revion's BRepSubtract
        // imports every cutter up front.
        let window1_bottom = vec![
            p(
                4.186_251_396_293_492_63,
                -3.160_553_589_522_334_23,
                8.999_999_761_581_420_90e-1,
            ),
            p(
                4.117_700_236_434_789_22,
                -3.365_387_034_117_671_48,
                8.999_999_761_581_420_90e-1,
            ),
            p(
                3.264_227_573_230_190_42,
                -3.079_757_208_939_679_64,
                8.999_999_761_581_420_90e-1,
            ),
            p(
                3.332_778_733_088_893_83,
                -2.874_923_764_344_342_38,
                8.999_999_761_581_420_90e-1,
            ),
        ];
        let window_height = 2.100_000_023_841_857_91 - 8.999_999_761_581_420_90e-1;
        let window1 = build_extruded(&mut store, window1_bottom, window_height);

        // Window 2 — second window placement captured verbatim from
        // the failing live session. iter=2 of the cascade subtracts
        // this from the post-(door + window1) wall and triggers
        // "zero-length vector" inside the boolean engine.
        let window2_bottom = vec![
            p(
                1.532_615_675_427_539_74,
                -2.591_498_101_432_384_35,
                8.999_999_761_581_420_90e-1,
            ),
            p(
                1.476_937_380_768_878_82,
                -2.800_198_677_737_773_87,
                8.999_999_761_581_420_90e-1,
            ),
            p(
                6.073_516_691_992_446_59e-1,
                -2.568_205_789_472_399_75,
                8.999_999_761_581_420_90e-1,
            ),
            p(
                6.630_299_638_579_055_80e-1,
                -2.359_505_213_167_010_23,
                8.999_999_761_581_420_90e-1,
            ),
        ];
        let window2 = build_extruded(&mut store, window2_bottom, window_height);

        // Now run the subtract chain. All three subtracts must succeed.
        let wall_minus_door = boolean_execute(&mut store, wall, door, BooleanOp::Subtract)
            .expect("Subtract(wall, door) must succeed (iter=0)");
        let wall_minus_door_and_w1 =
            boolean_execute(&mut store, wall_minus_door, window1, BooleanOp::Subtract)
                .expect("Subtract(wall_minus_door, window1) must succeed (iter=1)");
        let result = boolean_execute(
            &mut store,
            wall_minus_door_and_w1,
            window2,
            BooleanOp::Subtract,
        );
        assert!(
            result.is_ok(),
            "Subtract(wall_minus_door_and_w1, window2) must succeed (iter=2), got: {result:?}",
        );
    }

    fn build_extruded(store: &mut TopologyStore, bottom: Vec<Point3>, height: f64) -> SolidId {
        let wire = MakeWire::new(bottom, true).execute(store).unwrap();
        let face = MakeFace::new(wire, vec![]).execute(store).unwrap();
        Extrude::new(face, Vector3::new(0.0, 0.0, height))
            .execute(store)
            .unwrap()
    }

    #[test]
    fn subtract_preserves_existing_holes() {
        let mut store = TopologyStore::new();

        // plate: 10x10x4
        let plate = make_box(&mut store, 0.0, 0.0, 0.0, 10.0, 10.0, 4.0);
        // column: 2x2, extends beyond plate in z
        let column = make_box(&mut store, 4.0, 4.0, -0.5, 2.0, 2.0, 5.0);
        // plate - column → top and bottom each get a hole
        let holey = boolean_execute(&mut store, plate, column, BooleanOp::Subtract).unwrap();

        let count_faces_with_holes = |store: &TopologyStore, solid_id| {
            let solid = store.solid(solid_id).unwrap();
            let shell = store.shell(solid.outer_shell).unwrap();
            shell
                .faces
                .iter()
                .filter(|&&fid| {
                    let face = store.face(fid).unwrap();
                    !face.inner_wires.is_empty()
                })
                .count()
        };

        let holey_holes = count_faces_with_holes(&store, holey);
        assert_eq!(
            holey_holes, 2,
            "holey should have 2 faces with holes (top + bottom), got {holey_holes}"
        );

        // corner: cuts off x=7..12 strip, extends beyond plate in z
        let corner = make_box(&mut store, 7.0, 0.0, -0.5, 5.0, 10.0, 5.0);
        // holey - corner → right strip removed; left fragments (x=0..7) retain hole (centroid x=5<7)
        let result = boolean_execute(&mut store, holey, corner, BooleanOp::Subtract).unwrap();

        let result_holes = count_faces_with_holes(&store, result);
        assert_eq!(
            result_holes, 2,
            "result should still have 2 faces with holes, got {result_holes}"
        );
    }
}
