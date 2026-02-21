use crate::error::Result;
use crate::math::intersect_3d::{plane_plane_intersect, PlanePairRelation};
use crate::math::polygon_3d::clip_segment_to_polygon;
use crate::math::{Point3, TOLERANCE};
use crate::topology::{FaceId, FaceSurface, TopologyStore};

/// A segment where two faces intersect in 3D space.
#[derive(Debug, Clone)]
pub struct FaceFaceIntersection {
    pub start: Point3,
    pub end: Point3,
    pub face_a: FaceId,
    pub face_b: FaceId,
}

/// Computes the intersection segments between two planar faces.
///
/// Returns an empty list if the faces don't intersect (parallel, coincident,
/// or the intersection line doesn't pass through both faces).
///
/// # Errors
///
/// Returns an error if face or wire topology cannot be read.
pub fn intersect_face_face(
    store: &TopologyStore,
    face_a: FaceId,
    face_b: FaceId,
) -> Result<Vec<FaceFaceIntersection>> {
    let fa = store.face(face_a)?;
    let fb = store.face(face_b)?;

    let FaceSurface::Plane(ref plane_a) = fa.surface else {
        todo!("Face intersection for non-planar faces")
    };
    let FaceSurface::Plane(ref plane_b) = fb.surface else {
        todo!("Face intersection for non-planar faces")
    };

    // Step 1: intersect the two planes
    let relation = plane_plane_intersect(plane_a, plane_b);
    let (line_origin, line_dir) = match relation {
        PlanePairRelation::IntersectionLine { origin, direction } => (origin, direction),
        PlanePairRelation::Parallel { .. } | PlanePairRelation::Coincident => {
            return Ok(Vec::new());
        }
    };

    // Step 2: collect polygon vertices for both faces
    let poly_a = collect_face_polygon(store, face_a)?;
    let poly_b = collect_face_polygon(store, face_b)?;

    if poly_a.len() < 3 || poly_b.len() < 3 {
        return Ok(Vec::new());
    }

    // Step 3: create a large segment along the intersection line, then clip to each polygon
    // Find extent by projecting polygon points onto the line direction
    let all_points: Vec<&Point3> = poly_a.iter().chain(poly_b.iter()).collect();
    let (t_min, t_max) = line_extent(&line_origin, &line_dir, &all_points);
    let margin = 1.0; // small margin to avoid clipping issues at exact boundaries
    let seg_start = line_origin + line_dir * (t_min - margin);
    let seg_end = line_origin + line_dir * (t_max + margin);

    // Step 4: clip segment to face A polygon
    let intervals_a = clip_segment_to_polygon(&seg_start, &seg_end, &poly_a, plane_a);

    // Step 5: clip segment to face B polygon
    let intervals_b = clip_segment_to_polygon(&seg_start, &seg_end, &poly_b, plane_b);

    // Step 6: find overlap of 1D intervals
    let mut results = Vec::new();
    for &(a0, a1) in &intervals_a {
        for &(b0, b1) in &intervals_b {
            let overlap_start = a0.max(b0);
            let overlap_end = a1.min(b1);
            if overlap_end - overlap_start > TOLERANCE {
                let seg_dir = seg_end - seg_start;
                let start = seg_start + seg_dir * overlap_start;
                let end = seg_start + seg_dir * overlap_end;
                results.push(FaceFaceIntersection {
                    start,
                    end,
                    face_a,
                    face_b,
                });
            }
        }
    }

    Ok(results)
}

/// Collects the outer polygon vertices of a face.
pub(crate) fn collect_face_polygon(
    store: &TopologyStore,
    face_id: FaceId,
) -> Result<Vec<Point3>> {
    let face = store.face(face_id)?;
    let wire = store.wire(face.outer_wire)?;
    let mut polygon = Vec::with_capacity(wire.edges.len());
    for oe in &wire.edges {
        let edge = store.edge(oe.edge)?;
        let vid = if oe.forward { edge.start } else { edge.end };
        let vertex = store.vertex(vid)?;
        polygon.push(vertex.point);
    }
    Ok(polygon)
}

/// Computes the min/max projection of a set of points onto a line.
fn line_extent(origin: &Point3, dir: &crate::math::Vector3, points: &[&Point3]) -> (f64, f64) {
    let mut t_min = f64::INFINITY;
    let mut t_max = f64::NEG_INFINITY;
    for p in points {
        let t = (*p - origin).dot(dir);
        if t < t_min {
            t_min = t;
        }
        if t > t_max {
            t_max = t;
        }
    }
    (t_min, t_max)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::creation::{MakeFace, MakeWire};

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn make_face(store: &mut TopologyStore, points: Vec<Point3>) -> FaceId {
        let wire = MakeWire::new(points, true).execute(store).unwrap();
        MakeFace::new(wire, vec![]).execute(store).unwrap()
    }

    #[test]
    fn two_perpendicular_quads_crossing() {
        let mut store = TopologyStore::new();

        // Face A: XY-plane quad at z=0, from (−1,−1) to (1,1)
        let face_a = make_face(
            &mut store,
            vec![
                p(-1.0, -1.0, 0.0),
                p(1.0, -1.0, 0.0),
                p(1.0, 1.0, 0.0),
                p(-1.0, 1.0, 0.0),
            ],
        );

        // Face B: XZ-plane quad at y=0, from (−1,−1) to (1,1)
        let face_b = make_face(
            &mut store,
            vec![
                p(-1.0, 0.0, -1.0),
                p(1.0, 0.0, -1.0),
                p(1.0, 0.0, 1.0),
                p(-1.0, 0.0, 1.0),
            ],
        );

        let intersections = intersect_face_face(&store, face_a, face_b).unwrap();
        assert_eq!(intersections.len(), 1, "expected 1 intersection segment");

        let seg = &intersections[0];
        // Intersection should be along the X-axis at y=0, z=0, from x=-1 to x=1
        let dx = (seg.end - seg.start).norm();
        assert!(
            (dx - 2.0).abs() < 0.01,
            "intersection length should be ~2.0, got {dx}"
        );
    }

    #[test]
    fn parallel_quads_no_intersection() {
        let mut store = TopologyStore::new();

        let face_a = make_face(
            &mut store,
            vec![
                p(0.0, 0.0, 0.0),
                p(1.0, 0.0, 0.0),
                p(1.0, 1.0, 0.0),
                p(0.0, 1.0, 0.0),
            ],
        );

        let face_b = make_face(
            &mut store,
            vec![
                p(0.0, 0.0, 5.0),
                p(1.0, 0.0, 5.0),
                p(1.0, 1.0, 5.0),
                p(0.0, 1.0, 5.0),
            ],
        );

        let intersections = intersect_face_face(&store, face_a, face_b).unwrap();
        assert!(intersections.is_empty());
    }

    #[test]
    fn non_overlapping_quads_no_intersection() {
        let mut store = TopologyStore::new();

        // Two perpendicular quads that don't overlap
        let face_a = make_face(
            &mut store,
            vec![
                p(0.0, 0.0, 0.0),
                p(1.0, 0.0, 0.0),
                p(1.0, 1.0, 0.0),
                p(0.0, 1.0, 0.0),
            ],
        );

        // XZ quad at y=5 (far away from face_a)
        let face_b = make_face(
            &mut store,
            vec![
                p(0.0, 5.0, 0.0),
                p(1.0, 5.0, 0.0),
                p(1.0, 5.0, 1.0),
                p(0.0, 5.0, 1.0),
            ],
        );

        let intersections = intersect_face_face(&store, face_a, face_b).unwrap();
        assert!(intersections.is_empty());
    }
}
