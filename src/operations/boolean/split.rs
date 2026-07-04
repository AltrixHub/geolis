use crate::error::{OperationError, Result};
use crate::geometry::surface::Plane;
use crate::math::polygon_3d::{point_in_polygon_3d, polygon_area_3d};
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::topology::{FaceId, FaceSurface, TopologyStore};

use super::face_intersection::{collect_face_polygon, collect_inner_wire_polygons};

/// Which solid a face fragment originated from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolidSource {
    A,
    B,
}

/// A polygon fragment produced by splitting a face along intersection lines.
#[derive(Debug, Clone)]
pub struct FaceFragment {
    pub boundary: Vec<Point3>,
    pub inner_boundaries: Vec<Vec<Point3>>,
    pub plane: Plane,
    pub same_sense: bool,
    pub source_face: FaceId,
    pub source: SolidSource,
}

/// Splits a face along a set of cut segments.
///
/// Each cut is a `(start, end)` segment lying on the face's plane.
/// If no cuts intersect the face, returns a single fragment (the whole face).
///
/// Uses a 2D projection approach: project the polygon and cuts into the
/// face's UV space, split the polygon, then lift back to 3D.
pub fn split_face(
    store: &TopologyStore,
    face_id: FaceId,
    cuts: &[(Point3, Point3)],
    source: SolidSource,
) -> Result<Vec<FaceFragment>> {
    let face = store.face(face_id)?;
    let FaceSurface::Plane(ref plane) = face.surface else {
        if matches!(face.surface, FaceSurface::Nurbs(_)) {
            return Err(OperationError::Failed(
                "boolean operations on NURBS faces are not yet supported".into(),
            )
            .into());
        }
        todo!("Face splitting for non-planar faces")
    };

    let polygon = collect_face_polygon(store, face_id)?;
    let inner_polygons = collect_inner_wire_polygons(store, face_id)?;

    if polygon.len() < 3 {
        return Ok(Vec::new());
    }

    // Filter cuts that actually cross the polygon interior
    let relevant_cuts: Vec<&(Point3, Point3)> = cuts
        .iter()
        .filter(|(s, e)| {
            let len = (e - s).norm();
            len > TOLERANCE
        })
        .collect();

    if relevant_cuts.is_empty() {
        return Ok(vec![FaceFragment {
            boundary: polygon,
            inner_boundaries: inner_polygons,
            plane: plane.clone(),
            same_sense: face.same_sense,
            source_face: face_id,
            source,
        }]);
    }

    // For Phase 1 with planar faces, we support single straight cuts through a polygon.
    // Split the polygon by each cut line sequentially.
    let mut fragments = vec![polygon];

    for cut in &relevant_cuts {
        let mut next_fragments = Vec::new();
        for poly in &fragments {
            let split_result = split_polygon_by_line(poly, &cut.0, &cut.1, plane);
            next_fragments.extend(split_result);
        }
        fragments = next_fragments;
    }

    // Split inner polygons by the same cut lines
    let mut inner_fragments: Vec<Vec<Point3>> = Vec::new();
    for inner in &inner_polygons {
        let mut current = vec![inner.clone()];
        for cut in &relevant_cuts {
            let mut next = Vec::new();
            for poly in &current {
                let split_result = split_polygon_by_line(poly, &cut.0, &cut.1, plane);
                next.extend(split_result);
            }
            current = next;
        }
        // Filter degenerate inner fragments
        let normal = plane.plane_normal();
        let min_area = TOLERANCE * TOLERANCE;
        for frag in current {
            if frag.len() >= 3
                && polygon_area_3d(&frag, normal) > min_area
                && newell_normal_3d(&frag).norm() > TOLERANCE
            {
                inner_fragments.push(frag);
            }
        }
    }

    // Filter out degenerate outer fragments and associate inner fragments.
    // The Newell normal check rejects 3D-collinear slivers that may have a tiny
    // nonzero projected area but cannot define a plane downstream
    // (mirrors the guard in `MakeFace::compute_plane_from_points`).
    let normal = plane.plane_normal();
    let min_area = TOLERANCE * TOLERANCE;
    let result = fragments
        .into_iter()
        .filter(|f| {
            f.len() >= 3
                && polygon_area_3d(f, normal) > min_area
                && newell_normal_3d(f).norm() > TOLERANCE
        })
        .map(|boundary| {
            let inners = associate_inner_fragments(&boundary, &inner_fragments, plane);
            FaceFragment {
                boundary,
                inner_boundaries: inners,
                plane: plane.clone(),
                same_sense: face.same_sense,
                source_face: face_id,
                source,
            }
        })
        .collect();

    Ok(result)
}

/// Associates inner polygon fragments with an outer boundary using centroid containment.
fn associate_inner_fragments(
    outer_boundary: &[Point3],
    inner_fragments: &[Vec<Point3>],
    plane: &Plane,
) -> Vec<Vec<Point3>> {
    inner_fragments
        .iter()
        .filter(|inner| {
            let centroid = polygon_centroid(inner);
            point_in_polygon_3d(&centroid, outer_boundary, plane)
        })
        .cloned()
        .collect()
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

/// Computes the unnormalized Newell normal of a 3D polygon.
///
/// Returns a zero-magnitude vector when the points are collinear in 3D, which
/// is the same condition that `MakeFace::compute_plane_from_points` uses to
/// reject "all points collinear" inputs. Used to filter out sliver fragments
/// emitted by `split_polygon_by_line` that would later fail face construction.
pub(super) fn newell_normal_3d(points: &[Point3]) -> Vector3 {
    let n = points.len();
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    for i in 0..n {
        let curr = &points[i];
        let next = &points[(i + 1) % n];
        nx += (curr.y - next.y) * (curr.z + next.z);
        ny += (curr.z - next.z) * (curr.x + next.x);
        nz += (curr.x - next.x) * (curr.y + next.y);
    }
    Vector3::new(nx, ny, nz)
}

/// Splits a polygon by an infinite line defined by two points on the face plane.
///
/// Projects to UV space, splits, and returns the resulting polygon(s).
#[allow(clippy::many_single_char_names)]
fn split_polygon_by_line(
    polygon: &[Point3],
    line_p0: &Point3,
    line_p1: &Point3,
    plane: &Plane,
) -> Vec<Vec<Point3>> {
    let n = polygon.len();
    if n < 3 {
        return vec![polygon.to_vec()];
    }

    // Project everything to UV space
    let u_dir = plane.u_dir();
    let v_dir = plane.v_dir();
    let origin = plane.origin();

    let project = |p: &Point3| -> (f64, f64) {
        let diff = p - origin;
        (diff.dot(u_dir), diff.dot(v_dir))
    };

    let unproject = |u: f64, v: f64| -> Point3 { origin + u_dir * u + v_dir * v };

    let poly_uv: Vec<(f64, f64)> = polygon.iter().map(&project).collect();
    let lp0 = project(line_p0);
    let lp1 = project(line_p1);

    // Line direction in UV
    let lu = lp1.0 - lp0.0;
    let lv = lp1.1 - lp0.1;

    // Classify each polygon vertex relative to the line
    // sign = cross product of (line_dir) x (vertex - line_p0)
    let signs: Vec<f64> = poly_uv
        .iter()
        .map(|&(pu, pv)| {
            let du = pu - lp0.0;
            let dv = pv - lp0.1;
            lu * dv - lv * du // cross product
        })
        .collect();

    // Check if all vertices are on the same side
    let has_positive = signs.iter().any(|&s| s > TOLERANCE);
    let has_negative = signs.iter().any(|&s| s < -TOLERANCE);

    if !has_positive || !has_negative {
        // No split — all on one side (or on the line)
        return vec![polygon.to_vec()];
    }

    // Split the polygon into two halves
    let mut side_a: Vec<Point3> = Vec::new();
    let mut side_b: Vec<Point3> = Vec::new();

    for i in 0..n {
        let j = (i + 1) % n;
        let si = signs[i];
        let sj = signs[j];
        let pi = &poly_uv[i];
        let pj = &poly_uv[j];

        // Add current vertex to appropriate side(s)
        if si >= -TOLERANCE {
            side_a.push(polygon[i]);
        }
        if si <= TOLERANCE {
            side_b.push(polygon[i]);
        }

        // Check if edge crosses the line
        if (si > TOLERANCE && sj < -TOLERANCE) || (si < -TOLERANCE && sj > TOLERANCE) {
            // Compute intersection point
            let t = si / (si - sj);
            let u = pi.0 + t * (pj.0 - pi.0);
            let v = pi.1 + t * (pj.1 - pi.1);
            let crossing = unproject(u, v);
            side_a.push(crossing);
            side_b.push(crossing);
        }
    }

    let mut result = Vec::new();
    if side_a.len() >= 3 {
        result.push(side_a);
    }
    if side_b.len() >= 3 {
        result.push(side_b);
    }
    if result.is_empty() {
        result.push(polygon.to_vec());
    }
    result
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::creation::{MakeFace, MakeWire};

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn newell_normal_3d_zero_for_collinear_points() {
        // Four points collinear along the X axis. The Newell normal must be
        // zero, even though the surrounding code may compute a tiny non-zero
        // projected area when the supplied normal is mis-oriented.
        let pts = vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(2.0, 0.0, 0.0),
            p(3.0, 0.0, 0.0),
        ];
        let normal = newell_normal_3d(&pts);
        assert!(
            normal.norm() < TOLERANCE,
            "expected Newell normal to vanish for collinear points, got {normal:?}",
        );
    }

    #[test]
    fn newell_normal_3d_nonzero_for_proper_polygon() {
        // A unit square in the XY plane — should yield a normal pointing in +Z.
        let pts = vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 0.0),
            p(0.0, 1.0, 0.0),
        ];
        let normal = newell_normal_3d(&pts);
        assert!(
            normal.norm() > TOLERANCE,
            "expected non-degenerate normal, got {normal:?}",
        );
        // Newell sums an oriented area-times-two, so the unit square gives 2.
        assert!((normal.z - 2.0).abs() < 1e-9);
        assert!(normal.x.abs() < 1e-9);
        assert!(normal.y.abs() < 1e-9);
    }

    /// Demonstrates the failure mode the Newell check guards against: a
    /// 4-vertex sliver that lies on a single 3D line can still produce a
    /// non-zero projected area against a mis-oriented plane normal, so the
    /// area-only filter would accept it and `MakeFace::compute_plane_from_points`
    /// would then reject it as "all points collinear". The Newell normal of
    /// the sliver is exactly zero and lets us reject it upstream.
    #[test]
    fn collinear_sliver_caught_by_newell_but_not_by_projected_area() {
        // 4 vertices strictly along the X axis — collinear in 3D.
        let sliver = vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(2.0, 0.0, 0.0),
            p(3.0, 0.0, 0.0),
        ];
        let newell = newell_normal_3d(&sliver);
        assert!(
            newell.norm() < TOLERANCE,
            "Newell normal of a collinear sliver must vanish",
        );

        // MakeFace would reject this sliver — confirm the production guard
        // would also catch it.
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(sliver, true).execute(&mut store).unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store);
        assert!(
            face.is_err(),
            "MakeFace must reject a 4-vertex collinear sliver (got {face:?})",
        );
    }

    fn make_face(store: &mut TopologyStore, points: Vec<Point3>) -> FaceId {
        let wire = MakeWire::new(points, true).execute(store).unwrap();
        MakeFace::new(wire, vec![]).execute(store).unwrap()
    }

    #[test]
    fn no_cuts_returns_whole_face() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![
                p(0.0, 0.0, 0.0),
                p(4.0, 0.0, 0.0),
                p(4.0, 4.0, 0.0),
                p(0.0, 4.0, 0.0),
            ],
        );
        let fragments = split_face(&store, face, &[], SolidSource::A).unwrap();
        assert_eq!(fragments.len(), 1);
        assert_eq!(fragments[0].boundary.len(), 4);
    }

    #[test]
    fn cut_through_middle_produces_two_fragments() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![
                p(0.0, 0.0, 0.0),
                p(4.0, 0.0, 0.0),
                p(4.0, 4.0, 0.0),
                p(0.0, 4.0, 0.0),
            ],
        );
        // Horizontal cut at y=2
        let cuts = vec![(p(0.0, 2.0, 0.0), p(4.0, 2.0, 0.0))];
        let fragments = split_face(&store, face, &cuts, SolidSource::A).unwrap();
        assert_eq!(
            fragments.len(),
            2,
            "expected 2 fragments, got {}",
            fragments.len()
        );

        // Both fragments should be rectangles (4 vertices)
        for (i, f) in fragments.iter().enumerate() {
            assert!(
                f.boundary.len() >= 3,
                "fragment {i} has {} vertices",
                f.boundary.len()
            );
        }
    }

    #[test]
    fn cut_outside_face_returns_one_fragment() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![
                p(0.0, 0.0, 0.0),
                p(4.0, 0.0, 0.0),
                p(4.0, 4.0, 0.0),
                p(0.0, 4.0, 0.0),
            ],
        );
        // Cut far outside the face
        let cuts = vec![(p(10.0, 10.0, 0.0), p(20.0, 10.0, 0.0))];
        let fragments = split_face(&store, face, &cuts, SolidSource::A).unwrap();
        assert_eq!(fragments.len(), 1);
    }

    #[test]
    fn fragment_preserves_source_info() {
        let mut store = TopologyStore::new();
        let face = make_face(
            &mut store,
            vec![
                p(0.0, 0.0, 0.0),
                p(4.0, 0.0, 0.0),
                p(4.0, 4.0, 0.0),
                p(0.0, 4.0, 0.0),
            ],
        );
        let fragments = split_face(&store, face, &[], SolidSource::B).unwrap();
        assert_eq!(fragments[0].source, SolidSource::B);
        assert_eq!(fragments[0].source_face, face);
    }

    fn make_face_with_hole(store: &mut TopologyStore) -> FaceId {
        // outer: 10x10, inner: 2x2 hole at center (4..6, 4..6)
        let outer = MakeWire::new(
            vec![
                p(0.0, 0.0, 0.0),
                p(10.0, 0.0, 0.0),
                p(10.0, 10.0, 0.0),
                p(0.0, 10.0, 0.0),
            ],
            true,
        )
        .execute(store)
        .unwrap();
        let inner = MakeWire::new(
            vec![
                p(4.0, 4.0, 0.0),
                p(6.0, 4.0, 0.0),
                p(6.0, 6.0, 0.0),
                p(4.0, 6.0, 0.0),
            ],
            true,
        )
        .execute(store)
        .unwrap();
        MakeFace::new(outer, vec![inner]).execute(store).unwrap()
    }

    #[test]
    fn split_face_preserves_inner_wire_no_cuts() {
        let mut store = TopologyStore::new();
        let face = make_face_with_hole(&mut store);
        let fragments = split_face(&store, face, &[], SolidSource::A).unwrap();
        assert_eq!(fragments.len(), 1, "no cuts → 1 fragment");
        assert_eq!(
            fragments[0].inner_boundaries.len(),
            1,
            "inner wire must be preserved"
        );
        assert_eq!(
            fragments[0].inner_boundaries[0].len(),
            4,
            "inner wire has 4 vertices"
        );
    }

    #[test]
    fn split_face_inner_in_one_fragment() {
        // outer: 10x10, inner: 2x2 hole at bottom-left (1..3, 1..3)
        // cut at y=5 → hole stays in the lower fragment (y<5), centroid (2,2)
        let mut store = TopologyStore::new();
        let outer = MakeWire::new(
            vec![
                p(0.0, 0.0, 0.0),
                p(10.0, 0.0, 0.0),
                p(10.0, 10.0, 0.0),
                p(0.0, 10.0, 0.0),
            ],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let inner = MakeWire::new(
            vec![
                p(1.0, 1.0, 0.0),
                p(3.0, 1.0, 0.0),
                p(3.0, 3.0, 0.0),
                p(1.0, 3.0, 0.0),
            ],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(outer, vec![inner])
            .execute(&mut store)
            .unwrap();
        let cuts = vec![(p(0.0, 5.0, 0.0), p(10.0, 5.0, 0.0))];
        let fragments = split_face(&store, face, &cuts, SolidSource::A).unwrap();
        assert_eq!(fragments.len(), 2, "cut at y=5 → 2 fragments");
        let total_inner: usize = fragments.iter().map(|f| f.inner_boundaries.len()).sum();
        assert_eq!(total_inner, 1, "hole must appear in exactly one fragment");
        let frag_with_hole = fragments
            .iter()
            .find(|f| !f.inner_boundaries.is_empty())
            .unwrap();
        let centroid_y: f64 = frag_with_hole.inner_boundaries[0]
            .iter()
            .map(|pt| pt.y)
            .sum::<f64>()
            / frag_with_hole.inner_boundaries[0].len() as f64;
        assert!(
            centroid_y < 5.0,
            "hole centroid must be in the lower fragment (y<5), got {centroid_y}"
        );
    }

    #[test]
    fn split_face_splits_inner_wire() {
        // outer: 10x10, inner: 4x6 hole (3..7, 2..8)
        // cut at y=5 bisects the inner wire → each outer fragment gets one inner fragment
        let mut store = TopologyStore::new();
        let outer = MakeWire::new(
            vec![
                p(0.0, 0.0, 0.0),
                p(10.0, 0.0, 0.0),
                p(10.0, 10.0, 0.0),
                p(0.0, 10.0, 0.0),
            ],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let inner = MakeWire::new(
            vec![
                p(3.0, 2.0, 0.0),
                p(7.0, 2.0, 0.0),
                p(7.0, 8.0, 0.0),
                p(3.0, 8.0, 0.0),
            ],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(outer, vec![inner])
            .execute(&mut store)
            .unwrap();
        let cuts = vec![(p(0.0, 5.0, 0.0), p(10.0, 5.0, 0.0))];
        let fragments = split_face(&store, face, &cuts, SolidSource::A).unwrap();
        assert_eq!(fragments.len(), 2, "cut at y=5 → 2 fragments");
        let total_inner: usize = fragments.iter().map(|f| f.inner_boundaries.len()).sum();
        assert_eq!(
            total_inner, 2,
            "each fragment should carry one inner fragment"
        );
        for (i, frag) in fragments.iter().enumerate() {
            assert_eq!(
                frag.inner_boundaries.len(),
                1,
                "fragment {i} should have exactly 1 inner boundary"
            );
        }
    }
}
