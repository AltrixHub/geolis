use crate::error::Result;
use crate::geometry::surface::Plane;
use crate::math::polygon_3d::polygon_area_3d;
use crate::math::{Point3, TOLERANCE};
use crate::topology::{FaceId, FaceSurface, TopologyStore};

use super::face_intersection::collect_face_polygon;

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
    let FaceSurface::Plane(ref plane) = face.surface;

    let polygon = collect_face_polygon(store, face_id)?;
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

    // Filter out degenerate fragments
    let normal = plane.plane_normal();
    let min_area = TOLERANCE * TOLERANCE;
    let result = fragments
        .into_iter()
        .filter(|f| f.len() >= 3 && polygon_area_3d(f, normal) > min_area)
        .map(|boundary| FaceFragment {
            boundary,
            plane: plane.clone(),
            same_sense: face.same_sense,
            source_face: face_id,
            source,
        })
        .collect();

    Ok(result)
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

    let unproject = |u: f64, v: f64| -> Point3 {
        origin + u_dir * u + v_dir * v
    };

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
        assert_eq!(fragments.len(), 2, "expected 2 fragments, got {}", fragments.len());

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
}
