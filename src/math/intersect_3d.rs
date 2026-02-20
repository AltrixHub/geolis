use crate::geometry::surface::Plane;

use super::{Point3, Vector3, TOLERANCE};

/// Relationship between two planes.
#[derive(Debug)]
pub enum PlanePairRelation {
    /// Planes intersect along a line.
    IntersectionLine {
        origin: Point3,
        direction: Vector3,
    },
    /// Planes are parallel but not coincident.
    Parallel { distance: f64 },
    /// Planes are the same (coincident).
    Coincident,
}

/// Computes the intersection of two planes.
///
/// Returns an [`IntersectionLine`](PlanePairRelation::IntersectionLine) with a
/// unit-length `direction` when the planes cross, [`Parallel`](PlanePairRelation::Parallel)
/// when they don't, or [`Coincident`](PlanePairRelation::Coincident) when they overlap.
#[must_use]
pub fn plane_plane_intersect(a: &Plane, b: &Plane) -> PlanePairRelation {
    let na = a.plane_normal();
    let nb = b.plane_normal();

    let dir = na.cross(nb);
    let dir_len = dir.norm();

    if dir_len < TOLERANCE {
        // Normals are (anti-)parallel — planes are parallel or coincident.
        let diff = b.origin() - a.origin();
        let dist = diff.dot(na).abs();
        if dist < TOLERANCE {
            PlanePairRelation::Coincident
        } else {
            PlanePairRelation::Parallel { distance: dist }
        }
    } else {
        let dir = dir / dir_len;

        // Find a point on the intersection line.
        // Solve na.dot(p - oa) = 0 AND nb.dot(p - ob) = 0 simultaneously.
        // Write p = oa + s * na + t * nb + u * dir  (u is free, set u = 0).
        // na.dot(s*na + t*nb + oa - oa) = 0  =>  s + t*(na.nb) = 0
        // nb.dot(s*na + t*nb + oa - ob) = 0  =>  s*(na.nb) + t = nb.dot(ob - oa)
        let d1 = 0.0_f64; // na.dot(oa - oa)
        let d2 = nb.dot(&(b.origin() - a.origin()));
        let dot_nn = na.dot(nb);
        let denom = 1.0 - dot_nn * dot_nn;

        let origin = if denom.abs() < TOLERANCE {
            // Fallback: project a.origin() onto intersection line
            *a.origin()
        } else {
            let s = (d1 - dot_nn * d2) / denom; // d1 is 0
            let t = (d2 - dot_nn * d1) / denom;
            a.origin() + na * s + nb * t
        };

        PlanePairRelation::IntersectionLine { origin, direction: dir }
    }
}

/// Relationship of a line with a plane.
#[derive(Debug)]
pub enum LinePlaneRelation {
    /// Line intersects the plane at a single point.
    Point { point: Point3, t: f64 },
    /// Line is parallel to the plane (does not intersect).
    Parallel,
    /// Line lies entirely on the plane.
    OnPlane,
}

/// Computes the intersection of a line `origin + t * dir` with a plane.
#[must_use]
pub fn line_plane_intersect(
    origin: &Point3,
    dir: &Vector3,
    plane: &Plane,
) -> LinePlaneRelation {
    let normal = plane.plane_normal();
    let denom = normal.dot(dir);

    let diff = plane.origin() - origin;
    let numer = normal.dot(&diff);

    if denom.abs() < TOLERANCE {
        // Line is parallel to the plane
        if numer.abs() < TOLERANCE {
            LinePlaneRelation::OnPlane
        } else {
            LinePlaneRelation::Parallel
        }
    } else {
        let t = numer / denom;
        let point = origin + dir * t;
        LinePlaneRelation::Point { point, t }
    }
}

/// Classification of a point relative to a plane.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointPlaneClassification {
    /// Point is on the positive side (in the direction of the normal).
    Front,
    /// Point is on the negative side (opposite the normal).
    Back,
    /// Point lies on the plane (within tolerance).
    On,
}

/// Classifies a point relative to a plane.
#[must_use]
pub fn classify_point_plane(point: &Point3, plane: &Plane) -> PointPlaneClassification {
    let diff = point - plane.origin();
    let dist = plane.plane_normal().dot(&diff);

    if dist > TOLERANCE {
        PointPlaneClassification::Front
    } else if dist < -TOLERANCE {
        PointPlaneClassification::Back
    } else {
        PointPlaneClassification::On
    }
}

/// Signed distance from a point to a plane.
/// Positive = on the normal side, negative = opposite.
#[must_use]
pub fn signed_distance_to_plane(point: &Point3, plane: &Plane) -> f64 {
    let diff = point - plane.origin();
    plane.plane_normal().dot(&diff)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn v(x: f64, y: f64, z: f64) -> Vector3 {
        Vector3::new(x, y, z)
    }

    // ── plane_plane_intersect ──

    #[test]
    fn perpendicular_planes_intersect() {
        // XY-plane and XZ-plane should intersect along the X-axis
        let xy = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        let xz = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 1.0, 0.0)).unwrap();

        let result = plane_plane_intersect(&xy, &xz);
        match result {
            PlanePairRelation::IntersectionLine { direction, .. } => {
                // Direction should be along X-axis (±)
                assert!(
                    direction.x.abs() > 0.99,
                    "expected X-axis direction, got {direction:?}"
                );
            }
            other => panic!("expected IntersectionLine, got {other:?}"),
        }
    }

    #[test]
    fn oblique_planes_intersect() {
        // XY-plane and a tilted plane
        let xy = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        let tilted = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 1.0, 1.0)).unwrap();

        let result = plane_plane_intersect(&xy, &tilted);
        match result {
            PlanePairRelation::IntersectionLine { direction, .. } => {
                // Intersection should be along X-axis
                assert!(
                    direction.x.abs() > 0.99,
                    "expected X-axis direction, got {direction:?}"
                );
            }
            other => panic!("expected IntersectionLine, got {other:?}"),
        }
    }

    #[test]
    fn parallel_planes() {
        let a = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        let b = Plane::from_normal(p(0.0, 0.0, 5.0), v(0.0, 0.0, 1.0)).unwrap();

        let result = plane_plane_intersect(&a, &b);
        match result {
            PlanePairRelation::Parallel { distance } => {
                assert!((distance - 5.0).abs() < TOLERANCE);
            }
            other => panic!("expected Parallel, got {other:?}"),
        }
    }

    #[test]
    fn coincident_planes() {
        let a = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        let b = Plane::from_normal(p(1.0, 2.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();

        assert!(matches!(
            plane_plane_intersect(&a, &b),
            PlanePairRelation::Coincident
        ));
    }

    #[test]
    fn anti_parallel_planes_are_parallel() {
        let a = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        let b = Plane::from_normal(p(0.0, 0.0, 3.0), v(0.0, 0.0, -1.0)).unwrap();

        let result = plane_plane_intersect(&a, &b);
        match result {
            PlanePairRelation::Parallel { distance } => {
                assert!((distance - 3.0).abs() < TOLERANCE);
            }
            other => panic!("expected Parallel, got {other:?}"),
        }
    }

    #[test]
    fn intersection_point_lies_on_both_planes() {
        let a = Plane::from_normal(p(1.0, 0.0, 0.0), v(1.0, 0.0, 0.0)).unwrap();
        let b = Plane::from_normal(p(0.0, 2.0, 0.0), v(0.0, 1.0, 0.0)).unwrap();

        match plane_plane_intersect(&a, &b) {
            PlanePairRelation::IntersectionLine { origin, direction } => {
                // origin should lie on both planes
                let dist_a = (origin - p(1.0, 0.0, 0.0)).dot(&v(1.0, 0.0, 0.0));
                let dist_b = (origin - p(0.0, 2.0, 0.0)).dot(&v(0.0, 1.0, 0.0));
                assert!(
                    dist_a.abs() < TOLERANCE,
                    "origin not on plane A: dist = {dist_a}"
                );
                assert!(
                    dist_b.abs() < TOLERANCE,
                    "origin not on plane B: dist = {dist_b}"
                );
                // direction should be along Z-axis
                assert!(direction.z.abs() > 0.99);
            }
            other => panic!("expected IntersectionLine, got {other:?}"),
        }
    }

    // ── line_plane_intersect ──

    #[test]
    fn line_hits_plane() {
        let plane = Plane::from_normal(p(0.0, 0.0, 5.0), v(0.0, 0.0, 1.0)).unwrap();
        let result = line_plane_intersect(&p(0.0, 0.0, 0.0), &v(0.0, 0.0, 1.0), &plane);
        match result {
            LinePlaneRelation::Point { point, t } => {
                assert!((t - 5.0).abs() < TOLERANCE);
                assert!((point.z - 5.0).abs() < TOLERANCE);
            }
            other => panic!("expected Point, got {other:?}"),
        }
    }

    #[test]
    fn line_parallel_to_plane() {
        let plane = Plane::from_normal(p(0.0, 0.0, 5.0), v(0.0, 0.0, 1.0)).unwrap();
        let result = line_plane_intersect(&p(0.0, 0.0, 0.0), &v(1.0, 0.0, 0.0), &plane);
        assert!(matches!(result, LinePlaneRelation::Parallel));
    }

    #[test]
    fn line_on_plane() {
        let plane = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        let result = line_plane_intersect(&p(1.0, 2.0, 0.0), &v(1.0, 0.0, 0.0), &plane);
        assert!(matches!(result, LinePlaneRelation::OnPlane));
    }

    #[test]
    fn line_oblique_to_plane() {
        let plane = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        let result =
            line_plane_intersect(&p(0.0, 0.0, -3.0), &v(1.0, 1.0, 1.0), &plane);
        match result {
            LinePlaneRelation::Point { point, t } => {
                assert!((t - 3.0).abs() < TOLERANCE);
                assert!((point.z).abs() < TOLERANCE);
            }
            other => panic!("expected Point, got {other:?}"),
        }
    }

    // ── classify_point_plane ──

    #[test]
    fn point_in_front_of_plane() {
        let plane = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        assert_eq!(
            classify_point_plane(&p(0.0, 0.0, 1.0), &plane),
            PointPlaneClassification::Front
        );
    }

    #[test]
    fn point_behind_plane() {
        let plane = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        assert_eq!(
            classify_point_plane(&p(0.0, 0.0, -1.0), &plane),
            PointPlaneClassification::Back
        );
    }

    #[test]
    fn point_on_plane() {
        let plane = Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap();
        assert_eq!(
            classify_point_plane(&p(5.0, 3.0, 0.0), &plane),
            PointPlaneClassification::On
        );
    }
}
