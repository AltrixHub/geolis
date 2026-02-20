use crate::geometry::surface::Plane;

use super::{Point3, Vector3, TOLERANCE};

/// Projects a 3D point onto the UV coordinate system of a plane.
///
/// Returns `(u, v)` coordinates.
#[must_use]
fn project_to_uv(point: &Point3, plane: &Plane) -> (f64, f64) {
    let diff = point - plane.origin();
    let u = diff.dot(plane.u_dir());
    let v = diff.dot(plane.v_dir());
    (u, v)
}

/// Point-in-polygon test for a 3D point coplanar with the polygon.
///
/// Projects to the face's UV coordinate space and uses the winding number
/// algorithm. Returns `true` if the point is inside or on the boundary.
#[must_use]
pub fn point_in_polygon_3d(point: &Point3, polygon: &[Point3], plane: &Plane) -> bool {
    if polygon.len() < 3 {
        return false;
    }

    let (px, py) = project_to_uv(point, plane);
    let uvs: Vec<(f64, f64)> = polygon.iter().map(|p| project_to_uv(p, plane)).collect();

    winding_number_2d(px, py, &uvs) != 0
}

/// Winding number of point `(px, py)` with respect to polygon `verts`.
///
/// Non-zero => inside, zero => outside.
fn winding_number_2d(px: f64, py: f64, verts: &[(f64, f64)]) -> i32 {
    let n = verts.len();
    let mut winding = 0i32;
    for i in 0..n {
        let (x0, y0) = verts[i];
        let (x1, y1) = verts[(i + 1) % n];

        if y0 <= py {
            if y1 > py && cross_2d(x1 - x0, y1 - y0, px - x0, py - y0) > 0.0 {
                winding += 1;
            }
        } else if y1 <= py && cross_2d(x1 - x0, y1 - y0, px - x0, py - y0) < 0.0 {
            winding -= 1;
        }
    }
    winding
}

/// 2D cross product: `(ax * by - ay * bx)`.
#[inline]
fn cross_2d(ax: f64, ay: f64, bx: f64, by: f64) -> f64 {
    ax * by - ay * bx
}

/// Clips a line segment `(seg_start, seg_end)` to a convex polygon boundary.
///
/// Both the segment and polygon must be coplanar with the given plane.
/// Returns the sub-segments that lie inside the polygon as `(t_start, t_end)`
/// pairs where `t` is the parameter along the original segment `[0, 1]`.
///
/// For non-convex polygons, multiple sub-segments may be returned.
#[must_use]
pub fn clip_segment_to_polygon(
    seg_start: &Point3,
    seg_end: &Point3,
    polygon: &[Point3],
    plane: &Plane,
) -> Vec<(f64, f64)> {
    if polygon.len() < 3 {
        return Vec::new();
    }

    let (su, sv) = project_to_uv(seg_start, plane);
    let (eu, ev) = project_to_uv(seg_end, plane);
    let du = eu - su;
    let dv = ev - sv;

    let uvs: Vec<(f64, f64)> = polygon.iter().map(|p| project_to_uv(p, plane)).collect();

    // Collect all t-values where the segment crosses a polygon edge
    let n = uvs.len();
    let mut crossings: Vec<f64> = Vec::new();

    for i in 0..n {
        let (ex0, ey0) = uvs[i];
        let (ex1, ey1) = uvs[(i + 1) % n];
        let edx = ex1 - ex0;
        let edy = ey1 - ey0;

        let cross = du * edy - dv * edx;
        if cross.abs() < TOLERANCE {
            continue; // Parallel
        }

        let dx = ex0 - su;
        let dy = ey0 - sv;
        let t = (dx * edy - dy * edx) / cross;
        let u_edge = (dx * dv - dy * du) / cross;

        let eps = TOLERANCE;
        if t >= -eps && t <= 1.0 + eps && u_edge >= -eps && u_edge <= 1.0 + eps {
            crossings.push(t.clamp(0.0, 1.0));
        }
    }

    crossings.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    crossings.dedup_by(|a, b| (*a - *b).abs() < TOLERANCE);

    // Build candidate intervals from crossings + endpoints
    let mut sample_ts = vec![0.0];
    sample_ts.extend_from_slice(&crossings);
    sample_ts.push(1.0);
    sample_ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    sample_ts.dedup_by(|a, b| (*a - *b).abs() < TOLERANCE);

    // For each interval between consecutive t-values, test midpoint containment
    let mut result = Vec::new();
    for win in sample_ts.windows(2) {
        let t0 = win[0];
        let t1 = win[1];
        if (t1 - t0).abs() < TOLERANCE {
            continue;
        }
        let mid_t = (t0 + t1) * 0.5;
        let mid_u = su + du * mid_t;
        let mid_v = sv + dv * mid_t;
        if winding_number_2d(mid_u, mid_v, &uvs) != 0 {
            // Merge with previous interval if contiguous
            if let Some(last) = result.last_mut() {
                let (_, ref mut last_end): (f64, f64) = *last;
                if (t0 - *last_end).abs() < TOLERANCE {
                    *last_end = t1;
                    continue;
                }
            }
            result.push((t0, t1));
        }
    }

    result
}

/// Compute the 3D point along a segment at parameter `t`.
#[must_use]
pub fn segment_point_at(start: &Point3, end: &Point3, t: f64) -> Point3 {
    let dir = end - start;
    start + dir * t
}

/// Compute the area of a 3D polygon (coplanar points).
///
/// Uses the cross-product summation method projected along the polygon normal.
#[must_use]
pub fn polygon_area_3d(points: &[Point3], normal: &Vector3) -> f64 {
    if points.len() < 3 {
        return 0.0;
    }
    let n = points.len();
    let mut cross_sum = Vector3::new(0.0, 0.0, 0.0);
    let o = &points[0];
    for i in 1..n {
        let a = points[i] - o;
        let b = points[(i + 1) % n] - o;
        cross_sum += a.cross(&b);
    }
    0.5 * cross_sum.dot(normal).abs()
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

    fn xy_plane() -> Plane {
        Plane::from_normal(p(0.0, 0.0, 0.0), v(0.0, 0.0, 1.0)).unwrap()
    }

    fn unit_square() -> Vec<Point3> {
        vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 0.0),
            p(0.0, 1.0, 0.0),
        ]
    }

    // ── point_in_polygon_3d ──

    #[test]
    fn point_inside_square() {
        let plane = xy_plane();
        assert!(point_in_polygon_3d(&p(0.5, 0.5, 0.0), &unit_square(), &plane));
    }

    #[test]
    fn point_outside_square() {
        let plane = xy_plane();
        assert!(!point_in_polygon_3d(&p(2.0, 0.5, 0.0), &unit_square(), &plane));
    }

    #[test]
    fn point_inside_triangle() {
        let plane = xy_plane();
        let tri = vec![p(0.0, 0.0, 0.0), p(4.0, 0.0, 0.0), p(2.0, 3.0, 0.0)];
        assert!(point_in_polygon_3d(&p(2.0, 1.0, 0.0), &tri, &plane));
    }

    #[test]
    fn point_outside_triangle() {
        let plane = xy_plane();
        let tri = vec![p(0.0, 0.0, 0.0), p(4.0, 0.0, 0.0), p(2.0, 3.0, 0.0)];
        assert!(!point_in_polygon_3d(&p(5.0, 5.0, 0.0), &tri, &plane));
    }

    // ── clip_segment_to_polygon ──

    #[test]
    fn segment_fully_inside() {
        let plane = xy_plane();
        let sq = unit_square();
        let result = clip_segment_to_polygon(
            &p(0.2, 0.5, 0.0),
            &p(0.8, 0.5, 0.0),
            &sq,
            &plane,
        );
        assert_eq!(result.len(), 1);
        assert!((result[0].0).abs() < TOLERANCE);
        assert!((result[0].1 - 1.0).abs() < TOLERANCE);
    }

    #[test]
    fn segment_fully_outside() {
        let plane = xy_plane();
        let sq = unit_square();
        let result = clip_segment_to_polygon(
            &p(2.0, 0.5, 0.0),
            &p(3.0, 0.5, 0.0),
            &sq,
            &plane,
        );
        assert!(result.is_empty());
    }

    #[test]
    fn segment_partial_clip() {
        let plane = xy_plane();
        let sq = unit_square();
        // Segment from outside to inside
        let result = clip_segment_to_polygon(
            &p(-0.5, 0.5, 0.0),
            &p(0.5, 0.5, 0.0),
            &sq,
            &plane,
        );
        assert_eq!(result.len(), 1);
        // t=0.5 is where it enters the square (x=0.0)
        assert!((result[0].0 - 0.5).abs() < 0.01);
        assert!((result[0].1 - 1.0).abs() < 0.01);
    }

    #[test]
    fn segment_through_polygon() {
        let plane = xy_plane();
        let sq = unit_square();
        // Segment goes clean through the square from left to right
        let result = clip_segment_to_polygon(
            &p(-1.0, 0.5, 0.0),
            &p(2.0, 0.5, 0.0),
            &sq,
            &plane,
        );
        assert_eq!(result.len(), 1);
        // Enters at t = 1/3 (x=0), exits at t = 2/3 (x=1)
        let expected_t0 = 1.0 / 3.0;
        let expected_t1 = 2.0 / 3.0;
        assert!(
            (result[0].0 - expected_t0).abs() < 0.01,
            "t0 = {}, expected {}",
            result[0].0,
            expected_t0
        );
        assert!(
            (result[0].1 - expected_t1).abs() < 0.01,
            "t1 = {}, expected {}",
            result[0].1,
            expected_t1
        );
    }

    // ── polygon_area_3d ──

    #[test]
    fn unit_square_area() {
        let area = polygon_area_3d(&unit_square(), &v(0.0, 0.0, 1.0));
        assert!((area - 1.0).abs() < TOLERANCE);
    }

    #[test]
    fn triangle_area() {
        let tri = vec![p(0.0, 0.0, 0.0), p(4.0, 0.0, 0.0), p(0.0, 3.0, 0.0)];
        let area = polygon_area_3d(&tri, &v(0.0, 0.0, 1.0));
        assert!((area - 6.0).abs() < TOLERANCE);
    }
}
