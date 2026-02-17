use std::f64::consts::PI;

use super::{Point3, Vector3, TOLERANCE};

/// Parametric 2D line-line intersection.
///
/// Given lines `p1 + t * d1` and `p2 + u * d2`, returns `(t, u)` if not parallel.
#[must_use]
pub fn line_line_intersect_2d(
    p1: &Point3,
    d1: &Vector3,
    p2: &Point3,
    d2: &Vector3,
) -> Option<(f64, f64)> {
    let cross = d1.x * d2.y - d1.y * d2.x;
    if cross.abs() < TOLERANCE {
        return None;
    }
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let t = (dx * d2.y - dy * d2.x) / cross;
    let u = (dx * d1.y - dy * d1.x) / cross;
    Some((t, u))
}

/// Bounded segment-segment intersection in 2D.
///
/// Returns `(intersection_point, t, u)` where `t` and `u` are in `[0, 1]`.
#[must_use]
pub fn segment_segment_intersect_2d(
    a0: &Point3,
    a1: &Point3,
    b0: &Point3,
    b1: &Point3,
) -> Option<(Point3, f64, f64)> {
    let da = Vector3::new(a1.x - a0.x, a1.y - a0.y, 0.0);
    let db = Vector3::new(b1.x - b0.x, b1.y - b0.y, 0.0);

    let cross = da.x * db.y - da.y * db.x;
    if cross.abs() < TOLERANCE {
        return None;
    }

    let dx = b0.x - a0.x;
    let dy = b0.y - a0.y;
    let t = (dx * db.y - dy * db.x) / cross;
    let u = (dx * da.y - dy * da.x) / cross;

    // Use a small epsilon to include endpoints.
    let eps = TOLERANCE;
    if t >= -eps && t <= 1.0 + eps && u >= -eps && u <= 1.0 + eps {
        let t_clamped = t.clamp(0.0, 1.0);
        let pt = Point3::new(a0.x + da.x * t_clamped, a0.y + da.y * t_clamped, a0.z);
        Some((pt, t_clamped, u.clamp(0.0, 1.0)))
    } else {
        None
    }
}

/// Linear interpolation: `origin + dir * t`.
#[must_use]
pub fn point_at(origin: &Point3, dir: &Vector3, t: f64) -> Point3 {
    Point3::new(origin.x + dir.x * t, origin.y + dir.y * t, origin.z)
}

/// Intersection of a line segment with a circular arc in 2D.
///
/// The segment goes from `(ax0, ay0)` to `(ax1, ay1)`.
/// The arc has center `(cx, cy)`, `radius`, `start_angle`, and `sweep`.
///
/// Returns a vector of `((x, y), t_seg, t_arc)` where:
/// - `t_seg` is the parameter on the segment `[0, 1]`
/// - `t_arc` is the parameter on the arc `[0, 1]`
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn line_arc_intersect_2d(
    ax0: f64,
    ay0: f64,
    ax1: f64,
    ay1: f64,
    cx: f64,
    cy: f64,
    radius: f64,
    start_angle: f64,
    sweep: f64,
) -> Vec<((f64, f64), f64, f64)> {
    let mut results = Vec::new();
    if radius < TOLERANCE || sweep.abs() < TOLERANCE {
        return results;
    }

    let dx = ax1 - ax0;
    let dy = ay1 - ay0;
    let seg_len_sq = dx * dx + dy * dy;
    if seg_len_sq < TOLERANCE * TOLERANCE {
        return results;
    }

    // Substitute parametric line into circle equation:
    // (ax0 + t*dx - cx)² + (ay0 + t*dy - cy)² = r²
    let fx = ax0 - cx;
    let fy = ay0 - cy;
    let a = seg_len_sq;
    let b = 2.0 * (fx * dx + fy * dy);
    let c = fx * fx + fy * fy - radius * radius;
    let discriminant = b * b - 4.0 * a * c;

    if discriminant < -TOLERANCE {
        return results;
    }
    let disc_sqrt = discriminant.max(0.0).sqrt();

    let eps = TOLERANCE;
    let t_roots = if disc_sqrt < TOLERANCE * 100.0 {
        // Tangent case: single root.
        vec![-b / (2.0 * a)]
    } else {
        vec![(-b - disc_sqrt) / (2.0 * a), (-b + disc_sqrt) / (2.0 * a)]
    };

    for t_seg in t_roots {
        if t_seg < -eps || t_seg > 1.0 + eps {
            continue;
        }
        let t_seg = t_seg.clamp(0.0, 1.0);

        let px = ax0 + t_seg * dx;
        let py = ay0 + t_seg * dy;

        // Check if point is within the arc's angular range.
        let angle = (py - cy).atan2(px - cx);
        if let Some(t_arc) = angle_to_arc_param(angle, start_angle, sweep) {
            results.push(((px, py), t_seg, t_arc));
        }
    }

    results
}

/// Intersection of two circular arcs in 2D.
///
/// Arc 1: center `(c1x, c1y)`, `r1`, `start1`, `sweep1`.
/// Arc 2: center `(c2x, c2y)`, `r2`, `start2`, `sweep2`.
///
/// Returns a vector of `((x, y), t1, t2)` where `t1` and `t2` are arc parameters in `[0, 1]`.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn arc_arc_intersect_2d(
    c1x: f64,
    c1y: f64,
    r1: f64,
    start1: f64,
    sweep1: f64,
    c2x: f64,
    c2y: f64,
    r2: f64,
    start2: f64,
    sweep2: f64,
) -> Vec<((f64, f64), f64, f64)> {
    let mut results = Vec::new();
    if r1 < TOLERANCE || r2 < TOLERANCE {
        return results;
    }

    let dx = c2x - c1x;
    let dy = c2y - c1y;
    let dist_sq = dx * dx + dy * dy;
    let dist = dist_sq.sqrt();

    if dist < TOLERANCE {
        // Concentric circles — no intersection points (or infinite if same radius).
        return results;
    }

    // Check if circles intersect.
    let sum = r1 + r2;
    let diff = (r1 - r2).abs();
    if dist > sum + TOLERANCE || dist < diff - TOLERANCE {
        return results;
    }

    // Distance from c1 along the line c1→c2 to the radical line.
    let a = (r1 * r1 - r2 * r2 + dist_sq) / (2.0 * dist);
    let h_sq = r1 * r1 - a * a;
    if h_sq < -TOLERANCE {
        return results;
    }
    let h = h_sq.max(0.0).sqrt();

    // Midpoint on the radical line.
    let mx = c1x + a * dx / dist;
    let my = c1y + a * dy / dist;

    // Perpendicular direction.
    let px = -dy / dist;
    let py = dx / dist;

    // Two candidate intersection points (or one if tangent).
    let candidates = if h < TOLERANCE {
        vec![(mx, my)]
    } else {
        vec![(mx + h * px, my + h * py), (mx - h * px, my - h * py)]
    };

    let eps = TOLERANCE;
    for (ix, iy) in candidates {
        let angle1 = (iy - c1y).atan2(ix - c1x);
        let angle2 = (iy - c2y).atan2(ix - c2x);

        let t1 = angle_to_arc_param(angle1, start1, sweep1);
        let t2 = angle_to_arc_param(angle2, start2, sweep2);

        if let (Some(t1), Some(t2)) = (t1, t2) {
            // Verify the point is close to both arcs.
            let d1 = ((ix - c1x).powi(2) + (iy - c1y).powi(2)).sqrt();
            let d2 = ((ix - c2x).powi(2) + (iy - c2y).powi(2)).sqrt();
            if (d1 - r1).abs() < eps && (d2 - r2).abs() < eps {
                results.push(((ix, iy), t1, t2));
            }
        }
    }

    results
}

/// Converts an absolute angle to an arc parameter `t` in `[0, 1]`.
///
/// Returns `None` if the angle is not within the arc's angular range.
fn angle_to_arc_param(angle: f64, start_angle: f64, sweep: f64) -> Option<f64> {
    let eps = TOLERANCE * 100.0;

    // Compute the angular offset from start_angle to angle in the sweep direction.
    let mut delta = angle - start_angle;

    // Normalize delta to match the sweep direction.
    if sweep > 0.0 {
        while delta < -eps {
            delta += 2.0 * PI;
        }
        while delta > 2.0 * PI + eps {
            delta -= 2.0 * PI;
        }
    } else {
        while delta > eps {
            delta -= 2.0 * PI;
        }
        while delta < -2.0 * PI - eps {
            delta += 2.0 * PI;
        }
    }

    let t = delta / sweep;
    if t >= -eps && t <= 1.0 + eps {
        Some(t.clamp(0.0, 1.0))
    } else {
        None
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn line_line_perpendicular() {
        let p1 = Point3::new(0.0, 0.0, 0.0);
        let d1 = Vector3::new(1.0, 0.0, 0.0);
        let p2 = Point3::new(0.5, -1.0, 0.0);
        let d2 = Vector3::new(0.0, 1.0, 0.0);
        let (t, u) = line_line_intersect_2d(&p1, &d1, &p2, &d2).unwrap();
        assert!((t - 0.5).abs() < TOLERANCE);
        assert!((u - 1.0).abs() < TOLERANCE);
    }

    #[test]
    fn line_line_parallel_returns_none() {
        let p1 = Point3::new(0.0, 0.0, 0.0);
        let d1 = Vector3::new(1.0, 0.0, 0.0);
        let p2 = Point3::new(0.0, 1.0, 0.0);
        let d2 = Vector3::new(1.0, 0.0, 0.0);
        assert!(line_line_intersect_2d(&p1, &d1, &p2, &d2).is_none());
    }

    #[test]
    fn segment_segment_crossing() {
        let a0 = Point3::new(0.0, 0.0, 0.0);
        let a1 = Point3::new(2.0, 2.0, 0.0);
        let b0 = Point3::new(0.0, 2.0, 0.0);
        let b1 = Point3::new(2.0, 0.0, 0.0);
        let (pt, t, u) = segment_segment_intersect_2d(&a0, &a1, &b0, &b1).unwrap();
        assert!((pt.x - 1.0).abs() < TOLERANCE);
        assert!((pt.y - 1.0).abs() < TOLERANCE);
        assert!((t - 0.5).abs() < TOLERANCE);
        assert!((u - 0.5).abs() < TOLERANCE);
    }

    #[test]
    fn segment_segment_no_crossing() {
        let a0 = Point3::new(0.0, 0.0, 0.0);
        let a1 = Point3::new(1.0, 0.0, 0.0);
        let b0 = Point3::new(0.0, 1.0, 0.0);
        let b1 = Point3::new(1.0, 1.0, 0.0);
        assert!(segment_segment_intersect_2d(&a0, &a1, &b0, &b1).is_none());
    }

    #[test]
    fn point_at_interpolation() {
        let origin = Point3::new(1.0, 2.0, 3.0);
        let dir = Vector3::new(4.0, 6.0, 0.0);
        let pt = point_at(&origin, &dir, 0.5);
        assert!((pt.x - 3.0).abs() < TOLERANCE);
        assert!((pt.y - 5.0).abs() < TOLERANCE);
        assert!((pt.z - 3.0).abs() < TOLERANCE);
    }

    // ── line-arc intersection tests ──

    #[test]
    fn line_arc_two_crossings() {
        // Horizontal segment through unit circle at y=0.
        // Arc: full semicircle from angle 0 to π (CCW), center at origin, radius 1.
        let hits = line_arc_intersect_2d(
            -2.0, 0.0, 2.0, 0.0, // segment
            0.0, 0.0, 1.0,        // center, radius
            0.0, PI,               // start_angle=0, sweep=π
        );
        // Should hit at (1, 0) (t_arc=0) and (-1, 0) (t_arc=1).
        assert_eq!(hits.len(), 2, "expected 2 hits, got {}", hits.len());
    }

    #[test]
    fn line_arc_no_crossing() {
        // Segment entirely outside the circle.
        let hits = line_arc_intersect_2d(
            3.0, 0.0, 4.0, 0.0,
            0.0, 0.0, 1.0,
            0.0, PI,
        );
        assert!(hits.is_empty());
    }

    #[test]
    fn line_arc_tangent() {
        // Horizontal segment tangent to unit circle at (0, 1).
        // Arc: upper half, from π/2 to -π/2 (CW through top).
        // Actually use CCW arc from 0 to π which goes through π/2.
        let hits = line_arc_intersect_2d(
            -1.0, 1.0, 1.0, 1.0,  // horizontal line at y=1
            0.0, 0.0, 1.0,        // center, radius
            0.0, PI,               // CCW semicircle
        );
        // Should hit at (0, 1) which is at angle π/2 (t_arc = 0.5).
        assert_eq!(hits.len(), 1, "hits={hits:?}");
        assert!((hits[0].0 .0).abs() < 1e-6, "x={}", hits[0].0 .0);
        assert!((hits[0].0 .1 - 1.0).abs() < 1e-6, "y={}", hits[0].0 .1);
    }

    #[test]
    fn line_arc_miss_outside_arc_range() {
        // Segment crosses the circle but NOT within the arc's angular range.
        // Arc from 0 to π/2 (first quadrant only).
        let hits = line_arc_intersect_2d(
            -2.0, 0.0, 2.0, 0.0,  // horizontal at y=0
            0.0, 0.0, 1.0,
            PI / 4.0, PI / 4.0,   // arc from 45° to 90° only
        );
        // Circle intersects at x=±1, y=0. Angles 0 and π.
        // Neither is in range [π/4, π/2]. Should return empty.
        assert!(hits.is_empty(), "hits={hits:?}");
    }

    // ── arc-arc intersection tests ──

    #[test]
    fn arc_arc_two_crossings() {
        // Two unit circles, centers at (0,0) and (1,0).
        // Intersection points at (0.5, ±√3/2) ≈ (0.5, ±0.866).
        // Use large arcs that cover both intersection angles.
        // Arc1: center (0,0), covers [-π, π] (full circle sweep of 2π from -π)
        // Arc2: center (1,0), covers [0, 2π] (full circle)
        let hits = arc_arc_intersect_2d(
            0.0, 0.0, 1.0, -PI, 2.0 * PI,       // arc1: nearly full circle
            1.0, 0.0, 1.0, 0.0, 2.0 * PI,        // arc2: nearly full circle
        );
        assert_eq!(hits.len(), 2, "hits={hits:?}");
        // Verify intersection points.
        let sqrt3_2 = 3.0_f64.sqrt() / 2.0;
        let (mut y0, mut y1) = (hits[0].0 .1, hits[1].0 .1);
        if y0 > y1 {
            std::mem::swap(&mut y0, &mut y1);
        }
        assert!((y0 + sqrt3_2).abs() < 1e-6, "y0={y0}");
        assert!((y1 - sqrt3_2).abs() < 1e-6, "y1={y1}");
    }

    #[test]
    fn arc_arc_no_overlap() {
        // Two circles too far apart.
        let hits = arc_arc_intersect_2d(
            0.0, 0.0, 1.0, 0.0, PI,
            5.0, 0.0, 1.0, 0.0, PI,
        );
        assert!(hits.is_empty());
    }

    #[test]
    fn arc_arc_tangent() {
        // Two unit circles tangent externally at (1, 0).
        // Arc1 covers angle 0. Arc2 covers angle π.
        let hits = arc_arc_intersect_2d(
            0.0, 0.0, 1.0, -PI / 4.0, PI / 2.0,   // covers angle 0
            2.0, 0.0, 1.0, PI / 2.0, PI,            // covers angle π
        );
        assert_eq!(hits.len(), 1, "hits={hits:?}");
        assert!((hits[0].0 .0 - 1.0).abs() < 1e-6);
        assert!((hits[0].0 .1).abs() < 1e-6);
    }

    #[test]
    fn arc_arc_miss_outside_range() {
        // Circles overlap, but arcs don't cover the intersection angles.
        let hits = arc_arc_intersect_2d(
            0.0, 0.0, 1.0, 0.0, PI / 4.0,         // small arc near angle 0
            1.0, 0.0, 1.0, PI, PI / 4.0,            // small arc near angle π
        );
        // The intersection points of these circles are at y ≈ ±0.866, angles ≈ ±60°.
        // Arc1 only covers [0°, 45°], arc2 only covers [180°, 225°].
        // Neither contains the intersection angles.
        assert!(hits.is_empty(), "hits={hits:?}");
    }
}
