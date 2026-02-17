use std::f64::consts::PI;

/// Returns the minimum distance from point `(px, py)` to the line segment
/// from `(ax, ay)` to `(bx, by)`.
#[must_use]
pub fn point_to_segment_dist(px: f64, py: f64, ax: f64, ay: f64, bx: f64, by: f64) -> f64 {
    let dx = bx - ax;
    let dy = by - ay;
    let len_sq = dx * dx + dy * dy;

    if len_sq < 1e-20 {
        // Degenerate segment (zero length).
        return ((px - ax).powi(2) + (py - ay).powi(2)).sqrt();
    }

    // Project point onto the infinite line, clamp to [0, 1].
    let t = ((px - ax) * dx + (py - ay) * dy) / len_sq;
    let t = t.clamp(0.0, 1.0);

    let closest_x = ax + t * dx;
    let closest_y = ay + t * dy;

    ((px - closest_x).powi(2) + (py - closest_y).powi(2)).sqrt()
}

/// Returns the minimum distance from point `(px, py)` to a circular arc.
///
/// The arc is defined by center `(cx, cy)`, `radius`, `start_angle`, and `sweep`.
///
/// If the point's angle (relative to center) falls within the arc range,
/// the distance is `||point - center| - radius|`.
/// Otherwise, the distance is the minimum of the distances to the two arc endpoints.
#[must_use]
pub fn point_to_arc_dist(
    px: f64,
    py: f64,
    cx: f64,
    cy: f64,
    radius: f64,
    start_angle: f64,
    sweep: f64,
) -> f64 {
    let dx = px - cx;
    let dy = py - cy;
    let dist_to_center = (dx * dx + dy * dy).sqrt();

    // Check if the point's angle is within the arc range.
    let angle = dy.atan2(dx);
    if angle_in_arc_range(angle, start_angle, sweep) {
        return (dist_to_center - radius).abs();
    }

    // Point is outside the arc's angular range. Check distance to endpoints.
    let end_angle = start_angle + sweep;
    let ep0_x = cx + radius * start_angle.cos();
    let ep0_y = cy + radius * start_angle.sin();
    let ep1_x = cx + radius * end_angle.cos();
    let ep1_y = cy + radius * end_angle.sin();

    let d0 = ((px - ep0_x).powi(2) + (py - ep0_y).powi(2)).sqrt();
    let d1 = ((px - ep1_x).powi(2) + (py - ep1_y).powi(2)).sqrt();

    d0.min(d1)
}

/// Checks if an angle falls within an arc's angular range.
fn angle_in_arc_range(angle: f64, start_angle: f64, sweep: f64) -> bool {
    let eps = 1e-10;
    let mut delta = angle - start_angle;

    if sweep > 0.0 {
        while delta < -eps {
            delta += 2.0 * PI;
        }
        while delta > 2.0 * PI + eps {
            delta -= 2.0 * PI;
        }
        delta >= -eps && delta <= sweep + eps
    } else {
        while delta > eps {
            delta -= 2.0 * PI;
        }
        while delta < -2.0 * PI - eps {
            delta += 2.0 * PI;
        }
        delta <= eps && delta >= sweep - eps
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    // ── point_to_segment_dist tests ──

    #[test]
    fn segment_dist_perpendicular_projection() {
        // Point (1, 1) to segment (0,0)→(2,0). Closest at (1,0), dist = 1.
        let d = point_to_segment_dist(1.0, 1.0, 0.0, 0.0, 2.0, 0.0);
        assert!((d - 1.0).abs() < TOL, "d={d}");
    }

    #[test]
    fn segment_dist_endpoint_closest() {
        // Point (-1, 0) to segment (0,0)→(2,0). Closest at (0,0), dist = 1.
        let d = point_to_segment_dist(-1.0, 0.0, 0.0, 0.0, 2.0, 0.0);
        assert!((d - 1.0).abs() < TOL, "d={d}");
    }

    #[test]
    fn segment_dist_on_segment() {
        // Point on the segment itself.
        let d = point_to_segment_dist(1.0, 0.0, 0.0, 0.0, 2.0, 0.0);
        assert!(d.abs() < TOL, "d={d}");
    }

    #[test]
    fn segment_dist_degenerate() {
        // Zero-length segment: distance is point-to-point.
        let d = point_to_segment_dist(3.0, 4.0, 0.0, 0.0, 0.0, 0.0);
        assert!((d - 5.0).abs() < TOL, "d={d}");
    }

    // ── point_to_arc_dist tests ──

    #[test]
    fn arc_dist_in_range() {
        // Point at (0, 2) to CCW semicircle centered at origin, radius 1.
        // Angle of point = π/2, which is in [0, π]. Distance = |2 - 1| = 1.
        let d = point_to_arc_dist(0.0, 2.0, 0.0, 0.0, 1.0, 0.0, PI);
        assert!((d - 1.0).abs() < TOL, "d={d}");
    }

    #[test]
    fn arc_dist_outside_range() {
        // Point at (0, -2) to CCW semicircle from 0 to π (upper half).
        // Angle = -π/2, not in [0, π].
        // Arc endpoints: (1, 0) and (-1, 0).
        // Distance to (1,0) = √(1+4) = √5, to (-1,0) = √(1+4) = √5.
        let d = point_to_arc_dist(0.0, -2.0, 0.0, 0.0, 1.0, 0.0, PI);
        let expected = 5.0_f64.sqrt();
        assert!((d - expected).abs() < 1e-6, "d={d}");
    }

    #[test]
    fn arc_dist_on_arc() {
        // Point exactly on the arc.
        let d = point_to_arc_dist(0.0, 1.0, 0.0, 0.0, 1.0, 0.0, PI);
        assert!(d.abs() < TOL, "d={d}");
    }

    #[test]
    fn arc_dist_inside_arc() {
        // Point at center. In angular range, distance = radius.
        let d = point_to_arc_dist(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, PI);
        assert!((d - 1.0).abs() < TOL, "d={d}");
    }
}
