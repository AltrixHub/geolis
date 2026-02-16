use crate::error::{OperationError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

/// When `cos(angle between consecutive segments) < this`, use a flat cap
/// instead of a miter join. Only for near-180° reversals (> ~169°).
const FLAT_CAP_COS: f64 = -0.98;

/// Maximum miter distance as a multiple of `|distance|`. When the miter
/// extends further than this, a bevel (two points) is used instead.
/// A limit of 4.0 clips at ~30° turn angles (matches SVG default).
const MITER_LIMIT: f64 = 4.0;

/// Offsets a 2D polyline by a given distance with self-intersection trimming.
///
/// Works directly with `Vec<Point3>` (XY plane, Z ignored) for immediate
/// usability without requiring the topology store.
///
/// # Algorithm
///
/// 1. **Phase A**: Offset each segment perpendicular to its direction
/// 2. **Phase B**: Build raw offset polyline by intersecting consecutive offset lines.
///    Near-antiparallel segments (hairpins, cross arms) get a flat cap instead of
///    a divergent miter.
/// 3. **Phase C**: Walk the raw polyline, detect self-intersections, remove loops
///
/// # Sign Convention
///
/// - Positive distance: left offset (relative to walking direction)
/// - Negative distance: right offset
#[derive(Debug)]
pub struct PolylineOffset2D {
    points: Vec<Point3>,
    distance: f64,
    closed: bool,
}

impl PolylineOffset2D {
    /// Creates a new polyline offset operation.
    #[must_use]
    pub fn new(points: Vec<Point3>, distance: f64, closed: bool) -> Self {
        Self {
            points,
            distance,
            closed,
        }
    }

    /// Executes the offset operation.
    ///
    /// For **closed** polylines, produces an inset/outset polygon.
    /// For **open** polylines, produces a closed polygon tracing both sides of the
    /// path with flat caps at the endpoints (like a stroke/buffer outline).
    /// Both `d > 0` and `d < 0` produce the same outline (using `|d|`).
    ///
    /// # Errors
    ///
    /// - `OperationError::InvalidInput` if fewer than 2 points are provided
    /// - `OperationError::Failed` if the offset collapses to fewer than 2 points
    pub fn execute(&self) -> Result<Vec<Point3>> {
        let n = self.points.len();
        if n < 2 {
            return Err(OperationError::InvalidInput(
                "at least 2 points are required for polyline offset".to_owned(),
            )
            .into());
        }

        if self.distance.abs() < TOLERANCE {
            return Ok(self.points.clone());
        }

        if self.closed {
            self.execute_closed()
        } else {
            self.execute_open()
        }
    }

    /// Executes offset for a closed polyline (inset/outset polygon).
    fn execute_closed(&self) -> Result<Vec<Point3>> {
        // Phase A & B: Build raw offset polyline
        let (raw, _flat_cap_segments) = self.build_raw_offset()?;

        if raw.len() < 2 {
            return Err(OperationError::Failed(
                "offset collapsed to fewer than 2 points".to_owned(),
            )
            .into());
        }

        // Phase C: Remove self-intersections
        let winding_sign = signed_area_2d(&self.points).signum();
        let trimmed = trim_closed_loops(&raw, winding_sign);

        if trimmed.len() < 2 {
            return Err(OperationError::Failed(
                "offset collapsed to fewer than 2 points".to_owned(),
            )
            .into());
        }

        // Detect collapse: inward offset should reduce area.
        if trimmed.len() >= 3 {
            let original_area = signed_area_2d(&self.points);
            let result_area = signed_area_2d(&trimmed);
            let is_inward = original_area * self.distance > 0.0;
            if is_inward
                && original_area.abs() > TOLERANCE
                && result_area.abs() > original_area.abs()
            {
                return Err(OperationError::Failed(
                    "offset collapsed (passed through center)".to_owned(),
                )
                .into());
            }
        }

        Ok(trimmed)
    }

    /// Executes offset for an open polyline: produces a closed polygon tracing
    /// both sides with flat caps (buffer/stroke outline).
    fn execute_open(&self) -> Result<Vec<Point3>> {
        let d = self.distance.abs();

        // Forward offset: left side of original path at +|d|.
        let forward = build_one_side_offset(&self.points, d)?;

        // Backward offset: reverse input, left offset at +|d|.
        // For polylines with 180° reversals, the backward path retraces arm
        // segments already covered by the forward path's flat caps. Filter
        // these duplicates to avoid collinear overlapping edges that confuse
        // the self-intersection trimmer.
        let reversed_input: Vec<Point3> = self.points.iter().copied().rev().collect();
        let backward = build_one_side_offset(&reversed_input, d)?;

        let dup_tol_sq = TOLERANCE * TOLERANCE * 1e8; // ~1e-12
        let filtered_backward: Vec<Point3> = backward
            .iter()
            .filter(|bp| {
                !forward
                    .iter()
                    .any(|fp| (bp.x - fp.x).powi(2) + (bp.y - fp.y).powi(2) < dup_tol_sq)
            })
            .copied()
            .collect();

        // Stitch into a closed polygon: forward + filtered backward.
        let mut combined = forward;
        combined.extend(filtered_backward);

        if combined.len() < 3 {
            return Err(OperationError::Failed(
                "offset collapsed to fewer than 3 points".to_owned(),
            )
            .into());
        }

        // Trim self-intersections using the combined polygon's own winding.
        let winding = signed_area_2d(&combined).signum();
        let sign = if winding.abs() < 0.5 { 1.0 } else { winding };
        let trimmed = trim_closed_loops(&combined, sign);

        if trimmed.len() < 3 {
            return Err(OperationError::Failed(
                "offset collapsed to fewer than 3 points".to_owned(),
            )
            .into());
        }

        // Canonicalize: rotate to leftmost-bottommost vertex for deterministic output.
        Ok(rotate_to_canonical_start(&trimmed))
    }

    /// Phase A & B: Offset each segment and intersect consecutive offset lines
    /// to build the raw (untrimmed) offset polyline.
    ///
    /// Near-antiparallel segments (angle > ~154°) produce a **flat cap** (two
    /// points) instead of a single miter point that would diverge to infinity.
    ///
    /// Returns the raw offset points and a sorted list of flat-cap segment
    /// indices (segments connecting the two flat-cap points at each reversal).
    fn build_raw_offset(&self) -> Result<(Vec<Point3>, Vec<usize>)> {
        if self.closed {
            build_one_side_offset_closed(&self.points, self.distance)
        } else {
            let raw = build_one_side_offset(&self.points, self.distance)?;
            Ok((raw, Vec::new()))
        }
    }

}

/// Builds a one-sided offset for an open polyline (Phase A + B).
///
/// Takes raw points and a signed distance. Returns the offset path as a list
/// of points (no flat-cap tracking, used for both-sides stitching).
fn build_one_side_offset(points: &[Point3], distance: f64) -> Result<Vec<Point3>> {
    let n = points.len();
    let segment_count = n - 1;

    // Phase A: Compute offset segments.
    let mut offset_segments: Vec<(Point3, Point3)> = Vec::with_capacity(segment_count);
    let mut directions: Vec<Vector3> = Vec::with_capacity(segment_count);

    for i in 0..segment_count {
        let j = i + 1;
        let dir = segment_direction(&points[i], &points[j])?;
        let normal = left_normal(dir);
        let offset = normal * distance;

        let a = Point3::new(
            points[i].x + offset.x,
            points[i].y + offset.y,
            points[i].z,
        );
        let b = Point3::new(
            points[j].x + offset.x,
            points[j].y + offset.y,
            points[j].z,
        );
        offset_segments.push((a, b));
        directions.push(dir);
    }

    if offset_segments.is_empty() {
        return Err(OperationError::Failed("no valid segments to offset".to_owned()).into());
    }

    // Phase B: Build raw polyline by intersecting consecutive offset segments.
    let mut raw = Vec::with_capacity(n * 2);

    // First point: start of the first offset segment.
    raw.push(offset_segments[0].0);

    // Interior corners.
    for i in 1..n - 1 {
        push_corner(
            &mut raw,
            &offset_segments[i - 1],
            &offset_segments[i],
            &directions[i - 1],
            &directions[i],
            &points[i],
            distance,
        );
    }

    // Last point: end of the last offset segment.
    raw.push(offset_segments[segment_count - 1].1);

    Ok(raw)
}

/// Builds a one-sided offset for a closed polyline (Phase A + B) with flat-cap tracking.
fn build_one_side_offset_closed(
    points: &[Point3],
    distance: f64,
) -> Result<(Vec<Point3>, Vec<usize>)> {
    let n = points.len();
    let segment_count = n;

    // Phase A: Compute offset segments.
    let mut offset_segments: Vec<(Point3, Point3)> = Vec::with_capacity(segment_count);
    let mut directions: Vec<Vector3> = Vec::with_capacity(segment_count);

    for i in 0..segment_count {
        let j = (i + 1) % n;
        let dir = segment_direction(&points[i], &points[j])?;
        let normal = left_normal(dir);
        let offset = normal * distance;

        let a = Point3::new(
            points[i].x + offset.x,
            points[i].y + offset.y,
            points[i].z,
        );
        let b = Point3::new(
            points[j].x + offset.x,
            points[j].y + offset.y,
            points[j].z,
        );
        offset_segments.push((a, b));
        directions.push(dir);
    }

    if offset_segments.is_empty() {
        return Err(OperationError::Failed("no valid segments to offset".to_owned()).into());
    }

    // Phase B: Build raw polyline by intersecting consecutive offset segments.
    let mut raw = Vec::with_capacity(n * 2);
    let mut flat_cap_segments: Vec<usize> = Vec::new();

    for i in 0..segment_count {
        let prev = if i == 0 { segment_count - 1 } else { i - 1 };
        let pre_len = raw.len();
        let is_flat_cap = push_corner(
            &mut raw,
            &offset_segments[prev],
            &offset_segments[i],
            &directions[prev],
            &directions[i],
            &points[i],
            distance,
        );
        if is_flat_cap {
            flat_cap_segments.push(pre_len);
        }
    }

    Ok((raw, flat_cap_segments))
}

/// Checks whether segments i and j are adjacent in a closed polyline.
fn are_adjacent(i: usize, j: usize, n: usize) -> bool {
    let diff = i.abs_diff(j);
    diff == 1 || diff == n - 1
}

/// Finds the first self-intersection between non-adjacent segments in a closed polygon.
///
/// Skips endpoint-to-endpoint touches (both `t` and `u` near 0 or 1) which
/// occur when the combined outline revisits the same geometric point at
/// non-adjacent vertex positions. Only genuine crossings (at least one
/// parameter in the interior) are reported.
///
/// Returns `(i, j, intersection_point)` where `i < j` are segment indices.
fn find_first_self_intersection(points: &[Point3]) -> Option<(usize, usize, Point3)> {
    let n = points.len();
    if n < 4 {
        return None;
    }
    let eps = TOLERANCE * 100.0;
    for i in 0..n {
        let i_next = (i + 1) % n;
        for j in (i + 2)..n {
            if are_adjacent(i, j, n) {
                continue;
            }
            let j_next = (j + 1) % n;
            if let Some((pt, t, u)) = segment_segment_intersect_2d(
                &points[i],
                &points[i_next],
                &points[j],
                &points[j_next],
            ) {
                // Skip endpoint-to-endpoint touches: both parameters at endpoints.
                let t_at_end = t < eps || t > 1.0 - eps;
                let u_at_end = u < eps || u > 1.0 - eps;
                if t_at_end && u_at_end {
                    continue;
                }
                return Some((i, j, pt));
            }
        }
    }
    None
}

/// Splits a closed polygon at the intersection of segments `i` and `j` into two sub-polygons.
///
/// Assumes `i < j`. Returns two vertex lists representing the two loops created by the split:
/// - Sub-path A: `[intersection, P(i+1), ..., P(j)]`
/// - Sub-path B: `[intersection, P(j+1), ..., P(i)]` (wrapping around)
#[allow(clippy::many_single_char_names)]
fn split_at_intersection(
    points: &[Point3],
    seg_i: usize,
    seg_j: usize,
    intersection: Point3,
) -> (Vec<Point3>, Vec<Point3>) {
    let n = points.len();

    // Sub-path A: intersection, then vertices (seg_i+1) through seg_j inclusive.
    let mut a = Vec::with_capacity(seg_j - seg_i + 1);
    a.push(intersection);
    a.extend_from_slice(&points[(seg_i + 1)..=seg_j]);

    // Sub-path B: intersection, then vertices (seg_j+1)%n through seg_i inclusive (wrapping).
    let b_vertex_count = n - (seg_j - seg_i);
    let mut b = Vec::with_capacity(b_vertex_count + 1);
    b.push(intersection);
    let mut idx = (seg_j + 1) % n;
    loop {
        b.push(points[idx]);
        if idx == seg_i {
            break;
        }
        idx = (idx + 1) % n;
    }

    (a, b)
}

/// Removes degenerate vertices from a closed polygon: consecutive duplicates
/// and collinear (on-edge) points.
///
/// Essential after splitting at intersections that land on existing vertices
/// (e.g., cross-shape arm collapse produces endpoint intersections that create
/// duplicate and collinear vertices in sub-paths).
fn clean_polygon(points: &[Point3]) -> Vec<Point3> {
    if points.len() < 3 {
        return points.to_vec();
    }

    // Step 1: Remove consecutive near-duplicates.
    let tol_sq = TOLERANCE * TOLERANCE * 100.0; // 1e-18 — well within float noise
    let mut deduped: Vec<Point3> = Vec::with_capacity(points.len());
    for &pt in points {
        if let Some(&last) = deduped.last() {
            if (pt.x - last.x).powi(2) + (pt.y - last.y).powi(2) < tol_sq {
                continue;
            }
        }
        deduped.push(pt);
    }
    // Wrap-around: check last vs first.
    if deduped.len() > 1 {
        let first = deduped[0];
        let last = deduped[deduped.len() - 1];
        if (last.x - first.x).powi(2) + (last.y - first.y).powi(2) < tol_sq {
            deduped.pop();
        }
    }

    if deduped.len() < 3 {
        return deduped;
    }

    // Step 2: Remove collinear vertices (single pass, cross-product check).
    let n = deduped.len();
    let mut cleaned = Vec::with_capacity(n);
    for i in 0..n {
        let prev = if i == 0 { n - 1 } else { i - 1 };
        let next = (i + 1) % n;
        let cross = (deduped[i].x - deduped[prev].x) * (deduped[next].y - deduped[i].y)
            - (deduped[i].y - deduped[prev].y) * (deduped[next].x - deduped[i].x);
        if cross.abs() >= TOLERANCE {
            cleaned.push(deduped[i]);
        }
    }

    // Don't reduce below 3 vertices — return deduped as fallback.
    if cleaned.len() < 3 {
        return deduped;
    }
    cleaned
}

/// Recursively removes self-intersection loops from a closed polygon.
///
/// At each self-intersection, splits into two sub-polygons, recursively trims
/// both, then keeps the one whose winding matches the original polygon.
/// Junk loops from collapsed features wind oppositely and are discarded.
///
/// `winding_sign` is `+1.0` for CCW originals, `-1.0` for CW.
/// Convergence is guaranteed because each split strictly reduces vertex count.
fn trim_closed_loops(points: &[Point3], winding_sign: f64) -> Vec<Point3> {
    // Clean up degenerate vertices before intersection detection.
    // Arm-collapse offsets can produce duplicate/collinear vertices from
    // endpoint intersections that confuse the self-intersection scanner.
    let pts = clean_polygon(points);

    if pts.len() < 4 {
        return pts;
    }
    match find_first_self_intersection(&pts) {
        None => pts,
        Some((i, j, pt)) => {
            let (a, b) = split_at_intersection(&pts, i, j, pt);
            let trimmed_a = trim_closed_loops(&a, winding_sign);
            let trimmed_b = trim_closed_loops(&b, winding_sign);
            let area_a = signed_area_2d(&trimmed_a);
            let area_b = signed_area_2d(&trimmed_b);
            let a_correct = area_a * winding_sign > 0.0;
            let b_correct = area_b * winding_sign > 0.0;
            match (a_correct, b_correct) {
                (true, false) => trimmed_a,
                (false, true) => trimmed_b,
                _ => {
                    // Both match or neither: keep larger absolute area.
                    if area_a.abs() >= area_b.abs() {
                        trimmed_a
                    } else {
                        trimmed_b
                    }
                }
            }
        }
    }
}

/// Pushes corner point(s) into `raw`.
///
/// Returns `true` if a flat cap (180° reversal barrier) was inserted.
///
/// - Near-antiparallel segments: flat cap (two points) with barrier.
/// - Miter exceeding `MITER_LIMIT`: bevel (two points), no barrier.
/// - Normal corners: single miter intersection point.
fn push_corner(
    raw: &mut Vec<Point3>,
    seg_prev: &(Point3, Point3),
    seg_next: &(Point3, Point3),
    dir_prev: &Vector3,
    dir_next: &Vector3,
    original_corner: &Point3,
    distance: f64,
) -> bool {
    let cos_angle = dir_prev.x * dir_next.x + dir_prev.y * dir_next.y;

    if cos_angle < FLAT_CAP_COS {
        // Near-antiparallel (> ~169°): flat cap with barrier.
        raw.push(seg_prev.1);
        raw.push(seg_next.0);
        return true;
    }

    // Normal corner: try miter, fall back to bevel if too long.
    let corner = intersect_offset_lines(seg_prev, seg_next, original_corner, distance);
    let dx = corner.x - original_corner.x;
    let dy = corner.y - original_corner.y;
    let miter_dist_sq = dx * dx + dy * dy;
    let limit = MITER_LIMIT * distance.abs();

    if miter_dist_sq > limit * limit {
        // Miter too long: bevel (two points), no barrier.
        raw.push(seg_prev.1);
        raw.push(seg_next.0);
    } else {
        raw.push(corner);
    }
    false
}

/// Computes the normalized direction from point `a` to point `b`.
fn segment_direction(a: &Point3, b: &Point3) -> Result<Vector3> {
    let d = b - a;
    let len = (d.x * d.x + d.y * d.y).sqrt();
    if len < TOLERANCE {
        return Err(OperationError::InvalidInput(format!(
            "zero-length segment between ({}, {}) and ({}, {})",
            a.x, a.y, b.x, b.y
        ))
        .into());
    }
    Ok(Vector3::new(d.x / len, d.y / len, 0.0))
}

/// Returns the left-pointing normal of a direction vector in the XY plane.
fn left_normal(dir: Vector3) -> Vector3 {
    Vector3::new(-dir.y, dir.x, 0.0)
}

/// Intersects two offset lines and returns the corner point.
///
/// Falls back to the offset midpoint if lines are parallel.
fn intersect_offset_lines(
    seg_prev: &(Point3, Point3),
    seg_next: &(Point3, Point3),
    original_corner: &Point3,
    distance: f64,
) -> Point3 {
    let d_prev = Vector3::new(
        seg_prev.1.x - seg_prev.0.x,
        seg_prev.1.y - seg_prev.0.y,
        0.0,
    );
    let d_next = Vector3::new(
        seg_next.1.x - seg_next.0.x,
        seg_next.1.y - seg_next.0.y,
        0.0,
    );

    if let Some((t, _u)) = line_line_intersect_2d(&seg_prev.1, &d_prev, &seg_next.0, &d_next) {
        // Intersection point from the end of the previous segment.
        point_at(&seg_prev.1, &d_prev, t)
    } else {
        // Parallel segments: use the endpoint of the previous offset segment,
        // or fall back to shifting the original corner.
        let normal = left_normal(
            Vector3::new(d_prev.x, d_prev.y, 0.0)
                .try_normalize(TOLERANCE)
                .unwrap_or(Vector3::new(1.0, 0.0, 0.0)),
        );
        Point3::new(
            original_corner.x + normal.x * distance,
            original_corner.y + normal.y * distance,
            original_corner.z,
        )
    }
}

/// Parametric 2D line-line intersection.
///
/// Given lines `p1 + t * d1` and `p2 + u * d2`, returns `(t, u)` if not parallel.
fn line_line_intersect_2d(
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
fn segment_segment_intersect_2d(
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
fn point_at(origin: &Point3, dir: &Vector3, t: f64) -> Point3 {
    Point3::new(origin.x + dir.x * t, origin.y + dir.y * t, origin.z)
}

/// Rotates a closed polygon so it starts at the leftmost vertex (smallest x),
/// breaking ties by smallest y. Ensures deterministic output for tests.
fn rotate_to_canonical_start(points: &[Point3]) -> Vec<Point3> {
    if points.len() < 2 {
        return points.to_vec();
    }
    let mut best = 0;
    for (i, pt) in points.iter().enumerate().skip(1) {
        let b = &points[best];
        if pt.x < b.x - TOLERANCE || (pt.x - b.x).abs() < TOLERANCE && pt.y < b.y {
            best = i;
        }
    }
    if best == 0 {
        return points.to_vec();
    }
    let mut rotated = Vec::with_capacity(points.len());
    rotated.extend_from_slice(&points[best..]);
    rotated.extend_from_slice(&points[..best]);
    rotated
}

/// Computes the signed area of a polygon in the XY plane (shoelace formula).
///
/// Positive for counter-clockwise, negative for clockwise.
fn signed_area_2d(points: &[Point3]) -> f64 {
    let n = points.len();
    if n < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        sum += points[i].x * points[j].y - points[j].x * points[i].y;
    }
    sum * 0.5
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// Helper: asserts two points are approximately equal.
    fn assert_point_near(a: &Point3, b: &Point3, tol: f64, msg: &str) {
        let d = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt();
        assert!(
            d < tol,
            "{msg}: expected ({}, {}), got ({}, {}), dist={d}",
            b.x,
            b.y,
            a.x,
            a.y
        );
    }

    #[test]
    fn straight_line_offset() {
        // Horizontal line: both-sides offset produces a rectangle.
        let points = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0)];
        let op = PolylineOffset2D::new(points, 1.0, false);
        let result = op.execute().unwrap();

        assert_eq!(result.len(), 4);
        // Canonical start: leftmost-bottommost = (0, -1).
        assert_point_near(&result[0], &Point3::new(0.0, -1.0, 0.0), 1e-9, "BL");
        assert_point_near(&result[1], &Point3::new(0.0, 1.0, 0.0), 1e-9, "TL");
        assert_point_near(&result[2], &Point3::new(10.0, 1.0, 0.0), 1e-9, "TR");
        assert_point_near(&result[3], &Point3::new(10.0, -1.0, 0.0), 1e-9, "BR");
    }

    #[test]
    fn l_shape_90_degree() {
        // Right-going then up-going: both-sides offset → 6-vertex polygon.
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
        ];
        let op = PolylineOffset2D::new(points, 1.0, false);
        let result = op.execute().unwrap();

        assert_eq!(result.len(), 6);
        // Canonical start: leftmost-bottommost = (0, -1).
        assert_point_near(&result[0], &Point3::new(0.0, -1.0, 0.0), 1e-9, "v0");
        assert_point_near(&result[1], &Point3::new(0.0, 1.0, 0.0), 1e-9, "v1");
        assert_point_near(&result[2], &Point3::new(9.0, 1.0, 0.0), 1e-9, "v2");
        assert_point_near(&result[3], &Point3::new(9.0, 10.0, 0.0), 1e-9, "v3");
        assert_point_near(&result[4], &Point3::new(11.0, 10.0, 0.0), 1e-9, "v4");
        assert_point_near(&result[5], &Point3::new(11.0, -1.0, 0.0), 1e-9, "v5");
    }

    #[test]
    fn reversal_180_degree_flat_cap() {
        // Exact 180-degree reversal: go right then back left.
        // Should produce a flat cap (2 points) at the tip instead of diverging.
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
        ];
        let op = PolylineOffset2D::new(points, 1.0, false);
        let result = op.execute().unwrap();

        // Should produce a valid result (not divergent spike).
        assert!(result.len() >= 2, "should have at least 2 points");
        // All points should be within reasonable bounds.
        for pt in &result {
            assert!(
                pt.x.abs() < 15.0 && pt.y.abs() < 5.0,
                "point ({}, {}) out of bounds — flat cap failed",
                pt.x,
                pt.y
            );
        }
    }

    #[test]
    fn cross_pattern_no_divergence() {
        // Cross / plus shape: center with 4 arms, each arm reverses.
        let c = Point3::new(0.0, 0.0, 0.0);
        let points = vec![
            Point3::new(-5.0, 0.0, 0.0),
            c,
            Point3::new(0.0, 5.0, 0.0),
            c,
            Point3::new(5.0, 0.0, 0.0),
            c,
            Point3::new(0.0, -5.0, 0.0),
        ];
        let op = PolylineOffset2D::new(points, 0.5, false);
        let result = op.execute().unwrap();

        assert!(result.len() >= 2);
        // No point should diverge beyond reasonable bounds.
        for pt in &result {
            assert!(
                pt.x.abs() < 10.0 && pt.y.abs() < 10.0,
                "point ({}, {}) diverged",
                pt.x,
                pt.y
            );
        }
    }

    #[test]
    fn closed_square_inward_offset() {
        // CCW square: left normal points inward, so positive distance = inward.
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
            Point3::new(0.0, 10.0, 0.0),
        ];
        let op = PolylineOffset2D::new(points, 1.0, true);
        let result = op.execute().unwrap();

        assert_eq!(result.len(), 4);
        assert_point_near(&result[0], &Point3::new(1.0, 1.0, 0.0), 1e-9, "corner 0");
        assert_point_near(&result[1], &Point3::new(9.0, 1.0, 0.0), 1e-9, "corner 1");
        assert_point_near(&result[2], &Point3::new(9.0, 9.0, 0.0), 1e-9, "corner 2");
        assert_point_near(&result[3], &Point3::new(1.0, 9.0, 0.0), 1e-9, "corner 3");
    }

    #[test]
    fn closed_triangle_collapse_error() {
        // Small equilateral triangle with large inward offset -> collapse.
        // CCW triangle: positive distance = left = inward.
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(1.0, 1.732, 0.0),
        ];
        // Inradius ≈ 0.577, so offset by 5.0 inward should collapse.
        let op = PolylineOffset2D::new(points, 5.0, true);
        assert!(op.execute().is_err());
    }

    #[test]
    fn parallel_consecutive_segments() {
        // Collinear points: both-sides offset produces a rectangle
        // (collinear intermediate vertices cleaned by `clean_polygon`).
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
        ];
        let op = PolylineOffset2D::new(points, 1.0, false);
        let result = op.execute().unwrap();

        assert_eq!(result.len(), 4);
        assert_point_near(&result[0], &Point3::new(0.0, -1.0, 0.0), 1e-9, "BL");
        assert_point_near(&result[1], &Point3::new(0.0, 1.0, 0.0), 1e-9, "TL");
        assert_point_near(&result[2], &Point3::new(10.0, 1.0, 0.0), 1e-9, "TR");
        assert_point_near(&result[3], &Point3::new(10.0, -1.0, 0.0), 1e-9, "BR");
    }

    #[test]
    fn zero_distance_returns_copy() {
        let points = vec![Point3::new(1.0, 2.0, 0.0), Point3::new(3.0, 4.0, 0.0)];
        let op = PolylineOffset2D::new(points.clone(), 0.0, false);
        let result = op.execute().unwrap();

        assert_eq!(result.len(), 2);
        assert_point_near(&result[0], &points[0], 1e-12, "copy start");
        assert_point_near(&result[1], &points[1], 1e-12, "copy end");
    }

    #[test]
    fn fewer_than_two_points_error() {
        let points = vec![Point3::new(0.0, 0.0, 0.0)];
        let op = PolylineOffset2D::new(points, 1.0, false);
        assert!(op.execute().is_err());

        let empty: Vec<Point3> = vec![];
        let op = PolylineOffset2D::new(empty, 1.0, false);
        assert!(op.execute().is_err());
    }

    #[test]
    fn sharp_hairpin_trimming() {
        // Hairpin: goes right then back left, nearly overlapping.
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.5, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        // Offset by 2.0 (> half the 0.5 height separation).
        let op = PolylineOffset2D::new(points, 2.0, false);
        let result = op.execute().unwrap();

        // Should be trimmed due to self-intersection.
        assert!(result.len() >= 2, "should have at least 2 points");
    }

    #[test]
    fn negative_distance_right_offset() {
        // Horizontal line: d<0 produces the same both-sides outline as d>0 (uses |d|).
        let points = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0)];
        let op = PolylineOffset2D::new(points, -1.0, false);
        let result = op.execute().unwrap();

        assert_eq!(result.len(), 4);
        assert_point_near(&result[0], &Point3::new(0.0, -1.0, 0.0), 1e-9, "BL");
        assert_point_near(&result[1], &Point3::new(0.0, 1.0, 0.0), 1e-9, "TL");
        assert_point_near(&result[2], &Point3::new(10.0, 1.0, 0.0), 1e-9, "TR");
        assert_point_near(&result[3], &Point3::new(10.0, -1.0, 0.0), 1e-9, "BR");
    }

    #[test]
    fn closed_square_outward_offset() {
        // CCW square: negative distance (right) = outward.
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
            Point3::new(0.0, 10.0, 0.0),
        ];
        let op = PolylineOffset2D::new(points, -1.0, true);
        let result = op.execute().unwrap();

        assert_eq!(result.len(), 4);
        assert_point_near(&result[0], &Point3::new(-1.0, -1.0, 0.0), 1e-9, "corner 0");
        assert_point_near(&result[1], &Point3::new(11.0, -1.0, 0.0), 1e-9, "corner 1");
        assert_point_near(&result[2], &Point3::new(11.0, 11.0, 0.0), 1e-9, "corner 2");
        assert_point_near(&result[3], &Point3::new(-1.0, 11.0, 0.0), 1e-9, "corner 3");
    }

    // ── Failing pattern tests (from viewer) ────────────────────────

    /// Helper: assert no point diverges beyond `max_dist` from origin.
    fn assert_bounded(result: &[Point3], origin: &Point3, max_dist: f64, label: &str) {
        for (i, pt) in result.iter().enumerate() {
            let d = ((pt.x - origin.x).powi(2) + (pt.y - origin.y).powi(2)).sqrt();
            assert!(
                d < max_dist,
                "{label} point {i} ({}, {}) is {d:.2} from origin, exceeds {max_dist}",
                pt.x, pt.y
            );
        }
    }

    /// Helper: check that a closed offset polygon has no self-intersections.
    fn assert_no_self_intersection(result: &[Point3], label: &str) {
        let n = result.len();
        for i in 0..n {
            let i_next = (i + 1) % n;
            for j in (i + 2)..n {
                if i == 0 && j == n - 1 {
                    continue; // adjacent (first-last)
                }
                let j_next = (j + 1) % n;
                if let Some((_pt, _t, _u)) = segment_segment_intersect_2d(
                    &result[i],
                    &result[i_next],
                    &result[j],
                    &result[j_next],
                ) {
                    panic!(
                        "{label}: self-intersection between segments {i}-{i_next} and {j}-{j_next}"
                    );
                }
            }
        }
    }


    #[test]
    fn closed_arrow_outward_no_spike() {
        // Closed arrow / chevron (concave) — outward offset.
        let points = vec![
            Point3::new(1.0, -5.0, 0.0),
            Point3::new(5.0, -3.5, 0.0),
            Point3::new(1.0, -2.5, 0.0),
            Point3::new(2.5, -3.5, 0.0),
        ];
        let op = PolylineOffset2D::new(points.clone(), -0.2, true);
        let result = op.execute().unwrap();
        // All points should stay within reasonable bounds of the original shape.
        assert_bounded(&result, &Point3::new(3.0, -3.75, 0.0), 5.0, "arrow outward");
        assert!(result.len() >= 4, "should have at least 4 points");
    }

    #[test]
    fn closed_arrow_inward_collapses() {
        // Closed arrow / chevron is too thin at its acute vertices for any
        // meaningful inward offset — the self-intersection trimmer correctly
        // reduces it. Should not panic or produce divergent points.
        let points = vec![
            Point3::new(1.0, -5.0, 0.0),
            Point3::new(5.0, -3.5, 0.0),
            Point3::new(1.0, -2.5, 0.0),
            Point3::new(2.5, -3.5, 0.0),
        ];
        let op = PolylineOffset2D::new(points, 0.1, true);
        let result = op.execute().unwrap();
        // Result is trimmed but should not diverge.
        assert_bounded(&result, &Point3::new(3.0, -3.75, 0.0), 5.0, "arrow inward");
    }

    #[test]
    fn closed_cross_outward_no_self_intersection() {
        // Closed cross / plus outline (12 vertices, 8 concave corners).
        let points = vec![
            Point3::new(2.0, 3.5, 0.0),
            Point3::new(3.5, 3.5, 0.0),
            Point3::new(3.5, 2.0, 0.0),
            Point3::new(4.5, 2.0, 0.0),
            Point3::new(4.5, 3.5, 0.0),
            Point3::new(6.0, 3.5, 0.0),
            Point3::new(6.0, 4.5, 0.0),
            Point3::new(4.5, 4.5, 0.0),
            Point3::new(4.5, 6.0, 0.0),
            Point3::new(3.5, 6.0, 0.0),
            Point3::new(3.5, 4.5, 0.0),
            Point3::new(2.0, 4.5, 0.0),
        ];
        // Outward (negative for CCW).
        let op = PolylineOffset2D::new(points, -0.3, true);
        let result = op.execute().unwrap();
        assert_bounded(
            &result,
            &Point3::new(4.0, 4.0, 0.0),
            5.0,
            "cross outward",
        );
        assert_no_self_intersection(&result, "cross outward");
    }

    #[test]
    fn closed_t_shape_outward_no_self_intersection() {
        // Closed T-shape outline (8 vertices, concave).
        let points = vec![
            Point3::new(8.0, 2.0, 0.0),
            Point3::new(13.0, 2.0, 0.0),
            Point3::new(13.0, 3.0, 0.0),
            Point3::new(11.5, 3.0, 0.0),
            Point3::new(11.5, 6.0, 0.0),
            Point3::new(9.5, 6.0, 0.0),
            Point3::new(9.5, 3.0, 0.0),
            Point3::new(8.0, 3.0, 0.0),
        ];
        let op = PolylineOffset2D::new(points, -0.2, true);
        let result = op.execute().unwrap();
        assert_bounded(
            &result,
            &Point3::new(10.5, 4.0, 0.0),
            5.0,
            "T outward",
        );
        assert_no_self_intersection(&result, "T outward");
    }

    #[test]
    fn closed_star_no_spike() {
        // Closed star (5-pointed, sharp angles ~36°).
        let points = vec![
            Point3::new(9.5, -5.0, 0.0),
            Point3::new(10.1, -6.8, 0.0),
            Point3::new(12.0, -6.8, 0.0),
            Point3::new(10.5, -8.0, 0.0),
            Point3::new(11.1, -9.8, 0.0),
            Point3::new(9.5, -8.8, 0.0),
            Point3::new(7.9, -9.8, 0.0),
            Point3::new(8.5, -8.0, 0.0),
            Point3::new(7.0, -6.8, 0.0),
            Point3::new(8.9, -6.8, 0.0),
        ];
        let center = Point3::new(9.5, -7.5, 0.0);
        // Outward offset.
        let op = PolylineOffset2D::new(points.clone(), -0.15, true);
        let result = op.execute().unwrap();
        // Star tips should not spike far out.
        assert_bounded(&result, &center, 5.0, "star outward");
        // Inward offset.
        let op = PolylineOffset2D::new(points, 0.15, true);
        let result = op.execute().unwrap();
        assert_bounded(&result, &center, 5.0, "star inward");
    }

    #[test]
    fn open_cross_reversal_preserves_all_arms() {
        // Open cross: center with arms, 180° reversals.
        // Both +/- offsets should trace all arms without losing segments.
        let points = open_cross_points();
        for dist in &[0.3, -0.3] {
            let op = PolylineOffset2D::new(points.clone(), *dist, false);
            let result = op.execute().unwrap();
            // 6 original segments → 9 raw points (3 flat caps × 2 points each + start + end).
            // No trimming should remove any arm.
            assert!(
                result.len() >= 8,
                "dist={dist}: expected ≥8 points, got {}",
                result.len()
            );
            // No point should diverge.
            assert_bounded(
                &result,
                &Point3::new(0.0, 0.0, 0.0),
                5.0,
                &format!("cross reversal dist={dist}"),
            );
        }
    }

    #[test]
    fn open_t_junction_reversal_preserves_arms() {
        // T-junction: stem + crossbar with reversal at junction.
        // Both-sides offset produces a closed polygon. The self-intersection
        // trimmer keeps the largest valid sub-polygon (≥4 vertices).
        let points = vec![
            Point3::new(-4.0, 0.0, 0.0),
            Point3::new(-2.0, 0.0, 0.0),
            Point3::new(-2.0, -1.5, 0.0),
            Point3::new(-2.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
        ];
        for dist in &[0.25, -0.25] {
            let op = PolylineOffset2D::new(points.clone(), *dist, false);
            let result = op.execute().unwrap();
            assert!(
                result.len() >= 4,
                "dist={dist}: expected ≥4 points, got {}",
                result.len()
            );
            assert_bounded(
                &result,
                &Point3::new(-2.0, 0.0, 0.0),
                5.0,
                &format!("T-junction dist={dist}"),
            );
        }
    }

    // ── T-shape ground truth tests ─────────────────────────────────────
    //
    // T-shape: bar 10×1 (x:0..10, y:0..1), stem 4×5 (x:3..7, y:1..6).
    // CCW vertices: (0,0)→(10,0)→(10,1)→(7,1)→(7,6)→(3,6)→(3,1)→(0,1)
    //
    // For CCW polygon: positive distance = inward, negative = outward.
    //
    // Inward bar-collapse threshold: d > 0.5 (bar height=1, so 2d > 1).
    // When bar collapses, only the stem rectangle survives.
    //
    // Hand-computed expected vertices verified visually in the viewer
    // (examples/viewer/patterns/offset_intersection.rs).

    /// T-shape CCW vertices at origin.
    fn t_shape_points() -> Vec<Point3> {
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 1.0, 0.0),
            Point3::new(7.0, 1.0, 0.0),
            Point3::new(7.0, 6.0, 0.0),
            Point3::new(3.0, 6.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ]
    }

    /// Assert polygon vertices match expected (same count, same order).
    fn assert_polygon_eq(result: &[Point3], expected: &[Point3], tol: f64, label: &str) {
        assert_eq!(
            result.len(),
            expected.len(),
            "{label}: vertex count mismatch — got {}, expected {}",
            result.len(),
            expected.len()
        );
        for (i, (r, e)) in result.iter().zip(expected.iter()).enumerate() {
            assert_point_near(r, e, tol, &format!("{label} vertex {i}"));
        }
    }

    // ── Inward offset (positive distance, CCW) ────────────────────────

    #[test]
    fn t_shape_inward_d03_full_shape() {
        // d=0.3: bar and stem both survive → inset T-shape (8 vertices).
        let result = PolylineOffset2D::new(t_shape_points(), 0.3, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(0.3, 0.3, 0.0),
            Point3::new(9.7, 0.3, 0.0),
            Point3::new(9.7, 0.7, 0.0),
            Point3::new(6.7, 0.7, 0.0),
            Point3::new(6.7, 5.7, 0.0),
            Point3::new(3.3, 5.7, 0.0),
            Point3::new(3.3, 0.7, 0.0),
            Point3::new(0.3, 0.7, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "T inward d=0.3");
    }

    #[test]
    fn t_shape_inward_d06_bar_collapsed() {
        // d=0.6: bar collapses (2d=1.2 > bar height=1) → stem rectangle only.
        let result = PolylineOffset2D::new(t_shape_points(), 0.6, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(3.6, 0.6, 0.0),
            Point3::new(6.4, 0.6, 0.0),
            Point3::new(6.4, 5.4, 0.0),
            Point3::new(3.6, 5.4, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "T inward d=0.6");
    }

    #[test]
    fn t_shape_inward_d08_bar_collapsed() {
        // d=0.8: bar collapses → stem rectangle only.
        let result = PolylineOffset2D::new(t_shape_points(), 0.8, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(3.8, 0.8, 0.0),
            Point3::new(6.2, 0.8, 0.0),
            Point3::new(6.2, 5.2, 0.0),
            Point3::new(3.8, 5.2, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "T inward d=0.8");
    }

    #[test]
    fn t_shape_inward_d15_narrow_stem() {
        // d=1.5: bar collapsed, stem narrow (width=1.0, height=3.0).
        let result = PolylineOffset2D::new(t_shape_points(), 1.5, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(4.5, 1.5, 0.0),
            Point3::new(5.5, 1.5, 0.0),
            Point3::new(5.5, 4.5, 0.0),
            Point3::new(4.5, 4.5, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "T inward d=1.5");
    }

    // ── Outward offset (negative distance, CCW) ──────────────────────

    #[test]
    fn t_shape_outward_d03() {
        // d=0.3 outward: T-shape expands, 8 vertices, no self-intersection.
        let result = PolylineOffset2D::new(t_shape_points(), -0.3, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(-0.3, -0.3, 0.0),
            Point3::new(10.3, -0.3, 0.0),
            Point3::new(10.3, 1.3, 0.0),
            Point3::new(7.3, 1.3, 0.0),
            Point3::new(7.3, 6.3, 0.0),
            Point3::new(2.7, 6.3, 0.0),
            Point3::new(2.7, 1.3, 0.0),
            Point3::new(-0.3, 1.3, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "T outward d=0.3");
    }

    #[test]
    fn t_shape_outward_d06() {
        // d=0.6 outward: T-shape expands, 8 vertices.
        let result = PolylineOffset2D::new(t_shape_points(), -0.6, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(-0.6, -0.6, 0.0),
            Point3::new(10.6, -0.6, 0.0),
            Point3::new(10.6, 1.6, 0.0),
            Point3::new(7.6, 1.6, 0.0),
            Point3::new(7.6, 6.6, 0.0),
            Point3::new(2.4, 6.6, 0.0),
            Point3::new(2.4, 1.6, 0.0),
            Point3::new(-0.6, 1.6, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "T outward d=0.6");
    }

    #[test]
    fn t_shape_outward_d08() {
        // d=0.8 outward: T-shape expands, 8 vertices.
        let result = PolylineOffset2D::new(t_shape_points(), -0.8, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(-0.8, -0.8, 0.0),
            Point3::new(10.8, -0.8, 0.0),
            Point3::new(10.8, 1.8, 0.0),
            Point3::new(7.8, 1.8, 0.0),
            Point3::new(7.8, 6.8, 0.0),
            Point3::new(2.2, 6.8, 0.0),
            Point3::new(2.2, 1.8, 0.0),
            Point3::new(-0.8, 1.8, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "T outward d=0.8");
    }

    #[test]
    fn t_shape_outward_d15() {
        // d=1.5 outward: T-shape expands, 8 vertices.
        let result = PolylineOffset2D::new(t_shape_points(), -1.5, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(-1.5, -1.5, 0.0),
            Point3::new(11.5, -1.5, 0.0),
            Point3::new(11.5, 2.5, 0.0),
            Point3::new(8.5, 2.5, 0.0),
            Point3::new(8.5, 7.5, 0.0),
            Point3::new(1.5, 7.5, 0.0),
            Point3::new(1.5, 2.5, 0.0),
            Point3::new(-1.5, 2.5, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "T outward d=1.5");
    }

    // ── Cross-shape ground truth tests ────────────────────────────────
    //
    // Cross: horizontal arm (x:0..10, y:3..5, height=2),
    //        vertical arm   (x:3..7,  y:0..10, width=4).
    // 12 CCW vertices:
    //   (3,0)→(7,0)→(7,3)→(10,3)→(10,5)→(7,5)→(7,10)→(3,10)→(3,5)→(0,5)→(0,3)→(3,3)
    //
    // Horizontal arm collapses when d > 1.0 (arm height=2, so 2d > 2).
    // When arms collapse, only the vertical arm rectangle survives.

    /// Cross shape CCW vertices at origin.
    fn cross_shape_points() -> Vec<Point3> {
        vec![
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(7.0, 0.0, 0.0),
            Point3::new(7.0, 3.0, 0.0),
            Point3::new(10.0, 3.0, 0.0),
            Point3::new(10.0, 5.0, 0.0),
            Point3::new(7.0, 5.0, 0.0),
            Point3::new(7.0, 10.0, 0.0),
            Point3::new(3.0, 10.0, 0.0),
            Point3::new(3.0, 5.0, 0.0),
            Point3::new(0.0, 5.0, 0.0),
            Point3::new(0.0, 3.0, 0.0),
            Point3::new(3.0, 3.0, 0.0),
        ]
    }

    /// Open cross / plus — 4 arms from center with 180-degree reversals.
    ///
    /// Path: left → center → up → center → right → center → down.
    fn open_cross_points() -> Vec<Point3> {
        vec![
            Point3::new(-1.5, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.5, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.5, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, -1.5, 0.0),
        ]
    }

    // ── Inward offset (positive distance, CCW) ────────────────────────

    #[test]
    fn cross_inward_d05_full_shape() {
        // d=0.5: all arms survive → inset cross (12 vertices).
        let result = PolylineOffset2D::new(cross_shape_points(), 0.5, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(3.5, 0.5, 0.0),
            Point3::new(6.5, 0.5, 0.0),
            Point3::new(6.5, 3.5, 0.0),
            Point3::new(9.5, 3.5, 0.0),
            Point3::new(9.5, 4.5, 0.0),
            Point3::new(6.5, 4.5, 0.0),
            Point3::new(6.5, 9.5, 0.0),
            Point3::new(3.5, 9.5, 0.0),
            Point3::new(3.5, 4.5, 0.0),
            Point3::new(0.5, 4.5, 0.0),
            Point3::new(0.5, 3.5, 0.0),
            Point3::new(3.5, 3.5, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "cross inward d=0.5");
    }

    #[test]
    fn cross_inward_d15_arms_collapsed() {
        // d=1.5: horizontal arms collapse (2d=3 > height=2) → vertical arm rect.
        let result = PolylineOffset2D::new(cross_shape_points(), 1.5, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(4.5, 1.5, 0.0),
            Point3::new(5.5, 1.5, 0.0),
            Point3::new(5.5, 8.5, 0.0),
            Point3::new(4.5, 8.5, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "cross inward d=1.5");
    }

    // ── Outward offset (negative distance, CCW) ──────────────────────

    #[test]
    fn cross_outward_d05() {
        // d=0.5 outward: cross expands, 12 vertices.
        let result = PolylineOffset2D::new(cross_shape_points(), -0.5, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(2.5, -0.5, 0.0),
            Point3::new(7.5, -0.5, 0.0),
            Point3::new(7.5, 2.5, 0.0),
            Point3::new(10.5, 2.5, 0.0),
            Point3::new(10.5, 5.5, 0.0),
            Point3::new(7.5, 5.5, 0.0),
            Point3::new(7.5, 10.5, 0.0),
            Point3::new(2.5, 10.5, 0.0),
            Point3::new(2.5, 5.5, 0.0),
            Point3::new(-0.5, 5.5, 0.0),
            Point3::new(-0.5, 2.5, 0.0),
            Point3::new(2.5, 2.5, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "cross outward d=0.5");
    }

    #[test]
    fn cross_outward_d15() {
        // d=1.5 outward: cross expands, 12 vertices.
        let result = PolylineOffset2D::new(cross_shape_points(), -1.5, true)
            .execute()
            .unwrap();
        let expected = [
            Point3::new(1.5, -1.5, 0.0),
            Point3::new(8.5, -1.5, 0.0),
            Point3::new(8.5, 1.5, 0.0),
            Point3::new(11.5, 1.5, 0.0),
            Point3::new(11.5, 6.5, 0.0),
            Point3::new(8.5, 6.5, 0.0),
            Point3::new(8.5, 11.5, 0.0),
            Point3::new(1.5, 11.5, 0.0),
            Point3::new(1.5, 6.5, 0.0),
            Point3::new(-1.5, 6.5, 0.0),
            Point3::new(-1.5, 1.5, 0.0),
            Point3::new(1.5, 1.5, 0.0),
        ];
        assert_polygon_eq(&result, &expected, 1e-6, "cross outward d=1.5");
    }

    // ── Open cross with 180-degree reversals ──────────────────────────
    //
    // Open cross: left(-1.5,0) → center(0,0) → up(0,1.5) → center → right(1.5,0) → center → down(0,-1.5)
    // 6 segments with 3 flat caps (180° reversals at up/right/down tips).
    //
    // d>0 (left offset) traces the "upper/left" outline of all arms cleanly.
    // d<0 (right offset) traces the "lower/right" outline but currently has
    // self-crossings at the center due to flat-cap barriers preventing trimming.

    /// Expected closed 12-vertex cross outline at distance d from the open cross.
    ///
    /// Both d>0 and d<0 offsets should produce this same closed polygon.
    /// The offset of an open polyline traces both sides with flat caps,
    /// forming a closed outline.
    fn open_cross_expected(d: f64) -> Vec<Point3> {
        vec![
            Point3::new(-1.5, -d, 0.0),  //  0: left arm bottom-left
            Point3::new(-1.5, d, 0.0),    //  1: left arm top-left (cap)
            Point3::new(-d, d, 0.0),      //  2: center TL
            Point3::new(-d, 1.5, 0.0),    //  3: up arm top-left
            Point3::new(d, 1.5, 0.0),     //  4: up arm top-right (cap)
            Point3::new(d, d, 0.0),       //  5: center TR
            Point3::new(1.5, d, 0.0),     //  6: right arm top-right
            Point3::new(1.5, -d, 0.0),    //  7: right arm bottom-right (cap)
            Point3::new(d, -d, 0.0),      //  8: center BR
            Point3::new(d, -1.5, 0.0),    //  9: down arm bottom-right
            Point3::new(-d, -1.5, 0.0),   // 10: down arm bottom-left (cap)
            Point3::new(-d, -d, 0.0),     // 11: center BL
        ]
    }

    #[test]
    fn open_cross_d03_positive() {
        let result = PolylineOffset2D::new(open_cross_points(), 0.3, false)
            .execute()
            .unwrap();
        let expected = open_cross_expected(0.3);
        assert_polygon_eq(&result, &expected, 1e-6, "open cross d=+0.3");
    }

    #[test]
    fn open_cross_d03_negative() {
        let result = PolylineOffset2D::new(open_cross_points(), -0.3, false)
            .execute()
            .unwrap();
        let expected = open_cross_expected(0.3);
        assert_polygon_eq(&result, &expected, 1e-6, "open cross d=-0.3");
    }

    #[test]
    fn open_cross_d05_positive() {
        let result = PolylineOffset2D::new(open_cross_points(), 0.5, false)
            .execute()
            .unwrap();
        let expected = open_cross_expected(0.5);
        assert_polygon_eq(&result, &expected, 1e-6, "open cross d=+0.5");
    }

    #[test]
    fn open_cross_d05_negative() {
        let result = PolylineOffset2D::new(open_cross_points(), -0.5, false)
            .execute()
            .unwrap();
        let expected = open_cross_expected(0.5);
        assert_polygon_eq(&result, &expected, 1e-6, "open cross d=-0.5");
    }

}
