use super::{Point3, Vector3, TOLERANCE};
use crate::error::{OperationError, Result};

/// Computes the signed area of a polygon in the XY plane (shoelace formula).
///
/// Positive for counter-clockwise, negative for clockwise.
#[must_use]
pub fn signed_area_2d(points: &[Point3]) -> f64 {
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

/// Rotates a closed polygon so it starts at the leftmost vertex (smallest x),
/// breaking ties by smallest y. Ensures deterministic output for tests.
#[must_use]
pub fn rotate_to_canonical_start(points: &[Point3]) -> Vec<Point3> {
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

/// Returns the leftmost-bottommost vertex of a polygon (for tie-breaking in sort).
#[must_use]
pub fn leftmost_bottom(points: &[Point3]) -> Point3 {
    let mut best = points[0];
    for &pt in &points[1..] {
        if pt.x < best.x - TOLERANCE || ((pt.x - best.x).abs() < TOLERANCE && pt.y < best.y) {
            best = pt;
        }
    }
    best
}

/// Computes the normalized direction from point `a` to point `b`.
///
/// # Errors
///
/// Returns `OperationError::InvalidInput` if the segment has zero length.
pub fn segment_direction(a: &Point3, b: &Point3) -> Result<Vector3> {
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
#[must_use]
pub fn left_normal(dir: Vector3) -> Vector3 {
    Vector3::new(-dir.y, dir.x, 0.0)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn signed_area_ccw_square() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let area = signed_area_2d(&pts);
        assert!((area - 1.0).abs() < TOLERANCE);
    }

    #[test]
    fn signed_area_cw_square() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let area = signed_area_2d(&pts);
        assert!((area + 1.0).abs() < TOLERANCE);
    }

    #[test]
    fn signed_area_degenerate() {
        assert!((signed_area_2d(&[Point3::new(0.0, 0.0, 0.0)])).abs() < TOLERANCE);
        assert!((signed_area_2d(&[])).abs() < TOLERANCE);
    }

    #[test]
    fn canonical_start_already_leftmost() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ];
        let rotated = rotate_to_canonical_start(&pts);
        assert!((rotated[0].x).abs() < TOLERANCE);
        assert!((rotated[0].y).abs() < TOLERANCE);
    }

    #[test]
    fn canonical_start_rotation() {
        let pts = vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
        ];
        let rotated = rotate_to_canonical_start(&pts);
        assert!((rotated[0].x).abs() < TOLERANCE);
        assert!((rotated[0].y).abs() < TOLERANCE);
    }

    #[test]
    fn leftmost_bottom_basic() {
        let pts = vec![
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(0.5, 0.5, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ];
        let lb = leftmost_bottom(&pts);
        assert!((lb.x - 0.5).abs() < TOLERANCE);
        assert!((lb.y - 0.5).abs() < TOLERANCE);
    }

    #[test]
    fn segment_direction_basic() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(3.0, 4.0, 0.0);
        let dir = segment_direction(&a, &b).unwrap();
        assert!((dir.x - 0.6).abs() < TOLERANCE);
        assert!((dir.y - 0.8).abs() < TOLERANCE);
    }

    #[test]
    fn segment_direction_zero_length() {
        let a = Point3::new(1.0, 1.0, 0.0);
        let b = Point3::new(1.0, 1.0, 0.0);
        assert!(segment_direction(&a, &b).is_err());
    }

    #[test]
    fn left_normal_basic() {
        let dir = Vector3::new(1.0, 0.0, 0.0);
        let n = left_normal(dir);
        assert!((n.x).abs() < TOLERANCE);
        assert!((n.y - 1.0).abs() < TOLERANCE);
    }
}
