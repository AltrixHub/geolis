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

/// Triangulates a simple polygon by ear clipping its XY projection.
///
/// Returns one index triple per triangle — indices into `points`, each
/// triple wound counter-clockwise in XY whichever way the input winds.
/// **No Steiner points**: every triangle corner is an input vertex, and
/// the z of `points` is never read. That is what separates this from the
/// tessellation layer's constrained Delaunay pass, which answers in
/// positions rather than indices: a caller may lift each index back to
/// its own 3D position (a warped cap over a plan-projected ring) and
/// still derive shared topology from the index pairs, because the two
/// triangles meeting at an interior diagonal name it with the same pair.
///
/// Only strictly convex corners are clipped, so no triangle is
/// degenerate: the remaining polygon keeps positive area at every step,
/// and the final one is a non-degenerate triangle.
///
/// # Errors
///
/// Returns [`OperationError::InvalidInput`] when the polygon has fewer
/// than 3 vertices or its XY projection encloses no area, and
/// [`OperationError::Failed`] when no ear can be clipped — by the
/// two-ears theorem that means the projection is not a simple polygon
/// (self-intersecting, or pinched by a repeated vertex).
pub fn triangulate_polygon_xy(points: &[Point3]) -> Result<Vec<[usize; 3]>> {
    let n = points.len();
    if n < 3 {
        return Err(OperationError::InvalidInput(format!(
            "polygon triangulation needs at least 3 vertices, got {n}"
        ))
        .into());
    }
    let area = signed_area_2d(points);
    if area.abs() <= TOLERANCE {
        return Err(OperationError::InvalidInput(
            "polygon triangulation needs an XY projection that encloses area".into(),
        )
        .into());
    }

    // Walk counter-clockwise so a strictly positive turn is a convex
    // corner and every clipped ear comes out counter-clockwise too.
    let mut remaining: Vec<usize> = if area > 0.0 {
        (0..n).collect()
    } else {
        (0..n).rev().collect()
    };

    let mut triangles = Vec::with_capacity(n - 2);
    while remaining.len() > 3 {
        let m = remaining.len();
        let ear = find_ear(points, &remaining).ok_or_else(|| {
            OperationError::Failed(
                "polygon triangulation found no ear: the XY projection is not a simple polygon"
                    .into(),
            )
        })?;
        triangles.push([
            remaining[(ear + m - 1) % m],
            remaining[ear],
            remaining[(ear + 1) % m],
        ]);
        remaining.remove(ear);
    }
    triangles.push([remaining[0], remaining[1], remaining[2]]);
    Ok(triangles)
}

/// Finds the position in `remaining` of a clippable ear: a strictly
/// convex corner whose triangle holds no other remaining vertex. A
/// vertex ON the candidate triangle's boundary blocks the ear too —
/// clipping it would run the new diagonal through that vertex and leave
/// a T-junction where two triangles disagree about their shared edge.
fn find_ear(points: &[Point3], remaining: &[usize]) -> Option<usize> {
    let m = remaining.len();
    (0..m).find(|&k| {
        let prev = (k + m - 1) % m;
        let next = (k + 1) % m;
        let (a, b, c) = (
            points[remaining[prev]],
            points[remaining[k]],
            points[remaining[next]],
        );
        if orient_2d(&a, &b, &c) <= 0.0 {
            return false;
        }
        !remaining.iter().enumerate().any(|(j, &idx)| {
            j != k && j != prev && j != next && in_triangle_2d(&points[idx], &a, &b, &c)
        })
    })
}

/// Twice the signed area of triangle `a b c` in XY (positive when
/// counter-clockwise).
fn orient_2d(a: &Point3, b: &Point3, c: &Point3) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

/// Whether `p` lies inside or on the boundary of the counter-clockwise
/// triangle `a b c`, in XY.
fn in_triangle_2d(p: &Point3, a: &Point3, b: &Point3, c: &Point3) -> bool {
    orient_2d(a, b, p) >= 0.0 && orient_2d(b, c, p) >= 0.0 && orient_2d(c, a, p) >= 0.0
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

    // ── Ear-clipping triangulation ─────────────────────────────

    /// Sums the triangles' XY areas and asserts they tile the polygon:
    /// `n - 2` triangles, every one counter-clockwise, total area equal
    /// to the polygon's, and no index outside the ring.
    fn assert_tiles_polygon(points: &[Point3], triangles: &[[usize; 3]]) {
        assert_eq!(triangles.len(), points.len() - 2, "triangle count");
        let mut total = 0.0;
        for tri in triangles {
            for &i in tri {
                assert!(i < points.len(), "index {i} out of range");
            }
            assert_ne!(tri[0], tri[1]);
            assert_ne!(tri[1], tri[2]);
            assert_ne!(tri[2], tri[0]);
            let area2 = orient_2d(&points[tri[0]], &points[tri[1]], &points[tri[2]]);
            assert!(area2 > 0.0, "triangle {tri:?} is not counter-clockwise");
            total += 0.5 * area2;
        }
        let expected = signed_area_2d(points).abs();
        assert!(
            (total - expected).abs() < 1e-9,
            "triangles cover {total}, polygon is {expected}"
        );
    }

    fn ring(coords: &[(f64, f64)]) -> Vec<Point3> {
        // A varying z proves the triangulation reads the XY projection only.
        coords
            .iter()
            .enumerate()
            .map(|(i, &(x, y))| {
                #[allow(clippy::cast_precision_loss)]
                Point3::new(x, y, i as f64)
            })
            .collect()
    }

    #[test]
    fn triangulates_a_convex_ring_either_way_round() {
        let ccw = ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]);
        assert_tiles_polygon(&ccw, &triangulate_polygon_xy(&ccw).unwrap());

        let cw: Vec<Point3> = ccw.iter().rev().copied().collect();
        assert_tiles_polygon(&cw, &triangulate_polygon_xy(&cw).unwrap());
    }

    #[test]
    fn triangulates_a_non_convex_ring() {
        // L-shape: the reflex corner rules out a naive fan.
        let l = ring(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 2.0),
            (2.0, 2.0),
            (2.0, 4.0),
            (0.0, 4.0),
        ]);
        assert_tiles_polygon(&l, &triangulate_polygon_xy(&l).unwrap());
    }

    #[test]
    fn triangulates_a_thin_arc_band() {
        // Annular sector: the shape a curved ramp slab projects onto.
        let mut coords = Vec::new();
        for i in 0..=12 {
            let a = f64::from(i) * std::f64::consts::FRAC_PI_2 / 12.0;
            coords.push((3.0 * a.cos(), 3.0 * a.sin()));
        }
        for i in (0..=12).rev() {
            let a = f64::from(i) * std::f64::consts::FRAC_PI_2 / 12.0;
            coords.push((2.0 * a.cos(), 2.0 * a.sin()));
        }
        let band = ring(&coords);
        assert_tiles_polygon(&band, &triangulate_polygon_xy(&band).unwrap());
    }

    #[test]
    fn collinear_ring_vertices_never_produce_a_degenerate_triangle() {
        // Midpoints on every edge: they are never strictly convex, so
        // they may only be clipped once their neighbours are gone.
        let square = ring(&[
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (2.0, 2.0),
            (1.0, 2.0),
            (0.0, 2.0),
            (0.0, 1.0),
        ]);
        assert_tiles_polygon(&square, &triangulate_polygon_xy(&square).unwrap());
    }

    #[test]
    fn rejects_rings_that_are_not_simple_polygons() {
        let too_few = ring(&[(0.0, 0.0), (1.0, 0.0)]);
        assert!(triangulate_polygon_xy(&too_few).is_err());

        let collinear = ring(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]);
        assert!(matches!(
            triangulate_polygon_xy(&collinear),
            Err(crate::error::GeolisError::Operation(
                OperationError::InvalidInput(_)
            ))
        ));

        // Bow tie: the two lobes cancel, so the ring encloses no area.
        let bow_tie = ring(&[(0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 2.0)]);
        assert!(matches!(
            triangulate_polygon_xy(&bow_tie),
            Err(crate::error::GeolisError::Operation(
                OperationError::InvalidInput(_)
            ))
        ));

        // Pinched at a repeated vertex: area is positive but every
        // candidate ear is blocked by the pinch, so no ear exists.
        let pinched = ring(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (2.0, 2.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (2.0, 2.0),
        ]);
        assert!(matches!(
            triangulate_polygon_xy(&pinched),
            Err(crate::error::GeolisError::Operation(
                OperationError::Failed(_)
            ))
        ));
    }
}
