mod filter;
mod raw_offset;
mod self_intersect;
mod slice;
mod stitch;

use crate::error::{OperationError, Result};
use crate::geometry::pline::{Pline, PlineVertex};

/// Offsets a polyline (with potential arc segments) using the slice-and-filter
/// algorithm.
///
/// This is the new implementation that replaces the recursive split approach
/// in `PolylineOffset2D` with a more robust cavalier_contours-style algorithm.
#[derive(Debug)]
pub struct PlineOffset2D {
    pline: Pline,
    distance: f64,
}

impl PlineOffset2D {
    /// Creates a new polyline offset operation.
    #[must_use]
    pub fn new(pline: Pline, distance: f64) -> Self {
        Self { pline, distance }
    }

    /// Executes the offset, returning one or more result polylines.
    ///
    /// # Errors
    ///
    /// Returns `OperationError::InvalidInput` if the polyline has fewer than
    /// 2 vertices, or `OperationError::Failed` if the offset collapses entirely.
    pub fn execute(&self) -> Result<Vec<Pline>> {
        if self.pline.vertices.len() < 2 {
            return Err(OperationError::InvalidInput(
                "at least 2 vertices required for pline offset".to_owned(),
            )
            .into());
        }

        if self.distance.abs() < crate::math::TOLERANCE {
            return Ok(vec![self.pline.clone()]);
        }

        if self.pline.closed {
            self.execute_closed()
        } else {
            self.execute_open()
        }
    }

    /// Executes offset for closed polylines using the standard slice-and-filter
    /// pipeline.
    fn execute_closed(&self) -> Result<Vec<Pline>> {
        // Step 1: Build raw offset polyline.
        let raw = raw_offset::build(&self.pline, self.distance)?;

        // Step 2: Find all self-intersections.
        let intersections = self_intersect::find_all(&raw);
        if intersections.is_empty() {
            return Ok(vec![raw]);
        }

        // Step 3: Slice at intersection points.
        let seg_count = raw.segment_count();
        let slices = slice::build(&raw.vertices, seg_count, &intersections);

        // Step 4: Filter slices by distance to original.
        let valid = filter::apply(&slices, &self.pline, self.distance);

        // Step 5: Stitch valid slices into result polylines.
        let result = stitch::connect(&valid);

        if result.is_empty() {
            return Err(OperationError::Failed(
                "offset collapsed completely".to_owned(),
            )
            .into());
        }

        Ok(result)
    }

    /// Executes offset for open polylines using the both-sides approach.
    ///
    /// Creates a closed buffer polygon by:
    /// 1. Detecting spoke patterns → analytical outline (fast path)
    /// 2. Otherwise: forward + backward raw offset, combined into a closed polygon,
    ///    then slice-and-filter to extract the outer boundary
    fn execute_open(&self) -> Result<Vec<Pline>> {
        let abs_d = self.distance.abs();

        // Fast path: spoke patterns (center visited multiple times).
        if let Some(center) = find_spoke_center(&self.pline) {
            let tips = extract_arm_tips(&self.pline, &center);
            if tips.len() >= 2 {
                let outline = build_spoke_outline(&center, &tips, abs_d)?;
                return Ok(vec![outline]);
            }
        }

        // General path: wall outline algorithm (handles arbitrary centerline networks).
        let wall = super::wall_outline::WallOutline2D::new(self.pline.clone(), abs_d);
        wall.execute()
    }
}

#[cfg(test)]
fn clean_pline(pline: &Pline) -> Pline {
    let tol_sq = crate::math::TOLERANCE * crate::math::TOLERANCE * 100.0;
    let verts = &pline.vertices;
    if verts.len() < 3 {
        return pline.clone();
    }

    // Step 1: Remove consecutive near-duplicates.
    let mut deduped = Vec::with_capacity(verts.len());
    for v in verts {
        if let Some(last) = deduped.last() {
            let last: &PlineVertex = last;
            let dx = v.x - last.x;
            let dy = v.y - last.y;
            if dx * dx + dy * dy < tol_sq {
                continue;
            }
        }
        deduped.push(*v);
    }
    // Wrap-around: check last vs first.
    if pline.closed && deduped.len() > 1 {
        let first = &deduped[0];
        let last = &deduped[deduped.len() - 1];
        let dx = last.x - first.x;
        let dy = last.y - first.y;
        if dx * dx + dy * dy < tol_sq {
            deduped.pop();
        }
    }

    if deduped.len() < 3 {
        return Pline {
            vertices: deduped,
            closed: pline.closed,
        };
    }

    // Step 2: Remove collinear LINE vertices (single pass, cross-product check).
    // Arc vertices are always kept since arcs are not collinear even if endpoints align.
    let n = deduped.len();
    let mut cleaned = Vec::with_capacity(n);
    for i in 0..n {
        let prev = if i == 0 { n - 1 } else { i - 1 };
        let next = (i + 1) % n;

        // Keep vertices that participate in arc segments.
        if deduped[i].bulge.abs() >= 1e-12 || deduped[prev].bulge.abs() >= 1e-12 {
            cleaned.push(deduped[i]);
            continue;
        }

        let cross = (deduped[i].x - deduped[prev].x) * (deduped[next].y - deduped[i].y)
            - (deduped[i].y - deduped[prev].y) * (deduped[next].x - deduped[i].x);
        if cross.abs() >= crate::math::TOLERANCE {
            cleaned.push(deduped[i]);
        }
    }

    // Don't reduce below 3 vertices — return deduped as fallback.
    if cleaned.len() < 3 {
        return Pline {
            vertices: deduped,
            closed: pline.closed,
        };
    }

    Pline {
        vertices: cleaned,
        closed: pline.closed,
    }
}

#[cfg(test)]
fn combine_both_sides(forward: &Pline, backward: &Pline) -> Pline {
    let tol_sq = 1e-10;
    let mut verts = Vec::with_capacity(forward.vertices.len() + backward.vertices.len());

    // Forward path (all vertices keep their bulges for arcs/lines to next vertex).
    for (i, v) in forward.vertices.iter().enumerate() {
        if i == forward.vertices.len() - 1 {
            // Last forward vertex: its bulge should be 0 (straight line to first
            // backward vertex — the flat cap at the end of the original path).
            verts.push(PlineVertex::line(v.x, v.y));
        } else {
            verts.push(*v);
        }
    }

    // Backward path: skip first/last vertices if they duplicate the connection points.
    let back_len = backward.vertices.len();
    let skip_first = if let (Some(f_last), Some(b_first)) =
        (forward.vertices.last(), backward.vertices.first())
    {
        let dx = f_last.x - b_first.x;
        let dy = f_last.y - b_first.y;
        dx * dx + dy * dy < tol_sq
    } else {
        false
    };
    let skip_last = if let (Some(f_first), Some(b_last)) =
        (forward.vertices.first(), backward.vertices.last())
    {
        let dx = f_first.x - b_last.x;
        let dy = f_first.y - b_last.y;
        dx * dx + dy * dy < tol_sq
    } else {
        false
    };

    let start = if skip_first { 1 } else { 0 };
    let end = if skip_last { back_len - 1 } else { back_len };

    for i in start..end {
        let v = &backward.vertices[i];
        if i == end - 1 {
            // Last backward vertex: straight line back to first forward vertex (start cap).
            verts.push(PlineVertex::line(v.x, v.y));
        } else {
            verts.push(*v);
        }
    }

    Pline {
        vertices: verts,
        closed: true,
    }
}

// ── Spoke-based offset for open polylines with a center point ──────

/// Finds a "center" vertex that appears multiple times in the polyline.
///
/// Used to detect spoke patterns where arms radiate from a shared center
/// (e.g., open cross, X-cross, fork).
fn find_spoke_center(pline: &Pline) -> Option<(f64, f64)> {
    let tol_sq = 1e-8;
    let verts = &pline.vertices;
    for i in 0..verts.len() {
        for j in (i + 1)..verts.len() {
            let dx = verts[i].x - verts[j].x;
            let dy = verts[i].y - verts[j].y;
            if dx * dx + dy * dy < tol_sq {
                return Some((verts[i].x, verts[i].y));
            }
        }
    }
    None
}

/// Extracts unique arm tip positions (vertices not at center) from a spoke polyline.
fn extract_arm_tips(pline: &Pline, center: &(f64, f64)) -> Vec<(f64, f64)> {
    let tol_sq = 1e-8;
    let mut tips = Vec::new();
    for v in &pline.vertices {
        let dx = v.x - center.0;
        let dy = v.y - center.1;
        if dx * dx + dy * dy < tol_sq {
            continue;
        }
        let is_dup = tips
            .iter()
            .any(|t: &(f64, f64)| (v.x - t.0).powi(2) + (v.y - t.1).powi(2) < tol_sq);
        if !is_dup {
            tips.push((v.x, v.y));
        }
    }
    tips
}

/// Builds a closed polygon outline around spoke arms at distance `d`.
///
/// For each arm, emits left-tip and right-tip at the arm endpoint, then a
/// center vertex connecting to the next arm. Arms are sorted by angle
/// descending (clockwise).
fn build_spoke_outline(
    center: &(f64, f64),
    tips: &[(f64, f64)],
    d: f64,
) -> Result<Pline> {
    use crate::math::intersect_2d::{line_line_intersect_2d, point_at};
    use crate::math::{Point3, Vector3};

    // Compute arm directions and angles.
    let mut arms: Vec<((f64, f64), (f64, f64), f64)> = Vec::with_capacity(tips.len());
    for &tip in tips {
        let dx = tip.0 - center.0;
        let dy = tip.1 - center.1;
        let len = (dx * dx + dy * dy).sqrt();
        if len < 1e-12 {
            continue;
        }
        let dir = (dx / len, dy / len);
        let angle = dir.1.atan2(dir.0);
        arms.push((tip, dir, angle));
    }

    if arms.len() < 2 {
        return Err(OperationError::Failed("spoke needs at least 2 arms".to_owned()).into());
    }

    // Sort by angle descending for CW outline.
    arms.sort_by(|a, b| {
        b.2.partial_cmp(&a.2)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let n = arms.len();
    let mut outline = Vec::with_capacity(3 * n);

    for i in 0..n {
        let (tip, dir, _) = &arms[i];
        let (_, next_dir, _) = &arms[(i + 1) % n];

        // Left normal: (-dy, dx).
        let ln = (-dir.1, dir.0);

        // Left-tip (left side when facing center→tip direction).
        outline.push(PlineVertex::line(tip.0 + d * ln.0, tip.1 + d * ln.1));
        // Right-tip (right side).
        outline.push(PlineVertex::line(tip.0 - d * ln.0, tip.1 - d * ln.1));

        // Center vertex: intersection of arm_i's right buffer edge
        // and arm_{i+1}'s left buffer edge.
        let ln_next = (-next_dir.1, next_dir.0);
        let base_right = Point3::new(center.0 - d * ln.0, center.1 - d * ln.1, 0.0);
        let base_left = Point3::new(
            center.0 + d * ln_next.0,
            center.1 + d * ln_next.1,
            0.0,
        );
        let d_right = Vector3::new(dir.0, dir.1, 0.0);
        let d_left = Vector3::new(next_dir.0, next_dir.1, 0.0);

        if let Some((t, _)) = line_line_intersect_2d(&base_right, &d_right, &base_left, &d_left) {
            let pt = point_at(&base_right, &d_right, t);
            outline.push(PlineVertex::line(pt.x, pt.y));
        } else {
            // Parallel arms: use midpoint.
            outline.push(PlineVertex::line(
                (base_right.x + base_left.x) * 0.5,
                (base_right.y + base_left.y) * 0.5,
            ));
        }
    }

    Ok(Pline {
        vertices: outline,
        closed: true,
    })
}

#[cfg(test)]
fn signed_area_pline(pline: &Pline) -> f64 {
    let n = pline.vertices.len();
    if n < 3 {
        return 0.0;
    }
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += pline.vertices[i].x * pline.vertices[j].y;
        area -= pline.vertices[j].x * pline.vertices[i].y;
    }
    area * 0.5
}

#[cfg(test)]
fn trim_recursive(pline: &Pline, winding_sign: f64) -> Vec<Pline> {
    let cleaned = clean_pline(pline);
    if cleaned.vertices.len() < 4 {
        return vec![cleaned];
    }

    // Find first self-intersection.
    let intersections = self_intersect::find_all(&cleaned);
    if intersections.is_empty() {
        return vec![cleaned];
    }
    let ix = &intersections[0];

    // Split at the intersection.
    let (part_a, part_b) = split_pline_at_intersection(&cleaned, ix);

    // Recursively trim each part.
    let parts_a = trim_recursive(&part_a, winding_sign);
    let parts_b = trim_recursive(&part_b, winding_sign);

    // Keep only parts with correct winding.
    let mut result: Vec<Pline> = parts_a
        .into_iter()
        .chain(parts_b)
        .filter(|p| signed_area_pline(p) * winding_sign > 0.0)
        .collect();

    if result.is_empty() {
        // Fallback: recursively trim and keep largest by area.
        let fallback_a = trim_recursive(&part_a, winding_sign);
        let fallback_b = trim_recursive(&part_b, winding_sign);
        if let Some(best) = fallback_a
            .into_iter()
            .chain(fallback_b)
            .max_by(|x, y| {
                signed_area_pline(x)
                    .abs()
                    .partial_cmp(&signed_area_pline(y).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
        {
            result.push(best);
        }
    }

    result
}

#[cfg(test)]
fn split_pline_at_intersection(
    pline: &Pline,
    ix: &self_intersect::Intersection,
) -> (Pline, Pline) {
    let n = pline.vertices.len();
    let (pt_x, pt_y) = ix.point;

    // Sub-path A: intersection → v[seg_i+1] → ... → v[seg_j] → (close back to intersection)
    let mut a_verts = Vec::new();
    // Intersection vertex: line to next (sub-arc not yet supported).
    a_verts.push(PlineVertex::line(pt_x, pt_y));
    let mut idx = (ix.seg_i + 1) % n;
    loop {
        if idx == ix.seg_j {
            // Last vertex: its closing segment goes to intersection, not original next.
            a_verts.push(PlineVertex::line(
                pline.vertices[idx].x,
                pline.vertices[idx].y,
            ));
            break;
        }
        // Intermediate vertices: preserve original bulge.
        a_verts.push(pline.vertices[idx]);
        idx = (idx + 1) % n;
    }

    // Sub-path B: intersection → v[seg_j+1] → ... → v[seg_i] → (close back to intersection)
    let mut b_verts = Vec::new();
    b_verts.push(PlineVertex::line(pt_x, pt_y));
    idx = (ix.seg_j + 1) % n;
    loop {
        if idx == ix.seg_i {
            b_verts.push(PlineVertex::line(
                pline.vertices[idx].x,
                pline.vertices[idx].y,
            ));
            break;
        }
        b_verts.push(pline.vertices[idx]);
        idx = (idx + 1) % n;
    }

    let part_a = Pline {
        vertices: a_verts,
        closed: true,
    };
    let part_b = Pline {
        vertices: b_verts,
        closed: true,
    };

    (part_a, part_b)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn square_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(10.0, 10.0),
                PlineVertex::line(0.0, 10.0),
            ],
            closed: true,
        }
    }

    #[test]
    fn square_inward_offset() {
        let op = PlineOffset2D::new(square_pline(), 1.0);
        let result = op.execute().unwrap();
        assert!(!result.is_empty(), "should produce at least one result");
        // The inward offset of a 10x10 square by 1 should be an 8x8 square.
        let poly = &result[0];
        assert_eq!(poly.vertices.len(), 4, "expected 4 vertices");
    }

    #[test]
    fn square_outward_offset() {
        let op = PlineOffset2D::new(square_pline(), -1.0);
        let result = op.execute().unwrap();
        assert!(!result.is_empty());
        let poly = &result[0];
        assert_eq!(poly.vertices.len(), 4, "expected 4 vertices");
    }

    #[test]
    fn no_self_intersection_passthrough() {
        // A simple triangle offset inward — no self-intersections expected.
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(5.0, 8.66),
            ],
            closed: true,
        };
        let op = PlineOffset2D::new(pline, 0.5);
        let result = op.execute().unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vertices.len(), 3);
    }

    // ── Arc segment tests ──

    /// A rounded rectangle: two horizontal lines connected by semicircular arcs.
    /// (0,0)→(10,0) line, (10,0)→(10,4) CCW semicircle, (10,4)→(0,4) line, (0,4)→(0,0) CCW semicircle.
    fn rounded_rect_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),   // seg 0: line →(10,0)
                PlineVertex::new(10.0, 0.0, 1.0), // seg 1: CCW semicircle →(10,4)
                PlineVertex::line(10.0, 4.0),   // seg 2: line →(0,4)
                PlineVertex::new(0.0, 4.0, 1.0),  // seg 3: CCW semicircle →(0,0)
            ],
            closed: true,
        }
    }

    #[test]
    fn rounded_rect_inward_offset() {
        let op = PlineOffset2D::new(rounded_rect_pline(), 0.5);
        let result = op.execute().unwrap();
        assert!(!result.is_empty(), "should produce at least one result");
        let poly = &result[0];
        // Should have arc segments (non-zero bulge) in the result.
        let has_arcs = poly.vertices.iter().any(|v| v.bulge.abs() > 1e-6);
        assert!(has_arcs, "result should contain arc segments");
    }

    #[test]
    fn rounded_rect_outward_offset() {
        let op = PlineOffset2D::new(rounded_rect_pline(), -0.5);
        let result = op.execute().unwrap();
        assert!(!result.is_empty(), "should produce at least one result");
    }

    #[test]
    fn semicircle_arc_no_self_intersect() {
        // Open pline with a single semicircular arc — should pass through with no intersections.
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::new(5.0, 0.0, 1.0),
                PlineVertex::line(10.0, 0.0),
            ],
            closed: false,
        };
        let op = PlineOffset2D::new(pline, 0.5);
        let result = op.execute().unwrap();
        assert_eq!(result.len(), 1);
    }

    // ── Open polyline tests (both-sides buffer) ──

    /// Checks that every expected point appears somewhere in the result vertices.
    fn assert_vertices_match(result: &[PlineVertex], expected: &[(f64, f64)], tol: f64) {
        assert_eq!(
            result.len(),
            expected.len(),
            "vertex count mismatch: got {}, expected {}",
            result.len(),
            expected.len()
        );
        for &(ex, ey) in expected {
            let found = result
                .iter()
                .any(|v| (v.x - ex).abs() < tol && (v.y - ey).abs() < tol);
            assert!(
                found,
                "expected vertex ({ex:.4}, {ey:.4}) not found in result"
            );
        }
    }

    /// Open cross: 4 arms from center with 180° reversals.
    fn open_cross_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(-1.5, 0.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(0.0, 1.5),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(1.5, 0.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(0.0, -1.5),
            ],
            closed: false,
        }
    }

    /// Expected cross outline at distance d (12 vertices).
    fn open_cross_expected(d: f64) -> Vec<(f64, f64)> {
        vec![
            (-1.5, -d),
            (-1.5, d),
            (-d, d),
            (-d, 1.5),
            (d, 1.5),
            (d, d),
            (1.5, d),
            (1.5, -d),
            (d, -d),
            (d, -1.5),
            (-d, -1.5),
            (-d, -d),
        ]
    }

    #[test]
    fn open_cross_d03_positive() {
        let op = PlineOffset2D::new(open_cross_pline(), 0.3);
        let result = op.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 closed polygon");
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        assert_vertices_match(&poly.vertices, &open_cross_expected(0.3), 0.05);
    }

    #[test]
    fn open_cross_d03_negative() {
        let op = PlineOffset2D::new(open_cross_pline(), -0.3);
        let result = op.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 closed polygon");
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        assert_vertices_match(&poly.vertices, &open_cross_expected(0.3), 0.05);
    }

    #[test]
    fn open_cross_d05() {
        let op = PlineOffset2D::new(open_cross_pline(), 0.5);
        let result = op.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 closed polygon");
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        assert_vertices_match(&poly.vertices, &open_cross_expected(0.5), 0.05);
    }

    /// X-cross: 2 diagonal lines crossing at center with 180° reversals.
    fn x_cross_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(-3.0, -3.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(3.0, 3.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(-3.0, 3.0),
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(3.0, -3.0),
            ],
            closed: false,
        }
    }

    /// Expected X-cross outline at distance d (12 vertices).
    fn x_cross_expected(a: f64, d: f64) -> Vec<(f64, f64)> {
        let s2 = std::f64::consts::SQRT_2;
        let h = d * s2 / 2.0;
        let d2 = d * s2;
        vec![
            (-a - h, -a + h),
            (-a + h, -a - h),
            (0.0, -d2),
            (a - h, -a - h),
            (a + h, -a + h),
            (d2, 0.0),
            (a + h, a - h),
            (a - h, a + h),
            (0.0, d2),
            (-a + h, a + h),
            (-a - h, a - h),
            (-d2, 0.0),
        ]
    }

    #[test]
    fn x_cross_d05() {
        let op = PlineOffset2D::new(x_cross_pline(), 0.5);
        let result = op.execute().unwrap();
        let expected = x_cross_expected(3.0, 0.5);
        assert_eq!(result.len(), 1, "expected 1 closed polygon");
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        assert_vertices_match(&poly.vertices, &expected, 0.05);
    }

    /// Fork (Y-shape): stem + 2 branches with reversal at junction.
    fn fork_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(5.0, 4.0),
                PlineVertex::line(0.0, 9.0),
                PlineVertex::line(5.0, 4.0),
                PlineVertex::line(10.0, 9.0),
            ],
            closed: false,
        }
    }

    /// Expected fork outline at distance d (9 vertices).
    fn fork_expected(d: f64) -> Vec<(f64, f64)> {
        let s2 = std::f64::consts::SQRT_2;
        let h = d * s2 / 2.0;
        let jy = 4.0 + d * (1.0 - s2);
        vec![
            (5.0 - d, 0.0),
            (5.0 + d, 0.0),
            (5.0 + d, jy),
            (10.0 + h, 9.0 - h),
            (10.0 - h, 9.0 + h),
            (5.0, 4.0 + d * s2),
            (h, 9.0 + h),
            (-h, 9.0 - h),
            (5.0 - d, jy),
        ]
    }

    #[test]
    fn fork_d05() {
        let op = PlineOffset2D::new(fork_pline(), 0.5);
        let result = op.execute().unwrap();
        assert_eq!(result.len(), 1, "expected 1 closed polygon");
        let poly = &result[0];
        assert!(poly.closed, "result should be closed");
        assert_vertices_match(&poly.vertices, &fork_expected(0.5), 0.05);
    }

    /// Double-cross (井-shape): 4 crossing line segments.
    fn double_cross_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(3.0, 0.0),
                PlineVertex::line(3.0, 10.0),
                PlineVertex::line(3.0, 7.0),
                PlineVertex::line(0.0, 7.0),
                PlineVertex::line(10.0, 7.0),
                PlineVertex::line(7.0, 7.0),
                PlineVertex::line(7.0, 10.0),
                PlineVertex::line(7.0, 0.0),
                PlineVertex::line(7.0, 3.0),
                PlineVertex::line(10.0, 3.0),
                PlineVertex::line(0.0, 3.0),
            ],
            closed: false,
        }
    }

    /// Expected double-cross outline at distance d (28 vertices).
    fn double_cross_expected(d: f64) -> Vec<(f64, f64)> {
        vec![
            (3.0 - d, 0.0),
            (3.0 + d, 0.0),
            (3.0 + d, 3.0 - d),
            (7.0 - d, 3.0 - d),
            (7.0 - d, 0.0),
            (7.0 + d, 0.0),
            (7.0 + d, 3.0 - d),
            (10.0, 3.0 - d),
            (10.0, 3.0 + d),
            (7.0 + d, 3.0 + d),
            (7.0 + d, 7.0 - d),
            (10.0, 7.0 - d),
            (10.0, 7.0 + d),
            (7.0 + d, 7.0 + d),
            (7.0 + d, 10.0),
            (7.0 - d, 10.0),
            (7.0 - d, 7.0 + d),
            (3.0 + d, 7.0 + d),
            (3.0 + d, 10.0),
            (3.0 - d, 10.0),
            (3.0 - d, 7.0 + d),
            (0.0, 7.0 + d),
            (0.0, 7.0 - d),
            (3.0 - d, 7.0 - d),
            (3.0 - d, 3.0 + d),
            (0.0, 3.0 + d),
            (0.0, 3.0 - d),
            (3.0 - d, 3.0 - d),
        ]
    }

    #[test]
    fn double_cross_d03() {
        let op = PlineOffset2D::new(double_cross_pline(), 0.3);
        let result = op.execute().unwrap();
        assert!(!result.is_empty(), "expected at least 1 polygon");
        let outer = result
            .iter()
            .max_by_key(|p| p.vertices.len())
            .unwrap();
        assert!(outer.closed, "outer boundary should be closed");
        assert_vertices_match(&outer.vertices, &double_cross_expected(0.3), 0.05);
    }

    #[test]
    fn double_cross_d08() {
        let op = PlineOffset2D::new(double_cross_pline(), 0.8);
        let result = op.execute().unwrap();
        assert!(
            !result.is_empty(),
            "expected at least 1 polygon"
        );
        let outer = result
            .iter()
            .max_by_key(|p| p.vertices.len())
            .unwrap();
        assert!(outer.closed, "outer boundary should be closed");
        assert_vertices_match(&outer.vertices, &double_cross_expected(0.8), 0.05);
    }

    #[test]
    fn debug_semicircle_combined() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::new(5.0, 0.0, 1.0),
                PlineVertex::line(10.0, 0.0),
            ],
            closed: false,
        };
        let abs_d = 0.5_f64;
        let forward = raw_offset::build(&pline, abs_d).unwrap();
        let reversed_input = pline.reversed();
        let backward = raw_offset::build(&reversed_input, abs_d).unwrap();
        let combined = combine_both_sides(&forward, &backward);
        let cleaned = clean_pline(&combined);

        eprintln!("=== Semicircle d=0.5 ===");
        eprintln!("Forward ({} verts):", forward.vertices.len());
        for (i, v) in forward.vertices.iter().enumerate() {
            eprintln!("  [{i}] ({:.4}, {:.4}) bulge={:.6}", v.x, v.y, v.bulge);
        }
        eprintln!("Backward ({} verts):", backward.vertices.len());
        for (i, v) in backward.vertices.iter().enumerate() {
            eprintln!("  [{i}] ({:.4}, {:.4}) bulge={:.6}", v.x, v.y, v.bulge);
        }
        eprintln!("Combined ({} verts):", combined.vertices.len());
        for (i, v) in combined.vertices.iter().enumerate() {
            eprintln!("  [{i}] ({:.4}, {:.4}) bulge={:.6}", v.x, v.y, v.bulge);
        }
        eprintln!("Cleaned ({} verts):", cleaned.vertices.len());
        for (i, v) in cleaned.vertices.iter().enumerate() {
            eprintln!("  [{i}] ({:.4}, {:.4}) bulge={:.6}", v.x, v.y, v.bulge);
        }

        let intersections = self_intersect::find_all(&cleaned);
        eprintln!("Intersections ({}):", intersections.len());
        for ix in &intersections {
            eprintln!(
                "  seg {} x seg {}: t=({:.4}, {:.4}) pt=({:.4}, {:.4})",
                ix.seg_i, ix.seg_j, ix.t_i, ix.t_j, ix.point.0, ix.point.1
            );
        }

        let winding = signed_area_pline(&cleaned);
        eprintln!("Winding area: {winding:.4}");

        let result = trim_recursive(&cleaned, if winding.abs() < 1e-10 { 1.0 } else { winding.signum() });
        eprintln!("Trim result: {} polygons", result.len());
        for (i, p) in result.iter().enumerate() {
            eprintln!("  Poly {i}: {} verts, area={:.4}", p.vertices.len(), signed_area_pline(p));
        }
    }

    #[test]
    fn debug_fork_combined() {
        let pline = fork_pline();
        let abs_d = 0.5_f64;
        let forward = raw_offset::build(&pline, abs_d).unwrap();
        let reversed_input = pline.reversed();
        let backward = raw_offset::build(&reversed_input, abs_d).unwrap();
        let combined = combine_both_sides(&forward, &backward);
        let cleaned = clean_pline(&combined);

        eprintln!("=== Fork d=0.5 ===");
        eprintln!("Forward ({} verts):", forward.vertices.len());
        for (i, v) in forward.vertices.iter().enumerate() {
            eprintln!("  [{i}] ({:.4}, {:.4}) bulge={:.6}", v.x, v.y, v.bulge);
        }
        eprintln!("Backward ({} verts):", backward.vertices.len());
        for (i, v) in backward.vertices.iter().enumerate() {
            eprintln!("  [{i}] ({:.4}, {:.4}) bulge={:.6}", v.x, v.y, v.bulge);
        }
        eprintln!("Combined ({} verts):", combined.vertices.len());
        for (i, v) in combined.vertices.iter().enumerate() {
            eprintln!("  [{i}] ({:.4}, {:.4}) bulge={:.6}", v.x, v.y, v.bulge);
        }
        eprintln!("Cleaned ({} verts):", cleaned.vertices.len());
        for (i, v) in cleaned.vertices.iter().enumerate() {
            eprintln!("  [{i}] ({:.4}, {:.4}) bulge={:.6}", v.x, v.y, v.bulge);
        }

        let intersections = self_intersect::find_all(&cleaned);
        eprintln!("Intersections ({}):", intersections.len());
        for ix in &intersections {
            eprintln!(
                "  seg {} x seg {}: t=({:.4}, {:.4}) pt=({:.4}, {:.4})",
                ix.seg_i, ix.seg_j, ix.t_i, ix.t_j, ix.point.0, ix.point.1
            );
        }

        let winding = signed_area_pline(&cleaned);
        eprintln!("Winding area: {winding:.4}");
    }

    #[test]
    fn debug_x_cross_combined() {
        let pline = x_cross_pline();
        let abs_d = 0.5_f64;
        let forward = raw_offset::build(&pline, abs_d).unwrap();
        let reversed_input = pline.reversed();
        let backward = raw_offset::build(&reversed_input, abs_d).unwrap();
        let combined = combine_both_sides(&forward, &backward);
        let cleaned = clean_pline(&combined);

        eprintln!("=== X-cross d=0.5 ===");
        eprintln!("Combined ({} verts):", combined.vertices.len());
        for (i, v) in combined.vertices.iter().enumerate() {
            eprintln!("  [{i}] ({:.4}, {:.4}) bulge={:.6}", v.x, v.y, v.bulge);
        }
        eprintln!("Cleaned ({} verts):", cleaned.vertices.len());
        for (i, v) in cleaned.vertices.iter().enumerate() {
            eprintln!("  [{i}] ({:.4}, {:.4}) bulge={:.6}", v.x, v.y, v.bulge);
        }

        let intersections = self_intersect::find_all(&cleaned);
        eprintln!("Intersections ({}):", intersections.len());
        for ix in &intersections {
            eprintln!(
                "  seg {} x seg {}: t=({:.4}, {:.4}) pt=({:.4}, {:.4})",
                ix.seg_i, ix.seg_j, ix.t_i, ix.t_j, ix.point.0, ix.point.1
            );
        }
    }

    #[test]
    fn debug_double_cross_combined() {
        use crate::math::distance_2d::point_to_segment_dist;

        let pline = double_cross_pline();
        let abs_d = 0.3_f64;
        let forward = raw_offset::build(&pline, abs_d).unwrap();
        let reversed_input = pline.reversed();
        let backward = raw_offset::build(&reversed_input, abs_d).unwrap();
        let combined = combine_both_sides(&forward, &backward);
        let cleaned = clean_pline(&combined);

        eprintln!("=== Double-cross d=0.3 slice-and-filter ===");
        eprintln!("Cleaned ({} verts):", cleaned.vertices.len());
        for (i, v) in cleaned.vertices.iter().enumerate() {
            eprintln!("  [{i}] ({:.4}, {:.4})", v.x, v.y);
        }

        let intersections = self_intersect::find_all(&cleaned);
        eprintln!("Intersections ({}):", intersections.len());

        let seg_count = cleaned.segment_count();
        let slices = slice::build(&cleaned.vertices, seg_count, &intersections);
        eprintln!("Slices ({}):", slices.len());
        for (i, s) in slices.iter().enumerate() {
            let mid_idx = s.vertices.len() / 2;
            let mid = &s.vertices[mid_idx];
            // Min dist to original pline
            let mut min_d = f64::MAX;
            let orig = &pline;
            let n = orig.vertices.len();
            for si in 0..orig.segment_count() {
                let v0 = &orig.vertices[si];
                let v1 = &orig.vertices[(si + 1) % n];
                let d = point_to_segment_dist(mid.x, mid.y, v0.x, v0.y, v1.x, v1.y);
                if d < min_d { min_d = d; }
            }
            eprintln!("  Slice {i}: {} verts, mid=({:.4}, {:.4}), dist_to_orig={:.4}",
                s.vertices.len(), mid.x, mid.y, min_d);
        }

        let valid = filter::apply(&slices, &pline, abs_d);
        eprintln!("Valid slices: {}", valid.len());
        for (i, s) in valid.iter().enumerate() {
            eprintln!("  Valid {i}: {} verts, start=({:.4}, {:.4}), end=({:.4}, {:.4})",
                s.vertices.len(),
                s.vertices.first().map_or(0.0, |v| v.x),
                s.vertices.first().map_or(0.0, |v| v.y),
                s.vertices.last().map_or(0.0, |v| v.x),
                s.vertices.last().map_or(0.0, |v| v.y));
        }

        let result = stitch::connect(&valid);
        eprintln!("Stitched: {} polys", result.len());
        for (i, p) in result.iter().enumerate() {
            eprintln!("  Poly {i}: {} verts", p.vertices.len());
            for (j, v) in p.vertices.iter().enumerate() {
                eprintln!("    [{j}] ({:.4}, {:.4})", v.x, v.y);
            }
        }
    }

    #[test]
    fn mixed_line_arc_square_with_rounded_corner() {
        // Square with one rounded corner (quarter-circle arc).
        let bulge = std::f64::consts::FRAC_PI_8.tan(); // quarter circle
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::new(10.0, 10.0, bulge), // quarter arc corner
                PlineVertex::line(0.0, 10.0),
            ],
            closed: true,
        };
        let op = PlineOffset2D::new(pline, 0.5);
        let result = op.execute().unwrap();
        assert!(!result.is_empty(), "should produce at least one result");
    }
}
