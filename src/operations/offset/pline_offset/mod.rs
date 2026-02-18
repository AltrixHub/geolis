mod filter;
mod raw_offset;
mod self_intersect;
mod slice;
mod stitch;

use crate::error::{OperationError, Result};
use crate::geometry::pline::Pline;

/// Offsets a polyline (with potential arc segments) using the slice-and-filter
/// algorithm.
///
/// For closed polylines: positive distance = inward, negative = outward.
/// For open polylines: positive distance = left side, negative = right side.
/// Returns offset curve(s) without endpoint caps.
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
        let result = stitch::connect(&valid, true);

        if result.is_empty() {
            return Err(OperationError::Failed(
                "offset collapsed completely".to_owned(),
            )
            .into());
        }

        Ok(result)
    }

    /// Executes offset for open polylines using the slice-and-filter pipeline.
    ///
    /// Positive distance offsets to the left (when facing along the polyline
    /// direction), negative distance offsets to the right.  Returns open
    /// polyline(s) without endpoint caps.
    fn execute_open(&self) -> Result<Vec<Pline>> {
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
        let result = stitch::connect(&valid, false);

        if result.is_empty() {
            return Err(OperationError::Failed(
                "offset collapsed completely".to_owned(),
            )
            .into());
        }

        Ok(result)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::pline::PlineVertex;

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
        assert!(!result[0].closed, "open pline offset should stay open");
    }

    // ── Open polyline single-sided offset tests ──

    #[test]
    fn open_l_shape_left_offset() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(5.0, 5.0),
            ],
            closed: false,
        };
        let op = PlineOffset2D::new(pline, 0.5); // Left offset
        let result = op.execute().unwrap();
        assert_eq!(result.len(), 1);
        assert!(!result[0].closed, "result should be open");
        assert_eq!(result[0].vertices.len(), 3);
    }

    #[test]
    fn open_l_shape_right_offset() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(5.0, 5.0),
            ],
            closed: false,
        };
        let op = PlineOffset2D::new(pline, -0.5); // Right offset
        let result = op.execute().unwrap();
        assert_eq!(result.len(), 1);
        assert!(!result[0].closed, "result should be open");
    }

    #[test]
    fn open_arc_segment_offset() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::new(5.0, 0.0, 1.0), // CCW semicircle
                PlineVertex::line(10.0, 0.0),
            ],
            closed: false,
        };
        let op = PlineOffset2D::new(pline, 0.5);
        let result = op.execute().unwrap();
        assert_eq!(result.len(), 1);
        assert!(!result[0].closed);
        let has_arcs = result[0].vertices.iter().any(|v| v.bulge.abs() > 1e-6);
        assert!(has_arcs, "offset should preserve arc segments");
    }

    #[test]
    fn open_straight_line_offset() {
        // Simplest case: a single straight segment.
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
            ],
            closed: false,
        };
        let op = PlineOffset2D::new(pline, 1.0); // Left = +Y
        let result = op.execute().unwrap();
        assert_eq!(result.len(), 1);
        assert!(!result[0].closed, "result should be open");
        assert_eq!(result[0].vertices.len(), 2);
        // Offset left by 1 → y should be ~1.0
        for v in &result[0].vertices {
            assert!((v.y - 1.0).abs() < 0.01, "y should be ~1.0, got {}", v.y);
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
