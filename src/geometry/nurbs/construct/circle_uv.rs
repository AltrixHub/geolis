use std::f64::consts::FRAC_1_SQRT_2;

use crate::error::{GeometryError, Result};
use crate::math::{Point2, TOLERANCE};

use crate::geometry::nurbs::{KnotVector, NurbsCurve2D};

impl NurbsCurve2D {
    /// Builds an exact rational quadratic full circle in 2D (UV) space using the
    /// canonical nine-control-point construction (The NURBS Book, A7.1 over four
    /// quarter arcs).
    ///
    /// The circle is centred at `center` with the given `radius`, traversed
    /// counter-clockwise starting from `(center.x + radius, center.y)`.
    ///
    /// # Errors
    ///
    /// Returns [`GeometryError::Degenerate`] if `radius` is not positive.
    pub fn circle_uv(center: Point2, radius: f64) -> Result<NurbsCurve2D> {
        if radius < TOLERANCE {
            return Err(GeometryError::Degenerate("circle radius must be positive".into()).into());
        }

        let r = radius;
        let (cx, cy) = (center.x, center.y);
        // Nine control points around the square circumscribing the circle:
        // axis points (weight 1) alternate with corner points (weight 1/sqrt2).
        let control = vec![
            Point2::new(cx + r, cy),     // 0 deg
            Point2::new(cx + r, cy + r), // corner
            Point2::new(cx, cy + r),     // 90 deg
            Point2::new(cx - r, cy + r), // corner
            Point2::new(cx - r, cy),     // 180 deg
            Point2::new(cx - r, cy - r), // corner
            Point2::new(cx, cy - r),     // 270 deg
            Point2::new(cx + r, cy - r), // corner
            Point2::new(cx + r, cy),     // back to 0 deg
        ];
        let w = FRAC_1_SQRT_2;
        let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
        let knots = KnotVector::new(vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ])?;
        NurbsCurve2D::new(control, weights, knots, 2)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn circle_uv_lies_on_circle() {
        let c = NurbsCurve2D::circle_uv(Point2::new(1.0, 2.0), 3.0).unwrap();
        let (t0, t1) = c.parameter_domain();
        for i in 0..=64 {
            let t = t0 + (t1 - t0) * f64::from(i) / 64.0;
            let p = c.point_at(t).unwrap();
            let r = ((p.x - 1.0).powi(2) + (p.y - 2.0).powi(2)).sqrt();
            assert!((r - 3.0).abs() < 1e-9, "radius off at t={t}: r={r}");
        }
    }

    #[test]
    fn rejects_zero_radius() {
        assert!(NurbsCurve2D::circle_uv(Point2::new(0.0, 0.0), 0.0).is_err());
    }
}
