use std::f64::consts::{FRAC_PI_2, PI};

use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use crate::geometry::nurbs::NurbsCurve3D;

impl NurbsCurve3D {
    /// Builds a closed rounded-rectangle profile in the plane spanned by
    /// `(x_axis, y_axis)` through `center`.
    ///
    /// The outline is four straight edges joined by four exact quarter-circle
    /// corner arcs, concatenated C0 into a single degree-2 curve (the straight
    /// edges are degree-elevated to match the rational-quadratic arcs). `width`
    /// and `height` are the OUTER dimensions along `x_axis` and `y_axis`;
    /// `corner_radius` must satisfy `2 * corner_radius < min(width, height)` so
    /// every straight edge keeps a positive length.
    ///
    /// Both axes are normalized internally and must be perpendicular. The result
    /// is exactly planar and its corner arcs are exact circular quadrants.
    ///
    /// # Errors
    ///
    /// Returns an error if `width`, `height`, or `corner_radius` is non-positive,
    /// `2 * corner_radius >= min(width, height)`, either axis is zero-length or
    /// they are not perpendicular, or curve construction / concatenation fails.
    pub fn rounded_rectangle(
        center: Point3,
        x_axis: Vector3,
        y_axis: Vector3,
        width: f64,
        height: f64,
        corner_radius: f64,
    ) -> Result<Self> {
        if width <= TOLERANCE || height <= TOLERANCE || corner_radius <= TOLERANCE {
            return Err(GeometryError::Degenerate(
                "rounded rectangle width, height and corner radius must be positive".into(),
            )
            .into());
        }
        if 2.0 * corner_radius >= width.min(height) {
            return Err(GeometryError::Degenerate(
                "rounded rectangle corner radius too large (need 2r < min(width, height))".into(),
            )
            .into());
        }

        let x_len = x_axis.norm();
        let y_len = y_axis.norm();
        if x_len < TOLERANCE || y_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let x_axis = x_axis / x_len;
        let y_axis = y_axis / y_len;
        if x_axis.dot(&y_axis).abs() > TOLERANCE {
            return Err(GeometryError::Degenerate(
                "rounded rectangle axes must be perpendicular".into(),
            )
            .into());
        }
        let normal = x_axis.cross(&y_axis);

        let hw = 0.5 * width;
        let hh = 0.5 * height;
        let r = corner_radius;
        // Inset corner-arc centers (local plane coordinates).
        let ix = hw - r;
        let iy = hh - r;

        // Map local plane coordinates to a 3D point.
        let p = |lx: f64, ly: f64| -> Point3 { center + x_axis * lx + y_axis * ly };
        // A degree-1 straight edge between two local points.
        let edge = |ax: f64, ay: f64, bx: f64, by: f64| -> Result<Self> {
            NurbsCurve3D::polyline(&[p(ax, ay), p(bx, by)])
        };
        // A quarter arc about an inset corner center at local (cx, cy).
        let corner = |cx: f64, cy: f64, start: f64, end: f64| -> Result<Self> {
            NurbsCurve3D::arc(p(cx, cy), r, normal, x_axis, start, end)
        };

        // Walk the outline counter-clockwise starting at the bottom of the
        // right edge: edge, corner, edge, corner, ... (8 segments).
        let segments = [
            edge(hw, -iy, hw, iy)?,                     // right edge (upward)
            corner(ix, iy, 0.0, FRAC_PI_2)?,            // top-right arc
            edge(ix, hh, -ix, hh)?,                     // top edge (leftward)
            corner(-ix, iy, FRAC_PI_2, PI)?,            // top-left arc
            edge(-hw, iy, -hw, -iy)?,                   // left edge (downward)
            corner(-ix, -iy, PI, PI + FRAC_PI_2)?,      // bottom-left arc
            edge(-ix, -hh, ix, -hh)?,                   // bottom edge (rightward)
            corner(ix, -iy, PI + FRAC_PI_2, 2.0 * PI)?, // bottom-right arc
        ];

        Self::concatenate(&segments)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// Signed distance from a local plane point to the rounded-rectangle
    /// boundary (negative inside, zero on the profile). Standard rounded-box SDF.
    fn rounded_box_sdf(lx: f64, ly: f64, hw: f64, hh: f64, r: f64) -> f64 {
        let qx = lx.abs() - (hw - r);
        let qy = ly.abs() - (hh - r);
        let outside = (qx.max(0.0).powi(2) + qy.max(0.0).powi(2)).sqrt();
        outside + qx.max(qy).min(0.0) - r
    }

    #[test]
    fn rounded_rectangle_is_closed() {
        let c = NurbsCurve3D::rounded_rectangle(
            Point3::origin(),
            Vector3::x(),
            Vector3::y(),
            2.6,
            2.0,
            0.35,
        )
        .unwrap();
        assert!(c.is_endpoint_closed(), "rounded rectangle must close");
        assert_eq!(c.degree(), 2, "concatenated to degree 2");
    }

    #[test]
    fn sampled_points_lie_on_the_profile() {
        // Tilted plane to exercise the axis mapping (x + z diagonal, y = +Y).
        let center = Point3::new(1.0, 2.0, 3.0);
        let x_axis = Vector3::new(1.0, 0.0, 1.0);
        let y_axis = Vector3::new(0.0, 1.0, 0.0);
        let (width, height, r) = (2.6, 2.0, 0.35);
        let c = NurbsCurve3D::rounded_rectangle(center, x_axis, y_axis, width, height, r).unwrap();

        let xn = x_axis.normalize();
        let yn = y_axis.normalize();
        let (t0, t1) = c.parameter_domain();
        for i in 0..=200 {
            let t = t0 + (t1 - t0) * f64::from(i) / 200.0;
            let pt = c.point_at(t).unwrap();
            let d = pt - center;
            // Project onto the local plane axes.
            let lx = d.dot(&xn);
            let ly = d.dot(&yn);
            // In-plane (no component along the normal).
            let normal = xn.cross(&yn);
            assert!(d.dot(&normal).abs() < 1e-9, "point off plane at t={t}");
            let sdf = rounded_box_sdf(lx, ly, 0.5 * width, 0.5 * height, r);
            assert!(sdf.abs() < 1e-9, "point off profile at t={t}: sdf={sdf}");
        }
    }

    #[test]
    fn rejects_oversized_radius() {
        // 2r = 2.0 == min(width, height) = 2.0 → rejected.
        let result = NurbsCurve3D::rounded_rectangle(
            Point3::origin(),
            Vector3::x(),
            Vector3::y(),
            2.6,
            2.0,
            1.0,
        );
        assert!(result.is_err(), "2r >= min(width, height) must be rejected");
    }

    #[test]
    fn rejects_non_perpendicular_axes() {
        let result = NurbsCurve3D::rounded_rectangle(
            Point3::origin(),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            2.0,
            2.0,
            0.3,
        );
        assert!(result.is_err(), "non-perpendicular axes must be rejected");
    }
}
