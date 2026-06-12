use nalgebra::{Point, SVector};

use crate::error::{GeometryError, Result};
use crate::math::TOLERANCE;

use super::basis::basis_functions;
use super::knot::KnotVector;

/// A NURBS curve in `D`-dimensional space.
///
/// Control points and weights are stored separately; weights must be
/// strictly positive. Evaluation uses homogeneous accumulation
/// (The NURBS Book, A4.1).
#[derive(Debug, Clone)]
pub struct NurbsCurve<const D: usize> {
    control_points: Vec<Point<f64, D>>,
    weights: Vec<f64>,
    knots: KnotVector,
    degree: usize,
}

/// A NURBS curve in 2D (used for surface-parameter-space curves).
pub type NurbsCurve2D = NurbsCurve<2>;

/// A NURBS curve in 3D space.
pub type NurbsCurve3D = NurbsCurve<3>;

impl<const D: usize> NurbsCurve<D> {
    /// Creates a NURBS curve, validating the structural invariants so that
    /// every later internal call to [`KnotVector::find_span`] / [`basis_functions`]
    /// is consistent by construction.
    ///
    /// # Errors
    ///
    /// Returns an error if `degree < 1`, there are fewer than `degree + 1`
    /// control points, the weight count differs from the control-point count,
    /// any weight is not strictly positive, or the knot count is not
    /// `control_points.len() + degree + 1`.
    pub fn new(
        control_points: Vec<Point<f64, D>>,
        weights: Vec<f64>,
        knots: KnotVector,
        degree: usize,
    ) -> Result<Self> {
        if degree < 1 {
            return Err(GeometryError::Degenerate("curve degree must be >= 1".into()).into());
        }
        if control_points.len() < degree + 1 {
            return Err(GeometryError::Degenerate(format!(
                "need at least {} control points for degree {degree}, got {}",
                degree + 1,
                control_points.len()
            ))
            .into());
        }
        if weights.len() != control_points.len() {
            return Err(GeometryError::Degenerate(format!(
                "weight count {} does not match control-point count {}",
                weights.len(),
                control_points.len()
            ))
            .into());
        }
        if weights.iter().any(|&w| w < TOLERANCE) {
            return Err(
                GeometryError::Degenerate("weights must be strictly positive".into()).into(),
            );
        }
        let expected_knots = control_points.len() + degree + 1;
        if knots.len() != expected_knots {
            return Err(GeometryError::Degenerate(format!(
                "expected {expected_knots} knots, got {}",
                knots.len()
            ))
            .into());
        }
        Ok(Self {
            control_points,
            weights,
            knots,
            degree,
        })
    }

    /// Creates a NURBS curve with all weights set to `1.0` (a B-spline curve).
    ///
    /// # Errors
    ///
    /// Same validation as [`NurbsCurve::new`].
    pub fn from_unweighted(
        control_points: Vec<Point<f64, D>>,
        knots: KnotVector,
        degree: usize,
    ) -> Result<Self> {
        let weights = vec![1.0; control_points.len()];
        Self::new(control_points, weights, knots, degree)
    }

    /// Control points of the curve.
    #[must_use]
    pub fn control_points(&self) -> &[Point<f64, D>] {
        &self.control_points
    }

    /// Weights of the control points.
    #[must_use]
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// Knot vector of the curve.
    #[must_use]
    pub fn knots(&self) -> &KnotVector {
        &self.knots
    }

    /// Degree of the curve.
    #[must_use]
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Parameter domain `[t_min, t_max]` of the curve.
    #[must_use]
    pub fn parameter_domain(&self) -> (f64, f64) {
        self.knots.domain(self.degree)
    }

    /// Validates that `t` lies within the parameter domain (with tolerance).
    fn validate_parameter(&self, t: f64) -> Result<()> {
        let (min, max) = self.parameter_domain();
        if t < min - TOLERANCE || t > max + TOLERANCE {
            return Err(GeometryError::ParameterOutOfRange {
                parameter: "t",
                value: t,
                min,
                max,
            }
            .into());
        }
        Ok(())
    }

    /// Evaluates the curve at parameter `t` (The NURBS Book, A4.1).
    ///
    /// # Errors
    ///
    /// Returns an error if `t` is outside the parameter domain or the rational
    /// denominator vanishes.
    pub fn point_at(&self, t: f64) -> Result<Point<f64, D>> {
        self.validate_parameter(t)?;
        let span = self
            .knots
            .find_span(self.degree, self.control_points.len(), t);
        let basis = basis_functions(&self.knots, span, t, self.degree);

        let mut numerator = SVector::<f64, D>::zeros();
        let mut denominator = 0.0;
        for (j, &b) in basis.iter().enumerate() {
            let idx = span - self.degree + j;
            let w = b * self.weights[idx];
            numerator += self.control_points[idx].coords * w;
            denominator += w;
        }
        if denominator.abs() < TOLERANCE {
            return Err(GeometryError::Degenerate("zero rational denominator".into()).into());
        }
        Ok(Point::from(numerator / denominator))
    }

    /// Evaluates derivatives up to `order` at `t` (The NURBS Book, A4.2).
    ///
    /// Returns `d` where `d[0]` is the position vector and `d[k]` is the
    /// k-th derivative.
    ///
    /// # Errors
    ///
    /// Returns an error if `t` is outside the parameter domain.
    pub fn derivatives(&self, t: f64, order: usize) -> Result<Vec<SVector<f64, D>>> {
        use super::basis::{basis_function_derivatives, binomial};

        self.validate_parameter(t)?;
        let span = self
            .knots
            .find_span(self.degree, self.control_points.len(), t);
        let nders = basis_function_derivatives(&self.knots, span, t, self.degree, order);

        // Homogeneous derivatives: a[k] = sum N^{(k)} * w_i * P_i, w[k] = sum N^{(k)} * w_i
        let mut a = vec![SVector::<f64, D>::zeros(); order + 1];
        let mut w = vec![0.0; order + 1];
        for (k, row) in nders.iter().enumerate().take(order + 1) {
            for (j, nd) in row.iter().enumerate() {
                let idx = span - self.degree + j;
                a[k] += self.control_points[idx].coords * (nd * self.weights[idx]);
                w[k] += nd * self.weights[idx];
            }
        }
        if w[0].abs() < TOLERANCE {
            return Err(GeometryError::Degenerate("zero rational denominator".into()).into());
        }

        // Rational derivatives: C^{(k)} = (A^{(k)} - sum_{i=1..k} C(k,i) w^{(i)} C^{(k-i)}) / w
        let mut ders = vec![SVector::<f64, D>::zeros(); order + 1];
        for k in 0..=order {
            let mut v = a[k];
            for i in 1..=k {
                v -= ders[k - i] * (binomial(k, i) * w[i]);
            }
            ders[k] = v / w[0];
        }
        Ok(ders)
    }

    /// Whether the curve endpoints coincide.
    #[must_use]
    pub fn is_endpoint_closed(&self) -> bool {
        let (t_min, t_max) = self.parameter_domain();
        match (self.point_at(t_min), self.point_at(t_max)) {
            (Ok(a), Ok(b)) => (a - b).norm() < TOLERANCE,
            _ => false,
        }
    }
}

use crate::geometry::curve::{Curve, CurveDomain};
use crate::math::{Point3, Vector3};

impl Curve for NurbsCurve3D {
    fn evaluate(&self, t: f64) -> Result<Point3> {
        self.point_at(t)
    }

    fn tangent(&self, t: f64) -> Result<Vector3> {
        let ders = self.derivatives(t, 1)?;
        let len = ders[1].norm();
        if len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        Ok(ders[1] / len)
    }

    fn domain(&self) -> CurveDomain {
        let (t_min, t_max) = self.parameter_domain();
        CurveDomain::new(t_min, t_max)
    }

    fn is_closed(&self) -> bool {
        self.is_endpoint_closed()
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::{Point3, TOLERANCE};
    use std::f64::consts::FRAC_1_SQRT_2;

    fn line_curve() -> NurbsCurve3D {
        // Degree-1 curve between two points = straight line
        NurbsCurve3D::from_unweighted(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 4.0, 6.0)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap()
    }

    fn quarter_circle() -> NurbsCurve3D {
        // Exact rational quadratic quarter circle in the XY plane, radius 1
        NurbsCurve3D::new(
            vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![1.0, FRAC_1_SQRT_2, 1.0],
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap(),
            2,
        )
        .unwrap()
    }

    #[test]
    fn rejects_count_mismatch() {
        let result = NurbsCurve3D::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![1.0],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        );
        assert!(result.is_err());
    }

    #[test]
    fn rejects_wrong_knot_count() {
        let result = NurbsCurve3D::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![1.0, 1.0],
            KnotVector::new(vec![0.0, 0.0, 0.5, 1.0, 1.0]).unwrap(),
            1,
        );
        assert!(result.is_err());
    }

    #[test]
    fn rejects_non_positive_weight() {
        let result = NurbsCurve3D::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![1.0, 0.0],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        );
        assert!(result.is_err());
    }

    #[test]
    fn degree_one_is_linear_interpolation() {
        let c = line_curve();
        let p = c.point_at(0.25).unwrap();
        assert!((p - Point3::new(0.5, 1.0, 1.5)).norm() < TOLERANCE);
    }

    #[test]
    fn rejects_parameter_outside_domain() {
        let c = line_curve();
        assert!(c.point_at(1.5).is_err());
        assert!(c.point_at(-0.1).is_err());
    }

    #[test]
    fn quarter_circle_midpoint_is_exact() {
        let c = quarter_circle();
        let p = c.point_at(0.5).unwrap();
        assert!((p - Point3::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn quarter_circle_stays_on_unit_circle() {
        let c = quarter_circle();
        for i in 0..=20 {
            let t = f64::from(i) / 20.0;
            let p = c.point_at(t).unwrap();
            assert!((p.coords.norm() - 1.0).abs() < 1e-12, "t={t}");
        }
    }

    #[test]
    fn nurbs_curve_2d_evaluates() {
        use crate::math::Point2;
        let c = NurbsCurve2D::from_unweighted(
            vec![Point2::new(0.0, 0.0), Point2::new(2.0, 2.0)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap();
        let p = c.point_at(0.5).unwrap();
        assert!((p - Point2::new(1.0, 1.0)).norm() < TOLERANCE);
    }

    use crate::geometry::curve::Curve;
    use crate::math::Vector3;

    #[test]
    fn first_derivative_matches_central_difference() {
        let c = quarter_circle();
        let h = 1e-6;
        for i in 1..10 {
            let t = f64::from(i) / 10.0;
            let d = c.derivatives(t, 1).unwrap();
            let fd = (c.point_at(t + h).unwrap() - c.point_at(t - h).unwrap()) / (2.0 * h);
            assert!((d[1] - fd).norm() < 1e-5, "t={t}");
        }
    }

    #[test]
    fn second_derivative_matches_central_difference() {
        let c = quarter_circle();
        let h = 1e-4;
        let t = 0.5;
        let d = c.derivatives(t, 2).unwrap();
        let fd = (c.point_at(t + h).unwrap().coords - 2.0 * c.point_at(t).unwrap().coords
            + c.point_at(t - h).unwrap().coords)
            / (h * h);
        assert!((d[2] - fd).norm() < 1e-3);
    }

    #[test]
    fn derivative_order_zero_is_point() {
        let c = quarter_circle();
        let d = c.derivatives(0.3, 0).unwrap();
        let p = c.point_at(0.3).unwrap();
        assert!((d[0] - p.coords).norm() < TOLERANCE);
    }

    #[test]
    fn curve_trait_evaluate_and_tangent() {
        let c = quarter_circle();
        let p = Curve::evaluate(&c, 0.0).unwrap();
        assert!((p - Point3::new(1.0, 0.0, 0.0)).norm() < TOLERANCE);
        // Tangent at start of the quarter circle points in +Y
        let t = c.tangent(0.0).unwrap();
        assert!((t - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-9);
        assert!(
            (t.norm() - 1.0).abs() < 1e-12,
            "tangent must be unit length"
        );
    }

    #[test]
    fn curve_trait_domain_and_closed() {
        let c = quarter_circle();
        let d = c.domain();
        assert!((d.t_min - 0.0).abs() < TOLERANCE);
        assert!((d.t_max - 1.0).abs() < TOLERANCE);
        assert!(!c.is_closed());
    }
}
