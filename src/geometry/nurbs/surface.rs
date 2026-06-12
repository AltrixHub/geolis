use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use super::basis::{basis_function_derivatives, basis_functions, binomial};
use super::curve::NurbsCurve3D;
use super::knot::KnotVector;

/// A NURBS surface in 3D space.
///
/// Control points form an `nu x nv` grid stored u-major: the point at
/// u-index `i` and v-index `j` is `control_points[i * nv + j]`.
#[derive(Debug, Clone)]
pub struct NurbsSurface {
    control_points: Vec<Point3>,
    weights: Vec<f64>,
    nu: usize,
    nv: usize,
    knots_u: KnotVector,
    knots_v: KnotVector,
    degree_u: usize,
    degree_v: usize,
}

impl NurbsSurface {
    /// Creates a NURBS surface, validating the structural invariants so that
    /// every later internal call to [`KnotVector::find_span`] / [`basis_functions`]
    /// is consistent by construction.
    ///
    /// # Errors
    ///
    /// Returns an error if either degree is `< 1`, the grid is too small for
    /// the degrees, the control-point count does not equal `nu * nv`, the
    /// weight count differs, any weight is not strictly positive and finite, or
    /// either knot count is not `n + degree + 1`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        control_points: Vec<Point3>,
        weights: Vec<f64>,
        nu: usize,
        nv: usize,
        knots_u: KnotVector,
        knots_v: KnotVector,
        degree_u: usize,
        degree_v: usize,
    ) -> Result<Self> {
        if degree_u < 1 || degree_v < 1 {
            return Err(GeometryError::Degenerate("surface degrees must be >= 1".into()).into());
        }
        if nu < degree_u + 1 || nv < degree_v + 1 {
            return Err(GeometryError::Degenerate(format!(
                "grid {nu}x{nv} too small for degrees {degree_u},{degree_v} \
                 (need at least {}x{})",
                degree_u + 1,
                degree_v + 1
            ))
            .into());
        }
        if control_points.len() != nu * nv {
            return Err(GeometryError::Degenerate(format!(
                "control-point count {} does not match grid {nu}x{nv} = {}",
                control_points.len(),
                nu * nv
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
        if weights.iter().any(|&w| w <= 0.0 || !w.is_finite()) {
            return Err(GeometryError::Degenerate(
                "weights must be strictly positive and finite".into(),
            )
            .into());
        }
        let expected_u = nu + degree_u + 1;
        if knots_u.len() != expected_u {
            return Err(GeometryError::Degenerate(format!(
                "expected {expected_u} u-knots, got {}",
                knots_u.len()
            ))
            .into());
        }
        let expected_v = nv + degree_v + 1;
        if knots_v.len() != expected_v {
            return Err(GeometryError::Degenerate(format!(
                "expected {expected_v} v-knots, got {}",
                knots_v.len()
            ))
            .into());
        }
        Ok(Self {
            control_points,
            weights,
            nu,
            nv,
            knots_u,
            knots_v,
            degree_u,
            degree_v,
        })
    }

    /// Creates a NURBS surface with all weights set to `1.0` (a B-spline
    /// surface).
    ///
    /// # Errors
    ///
    /// Same validation as [`NurbsSurface::new`].
    #[allow(clippy::too_many_arguments)]
    pub fn from_unweighted(
        control_points: Vec<Point3>,
        nu: usize,
        nv: usize,
        knots_u: KnotVector,
        knots_v: KnotVector,
        degree_u: usize,
        degree_v: usize,
    ) -> Result<Self> {
        let weights = vec![1.0; control_points.len()];
        Self::new(
            control_points,
            weights,
            nu,
            nv,
            knots_u,
            knots_v,
            degree_u,
            degree_v,
        )
    }

    /// Control point at u-index `i`, v-index `j`.
    #[must_use]
    pub fn control_point(&self, i: usize, j: usize) -> &Point3 {
        &self.control_points[i * self.nv + j]
    }

    /// Weight at u-index `i`, v-index `j`.
    #[must_use]
    pub fn weight(&self, i: usize, j: usize) -> f64 {
        self.weights[i * self.nv + j]
    }

    /// Grid size `(nu, nv)`.
    #[must_use]
    pub fn grid_size(&self) -> (usize, usize) {
        (self.nu, self.nv)
    }

    /// Degrees `(degree_u, degree_v)`.
    #[must_use]
    pub fn degrees(&self) -> (usize, usize) {
        (self.degree_u, self.degree_v)
    }

    /// Knot vector in the u direction.
    #[must_use]
    pub fn knots_u(&self) -> &KnotVector {
        &self.knots_u
    }

    /// Knot vector in the v direction.
    #[must_use]
    pub fn knots_v(&self) -> &KnotVector {
        &self.knots_v
    }

    /// Parameter domain `((u_min, u_max), (v_min, v_max))`.
    #[must_use]
    pub fn parameter_domain(&self) -> ((f64, f64), (f64, f64)) {
        (
            self.knots_u.domain(self.degree_u),
            self.knots_v.domain(self.degree_v),
        )
    }

    /// Validates that `(u, v)` lies within the parameter domain (with
    /// tolerance), reporting `u` then `v`.
    fn validate_parameters(&self, u: f64, v: f64) -> Result<()> {
        let ((u_min, u_max), (v_min, v_max)) = self.parameter_domain();
        if u < u_min - TOLERANCE || u > u_max + TOLERANCE {
            return Err(GeometryError::ParameterOutOfRange {
                parameter: "u",
                value: u,
                min: u_min,
                max: u_max,
            }
            .into());
        }
        if v < v_min - TOLERANCE || v > v_max + TOLERANCE {
            return Err(GeometryError::ParameterOutOfRange {
                parameter: "v",
                value: v,
                min: v_min,
                max: v_max,
            }
            .into());
        }
        Ok(())
    }

    /// Evaluates the surface at parameters `(u, v)` (The NURBS Book, A4.3).
    ///
    /// # Errors
    ///
    /// Returns an error if the parameters are outside the domain or the
    /// rational denominator vanishes.
    pub fn point_at(&self, u: f64, v: f64) -> Result<Point3> {
        self.validate_parameters(u, v)?;
        let span_u = self.knots_u.find_span(self.degree_u, self.nu, u);
        let span_v = self.knots_v.find_span(self.degree_v, self.nv, v);
        let basis_u = basis_functions(&self.knots_u, span_u, u, self.degree_u);
        let basis_v = basis_functions(&self.knots_v, span_v, v, self.degree_v);

        let mut numerator = Vector3::zeros();
        let mut denominator = 0.0;
        for (ru, &bu) in basis_u.iter().enumerate() {
            let i = span_u - self.degree_u + ru;
            for (rv, &bv) in basis_v.iter().enumerate() {
                let j = span_v - self.degree_v + rv;
                let w = bu * bv * self.weight(i, j);
                numerator += self.control_point(i, j).coords * w;
                denominator += w;
            }
        }
        if denominator.abs() < TOLERANCE {
            return Err(GeometryError::Degenerate("zero rational denominator".into()).into());
        }
        Ok(Point3::from(numerator / denominator))
    }

    /// Computes partial derivatives of the surface up to total order `order`
    /// at `(u, v)` (The NURBS Book, A3.6 for the homogeneous derivatives,
    /// A4.4 / eq 4.20 for the rational correction).
    ///
    /// Returns `skl` where `skl[k][l]` = ∂^{k+l}S/∂u^k∂v^l and `skl[0][0]` is
    /// the position vector.
    ///
    /// # Errors
    ///
    /// Returns an error if the parameters are outside the domain or the
    /// rational denominator vanishes.
    pub fn derivatives(&self, u: f64, v: f64, order: usize) -> Result<Vec<Vec<Vector3>>> {
        self.validate_parameters(u, v)?;
        let span_u = self.knots_u.find_span(self.degree_u, self.nu, u);
        let span_v = self.knots_v.find_span(self.degree_v, self.nv, v);
        let nders_u = basis_function_derivatives(&self.knots_u, span_u, u, self.degree_u, order);
        let nders_v = basis_function_derivatives(&self.knots_v, span_v, v, self.degree_v, order);

        // Homogeneous derivatives: aders[k][l] for the weighted points,
        // wders[k][l] for the weights.
        let mut aders = vec![vec![Vector3::zeros(); order + 1]; order + 1];
        let mut wders = vec![vec![0.0_f64; order + 1]; order + 1];
        let du = order.min(self.degree_u);
        let dv = order.min(self.degree_v);
        for k in 0..=du {
            // temp[s] = sum_r nders_u[k][r] * (w_ij * P_ij, w_ij)
            let mut temp = vec![(Vector3::zeros(), 0.0_f64); self.degree_v + 1];
            for s in 0..=self.degree_v {
                let j = span_v - self.degree_v + s;
                let mut acc_p = Vector3::zeros();
                let mut acc_w = 0.0;
                for r in 0..=self.degree_u {
                    let i = span_u - self.degree_u + r;
                    let nd = nders_u[k][r];
                    let w = self.weight(i, j);
                    acc_p += self.control_point(i, j).coords * (nd * w);
                    acc_w += nd * w;
                }
                temp[s] = (acc_p, acc_w);
            }
            let l_max = dv.min(order - k);
            for l in 0..=l_max {
                let mut acc_p = Vector3::zeros();
                let mut acc_w = 0.0;
                for s in 0..=self.degree_v {
                    let nd = nders_v[l][s];
                    acc_p += temp[s].0 * nd;
                    acc_w += temp[s].1 * nd;
                }
                aders[k][l] = acc_p;
                wders[k][l] = acc_w;
            }
        }

        if wders[0][0].abs() < TOLERANCE {
            return Err(GeometryError::Degenerate("zero rational denominator".into()).into());
        }

        // Rational correction (A4.4 / eq 4.20).
        let mut skl = vec![vec![Vector3::zeros(); order + 1]; order + 1];
        for k in 0..=order {
            for l in 0..=(order - k) {
                let mut value = aders[k][l];
                for j in 1..=l {
                    value -= skl[k][l - j] * (binomial(l, j) * wders[0][j]);
                }
                for i in 1..=k {
                    value -= skl[k - i][l] * (binomial(k, i) * wders[i][0]);
                    let mut inner = Vector3::zeros();
                    for j in 1..=l {
                        inner += skl[k - i][l - j] * (binomial(l, j) * wders[i][j]);
                    }
                    value -= inner * binomial(k, i);
                }
                skl[k][l] = value / wders[0][0];
            }
        }
        Ok(skl)
    }

    /// Convenience accessor returning `(point, dS/du, dS/dv)`.
    ///
    /// # Errors
    ///
    /// Returns an error if the parameters are outside the domain or the
    /// rational denominator vanishes.
    pub fn partials(&self, u: f64, v: f64) -> Result<(Point3, Vector3, Vector3)> {
        let skl = self.derivatives(u, v, 1)?;
        Ok((Point3::from(skl[0][0]), skl[1][0], skl[0][1]))
    }
}

use crate::geometry::surface::{Surface, SurfaceDomain};

impl Surface for NurbsSurface {
    fn evaluate(&self, u: f64, v: f64) -> Result<Point3> {
        self.point_at(u, v)
    }

    fn normal(&self, u: f64, v: f64) -> Result<Vector3> {
        let (_, su, sv) = self.partials(u, v)?;
        let n = su.cross(&sv);
        let len = n.norm();
        if len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        Ok(n / len)
    }

    fn domain(&self) -> SurfaceDomain {
        let ((u_min, u_max), (v_min, v_max)) = self.parameter_domain();
        SurfaceDomain::new(u_min, u_max, v_min, v_max)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::surface::Surface;
    use crate::math::{Point3, Vector3, TOLERANCE};

    /// 2x2 bilinear patch spanning [0,2]x[0,2] in the XY plane.
    fn bilinear_patch() -> NurbsSurface {
        NurbsSurface::from_unweighted(
            vec![
                // u-major: index = i * nv + j (i: u index, j: v index)
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(2.0, 2.0, 0.0),
            ],
            2,
            2,
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
            1,
        )
        .unwrap()
    }

    /// Quadratic-in-u patch with a z-lift (exact polynomial surface for
    /// derivative checks).
    fn parabolic_patch() -> NurbsSurface {
        NurbsSurface::from_unweighted(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(2.0, 0.0, 2.0),
                Point3::new(2.0, 1.0, 2.0),
            ],
            3,
            2,
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            2,
            1,
        )
        .unwrap()
    }

    #[test]
    fn rejects_grid_count_mismatch() {
        let result = NurbsSurface::from_unweighted(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(2.0, 2.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            2,
            2,
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
            1,
        );
        assert!(result.is_err());
    }

    #[test]
    fn bilinear_patch_interpolates_corners_and_center() {
        let s = bilinear_patch();
        assert!((s.point_at(0.0, 0.0).unwrap() - Point3::new(0.0, 0.0, 0.0)).norm() < TOLERANCE);
        assert!((s.point_at(1.0, 1.0).unwrap() - Point3::new(2.0, 2.0, 0.0)).norm() < TOLERANCE);
        assert!((s.point_at(0.5, 0.5).unwrap() - Point3::new(1.0, 1.0, 0.0)).norm() < TOLERANCE);
    }

    #[test]
    fn rejects_out_of_domain_parameters() {
        let s = bilinear_patch();
        assert!(s.point_at(1.5, 0.5).is_err());
        assert!(s.point_at(0.5, -0.5).is_err());
    }

    #[test]
    fn parabolic_patch_midpoint() {
        let s = parabolic_patch();
        let p = s.point_at(0.5, 0.5).unwrap();
        assert!((p - Point3::new(1.0, 0.5, 0.5)).norm() < 1e-12);
    }

    #[test]
    fn partials_match_central_differences() {
        let s = parabolic_patch();
        let h = 1e-6;
        for &(u, v) in &[(0.3, 0.4), (0.5, 0.5), (0.7, 0.2)] {
            let d = s.derivatives(u, v, 1).unwrap();
            let du_fd = (s.point_at(u + h, v).unwrap() - s.point_at(u - h, v).unwrap()) / (2.0 * h);
            let dv_fd = (s.point_at(u, v + h).unwrap() - s.point_at(u, v - h).unwrap()) / (2.0 * h);
            assert!((d[1][0] - du_fd).norm() < 1e-5, "du at ({u},{v})");
            assert!((d[0][1] - dv_fd).norm() < 1e-5, "dv at ({u},{v})");
        }
    }

    #[test]
    fn second_partials_match_central_differences() {
        let s = parabolic_patch();
        let h = 1e-4;
        let (u, v) = (0.5, 0.5);
        let d = s.derivatives(u, v, 2).unwrap();
        let duu_fd = (s.point_at(u + h, v).unwrap().coords
            - 2.0 * s.point_at(u, v).unwrap().coords
            + s.point_at(u - h, v).unwrap().coords)
            / (h * h);
        assert!((d[2][0] - duu_fd).norm() < 1e-3);
    }

    #[test]
    fn derivative_order_zero_is_point() {
        let s = parabolic_patch();
        let d = s.derivatives(0.3, 0.6, 0).unwrap();
        let p = s.point_at(0.3, 0.6).unwrap();
        assert!((d[0][0] - p.coords).norm() < TOLERANCE);
    }

    #[test]
    fn surface_trait_normal_of_planar_patch_is_z() {
        let s = bilinear_patch();
        let n = Surface::normal(&s, 0.5, 0.5).unwrap();
        assert!((n - Vector3::new(0.0, 0.0, 1.0)).norm() < TOLERANCE);
        assert!((n.norm() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn surface_trait_domain() {
        let s = bilinear_patch();
        let d = Surface::domain(&s);
        assert!((d.u_min - 0.0).abs() < TOLERANCE);
        assert!((d.u_max - 1.0).abs() < TOLERANCE);
        assert!((d.v_min - 0.0).abs() < TOLERANCE);
        assert!((d.v_max - 1.0).abs() < TOLERANCE);
    }
}
