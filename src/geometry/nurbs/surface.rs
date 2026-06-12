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

/// Validates that `value` lies within `[min, max]` (with tolerance),
/// reporting the named `parameter` axis on failure.
fn check_axis(parameter: &'static str, value: f64, min: f64, max: f64) -> Result<()> {
    if value < min - TOLERANCE || value > max + TOLERANCE {
        return Err(GeometryError::ParameterOutOfRange {
            parameter,
            value,
            min,
            max,
        }
        .into());
    }
    Ok(())
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
        check_axis("u", u, u_min, u_max)?;
        check_axis("v", v, v_min, v_max)?;
        Ok(())
    }

    /// Evaluates the surface at parameters `(u, v)` (The NURBS Book, A4.3).
    ///
    /// # Errors
    ///
    /// Returns an error if the parameters are outside the domain or the
    /// rational denominator vanishes.
    // A4.3 single-char bindings follow The NURBS Book notation.
    #[allow(clippy::many_single_char_names)]
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
    // A3.6 / A4.4 single-char bindings and index-driven loops follow The NURBS Book.
    #[allow(clippy::many_single_char_names, clippy::needless_range_loop)]
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

    /// Extracts the isoparametric curve at fixed `u` (a curve in the v
    /// direction) by homogeneously contracting the u-basis.
    ///
    /// # Errors
    ///
    /// Returns an error if `u` is outside the domain, an isocurve weight
    /// vanishes, or the resulting curve fails construction.
    pub fn isocurve_u(&self, u: f64) -> Result<NurbsCurve3D> {
        let (u_min, u_max) = self.knots_u.domain(self.degree_u);
        check_axis("u", u, u_min, u_max)?;
        let span_u = self.knots_u.find_span(self.degree_u, self.nu, u);
        let basis_u = basis_functions(&self.knots_u, span_u, u, self.degree_u);

        let mut points = Vec::with_capacity(self.nv);
        let mut weights = Vec::with_capacity(self.nv);
        for j in 0..self.nv {
            let mut wp = Vector3::zeros();
            let mut w = 0.0;
            for (r, &bu) in basis_u.iter().enumerate() {
                let i = span_u - self.degree_u + r;
                let wij = self.weight(i, j);
                wp += self.control_point(i, j).coords * (bu * wij);
                w += bu * wij;
            }
            if w.abs() < TOLERANCE {
                return Err(GeometryError::Degenerate("zero isocurve weight".into()).into());
            }
            points.push(Point3::from(wp / w));
            weights.push(w);
        }
        NurbsCurve3D::new(points, weights, self.knots_v.clone(), self.degree_v)
    }

    /// Extracts the isoparametric curve at fixed `v` (a curve in the u
    /// direction) by homogeneously contracting the v-basis.
    ///
    /// # Errors
    ///
    /// Returns an error if `v` is outside the domain, an isocurve weight
    /// vanishes, or the resulting curve fails construction.
    pub fn isocurve_v(&self, v: f64) -> Result<NurbsCurve3D> {
        let (v_min, v_max) = self.knots_v.domain(self.degree_v);
        check_axis("v", v, v_min, v_max)?;
        let span_v = self.knots_v.find_span(self.degree_v, self.nv, v);
        let basis_v = basis_functions(&self.knots_v, span_v, v, self.degree_v);

        let mut points = Vec::with_capacity(self.nu);
        let mut weights = Vec::with_capacity(self.nu);
        for i in 0..self.nu {
            let mut wp = Vector3::zeros();
            let mut w = 0.0;
            for (r, &bv) in basis_v.iter().enumerate() {
                let j = span_v - self.degree_v + r;
                let wij = self.weight(i, j);
                wp += self.control_point(i, j).coords * (bv * wij);
                w += bv * wij;
            }
            if w.abs() < TOLERANCE {
                return Err(GeometryError::Degenerate("zero isocurve weight".into()).into());
            }
            points.push(Point3::from(wp / w));
            weights.push(w);
        }
        NurbsCurve3D::new(points, weights, self.knots_u.clone(), self.degree_u)
    }

    /// Extracts the four boundary curves `[u_min edge, u_max edge, v_min edge,
    /// v_max edge]`.
    ///
    /// # Errors
    ///
    /// Returns an error if any isocurve extraction fails.
    pub fn boundary_curves(&self) -> Result<[NurbsCurve3D; 4]> {
        let ((u_min, u_max), (v_min, v_max)) = self.parameter_domain();
        Ok([
            self.isocurve_u(u_min)?,
            self.isocurve_u(u_max)?,
            self.isocurve_v(v_min)?,
            self.isocurve_v(v_max)?,
        ])
    }

    /// Finds the closest point on the surface to `query` via a coarse seed
    /// grid followed by a clamped Newton iteration (The NURBS Book, §6.1).
    ///
    /// Non-convergence is not an error: the best parameters found are returned.
    ///
    /// # Errors
    ///
    /// Returns an error only if evaluating the seed grid or the final point
    /// fails (e.g. a vanishing rational denominator).
    // §6.1 Newton iteration: su/sv/f/g/jNN bindings and the exact grid-index to
    // f64 conversions follow The NURBS Book notation.
    #[allow(
        clippy::many_single_char_names,
        clippy::similar_names,
        clippy::cast_precision_loss
    )]
    pub fn closest_point(
        &self,
        query: &Point3,
        options: &InversionOptions,
    ) -> Result<SurfaceInversion> {
        let ((u_min, u_max), (v_min, v_max)) = self.parameter_domain();
        let samples = options.seed_samples.max(2);

        // Coarse seed: minimize squared distance over the parameter grid.
        let mut best_u = u_min;
        let mut best_v = v_min;
        let mut best_dist_sq = f64::INFINITY;
        for iu in 0..=samples {
            let u = u_min + (u_max - u_min) * (iu as f64) / (samples as f64);
            for iv in 0..=samples {
                let v = v_min + (v_max - v_min) * (iv as f64) / (samples as f64);
                let p = self.point_at(u, v)?;
                let d_sq = (p - query).norm_squared();
                if d_sq < best_dist_sq {
                    best_dist_sq = d_sq;
                    best_u = u;
                    best_v = v;
                }
            }
        }

        let mut u = best_u;
        let mut v = best_v;
        for _ in 0..options.max_iterations {
            let skl = self.derivatives(u, v, 2)?;
            let su = skl[1][0];
            let sv = skl[0][1];
            let r = skl[0][0] - query.coords;
            let r_norm = r.norm();
            if r_norm < options.tolerance {
                break;
            }
            let f = r.dot(&su);
            let g = r.dot(&sv);
            let su_norm = su.norm();
            let sv_norm = sv.norm();
            // Zero cosine (orthogonality) convergence test.
            if su_norm > TOLERANCE && sv_norm > TOLERANCE {
                let cos_u = f.abs() / (su_norm * r_norm);
                let cos_v = g.abs() / (sv_norm * r_norm);
                if cos_u < 1e-10 && cos_v < 1e-10 {
                    break;
                }
            }
            let j00 = su_norm * su_norm + r.dot(&skl[2][0]);
            let j01 = su.dot(&sv) + r.dot(&skl[1][1]);
            let j11 = sv_norm * sv_norm + r.dot(&skl[0][2]);
            let det = j00 * j11 - j01 * j01;
            if det.abs() < TOLERANCE {
                break;
            }
            let du = (-f * j11 + g * j01) / det;
            let dv = (f * j01 - g * j00) / det;
            let new_u = (u + du).clamp(u_min, u_max);
            let new_v = (v + dv).clamp(v_min, v_max);
            let step = (new_u - u).abs() + (new_v - v).abs();
            u = new_u;
            v = new_v;
            if step < options.tolerance {
                break;
            }
        }

        let point = self.point_at(u, v)?;
        let distance = (point - query).norm();
        Ok(SurfaceInversion {
            u,
            v,
            point,
            distance,
        })
    }

    /// Axis-aligned bounding box of the control points.
    ///
    /// By the convex-hull property the surface lies entirely inside this box, so
    /// it is a conservative bound used for intersection-candidate pruning.
    /// Returns `(min_corner, max_corner)`.
    #[must_use]
    pub fn bounding_box(&self) -> (Point3, Point3) {
        let mut min = self.control_points[0].coords;
        let mut max = self.control_points[0].coords;
        for p in &self.control_points[1..] {
            min = min.inf(&p.coords);
            max = max.sup(&p.coords);
        }
        (Point3::from(min), Point3::from(max))
    }

    /// Inserts `u` into the u-knot vector to full multiplicity `degree_u`, then
    /// splits into two surfaces covering `[u_min, u]` and `[u, u_max]` in u
    /// (same v-domain). Mirrors [`NurbsCurve::split`] applied column-wise on the
    /// homogeneous grid.
    ///
    /// # Errors
    ///
    /// Returns an error if `u` is not strictly inside the u-domain interior or a
    /// sub-surface fails construction.
    pub fn split_u(&self, u: f64) -> Result<(Self, Self)> {
        let ((u_min, u_max), _) = self.parameter_domain();
        if u <= u_min + TOLERANCE || u >= u_max - TOLERANCE {
            return Err(GeometryError::ParameterOutOfRange {
                parameter: "u",
                value: u,
                min: u_min,
                max: u_max,
            }
            .into());
        }
        let p = self.degree_u;
        let s = self.knots_u.multiplicity(u);
        let (knots, grid, new_nu) = refine_grid_u(self, u, p - s);

        let first = knots
            .iter()
            .position(|&x| (x - u).abs() < TOLERANCE)
            .ok_or_else(|| GeometryError::Degenerate("split knot not found".into()))?;

        // Left: u-rows [0, first), knots [0..first+p] ++ [u].
        let mut left_knots = knots[..first + p].to_vec();
        left_knots.push(u);
        let left = Self::from_homogeneous_u(
            &grid,
            new_nu,
            self.nv,
            0,
            first,
            KnotVector::new(left_knots)?,
            self.knots_v.clone(),
            p,
            self.degree_v,
        )?;

        // Right: u-rows [first-1, ..), knots [u; p+1] ++ knots[first+p..].
        let mut right_knots = vec![u; p + 1];
        right_knots.extend_from_slice(&knots[first + p..]);
        let right = Self::from_homogeneous_u(
            &grid,
            new_nu,
            self.nv,
            first - 1,
            new_nu,
            KnotVector::new(right_knots)?,
            self.knots_v.clone(),
            p,
            self.degree_v,
        )?;
        Ok((left, right))
    }

    /// Inserts `v` into the v-knot vector to full multiplicity `degree_v`, then
    /// splits into two surfaces covering `[v_min, v]` and `[v, v_max]` in v
    /// (same u-domain). Mirrors [`NurbsCurve::split`] applied row-wise on the
    /// homogeneous grid.
    ///
    /// # Errors
    ///
    /// Returns an error if `v` is not strictly inside the v-domain interior or a
    /// sub-surface fails construction.
    pub fn split_v(&self, v: f64) -> Result<(Self, Self)> {
        let (_, (v_min, v_max)) = self.parameter_domain();
        if v <= v_min + TOLERANCE || v >= v_max - TOLERANCE {
            return Err(GeometryError::ParameterOutOfRange {
                parameter: "v",
                value: v,
                min: v_min,
                max: v_max,
            }
            .into());
        }
        let p = self.degree_v;
        let s = self.knots_v.multiplicity(v);
        let (knots, grid, new_nv) = refine_grid_v(self, v, p - s);

        let first = knots
            .iter()
            .position(|&x| (x - v).abs() < TOLERANCE)
            .ok_or_else(|| GeometryError::Degenerate("split knot not found".into()))?;

        let mut left_knots = knots[..first + p].to_vec();
        left_knots.push(v);
        let left = Self::from_homogeneous_v(
            &grid,
            self.nu,
            new_nv,
            0,
            first,
            self.knots_u.clone(),
            KnotVector::new(left_knots)?,
            self.degree_u,
            p,
        )?;

        let mut right_knots = vec![v; p + 1];
        right_knots.extend_from_slice(&knots[first + p..]);
        let right = Self::from_homogeneous_v(
            &grid,
            self.nu,
            new_nv,
            first - 1,
            new_nv,
            self.knots_u.clone(),
            KnotVector::new(right_knots)?,
            self.degree_u,
            p,
        )?;
        Ok((left, right))
    }

    /// Rebuilds a surface from a homogeneous `(w*P, w)` u-major grid keeping
    /// only u-rows `[u_lo, u_hi)`.
    #[allow(clippy::too_many_arguments)]
    fn from_homogeneous_u(
        grid: &[(Vector3, f64)],
        grid_nu: usize,
        nv: usize,
        u_lo: usize,
        u_hi: usize,
        knots_u: KnotVector,
        knots_v: KnotVector,
        degree_u: usize,
        degree_v: usize,
    ) -> Result<Self> {
        debug_assert!(u_hi <= grid_nu);
        let mut points = Vec::with_capacity((u_hi - u_lo) * nv);
        let mut weights = Vec::with_capacity((u_hi - u_lo) * nv);
        for i in u_lo..u_hi {
            for j in 0..nv {
                let (wp, w) = grid[i * nv + j];
                if w <= 0.0 {
                    return Err(GeometryError::Degenerate(
                        "non-positive homogeneous weight".into(),
                    )
                    .into());
                }
                points.push(Point3::from(wp / w));
                weights.push(w);
            }
        }
        Self::new(
            points,
            weights,
            u_hi - u_lo,
            nv,
            knots_u,
            knots_v,
            degree_u,
            degree_v,
        )
    }

    /// Rebuilds a surface from a homogeneous `(w*P, w)` u-major grid keeping
    /// only v-columns `[v_lo, v_hi)`.
    #[allow(clippy::too_many_arguments)]
    fn from_homogeneous_v(
        grid: &[(Vector3, f64)],
        nu: usize,
        grid_nv: usize,
        v_lo: usize,
        v_hi: usize,
        knots_u: KnotVector,
        knots_v: KnotVector,
        degree_u: usize,
        degree_v: usize,
    ) -> Result<Self> {
        let mut points = Vec::with_capacity(nu * (v_hi - v_lo));
        let mut weights = Vec::with_capacity(nu * (v_hi - v_lo));
        for i in 0..nu {
            for j in v_lo..v_hi {
                let (wp, w) = grid[i * grid_nv + j];
                if w <= 0.0 {
                    return Err(GeometryError::Degenerate(
                        "non-positive homogeneous weight".into(),
                    )
                    .into());
                }
                points.push(Point3::from(wp / w));
                weights.push(w);
            }
        }
        Self::new(
            points,
            weights,
            nu,
            v_hi - v_lo,
            knots_u,
            knots_v,
            degree_u,
            degree_v,
        )
    }
}

/// A refined knot vector, the homogeneous `(w*P, w)` u-major grid it produced,
/// and the new control-point count along the refined axis.
type RefinedGrid = (Vec<f64>, Vec<(Vector3, f64)>, usize);

/// Per-column knot insertion in u: refines every v-column (which shares
/// `knots_u`/`degree_u`) by inserting `u` `times` times on the homogeneous
/// grid. Returns the refined u-knot vector, the new u-major homogeneous grid,
/// and the new u-row count. Mirrors [`NurbsCurve::insert_knot`] index
/// arithmetic, applied column-by-column.
// A5.1-style single-char bindings and index-driven copies follow The NURBS Book.
#[allow(
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    clippy::manual_memcpy
)]
fn refine_grid_u(s: &NurbsSurface, u: f64, times: usize) -> RefinedGrid {
    let p = s.degree_u;
    let np = s.nu;
    let nv = s.nv;
    let r = times;
    let mult = s.knots_u.multiplicity(u);
    let k = s.knots_u.find_span(p, np, u);
    let knots = s.knots_u.as_slice();
    let new_nu = np + r;

    let hw = |i: usize, j: usize| -> (Vector3, f64) {
        let w = s.weight(i, j);
        (s.control_point(i, j).coords * w, w)
    };

    let mut grid = vec![(Vector3::zeros(), 0.0); new_nu * nv];
    if r == 0 {
        for i in 0..np {
            for j in 0..nv {
                grid[i * nv + j] = hw(i, j);
            }
        }
        return (knots.to_vec(), grid, new_nu);
    }

    for j in 0..nv {
        // Front copy: leading unaffected rows.
        for i in 0..=(k - p) {
            grid[i * nv + j] = hw(i, j);
        }
        // Back copy: trailing unaffected rows.
        for i in (k - mult)..np {
            grid[(i + r) * nv + j] = hw(i, j);
        }
        // Affected window temp array.
        let mut work: Vec<(Vector3, f64)> = (0..=(p - mult)).map(|i| hw(k - p + i, j)).collect();
        for jj in 1..=r {
            let l = k - p + jj;
            for i in 0..=(p - jj - mult) {
                let alpha = (u - knots[l + i]) / (knots[i + k + 1] - knots[l + i]);
                work[i] = (
                    work[i + 1].0 * alpha + work[i].0 * (1.0 - alpha),
                    work[i + 1].1 * alpha + work[i].1 * (1.0 - alpha),
                );
            }
            grid[l * nv + j] = work[0];
            grid[(k + r - jj - mult) * nv + j] = work[p - jj - mult];
        }
        let l = k - p + r;
        for i in (l + 1)..(k - mult) {
            grid[i * nv + j] = work[i - l];
        }
    }

    let mut new_knots = Vec::with_capacity(knots.len() + r);
    new_knots.extend_from_slice(&knots[..=k]);
    new_knots.extend(std::iter::repeat_n(u, r));
    new_knots.extend_from_slice(&knots[k + 1..]);
    (new_knots, grid, new_nu)
}

/// Per-row knot insertion in v: refines every u-row (sharing `knots_v`/
/// `degree_v`) by inserting `v` `times` times on the homogeneous grid. Returns
/// the refined v-knot vector, the new u-major homogeneous grid, and the new
/// v-column count.
#[allow(
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    clippy::manual_memcpy
)]
fn refine_grid_v(s: &NurbsSurface, v: f64, times: usize) -> RefinedGrid {
    let p = s.degree_v;
    let np = s.nv;
    let nu = s.nu;
    let r = times;
    let mult = s.knots_v.multiplicity(v);
    let k = s.knots_v.find_span(p, np, v);
    let knots = s.knots_v.as_slice();
    let new_nv = np + r;

    let hw = |i: usize, j: usize| -> (Vector3, f64) {
        let w = s.weight(i, j);
        (s.control_point(i, j).coords * w, w)
    };

    let mut grid = vec![(Vector3::zeros(), 0.0); nu * new_nv];
    if r == 0 {
        for i in 0..nu {
            for j in 0..np {
                grid[i * new_nv + j] = hw(i, j);
            }
        }
        return (knots.to_vec(), grid, new_nv);
    }

    for i in 0..nu {
        for j in 0..=(k - p) {
            grid[i * new_nv + j] = hw(i, j);
        }
        for j in (k - mult)..np {
            grid[i * new_nv + (j + r)] = hw(i, j);
        }
        let mut work: Vec<(Vector3, f64)> = (0..=(p - mult)).map(|jj| hw(i, k - p + jj)).collect();
        for jj in 1..=r {
            let l = k - p + jj;
            for c in 0..=(p - jj - mult) {
                let alpha = (v - knots[l + c]) / (knots[c + k + 1] - knots[l + c]);
                work[c] = (
                    work[c + 1].0 * alpha + work[c].0 * (1.0 - alpha),
                    work[c + 1].1 * alpha + work[c].1 * (1.0 - alpha),
                );
            }
            grid[i * new_nv + l] = work[0];
            grid[i * new_nv + (k + r - jj - mult)] = work[p - jj - mult];
        }
        let l = k - p + r;
        for c in (l + 1)..(k - mult) {
            grid[i * new_nv + c] = work[c - l];
        }
    }

    let mut new_knots = Vec::with_capacity(knots.len() + r);
    new_knots.extend_from_slice(&knots[..=k]);
    new_knots.extend(std::iter::repeat_n(v, r));
    new_knots.extend_from_slice(&knots[k + 1..]);
    (new_knots, grid, new_nv)
}

/// Options controlling Newton point inversion on a surface.
#[derive(Debug, Clone, Copy)]
pub struct InversionOptions {
    /// Maximum Newton iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on the Euclidean residual and parameter step.
    pub tolerance: f64,
    /// Seed-grid samples per parametric direction. For oscillatory or
    /// high-curvature surfaces, raise this to avoid seeding the wrong basin.
    pub seed_samples: usize,
}

impl Default for InversionOptions {
    fn default() -> Self {
        Self {
            max_iterations: 64,
            tolerance: 1e-12,
            seed_samples: 16,
        }
    }
}

/// Result of a closest-point (inversion) query on a NURBS surface.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceInversion {
    /// U parameter of the closest point.
    pub u: f64,
    /// V parameter of the closest point.
    pub v: f64,
    /// The closest point on the surface.
    pub point: Point3,
    /// Distance from the query point.
    pub distance: f64,
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
// Tests reuse The NURBS Book single-char / su,sv-style notation.
#[allow(
    clippy::unwrap_used,
    clippy::many_single_char_names,
    clippy::similar_names
)]
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

    /// Quarter-cylinder shell: exact rational quadratic quarter circle in u
    /// (XY plane, radius 1, weights [1, 1/sqrt(2), 1]) extruded linearly in
    /// v along +Z by 2. A genuinely rational surface: interior weight != 1.
    fn quarter_cylinder_patch() -> NurbsSurface {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsSurface::new(
            vec![
                // u-major, nv = 2: (i, j) = i * 2 + j
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 2.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 2.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 2.0),
            ],
            vec![1.0, 1.0, w, w, 1.0, 1.0],
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

    #[test]
    fn isocurve_u_matches_surface() {
        let s = parabolic_patch();
        let c = s.isocurve_u(0.3).unwrap();
        for i in 0..=20 {
            let v = f64::from(i) / 20.0;
            assert!((c.point_at(v).unwrap() - s.point_at(0.3, v).unwrap()).norm() < 1e-12);
        }
    }

    #[test]
    fn isocurve_v_matches_surface() {
        let s = parabolic_patch();
        let c = s.isocurve_v(0.7).unwrap();
        for i in 0..=20 {
            let u = f64::from(i) / 20.0;
            assert!((c.point_at(u).unwrap() - s.point_at(u, 0.7).unwrap()).norm() < 1e-12);
        }
    }

    #[test]
    fn boundary_curves_trace_patch_edges() {
        let s = bilinear_patch();
        let [u0, u1, v0, v1] = s.boundary_curves().unwrap();
        assert!((u0.point_at(0.5).unwrap() - s.point_at(0.0, 0.5).unwrap()).norm() < 1e-12);
        assert!((u1.point_at(0.5).unwrap() - s.point_at(1.0, 0.5).unwrap()).norm() < 1e-12);
        assert!((v0.point_at(0.5).unwrap() - s.point_at(0.5, 0.0).unwrap()).norm() < 1e-12);
        assert!((v1.point_at(0.5).unwrap() - s.point_at(0.5, 1.0).unwrap()).norm() < 1e-12);
    }

    #[test]
    fn inversion_round_trips_surface_points() {
        let s = parabolic_patch();
        let options = InversionOptions::default();
        for &(u, v) in &[(0.2, 0.3), (0.5, 0.5), (0.85, 0.1)] {
            let p = s.point_at(u, v).unwrap();
            let result = s.closest_point(&p, &options).unwrap();
            assert!(
                result.distance < 1e-9,
                "distance {} at ({u},{v})",
                result.distance
            );
            let q = s.point_at(result.u, result.v).unwrap();
            assert!((q - p).norm() < 1e-9);
        }
    }

    #[test]
    fn inversion_projects_off_surface_point() {
        // For the planar bilinear patch, projection of an elevated point is
        // straight down and parameters equal the (scaled) XY coordinates.
        let s = bilinear_patch();
        let result = s
            .closest_point(&Point3::new(1.0, 1.0, 5.0), &InversionOptions::default())
            .unwrap();
        assert!((result.distance - 5.0).abs() < 1e-9);
        assert!((result.u - 0.5).abs() < 1e-9);
        assert!((result.v - 0.5).abs() < 1e-9);
    }

    #[test]
    fn inversion_clamps_to_domain_for_outside_points() {
        let s = bilinear_patch();
        let result = s
            .closest_point(&Point3::new(5.0, 1.0, 0.0), &InversionOptions::default())
            .unwrap();
        // Closest point is on the u_max edge
        assert!((result.u - 1.0).abs() < 1e-9);
        let p = s.point_at(result.u, result.v).unwrap();
        assert!((p - Point3::new(2.0, 1.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn rejects_wrong_knot_count() {
        let result = NurbsSurface::from_unweighted(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(2.0, 2.0, 0.0),
            ],
            2,
            2,
            KnotVector::new(vec![0.0, 0.0, 0.5, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
            1,
        );
        assert!(result.is_err());
    }

    #[test]
    fn rejects_non_positive_weight() {
        let result = NurbsSurface::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(2.0, 2.0, 0.0),
            ],
            vec![1.0, 1.0, 0.0, 1.0],
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
    fn rational_patch_lies_on_cylinder_and_derivatives_match_differences() {
        let s = quarter_cylinder_patch();
        let h = 1e-6;
        for &(u, v) in &[(0.2, 0.3), (0.5, 0.5), (0.8, 0.7)] {
            // Every point must lie on the unit cylinder x^2 + y^2 = 1
            let p = s.point_at(u, v).unwrap();
            let radial = (p.x * p.x + p.y * p.y).sqrt();
            assert!((radial - 1.0).abs() < 1e-12, "off cylinder at ({u},{v})");

            // First partials vs central differences (rational correction active)
            let d = s.derivatives(u, v, 2).unwrap();
            let du_fd = (s.point_at(u + h, v).unwrap() - s.point_at(u - h, v).unwrap()) / (2.0 * h);
            let dv_fd = (s.point_at(u, v + h).unwrap() - s.point_at(u, v - h).unwrap()) / (2.0 * h);
            assert!((d[1][0] - du_fd).norm() < 1e-5, "du at ({u},{v})");
            assert!((d[0][1] - dv_fd).norm() < 1e-5, "dv at ({u},{v})");

            // Second u-partial vs central differences
            let h2 = 1e-4;
            let duu_fd = (s.point_at(u + h2, v).unwrap().coords
                - 2.0 * s.point_at(u, v).unwrap().coords
                + s.point_at(u - h2, v).unwrap().coords)
                / (h2 * h2);
            assert!((d[2][0] - duu_fd).norm() < 1e-3, "duu at ({u},{v})");
        }
    }

    #[test]
    fn bounding_box_contains_sampled_points() {
        let s = quarter_cylinder_patch();
        let (min, max) = s.bounding_box();
        for iu in 0..=20 {
            let u = f64::from(iu) / 20.0;
            for iv in 0..=20 {
                let v = f64::from(iv) / 20.0;
                let p = s.point_at(u, v).unwrap();
                for k in 0..3 {
                    assert!(
                        p.coords[k] >= min.coords[k] - 1e-12
                            && p.coords[k] <= max.coords[k] + 1e-12,
                        "point ({u},{v}) axis {k} escapes box"
                    );
                }
            }
        }
    }

    /// Both sub-surfaces must agree with the original over their sub-domains.
    fn assert_split_u_matches(s: &NurbsSurface, u: f64) {
        let (left, right) = s.split_u(u).unwrap();
        let ((u_min, u_max), (v_min, v_max)) = s.parameter_domain();
        for iu in 0..=50 {
            let uu = u_min + (u - u_min) * f64::from(iu) / 50.0;
            for iv in 0..=10 {
                let vv = v_min + (v_max - v_min) * f64::from(iv) / 10.0;
                assert!(
                    (left.point_at(uu, vv).unwrap() - s.point_at(uu, vv).unwrap()).norm() < 1e-12,
                    "left mismatch at ({uu},{vv})"
                );
            }
        }
        for iu in 0..=50 {
            let uu = u + (u_max - u) * f64::from(iu) / 50.0;
            for iv in 0..=10 {
                let vv = v_min + (v_max - v_min) * f64::from(iv) / 10.0;
                assert!(
                    (right.point_at(uu, vv).unwrap() - s.point_at(uu, vv).unwrap()).norm() < 1e-12,
                    "right mismatch at ({uu},{vv})"
                );
            }
        }
    }

    fn assert_split_v_matches(s: &NurbsSurface, v: f64) {
        let (left, right) = s.split_v(v).unwrap();
        let ((u_min, u_max), (v_min, v_max)) = s.parameter_domain();
        for iv in 0..=50 {
            let vv = v_min + (v - v_min) * f64::from(iv) / 50.0;
            for iu in 0..=10 {
                let uu = u_min + (u_max - u_min) * f64::from(iu) / 10.0;
                assert!(
                    (left.point_at(uu, vv).unwrap() - s.point_at(uu, vv).unwrap()).norm() < 1e-12,
                    "left mismatch at ({uu},{vv})"
                );
            }
        }
        for iv in 0..=50 {
            let vv = v + (v_max - v) * f64::from(iv) / 50.0;
            for iu in 0..=10 {
                let uu = u_min + (u_max - u_min) * f64::from(iu) / 10.0;
                assert!(
                    (right.point_at(uu, vv).unwrap() - s.point_at(uu, vv).unwrap()).norm() < 1e-12,
                    "right mismatch at ({uu},{vv})"
                );
            }
        }
    }

    #[test]
    fn split_u_halves_match_bilinear() {
        assert_split_u_matches(&bilinear_patch(), 0.4);
    }

    #[test]
    fn split_v_halves_match_bilinear() {
        assert_split_v_matches(&bilinear_patch(), 0.6);
    }

    #[test]
    fn split_u_halves_match_rational_cylinder() {
        assert_split_u_matches(&quarter_cylinder_patch(), 0.5);
        assert_split_u_matches(&quarter_cylinder_patch(), 0.3);
    }

    #[test]
    fn split_v_halves_match_parabolic() {
        assert_split_v_matches(&parabolic_patch(), 0.5);
    }

    #[test]
    fn split_u_rejects_boundary() {
        assert!(bilinear_patch().split_u(0.0).is_err());
        assert!(bilinear_patch().split_u(1.0).is_err());
    }

    #[test]
    fn rational_patch_inversion_round_trips() {
        let s = quarter_cylinder_patch();
        let options = InversionOptions::default();
        for &(u, v) in &[(0.25, 0.4), (0.6, 0.8)] {
            let p = s.point_at(u, v).unwrap();
            let result = s.closest_point(&p, &options).unwrap();
            assert!(
                result.distance < 1e-9,
                "distance {} at ({u},{v})",
                result.distance
            );
        }
    }
}
