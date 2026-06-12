use nalgebra::{Point, SVector};

use crate::error::{GeometryError, Result};
use crate::math::TOLERANCE;

use super::basis::{basis_functions, binomial};
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
    /// any weight is not strictly positive and finite, or the knot count is not
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
        if weights.iter().any(|&w| w <= 0.0 || !w.is_finite()) {
            return Err(GeometryError::Degenerate(
                "weights must be strictly positive and finite".into(),
            )
            .into());
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
    // TODO(perf): hot-path callers (intersection marching) may want an allocation-free variant.
    pub fn derivatives(&self, t: f64, order: usize) -> Result<Vec<SVector<f64, D>>> {
        use super::basis::basis_function_derivatives;

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

    /// Returns a new curve with `u` inserted `times` times into the knot
    /// vector without changing the curve shape (The NURBS Book, A5.1).
    ///
    /// # Errors
    ///
    /// Returns an error if `u` is outside the domain or the resulting knot
    /// multiplicity would exceed the degree.
    // A5.1 single-char bindings and index-driven copies follow The NURBS Book.
    #[allow(
        clippy::many_single_char_names,
        clippy::needless_range_loop,
        clippy::manual_memcpy
    )]
    pub fn insert_knot(&self, u: f64, times: usize) -> Result<Self> {
        self.validate_parameter(u)?;
        let p = self.degree;
        let s = self.knots.multiplicity(u);
        if times == 0 {
            return Ok(self.clone());
        }
        if s + times > p {
            return Err(GeometryError::Degenerate(format!(
                "inserting {times} knots at u={u} would exceed degree {p} (multiplicity {s})"
            ))
            .into());
        }
        let np = self.control_points.len();
        let k = self.knots.find_span(p, np, u);
        let knots = self.knots.as_slice();
        let r = times;

        // Homogeneous control points (w * P, w)
        let hw = |i: usize| -> (SVector<f64, D>, f64) {
            (
                self.control_points[i].coords * self.weights[i],
                self.weights[i],
            )
        };

        let mut new_hp = vec![(SVector::<f64, D>::zeros(), 0.0); np + r];
        // Front copy: unaffected leading points new_hp[0..=k-p].
        for i in 0..=(k - p) {
            new_hp[i] = hw(i);
        }
        // Back copy: unaffected trailing points new_hp[i+r] for i in (k-s)..np.
        for i in (k - s)..np {
            new_hp[i + r] = hw(i);
        }

        // work: the book's Rw temp array over the affected p-s+1 points.
        let mut work: Vec<(SVector<f64, D>, f64)> = (0..=(p - s)).map(|i| hw(k - p + i)).collect();
        for j in 1..=r {
            let l = k - p + j;
            // Per-insertion alpha blend over the affected window.
            for i in 0..=(p - j - s) {
                let alpha = (u - knots[l + i]) / (knots[i + k + 1] - knots[l + i]);
                work[i] = (
                    work[i + 1].0 * alpha + work[i].0 * (1.0 - alpha),
                    work[i + 1].1 * alpha + work[i].1 * (1.0 - alpha),
                );
            }
            // Two boundary writes bracketing this insertion.
            new_hp[l] = work[0];
            new_hp[k + r - j - s] = work[p - j - s];
        }
        // Remaining-middle copy of the still-untouched work entries.
        let l = k - p + r;
        for i in (l + 1)..(k - s) {
            new_hp[i] = work[i - l];
        }

        // New knot vector: knots[0..=k] ++ [u; r] ++ knots[k+1..]
        let mut new_knots = Vec::with_capacity(knots.len() + r);
        new_knots.extend_from_slice(&knots[..=k]);
        new_knots.extend(std::iter::repeat_n(u, r));
        new_knots.extend_from_slice(&knots[k + 1..]);

        let (points, weights) = dehomogenize(&new_hp)?;
        Self::new(points, weights, KnotVector::new(new_knots)?, p)
    }

    /// Splits the curve at `u` into two curves covering `[t_min, u]` and
    /// `[u, t_max]`.
    ///
    /// # Errors
    ///
    /// Returns an error if `u` is outside the open domain interior.
    pub fn split(&self, u: f64) -> Result<(Self, Self)> {
        let (t_min, t_max) = self.parameter_domain();
        if u <= t_min + TOLERANCE || u >= t_max - TOLERANCE {
            return Err(GeometryError::ParameterOutOfRange {
                parameter: "u",
                value: u,
                min: t_min,
                max: t_max,
            }
            .into());
        }
        let p = self.degree;
        let s = self.knots.multiplicity(u);
        let refined = self.insert_knot(u, p - s)?;

        let knots = refined.knots.as_slice();
        // Index of the first knot equal to u (multiplicity is now exactly p).
        let first = knots
            .iter()
            .position(|&x| (x - u).abs() < TOLERANCE)
            .ok_or_else(|| GeometryError::Degenerate("split knot not found".into()))?;

        // Left: control points [0, first), knots [0..first+p] ++ [u]
        let mut left_knots = knots[..first + p].to_vec();
        left_knots.push(u);
        let left = Self::new(
            refined.control_points[..first].to_vec(),
            refined.weights[..first].to_vec(),
            KnotVector::new(left_knots)?,
            p,
        )?;

        // Right: control points [first-1, ..], knots [u; p+1] ++ knots[first+p..]
        // The domain-interior check guarantees `u` is an interior knot after refinement,
        // so first >= p + 1 >= 2 and `first - 1` cannot underflow.
        let mut right_knots = vec![u; p + 1];
        right_knots.extend_from_slice(&knots[first + p..]);
        let right = Self::new(
            refined.control_points[first - 1..].to_vec(),
            refined.weights[first - 1..].to_vec(),
            KnotVector::new(right_knots)?,
            p,
        )?;

        Ok((left, right))
    }

    /// Returns the curve with reversed parameterization (same shape).
    ///
    /// # Errors
    ///
    /// Returns an error if reconstructing the knot vector or curve fails. This
    /// path is unreachable for a valid curve (the reflected-reversed knots are
    /// non-decreasing and the structural invariants are preserved), but the
    /// error is propagated rather than masked by a silent fallback.
    pub fn reverse(&self) -> Result<Self> {
        let knots = self.knots.as_slice();
        let (a, b) = self.parameter_domain();
        // With `sum` computed once, `sum - k` is monotone under IEEE subtraction:
        // reversing a non-decreasing sequence yields exactly non-decreasing knots,
        // so no drift-clamp is needed.
        let sum = a + b;
        let new_knots: Vec<f64> = knots.iter().rev().map(|&k| sum - k).collect();
        let control_points: Vec<_> = self.control_points.iter().rev().copied().collect();
        let weights: Vec<_> = self.weights.iter().rev().copied().collect();
        Self::new(
            control_points,
            weights,
            KnotVector::new(new_knots)?,
            self.degree,
        )
    }

    /// Raises the curve degree by `t` without changing its shape
    /// (The NURBS Book, A5.9).
    ///
    /// Degree elevation operates on the homogeneous control points
    /// `(w * P, w)` so that rational curves (non-unit weights) are elevated
    /// exactly, then dehomogenizes the result.
    ///
    /// # Errors
    ///
    /// Returns an error if reconstructing the elevated knot vector or curve
    /// fails. This path is unreachable for a valid curve; the error is
    /// propagated rather than masked by a silent fallback. `t == 0` returns a
    /// clone.
    // A5.9: single-char bindings (a, b, ph, mh, kind, oldr, ...) and the
    // index-driven Bezier-segment loops follow The NURBS Book notation.
    #[allow(
        clippy::many_single_char_names,
        clippy::needless_range_loop,
        clippy::too_many_lines,
        clippy::similar_names
    )]
    pub fn elevate_degree(&self, t: usize) -> Result<Self> {
        if t == 0 {
            return Ok(self.clone());
        }
        let p = self.degree;
        let ph = p + t;
        let ph2 = ph / 2;

        // Homogeneous control points (w * P, w).
        type H<const D: usize> = (SVector<f64, D>, f64);
        let hw = |i: usize| -> H<D> {
            (
                self.control_points[i].coords * self.weights[i],
                self.weights[i],
            )
        };
        let zero_h = || -> H<D> { (SVector::<f64, D>::zeros(), 0.0) };
        let scale = |h: H<D>, s: f64| -> H<D> { (h.0 * s, h.1 * s) };
        let add = |a: H<D>, b: H<D>| -> H<D> { (a.0 + b.0, a.1 + b.1) };

        let u = self.knots.as_slice();
        let m = u.len() - 1;
        let n = self.control_points.len() - 1;

        // Bezout coefficients for degree-elevation of one Bezier segment
        // (eq 5.36): bezalfs[i][j].
        let mut bezalfs = vec![vec![0.0_f64; p + 1]; ph + 1];
        bezalfs[0][0] = 1.0;
        bezalfs[ph][p] = 1.0;
        for i in 1..=ph2 {
            let inv = 1.0 / binomial(ph, i);
            let mpi = p.min(i);
            for j in (i.saturating_sub(t))..=mpi {
                bezalfs[i][j] = inv * binomial(p, j) * binomial(t, i - j);
            }
        }
        for i in (ph2 + 1)..ph {
            let mpi = p.min(i);
            for j in (i.saturating_sub(t))..=mpi {
                bezalfs[i][j] = bezalfs[ph - i][p - j];
            }
        }

        // Output buffers (over-allocated; truncated to exact length at the end).
        let cap_q = (n + 1) * (t + 1) + ph + 1;
        let mut qw: Vec<H<D>> = vec![zero_h(); cap_q];
        let mut uh: Vec<f64> = vec![0.0; cap_q + ph + 1];

        let mut bpts: Vec<H<D>> = vec![zero_h(); p + 1];
        let mut ebpts: Vec<H<D>> = vec![zero_h(); ph + 1];
        let mut next_bpts: Vec<H<D>> = vec![zero_h(); p.saturating_sub(1)];
        let mut alfs: Vec<f64> = vec![0.0; p.saturating_sub(1).max(1)];

        let mut mh = ph;
        let mut kind = ph + 1;
        let mut r: isize = -1;
        let mut a = p;
        let mut b = p + 1;
        let mut cind = 1;
        let mut ua = u[0];

        qw[0] = hw(0);
        for i in 0..=ph {
            uh[i] = ua;
        }

        // Initialise the first Bezier segment.
        for i in 0..=p {
            bpts[i] = hw(i);
        }

        while b < m {
            let i_start = b;
            while b < m && (u[b] - u[b + 1]).abs() < TOLERANCE {
                b += 1;
            }
            let mul = b - i_start + 1;
            mh += mul + t;
            let ub = u[b];
            let oldr = r;
            r = p as isize - mul as isize;

            // Insert knot u[b] r times to make the next segment.
            let lbz = if oldr > 0 { (oldr + 2) as usize / 2 } else { 1 };
            let rbz = if r > 0 { ph - (r as usize + 1) / 2 } else { ph };

            if r > 0 {
                let r_us = r as usize;
                let numer = ub - ua;
                for k in ((mul + 1)..=p).rev() {
                    alfs[k - mul - 1] = numer / (u[a + k] - ua);
                }
                for j in 1..=r_us {
                    let save = r_us - j;
                    let s = mul + j;
                    for k in (s..=p).rev() {
                        let af = alfs[k - s];
                        bpts[k] = add(scale(bpts[k], af), scale(bpts[k - 1], 1.0 - af));
                    }
                    next_bpts[save] = bpts[p];
                }
            }

            // Degree-elevate the Bezier segment: bpts -> ebpts.
            for i in lbz..=ph {
                ebpts[i] = zero_h();
                let mpi = p.min(i);
                for j in (i.saturating_sub(t))..=mpi {
                    ebpts[i] = add(ebpts[i], scale(bpts[j], bezalfs[i][j]));
                }
            }

            // Remove knot u[a] oldr times (knot-removal blend of the new seg).
            if oldr > 1 {
                let oldr_us = oldr as usize;
                let mut first = kind - 2;
                let mut last = kind;
                let den = ub - ua;
                let mut bet = (ub - uh[kind - 1]) / den;
                for tr in 1..oldr_us {
                    let mut i = first;
                    let mut j = last;
                    let mut kj = j - kind + 1;
                    while j as isize - i as isize > tr as isize {
                        if i < cind {
                            let alf = (ub - uh[i]) / (ua - uh[i]);
                            qw[i] = add(scale(qw[i], alf), scale(qw[i - 1], 1.0 - alf));
                        }
                        if j >= lbz {
                            if (j as isize - tr as isize)
                                <= (kind as isize - ph as isize + last as isize)
                            {
                                let gam = (ub - uh[j - tr]) / den;
                                ebpts[kj] =
                                    add(scale(ebpts[kj], gam), scale(ebpts[kj + 1], 1.0 - gam));
                            } else {
                                ebpts[kj] =
                                    add(scale(ebpts[kj], bet), scale(ebpts[kj + 1], 1.0 - bet));
                            }
                        }
                        i += 1;
                        j -= 1;
                        kj = kj.wrapping_sub(1);
                    }
                    first = first.wrapping_sub(1);
                    last += 1;
                    bet = (ub - uh[kind - 1]) / den;
                }
            }

            // Load the knot ua.
            if a != p {
                for _ in 0..(ph as isize - oldr) {
                    uh[kind] = ua;
                    kind += 1;
                }
            }
            // Load control points into qw.
            for j in lbz..=rbz {
                qw[cind] = ebpts[j];
                cind += 1;
            }

            // Set up for the next pass through the loop.
            if b < m {
                for j in 0..(r as usize) {
                    bpts[j] = next_bpts[j];
                }
                for j in (r as usize)..=p {
                    bpts[j] = hw(b - p + j);
                }
                a = b;
                b += 1;
                ua = ub;
            } else {
                for i in 0..=ph {
                    uh[kind + i] = ub;
                }
            }
        }

        let nh = mh - ph - 1;
        let new_hp = &qw[..=nh];
        let new_knots = uh[..(nh + ph + 2)].to_vec();

        let (points, weights) = dehomogenize(new_hp)?;
        Self::new(points, weights, KnotVector::new(new_knots)?, ph)
    }
}

/// Converts homogeneous (w*P, w) pairs back to points and weights.
fn dehomogenize<const D: usize>(
    hp: &[(SVector<f64, D>, f64)],
) -> Result<(Vec<Point<f64, D>>, Vec<f64>)> {
    let mut points = Vec::with_capacity(hp.len());
    let mut weights = Vec::with_capacity(hp.len());
    for (wp, w) in hp {
        // Convex combination of strictly positive weights is strictly positive,
        // so this is a defensive guard rather than a tolerance question.
        if *w <= 0.0 {
            return Err(GeometryError::Degenerate("non-positive homogeneous weight".into()).into());
        }
        points.push(Point::from(wp / *w));
        weights.push(*w);
    }
    Ok((points, weights))
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
    fn tiny_weight_constructs_and_evaluates_finite() {
        let c = NurbsCurve3D::new(
            vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![1.0, 1e-12, 1.0],
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap(),
            2,
        )
        .unwrap();
        for i in 0..=10 {
            let t = f64::from(i) / 10.0;
            let p = c.point_at(t).unwrap();
            assert!(p.coords.iter().all(|c| c.is_finite()), "t={t}");
        }
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

    /// Sample both curves at shared parameters and assert max deviation.
    fn assert_same_shape(a: &NurbsCurve3D, b: &NurbsCurve3D, tol: f64) {
        let (t0, t1) = a.parameter_domain();
        for i in 0..=50 {
            let t = t0 + (t1 - t0) * f64::from(i) / 50.0;
            let pa = a.point_at(t).unwrap();
            let pb = b.point_at(t).unwrap();
            assert!((pa - pb).norm() < tol, "deviation at t={t}");
        }
    }

    #[test]
    fn knot_insertion_preserves_shape() {
        let c = quarter_circle();
        let refined = c.insert_knot(0.3, 1).unwrap();
        assert_eq!(refined.control_points().len(), 4);
        assert_eq!(refined.knots().len(), 7);
        assert_same_shape(&c, &refined, 1e-12);
    }

    #[test]
    fn knot_insertion_to_full_multiplicity_preserves_shape() {
        let c = quarter_circle();
        let refined = c.insert_knot(0.5, 2).unwrap();
        assert_eq!(refined.knots().multiplicity(0.5), 2);
        assert_same_shape(&c, &refined, 1e-12);
    }

    #[test]
    fn knot_insertion_beyond_degree_fails() {
        let c = quarter_circle();
        assert!(c.insert_knot(0.5, 3).is_err());
    }

    #[test]
    fn split_produces_matching_halves() {
        let c = quarter_circle();
        let (left, right) = c.split(0.4).unwrap();
        for i in 0..=20 {
            let t = 0.4 * f64::from(i) / 20.0;
            assert!((left.point_at(t).unwrap() - c.point_at(t).unwrap()).norm() < 1e-12);
        }
        for i in 0..=20 {
            let t = 0.4 + (1.0 - 0.4) * f64::from(i) / 20.0;
            assert!((right.point_at(t).unwrap() - c.point_at(t).unwrap()).norm() < 1e-12);
        }
    }

    #[test]
    fn reverse_traverses_backwards() {
        let c = quarter_circle();
        let r = c.reverse().unwrap();
        for i in 0..=20 {
            let t = f64::from(i) / 20.0;
            assert!((r.point_at(1.0 - t).unwrap() - c.point_at(t).unwrap()).norm() < 1e-12);
        }
    }

    #[test]
    fn elevate_degree_zero_is_clone() {
        let c = quarter_circle();
        let e = c.elevate_degree(0).unwrap();
        assert_eq!(e.degree(), c.degree());
        assert_eq!(e.control_points().len(), c.control_points().len());
        assert_same_shape(&c, &e, 1e-12);
    }

    #[test]
    fn elevate_line_to_quadratic_preserves_shape() {
        let c = line_curve();
        let e = c.elevate_degree(1).unwrap();
        assert_eq!(e.degree(), 2);
        // n + t rule: a single Bezier segment of degree 1 has 2 control points;
        // raising by 1 yields 3 control points and 6 knots.
        assert_eq!(e.control_points().len(), 3);
        assert_eq!(e.knots().len(), 6);
        assert_same_shape(&c, &e, 1e-12);
    }

    #[test]
    fn elevate_quarter_circle_by_one_preserves_shape_and_weights() {
        let c = quarter_circle();
        let e = c.elevate_degree(1).unwrap();
        assert_eq!(e.degree(), 3);
        assert!(e.weights().iter().all(|&w| w > 0.0));
        assert_same_shape(&c, &e, 1e-12);
    }

    #[test]
    fn elevate_quarter_circle_by_two_preserves_shape_and_weights() {
        let c = quarter_circle();
        let e = c.elevate_degree(2).unwrap();
        assert_eq!(e.degree(), 4);
        assert!(e.weights().iter().all(|&w| w > 0.0));
        assert_same_shape(&c, &e, 1e-12);
    }
}
