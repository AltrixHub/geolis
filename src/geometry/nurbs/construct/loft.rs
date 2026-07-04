use nalgebra::{DMatrix, DVector};

use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use crate::geometry::nurbs::{basis_functions, KnotVector, NurbsCurve3D, NurbsSurface};

impl NurbsSurface {
    /// Lofts (skins) a NURBS surface through a sequence of section curves
    /// (The NURBS Book, §10.3).
    ///
    /// The sections are made compatible internally (degree elevation to the
    /// common maximum degree, then a knot-vector union by knot insertion) so
    /// that they share a u direction. The v direction interpolates the
    /// resulting control net: for each u control index the homogeneous control
    /// points of the sections are interpolated by a degree-`degree_v` global
    /// interpolation, so every input section is reproduced exactly at its v
    /// parameter. `degree_v` defaults to `min(3, sections - 1)`.
    ///
    /// # Errors
    ///
    /// Returns an error if fewer than 2 sections are given, the requested
    /// `degree_v` is incompatible with the section count, compatibility
    /// processing fails, the interpolation system is singular, or the surface
    /// fails construction.
    #[allow(clippy::cast_precision_loss, clippy::needless_range_loop)]
    pub fn loft(sections: &[NurbsCurve3D], degree_v: Option<usize>) -> Result<Self> {
        if sections.len() < 2 {
            return Err(
                GeometryError::Degenerate("loft needs at least 2 section curves".into()).into(),
            );
        }
        let k = sections.len();
        let degree_v = degree_v.unwrap_or_else(|| 3.min(k - 1));
        if degree_v < 1 || degree_v > k - 1 {
            return Err(GeometryError::Degenerate(format!(
                "degree_v {degree_v} incompatible with {k} sections (need 1..={})",
                k - 1
            ))
            .into());
        }

        let compatible = make_compatible(sections)?;
        let degree_u = compatible[0].degree();
        let knots_u = compatible[0].knots().clone();
        let nu = compatible[0].control_points().len();

        // Homogeneous control nets: hp[k_idx][i] = (w * P, w).
        let nets: Vec<Vec<(Vector3, f64)>> = compatible
            .iter()
            .map(|c| {
                c.control_points()
                    .iter()
                    .zip(c.weights())
                    .map(|(p, &w)| (p.coords * w, w))
                    .collect()
            })
            .collect();

        // Section v parameters by chord length on the control-net centroids
        // (§10.3): average each section's control points, then chord-length
        // parameterize the resulting polyline.
        let centroids: Vec<Point3> = nets
            .iter()
            .map(|net| {
                let mut sum = Vector3::zeros();
                for &(wp, w) in net {
                    sum += wp / w;
                }
                Point3::from(sum / nu as f64)
            })
            .collect();
        let v_params = chord_length_params(&centroids)?;

        // Averaged v-knots (eq 9.8) shared by every u-index interpolation.
        let knots_v = averaged_knots(&v_params, degree_v)?;

        // Shared interpolation matrix: row m holds the degree-`degree_v` basis
        // values at v_params[m]. Identical for every u-index, so factor once.
        let mut mat = DMatrix::<f64>::zeros(k, k);
        for (m, &t) in v_params.iter().enumerate() {
            let span = knots_v.find_span(degree_v, k, t);
            let basis = basis_functions(&knots_v, span, t, degree_v);
            for (j, b) in basis.iter().enumerate() {
                mat[(m, span - degree_v + j)] = *b;
            }
        }
        let lu = mat.lu();

        // Solve the interpolation for each u-index over the 4 homogeneous
        // coordinates (w*x, w*y, w*z, w), forming the surface control net.
        let mut control_points = vec![Point3::origin(); nu * k];
        let mut weights = vec![0.0; nu * k];
        for i in 0..nu {
            let mut solved: [DVector<f64>; 4] = core::array::from_fn(|_| DVector::<f64>::zeros(k));
            for axis in 0..4 {
                let rhs = DVector::from_iterator(
                    k,
                    nets.iter().map(|net| {
                        let (wp, w) = net[i];
                        if axis < 3 {
                            wp[axis]
                        } else {
                            w
                        }
                    }),
                );
                solved[axis] = lu.solve(&rhs).ok_or_else(|| {
                    GeometryError::Degenerate("singular loft interpolation system".into())
                })?;
            }
            for j in 0..k {
                let w = solved[3][j];
                if w <= TOLERANCE {
                    return Err(GeometryError::Degenerate(
                        "non-positive interpolated loft weight".into(),
                    )
                    .into());
                }
                // u-major grid index = i * nv + j with nv = k.
                let idx = i * k + j;
                control_points[idx] =
                    Point3::new(solved[0][j] / w, solved[1][j] / w, solved[2][j] / w);
                weights[idx] = w;
            }
        }

        NurbsSurface::new(
            control_points,
            weights,
            nu,
            k,
            knots_u,
            knots_v,
            degree_u,
            degree_v,
        )
    }
}

/// Makes the section curves compatible: a common degree (the maximum, reached
/// by degree elevation) and a common knot vector (the union of all section
/// knot vectors over the normalized domain `[0, 1]`, reached by knot
/// insertion).
///
/// # Errors
///
/// Returns an error if degree elevation, reparameterization, or knot insertion
/// fails.
fn make_compatible(sections: &[NurbsCurve3D]) -> Result<Vec<NurbsCurve3D>> {
    // 1. Elevate all sections to the common maximum degree.
    let max_degree = sections.iter().map(NurbsCurve3D::degree).max().unwrap_or(1);
    let mut elevated: Vec<NurbsCurve3D> = Vec::with_capacity(sections.len());
    for c in sections {
        let normalized = reparameterize_unit(c)?;
        elevated.push(normalized.elevate_degree(max_degree - normalized.degree())?);
    }

    // 2. Collect the union of all interior knots (multiplicity-aware) over the
    //    shared domain [0, 1].
    let union = knot_union(&elevated, max_degree);

    // 3. Insert the missing knots into each section so all share the union.
    let mut compatible = Vec::with_capacity(elevated.len());
    for c in elevated {
        let mut curve = c;
        for &(u, target_mult) in &union {
            let have = curve.knots().multiplicity(u);
            if target_mult > have {
                curve = curve.insert_knot(u, target_mult - have)?;
            }
        }
        compatible.push(curve);
    }
    Ok(compatible)
}

/// Reparameterizes a curve's knot vector to the domain `[0, 1]` by an affine
/// map, leaving the geometry unchanged.
fn reparameterize_unit(curve: &NurbsCurve3D) -> Result<NurbsCurve3D> {
    let (a, b) = curve.parameter_domain();
    if (a - 0.0).abs() < TOLERANCE && (b - 1.0).abs() < TOLERANCE {
        return Ok(curve.clone());
    }
    let span = b - a;
    if span < TOLERANCE {
        return Err(GeometryError::Degenerate("degenerate curve domain".into()).into());
    }
    let new_knots: Vec<f64> = curve
        .knots()
        .as_slice()
        .iter()
        .map(|&t| ((t - a) / span).clamp(0.0, 1.0))
        .collect();
    NurbsCurve3D::new(
        curve.control_points().to_vec(),
        curve.weights().to_vec(),
        KnotVector::new(new_knots)?,
        curve.degree(),
    )
}

/// Computes the union of interior knots across all curves on `[0, 1]`, as
/// `(value, multiplicity)` pairs where multiplicity is the maximum seen.
fn knot_union(curves: &[NurbsCurve3D], degree: usize) -> Vec<(f64, usize)> {
    let mut union: Vec<(f64, usize)> = Vec::new();
    for curve in curves {
        let knots = curve.knots().as_slice();
        // Interior knots only: skip the clamped ends.
        for &u in &knots[degree + 1..knots.len() - degree - 1] {
            let mult = curve.knots().multiplicity(u);
            match union.iter_mut().find(|(v, _)| (*v - u).abs() < TOLERANCE) {
                Some(entry) => entry.1 = entry.1.max(mult),
                None => union.push((u, mult)),
            }
        }
    }
    union
}

/// Chord-length parameters (eq 9.5) over the polyline through `points`,
/// normalized to `[0, 1]`.
fn chord_length_params(points: &[Point3]) -> Result<Vec<f64>> {
    let n = points.len();
    let chords: Vec<f64> = points.windows(2).map(|w| (w[1] - w[0]).norm()).collect();
    let total: f64 = chords.iter().sum();
    if total < TOLERANCE {
        return Err(GeometryError::Degenerate("coincident loft sections".into()).into());
    }
    let mut params = Vec::with_capacity(n);
    params.push(0.0);
    let mut acc = 0.0;
    for chord in &chords {
        acc += chord;
        params.push(acc / total);
    }
    params[n - 1] = 1.0;
    Ok(params)
}

/// Averaged knot vector (eq 9.8) for a degree-`degree` global interpolation
/// through points at the given `params`.
#[allow(clippy::cast_precision_loss)]
fn averaged_knots(params: &[f64], degree: usize) -> Result<KnotVector> {
    let n = params.len();
    let mut knots = vec![0.0; degree + 1];
    for j in 1..(n - degree) {
        let avg: f64 = params[j..j + degree].iter().sum::<f64>() / degree as f64;
        knots.push(avg);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    KnotVector::new(knots)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn line(p0: Point3, p1: Point3) -> NurbsCurve3D {
        NurbsCurve3D::from_unweighted(
            vec![p0, p1],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap()
    }

    fn circle_at(z: f64, r: f64) -> NurbsCurve3D {
        NurbsCurve3D::circle(Point3::new(0.0, 0.0, z), r, Vector3::z(), Vector3::x()).unwrap()
    }

    #[test]
    fn loft_two_lines_is_bilinear_patch() {
        // Section 0: x-axis line at y=0; section 1: x-axis line at y=2, z=1.
        let s0 = line(Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0));
        let s1 = line(Point3::new(0.0, 2.0, 1.0), Point3::new(2.0, 2.0, 1.0));
        let surf = NurbsSurface::loft(&[s0.clone(), s1.clone()], None).unwrap();
        let ((u0, u1), (v0, v1)) = surf.parameter_domain();
        for i in 0..=10 {
            let u = u0 + (u1 - u0) * f64::from(i) / 10.0;
            let a = s0.point_at(u).unwrap();
            let b = s1.point_at(u).unwrap();
            for j in 0..=10 {
                let v = v0 + (v1 - v0) * f64::from(j) / 10.0;
                let vt = (v - v0) / (v1 - v0);
                let expected = a + (b - a) * vt;
                let got = surf.point_at(u, v).unwrap();
                assert!((got - expected).norm() < 1e-12, "bilinear off at ({u},{v})");
            }
        }
    }

    /// Sections must be reproduced exactly at their assigned v parameters.
    fn assert_sections_reproduced(
        surf: &NurbsSurface,
        sections: &[NurbsCurve3D],
        v_params: &[f64],
        tol: f64,
    ) {
        let ((u0, u1), _) = surf.parameter_domain();
        for (sec, &v) in sections.iter().zip(v_params) {
            for i in 0..=20 {
                let u = u0 + (u1 - u0) * f64::from(i) / 20.0;
                let got = surf.point_at(u, v).unwrap();
                let want = sec.point_at(u).unwrap();
                assert!((got - want).norm() < tol, "section off at (u={u}, v={v})");
            }
        }
    }

    #[test]
    fn loft_three_circles_reproduces_each_section() {
        let sections = vec![
            circle_at(0.0, 1.0),
            circle_at(1.0, 2.0),
            circle_at(2.0, 1.0),
        ];
        let surf = NurbsSurface::loft(&sections, None).unwrap();
        // v parameters are chord-length on centroids; all centroids share x=y=0
        // so they are uniform in z -> 0, 0.5, 1.
        assert_sections_reproduced(&surf, &sections, &[0.0, 0.5, 1.0], 1e-9);
    }

    #[test]
    fn loft_mixed_degrees_exercises_compatibility() {
        // Section 0: degree-1 line. Section 1: degree-2 quarter circle. Both
        // run roughly along +X then +Y; compatibility must elevate the line.
        let s0 = line(Point3::new(1.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0));
        let s1 = NurbsCurve3D::arc(
            Point3::origin(),
            1.0,
            Vector3::z(),
            Vector3::x(),
            0.0,
            std::f64::consts::FRAC_PI_2,
        )
        .unwrap();
        // Lift section 1 in z so the sections are distinct.
        let s1 = NurbsCurve3D::new(
            s1.control_points()
                .iter()
                .map(|p| Point3::new(p.x, p.y, 1.0))
                .collect(),
            s1.weights().to_vec(),
            s1.knots().clone(),
            s1.degree(),
        )
        .unwrap();
        let sections = vec![s0, s1];
        let surf = NurbsSurface::loft(&sections, None).unwrap();
        // Two sections -> v params 0 and 1.
        assert_sections_reproduced(&surf, &sections, &[0.0, 1.0], 1e-9);
    }

    #[test]
    fn rejects_single_section() {
        let s0 = line(Point3::origin(), Point3::new(1.0, 0.0, 0.0));
        assert!(NurbsSurface::loft(&[s0], None).is_err());
    }
}
