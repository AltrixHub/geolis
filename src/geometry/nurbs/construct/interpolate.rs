use nalgebra::{DMatrix, DVector};

use crate::error::{GeometryError, Result};
use crate::math::{Point3, TOLERANCE};

use crate::geometry::nurbs::{basis_functions, KnotVector, NurbsCurve3D};

/// Globally interpolates a NURBS curve of the given degree through the
/// points (The NURBS Book, A9.1) using chord-length parameterization and
/// knot averaging. Returns the curve and the parameter assigned to each
/// input point.
///
/// # Errors
///
/// Returns an error if the degree is < 1, fewer than `degree + 1` points
/// are given, the total chord length is degenerate, or the interpolation
/// system is singular.
pub fn interpolate_points(points: &[Point3], degree: usize) -> Result<(NurbsCurve3D, Vec<f64>)> {
    let n = points.len();
    if degree < 1 {
        return Err(GeometryError::Degenerate("degree must be >= 1".into()).into());
    }
    if n < degree + 1 {
        return Err(GeometryError::Degenerate(format!(
            "degree {degree} interpolation needs at least {} points, got {n}",
            degree + 1
        ))
        .into());
    }

    // Chord-length parameters (eq 9.5)
    let chords: Vec<f64> = points.windows(2).map(|w| (w[1] - w[0]).norm()).collect();
    let total: f64 = chords.iter().sum();
    if total < TOLERANCE {
        return Err(GeometryError::Degenerate("degenerate point set".into()).into());
    }
    let mut params = Vec::with_capacity(n);
    params.push(0.0);
    let mut acc = 0.0;
    for chord in &chords {
        acc += chord;
        params.push(acc / total);
    }
    params[n - 1] = 1.0;

    // Averaged knots (eq 9.8)
    let mut knots = vec![0.0; degree + 1];
    for j in 1..(n - degree) {
        let avg: f64 = params[j..j + degree].iter().sum::<f64>() / degree as f64;
        knots.push(avg);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    let knot_vector = KnotVector::new(knots)?;

    // Coefficient matrix: row i holds the basis values at params[i]
    let mut mat = DMatrix::<f64>::zeros(n, n);
    for (i, &t) in params.iter().enumerate() {
        let span = knot_vector.find_span(degree, n, t);
        let basis = basis_functions(&knot_vector, span, t, degree);
        for (j, b) in basis.iter().enumerate() {
            mat[(i, span - degree + j)] = *b;
        }
    }

    let lu = mat.lu();
    let mut coords = vec![[0.0; 3]; n];
    for axis in 0..3 {
        let rhs = DVector::from_iterator(n, points.iter().map(|p| p.coords[axis]));
        let solution = lu
            .solve(&rhs)
            .ok_or_else(|| GeometryError::Degenerate("singular interpolation system".into()))?;
        for (i, value) in solution.iter().enumerate() {
            coords[i][axis] = *value;
        }
    }
    let control_points: Vec<Point3> = coords
        .into_iter()
        .map(|[x, y, z]| Point3::new(x, y, z))
        .collect();

    let curve = NurbsCurve3D::from_unweighted(control_points, knot_vector, degree)?;
    Ok((curve, params))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn interpolated_curve_passes_through_all_points() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(3.0, 1.0, 1.0),
            Point3::new(4.0, 4.0, 0.0),
            Point3::new(6.0, 3.0, 2.0),
        ];
        let (curve, params) = interpolate_points(&pts, 3).unwrap();
        assert_eq!(params.len(), pts.len());
        for (p, &t) in pts.iter().zip(&params) {
            let q = curve.point_at(t).unwrap();
            assert!((q - p).norm() < 1e-9, "missed point at t={t}");
        }
    }

    #[test]
    fn interpolation_with_exactly_degree_plus_one_points() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
        ];
        let (curve, params) = interpolate_points(&pts, 3).unwrap();
        for (p, &t) in pts.iter().zip(&params) {
            let q = curve.point_at(t).unwrap();
            assert!((q - p).norm() < 1e-9);
        }
    }

    #[test]
    fn rejects_too_few_points() {
        let pts = vec![Point3::origin(), Point3::new(1.0, 0.0, 0.0)];
        assert!(interpolate_points(&pts, 3).is_err());
    }
}
