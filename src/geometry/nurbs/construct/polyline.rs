use crate::error::{GeometryError, Result};
use crate::math::{Point3, TOLERANCE};

use crate::geometry::nurbs::{KnotVector, NurbsCurve3D};

/// Builds a degree-1 NURBS curve through the given points with chord-length
/// parameterization normalized to `[0, 1]`. The curve passes through every
/// point exactly.
///
/// # Errors
///
/// Returns an error if fewer than 2 points are given or consecutive points
/// coincide.
pub fn nurbs_polyline(points: &[Point3]) -> Result<NurbsCurve3D> {
    if points.len() < 2 {
        return Err(GeometryError::Degenerate("polyline needs at least 2 points".into()).into());
    }
    let chords: Vec<f64> = points.windows(2).map(|w| (w[1] - w[0]).norm()).collect();
    if chords.iter().any(|&c| c < TOLERANCE) {
        return Err(GeometryError::Degenerate(
            "consecutive polyline points must not coincide".into(),
        )
        .into());
    }
    let total: f64 = chords.iter().sum();

    // Knots: [0, 0, t_1, ..., t_{n-2}, 1, 1] (degree 1, clamped)
    let mut knots = Vec::with_capacity(points.len() + 2);
    knots.push(0.0);
    knots.push(0.0);
    let mut acc = 0.0;
    for chord in &chords[..chords.len() - 1] {
        acc += chord;
        knots.push(acc / total);
    }
    knots.push(1.0);
    knots.push(1.0);

    NurbsCurve3D::from_unweighted(points.to_vec(), KnotVector::new(knots)?, 1)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::curve::Curve;
    use crate::math::TOLERANCE;

    #[test]
    fn polyline_passes_through_vertices() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(3.0, 2.0, 1.0),
        ];
        let c = nurbs_polyline(&pts).unwrap();
        // Vertices sit at normalized chord-length parameters.
        let chords = [1.0, 2.0, 5.0_f64.sqrt()];
        let total: f64 = chords.iter().sum();
        let params = [0.0, chords[0] / total, (chords[0] + chords[1]) / total, 1.0];
        for (p, t) in pts.iter().zip(params) {
            let q = c.point_at(t).unwrap();
            assert!((q - p).norm() < TOLERANCE, "vertex at t={t}");
        }
    }

    #[test]
    fn polyline_segments_are_straight() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 2.0, 0.0),
        ];
        let c = nurbs_polyline(&pts).unwrap();
        // Midpoint of first segment (chord param 0.25 = half of first chord)
        let p = c.point_at(0.25).unwrap();
        assert!((p - Point3::new(1.0, 0.0, 0.0)).norm() < TOLERANCE);
        assert!(!c.is_closed());
    }

    #[test]
    fn rejects_single_point() {
        assert!(nurbs_polyline(&[Point3::origin()]).is_err());
    }

    #[test]
    fn rejects_coincident_consecutive_points() {
        let r = nurbs_polyline(&[Point3::origin(), Point3::origin()]);
        assert!(r.is_err());
    }
}
