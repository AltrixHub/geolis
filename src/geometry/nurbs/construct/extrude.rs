use crate::error::{GeometryError, Result};
use crate::math::{Vector3, TOLERANCE};

use crate::geometry::nurbs::{KnotVector, NurbsCurve3D, NurbsSurface};

impl NurbsSurface {
    /// Builds the translational extrusion of a profile curve along `direction`
    /// (The NURBS Book, §8.3).
    ///
    /// The u direction inherits the profile's degree and knot vector; the v
    /// direction is linear (degree 1) between the profile and its translate.
    /// The extrusion is exact for rational profiles: the v=0 row reproduces the
    /// profile and the v=1 row reproduces the profile translated by
    /// `direction`, with matching weights so every surface point equals
    /// `profile(u) + v * direction`.
    ///
    /// # Errors
    ///
    /// Returns an error if `direction` is zero-length or the surface fails
    /// construction.
    pub fn extrude(profile: &NurbsCurve3D, direction: Vector3) -> Result<Self> {
        if direction.norm() < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }

        let nu = profile.control_points().len();
        let nv = 2;

        // Control grid (u-major, index = i * nv + j): row j=0 is the profile,
        // row j=1 is the profile translated by `direction`. Weights repeat the
        // profile weight along v so the linear blend stays rational-exact.
        let mut control_points = Vec::with_capacity(nu * nv);
        let mut weights = Vec::with_capacity(nu * nv);
        for (point, &w) in profile.control_points().iter().zip(profile.weights()) {
            control_points.push(*point);
            control_points.push(*point + direction);
            weights.push(w);
            weights.push(w);
        }

        let knots_v = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0])?;

        NurbsSurface::new(
            control_points,
            weights,
            nu,
            nv,
            profile.knots().clone(),
            knots_v,
            profile.degree(),
            1,
        )
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::surface::Surface;
    use crate::math::Point3;

    fn quarter_circle() -> NurbsCurve3D {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsCurve3D::new(
            vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![1.0, w, 1.0],
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap(),
            2,
        )
        .unwrap()
    }

    #[test]
    fn extruded_rational_curve_is_exact_translation() {
        let profile = quarter_circle();
        let dir = Vector3::new(0.0, 0.0, 2.0);
        let s = NurbsSurface::extrude(&profile, dir).unwrap();
        let ((u0, u1), (v0, v1)) = s.parameter_domain();
        for i in 0..=10 {
            let u = u0 + (u1 - u0) * f64::from(i) / 10.0;
            let profile_pt = profile.point_at(u).unwrap();
            for j in 0..=10 {
                let v = v0 + (v1 - v0) * f64::from(j) / 10.0;
                let expected = profile_pt + dir * v;
                let got = s.point_at(u, v).unwrap();
                assert!((got - expected).norm() < 1e-12, "mismatch at ({u},{v})");
            }
        }
    }

    #[test]
    fn extruded_line_is_planar_with_constant_normal() {
        let line = NurbsCurve3D::from_unweighted(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap();
        let s = NurbsSurface::extrude(&line, Vector3::new(0.0, 0.0, 3.0)).unwrap();
        let n0 = Surface::normal(&s, 0.5, 0.5).unwrap();
        for &(u, v) in &[(0.1, 0.2), (0.9, 0.8), (0.5, 0.0), (0.5, 1.0)] {
            let n = Surface::normal(&s, u, v).unwrap();
            assert!((n - n0).norm() < 1e-12, "normal varies at ({u},{v})");
        }
    }

    #[test]
    fn rejects_zero_direction() {
        let profile = quarter_circle();
        assert!(NurbsSurface::extrude(&profile, Vector3::zeros()).is_err());
    }
}
