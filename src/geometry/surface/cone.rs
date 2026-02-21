use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use super::{Surface, SurfaceDomain};

/// A conical surface in 3D space.
///
/// Defined by an apex point, an axis direction, a half-angle, and a
/// reference direction for u=0.
///
/// `P(u, v) = apex + v * (cos(alpha) * axis + sin(alpha) * (cos(u) * ref_dir + sin(u) * binormal))`
/// where `binormal = axis x ref_dir` and `alpha` is the half-angle.
///
/// The parameter `v >= 0` measures distance along the generator from the apex.
#[derive(Debug, Clone)]
pub struct Cone {
    apex: Point3,
    axis: Vector3,
    half_angle: f64,
    ref_dir: Vector3,
}

impl Cone {
    /// Creates a new cone.
    ///
    /// # Arguments
    ///
    /// * `apex` - The apex (tip) of the cone
    /// * `axis` - Axis direction from apex outward (will be normalized)
    /// * `half_angle` - Half-angle in radians (must be in `(0, pi/2)`)
    /// * `ref_dir` - Reference direction for u=0 (must be perpendicular to axis)
    ///
    /// # Errors
    ///
    /// Returns an error if the half-angle is out of range, axis is zero-length,
    /// or the reference direction is not perpendicular to the axis.
    pub fn new(
        apex: Point3,
        axis: Vector3,
        half_angle: f64,
        ref_dir: Vector3,
    ) -> Result<Self> {
        if half_angle <= TOLERANCE || half_angle >= std::f64::consts::FRAC_PI_2 - TOLERANCE {
            return Err(GeometryError::Degenerate(
                "cone half-angle must be in (0, pi/2)".into(),
            )
            .into());
        }

        let axis_len = axis.norm();
        if axis_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let axis = axis / axis_len;

        let ref_len = ref_dir.norm();
        if ref_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let ref_dir = ref_dir / ref_len;

        if axis.dot(&ref_dir).abs() > TOLERANCE {
            return Err(GeometryError::Degenerate(
                "reference direction must be perpendicular to axis".into(),
            )
            .into());
        }

        Ok(Self {
            apex,
            axis,
            half_angle,
            ref_dir,
        })
    }

    /// Returns the apex point.
    #[must_use]
    pub fn apex(&self) -> &Point3 {
        &self.apex
    }

    /// Returns the axis direction (unit vector).
    #[must_use]
    pub fn axis(&self) -> &Vector3 {
        &self.axis
    }

    /// Returns the half-angle in radians.
    #[must_use]
    pub fn half_angle(&self) -> f64 {
        self.half_angle
    }

    /// Returns the reference direction (u=0).
    #[must_use]
    pub fn ref_dir(&self) -> &Vector3 {
        &self.ref_dir
    }

    /// Computes the binormal direction (`axis x ref_dir`).
    fn binormal(&self) -> Vector3 {
        self.axis.cross(&self.ref_dir)
    }

    /// Computes the (u, v) parameters for a given 3D point on the surface.
    ///
    /// - `u` = angle around axis in `(-pi, pi]`
    /// - `v` = distance along the generator from the apex (>= 0)
    #[must_use]
    pub fn inverse(&self, point: &Point3) -> (f64, f64) {
        let dp = point - self.apex;
        let binormal = self.binormal();
        // v = distance along the generator
        let v = dp.norm();
        // For u, project dp onto the radial plane perpendicular to axis
        let radial_x = dp.dot(&self.ref_dir);
        let radial_y = dp.dot(&binormal);
        let u = radial_y.atan2(radial_x);
        (u, v)
    }
}

impl Surface for Cone {
    fn evaluate(&self, u: f64, v: f64) -> Result<Point3> {
        let binormal = self.binormal();
        let ca = self.half_angle.cos();
        let sa = self.half_angle.sin();
        let radial = self.ref_dir * u.cos() + binormal * u.sin();
        let dir = self.axis * ca + radial * sa;
        Ok(self.apex + dir * v)
    }

    fn normal(&self, u: f64, v: f64) -> Result<Vector3> {
        if v.abs() < TOLERANCE {
            return Err(GeometryError::Degenerate(
                "cone normal is degenerate at apex".into(),
            )
            .into());
        }
        let binormal = self.binormal();
        let ca = self.half_angle.cos();
        let sa = self.half_angle.sin();
        let radial = self.ref_dir * u.cos() + binormal * u.sin();
        // Outward normal: perpendicular to generator, pointing away from axis
        // N = sin(alpha) * axis_component - cos(alpha) * radial_component... actually:
        // The surface tangent in u-direction: dP/du = v * sin(alpha) * (-sin(u)*ref + cos(u)*binorm)
        // The surface tangent in v-direction: dP/dv = cos(alpha)*axis + sin(alpha)*radial
        // Normal = dP/du x dP/dv (outward)
        let du = (-self.ref_dir * u.sin() + binormal * u.cos()) * (v * sa);
        let dv = self.axis * ca + radial * sa;
        let n = du.cross(&dv);
        let len = n.norm();
        if len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        Ok(n / len)
    }

    fn domain(&self) -> SurfaceDomain {
        SurfaceDomain::new(0.0, std::f64::consts::TAU, 0.0, f64::INFINITY)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_4, FRAC_PI_2, TAU};

    fn z_cone_45() -> Cone {
        Cone::new(
            Point3::origin(),
            Vector3::z(),
            FRAC_PI_4,
            Vector3::x(),
        )
        .unwrap()
    }

    #[test]
    fn evaluate_at_apex() {
        let c = z_cone_45();
        let p = c.evaluate(0.0, 0.0).unwrap();
        assert!((p - Point3::origin()).norm() < TOLERANCE);
    }

    #[test]
    fn evaluate_along_generator() {
        let c = z_cone_45();
        // At u=0, v=1: P = (sin(45), 0, cos(45))
        let p = c.evaluate(0.0, 1.0).unwrap();
        let s = FRAC_PI_4.sin();
        let co = FRAC_PI_4.cos();
        assert!((p - Point3::new(s, 0.0, co)).norm() < 1e-9);
    }

    #[test]
    fn evaluate_at_pi_over_2() {
        let c = z_cone_45();
        let p = c.evaluate(FRAC_PI_2, 1.0).unwrap();
        let s = FRAC_PI_4.sin();
        let co = FRAC_PI_4.cos();
        assert!((p - Point3::new(0.0, s, co)).norm() < 1e-9);
    }

    #[test]
    fn normal_degenerate_at_apex() {
        let c = z_cone_45();
        let r = c.normal(0.0, 0.0);
        assert!(r.is_err());
    }

    #[test]
    fn normal_outward() {
        let c = z_cone_45();
        let n = c.normal(0.0, 1.0).unwrap();
        // Outward normal should have positive x-component at u=0
        assert!(n.x > 0.0);
    }

    #[test]
    fn domain_ranges() {
        let c = z_cone_45();
        let d = c.domain();
        assert!((d.u_min).abs() < TOLERANCE);
        assert!((d.u_max - TAU).abs() < TOLERANCE);
        assert!((d.v_min).abs() < TOLERANCE);
        assert!(d.v_max.is_infinite());
    }

    #[test]
    fn invalid_half_angle_zero() {
        let r = Cone::new(Point3::origin(), Vector3::z(), 0.0, Vector3::x());
        assert!(r.is_err());
    }

    #[test]
    fn invalid_half_angle_90() {
        let r = Cone::new(Point3::origin(), Vector3::z(), FRAC_PI_2, Vector3::x());
        assert!(r.is_err());
    }

    #[test]
    fn inverse_roundtrip() {
        let c = z_cone_45();
        for &(u, v) in &[(0.0, 1.0), (FRAC_PI_2, 2.0), (1.0, 3.0), (TAU * 0.75, 0.5)] {
            let p = c.evaluate(u, v).unwrap();
            let (u2, v2) = c.inverse(&p);
            let p2 = c.evaluate(u2, v2).unwrap();
            assert!((p - p2).norm() < 1e-9, "roundtrip failed for u={u}, v={v}");
        }
    }
}
