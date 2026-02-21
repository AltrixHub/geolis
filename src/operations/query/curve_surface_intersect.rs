use crate::error::Result;
use crate::geometry::curve::{Curve, Line};
use crate::geometry::surface::{Cone, Cylinder, Plane, Sphere};
use crate::math::{Point3, TOLERANCE};

/// A single intersection point between a curve and a surface.
#[derive(Debug, Clone, Copy)]
pub struct CurveSurfaceHit {
    /// Parameter on the curve.
    pub t: f64,
    /// U parameter on the surface.
    pub u: f64,
    /// V parameter on the surface.
    pub v: f64,
    /// 3D intersection point.
    pub point: Point3,
}

/// Computes intersection points between a line and a surface.
///
/// Uses analytic solutions for Line-Plane, Line-Cylinder, Line-Sphere,
/// and Line-Cone intersections.
pub struct LineSurfaceIntersect {
    line: Line,
    t_min: f64,
    t_max: f64,
}

impl LineSurfaceIntersect {
    /// Creates a new intersection query.
    ///
    /// `t_min` and `t_max` define the line segment to check.
    #[must_use]
    pub fn new(line: Line, t_min: f64, t_max: f64) -> Self {
        Self { line, t_min, t_max }
    }

    /// Intersects the line with a Plane.
    pub fn with_plane(&self, plane: &Plane) -> Result<Vec<CurveSurfaceHit>> {
        intersect_line_plane(&self.line, self.t_min, self.t_max, plane)
    }

    /// Intersects the line with a Cylinder.
    pub fn with_cylinder(&self, cyl: &Cylinder) -> Result<Vec<CurveSurfaceHit>> {
        intersect_line_cylinder(&self.line, self.t_min, self.t_max, cyl)
    }

    /// Intersects the line with a Sphere.
    pub fn with_sphere(&self, sph: &Sphere) -> Result<Vec<CurveSurfaceHit>> {
        intersect_line_sphere(&self.line, self.t_min, self.t_max, sph)
    }

    /// Intersects the line with a Cone.
    pub fn with_cone(&self, cone: &Cone) -> Result<Vec<CurveSurfaceHit>> {
        intersect_line_cone(&self.line, self.t_min, self.t_max, cone)
    }
}

/// Line-Plane intersection: solve `(O + t*D - P0) . N = 0`.
fn intersect_line_plane(
    line: &Line,
    t_min: f64,
    t_max: f64,
    plane: &Plane,
) -> Result<Vec<CurveSurfaceHit>> {
    let origin = line.evaluate(0.0)?;
    let dir = line.tangent(0.0)?;
    let normal = *plane.plane_normal();

    let denom = dir.dot(&normal);
    if denom.abs() < TOLERANCE {
        // Line is parallel to plane
        return Ok(vec![]);
    }

    let t = (plane.origin() - origin).dot(&normal) / denom;
    if t < t_min - TOLERANCE || t > t_max + TOLERANCE {
        return Ok(vec![]);
    }

    let point = line.evaluate(t)?;
    let dp = point - plane.origin();
    let u = dp.dot(plane.u_dir());
    let v = dp.dot(plane.v_dir());

    Ok(vec![CurveSurfaceHit { t, u, v, point }])
}

/// Line-Cylinder intersection: reduce to 2D circle-line problem.
fn intersect_line_cylinder(
    line: &Line,
    t_min: f64,
    t_max: f64,
    cyl: &Cylinder,
) -> Result<Vec<CurveSurfaceHit>> {
    let origin = line.evaluate(0.0)?;
    let dir = line.tangent(0.0)?;

    let axis = cyl.axis();
    let center = cyl.center();
    let r = cyl.radius();

    // Project everything onto the plane perpendicular to the cylinder axis
    let dp = origin - center;
    let dp_perp = dp - *axis * dp.dot(axis);
    let dir_perp = dir - *axis * dir.dot(axis);

    // Solve |dp_perp + t * dir_perp|^2 = r^2
    let a = dir_perp.dot(&dir_perp);
    let b = 2.0 * dp_perp.dot(&dir_perp);
    let c = dp_perp.dot(&dp_perp) - r * r;

    solve_quadratic_hits(line, t_min, t_max, a, b, c, |pt| cyl.inverse(pt))
}

/// Line-Sphere intersection: solve `|O + t*D - C|^2 = r^2`.
fn intersect_line_sphere(
    line: &Line,
    t_min: f64,
    t_max: f64,
    sph: &Sphere,
) -> Result<Vec<CurveSurfaceHit>> {
    let origin = line.evaluate(0.0)?;
    let dir = line.tangent(0.0)?;

    let dp = origin - sph.center();
    let r = sph.radius();

    let a = dir.dot(&dir);
    let b = 2.0 * dp.dot(&dir);
    let c = dp.dot(&dp) - r * r;

    solve_quadratic_hits(line, t_min, t_max, a, b, c, |pt| sph.inverse(pt))
}

/// Line-Cone intersection: solve the cone equation.
fn intersect_line_cone(
    line: &Line,
    t_min: f64,
    t_max: f64,
    cone: &Cone,
) -> Result<Vec<CurveSurfaceHit>> {
    let origin = line.evaluate(0.0)?;
    let dir = line.tangent(0.0)?;

    let apex = cone.apex();
    let axis = cone.axis();
    let ca = cone.half_angle().cos();
    let cos2 = ca * ca;

    let dp = origin - apex;

    // Cone equation: (D.axis)^2 * cos^2 = D.D * cos^2 + ...
    // Expanding: |P-apex|^2 * cos^2(alpha) = ((P-apex).axis)^2
    // Substitute P = O + t*dir:
    let d_dot_a = dir.dot(axis);
    let dp_dot_a = dp.dot(axis);

    let a = d_dot_a * d_dot_a - dir.dot(&dir) * cos2;
    let b = 2.0 * (d_dot_a * dp_dot_a - dp.dot(&dir) * cos2);
    let c = dp_dot_a * dp_dot_a - dp.dot(&dp) * cos2;

    // The cone equation gives both nappes. Filter to v > 0 (forward nappe).
    let mut hits = Vec::new();
    let disc = b * b - 4.0 * a * c;
    if disc < -TOLERANCE {
        return Ok(hits);
    }
    let disc = disc.max(0.0).sqrt();

    let candidates = if a.abs() < TOLERANCE {
        // Linear case: single solution
        if b.abs() < TOLERANCE {
            vec![]
        } else {
            vec![-c / b]
        }
    } else {
        vec![(-b - disc) / (2.0 * a), (-b + disc) / (2.0 * a)]
    };

    for t in candidates {
        if t < t_min - TOLERANCE || t > t_max + TOLERANCE {
            continue;
        }
        let point = line.evaluate(t)?;
        let (u, v) = cone.inverse(&point);
        // Only forward nappe
        if v > -TOLERANCE {
            hits.push(CurveSurfaceHit { t, u, v, point });
        }
    }

    // Sort by t and remove near-duplicates
    hits.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
    if hits.len() == 2 && (hits[0].t - hits[1].t).abs() < TOLERANCE {
        hits.pop();
    }

    Ok(hits)
}

/// Solves a quadratic `a*t^2 + b*t + c = 0` and returns intersection hits.
fn solve_quadratic_hits(
    line: &Line,
    t_min: f64,
    t_max: f64,
    a: f64,
    b: f64,
    c: f64,
    inverse: impl Fn(&Point3) -> (f64, f64),
) -> Result<Vec<CurveSurfaceHit>> {
    let mut hits = Vec::new();
    let disc = b * b - 4.0 * a * c;

    if disc < -TOLERANCE {
        return Ok(hits);
    }

    if a.abs() < TOLERANCE {
        // Degenerate: linear equation
        if b.abs() < TOLERANCE {
            return Ok(hits);
        }
        let t = -c / b;
        if t >= t_min - TOLERANCE && t <= t_max + TOLERANCE {
            let point = line.evaluate(t)?;
            let (u, v) = inverse(&point);
            hits.push(CurveSurfaceHit { t, u, v, point });
        }
        return Ok(hits);
    }

    let disc = disc.max(0.0).sqrt();
    let t1 = (-b - disc) / (2.0 * a);
    let t2 = (-b + disc) / (2.0 * a);

    for t in [t1, t2] {
        if t >= t_min - TOLERANCE && t <= t_max + TOLERANCE {
            let point = line.evaluate(t)?;
            let (u, v) = inverse(&point);
            hits.push(CurveSurfaceHit { t, u, v, point });
        }
    }

    // Sort by t and remove near-duplicate hits (tangent case)
    hits.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
    if hits.len() == 2 && (hits[0].t - hits[1].t).abs() < TOLERANCE {
        hits.pop();
    }

    Ok(hits)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::curve::Line;
    use crate::math::Vector3;

    #[test]
    fn line_through_sphere_gives_2_hits() {
        let sph = Sphere::new(Point3::origin(), 5.0, Vector3::z(), Vector3::x()).unwrap();
        let line = Line::new(Point3::new(-10.0, 0.0, 0.0), Vector3::x()).unwrap();

        let hits = LineSurfaceIntersect::new(line, 0.0, 20.0)
            .with_sphere(&sph)
            .unwrap();

        assert_eq!(hits.len(), 2);
        // Entry at t=5, exit at t=15
        assert!((hits[0].t - 5.0).abs() < 1e-6, "t0 = {}", hits[0].t);
        assert!((hits[1].t - 15.0).abs() < 1e-6, "t1 = {}", hits[1].t);
    }

    #[test]
    fn line_tangent_to_sphere_gives_1_hit() {
        let sph = Sphere::new(Point3::origin(), 5.0, Vector3::z(), Vector3::x()).unwrap();
        let line = Line::new(Point3::new(-10.0, 5.0, 0.0), Vector3::x()).unwrap();

        let hits = LineSurfaceIntersect::new(line, 0.0, 20.0)
            .with_sphere(&sph)
            .unwrap();

        assert_eq!(hits.len(), 1);
        assert!((hits[0].point.y - 5.0).abs() < 1e-6);
    }

    #[test]
    fn line_misses_sphere_gives_0_hits() {
        let sph = Sphere::new(Point3::origin(), 5.0, Vector3::z(), Vector3::x()).unwrap();
        let line = Line::new(Point3::new(-10.0, 6.0, 0.0), Vector3::x()).unwrap();

        let hits = LineSurfaceIntersect::new(line, 0.0, 20.0)
            .with_sphere(&sph)
            .unwrap();

        assert_eq!(hits.len(), 0);
    }

    #[test]
    fn line_through_cylinder_gives_2_hits() {
        let cyl = Cylinder::new(Point3::origin(), 3.0, Vector3::z(), Vector3::x()).unwrap();
        let line = Line::new(Point3::new(-10.0, 0.0, 2.0), Vector3::x()).unwrap();

        let hits = LineSurfaceIntersect::new(line, 0.0, 20.0)
            .with_cylinder(&cyl)
            .unwrap();

        assert_eq!(hits.len(), 2);
        // Entry at x=-3, exit at x=3
        assert!((hits[0].point.x - (-3.0)).abs() < 1e-6);
        assert!((hits[1].point.x - 3.0).abs() < 1e-6);
    }

    #[test]
    fn line_plane_intersection() {
        let plane = Plane::from_normal(Point3::new(0.0, 0.0, 5.0), Vector3::z()).unwrap();
        let line = Line::new(Point3::origin(), Vector3::z()).unwrap();

        let hits = LineSurfaceIntersect::new(line, 0.0, 10.0)
            .with_plane(&plane)
            .unwrap();

        assert_eq!(hits.len(), 1);
        assert!((hits[0].t - 5.0).abs() < 1e-6);
        assert!((hits[0].point.z - 5.0).abs() < 1e-6);
    }

    #[test]
    fn line_parallel_to_plane_gives_0_hits() {
        let plane = Plane::from_normal(Point3::new(0.0, 0.0, 5.0), Vector3::z()).unwrap();
        let line = Line::new(Point3::origin(), Vector3::x()).unwrap();

        let hits = LineSurfaceIntersect::new(line, 0.0, 10.0)
            .with_plane(&plane)
            .unwrap();

        assert_eq!(hits.len(), 0);
    }

    #[test]
    fn line_through_cone_gives_2_hits() {
        let cone =
            Cone::new(Point3::origin(), Vector3::z(), std::f64::consts::FRAC_PI_4, Vector3::x())
                .unwrap();
        let line = Line::new(Point3::new(-10.0, 0.0, 5.0), Vector3::x()).unwrap();

        let hits = LineSurfaceIntersect::new(line, 0.0, 20.0)
            .with_cone(&cone)
            .unwrap();

        // At z=5, cone radius = 5*tan(45°) = 5
        assert_eq!(hits.len(), 2);
        assert!((hits[0].point.x - (-5.0)).abs() < 1e-4);
        assert!((hits[1].point.x - 5.0).abs() < 1e-4);
    }

    #[test]
    fn line_outside_t_range_gives_0_hits() {
        let sph = Sphere::new(Point3::origin(), 5.0, Vector3::z(), Vector3::x()).unwrap();
        let line = Line::new(Point3::new(-10.0, 0.0, 0.0), Vector3::x()).unwrap();

        // Intersection at t=5 and t=15, but range is [0, 3]
        let hits = LineSurfaceIntersect::new(line, 0.0, 3.0)
            .with_sphere(&sph)
            .unwrap();

        assert_eq!(hits.len(), 0);
    }
}
