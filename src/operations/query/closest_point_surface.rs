use crate::error::Result;
use crate::geometry::surface::Surface;
use crate::math::{Point3, TOLERANCE};
use crate::topology::{FaceId, FaceSurface, TopologyStore};

/// Result of a closest-point-on-surface query.
#[derive(Debug, Clone, Copy)]
pub struct SurfacePoint {
    /// U parameter on the surface.
    pub u: f64,
    /// V parameter on the surface.
    pub v: f64,
    /// 3D point on the surface.
    pub point: Point3,
    /// Distance from the query point to the surface point.
    pub distance: f64,
}

/// Finds the closest point on a surface to a given query point.
///
/// Uses analytic solutions for Plane, Cylinder, Sphere, and Cone.
pub struct ClosestPointOnSurface {
    face: FaceId,
    query: Point3,
}

impl ClosestPointOnSurface {
    /// Creates a new query.
    #[must_use]
    pub fn new(face: FaceId, query: Point3) -> Self {
        Self { face, query }
    }

    /// Executes the query.
    ///
    /// # Errors
    ///
    /// Returns an error if the face is not found or the surface type
    /// is unsupported.
    pub fn execute(&self, store: &TopologyStore) -> Result<SurfacePoint> {
        let face = store.face(self.face)?;
        match &face.surface {
            FaceSurface::Plane(plane) => closest_on_plane(plane, &self.query),
            FaceSurface::Cylinder(cyl) => closest_on_cylinder(cyl, &self.query),
            FaceSurface::Sphere(sph) => closest_on_sphere(sph, &self.query),
            FaceSurface::Cone(cone) => closest_on_cone(cone, &self.query),
            FaceSurface::Torus(torus) => closest_on_torus(torus, &self.query),
        }
    }
}

fn closest_on_plane(
    plane: &crate::geometry::surface::Plane,
    query: &Point3,
) -> Result<SurfacePoint> {
    let dp = query - plane.origin();
    let u = dp.dot(plane.u_dir());
    let v = dp.dot(plane.v_dir());
    let point = plane.evaluate(u, v)?;
    let distance = (query - point).norm();
    Ok(SurfacePoint {
        u,
        v,
        point,
        distance,
    })
}

fn closest_on_cylinder(
    cyl: &crate::geometry::surface::Cylinder,
    query: &Point3,
) -> Result<SurfacePoint> {
    let dp = query - cyl.center();
    let v = dp.dot(cyl.axis());
    let foot = cyl.center() + cyl.axis() * v;
    let radial = query - foot;
    let radial_len = radial.norm();

    let point = if radial_len < TOLERANCE {
        // Query is on the axis; pick the ref_dir direction
        foot + cyl.ref_dir() * cyl.radius()
    } else {
        foot + radial * (cyl.radius() / radial_len)
    };

    let (u, v_param) = cyl.inverse(&point);
    let distance = (query - point).norm();
    Ok(SurfacePoint {
        u,
        v: v_param,
        point,
        distance,
    })
}

fn closest_on_sphere(
    sph: &crate::geometry::surface::Sphere,
    query: &Point3,
) -> Result<SurfacePoint> {
    let dp = query - sph.center();
    let dp_len = dp.norm();

    let point = if dp_len < TOLERANCE {
        // Query is at center; pick the ref_dir direction
        *sph.center() + *sph.ref_dir() * sph.radius()
    } else {
        *sph.center() + dp * (sph.radius() / dp_len)
    };

    let (u, v) = sph.inverse(&point);
    let distance = (query - point).norm();
    Ok(SurfacePoint {
        u,
        v,
        point,
        distance,
    })
}

fn closest_on_cone(
    cone: &crate::geometry::surface::Cone,
    query: &Point3,
) -> Result<SurfacePoint> {
    let dp = query - cone.apex();
    let axis_proj = dp.dot(cone.axis());
    let radial = dp - *cone.axis() * axis_proj;
    let radial_len = radial.norm();

    // Project onto the cone surface: find the closest point on the generator line
    let sa = cone.half_angle().sin();
    let ca = cone.half_angle().cos();

    // The generator direction at the query's azimuthal angle
    let (u, radial_dir) = if radial_len < TOLERANCE {
        (0.0, *cone.ref_dir())
    } else {
        let rd = radial / radial_len;
        let binormal = cone.axis().cross(cone.ref_dir());
        let u = dp.dot(&binormal).atan2(dp.dot(cone.ref_dir()));
        (u, rd)
    };

    // Generator direction: cos(α)*axis + sin(α)*radial_dir
    let gen_dir = *cone.axis() * ca + radial_dir * sa;

    // Project dp onto the generator direction to find v
    let v = dp.dot(&gen_dir).max(0.0);
    let point = *cone.apex() + gen_dir * v;

    let distance = (query - point).norm();
    Ok(SurfacePoint {
        u,
        v,
        point,
        distance,
    })
}

fn closest_on_torus(
    torus: &crate::geometry::surface::Torus,
    query: &Point3,
) -> Result<SurfacePoint> {
    // Use inverse() as initial estimate, then evaluate
    let (u, v) = torus.inverse(query);
    let point = torus.evaluate(u, v)?;
    let distance = (query - point).norm();
    Ok(SurfacePoint {
        u,
        v,
        point,
        distance,
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::surface::{Cylinder, Plane, Sphere};
    use crate::math::Vector3;
    use crate::topology::{FaceData, FaceSurface, TopologyStore, VertexData, WireData};

    fn make_plane_face(store: &mut TopologyStore) -> FaceId {
        let plane = Plane::from_normal(Point3::origin(), Vector3::z()).unwrap();
        let v0 = store.add_vertex(VertexData::new(Point3::new(-10.0, -10.0, 0.0)));
        let wire = store.add_wire(WireData {
            edges: vec![],
            is_closed: true,
        });
        // Add dummy vertex to avoid empty wire issues in other contexts
        let _ = v0;
        store.add_face(FaceData {
            surface: FaceSurface::Plane(plane),
            outer_wire: wire,
            inner_wires: vec![],
            same_sense: true,
        })
    }

    fn make_cylinder_face(store: &mut TopologyStore) -> FaceId {
        let cyl = Cylinder::new(Point3::origin(), 5.0, Vector3::z(), Vector3::x()).unwrap();
        let wire = store.add_wire(WireData {
            edges: vec![],
            is_closed: true,
        });
        store.add_face(FaceData {
            surface: FaceSurface::Cylinder(cyl),
            outer_wire: wire,
            inner_wires: vec![],
            same_sense: true,
        })
    }

    fn make_sphere_face(store: &mut TopologyStore) -> FaceId {
        let sph = Sphere::new(Point3::origin(), 5.0, Vector3::z(), Vector3::x()).unwrap();
        let wire = store.add_wire(WireData {
            edges: vec![],
            is_closed: true,
        });
        store.add_face(FaceData {
            surface: FaceSurface::Sphere(sph),
            outer_wire: wire,
            inner_wires: vec![],
            same_sense: true,
        })
    }

    #[test]
    fn plane_closest_point() {
        let mut store = TopologyStore::new();
        let face = make_plane_face(&mut store);

        let result = ClosestPointOnSurface::new(face, Point3::new(3.0, 4.0, 7.0))
            .execute(&store)
            .unwrap();

        assert!((result.point.z).abs() < 1e-10);
        assert!((result.distance - 7.0).abs() < 1e-10);
    }

    #[test]
    fn cylinder_closest_from_outside() {
        let mut store = TopologyStore::new();
        let face = make_cylinder_face(&mut store);

        // Point at (10, 0, 3) → closest on cylinder is (5, 0, 3)
        let result = ClosestPointOnSurface::new(face, Point3::new(10.0, 0.0, 3.0))
            .execute(&store)
            .unwrap();

        assert!((result.point.x - 5.0).abs() < 1e-6);
        assert!((result.point.y).abs() < 1e-6);
        assert!((result.point.z - 3.0).abs() < 1e-6);
        assert!((result.distance - 5.0).abs() < 1e-6);
    }

    #[test]
    fn sphere_closest_from_center() {
        let mut store = TopologyStore::new();
        let face = make_sphere_face(&mut store);

        // Point at center → closest is along ref_dir
        let result = ClosestPointOnSurface::new(face, Point3::origin())
            .execute(&store)
            .unwrap();

        assert!((result.distance - 5.0).abs() < 1e-6);
    }

    #[test]
    fn sphere_closest_from_outside() {
        let mut store = TopologyStore::new();
        let face = make_sphere_face(&mut store);

        // Point at (10, 0, 0) → closest on sphere is (5, 0, 0)
        let result = ClosestPointOnSurface::new(face, Point3::new(10.0, 0.0, 0.0))
            .execute(&store)
            .unwrap();

        assert!((result.point.x - 5.0).abs() < 1e-6);
        assert!((result.distance - 5.0).abs() < 1e-6);
    }
}
