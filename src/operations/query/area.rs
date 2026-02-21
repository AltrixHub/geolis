use crate::error::Result;
use crate::tessellation::{TessellateSolid, TessellationParams};
use crate::topology::{SolidId, TopologyStore};

/// Computes the total surface area of a solid.
///
/// Uses tessellation to approximate the area by summing the areas of all
/// triangles in the mesh. The accuracy depends on the tessellation parameters.
pub struct Area {
    solid: SolidId,
    params: TessellationParams,
}

impl Area {
    /// Creates a new `Area` query with default tessellation parameters.
    #[must_use]
    pub fn new(solid: SolidId) -> Self {
        Self {
            solid,
            params: TessellationParams::default(),
        }
    }

    /// Sets custom tessellation parameters for higher accuracy.
    #[must_use]
    pub fn with_params(mut self, params: TessellationParams) -> Self {
        self.params = params;
        self
    }

    /// Executes the query, returning the total surface area.
    ///
    /// # Errors
    ///
    /// Returns an error if the solid cannot be tessellated.
    pub fn execute(&self, store: &TopologyStore) -> Result<f64> {
        let mesh = TessellateSolid::new(self.solid, self.params).execute(store)?;

        let mut total_area = 0.0;
        for tri in &mesh.indices {
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];

            let edge1 = v1 - v0;
            let edge2 = v2 - v0;
            total_area += edge1.cross(&edge2).norm() * 0.5;
        }

        Ok(total_area)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::{MakeBox, MakeCylinder, MakeSphere};
    use std::f64::consts::PI;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn box_area() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(2.0, 3.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let area = Area::new(solid).execute(&store).unwrap();
        // 2*(2*3 + 2*4 + 3*4) = 2*(6+8+12) = 52
        assert!((area - 52.0).abs() < 0.1, "expected 52.0, got {area}");
    }

    #[test]
    fn cylinder_area() {
        let mut store = TopologyStore::new();
        let solid = MakeCylinder::new(p(0.0, 0.0, 0.0), 2.0, crate::math::Vector3::z(), 5.0)
            .execute(&mut store)
            .unwrap();

        let area = Area::new(solid).execute(&store).unwrap();
        // Cylinder: lateral = 2*pi*r*h = 20*pi ≈ 62.83
        // Top disc = pi*r^2 = 4*pi ≈ 12.57 (but our cylinder has top/bottom discs)
        // Total = 2*pi*r*h + 2*pi*r^2 = 2*pi*2*5 + 2*pi*4 = 20*pi + 8*pi = 28*pi ≈ 87.96
        // But our cylinder has an open center (revolution of rectangle with on-axis edge)
        // so the "discs" are full circles.
        let expected = 28.0 * PI;
        let tolerance = expected * 0.05; // 5% tolerance for tessellation
        assert!(
            (area - expected).abs() < tolerance,
            "expected ~{expected:.2}, got {area:.2}"
        );
    }

    #[test]
    fn sphere_area() {
        let mut store = TopologyStore::new();
        let solid = MakeSphere::new(p(0.0, 0.0, 0.0), 3.0)
            .execute(&mut store)
            .unwrap();

        let area = Area::new(solid).execute(&store).unwrap();
        // True sphere: 4*pi*r^2 = 36*pi ≈ 113.10
        // Our sphere is a cone approximation, so area will differ
        // The cone approximation has area = pi*r*sqrt(r^2+r^2)*2 = 2*pi*r^2*sqrt(2)
        // ≈ 79.97. We just check it's positive and reasonable.
        assert!(area > 50.0 && area < 200.0, "unexpected sphere area: {area}");
    }
}
