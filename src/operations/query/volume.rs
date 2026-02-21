use crate::error::Result;
use crate::tessellation::{TessellateSolid, TessellationParams};
use crate::topology::{SolidId, TopologyStore};

/// Computes the volume of a solid.
///
/// Uses tessellation and the signed tetrahedron method. For each triangle,
/// computes `(1/6) * v0 . (v1 x v2)` and sums over all triangles.
///
/// The mesh normals are used to correct for any winding inconsistencies
/// between faces, making the computation robust even when face tessellations
/// have mixed winding orders.
pub struct Volume {
    solid: SolidId,
    params: TessellationParams,
}

impl Volume {
    /// Creates a new `Volume` query with default tessellation parameters.
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

    /// Executes the query, returning the volume (absolute value).
    ///
    /// # Errors
    ///
    /// Returns an error if the solid cannot be tessellated.
    pub fn execute(&self, store: &TopologyStore) -> Result<f64> {
        let mesh = TessellateSolid::new(self.solid, self.params).execute(store)?;

        let mut signed_volume = 0.0;
        for tri in &mesh.indices {
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];

            // Signed volume of tetrahedron formed by triangle and origin
            let cross = (v1 - v0).cross(&(v2 - v0));
            let det = v0.coords.dot(&v1.coords.cross(&v2.coords));

            // Use mesh normals to determine the correct sign.
            // If the triangle's geometric normal (from winding) disagrees with
            // the stored mesh normal, flip the contribution sign.
            let avg_normal = mesh.normals[tri[0] as usize]
                + mesh.normals[tri[1] as usize]
                + mesh.normals[tri[2] as usize];

            if avg_normal.dot(&cross) >= 0.0 {
                signed_volume += det;
            } else {
                signed_volume -= det;
            }
        }

        Ok(signed_volume.abs() / 6.0)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::{MakeBox, MakeCone, MakeCylinder, MakeSphere};
    use std::f64::consts::PI;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn box_volume() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(2.0, 3.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let volume = Volume::new(solid).execute(&store).unwrap();
        assert!((volume - 24.0).abs() < 0.1, "expected 24.0, got {volume}");
    }

    #[test]
    fn cylinder_volume() {
        let mut store = TopologyStore::new();
        let solid = MakeCylinder::new(p(0.0, 0.0, 0.0), 2.0, crate::math::Vector3::z(), 5.0)
            .execute(&mut store)
            .unwrap();

        let volume = Volume::new(solid).execute(&store).unwrap();
        // pi * r^2 * h = pi * 4 * 5 = 20*pi ≈ 62.83
        let expected = 20.0 * PI;
        let tolerance = expected * 0.05;
        assert!(
            (volume - expected).abs() < tolerance,
            "expected ~{expected:.2}, got {volume:.2}"
        );
    }

    #[test]
    fn cone_volume() {
        let mut store = TopologyStore::new();
        let solid = MakeCone::new(p(0.0, 0.0, 0.0), 3.0, 0.0, crate::math::Vector3::z(), 4.0)
            .execute(&mut store)
            .unwrap();

        let volume = Volume::new(solid).execute(&store).unwrap();
        // (1/3) * pi * r^2 * h = (1/3) * pi * 9 * 4 = 12*pi ≈ 37.70
        let expected = 12.0 * PI;
        let tolerance = expected * 0.05;
        assert!(
            (volume - expected).abs() < tolerance,
            "expected ~{expected:.2}, got {volume:.2}"
        );
    }

    #[test]
    fn sphere_volume() {
        let mut store = TopologyStore::new();
        let solid = MakeSphere::new(p(0.0, 0.0, 0.0), 3.0)
            .execute(&mut store)
            .unwrap();

        let volume = Volume::new(solid).execute(&store).unwrap();
        // Our sphere is a cone approximation (two cones back-to-back):
        // Each cone: (1/3)*pi*r^2*r = (1/3)*pi*27 → total = 2*(pi*27/3) = 18*pi ≈ 56.55
        let expected = 18.0 * PI;
        let tolerance = expected * 0.10; // looser tolerance for approximated shape
        assert!(
            (volume - expected).abs() < tolerance,
            "expected ~{expected:.2}, got {volume:.2}"
        );
    }

    #[test]
    fn offset_box_volume() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(1.0, 2.0, 3.0), p(3.0, 5.0, 7.0))
            .execute(&mut store)
            .unwrap();

        let volume = Volume::new(solid).execute(&store).unwrap();
        // 2 * 3 * 4 = 24
        assert!((volume - 24.0).abs() < 0.1, "expected 24.0, got {volume}");
    }
}
