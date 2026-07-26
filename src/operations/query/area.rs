use crate::error::{OperationError, Result};
use crate::math::Vector3;
use crate::tessellation::{TessellateSolid, TessellationParams};
use crate::topology::{SolidId, TopologyStore};

/// Computes the surface area of a solid, optionally restricted to the
/// surface facing a given direction.
///
/// Uses tessellation to approximate the area by summing the areas of all
/// triangles in the mesh. The accuracy depends on the tessellation parameters.
pub struct Area {
    solid: SolidId,
    params: TessellationParams,
    facing: Option<Vector3>,
}

impl Area {
    /// Creates a new `Area` query with default tessellation parameters.
    #[must_use]
    pub fn new(solid: SolidId) -> Self {
        Self {
            solid,
            params: TessellationParams::default(),
            facing: None,
        }
    }

    /// Sets custom tessellation parameters for higher accuracy.
    #[must_use]
    pub fn with_params(mut self, params: TessellationParams) -> Self {
        self.params = params;
        self
    }

    /// Restricts the sum to the surface whose outward normal points into
    /// the same open hemisphere as `direction` — the measurement behind
    /// questions like "how much roofing does this solid need" (`+Z`) or
    /// "how much of it is glazed on this elevation".
    ///
    /// Surface exactly perpendicular to `direction` contributes nothing, so
    /// the up-facing and down-facing areas of a closed solid never
    /// double-count its vertical sides. Only the direction matters, not its
    /// magnitude.
    #[must_use]
    pub fn with_facing(mut self, direction: Vector3) -> Self {
        self.facing = Some(direction);
        self
    }

    /// Executes the query, returning the surface area (total, or only the
    /// part facing the [`Area::with_facing`] direction).
    ///
    /// # Errors
    ///
    /// Returns an error if the solid cannot be tessellated, or if a facing
    /// direction was given that is non-finite or degenerate (zero length) —
    /// such a direction selects no hemisphere at all, and silently
    /// answering `0.0` would read as a solid with no surface.
    pub fn execute(&self, store: &TopologyStore) -> Result<f64> {
        if let Some(direction) = self.facing {
            let norm = direction.norm();
            if !norm.is_finite() || norm < crate::math::TOLERANCE {
                return Err(OperationError::InvalidInput(format!(
                    "area facing direction must be finite and non-degenerate, got {direction:?}"
                ))
                .into());
            }
        }
        let mesh = TessellateSolid::new(self.solid, self.params).execute(store)?;

        let mut total_area = 0.0;
        for tri in &mesh.indices {
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];

            let edge1 = v1 - v0;
            let edge2 = v2 - v0;
            // `cross` is the outward normal scaled by twice the triangle
            // area, so one product answers both the facing test and the
            // area contribution.
            let cross = edge1.cross(&edge2);
            if self
                .facing
                .is_some_and(|direction| cross.dot(&direction) <= 0.0)
            {
                continue;
            }
            total_area += cross.norm() * 0.5;
        }

        Ok(total_area)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::{Point3, Vector3};
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

    /// A box's up-facing surface is exactly its top face: the four sides
    /// are perpendicular to `+Z` and the bottom faces away.
    #[test]
    fn facing_up_selects_only_the_top_face() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(2.0, 3.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let up = Area::new(solid)
            .with_facing(Vector3::z())
            .execute(&store)
            .unwrap();
        assert!((up - 6.0).abs() < 1e-9, "top face is 2*3 = 6, got {up}");
        let down = Area::new(solid)
            .with_facing(-Vector3::z())
            .execute(&store)
            .unwrap();
        assert!((down - 6.0).abs() < 1e-9, "bottom face is 6, got {down}");
        let side = Area::new(solid)
            .with_facing(Vector3::x())
            .execute(&store)
            .unwrap();
        assert!(
            (side - 12.0).abs() < 1e-9,
            "+X face is 3*4 = 12, got {side}"
        );
    }

    /// Only the direction matters — scaling it must not scale the answer.
    #[test]
    fn facing_direction_magnitude_is_irrelevant() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(2.0, 3.0, 4.0))
            .execute(&mut store)
            .unwrap();
        let unit = Area::new(solid)
            .with_facing(Vector3::z())
            .execute(&store)
            .unwrap();
        let scaled = Area::new(solid)
            .with_facing(Vector3::new(0.0, 0.0, 100.0))
            .execute(&store)
            .unwrap();
        assert!((unit - scaled).abs() < 1e-12);
    }

    /// The facing partition is exhaustive for an axis-aligned box with no
    /// surface perpendicular to the probe direction.
    #[test]
    fn opposite_facings_partition_the_surface() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(2.0, 3.0, 4.0))
            .execute(&mut store)
            .unwrap();
        let diagonal = Vector3::new(1.0, 1.0, 1.0);
        let front = Area::new(solid)
            .with_facing(diagonal)
            .execute(&store)
            .unwrap();
        let back = Area::new(solid)
            .with_facing(-diagonal)
            .execute(&store)
            .unwrap();
        let total = Area::new(solid).execute(&store).unwrap();
        assert!(
            (front + back - total).abs() < 1e-9,
            "{front} + {back} != {total}"
        );
    }

    /// A degenerate facing direction selects no hemisphere; answering
    /// `0.0` would read as a solid with no surface, so it must error.
    #[test]
    fn rejects_degenerate_facing_direction() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(1.0, 1.0, 1.0))
            .execute(&mut store)
            .unwrap();
        for bad in [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(f64::NAN, 0.0, 1.0),
            Vector3::new(f64::INFINITY, 0.0, 0.0),
        ] {
            assert!(
                Area::new(solid).with_facing(bad).execute(&store).is_err(),
                "facing {bad:?} must be rejected"
            );
        }
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
        assert!(
            area > 50.0 && area < 200.0,
            "unexpected sphere area: {area}"
        );
    }
}
