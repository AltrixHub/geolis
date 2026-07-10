//! Faceted loft between two planar profiles.
//!
//! Builds a watertight solid from two closed planar profiles with the
//! same vertex count: planar caps plus triangulated side faces (each
//! side quad is split along a diagonal, so twisted correspondences
//! stay planar per face). The NURBS-surface loft
//! ([`crate::geometry::NurbsSurface::loft`]) remains the smooth
//! counterpart; sealing it into a solid is follow-up kernel work.

use crate::error::{OperationError, Result};
use crate::math::Point3;
use crate::operations::creation::{MakeFace, MakeSolid};
use crate::topology::{OrientedEdge, ShellData, SolidId, TopologyStore};

use super::extrude::{create_closed_wire, create_line_edge, create_loop_edges, newell_normal};

/// Loft two closed planar profiles (equal vertex counts, index-matched
/// correspondence) into a faceted solid.
pub struct MakeLoft {
    bottom: Vec<Point3>,
    top: Vec<Point3>,
}

impl MakeLoft {
    #[must_use]
    pub fn new(bottom: Vec<Point3>, top: Vec<Point3>) -> Self {
        Self { bottom, top }
    }

    /// Executes the loft, creating a closed solid in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error when the profiles have fewer than 3 vertices or
    /// mismatched counts, are degenerate (zero Newell normal), or when
    /// the cap faces cannot be built (non-planar profiles).
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        let n = self.bottom.len();
        if n < 3 || self.top.len() != n {
            return Err(OperationError::InvalidInput(format!(
                "loft profiles need matching vertex counts >= 3, got {} and {}",
                n,
                self.top.len()
            ))
            .into());
        }

        // Orient both profiles so the bottom's Newell normal points away
        // from the top (outward below); reverse BOTH to keep the
        // index correspondence intact.
        let bottom_normal = newell_normal(&self.bottom)?;
        let centroid = |pts: &[Point3]| {
            let mut c = Point3::new(0.0, 0.0, 0.0);
            #[expect(clippy::cast_precision_loss, reason = "profile counts are small")]
            let inv = 1.0 / pts.len() as f64;
            for p in pts {
                c.x += p.x * inv;
                c.y += p.y * inv;
                c.z += p.z * inv;
            }
            c
        };
        let up = centroid(&self.top) - centroid(&self.bottom);
        let should_reverse = bottom_normal.dot(&up) < 0.0;
        let (bottom_points, top_points): (Vec<Point3>, Vec<Point3>) = if should_reverse {
            (
                self.bottom.iter().rev().copied().collect(),
                self.top.iter().rev().copied().collect(),
            )
        } else {
            (self.bottom.clone(), self.top.clone())
        };

        let bottom_verts: Vec<_> = bottom_points
            .iter()
            .map(|p| store.add_vertex(crate::topology::VertexData::new(*p)))
            .collect();
        let top_verts: Vec<_> = top_points
            .iter()
            .map(|p| store.add_vertex(crate::topology::VertexData::new(*p)))
            .collect();

        let bottom_edges = create_loop_edges(store, &bottom_verts, &bottom_points)?;
        let top_edges = create_loop_edges(store, &top_verts, &top_points)?;
        // Vertical edges i → i and diagonal edges bottom[i] → top[j].
        let mut vert_edges = Vec::with_capacity(n);
        let mut diag_edges = Vec::with_capacity(n);
        for i in 0..n {
            let j = (i + 1) % n;
            vert_edges.push(create_line_edge(
                store,
                bottom_verts[i],
                top_verts[i],
                bottom_points[i],
                top_points[i],
            )?);
            diag_edges.push(create_line_edge(
                store,
                bottom_verts[i],
                top_verts[j],
                bottom_points[i],
                top_points[j],
            )?);
        }

        let mut all_faces = Vec::with_capacity(2 * n + 2);

        // Bottom cap: reversed winding (outward below).
        let bottom_wire_edges: Vec<OrientedEdge> = (0..n)
            .rev()
            .map(|i| OrientedEdge::new(bottom_edges[i], false))
            .collect();
        let bottom_wire = create_closed_wire(store, bottom_wire_edges);
        all_faces.push(MakeFace::new(bottom_wire, vec![]).execute(store)?);

        // Top cap: forward winding (outward above).
        let top_wire_edges: Vec<OrientedEdge> = (0..n)
            .map(|i| OrientedEdge::new(top_edges[i], true))
            .collect();
        let top_wire = create_closed_wire(store, top_wire_edges);
        all_faces.push(MakeFace::new(top_wire, vec![]).execute(store)?);

        // Side triangles: quad (b_i, b_j, t_j, t_i) split along b_i → t_j.
        for i in 0..n {
            let j = (i + 1) % n;
            let tri1 = vec![
                OrientedEdge::new(bottom_edges[i], true),
                OrientedEdge::new(vert_edges[j], true),
                OrientedEdge::new(diag_edges[i], false),
            ];
            let wire1 = create_closed_wire(store, tri1);
            all_faces.push(MakeFace::new(wire1, vec![]).execute(store)?);

            let tri2 = vec![
                OrientedEdge::new(diag_edges[i], true),
                OrientedEdge::new(top_edges[i], false),
                OrientedEdge::new(vert_edges[i], false),
            ];
            let wire2 = create_closed_wire(store, tri2);
            all_faces.push(MakeFace::new(wire2, vec![]).execute(store)?);
        }

        let shell_id = store.add_shell(ShellData {
            faces: all_faces,
            is_closed: true,
        });
        MakeSolid::new(shell_id, vec![]).execute(store)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::query::Volume;

    fn square(cx: f64, cy: f64, half: f64, z: f64) -> Vec<Point3> {
        vec![
            Point3::new(cx - half, cy - half, z),
            Point3::new(cx + half, cy - half, z),
            Point3::new(cx + half, cy + half, z),
            Point3::new(cx - half, cy + half, z),
        ]
    }

    #[test]
    fn straight_loft_matches_the_extruded_box() {
        let mut store = TopologyStore::new();
        let solid = MakeLoft::new(square(0.0, 0.0, 2.0, 0.0), square(0.0, 0.0, 2.0, 3.0))
            .execute(&mut store)
            .unwrap();
        let volume = Volume::new(solid).execute(&store).unwrap();
        assert!((volume - 48.0).abs() < 1e-9, "volume = {volume}");
    }

    #[test]
    fn tapered_loft_forms_a_frustum() {
        // Square frustum: bottom 4x4, top 2x2, height 3.
        // V = h/3 * (A1 + A2 + sqrt(A1*A2)) = 1 * (16 + 4 + 8) = 28.
        let mut store = TopologyStore::new();
        let solid = MakeLoft::new(square(0.0, 0.0, 2.0, 0.0), square(0.0, 0.0, 1.0, 3.0))
            .execute(&mut store)
            .unwrap();
        let volume = Volume::new(solid).execute(&store).unwrap();
        assert!((volume - 28.0).abs() < 1e-9, "volume = {volume}");
    }

    #[test]
    fn offset_loft_shears_without_losing_volume() {
        // Sheared prism keeps the base-area * height volume.
        let mut store = TopologyStore::new();
        let solid = MakeLoft::new(square(0.0, 0.0, 2.0, 0.0), square(3.0, 1.0, 2.0, 3.0))
            .execute(&mut store)
            .unwrap();
        let volume = Volume::new(solid).execute(&store).unwrap();
        assert!((volume - 48.0).abs() < 1e-9, "volume = {volume}");
    }

    #[test]
    fn mismatched_profiles_are_rejected() {
        let mut store = TopologyStore::new();
        let mut top = square(0.0, 0.0, 1.0, 3.0);
        top.pop();
        assert!(MakeLoft::new(square(0.0, 0.0, 2.0, 0.0), top)
            .execute(&mut store)
            .is_err());
    }
}
