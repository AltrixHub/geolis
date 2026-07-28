//! Tangent-plane angle between the two faces that share an edge — the
//! "is this edge a visible crease?" query.
//!
//! A wireframe must draw the edges a viewer perceives as creases and skip
//! the ones that are only a topological seam. Topology alone cannot tell
//! them apart: a solid is free to split one smooth surface into several
//! faces (an extruded full circle becomes several cylindrical patches
//! because a segmented prism needs at least three profile segments; a
//! union inserts a split vertex on a flat side). Those internal boundaries
//! carry no geometric discontinuity — the surface normal is the SAME on
//! both sides — while a real corner turns the normal.
//!
//! [`SurfaceAngleAcrossEdge`] answers that with geometry: it evaluates
//! both faces' surface normals at one point of the shared edge and returns
//! the unsigned angle between the tangent planes. `0` means the faces meet
//! smoothly (no crease); `π/2` means they meet at a right angle. A caller
//! filters on its own threshold.
//!
//! Face orientation (`same_sense`, ring winding) is deliberately NOT
//! applied: a crease is a property of the surfaces, not of which side the
//! material is on, so the answer is `acos(|n₀ · n₁|)` and lands in
//! `[0, π/2]`.
//!
//! # Where the parameters come from
//!
//! | Face surface | UV of the sample point |
//! |---|---|
//! | `Plane` | not needed — the normal is constant |
//! | any, with a pcurve for the edge | the pcurve at the edge's mid parameter (exact, same-parameter convention) |
//! | any, without a pcurve | [`ClosestPointOnSurface`] of the edge's mid point |
//!
//! The pcurve path is exact and closed-form, so the common case (solids
//! built with shared-edge topology) costs no surface inversion.

use crate::error::{OperationError, Result};
use crate::geometry::surface::Surface;
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::topology::{EdgeId, FaceId, FaceSurface, TopologyStore};

use super::{ClosestPointOnSurface, PointOnCurve};

/// The unsigned angle (radians, in `[0, π/2]`) between the tangent planes
/// of two faces at the edge they share.
pub struct SurfaceAngleAcrossEdge {
    edge: EdgeId,
    face_a: FaceId,
    face_b: FaceId,
}

impl SurfaceAngleAcrossEdge {
    /// Creates the query for `edge` shared by `face_a` and `face_b`.
    #[must_use]
    pub fn new(edge: EdgeId, face_a: FaceId, face_b: FaceId) -> Self {
        Self {
            edge,
            face_a,
            face_b,
        }
    }

    /// Executes the query.
    ///
    /// # Errors
    ///
    /// Returns an error when the edge or either face is missing, when a
    /// face's surface cannot be evaluated at the shared edge, or when
    /// either normal is degenerate.
    pub fn execute(&self, store: &TopologyStore) -> Result<f64> {
        let edge = store.edge(self.edge)?;
        let t = 0.5 * (edge.t_start + edge.t_end);
        let normal_a = face_normal_at_edge(store, self.face_a, self.edge, t)?;
        let normal_b = face_normal_at_edge(store, self.face_b, self.edge, t)?;
        let cos = normal_a.dot(&normal_b).abs().clamp(0.0, 1.0);
        Ok(cos.acos())
    }
}

/// The unit surface normal of `face` at parameter `t` along its boundary
/// `edge`.
fn face_normal_at_edge(
    store: &TopologyStore,
    face: FaceId,
    edge: EdgeId,
    t: f64,
) -> Result<Vector3> {
    let face_data = store.face(face)?;
    // A plane's normal is constant, so no sample point is needed at all —
    // which also covers the planar faces built before shared-edge topology
    // (no pcurves recorded).
    if let FaceSurface::Plane(plane) = &face_data.surface {
        return unit(*plane.plane_normal());
    }
    // Same-parameter convention: a recorded pcurve is parameterized
    // identically to the edge curve, so `t` maps straight through — exact
    // and closed-form. Without one, invert the edge's mid point onto the
    // surface instead.
    let (u, v) = if let Some(pcurve) = face_data.pcurve_for(edge) {
        let uv = pcurve.point_at(t)?;
        (uv.x, uv.y)
    } else {
        let point: Point3 = PointOnCurve::new(edge, t).execute(store)?;
        let found = ClosestPointOnSurface::new(face, point).execute(store)?;
        (found.u, found.v)
    };
    surface_normal(&face_data.surface, u, v)
}

/// The surface normal at `(u, v)`, for every surface class.
fn surface_normal(surface: &FaceSurface, u: f64, v: f64) -> Result<Vector3> {
    match surface {
        FaceSurface::Plane(s) => unit(*s.plane_normal()),
        FaceSurface::Cylinder(s) => unit(s.normal(u, v)?),
        FaceSurface::Cone(s) => unit(s.normal(u, v)?),
        FaceSurface::Sphere(s) => unit(s.normal(u, v)?),
        FaceSurface::Torus(s) => unit(s.normal(u, v)?),
        FaceSurface::Nurbs(s) => {
            let (_, du, dv) = s.partials(u, v)?;
            unit(du.cross(&dv))
        }
    }
}

/// Normalizes `v`, reporting a degenerate (zero-length) normal.
fn unit(v: Vector3) -> Result<Vector3> {
    let norm = v.norm();
    if !norm.is_finite() || norm < TOLERANCE {
        return Err(OperationError::InvalidInput(format!(
            "surface normal is degenerate (norm = {norm})"
        ))
        .into());
    }
    Ok(v / norm)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Vector3 as V3;
    use crate::operations::creation::{MakeFace, MakeSegmentedPrism, MakeWire, ProfileSegment};
    use crate::operations::shaping::Extrude;
    use crate::topology::{FaceId, ShellId};
    use std::collections::HashMap;
    use std::f64::consts::{FRAC_PI_2, TAU};

    /// Every edge of `shell` that exactly two faces share, with that pair.
    fn shared_edges(store: &TopologyStore, shell: ShellId) -> Vec<(EdgeId, FaceId, FaceId)> {
        let mut by_edge: HashMap<EdgeId, Vec<FaceId>> = HashMap::new();
        for &face_id in &store.shell(shell).unwrap().faces {
            let face = store.face(face_id).unwrap();
            let wires = std::iter::once(face.outer_wire).chain(face.inner_wires.iter().copied());
            for wire_id in wires {
                for oe in &store.wire(wire_id).unwrap().edges {
                    by_edge.entry(oe.edge).or_default().push(face_id);
                }
            }
        }
        by_edge
            .into_iter()
            .filter(|(_, faces)| faces.len() == 2)
            .map(|(edge, faces)| (edge, faces[0], faces[1]))
            .collect()
    }

    /// A full circle of radius `r` in the z = 0 plane, as three arc
    /// segments (the minimum a segmented prism accepts).
    fn circle_profile(r: f64) -> Vec<ProfileSegment> {
        let center = Point3::new(0.0, 0.0, 0.0);
        (0..3)
            .map(|k| {
                let a0 = f64::from(k) * TAU / 3.0;
                ProfileSegment::Arc {
                    center,
                    radius: r,
                    normal: V3::z(),
                    ref_dir: V3::new(a0.cos(), a0.sin(), 0.0),
                    start_angle: 0.0,
                    end_angle: TAU / 3.0,
                }
            })
            .collect()
    }

    /// The three lateral patches of an extruded circle are ONE cylinder:
    /// every vertical seam between them reads as a zero angle, while the
    /// rim edges against the planar caps read as a right angle.
    #[test]
    fn extruded_circle_seams_are_smooth_and_rims_are_right_angles() {
        let mut store = TopologyStore::new();
        let solid = MakeSegmentedPrism::new(circle_profile(0.2), V3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();
        let shell = store.solid(solid).unwrap().outer_shell;

        let mut smooth = 0;
        let mut right = 0;
        for (edge, a, b) in shared_edges(&store, shell) {
            let angle = SurfaceAngleAcrossEdge::new(edge, a, b)
                .execute(&store)
                .unwrap();
            if angle < 1e-6 {
                smooth += 1;
            } else {
                assert!(
                    (angle - FRAC_PI_2).abs() < 1e-9,
                    "a rim edge must be a right angle; got {angle}"
                );
                right += 1;
            }
        }
        assert_eq!(smooth, 3, "one smooth seam per lateral patch joint");
        assert_eq!(right, 6, "three bottom + three top rim edges");
    }

    /// A chord-faceted circle is NOT smooth: each facet joint turns the
    /// normal by the facet angle, so faceting can never be mistaken for a
    /// tangent-continuous seam.
    #[test]
    fn chord_facet_joints_report_the_facet_angle() {
        const N: usize = 24;
        let mut store = TopologyStore::new();
        let points: Vec<Point3> = (0..N)
            .map(|k| {
                #[allow(clippy::cast_precision_loss)]
                let a = k as f64 * TAU / N as f64;
                Point3::new(0.2 * a.cos(), 0.2 * a.sin(), 0.0)
            })
            .collect();
        let wire = MakeWire::new(points, true).execute(&mut store).unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, V3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();
        let shell = store.solid(solid).unwrap().outer_shell;

        #[allow(clippy::cast_precision_loss)]
        let facet_angle = TAU / N as f64;
        let mut facets = 0;
        for (edge, a, b) in shared_edges(&store, shell) {
            let angle = SurfaceAngleAcrossEdge::new(edge, a, b)
                .execute(&store)
                .unwrap();
            if (angle - facet_angle).abs() < 1e-9 {
                facets += 1;
            }
        }
        assert_eq!(
            facets, N,
            "every one of the {N} facet joints must report the {facet_angle} rad turn",
        );
    }

    /// A box: every edge is a right angle, planar face pair included.
    #[test]
    fn box_edges_are_right_angles() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, V3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();
        let shell = store.solid(solid).unwrap().outer_shell;

        let edges = shared_edges(&store, shell);
        assert_eq!(edges.len(), 12, "a box has 12 shared edges");
        for (edge, a, b) in edges {
            let angle = SurfaceAngleAcrossEdge::new(edge, a, b)
                .execute(&store)
                .unwrap();
            assert!(
                (angle - FRAC_PI_2).abs() < 1e-9,
                "box edge angle must be π/2; got {angle}"
            );
        }
    }
}
