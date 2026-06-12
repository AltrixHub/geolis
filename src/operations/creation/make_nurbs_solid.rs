//! Closed NURBS-faced solid builders for the through-cut boolean (P6/P7).
//!
//! These constructors exist to produce closed solids whose boundary contains
//! genuine NURBS faces, so the through-cut subtract and its demo have real
//! curved input to operate on. They are deliberately minimal.

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{NurbsCurve3D, NurbsSurface};
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::topology::{FaceId, ShellData, SolidData, SolidId, TopologyStore};

use super::{MakeFace, MakeNurbsFace, MakeWire};

/// A vertical NURBS tube: a circular profile extruded along `+Z`, capped by two
/// planar disks.
///
/// The side is a single closed `NurbsSurface` (`NurbsSurface::extrude` of a
/// rational circle) whose `u` parameter runs around the tube (closed) and whose
/// `v` parameter runs along the axis — exactly the band UV topology the
/// through-cut subtract requires.
pub struct MakeNurbsTube {
    center: Point3,
    radius: f64,
    height: f64,
}

impl MakeNurbsTube {
    /// Creates a tube of `radius` rising `height` along `+Z` from `center`
    /// (the center of the bottom cap circle).
    #[must_use]
    pub fn new(center: Point3, radius: f64, height: f64) -> Self {
        Self {
            center,
            radius,
            height,
        }
    }

    /// Builds the tube solid in the store.
    ///
    /// # Errors
    ///
    /// Returns an error if `radius` or `height` is non-positive, or if any
    /// underlying surface / face construction fails.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        if self.radius <= TOLERANCE || self.height <= TOLERANCE {
            return Err(OperationError::InvalidInput(
                "tube radius and height must be positive".into(),
            )
            .into());
        }

        let axis = Vector3::new(0.0, 0.0, self.height);
        let bottom_circle =
            NurbsCurve3D::circle(self.center, self.radius, Vector3::z(), Vector3::x())?;
        let side_surface = NurbsSurface::extrude(&bottom_circle, axis)?;
        let side = MakeNurbsFace::new(side_surface).execute(store)?;

        let top_center = self.center + axis;
        let bottom = self.cap_face(store, self.center, false)?;
        let top = self.cap_face(store, top_center, true)?;

        Ok(finish_solid(store, vec![side, bottom, top]))
    }

    /// Builds a planar circular cap as a polygonal disk. `upward` orients the
    /// cap so its outward normal points away from the tube body.
    fn cap_face(&self, store: &mut TopologyStore, center: Point3, upward: bool) -> Result<FaceId> {
        const CAP_SEGMENTS: usize = 48;
        let mut pts = Vec::with_capacity(CAP_SEGMENTS);
        for i in 0..CAP_SEGMENTS {
            #[allow(clippy::cast_precision_loss)]
            let angle = std::f64::consts::TAU * (i as f64) / (CAP_SEGMENTS as f64);
            // Bottom cap winds clockwise (normal -Z), top cap counter-clockwise
            // (normal +Z) when the wire is built in this order.
            let a = if upward { angle } else { -angle };
            pts.push(Point3::new(
                center.x + self.radius * a.cos(),
                center.y + self.radius * a.sin(),
                center.z,
            ));
        }
        let wire = MakeWire::new(pts, true).execute(store)?;
        MakeFace::new(wire, vec![]).execute(store)
    }
}

/// A curved slab: a NURBS sheet (`front`) thickened along `+Z` by `thickness`.
///
/// The front face is a bicubic-ish NURBS patch with a gentle central rise; the
/// back face is the front control net translated down by `thickness` (a
/// control-net offset — an approximation of a true offset surface, exact here
/// because the translation is rigid). Four planar side faces close the slab
/// from the front/back boundary curves.
///
/// The slab spans `[0, size] x [0, size]` in XY; the front patch sits at
/// `z = base_z + bulge * shape(u,v)` and the back patch `thickness` below it.
pub struct MakeCurvedSlab {
    size: f64,
    base_z: f64,
    bulge: f64,
    thickness: f64,
}

impl MakeCurvedSlab {
    /// Creates a curved slab spanning `[0,size]^2` in XY, front patch peaking
    /// `bulge` above `base_z`, thickened `thickness` downward.
    #[must_use]
    pub fn new(size: f64, base_z: f64, bulge: f64, thickness: f64) -> Self {
        Self {
            size,
            base_z,
            bulge,
            thickness,
        }
    }

    /// Builds the slab solid in the store.
    ///
    /// # Errors
    ///
    /// Returns an error if `size` or `thickness` is non-positive, or any face
    /// construction fails.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        if self.size <= TOLERANCE || self.thickness <= TOLERANCE {
            return Err(OperationError::InvalidInput(
                "slab size and thickness must be positive".into(),
            )
            .into());
        }

        let front_surface = self.sheet_surface(0.0);
        let back_surface = self.sheet_surface(-self.thickness);
        let front_surface = front_surface?;
        let back_surface = back_surface?;

        // same_sense defaults to true (normal = +Z-ish, outward for the top).
        let front = MakeNurbsFace::new(front_surface.clone()).execute(store)?;
        // Back face normal must point downward (outward). Flip same_sense.
        let back = MakeNurbsFace::new(back_surface.clone()).execute(store)?;
        store.face_mut(back)?.same_sense = false;

        let sides = side_faces(store, &front_surface, &back_surface)?;

        let mut faces = vec![front, back];
        faces.extend(sides);
        Ok(finish_solid(store, faces))
    }

    /// Builds the 4×4 control-net NURBS sheet at vertical offset `dz` from the
    /// front patch's nominal height.
    fn sheet_surface(&self, dz: f64) -> Result<NurbsSurface> {
        use crate::geometry::nurbs::KnotVector;
        // Normalized central-rise profile over the 4x4 control net.
        let shape = [
            [0.0_f64, 0.3, 0.3, 0.0],
            [0.3, 1.0, 1.0, 0.3],
            [0.3, 1.0, 1.0, 0.3],
            [0.0, 0.3, 0.3, 0.0],
        ];
        let step = self.size / 3.0;
        let mut control = Vec::with_capacity(16);
        for (i, row) in shape.iter().enumerate() {
            for (j, &s) in row.iter().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let x = i as f64 * step;
                #[allow(clippy::cast_precision_loss)]
                let y = j as f64 * step;
                let z = self.base_z + self.bulge * s + dz;
                control.push(Point3::new(x, y, z));
            }
        }
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])?;
        NurbsSurface::from_unweighted(control, 4, 4, knots.clone(), knots, 3, 3)
    }
}

/// Builds the four planar side faces from the front/back boundary curves,
/// sampling each boundary into a polyline ribbon and closing it as a quad strip
/// face (planar per side because both boundaries share an XY edge).
fn side_faces(
    store: &mut TopologyStore,
    front: &NurbsSurface,
    back: &NurbsSurface,
) -> Result<Vec<FaceId>> {
    const N: usize = 16;
    let ((u0, u1), (v0, v1)) = front.parameter_domain();
    // Four domain edges as (u,v) parameter walks. Order chosen so each ribbon's
    // front-edge + back-edge + verticals wind consistently.
    let edges: [[(f64, f64); 2]; 4] = [
        [(u0, v0), (u1, v0)], // y = 0 side
        [(u1, v0), (u1, v1)], // x = size side
        [(u1, v1), (u0, v1)], // y = size side
        [(u0, v1), (u0, v0)], // x = 0 side
    ];

    let mut faces = Vec::with_capacity(4);
    for edge in edges {
        let face = side_ribbon(store, front, back, edge, N)?;
        faces.push(face);
    }
    Ok(faces)
}

/// Builds one planar side face. The front and back patches share the same XY
/// boundary (the thickening is a pure vertical translation of the control net),
/// so each side lies in a vertical plane and is a valid planar quad strip: front
/// edge forward, then down to the back edge, then back edge reversed.
fn side_ribbon(
    store: &mut TopologyStore,
    front: &NurbsSurface,
    back: &NurbsSurface,
    edge: [(f64, f64); 2],
    n: usize,
) -> Result<FaceId> {
    let [start, end] = edge;
    let mut ring: Vec<Point3> = Vec::with_capacity(2 * (n + 1));
    // Front edge start -> end.
    for i in 0..=n {
        #[allow(clippy::cast_precision_loss)]
        let t = i as f64 / n as f64;
        let u = start.0 + (end.0 - start.0) * t;
        let v = start.1 + (end.1 - start.1) * t;
        ring.push(front.point_at(u, v)?);
    }
    // Back edge end -> start.
    for i in 0..=n {
        #[allow(clippy::cast_precision_loss)]
        let t = i as f64 / n as f64;
        let u = end.0 + (start.0 - end.0) * t;
        let v = end.1 + (start.1 - end.1) * t;
        ring.push(back.point_at(u, v)?);
    }
    // Drop duplicate seam points so MakeWire does not see coincident pairs.
    dedup_ring(&mut ring);
    let wire = MakeWire::new(ring, true).execute(store)?;
    MakeFace::new(wire, vec![]).execute(store)
}

/// Removes consecutive coincident points and a coincident wrap-around point.
fn dedup_ring(pts: &mut Vec<Point3>) {
    pts.dedup_by(|a, b| (*a - *b).norm() < TOLERANCE);
    while pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).norm() < TOLERANCE {
        pts.pop();
    }
}

/// A closed solid of revolution: a planar profile revolved a full turn about an
/// axis-parallel line, closed by two planar caps at the profile ends.
///
/// The profile is given as `(radius, z)` samples; it is revolved about the `+Z`
/// axis through the origin. The revolved wall is a single closed NURBS surface
/// (`NurbsSurface::revolve`); the two caps are planar annular/disk faces at the
/// first and last profile heights.
pub struct MakeRevolvedSolid {
    profile: Vec<(f64, f64)>,
}

impl MakeRevolvedSolid {
    /// Creates a revolved solid from `(radius, z)` profile samples (at least 2).
    #[must_use]
    pub fn new(profile: Vec<(f64, f64)>) -> Self {
        Self { profile }
    }

    /// Builds the revolved solid in the store.
    ///
    /// # Errors
    ///
    /// Returns an error if fewer than 2 profile samples are given, any radius is
    /// non-positive, or surface / face construction fails.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        if self.profile.len() < 2 {
            return Err(OperationError::InvalidInput(
                "revolved solid needs at least 2 profile samples".into(),
            )
            .into());
        }
        if self.profile.iter().any(|&(r, _)| r <= TOLERANCE) {
            return Err(OperationError::InvalidInput(
                "revolved profile radii must be positive".into(),
            )
            .into());
        }

        let profile_pts: Vec<Point3> = self
            .profile
            .iter()
            .map(|&(r, z)| Point3::new(r, 0.0, z))
            .collect();
        let (profile_curve, _) =
            NurbsCurve3D::interpolate(&profile_pts, 3.min(profile_pts.len() - 1))?;
        let wall_surface = NurbsSurface::revolve(
            &profile_curve,
            Point3::origin(),
            Vector3::z(),
            std::f64::consts::TAU,
        )?;
        let wall = MakeNurbsFace::new(wall_surface).execute(store)?;

        let (r0, z0) = self.profile[0];
        let (r1, z1) = self.profile[self.profile.len() - 1];
        let bottom = disk_face(store, Point3::new(0.0, 0.0, z0), r0, false)?;
        let top = disk_face(store, Point3::new(0.0, 0.0, z1), r1, true)?;

        Ok(finish_solid(store, vec![wall, bottom, top]))
    }
}

/// Builds a planar polygonal disk face centered at `center` with `radius`.
/// `upward` orients the disk so its normal points away from the body.
fn disk_face(
    store: &mut TopologyStore,
    center: Point3,
    radius: f64,
    upward: bool,
) -> Result<FaceId> {
    const SEGMENTS: usize = 48;
    let mut pts = Vec::with_capacity(SEGMENTS);
    for i in 0..SEGMENTS {
        #[allow(clippy::cast_precision_loss)]
        let angle = std::f64::consts::TAU * (i as f64) / (SEGMENTS as f64);
        let a = if upward { angle } else { -angle };
        pts.push(Point3::new(
            center.x + radius * a.cos(),
            center.y + radius * a.sin(),
            center.z,
        ));
    }
    let wire = MakeWire::new(pts, true).execute(store)?;
    MakeFace::new(wire, vec![]).execute(store)
}

/// Wraps a face list into a closed shell + solid.
pub(crate) fn finish_solid(store: &mut TopologyStore, faces: Vec<FaceId>) -> SolidId {
    let shell = store.add_shell(ShellData {
        faces,
        is_closed: true,
    });
    store.add_solid(SolidData {
        outer_shell: shell,
        inner_shells: vec![],
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::tessellation::{TessellateSolid, TessellationParams};
    use crate::topology::FaceSurface;
    use std::collections::HashMap;

    /// Asserts the tessellated solid mesh is edge-manifold: every undirected
    /// triangle edge is shared by 1 or 2 triangles (P5 manifold pattern).
    fn assert_manifold(store: &TopologyStore, solid: SolidId) {
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(store)
            .unwrap();
        assert!(!mesh.indices.is_empty(), "empty mesh");
        let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
        for tri in &mesh.indices {
            for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }
        for (&(a, b), &c) in &counts {
            assert!(c == 1 || c == 2, "edge ({a},{b}) used {c} times");
        }
    }

    #[test]
    fn tube_builds_three_faces() {
        let mut store = TopologyStore::new();
        let solid = MakeNurbsTube::new(Point3::new(2.0, 2.0, -1.0), 0.8, 4.0)
            .execute(&mut store)
            .unwrap();
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        assert_eq!(shell.faces.len(), 3, "tube = side + 2 caps");
        // Side face is a NURBS surface.
        let side = store.face(shell.faces[0]).unwrap();
        assert!(matches!(side.surface, FaceSurface::Nurbs(_)));
    }

    #[test]
    fn tube_tessellates_manifold() {
        let mut store = TopologyStore::new();
        let solid = MakeNurbsTube::new(Point3::new(0.0, 0.0, 0.0), 1.0, 3.0)
            .execute(&mut store)
            .unwrap();
        assert_manifold(&store, solid);
    }

    #[test]
    fn slab_builds_six_faces_two_nurbs() {
        let mut store = TopologyStore::new();
        let solid = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        assert_eq!(shell.faces.len(), 6, "slab = front + back + 4 sides");
        let nurbs_faces = shell
            .faces
            .iter()
            .filter(|&&f| matches!(store.face(f).unwrap().surface, FaceSurface::Nurbs(_)))
            .count();
        assert_eq!(nurbs_faces, 2, "front + back are NURBS");
    }

    #[test]
    fn slab_tessellates_manifold() {
        let mut store = TopologyStore::new();
        let solid = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        assert_manifold(&store, solid);
    }

    #[test]
    fn revolved_solid_builds_and_tessellates() {
        let mut store = TopologyStore::new();
        let profile = vec![(1.0, 0.0), (1.6, 1.0), (1.2, 2.0), (1.8, 3.0)];
        let solid = MakeRevolvedSolid::new(profile).execute(&mut store).unwrap();
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        assert_eq!(shell.faces.len(), 3, "wall + 2 caps");
        assert_manifold(&store, solid);
    }

    #[test]
    fn rejects_degenerate_inputs() {
        let mut store = TopologyStore::new();
        assert!(MakeNurbsTube::new(Point3::origin(), 0.0, 1.0)
            .execute(&mut store)
            .is_err());
        assert!(MakeCurvedSlab::new(1.0, 0.0, 0.0, 0.0)
            .execute(&mut store)
            .is_err());
        assert!(MakeRevolvedSolid::new(vec![(1.0, 0.0)])
            .execute(&mut store)
            .is_err());
    }
}
