//! Closed NURBS-faced solid builders for the through-cut boolean (P6/P7).
//!
//! These constructors exist to produce closed solids whose boundary contains
//! genuine NURBS faces, so the through-cut subtract and its demo have real
//! curved input to operate on. They are deliberately minimal.

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{NurbsCurve3D, NurbsSurface};
use crate::geometry::surface::Surface;
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::topology::{FaceId, ShellData, SolidData, SolidId, TopologyStore};

use super::{MakeFace, MakeNurbsFace, MakeWire};

/// A NURBS prism: a closed planar profile extruded along `direction`, capped by
/// two planar polygonal faces sampled from the profile.
///
/// The side is a single `NurbsSurface` (`NurbsSurface::extrude` of the profile)
/// whose `u` parameter runs around the profile (closed for a closed profile) and
/// whose `v` parameter runs along `direction` — exactly the band UV topology the
/// through-cut boolean expects. The bottom cap sits at the profile and the top
/// cap at the profile translated by `direction`; each cap is oriented so its
/// outward normal points away from the prism body.
///
/// The profile must be closed and planar; the caps are planar polygons sampled
/// from it.
pub struct MakeNurbsPrism {
    profile: NurbsCurve3D,
    direction: Vector3,
}

impl MakeNurbsPrism {
    /// Number of profile samples per cap polygon.
    const CAP_SEGMENTS: usize = 64;

    /// Creates a prism extruding `profile` along `direction`.
    #[must_use]
    pub fn new(profile: NurbsCurve3D, direction: Vector3) -> Self {
        Self { profile, direction }
    }

    /// Builds the prism solid in the store.
    ///
    /// # Errors
    ///
    /// Returns an error if `direction` is zero-length, the profile is not closed,
    /// or any underlying surface / face construction fails.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        if self.direction.norm() <= TOLERANCE {
            return Err(
                OperationError::InvalidInput("prism direction must be non-zero".into()).into(),
            );
        }
        if !self.profile.is_endpoint_closed() {
            return Err(OperationError::InvalidInput("prism profile must be closed".into()).into());
        }

        let side_surface = NurbsSurface::extrude(&self.profile, self.direction)?;
        let side = MakeNurbsFace::new(side_surface).execute(store)?;

        // Sample the profile once; the bottom cap uses these points, the top cap
        // the same points translated by `direction`.
        let base = self.sample_profile()?;
        let newell = newell_normal(&base);
        let along = newell.dot(&self.direction);

        // Bottom cap: outward normal opposes `direction`. Top cap: along it.
        let bottom = cap_face(store, &base, along > 0.0)?;
        let top_pts: Vec<Point3> = base.iter().map(|p| p + self.direction).collect();
        let top = cap_face(store, &top_pts, along < 0.0)?;

        Ok(finish_solid(store, vec![side, bottom, top]))
    }

    /// Samples the profile into `CAP_SEGMENTS` distinct points (excludes the
    /// closing duplicate at the domain end).
    fn sample_profile(&self) -> Result<Vec<Point3>> {
        let (t0, t1) = self.profile.parameter_domain();
        let mut pts = Vec::with_capacity(Self::CAP_SEGMENTS);
        for i in 0..Self::CAP_SEGMENTS {
            #[allow(clippy::cast_precision_loss)]
            let frac = i as f64 / Self::CAP_SEGMENTS as f64;
            pts.push(self.profile.point_at(t0 + (t1 - t0) * frac)?);
        }
        Ok(pts)
    }
}

/// A vertical NURBS tube: a circular profile extruded along `+Z`, capped by two
/// planar disks. Delegates to [`MakeNurbsPrism`] with a rational-circle profile.
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
        let profile = NurbsCurve3D::circle(self.center, self.radius, Vector3::z(), Vector3::x())?;
        let axis = Vector3::new(0.0, 0.0, self.height);
        MakeNurbsPrism::new(profile, axis).execute(store)
    }
}

/// Builds a planar polygonal cap face from an ordered ring of coplanar points.
/// When `reverse` is set the point order is flipped so the face's right-hand
/// normal points away from the prism body.
fn cap_face(store: &mut TopologyStore, points: &[Point3], reverse: bool) -> Result<FaceId> {
    let mut pts = points.to_vec();
    if reverse {
        pts.reverse();
    }
    let wire = MakeWire::new(pts, true).execute(store)?;
    MakeFace::new(wire, vec![]).execute(store)
}

/// Newell's method normal of a (possibly non-convex) planar polygon. The
/// magnitude is twice the polygon area; only its direction is used.
fn newell_normal(pts: &[Point3]) -> Vector3 {
    let mut n = Vector3::zeros();
    let m = pts.len();
    for i in 0..m {
        let a = pts[i];
        let b = pts[(i + 1) % m];
        n.x += (a.y - b.y) * (a.z + b.z);
        n.y += (a.z - b.z) * (a.x + b.x);
        n.z += (a.x - b.x) * (a.y + b.y);
    }
    n
}

/// A curved slab: a NURBS sheet (`front`) thickened along `+Z` by `thickness`.
///
/// The front face is a bicubic-ish NURBS patch with a gentle central rise; the
/// back face is the front control net translated down by `thickness` (a
/// control-net offset — an approximation of a true offset surface, exact here
/// because the translation is rigid). Four exact ruled NURBS side faces close
/// the slab: each is `NurbsSurface::extrude` of a front boundary isocurve along
/// `(0, 0, -thickness)`, so every side boundary is geometrically identical to
/// the front/back curved boundaries (no lens-shaped gaps).
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

        let sides = side_faces(store, &front_surface, self.thickness)?;

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

/// Builds the four exact ruled NURBS side faces of the slab.
///
/// The back sheet is the front control net translated by `(0, 0, -thickness)`,
/// so each side wall is exactly the ruled surface between a front boundary
/// isocurve and that same curve translated down by `thickness`. That ruled
/// surface is `NurbsSurface::extrude(front_boundary, (0, 0, -thickness))`
/// (`degree_v = 1`, exact): its `v = 0` boundary isocurve reproduces the front
/// face's corresponding boundary curve and its `v = 1` isocurve the back face's,
/// so the side faces meet the curved front/back faces with no lens-shaped gap.
fn side_faces(
    store: &mut TopologyStore,
    front: &NurbsSurface,
    thickness: f64,
) -> Result<Vec<FaceId>> {
    let extrude_dir = Vector3::new(0.0, 0.0, -thickness);
    // XY center of the slab, used to orient each side normal outward.
    let center = front.point_at(0.5, 0.5)?;

    let mut faces = Vec::with_capacity(4);
    for boundary in front.boundary_curves()? {
        let side_surface = NurbsSurface::extrude(&boundary, extrude_dir)?;
        let face = MakeNurbsFace::new(side_surface.clone()).execute(store)?;
        // Flip the sense on the sides whose natural normal points inward, so
        // every side face's normal points away from the slab body.
        if !side_normal_points_outward(&side_surface, &center)? {
            store.face_mut(face)?.same_sense = false;
        }
        faces.push(face);
    }
    Ok(faces)
}

/// Reports whether the ruled side surface's natural normal (at its midpoint)
/// points away from the slab body, measured radially in XY from `center`.
fn side_normal_points_outward(side: &NurbsSurface, center: &Point3) -> Result<bool> {
    let ((u_min, u_max), (v_min, v_max)) = side.parameter_domain();
    let u = 0.5 * (u_min + u_max);
    let v = 0.5 * (v_min + v_max);
    let mid = side.point_at(u, v)?;
    let normal = Surface::normal(side, u, v)?;
    // Radial outward direction in the XY plane (side walls are vertical-ish).
    let outward = Vector3::new(mid.x - center.x, mid.y - center.y, 0.0);
    Ok(normal.dot(&outward) >= 0.0)
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
    fn prism_builds_side_and_two_caps() {
        let mut store = TopologyStore::new();
        let profile = NurbsCurve3D::rounded_rectangle(
            Point3::new(0.0, 0.0, 0.0),
            Vector3::x(),
            Vector3::y(),
            2.6,
            2.0,
            0.35,
        )
        .unwrap();
        let solid = MakeNurbsPrism::new(profile, Vector3::new(0.0, 0.0, 1.2))
            .execute(&mut store)
            .unwrap();
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        assert_eq!(shell.faces.len(), 3, "prism = side + 2 caps");
        assert!(matches!(
            store.face(shell.faces[0]).unwrap().surface,
            FaceSurface::Nurbs(_)
        ));
    }

    #[test]
    fn prism_tessellates_manifold() {
        let mut store = TopologyStore::new();
        let profile = NurbsCurve3D::rounded_rectangle(
            Point3::new(0.0, 0.0, 0.0),
            Vector3::x(),
            Vector3::y(),
            2.6,
            2.0,
            0.35,
        )
        .unwrap();
        let solid = MakeNurbsPrism::new(profile, Vector3::new(0.0, 0.0, 1.2))
            .execute(&mut store)
            .unwrap();
        assert_manifold(&store, solid);
    }

    #[test]
    fn prism_side_vertices_lie_on_extruded_surface() {
        use crate::geometry::nurbs::InversionOptions;
        let mut store = TopologyStore::new();
        let profile =
            NurbsCurve3D::circle(Point3::new(1.0, 1.0, 0.0), 0.8, Vector3::z(), Vector3::x())
                .unwrap();
        let direction = Vector3::new(0.3, 0.0, 2.0);
        let expected = NurbsSurface::extrude(&profile, direction).unwrap();
        let solid = MakeNurbsPrism::new(profile, direction)
            .execute(&mut store)
            .unwrap();
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        let side = nurbs_of(&store, shell.faces[0]);
        let opts = InversionOptions::default();
        let ((u0, u1), (v0, v1)) = side.parameter_domain();
        for i in 0..=8 {
            for j in 0..=4 {
                let u = u0 + (u1 - u0) * f64::from(i) / 8.0;
                let v = v0 + (v1 - v0) * f64::from(j) / 4.0;
                let p = side.point_at(u, v).unwrap();
                let inv = expected.closest_point(&p, &opts).unwrap();
                assert!(inv.distance < 1e-9, "side vertex off extruded surface");
            }
        }
    }

    #[test]
    fn prism_rejects_open_profile() {
        let mut store = TopologyStore::new();
        let open = NurbsCurve3D::polyline(&[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ])
        .unwrap();
        assert!(MakeNurbsPrism::new(open, Vector3::z())
            .execute(&mut store)
            .is_err());
    }

    #[test]
    fn slab_builds_six_nurbs_faces() {
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
        assert_eq!(
            nurbs_faces, 6,
            "front + back + 4 exact ruled side faces are all NURBS"
        );
    }

    /// Extracts the NURBS surface backing a face, panicking otherwise.
    fn nurbs_of(store: &TopologyStore, face: FaceId) -> NurbsSurface {
        match &store.face(face).unwrap().surface {
            FaceSurface::Nurbs(n) => n.clone(),
            other => panic!("expected a NURBS face, got {other:?}"),
        }
    }

    /// Asserts two 3D NURBS curves coincide geometrically: sampled at 50 shared
    /// parameter fractions, every pair of points is within 1e-9.
    fn assert_curves_coincide(a: &NurbsCurve3D, b: &NurbsCurve3D, label: &str) {
        let (a0, a1) = a.parameter_domain();
        let (b0, b1) = b.parameter_domain();
        for i in 0..=50 {
            let f = f64::from(i) / 50.0;
            let pa = a.point_at(a0 + (a1 - a0) * f).unwrap();
            let pb = b.point_at(b0 + (b1 - b0) * f).unwrap();
            let d = (pa - pb).norm();
            assert!(
                d < 1e-9,
                "{label}: boundary mismatch at f={f}: distance {d}"
            );
        }
    }

    /// The exact ruled side faces close the gaps: each side face's top (`v=0`)
    /// boundary isocurve coincides with the front face's corresponding boundary
    /// curve, and its bottom (`v=1`) isocurve with the back face's — geometric
    /// coincidence, not merely shared corner points. The old planar sides (chord
    /// through the corners) failed this by design.
    #[test]
    fn slab_side_boundaries_coincide_with_faces() {
        let mut store = TopologyStore::new();
        let solid = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        // Faces are stored as [front, back, side_0 .. side_3], where side_i is
        // extruded from front.boundary_curves()[i] (same isocurve order for the
        // back), so side_i pairs with front/back boundary curve i.
        let front = nurbs_of(&store, shell.faces[0]);
        let back = nurbs_of(&store, shell.faces[1]);
        let front_boundaries = front.boundary_curves().unwrap();
        let back_boundaries = back.boundary_curves().unwrap();

        for (i, &face) in shell.faces[2..6].iter().enumerate() {
            let side = nurbs_of(&store, face);
            let ((_, _), (v_min, v_max)) = side.parameter_domain();
            let top = side.isocurve_v(v_min).unwrap();
            let bottom = side.isocurve_v(v_max).unwrap();
            assert_curves_coincide(
                &top,
                &front_boundaries[i],
                &format!("side {i} top vs front"),
            );
            assert_curves_coincide(
                &bottom,
                &back_boundaries[i],
                &format!("side {i} bottom vs back"),
            );
        }
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
