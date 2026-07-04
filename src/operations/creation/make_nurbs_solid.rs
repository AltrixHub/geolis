//! Closed NURBS-faced solid builders for the through-cut boolean (P6/P7).
//!
//! These constructors exist to produce closed solids whose boundary contains
//! genuine NURBS faces, so the through-cut subtract and its demo have real
//! curved input to operate on. They are deliberately minimal.

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{KnotVector, NurbsCurve2D, NurbsCurve3D, NurbsSurface};
use crate::geometry::surface::Surface;
use crate::math::{Point2, Point3, Vector3, TOLERANCE};
use crate::topology::{
    EdgeCurve, EdgeData, EdgeId, EdgeName, EdgeRole, FaceId, FaceName, FacePcurve, FaceRole,
    FaceSurface, OpId, OrientedEdge, ShellData, SolidData, SolidId, TopologyStore, VertexData,
    WireData,
};

use super::{MakeFace, MakeNurbsFace};

/// A NURBS prism: a closed planar profile extruded along `direction`, capped by
/// two exact planar disks on the profile's boundary curves.
///
/// The side is a single `NurbsSurface` (`NurbsSurface::extrude` of the profile)
/// whose `u` parameter runs around the profile (closed for a closed profile) and
/// whose `v` parameter runs along `direction` — exactly the band UV topology the
/// through-cut boolean expects. The bottom and top ring curves are **shared
/// edges**: the side face's wire and each cap's wire reference the same
/// [`EdgeId`], so the solid tessellates boundary-conformally by construction
/// (both faces consume the same per-edge samples).
pub struct MakeNurbsPrism {
    profile: NurbsCurve3D,
    direction: Vector3,
    op_id: Option<OpId>,
}

impl MakeNurbsPrism {
    /// Creates a prism extruding `profile` along `direction`.
    #[must_use]
    pub fn new(profile: NurbsCurve3D, direction: Vector3) -> Self {
        Self {
            profile,
            direction,
            op_id: None,
        }
    }

    /// Registers persistent names for the prism's faces and ring edges under
    /// the caller-supplied operation identity.
    #[must_use]
    pub fn with_op_id(mut self, op: OpId) -> Self {
        self.op_id = Some(op);
        self
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
        let ((u0, u1), (v0, v1)) = side_surface.parameter_domain();

        // Shared ring edges: the side surface's exact v-boundary isocurves (the
        // bottom one IS the profile, with identical knots), each a closed edge
        // with a single shared vertex.
        let bottom_edge = closed_ring_edge(store, side_surface.isocurve_v(v0)?)?;
        let top_edge = closed_ring_edge(store, side_surface.isocurve_v(v1)?)?;

        // Side face: one wire referencing both rings; the pcurves map the ring
        // parameter straight onto the u axis at the fixed v (same-parameter by
        // construction: extrude preserves the profile knots in u).
        let side_wire = store.add_wire(WireData {
            edges: vec![
                OrientedEdge::new(bottom_edge, true),
                OrientedEdge::new(top_edge, false),
            ],
            is_closed: true,
        });
        let pcurves = vec![
            FacePcurve {
                edge: bottom_edge,
                curve: iso_pcurve_u(u0, u1, v0)?,
            },
            FacePcurve {
                edge: top_edge,
                curve: iso_pcurve_u(u0, u1, v1)?,
            },
        ];
        let side = MakeNurbsFace::new(side_surface)
            .with_boundary(side_wire, pcurves)
            .execute(store)?;

        // Caps: planar faces whose outer wire is the SAME ring edge. Bottom
        // outward normal opposes `direction`; top points along it.
        let bottom = cap_from_ring(store, bottom_edge, -self.direction)?;
        let top = cap_from_ring(store, top_edge, self.direction)?;

        if let Some(op) = &self.op_id {
            bind_created_face(store, side, op, FaceRole::Side(0));
            bind_created_face(store, bottom, op, FaceRole::CapStart);
            bind_created_face(store, top, op, FaceRole::CapEnd);
            bind_created_edge(store, bottom_edge, op, EdgeRole::RingStart);
            bind_created_edge(store, top_edge, op, EdgeRole::RingEnd);
        }

        Ok(finish_solid(store, vec![side, bottom, top]))
    }
}

/// Registers a creation-op face name.
pub(crate) fn bind_created_face(
    store: &mut TopologyStore,
    face: FaceId,
    op: &OpId,
    role: FaceRole,
) {
    store.names_mut().bind_face(
        face,
        FaceName::Created {
            op: op.clone(),
            role,
        },
    );
}

/// Registers a creation-op edge name.
fn bind_created_edge(store: &mut TopologyStore, edge: EdgeId, op: &OpId, role: EdgeRole) {
    store.names_mut().bind_edge(
        edge,
        EdgeName::Created {
            op: op.clone(),
            role,
        },
    );
}

/// Adds a closed ring edge (start == end vertex) for `curve`.
fn closed_ring_edge(store: &mut TopologyStore, curve: NurbsCurve3D) -> Result<EdgeId> {
    let (t0, t1) = curve.parameter_domain();
    let start = curve.point_at(t0)?;
    let vertex = store.add_vertex(VertexData { point: start });
    Ok(store.add_edge(EdgeData {
        start: vertex,
        end: vertex,
        curve: EdgeCurve::Nurbs(curve),
        t_start: t0,
        t_end: t1,
    }))
}

/// Degree-1 pcurve mapping a ring edge's parameter onto the `u` axis at a
/// fixed `v` (`t → (t, v)` over `[u0, u1]`).
pub(crate) fn iso_pcurve_u(u0: f64, u1: f64, v: f64) -> Result<NurbsCurve2D> {
    NurbsCurve2D::from_unweighted(
        vec![Point2::new(u0, v), Point2::new(u1, v)],
        KnotVector::new(vec![u0, u0, u1, u1])?,
        1,
    )
}

/// Planar cap face on a shared ring edge, oriented so its stored normal points
/// along `outward`.
fn cap_from_ring(store: &mut TopologyStore, ring: EdgeId, outward: Vector3) -> Result<FaceId> {
    let wire = store.add_wire(WireData {
        edges: vec![OrientedEdge::new(ring, true)],
        is_closed: true,
    });
    let face = MakeFace::new(wire, vec![]).execute(store)?;
    let flip = match &store.face(face)?.surface {
        FaceSurface::Plane(plane) => plane.plane_normal().dot(&outward) < 0.0,
        _ => false,
    };
    if flip {
        store.face_mut(face)?.same_sense = false;
    }
    Ok(face)
}

/// A vertical NURBS tube: a circular profile extruded along `+Z`, capped by two
/// planar disks. Delegates to [`MakeNurbsPrism`] with a rational-circle profile.
pub struct MakeNurbsTube {
    center: Point3,
    radius: f64,
    height: f64,
    op_id: Option<OpId>,
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
            op_id: None,
        }
    }

    /// Registers persistent names under the caller-supplied operation
    /// identity (delegated to the underlying prism).
    #[must_use]
    pub fn with_op_id(mut self, op: OpId) -> Self {
        self.op_id = Some(op);
        self
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
        let mut prism = MakeNurbsPrism::new(profile, axis);
        if let Some(op) = &self.op_id {
            prism = prism.with_op_id(op.clone());
        }
        prism.execute(store)
    }
}

/// A curved (plan-arc) wall: a vertical prism whose plan footprint is a circular
/// annular sector.
///
/// The wall bends along a circular arc in plan (about `arc_center`, mid-surface
/// `radius`, from `start_angle` to `end_angle` measured from `+X` toward `+Y`),
/// rises `height` along `+Z`, and is `thickness` thick radially. Its inner and
/// outer faces are exact concentric-arc extrusions at `radius ∓ thickness / 2`
/// (the concentric arcs are the exact radial offset — no approximation); the
/// top/bottom are exact ruled sectors between the inner and outer arcs and the
/// two ends are exact ruled radial-vertical rectangles. All six faces are NURBS
/// and share their boundary curves exactly, so the solid tessellates
/// boundary-conformally.
pub struct MakeCurvedWall {
    arc_center: Point3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
    height: f64,
    thickness: f64,
    op_id: Option<OpId>,
}

impl MakeCurvedWall {
    /// Creates a curved wall about `arc_center` (base plane at `arc_center.z`).
    #[must_use]
    pub fn new(
        arc_center: Point3,
        radius: f64,
        start_angle: f64,
        end_angle: f64,
        height: f64,
        thickness: f64,
    ) -> Self {
        Self {
            arc_center,
            radius,
            start_angle,
            end_angle,
            height,
            thickness,
            op_id: None,
        }
    }

    /// Registers persistent names for the wall's six faces under the
    /// caller-supplied operation identity (0 = inner, 1 = outer, 2 = start
    /// end, 3 = end end, plus Top / Bottom).
    #[must_use]
    pub fn with_op_id(mut self, op: OpId) -> Self {
        self.op_id = Some(op);
        self
    }

    /// Builds the curved wall solid in the store.
    ///
    /// # Errors
    ///
    /// Returns an error if `height` or `thickness` is non-positive, the inner
    /// radius `radius - thickness / 2` is non-positive, or any arc / surface /
    /// face construction fails (including an out-of-range angular sweep).
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        if self.height <= TOLERANCE || self.thickness <= TOLERANCE {
            return Err(OperationError::InvalidInput(
                "curved wall height and thickness must be positive".into(),
            )
            .into());
        }
        let half = 0.5 * self.thickness;
        if self.radius - half <= TOLERANCE {
            return Err(OperationError::InvalidInput(
                "curved wall inner radius (radius - thickness/2) must be positive".into(),
            )
            .into());
        }

        let normal = Vector3::z();
        let ref_dir = Vector3::x();
        let base = self.arc_center;
        let top = self.arc_center + Vector3::new(0.0, 0.0, self.height);
        let up = Vector3::new(0.0, 0.0, self.height);

        let arc = |c: Point3, r: f64| -> Result<NurbsCurve3D> {
            NurbsCurve3D::arc(c, r, normal, ref_dir, self.start_angle, self.end_angle)
        };
        let inner_base = arc(base, self.radius - half)?;
        let outer_base = arc(base, self.radius + half)?;
        let inner_top = arc(top, self.radius - half)?;
        let outer_top = arc(top, self.radius + half)?;

        // Radial bottom edges at the two angular ends (inner -> outer).
        let (ti0, ti1) = inner_base.parameter_domain();
        let start_radial =
            NurbsCurve3D::polyline(&[inner_base.point_at(ti0)?, outer_base.point_at(ti0)?])?;
        let end_radial =
            NurbsCurve3D::polyline(&[inner_base.point_at(ti1)?, outer_base.point_at(ti1)?])?;

        // Six exact NURBS faces.
        let inner_surf = NurbsSurface::extrude(&inner_base, up)?;
        let outer_surf = NurbsSurface::extrude(&outer_base, up)?;
        let bottom_surf = ruled_surface(&inner_base, &outer_base)?;
        let top_surf = ruled_surface(&inner_top, &outer_top)?;
        let start_surf = NurbsSurface::extrude(&start_radial, up)?;
        let end_surf = NurbsSurface::extrude(&end_radial, up)?;

        // Body centroid, used to orient every face normal outward.
        let mid_angle = 0.5 * (self.start_angle + self.end_angle);
        let radial = Vector3::new(mid_angle.cos(), mid_angle.sin(), 0.0);
        let centroid =
            self.arc_center + radial * self.radius + Vector3::new(0.0, 0.0, 0.5 * self.height);

        let mut faces = Vec::with_capacity(6);
        for surf in [
            inner_surf,
            outer_surf,
            bottom_surf,
            top_surf,
            start_surf,
            end_surf,
        ] {
            faces.push(make_outward_face(store, surf, &centroid)?);
        }

        if let Some(op) = &self.op_id {
            let roles = [
                FaceRole::Side(0),
                FaceRole::Side(1),
                FaceRole::Bottom,
                FaceRole::Top,
                FaceRole::Side(2),
                FaceRole::Side(3),
            ];
            for (&face, role) in faces.iter().zip(roles) {
                bind_created_face(store, face, op, role);
            }
        }

        Ok(finish_solid(store, faces))
    }
}

/// Builds a NURBS face from `surface` and flips its sense when the natural
/// midpoint normal points toward `centroid` (so the stored normal faces
/// outward from the body).
fn make_outward_face(
    store: &mut TopologyStore,
    surface: NurbsSurface,
    centroid: &Point3,
) -> Result<FaceId> {
    let ((u0, u1), (v0, v1)) = surface.parameter_domain();
    let u = 0.5 * (u0 + u1);
    let v = 0.5 * (v0 + v1);
    let mid = surface.point_at(u, v)?;
    let normal = Surface::normal(&surface, u, v)?;
    let face = MakeNurbsFace::new(surface).execute(store)?;
    if normal.dot(&(mid - centroid)) < 0.0 {
        store.face_mut(face)?.same_sense = false;
    }
    Ok(face)
}

/// Builds the exact ruled NURBS surface between two curves that share a degree,
/// knot vector, and per-index weights (row `v = 0` is `a`, row `v = 1` is `b`).
///
/// With matching endpoint weights each `v`-isoline is a straight segment, so for
/// two concentric coplanar arcs the ruled surface is the exact planar annular
/// sector between them.
fn ruled_surface(a: &NurbsCurve3D, b: &NurbsCurve3D) -> Result<NurbsSurface> {
    use crate::geometry::nurbs::KnotVector;
    let nu = a.control_points().len();
    let mut control_points = Vec::with_capacity(nu * 2);
    let mut weights = Vec::with_capacity(nu * 2);
    for i in 0..nu {
        control_points.push(a.control_points()[i]);
        control_points.push(b.control_points()[i]);
        weights.push(a.weights()[i]);
        weights.push(b.weights()[i]);
    }
    let knots_v = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0])?;
    NurbsSurface::new(
        control_points,
        weights,
        nu,
        2,
        a.knots().clone(),
        knots_v,
        a.degree(),
        1,
    )
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
    op_id: Option<OpId>,
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
            op_id: None,
            size,
            base_z,
            bulge,
            thickness,
        }
    }

    /// Registers persistent names for the slab's six faces under the
    /// caller-supplied operation identity (Top = front, Bottom = back,
    /// Side(0..=3) = the ruled side walls in build order).
    #[must_use]
    pub fn with_op_id(mut self, op: OpId) -> Self {
        self.op_id = Some(op);
        self
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

        if let Some(op) = &self.op_id {
            bind_created_face(store, faces[0], op, FaceRole::Top);
            bind_created_face(store, faces[1], op, FaceRole::Bottom);
            for (k, &side) in faces[2..].iter().enumerate() {
                #[allow(clippy::cast_possible_truncation)]
                bind_created_face(store, side, op, FaceRole::Side(k as u8));
            }
        }

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
    op_id: Option<OpId>,
}

impl MakeRevolvedSolid {
    /// Creates a revolved solid from `(radius, z)` profile samples (at least 2).
    #[must_use]
    pub fn new(profile: Vec<(f64, f64)>) -> Self {
        Self {
            profile,
            op_id: None,
        }
    }

    /// Registers persistent names (`Wall`, `CapStart` / `CapEnd`, ring edges)
    /// under the caller-supplied operation identity.
    #[must_use]
    pub fn with_op_id(mut self, op: OpId) -> Self {
        self.op_id = Some(op);
        self
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
        let ((u0, u1), (v0, v1)) = wall_surface.parameter_domain();

        // Shared ring edges: the wall's exact u-boundary isocurves (rational
        // circles at the first / last profile heights; the revolve puts the
        // profile in u and the azimuth in v).
        let bottom_edge = closed_ring_edge(store, wall_surface.isocurve_u(u0)?)?;
        let top_edge = closed_ring_edge(store, wall_surface.isocurve_u(u1)?)?;

        let wall_wire = store.add_wire(WireData {
            edges: vec![
                OrientedEdge::new(bottom_edge, true),
                OrientedEdge::new(top_edge, false),
            ],
            is_closed: true,
        });
        let pcurves = vec![
            FacePcurve {
                edge: bottom_edge,
                curve: iso_pcurve_v(v0, v1, u0)?,
            },
            FacePcurve {
                edge: top_edge,
                curve: iso_pcurve_v(v0, v1, u1)?,
            },
        ];
        let wall = MakeNurbsFace::new(wall_surface)
            .with_boundary(wall_wire, pcurves)
            .execute(store)?;

        // Caps: exact disks on the SAME ring edges, normals away from the body.
        let bottom = cap_from_ring(store, bottom_edge, -Vector3::z())?;
        let top = cap_from_ring(store, top_edge, Vector3::z())?;

        if let Some(op) = &self.op_id {
            bind_created_face(store, wall, op, FaceRole::Wall);
            bind_created_face(store, bottom, op, FaceRole::CapStart);
            bind_created_face(store, top, op, FaceRole::CapEnd);
            bind_created_edge(store, bottom_edge, op, EdgeRole::RingStart);
            bind_created_edge(store, top_edge, op, EdgeRole::RingEnd);
        }

        Ok(finish_solid(store, vec![wall, bottom, top]))
    }
}

/// Degree-1 pcurve mapping a ring edge's parameter onto the `v` axis at a
/// fixed `u` (`t → (u, t)` over `[v0, v1]`).
fn iso_pcurve_v(v0: f64, v1: f64, u: f64) -> Result<NurbsCurve2D> {
    NurbsCurve2D::from_unweighted(
        vec![Point2::new(u, v0), Point2::new(u, v1)],
        KnotVector::new(vec![v0, v0, v1, v1])?,
        1,
    )
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
#[allow(clippy::unwrap_used, clippy::expect_used)]
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

    /// F4a determinism: rebuilding the same op with the same op id into a
    /// fresh store resolves every persistent name, and the resolved faces
    /// carry identical geometry.
    #[test]
    fn creation_names_are_rebuild_stable() {
        use crate::topology::{EdgeName, EdgeRole, FaceName, FaceRole, OpId};

        let build = || {
            let mut store = TopologyStore::new();
            let solid = MakeNurbsTube::new(Point3::new(2.0, 2.0, -1.0), 0.8, 4.0)
                .with_op_id(OpId::new("tube1"))
                .execute(&mut store)
                .unwrap();
            (store, solid)
        };
        let (store_a, _) = build();
        let (store_b, _) = build();

        let side = FaceName::Created {
            op: OpId::new("tube1"),
            role: FaceRole::Side(0),
        };
        let cap_end = FaceName::Created {
            op: OpId::new("tube1"),
            role: FaceRole::CapEnd,
        };
        let ring_start = EdgeName::Created {
            op: OpId::new("tube1"),
            role: EdgeRole::RingStart,
        };

        for name in [&side, &cap_end] {
            let fa = store_a.names().face(name).expect("resolves in build A");
            let fb = store_b.names().face(name).expect("resolves in build B");
            // Same geometry across rebuilds: compare a surface sample.
            let sample = |store: &TopologyStore, f| match &store.face(f).unwrap().surface {
                FaceSurface::Nurbs(s) => s.point_at(0.3, 0.6).unwrap(),
                FaceSurface::Plane(p) => *p.origin(),
                _ => panic!("unexpected surface"),
            };
            let pa = sample(&store_a, fa);
            let pb = sample(&store_b, fb);
            assert!(
                (pa - pb).norm() < 1e-12,
                "rebuilt face for {name:?} moved: {pa:?} vs {pb:?}"
            );
        }
        assert!(store_a.names().edge(&ring_start).is_some());
        assert!(store_b.names().edge(&ring_start).is_some());

        // Without an op id nothing is registered.
        let mut store = TopologyStore::new();
        MakeNurbsTube::new(Point3::new(2.0, 2.0, -1.0), 0.8, 4.0)
            .execute(&mut store)
            .unwrap();
        assert!(store.names().face(&side).is_none());
    }

    /// F4a: the slab and wall builders name all six faces deterministically.
    #[test]
    fn slab_and_wall_faces_are_named() {
        use crate::topology::{FaceName, FaceRole, OpId};

        let mut store = TopologyStore::new();
        MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .with_op_id(OpId::new("slab1"))
            .execute(&mut store)
            .unwrap();
        for role in [
            FaceRole::Top,
            FaceRole::Bottom,
            FaceRole::Side(0),
            FaceRole::Side(1),
            FaceRole::Side(2),
            FaceRole::Side(3),
        ] {
            let name = FaceName::Created {
                op: OpId::new("slab1"),
                role: role.clone(),
            };
            assert!(
                store.names().face(&name).is_some(),
                "slab face {role:?} must be named"
            );
        }

        let mut store = TopologyStore::new();
        MakeCurvedWall::new(Point3::origin(), 5.0, 0.2, 1.4, 3.0, 0.4)
            .with_op_id(OpId::new("wall1"))
            .execute(&mut store)
            .unwrap();
        for role in [
            FaceRole::Side(0),
            FaceRole::Side(1),
            FaceRole::Bottom,
            FaceRole::Top,
            FaceRole::Side(2),
            FaceRole::Side(3),
        ] {
            let name = FaceName::Created {
                op: OpId::new("wall1"),
                role: role.clone(),
            };
            assert!(
                store.names().face(&name).is_some(),
                "wall face {role:?} must be named"
            );
        }
    }

    /// F4a: the revolved solid names its wall, caps, and ring edges.
    #[test]
    fn revolved_solid_faces_are_named() {
        use crate::topology::{EdgeName, EdgeRole, FaceName, FaceRole, OpId};

        let mut store = TopologyStore::new();
        MakeRevolvedSolid::new(vec![(2.0, 0.0), (2.4, 1.2), (2.1, 2.4)])
            .with_op_id(OpId::new("vase1"))
            .execute(&mut store)
            .unwrap();
        for role in [FaceRole::Wall, FaceRole::CapStart, FaceRole::CapEnd] {
            let name = FaceName::Created {
                op: OpId::new("vase1"),
                role: role.clone(),
            };
            assert!(store.names().face(&name).is_some(), "{role:?} named");
        }
        for role in [EdgeRole::RingStart, EdgeRole::RingEnd] {
            let name = EdgeName::Created {
                op: OpId::new("vase1"),
                role,
            };
            assert!(store.names().edge(&name).is_some(), "{role:?} named");
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

    /// A representative curved wall matching the P8 pattern geometry.
    fn sample_wall(store: &mut TopologyStore) -> SolidId {
        use std::f64::consts::PI;
        MakeCurvedWall::new(
            Point3::origin(),
            8.0,
            55.0 * PI / 180.0,
            125.0 * PI / 180.0,
            6.0,
            0.4,
        )
        .execute(store)
        .unwrap()
    }

    #[test]
    fn curved_wall_builds_six_nurbs_faces() {
        let mut store = TopologyStore::new();
        let solid = sample_wall(&mut store);
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        assert_eq!(shell.faces.len(), 6, "inner+outer+top+bottom+2 ends");
        let nurbs = shell
            .faces
            .iter()
            .filter(|&&f| matches!(store.face(f).unwrap().surface, FaceSurface::Nurbs(_)))
            .count();
        assert_eq!(nurbs, 6, "all curved-wall faces are NURBS");
    }

    #[test]
    fn curved_wall_inner_face_sits_at_inner_radius() {
        let mut store = TopologyStore::new();
        let solid = sample_wall(&mut store);
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        // faces[0] is the inner face by construction order.
        let inner = nurbs_of(&store, shell.faces[0]);
        let inner_radius = 8.0 - 0.2; // radius - thickness/2
        let ((u0, u1), (v0, v1)) = inner.parameter_domain();
        for i in 0..=20 {
            for j in 0..=4 {
                let u = u0 + (u1 - u0) * f64::from(i) / 20.0;
                let v = v0 + (v1 - v0) * f64::from(j) / 4.0;
                let p = inner.point_at(u, v).unwrap();
                let dxy = (p.x * p.x + p.y * p.y).sqrt();
                assert!(
                    (dxy - inner_radius).abs() < 1e-9,
                    "inner sample off radius: {dxy} vs {inner_radius}"
                );
            }
        }
    }

    #[test]
    fn curved_wall_tessellates_manifold() {
        let mut store = TopologyStore::new();
        let solid = sample_wall(&mut store);
        assert_manifold(&store, solid);
    }

    #[test]
    fn curved_wall_boundaries_conform() {
        use crate::tessellation::max_adjacent_boundary_deviation;
        let mut store = TopologyStore::new();
        let solid = sample_wall(&mut store);
        let dev = max_adjacent_boundary_deviation(&store, solid);
        assert!(
            dev < 1e-6,
            "curved wall boundary deviation {dev} exceeds 1e-6"
        );
    }

    #[test]
    fn curved_wall_rejects_degenerate_inputs() {
        let mut store = TopologyStore::new();
        // thickness larger than 2*radius → inner radius non-positive.
        assert!(
            MakeCurvedWall::new(Point3::origin(), 0.1, 0.0, 1.0, 2.0, 1.0)
                .execute(&mut store)
                .is_err()
        );
        // zero height.
        assert!(
            MakeCurvedWall::new(Point3::origin(), 8.0, 0.0, 1.0, 0.0, 0.4)
                .execute(&mut store)
                .is_err()
        );
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
