use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{NurbsCurve2D, NurbsCurve3D, NurbsSurface};
use crate::math::{Point2, Point3, TOLERANCE};
use crate::topology::{
    EdgeCurve, EdgeData, FaceData, FaceId, FacePcurve, FaceSurface, FaceTrim, OrientedEdge,
    TopologyStore, TrimLoop, VertexData, VertexId, WireData, WireId,
};

/// Number of samples taken along each pcurve when approximating its 3D image.
const PCURVE_SAMPLES: usize = 24;

/// Builds a NURBS face: full-domain by default, optionally trimmed.
///
/// An untrimmed face takes its boundary wire from the surface's four exact
/// boundary isocurves. A trimmed face additionally stores UV-space trim loops
/// and builds 3D hole edges by mapping each pcurve through the surface.
/// A solid builder that owns shared boundary edges can instead hand the face a
/// pre-built wire plus its per-edge pcurves via [`Self::with_boundary`].
pub struct MakeNurbsFace {
    surface: NurbsSurface,
    trim: Option<FaceTrim>,
    boundary: Option<(WireId, Vec<FacePcurve>)>,
}

impl MakeNurbsFace {
    /// Creates a new `MakeNurbsFace` for the full surface domain.
    #[must_use]
    pub fn new(surface: NurbsSurface) -> Self {
        Self {
            surface,
            trim: None,
            boundary: None,
        }
    }

    /// Attaches trim loops to the face.
    #[must_use]
    pub fn with_trim(mut self, trim: FaceTrim) -> Self {
        self.trim = Some(trim);
        self
    }

    /// Uses a pre-built outer wire (with edges shared across adjacent faces)
    /// and this face's UV images of those edges instead of fabricating a
    /// private four-isocurve boundary.
    #[must_use]
    pub fn with_boundary(mut self, wire: WireId, pcurves: Vec<FacePcurve>) -> Self {
        self.boundary = Some((wire, pcurves));
        self
    }

    /// Executes the operation, creating the NURBS face in the store.
    ///
    /// # Errors
    ///
    /// Returns an error if a boundary isocurve cannot be extracted, the
    /// boundary degenerates to fewer than one edge, or (when trimmed) any loop
    /// fails validation (open loop, control point outside the parameter domain,
    /// or wrong winding).
    pub fn execute(&self, store: &mut TopologyStore) -> Result<FaceId> {
        let (outer_wire, pcurves) = match &self.boundary {
            Some((wire, pcurves)) => (*wire, pcurves.clone()),
            None => (self.build_boundary_wire(store)?, Vec::new()),
        };

        if let Some(trim) = &self.trim {
            self.validate_trim(trim)?;
        }

        // Hole edges (3D images of the hole pcurves) become inner wires so the
        // 3D topology stays consistent with the UV trim.
        let mut inner_wires = Vec::new();
        if let Some(trim) = &self.trim {
            for hole in &trim.holes {
                inner_wires.push(self.build_hole_wire(store, hole)?);
            }
        }

        let face_id = store.add_face(FaceData {
            surface: FaceSurface::Nurbs(self.surface.clone()),
            outer_wire,
            inner_wires,
            same_sense: true,
            trim: self.trim.clone(),
            pcurves,
        });
        Ok(face_id)
    }

    /// Builds the four-isocurve boundary wire, sharing corner vertices when
    /// corners coincide (closed seams) and skipping zero-length boundary edges.
    fn build_boundary_wire(&self, store: &mut TopologyStore) -> Result<crate::topology::WireId> {
        let [u_min_edge, u_max_edge, v_min_edge, v_max_edge] = self.surface.boundary_curves()?;
        let ((u_min, u_max), (v_min, v_max)) = self.surface.parameter_domain();

        // Corner points in traversal order: (u_min,v_min) -> (u_max,v_min)
        // -> (u_max,v_max) -> (u_min,v_max).
        let c00 = self.surface.point_at(u_min, v_min)?;
        let c10 = self.surface.point_at(u_max, v_min)?;
        let c11 = self.surface.point_at(u_max, v_max)?;
        let c01 = self.surface.point_at(u_min, v_max)?;

        // Share coincident corner vertices.
        let mut corners: Vec<(Point3, VertexId)> = Vec::new();
        let mut vertex_for = |store: &mut TopologyStore, p: Point3| -> VertexId {
            for (cp, id) in &corners {
                if (cp - p).norm() < TOLERANCE {
                    return *id;
                }
            }
            let id = store.add_vertex(VertexData::new(p));
            corners.push((p, id));
            id
        };

        let v00 = vertex_for(store, c00);
        let v10 = vertex_for(store, c10);
        let v11 = vertex_for(store, c11);
        let v01 = vertex_for(store, c01);

        // The four boundary curves stored in their natural isocurve direction,
        // each paired with the wire-traversal orientation needed to walk the
        // rectangle v00 -> v10 -> v11 -> v01 -> v00.
        //   v_min edge: curve in u, natural (u_min,v_min)->(u_max,v_min) = v00->v10, forward.
        //   u_max edge: curve in v, natural (u_max,v_min)->(u_max,v_max) = v10->v11, forward.
        //   v_max edge: curve in u, natural (u_min,v_max)->(u_max,v_max) = v01->v11, reversed.
        //   u_min edge: curve in v, natural (u_min,v_min)->(u_min,v_max) = v00->v01, reversed.
        let segments = [
            (v_min_edge, v00, v10, true),
            (u_max_edge, v10, v11, true),
            (v_max_edge, v01, v11, false),
            (u_min_edge, v00, v01, false),
        ];

        let mut oriented_edges = Vec::with_capacity(4);
        for (curve, natural_start, natural_end, forward) in segments {
            // Skip zero-length boundary segments (collapsed seams).
            if natural_start == natural_end {
                let (t_min, t_max) = curve.parameter_domain();
                let len = (curve.point_at(t_max)? - curve.point_at(t_min)?).norm();
                if len < TOLERANCE {
                    continue;
                }
            }
            let (t_min, t_max) = curve.parameter_domain();
            let edge_id = store.add_edge(EdgeData {
                start: natural_start,
                end: natural_end,
                curve: EdgeCurve::Nurbs(curve),
                t_start: t_min,
                t_end: t_max,
            });
            oriented_edges.push(OrientedEdge::new(edge_id, forward));
        }

        if oriented_edges.is_empty() {
            return Err(OperationError::Failed(
                "NURBS face boundary collapsed to zero edges".into(),
            )
            .into());
        }

        Ok(store.add_wire(WireData {
            edges: oriented_edges,
            is_closed: true,
        }))
    }

    /// Builds a 3D inner wire approximating a hole loop by sampling each pcurve
    /// in UV, mapping samples through the surface, and interpolating a degree-3
    /// NURBS curve (documented approximation).
    fn build_hole_wire(
        &self,
        store: &mut TopologyStore,
        hole: &TrimLoop,
    ) -> Result<crate::topology::WireId> {
        let mut oriented_edges = Vec::with_capacity(hole.curves.len());
        let mut prev_vertex: Option<VertexId> = None;
        let mut first_vertex: Option<VertexId> = None;

        for pcurve in &hole.curves {
            let pts3d = self.sample_pcurve_3d(pcurve)?;
            let (curve, _params) = NurbsCurve3D::interpolate(&pts3d, 3)?;

            let start_point = *pts3d
                .first()
                .ok_or_else(|| OperationError::Failed("empty pcurve sample".to_owned()))?;
            let end_point = *pts3d
                .last()
                .ok_or_else(|| OperationError::Failed("empty pcurve sample".to_owned()))?;

            let start_v = if let Some(v) = prev_vertex {
                v
            } else {
                let v = store.add_vertex(VertexData::new(start_point));
                first_vertex = Some(v);
                v
            };
            // Close back to the first vertex on the final segment when the loop
            // returns to its start.
            let end_v = if let Some(fv) = first_vertex {
                if (end_point - store.vertex(fv)?.point).norm() < TOLERANCE {
                    fv
                } else {
                    store.add_vertex(VertexData::new(end_point))
                }
            } else {
                store.add_vertex(VertexData::new(end_point))
            };

            let (t_min, t_max) = curve.parameter_domain();
            let edge_id = store.add_edge(EdgeData {
                start: start_v,
                end: end_v,
                curve: EdgeCurve::Nurbs(curve),
                t_start: t_min,
                t_end: t_max,
            });
            oriented_edges.push(OrientedEdge::new(edge_id, true));
            prev_vertex = Some(end_v);
        }

        Ok(store.add_wire(WireData {
            edges: oriented_edges,
            is_closed: true,
        }))
    }

    /// Samples a UV pcurve and maps each sample through the surface to 3D.
    fn sample_pcurve_3d(&self, pcurve: &NurbsCurve2D) -> Result<Vec<Point3>> {
        let (t_min, t_max) = pcurve.parameter_domain();
        let mut pts = Vec::with_capacity(PCURVE_SAMPLES + 1);
        for i in 0..=PCURVE_SAMPLES {
            #[allow(clippy::cast_precision_loss)]
            let frac = i as f64 / PCURVE_SAMPLES as f64;
            let t = t_min + frac * (t_max - t_min);
            let uv = pcurve.point_at(t)?;
            pts.push(self.surface.point_at(uv.x, uv.y)?);
        }
        Ok(pts)
    }

    /// Validates trim loops: each loop chains head-to-tail and closes in UV,
    /// all control points lie inside the parameter domain (padded by tolerance),
    /// and windings are correct (outer CCW, holes CW).
    fn validate_trim(&self, trim: &FaceTrim) -> Result<()> {
        let ((u_min, u_max), (v_min, v_max)) = self.surface.parameter_domain();
        let domain = (u_min, u_max, v_min, v_max);

        validate_loop_chain(&trim.outer)?;
        validate_loop_domain(&trim.outer, domain)?;
        if signed_area(&trim.outer)? <= 0.0 {
            return Err(OperationError::InvalidInput(
                "outer trim loop must wind counter-clockwise".into(),
            )
            .into());
        }

        for hole in &trim.holes {
            validate_loop_chain(hole)?;
            validate_loop_domain(hole, domain)?;
            if signed_area(hole)? >= 0.0 {
                return Err(OperationError::InvalidInput(
                    "hole trim loop must wind clockwise".into(),
                )
                .into());
            }
        }
        Ok(())
    }
}

/// Checks that consecutive curves in the loop join head-to-tail and that the
/// loop closes (last end coincides with first start) within tolerance.
fn validate_loop_chain(loop_: &TrimLoop) -> Result<()> {
    if loop_.curves.is_empty() {
        return Err(OperationError::InvalidInput("trim loop has no curves".into()).into());
    }
    let endpoints: Vec<(Point2, Point2)> = loop_
        .curves
        .iter()
        .map(|c| {
            let (t0, t1) = c.parameter_domain();
            Ok((c.point_at(t0)?, c.point_at(t1)?))
        })
        .collect::<Result<_>>()?;

    for i in 0..endpoints.len() {
        let end = endpoints[i].1;
        let next_start = endpoints[(i + 1) % endpoints.len()].0;
        if (end - next_start).norm() > TOLERANCE {
            return Err(OperationError::InvalidInput(format!(
                "trim loop is not closed between curve {i} and {}",
                (i + 1) % endpoints.len()
            ))
            .into());
        }
    }
    Ok(())
}

/// Checks that every control point lies inside the parameter domain padded by
/// `TOLERANCE`.
fn validate_loop_domain(loop_: &TrimLoop, domain: (f64, f64, f64, f64)) -> Result<()> {
    let (u_min, u_max, v_min, v_max) = domain;
    for curve in &loop_.curves {
        for cp in curve.control_points() {
            if cp.x < u_min - TOLERANCE
                || cp.x > u_max + TOLERANCE
                || cp.y < v_min - TOLERANCE
                || cp.y > v_max + TOLERANCE
            {
                return Err(OperationError::InvalidInput(format!(
                    "trim control point ({}, {}) is outside the parameter domain",
                    cp.x, cp.y
                ))
                .into());
            }
        }
    }
    Ok(())
}

/// Signed area of the polygon formed by sampling the loop's pcurves (shoelace).
/// Positive = counter-clockwise.
fn signed_area(loop_: &TrimLoop) -> Result<f64> {
    let mut poly: Vec<Point2> = Vec::new();
    for curve in &loop_.curves {
        let (t0, t1) = curve.parameter_domain();
        // Sample the interior of each segment; drop the tail to avoid the
        // duplicate shared with the next segment's head.
        for i in 0..PCURVE_SAMPLES {
            #[allow(clippy::cast_precision_loss)]
            let frac = i as f64 / PCURVE_SAMPLES as f64;
            let t = t0 + frac * (t1 - t0);
            poly.push(curve.point_at(t)?);
        }
    }
    let n = poly.len();
    let mut area2 = 0.0;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        area2 += a.x * b.y - b.x * a.y;
    }
    Ok(0.5 * area2)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::nurbs::KnotVector;

    /// Bilinear planar patch over [0,4]x[0,4] in the z=0 plane.
    fn planar_patch() -> NurbsSurface {
        NurbsSurface::from_unweighted(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 4.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
                Point3::new(4.0, 4.0, 0.0),
            ],
            2,
            2,
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
            1,
        )
        .unwrap()
    }

    /// A degree-1 polyline pcurve loop (closed rectangle) in UV.
    fn rect_loop(u0: f64, v0: f64, u1: f64, v1: f64, ccw: bool) -> TrimLoop {
        let corners = if ccw {
            [(u0, v0), (u1, v0), (u1, v1), (u0, v1)]
        } else {
            [(u0, v0), (u0, v1), (u1, v1), (u1, v0)]
        };
        let mut curves = Vec::new();
        for i in 0..4 {
            let a = corners[i];
            let b = corners[(i + 1) % 4];
            curves.push(
                NurbsCurve2D::from_unweighted(
                    vec![Point2::new(a.0, a.1), Point2::new(b.0, b.1)],
                    KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
                    1,
                )
                .unwrap(),
            );
        }
        TrimLoop::new(curves)
    }

    #[test]
    fn untrimmed_face_has_four_boundary_edges() {
        let mut store = TopologyStore::new();
        let face_id = MakeNurbsFace::new(planar_patch())
            .execute(&mut store)
            .unwrap();
        let face = store.face(face_id).unwrap();
        let wire = store.wire(face.outer_wire).unwrap();
        assert_eq!(wire.edges.len(), 4, "expected 4 boundary edges");
        assert!(matches!(face.surface, FaceSurface::Nurbs(_)));
        assert!(face.trim.is_none());
    }

    #[test]
    fn untrimmed_face_shares_corner_vertices() {
        let mut store = TopologyStore::new();
        let face_id = MakeNurbsFace::new(planar_patch())
            .execute(&mut store)
            .unwrap();
        let face = store.face(face_id).unwrap();
        let wire = store.wire(face.outer_wire).unwrap();
        // Collect the traversal-start vertex of each oriented edge: a closed
        // 4-edge rectangle uses exactly 4 distinct vertices.
        let mut starts = std::collections::HashSet::new();
        for oe in &wire.edges {
            let edge = store.edge(oe.edge).unwrap();
            let sv = if oe.forward { edge.start } else { edge.end };
            starts.insert(sv);
        }
        assert_eq!(starts.len(), 4, "rectangle boundary must share 4 corners");
    }

    #[test]
    fn trimmed_face_stores_trim_and_hole_wire() {
        let mut store = TopologyStore::new();
        let outer = rect_loop(0.0, 0.0, 1.0, 1.0, true);
        let hole = rect_loop(0.4, 0.4, 0.6, 0.6, false);
        let trim = FaceTrim::new(outer, vec![hole]);
        let face_id = MakeNurbsFace::new(planar_patch())
            .with_trim(trim)
            .execute(&mut store)
            .unwrap();
        let face = store.face(face_id).unwrap();
        assert!(face.trim.is_some());
        assert_eq!(face.inner_wires.len(), 1, "expected one hole wire");
        let hole_wire = store.wire(face.inner_wires[0]).unwrap();
        assert_eq!(hole_wire.edges.len(), 4, "hole loop has 4 pcurves");
    }

    #[test]
    fn open_outer_loop_is_rejected() {
        let mut store = TopologyStore::new();
        let mut outer = rect_loop(0.0, 0.0, 1.0, 1.0, true);
        outer.curves.pop(); // break closure
        let trim = FaceTrim::new(outer, vec![]);
        let result = MakeNurbsFace::new(planar_patch())
            .with_trim(trim)
            .execute(&mut store);
        assert!(result.is_err(), "open loop must be rejected");
    }

    #[test]
    fn out_of_domain_loop_is_rejected() {
        let mut store = TopologyStore::new();
        // Domain is [0,1]x[0,1]; this loop extends to 2.0.
        let outer = rect_loop(0.0, 0.0, 2.0, 1.0, true);
        let trim = FaceTrim::new(outer, vec![]);
        let result = MakeNurbsFace::new(planar_patch())
            .with_trim(trim)
            .execute(&mut store);
        assert!(result.is_err(), "out-of-domain loop must be rejected");
    }

    #[test]
    fn wrong_winding_outer_loop_is_rejected() {
        let mut store = TopologyStore::new();
        // Outer loop wound clockwise (ccw = false) must be rejected.
        let outer = rect_loop(0.0, 0.0, 1.0, 1.0, false);
        let trim = FaceTrim::new(outer, vec![]);
        let result = MakeNurbsFace::new(planar_patch())
            .with_trim(trim)
            .execute(&mut store);
        assert!(result.is_err(), "clockwise outer loop must be rejected");
    }

    #[test]
    fn wrong_winding_hole_loop_is_rejected() {
        let mut store = TopologyStore::new();
        let outer = rect_loop(0.0, 0.0, 1.0, 1.0, true);
        // Hole wound counter-clockwise (ccw = true) must be rejected.
        let hole = rect_loop(0.4, 0.4, 0.6, 0.6, true);
        let trim = FaceTrim::new(outer, vec![hole]);
        let result = MakeNurbsFace::new(planar_patch())
            .with_trim(trim)
            .execute(&mut store);
        assert!(result.is_err(), "ccw hole loop must be rejected");
    }
}
