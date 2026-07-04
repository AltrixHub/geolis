//! Pocket (blind-cut) subtract: the tool enters the target through one face
//! and ends inside it.
//!
//! Result faces for one pocket tool face:
//! - the punched entry face (a trim hole where the tool enters),
//! - a band down the tool side surface from the entry loop to the tool's
//!   buried end,
//! - the buried tool cap, sense-flipped, as the pocket floor.
//!
//! The buried end is resolved locally from the entry geometry (no
//! point-in-solid classification): the tool side direction that goes AGAINST
//! the target face's outward normal at the entry loop leads into the
//! material. Grazing / ambiguous entries are typed errors.
//!
//! Shared-edge topology (F2) makes the floor adjacency structural: the tool's
//! side face and its caps share ring `EdgeId`s, so the band's buried boundary
//! reuses the cap's ring wire, and the band's trim samples the buried ring at
//! the same chord-adaptive parameters the tessellation cache will use — the
//! band and the floor emit identical rim vertices.

use crate::error::{OperationError, Result};
use crate::geometry::surface::Surface;
use crate::math::{Point2, TOLERANCE};
use crate::tessellation::{tessellate_nurbs_curve_params, CurveTessellationOptions};
use crate::topology::{EdgeCurve, EdgeId, FaceId, FaceSurface, TopologyStore, WireId};

use super::loops::CutLoop;
use super::stitch::CutChain;

/// The buried end of a pocket tool.
pub(crate) struct BuriedEnd {
    /// The tool side surface's v-boundary that lies inside the target.
    pub v_boundary: f64,
    /// The shared ring wire at that boundary (the buried cap's outer wire).
    pub ring_wire: WireId,
    /// The buried cap face.
    pub cap_face: FaceId,
}

/// The buried end of a multi-face pocket tool: the shared bottom ring
/// resolved across every side face the entry chain crosses.
pub(crate) struct BuriedChainEnd {
    /// The buried cap face (its flipped copy becomes the pocket floor).
    pub cap_face: FaceId,
    /// Per chain segment, in chain order: that side face's shared ring edge
    /// at the buried end and the face's buried `v` bound.
    pub rings: Vec<(EdgeId, f64)>,
}

/// Which parametric end of a tool side face is buried inside the target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BuriedSide {
    /// The `v0` end.
    Low,
    /// The `v1` end.
    High,
}

/// Resolves which end of the pocket tool is buried inside the target.
///
/// # Errors
///
/// Returns a typed error when the entry is grazing/ambiguous, the tool side
/// face carries no pcurves (not built by the shared-edge creation ops), or no
/// cap shares the buried ring edge.
pub(crate) fn resolve_buried_end(
    store: &TopologyStore,
    entry: &CutLoop,
    tool_faces: &[FaceId],
) -> Result<BuriedEnd> {
    let side = buried_side(store, entry)?;
    let v_boundary = buried_v_bound(store, entry.tool_face, side)?;
    let ring_edge = side_ring_edge(store, entry.tool_face, v_boundary)?;
    let cap_face = cap_sharing_ring_edge(store, tool_faces, entry.tool_face, ring_edge)?;
    Ok(BuriedEnd {
        v_boundary,
        ring_wire: store.face(cap_face)?.outer_wire,
        cap_face,
    })
}

/// Resolves the buried end of a multi-face pocket tool: the buried side is
/// probed on the chain's first segment, then every crossed side face
/// contributes its own shared ring edge at that end, and the buried cap must
/// carry ALL of those edges on its outer wire (the shared bottom ring).
///
/// # Errors
///
/// Returns a typed error when the entry is grazing/ambiguous, a side face
/// lacks pcurves or a boundary edge at the buried end, or no single cap
/// shares the whole buried ring.
pub(crate) fn resolve_buried_chain_end(
    store: &TopologyStore,
    entry: &CutChain,
    tool_faces: &[FaceId],
) -> Result<BuriedChainEnd> {
    let first = entry
        .segments
        .first()
        .ok_or_else(|| OperationError::Failed("empty chained entry loop".into()))?;
    let side = buried_side(store, first)?;

    let mut rings = Vec::with_capacity(entry.segments.len());
    for seg in &entry.segments {
        let v_boundary = buried_v_bound(store, seg.tool_face, side)?;
        let ring_edge = side_ring_edge(store, seg.tool_face, v_boundary)?;
        rings.push((ring_edge, v_boundary));
    }

    let cap_face = cap_sharing_ring_edge(store, tool_faces, first.tool_face, rings[0].0)?;
    let cap_wire = store.wire(store.face(cap_face)?.outer_wire)?;
    for &(edge, _) in &rings {
        if !cap_wire.edges.iter().any(|oe| oe.edge == edge) {
            return Err(OperationError::Failed(
                "pocket subtract found no single tool cap sharing the whole \
                 buried ring"
                    .into(),
            )
            .into());
        }
    }
    Ok(BuriedChainEnd { cap_face, rings })
}

/// Probes which side of the entry loop's tool face points INTO the target
/// material: the direction against the target's outward normal.
fn buried_side(store: &TopologyStore, entry: &CutLoop) -> Result<BuriedSide> {
    // Target outward normal at the entry loop.
    let target = store.face(entry.target_face)?;
    let FaceSurface::Nurbs(target_surf) = &target.surface else {
        return Err(
            OperationError::Failed("pocket subtract requires a NURBS entry face".into()).into(),
        );
    };
    let (target_u, target_v) = mean_uv(&entry.branch.uv_a);
    let mut normal = Surface::normal(target_surf, target_u, target_v)?;
    if !target.same_sense {
        normal = -normal;
    }
    let centroid = target_surf.point_at(target_u, target_v)?;

    // Tool side surface probed a small step toward each v end.
    let side = store.face(entry.tool_face)?;
    let FaceSurface::Nurbs(side_surf) = &side.surface else {
        return Err(OperationError::Failed(
            "pocket subtract requires a NURBS tool side face".into(),
        )
        .into());
    };
    let (side_u, side_v) = mean_uv(&entry.branch.uv_b);
    let (_, (v0, v1)) = side_surf.parameter_domain();
    let delta = 0.05 * (v1 - v0);
    let dot_toward = |v: f64| -> Result<f64> {
        let p = side_surf.point_at(side_u, v.clamp(v0, v1))?;
        Ok((p - centroid).dot(&normal))
    };
    let dot_lo = dot_toward(side_v - delta)?;
    let dot_hi = dot_toward(side_v + delta)?;

    match (dot_lo < -TOLERANCE, dot_hi < -TOLERANCE) {
        (true, false) => Ok(BuriedSide::Low),
        (false, true) => Ok(BuriedSide::High),
        _ => Err(OperationError::Failed(
            "pocket subtract could not resolve the buried tool end \
             (grazing or ambiguous entry)"
                .into(),
        )
        .into()),
    }
}

/// The buried `v` bound of one tool side face for the given buried side.
fn buried_v_bound(store: &TopologyStore, side_face: FaceId, side: BuriedSide) -> Result<f64> {
    let face = store.face(side_face)?;
    let FaceSurface::Nurbs(surf) = &face.surface else {
        return Err(OperationError::Failed(
            "pocket subtract requires a NURBS tool side face".into(),
        )
        .into());
    };
    let (_, (v0, v1)) = surf.parameter_domain();
    Ok(match side {
        BuriedSide::Low => v0,
        BuriedSide::High => v1,
    })
}

/// The shared ring edge of a tool side face at the buried `v` boundary, found
/// via the face's pcurves.
fn side_ring_edge(store: &TopologyStore, side_face: FaceId, v_boundary: f64) -> Result<EdgeId> {
    let face = store.face(side_face)?;
    let side_wire = store.wire(face.outer_wire)?;
    let mut ring_edge = None;
    for oe in &side_wire.edges {
        let Some(pcurve) = face.pcurve_for(oe.edge) else {
            return Err(OperationError::Failed(
                "pocket subtract requires a shared-edge tool (side face \
                 without pcurves)"
                    .into(),
            )
            .into());
        };
        let (t0, t1) = pcurve.parameter_domain();
        let mid = pcurve.point_at(0.5 * (t0 + t1))?;
        if (mid.y - v_boundary).abs() < 1e-9 {
            ring_edge = Some(oe.edge);
        }
    }
    ring_edge.ok_or_else(|| {
        OperationError::Failed(
            "pocket subtract found no side-face boundary edge at the buried \
             end"
            .into(),
        )
        .into()
    })
}

/// The tool cap face whose outer wire contains the given buried ring edge.
fn cap_sharing_ring_edge(
    store: &TopologyStore,
    tool_faces: &[FaceId],
    entry_side_face: FaceId,
    ring_edge: EdgeId,
) -> Result<FaceId> {
    for &f in tool_faces {
        if f == entry_side_face {
            continue;
        }
        let face = store.face(f)?;
        let Ok(wire) = store.wire(face.outer_wire) else {
            continue;
        };
        if wire.edges.iter().any(|oe| oe.edge == ring_edge) {
            return Ok(f);
        }
    }
    Err(OperationError::Failed(
        "pocket subtract found no tool cap sharing the buried ring edge".into(),
    )
    .into())
}

/// UV samples of the buried ring in the tool side surface's parameter space:
/// the straight `v = v_boundary` line at the ring curve's chord-adaptive
/// parameters — the SAME parameters the tessellation cache computes for the
/// shared ring edge, so the band rim and the pocket floor rim coincide.
///
/// # Errors
///
/// Propagates store lookups and curve sampling errors.
pub(crate) fn buried_ring_uv(
    store: &TopologyStore,
    ring_wire: WireId,
    v_boundary: f64,
) -> Result<Vec<Point2>> {
    let wire = store.wire(ring_wire)?;
    let [ring] = wire.edges.as_slice() else {
        return Err(OperationError::Failed(
            "pocket buried ring wire must consist of the single shared ring \
             edge"
                .into(),
        )
        .into());
    };
    buried_edge_uv(store, ring.edge, v_boundary)
}

/// UV samples of ONE buried ring edge in its side face's parameter space: the
/// straight `v = v_boundary` line at the edge curve's chord-adaptive
/// parameters (same-parameter convention: the shared ring edge's curve
/// parameter equals the side surface's `u`), so the band fragment rim and the
/// pocket floor rim coincide.
///
/// # Errors
///
/// Propagates store lookups and curve sampling errors.
pub(crate) fn buried_edge_uv(
    store: &TopologyStore,
    ring_edge: EdgeId,
    v_boundary: f64,
) -> Result<Vec<Point2>> {
    let edge = store.edge(ring_edge)?;
    let EdgeCurve::Nurbs(curve) = &edge.curve else {
        return Err(
            OperationError::Failed("pocket buried ring edge must be a NURBS ring".into()).into(),
        );
    };
    let params = tessellate_nurbs_curve_params(
        curve,
        &CurveTessellationOptions {
            chord_tolerance: 1e-3,
            max_depth: 16,
        },
    )?;
    Ok(params
        .into_iter()
        .map(|t| Point2::new(t, v_boundary))
        .collect())
}

/// Mean of a UV trace.
fn mean_uv(uv: &[Point2]) -> (f64, f64) {
    if uv.is_empty() {
        return (0.0, 0.0);
    }
    #[allow(clippy::cast_precision_loss)]
    let inv = 1.0 / uv.len() as f64;
    (
        uv.iter().map(|p| p.x).sum::<f64>() * inv,
        uv.iter().map(|p| p.y).sum::<f64>() * inv,
    )
}

/// A pocket floor face: a sense-flipped copy of the buried tool cap. Its
/// outward normal pointed away from the tool body; the result solid's floor
/// must face INTO the cavity instead.
pub(crate) fn pocket_floor(store: &mut TopologyStore, cap: FaceId) -> Result<FaceId> {
    let copy = super::assemble::copy_face(store, cap)?;
    let face = store.face_mut(copy)?;
    face.same_sense = !face.same_sense;
    Ok(copy)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::boolean::nurbs::loops::{
        collect_nurbs_faces, extract_cut_loops, ToolFaceCut,
    };
    use crate::operations::creation::{MakeCurvedSlab, MakeNurbsTube};
    use crate::topology::SolidId;

    fn solid_faces(store: &TopologyStore, solid: SolidId) -> Vec<FaceId> {
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        shell.faces.clone()
    }

    /// Slab + tube ending inside it; returns the pocket entry loop and the
    /// tool's face list.
    fn pocket_fixture() -> (TopologyStore, CutLoop, Vec<FaceId>) {
        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        // Back face z ≈ 0 at the tube footprint; tube from below ends at
        // z = 0.5, inside the slab (front ≈ 0.98 there).
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -3.0), 0.7, 3.5)
            .execute(&mut store)
            .unwrap();
        let target = collect_nurbs_faces(&store, &solid_faces(&store, slab));
        let tool = collect_nurbs_faces(&store, &solid_faces(&store, tube));
        let cuts = extract_cut_loops(&target, &tool).unwrap();
        assert_eq!(cuts.len(), 1);
        let ToolFaceCut::Pocket { entry, .. } = &cuts[0] else {
            panic!("expected a pocket cut, got {:?}", cuts[0]);
        };
        let tool_faces = solid_faces(&store, tube);
        (store, entry.clone(), tool_faces)
    }

    #[test]
    fn buried_end_is_the_upper_cap() {
        let (store, entry, tool_faces) = pocket_fixture();
        let buried = resolve_buried_end(&store, &entry, &tool_faces).unwrap();
        // The tube rises from below; its buried end is the TOP (v = v1).
        let side = store.face(entry.tool_face).unwrap();
        let crate::topology::FaceSurface::Nurbs(surf) = &side.surface else {
            panic!("tool side must be NURBS");
        };
        let (_, (_, v1)) = surf.parameter_domain();
        assert!(
            (buried.v_boundary - v1).abs() < 1e-12,
            "expected the v1 (top) end, got {}",
            buried.v_boundary
        );
        // The cap plane sits at the tube's top height z = 0.5.
        let cap = store.face(buried.cap_face).unwrap();
        let crate::topology::FaceSurface::Plane(plane) = &cap.surface else {
            panic!("buried cap must be planar");
        };
        assert!((plane.origin().z - 0.5).abs() < 1e-9);
        // The ring wire's edge is one of the side face's wire edges (shared).
        let ring_edge = store.wire(buried.ring_wire).unwrap().edges[0].edge;
        let side_wire = store.wire(side.outer_wire).unwrap();
        assert!(
            side_wire.edges.iter().any(|oe| oe.edge == ring_edge),
            "buried ring edge must be shared with the side face"
        );
    }

    #[test]
    fn buried_ring_uv_spans_full_u_at_v_boundary() {
        let (store, entry, tool_faces) = pocket_fixture();
        let buried = resolve_buried_end(&store, &entry, &tool_faces).unwrap();
        let uv = buried_ring_uv(&store, buried.ring_wire, buried.v_boundary).unwrap();
        assert!(uv.len() >= 4);
        let side = store.face(entry.tool_face).unwrap();
        let crate::topology::FaceSurface::Nurbs(surf) = &side.surface else {
            panic!("tool side must be NURBS");
        };
        let ((u0, u1), _) = surf.parameter_domain();
        assert!((uv.first().unwrap().x - u0).abs() < 1e-12);
        assert!((uv.last().unwrap().x - u1).abs() < 1e-12);
        assert!(uv.iter().all(|p| (p.y - buried.v_boundary).abs() < 1e-12));
    }
}
