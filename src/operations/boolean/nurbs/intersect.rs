//! Keep-inside through-cut Intersect for NURBS-faced solids (P8).
//!
//! The counterpart of the through-cut [`subtract_through_cut`]. Where subtract
//! keeps the part of the target OUTSIDE the tool (punching SSI loops as trim
//! holes), intersect keeps the part INSIDE the tool: each SSI loop becomes the
//! **outer** trim boundary of its target face (keep the disc/patch inside the
//! loop), and the tool band strip between the entry and exit loops becomes the
//! plug's side wall (`same_sense = true`, normals pointing outward from the
//! plug rather than into a hole).
//!
//! The loops / seam-fill / band machinery is shared verbatim with subtract —
//! only the trim polarity and the band orientation differ.
//!
//! ## v1 topology assumptions
//!
//! - Each target face receives at most one SSI loop (the tool passes cleanly
//!   through, so its single side surface cuts each spanned target face once).
//! - Target faces with no loop are kept whole iff a representative interior
//!   point classifies as [`PointClassification::Inside`] the tool; faces
//!   entirely outside the tool are dropped. For the demo topologies (a blank
//!   fully spanned by the tool in one direction) every no-loop face lies outside
//!   and is dropped.
//!
//! [`subtract_through_cut`]: super::assemble::subtract_through_cut

use std::collections::HashMap;

use nalgebra::Vector3;

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{
    intersect_curve_surface, IntersectionOptions, NurbsCurve3D, NurbsSurface,
};
use crate::math::Point3;
use crate::topology::{FaceId, FaceSurface, FaceTrim, SolidId, TopologyStore, WireId};

use super::assemble::{assert_no_cap_intersection, copy_face, finish_solid, solid_faces};
use super::band::{build_band_face_oriented, BandRingWires};
use super::loops::{collect_nurbs_faces, extract_cut_loops, CutLoop};
use super::punch::{build_ring_wire, ssi_trim_loop};

/// Executes the keep-inside through-cut `target ∩ tool`.
///
/// # Errors
///
/// Returns [`OperationError::Failed`] naming the unsupported case when the
/// target has no NURBS faces, a tool cap intersects the target, the loop
/// preconditions are violated, a target face receives more than one loop, or any
/// geometric sub-step fails.
pub(crate) fn intersect_through_cut(
    store: &mut TopologyStore,
    target: SolidId,
    tool: SolidId,
) -> Result<SolidId> {
    let target_faces = solid_faces(store, target)?;
    let tool_faces = solid_faces(store, tool)?;

    let target_nurbs = collect_nurbs_faces(store, &target_faces);
    let tool_nurbs = collect_nurbs_faces(store, &tool_faces);
    let tool_surfaces: Vec<NurbsSurface> = tool_nurbs.iter().map(|(_, s)| s.clone()).collect();

    if target_nurbs.is_empty() {
        return Err(OperationError::Failed(
            "keep-inside intersect requires a NURBS-faced target".into(),
        )
        .into());
    }

    assert_no_cap_intersection(store, &target_faces, &tool_faces)?;

    let cuts = extract_cut_loops(&target_nurbs, &tool_nurbs)?;

    // Which target faces receive a loop (kept and trimmed to the inside). v1
    // intersect requires at most one loop per target face.
    let mut faces_with_loops: HashMap<FaceId, ()> = HashMap::new();
    for cut in &cuts {
        for l in &cut.loops {
            if faces_with_loops.insert(l.target_face, ()).is_some() {
                return Err(OperationError::Failed(
                    "keep-inside intersect requires at most one loop per target \
                     face (v1 topology)"
                        .into(),
                )
                .into());
            }
        }
    }

    // Copy the faces we keep, recording original -> copy for loop remapping.
    // Faces with a loop are always kept (trimmed to the inside); a no-loop face
    // is kept whole only if it lies inside the tool.
    let mut id_map: HashMap<FaceId, FaceId> = HashMap::new();
    let mut result_faces: Vec<FaceId> = Vec::new();
    for &fid in &target_faces {
        if faces_with_loops.contains_key(&fid) {
            let copy = copy_face(store, fid)?;
            id_map.insert(fid, copy);
            result_faces.push(copy);
        } else if face_is_inside_tool(store, fid, &tool_surfaces)? {
            result_faces.push(copy_face(store, fid)?);
        }
        // else: entirely outside the tool — dropped.
    }

    // Trim each loop's target-face copy to the inside, then build the band
    // (plug wall) sharing the exact ring wires. `same_sense = true` orients the
    // band outward from the kept plug (opposite the subtract hole wall).
    for cut in &cuts {
        let entry = punch_inside_onto_copy(store, &cut.loops[0], &id_map)?;
        let exit = punch_inside_onto_copy(store, &cut.loops[1], &id_map)?;
        let band = build_band_face_oriented(store, cut, BandRingWires { entry, exit }, true)?;
        result_faces.push(band);
    }

    Ok(finish_solid(store, result_faces))
}

/// Trims the COPIED target face of `loop_` to the disc/patch inside the loop and
/// returns the kept-boundary ring [`WireId`] so the band can share it.
fn punch_inside_onto_copy(
    store: &mut TopologyStore,
    loop_: &CutLoop,
    id_map: &HashMap<FaceId, FaceId>,
) -> Result<WireId> {
    let copied_target = *id_map.get(&loop_.target_face).ok_or_else(|| {
        OperationError::Failed("cut loop references an unknown target face".into())
    })?;
    punch_inside(store, copied_target, loop_)
}

/// Sets a target face's trim to keep only the region inside the SSI loop.
///
/// The loop's `uv_a` trace becomes the CCW **outer** trim boundary (holes
/// cleared), and the 3D ring wire replaces the face's boundary wire so the plug
/// disc is bounded by exactly the SSI ring shared with the band.
fn punch_inside(store: &mut TopologyStore, face_id: FaceId, loop_: &CutLoop) -> Result<WireId> {
    if !matches!(store.face(face_id)?.surface, FaceSurface::Nurbs(_)) {
        return Err(OperationError::Failed(
            "keep-inside intersect can only trim NURBS target faces".into(),
        )
        .into());
    }
    let outer = ssi_trim_loop(&loop_.branch.uv_a, false)?;
    let ring = build_ring_wire(store, &loop_.branch.points)?;

    let face = store.face_mut(face_id)?;
    face.trim = Some(FaceTrim::new(outer, Vec::new()));
    face.outer_wire = ring;
    face.inner_wires = Vec::new();
    Ok(ring)
}

/// Whether a representative interior point of `face` lies inside the tool.
///
/// v1 classification (per the plan): cast a long in-plane ray from the face
/// centroid and count crossings with the tool's NURBS side surface(s). An odd
/// crossing count means the probe is enclosed by the tool wall (inside); an even
/// count (0 or 2) means it is outside. This is exact for the through-cut tool
/// topology (a closed extruded profile whose axis is transverse to the ray).
///
/// The ray direction is oblique in the XY plane (and flat in Z, so parity across
/// the side wall is well defined) to avoid grazing the tool's parametric u-seam,
/// which lies along the profile's reference direction. A caps-only or
/// ray-parallel tool is out of scope for v1 and is rejected upstream by the loop
/// preconditions.
fn face_is_inside_tool(
    store: &TopologyStore,
    face: FaceId,
    tool_surfaces: &[NurbsSurface],
) -> Result<bool> {
    let probe = face_centroid(store, face)?;
    // Oblique XY direction (non-axis-aligned slope) so the ray never runs along
    // an axis-aligned tool seam; flat in Z so it stays transverse to the wall.
    let dir = Vector3::new(1.0, 0.436, 0.0);
    let ray = NurbsCurve3D::polyline(&[probe, probe + dir * 1.0e4])?;
    let options = IntersectionOptions::default();
    let mut crossings = 0usize;
    for surface in tool_surfaces {
        crossings += intersect_curve_surface(&ray, surface, &options)?.len();
    }
    Ok(crossings % 2 == 1)
}

/// A representative interior point of a face: the average of its outer boundary
/// wire's edge vertices. For the convex planar side faces of the demo blanks
/// this is the face centroid.
fn face_centroid(store: &TopologyStore, face: FaceId) -> Result<Point3> {
    let outer = store.face(face)?.outer_wire;
    let wire = store.wire(outer)?;
    let mut sum: Vector3<f64> = Vector3::zeros();
    let mut count = 0usize;
    for oe in &wire.edges {
        let edge = store.edge(oe.edge)?;
        for v in [edge.start, edge.end] {
            sum += store.vertex(v)?.point.coords;
            count += 1;
        }
    }
    if count == 0 {
        return Err(OperationError::Failed("face has an empty boundary wire".into()).into());
    }
    #[allow(clippy::cast_precision_loss)]
    Ok(Point3::from(sum / count as f64))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::creation::{MakeCurvedSlab, MakeNurbsTube};
    use crate::tessellation::{TessellateSolid, TessellationParams};
    use std::collections::HashMap as StdHashMap;

    /// Builds slab ∩ tube and returns (store, plug solid).
    fn slab_cap_tube(radius: f64) -> (TopologyStore, SolidId) {
        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -1.5), radius, 5.0)
            .execute(&mut store)
            .unwrap();
        let plug = intersect_through_cut(&mut store, slab, tube).unwrap();
        (store, plug)
    }

    #[test]
    fn plug_has_two_trimmed_discs_and_one_band() {
        let (store, plug) = slab_cap_tube(0.7);
        let shell = store.shell(store.solid(plug).unwrap().outer_shell).unwrap();
        // Two trimmed target discs (front + back inside the loops) + 1 band.
        assert_eq!(shell.faces.len(), 3, "plug = 2 discs + 1 band");
        // Every kept face carries a trim whose outer is the SSI loop and no holes.
        for &f in &shell.faces {
            let face = store.face(f).unwrap();
            let trim = face.trim.as_ref().expect("plug face must be trimmed");
            assert!(trim.holes.is_empty(), "keep-inside face has no holes");
        }
    }

    #[test]
    fn plug_tessellates_manifold() {
        let (store, plug) = slab_cap_tube(0.7);
        let mesh = TessellateSolid::new(plug, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty(), "empty plug mesh");
        let mut counts: StdHashMap<(u32, u32), usize> = StdHashMap::new();
        for tri in &mesh.indices {
            for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }
        for (&(a, b), &c) in &counts {
            assert!(c == 1 || c == 2, "plug edge ({a},{b}) used {c} times");
        }
    }

    #[test]
    fn plug_lies_within_tube_radius_and_spans_thickness() {
        const RADIUS: f64 = 0.7;
        let (store, plug) = slab_cap_tube(RADIUS);
        let mesh = TessellateSolid::new(plug, TessellationParams::default())
            .execute(&store)
            .unwrap();
        let axis = Point3::new(3.0, 3.0, 0.0);
        let mut zmin = f64::INFINITY;
        let mut zmax = f64::NEG_INFINITY;
        for v in &mesh.vertices {
            let dxy = ((v.x - axis.x).powi(2) + (v.y - axis.y).powi(2)).sqrt();
            assert!(
                dxy <= RADIUS + 1e-6,
                "plug vertex ({:.3},{:.3},{:.3}) escapes the tube radius (dxy={dxy})",
                v.x,
                v.y,
                v.z
            );
            zmin = zmin.min(v.z);
            zmax = zmax.max(v.z);
        }
        // The slab body spans roughly z in [-1, 1.5]; the plug fills the tube's
        // slice through it, so its z-extent covers a meaningful thickness.
        assert!(
            zmax - zmin > 1.0,
            "plug z-extent {} too small (zmin={zmin}, zmax={zmax})",
            zmax - zmin
        );
    }

    #[test]
    fn inputs_are_preserved() {
        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -1.5), 0.7, 5.0)
            .execute(&mut store)
            .unwrap();
        let slab_faces = store
            .shell(store.solid(slab).unwrap().outer_shell)
            .unwrap()
            .faces
            .clone();
        let _ = intersect_through_cut(&mut store, slab, tube).unwrap();
        for f in slab_faces {
            let face = store.face(f).unwrap();
            assert!(face.trim.is_none(), "input slab face must stay untrimmed");
        }
    }
}
