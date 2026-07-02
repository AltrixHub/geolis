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
//! - A target face with NO loop is dropped. In the through-cut topology the tool
//!   passes fully through the target, so the cut faces (trimmed to their SSI
//!   loops) already bound the kept plug; every un-cut target face lies outside
//!   the tool and is discarded. (A point-in-solid classification against the
//!   NURBS-faced tool is deliberately NOT attempted: a robust one is out of
//!   scope for v1 and the topological invariant makes it unnecessary — the kept
//!   plug is exactly the cut discs plus the tool band.)
//!
//! [`subtract_through_cut`]: super::assemble::subtract_through_cut

use std::collections::HashMap;

use crate::error::{OperationError, Result};
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

    // Copy the cut faces (each trimmed to the inside below), recording
    // original -> copy for loop remapping. Un-cut target faces lie outside the
    // through tool and are dropped (see the module-level topology invariant).
    let mut id_map: HashMap<FaceId, FaceId> = HashMap::new();
    let mut result_faces: Vec<FaceId> = Vec::new();
    for &fid in &target_faces {
        if faces_with_loops.contains_key(&fid) {
            let copy = copy_face(store, fid)?;
            id_map.insert(fid, copy);
            result_faces.push(copy);
        }
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

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
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

    /// The chained frame boolean: `Subtract(Intersect(blank, outer_prism),
    /// inner_prism)`. The second (subtract) cut punches a hole into a face whose
    /// `trim.outer` is ALREADY the intersect's SSI window loop, so the frame ring
    /// must come out manifold with two loop-outer discs each carrying one hole.
    #[test]
    fn chained_frame_ring_is_manifold_with_loop_outer_holes() {
        use crate::geometry::nurbs::NurbsCurve3D;
        use crate::math::Vector3;
        use crate::operations::boolean::{Intersect, Subtract};
        use crate::operations::creation::{MakeCurvedWall, MakeNurbsPrism};
        use crate::topology::FaceSurface;
        use std::f64::consts::PI;

        let mut store = TopologyStore::new();
        let deg = |d: f64| d * PI / 180.0;

        // Frame blank: a slightly thicker curved wall spanning just past the
        // window angular extent.
        let blank = MakeCurvedWall::new(Point3::origin(), 8.0, deg(78.0), deg(102.0), 6.0, 0.52)
            .execute(&mut store)
            .unwrap();

        // The window box: a radial rounded-rectangle prism through the wall.
        let window_center = Point3::new(0.0, 7.4, 3.0);
        let tangential = Vector3::x();
        let vertical = Vector3::z();
        let radial = Vector3::new(0.0, 1.2, 0.0);

        let outer_profile =
            NurbsCurve3D::rounded_rectangle(window_center, tangential, vertical, 2.6, 2.0, 0.35)
                .unwrap();
        let outer_prism = MakeNurbsPrism::new(outer_profile, radial)
            .execute(&mut store)
            .unwrap();

        let plate = Intersect::new(blank, outer_prism)
            .execute(&mut store)
            .unwrap();

        // Inner corner radius 0.3: the SSI marcher cannot follow the tighter
        // 0.25 corner cleanly (it fragments the loop into open branches), so the
        // inner opening uses the slightly rounder 0.3 the marcher handles.
        let inner_profile =
            NurbsCurve3D::rounded_rectangle(window_center, tangential, vertical, 2.1, 1.5, 0.3)
                .unwrap();
        let inner_prism = MakeNurbsPrism::new(inner_profile, radial)
            .execute(&mut store)
            .unwrap();

        let frame = Subtract::new(plate, inner_prism)
            .execute(&mut store)
            .unwrap();

        // Manifold.
        let mesh = TessellateSolid::new(frame, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty(), "empty frame mesh");
        let mut counts: StdHashMap<(u32, u32), usize> = StdHashMap::new();
        for tri in &mesh.indices {
            for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }
        for (&(a, b), &c) in &counts {
            assert!(c == 1 || c == 2, "frame edge ({a},{b}) used {c} times");
        }

        // Two disc faces whose trim.outer is a real SSI loop (many segments, not
        // the 4-corner full-domain rectangle) AND carrying exactly one hole — the
        // subtract cut correctly appended a hole to the intersect's loop-outer.
        let shell = store
            .shell(store.solid(frame).unwrap().outer_shell)
            .unwrap();
        let loop_outer_with_hole = shell
            .faces
            .iter()
            .filter(|&&f| {
                let face = store.face(f).unwrap();
                let Some(trim) = &face.trim else { return false };
                matches!(face.surface, FaceSurface::Nurbs(_))
                    && trim.outer.curves.len() > 4
                    && trim.holes.len() == 1
            })
            .count();
        assert_eq!(
            loop_outer_with_hole, 2,
            "frame must have two loop-outer discs each with one hole"
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
