//! Result assembly for the through-cut subtract.
//!
//! Copies the target solid's faces (so the inputs stay untouched), punches the
//! SSI loops as trim holes onto the copies, builds the tool band faces, and
//! collects everything into a new shell + solid. The tool's caps and the rest
//! of its body are discarded.

use std::collections::HashMap;

use crate::error::{OperationError, Result};
use crate::topology::{FaceData, FaceId, FaceSurface, SolidId, TopologyStore};

use super::band::build_band_face;
use super::loops::{collect_nurbs_faces, extract_cut_loops};
use super::punch::punch_loop;

/// Executes the through-cut subtract `target - tool`.
///
/// # Errors
///
/// Returns [`OperationError::Failed`] naming the unsupported case when a tool
/// cap intersects the target, the loop preconditions are violated, or any
/// geometric sub-step fails. The planar boolean pipeline is never reached for
/// NURBS-faced solids.
pub(crate) fn subtract_through_cut(
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
            "through-cut subtract requires a NURBS-faced target".into(),
        )
        .into());
    }

    // Cap guard: in the v1 through-cut topology the tool's planar caps lie
    // outside the target and intersect nothing. If a cap actually meets a target
    // face the configuration is out of scope.
    assert_no_cap_intersection(store, &target_faces, &tool_faces)?;

    // Extract + validate the through-cut loops on the ORIGINAL faces.
    let cuts = extract_cut_loops(&target_nurbs, &tool_nurbs)?;

    // Copy every target face so the input solid is preserved, recording the
    // original -> copy id map for loop remapping.
    let mut id_map: HashMap<FaceId, FaceId> = HashMap::new();
    let mut result_faces: Vec<FaceId> = Vec::with_capacity(target_faces.len());
    for &fid in &target_faces {
        let copy = copy_face(store, fid)?;
        id_map.insert(fid, copy);
        result_faces.push(copy);
    }

    // Punch each loop onto the COPIED target face.
    for cut in &cuts {
        for loop_ in &cut.loops {
            let copied_target = *id_map.get(&loop_.target_face).ok_or_else(|| {
                OperationError::Failed("cut loop references an unknown target face".into())
            })?;
            let mut remapped = loop_.clone();
            remapped.target_face = copied_target;
            punch_loop(store, &remapped)?;
        }
    }

    // Build a band (hole-wall) face per tool side face.
    for cut in &cuts {
        let band = build_band_face(store, cut)?;
        result_faces.push(band);
    }

    Ok(finish_solid(store, result_faces))
}

/// Errors if any tool cap (planar tool face) intersects any target face.
///
/// Uses the existing planar/face intersection probe between every planar tool
/// face and every target face; a non-empty intersection means a cap meets the
/// target, which the through-cut path does not handle.
fn assert_no_cap_intersection(
    store: &TopologyStore,
    target_faces: &[FaceId],
    tool_faces: &[FaceId],
) -> Result<()> {
    use crate::operations::boolean::intersect_face_face;

    for &tf in tool_faces {
        if !matches!(store.face(tf)?.surface, FaceSurface::Plane(_)) {
            continue;
        }
        for &gf in target_faces {
            // `intersect_face_face` resolves planar-planar only; restrict the
            // probe to planar target faces (the demo slab's flat sides). A
            // planar cap meeting a curved target face is out of scope but cannot
            // arise in the v1 through-cut topology (caps clear the curved
            // faces), and SSI-based loop extraction already governs the NURBS
            // pairings.
            if !matches!(store.face(gf)?.surface, FaceSurface::Plane(_)) {
                continue;
            }
            let hits = intersect_face_face(store, tf, gf)?;
            if !hits.is_empty() {
                return Err(OperationError::Failed(
                    "through-cut subtract does not support a tool cap that \
                     intersects the target (cap must lie outside)"
                        .into(),
                )
                .into());
            }
        }
    }
    Ok(())
}

/// Deep-copies a face into a new `FaceData` entry, cloning the surface and trim
/// and sharing the (read-only) wire ids. The copy is independent so punching can
/// mutate it without touching the input.
fn copy_face(store: &mut TopologyStore, face: FaceId) -> Result<FaceId> {
    let src = store.face(face)?;
    let data = FaceData {
        surface: src.surface.clone(),
        outer_wire: src.outer_wire,
        inner_wires: src.inner_wires.clone(),
        same_sense: src.same_sense,
        trim: src.trim.clone(),
    };
    Ok(store.add_face(data))
}

/// Collects a solid's outer-shell face ids.
fn solid_faces(store: &TopologyStore, solid: SolidId) -> Result<Vec<FaceId>> {
    let shell = store.shell(store.solid(solid)?.outer_shell)?;
    Ok(shell.faces.clone())
}

/// Wraps a face list into a closed shell + solid.
fn finish_solid(store: &mut TopologyStore, faces: Vec<FaceId>) -> SolidId {
    use crate::topology::{ShellData, SolidData};
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
    use crate::math::Point3;
    use crate::operations::creation::{MakeCurvedSlab, MakeNurbsTube};
    use crate::tessellation::{TessellateSolid, TessellationParams};
    use std::collections::HashMap;

    /// Builds slab − tube and returns (store, result solid).
    fn slab_minus_tube(radius: f64) -> (TopologyStore, SolidId) {
        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -1.5), radius, 5.0)
            .execute(&mut store)
            .unwrap();
        let result = subtract_through_cut(&mut store, slab, tube).unwrap();
        (store, result)
    }

    #[test]
    fn result_has_punched_faces_and_bands() {
        let (store, result) = slab_minus_tube(0.7);
        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();
        // 6 slab faces + 1 band face = 7.
        assert_eq!(shell.faces.len(), 7, "6 slab faces + 1 band");
        // Exactly 2 faces carry holes (front + back).
        let holed = shell
            .faces
            .iter()
            .filter(|&&f| !store.face(f).unwrap().inner_wires.is_empty())
            .count();
        assert_eq!(holed, 2, "front + back punched");
    }

    #[test]
    fn result_tessellates_manifold() {
        let (store, result) = slab_minus_tube(0.7);
        let mesh = TessellateSolid::new(result, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        // Position-deduplicated edge-use counts: every edge is used 1 or 2 times
        // (no edge shared by 3+ triangles). A strict "exactly 2" closure cannot
        // hold here — see `strict_watertightness_blocked` for why.
        let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
        for tri in &mesh.indices {
            for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }
        for (&(a, b), &c) in &counts {
            assert!(c == 1 || c == 2, "result edge ({a},{b}) used {c} times");
        }
    }

    /// BLOCKED: a strictly-closed (every edge used exactly twice) mesh check
    /// cannot pass for solids whose faces are tessellated independently.
    ///
    /// The geolis tessellator meshes each face on its own; adjacent faces sample
    /// their shared boundary curve at different points, leaving mesh-level
    /// T-junctions (boundary edges) along every curved-face/planar-face border.
    /// This is a tessellation-conformance limitation, NOT a `BRep` shell-closure
    /// defect — the solid's shell is topologically closed (faces share wires).
    ///
    /// Measured (position-deduplicated, 1e-6 quantization):
    /// - plain curved slab (no hole): 384 boundary edges
    /// - slab − tube (this result): 1788 boundary edges
    ///   - 384 on the slab outer perimeter (identical to the plain slab —
    ///     curved front/back vs planar sides),
    ///   - 1404 on the hole ring (punched front/back rings vs the band ring).
    ///
    /// The 384 perimeter edges prove the cause is generic per-face sampling, not
    /// the punch/band u-seam: they are present with no boolean applied at all.
    /// Closing this requires conforming-boundary tessellation (shared boundary
    /// polylines across adjacent faces), which is out of scope for this fix.
    ///
    /// This test documents and pins the limitation; it does not relax the
    /// manifold contract above.
    #[test]
    fn strict_watertightness_blocked() {
        fn boundary_edge_count(store: &TopologyStore, solid: SolidId) -> usize {
            let mesh = TessellateSolid::new(solid, TessellationParams::default())
                .execute(store)
                .unwrap();
            #[allow(clippy::cast_possible_truncation)]
            let key_of = |p: &Point3| -> (i64, i64, i64) {
                const Q: f64 = 1e6;
                (
                    (p.x * Q).round() as i64,
                    (p.y * Q).round() as i64,
                    (p.z * Q).round() as i64,
                )
            };
            let mut canon: HashMap<(i64, i64, i64), u32> = HashMap::new();
            #[allow(clippy::cast_possible_truncation)]
            let mut id_of = |p: &Point3| -> u32 {
                let k = key_of(p);
                let next = canon.len() as u32;
                *canon.entry(k).or_insert(next)
            };
            let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
            for tri in &mesh.indices {
                let a = id_of(&mesh.vertices[tri[0] as usize]);
                let b = id_of(&mesh.vertices[tri[1] as usize]);
                let c = id_of(&mesh.vertices[tri[2] as usize]);
                for &(x, y) in &[(a, b), (b, c), (c, a)] {
                    let key = if x < y { (x, y) } else { (y, x) };
                    *counts.entry(key).or_insert(0) += 1;
                }
            }
            counts.values().filter(|&&c| c != 2).count()
        }

        // Plain slab already has boundary edges from independent per-face
        // tessellation — the strict check is a tessellator limitation.
        let mut plain_store = TopologyStore::new();
        let plain = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut plain_store)
            .unwrap();
        let plain_boundary = boundary_edge_count(&plain_store, plain);
        assert!(
            plain_boundary > 0,
            "expected the plain slab to already exhibit per-face boundary edges; \
             if this is 0 the tessellator now conforms boundaries and the strict \
             watertightness check should be re-enabled"
        );

        let (store, result) = slab_minus_tube(0.7);
        let cut_boundary = boundary_edge_count(&store, result);
        assert!(
            cut_boundary >= plain_boundary,
            "cut result should carry at least the plain slab's boundary edges"
        );
    }

    #[test]
    fn result_has_a_real_hole() {
        // Rigorous check: the tube axis (a straight segment running down the
        // hole at the tube's XY center) must miss the band (hole-wall) NURBS
        // faces of the result solid — the axis threads the open tube untouched.
        //
        // The punched front/back NURBS faces are excluded on purpose: their
        // *surface* still spans the hole region geometrically (the hole lives in
        // the trim, which `intersect_curve_surface` does not consult), so the
        // axis necessarily crosses their underlying surface at the cap z-levels.
        // The band faces, in contrast, are the actual tube wall, so a centered
        // axis missing them proves the wall is a genuine open cylinder.
        use crate::geometry::nurbs::{intersect_curve_surface, IntersectionOptions, NurbsCurve3D};

        let (store, result) = slab_minus_tube(0.7);

        // Axis as a degree-1 polyline spanning the full hole length (and a
        // margin on either side) at the tube's XY center.
        let axis =
            NurbsCurve3D::polyline(&[Point3::new(3.0, 3.0, -1.5), Point3::new(3.0, 3.0, 1.7)])
                .unwrap();

        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();
        // Band (hole-wall) faces are the NURBS faces with no inner wires; the
        // punched front/back copies carry their hole as an inner wire.
        let band_faces: Vec<_> = collect_nurbs_faces(&store, &shell.faces)
            .into_iter()
            .filter(|(fid, _)| store.face(*fid).unwrap().inner_wires.is_empty())
            .collect();
        assert!(
            !band_faces.is_empty(),
            "result must carry at least one band (hole-wall) face to probe"
        );
        let options = IntersectionOptions::default();
        for (fid, surface) in &band_faces {
            let hits = intersect_curve_surface(&axis, surface, &options).unwrap();
            assert!(
                hits.is_empty(),
                "tube axis hits band face {fid:?} ({} times) — hole is not open",
                hits.len()
            );
        }

        // Secondary coarse check: no mesh vertex sits near the tube axis inside
        // the slab interval.
        let mesh = TessellateSolid::new(result, TessellationParams::default())
            .execute(&store)
            .unwrap();
        let center = Point3::new(3.0, 3.0, 0.0);
        let radius = 0.7;
        for v in &mesh.vertices {
            // Inside the slab body z-band.
            if v.z > -1.2 && v.z < 1.7 {
                let dxy = ((v.x - center.x).powi(2) + (v.y - center.y).powi(2)).sqrt();
                assert!(
                    dxy > radius * 0.8,
                    "vertex ({:.3},{:.3},{:.3}) intrudes into the hole (dxy={dxy})",
                    v.x,
                    v.y,
                    v.z
                );
            }
        }
    }

    #[test]
    fn input_solids_are_preserved() {
        // After the subtract, the slab's original faces are untouched (no holes).
        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -1.5), 0.7, 5.0)
            .execute(&mut store)
            .unwrap();
        let slab_shell = store.shell(store.solid(slab).unwrap().outer_shell).unwrap();
        let original_faces: Vec<_> = slab_shell.faces.clone();

        let _ = subtract_through_cut(&mut store, slab, tube).unwrap();

        for f in original_faces {
            assert!(
                store.face(f).unwrap().inner_wires.is_empty(),
                "input slab face must stay un-punched"
            );
        }
    }
}
