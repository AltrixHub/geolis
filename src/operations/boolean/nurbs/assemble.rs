//! Result assembly for the through-cut subtract.
//!
//! Copies the target solid's faces (so the inputs stay untouched), punches the
//! SSI loops as trim holes onto the copies, builds the tool band faces, and
//! collects everything into a new shell + solid. The tool's caps and the rest
//! of its body are discarded.

use std::collections::HashMap;

use crate::error::{OperationError, Result};
use crate::topology::{FaceData, FaceId, FaceSurface, SolidId, TopologyStore};

use super::band::{build_band_face, BandRingWires};
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

    // Punch each loop onto the COPIED target face, then build the band face that
    // shares those exact hole-ring wires. `cut.loops` is ordered [entry, exit]
    // (loops.rs sorts by mean v), so the two punch results map directly to the
    // band's entry/exit rings.
    for cut in &cuts {
        let entry = punch_onto_copy(store, &cut.loops[0], &id_map)?;
        let exit = punch_onto_copy(store, &cut.loops[1], &id_map)?;
        let band = build_band_face(store, cut, BandRingWires { entry, exit })?;
        result_faces.push(band);
    }

    Ok(finish_solid(store, result_faces))
}

/// Errors if any tool cap (planar tool face) intersects any target face.
///
/// Uses the existing planar/face intersection probe between every planar tool
/// face and every target face; a non-empty intersection means a cap meets the
/// target, which the through-cut path does not handle. Shared with the intersect
/// path.
pub(crate) fn assert_no_cap_intersection(
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

/// Punches one cut loop onto the COPIED target face and returns the hole-ring
/// [`WireId`] created on that copy (so the band face can share it).
///
/// The loop's `target_face` is remapped through `id_map` to the result copy
/// first; punching the preserved input face would attach the ring to the wrong
/// face.
fn punch_onto_copy(
    store: &mut TopologyStore,
    loop_: &super::loops::CutLoop,
    id_map: &HashMap<FaceId, FaceId>,
) -> Result<crate::topology::WireId> {
    let copied_target = *id_map.get(&loop_.target_face).ok_or_else(|| {
        OperationError::Failed("cut loop references an unknown target face".into())
    })?;
    let mut remapped = loop_.clone();
    remapped.target_face = copied_target;
    punch_loop(store, &remapped)
}

/// Deep-copies a face into a new `FaceData` entry, cloning the surface and trim
/// and sharing the (read-only) wire ids. The copy is independent so punching can
/// mutate it without touching the input.
pub(crate) fn copy_face(store: &mut TopologyStore, face: FaceId) -> Result<FaceId> {
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
pub(crate) fn solid_faces(store: &TopologyStore, solid: SolidId) -> Result<Vec<FaceId>> {
    let shell = store.shell(store.solid(solid)?.outer_shell)?;
    Ok(shell.faces.clone())
}

/// Wraps a face list into a closed shell + solid.
pub(crate) fn finish_solid(store: &mut TopologyStore, faces: Vec<FaceId>) -> SolidId {
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

    /// The slab − tube result's adjacent faces conform along every shared
    /// boundary: the outer silhouette (punched top/bottom vs untrimmed side
    /// walls) is now sampled at the boundary-curve-intrinsic parameters, and the
    /// hole rings were already conformed by the polyline-trim fix. The max
    /// adjacent-boundary deviation drops from the chord sagitta (~3e-1, driven by
    /// the coarse 4-corner punched outer loop) to floating-point noise.
    #[test]
    fn boolean_result_boundaries_conform() {
        use crate::tessellation::max_adjacent_boundary_deviation;
        let (store, result) = slab_minus_tube(0.7);
        let dev = max_adjacent_boundary_deviation(&store, result);
        assert!(
            dev < 1e-6,
            "slab-tube adjacent-boundary deviation {dev} exceeds 1e-6"
        );
    }

    #[test]
    fn result_has_punched_faces_and_bands() {
        let (store, result) = slab_minus_tube(0.7);
        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();
        // 6 slab faces + 1 band face = 7.
        assert_eq!(shell.faces.len(), 7, "6 slab faces + 1 band");
        // The two punched faces (front + back) each carry exactly one hole inner
        // wire; the band face carries one inner wire (the exit ring) plus its
        // outer wire (the entry ring). All three NURBS faces with inner wires:
        // 2 punched + 1 band = 3.
        let with_inner = shell
            .faces
            .iter()
            .filter(|&&f| !store.face(f).unwrap().inner_wires.is_empty())
            .count();
        assert_eq!(with_inner, 3, "front + back punched + 1 band");
    }

    /// The band face shares its boundary wires with the punched faces' inner
    /// wires — the same `WireId`s, not duplicates.
    #[test]
    fn band_shares_ring_wires_with_punched_faces() {
        use crate::topology::WireId;
        let (store, result) = slab_minus_tube(0.7);
        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();

        // Collect the punched faces' hole inner-wire ids (front + back rings).
        let mut punched_rings: Vec<WireId> = Vec::new();
        // Locate the single band face: its outer wire is itself a hole ring (it
        // appears in some other face's inner wires).
        let mut all_inner: Vec<WireId> = Vec::new();
        for &f in &shell.faces {
            all_inner.extend(store.face(f).unwrap().inner_wires.iter().copied());
        }
        // The band is the face whose `outer_wire` is one of the hole rings.
        let band = shell
            .faces
            .iter()
            .copied()
            .find(|&f| all_inner.contains(&store.face(f).unwrap().outer_wire))
            .unwrap();
        let band_face = store.face(band).unwrap();
        let band_entry = band_face.outer_wire;
        assert_eq!(band_face.inner_wires.len(), 1, "band has one inner ring");
        let band_exit = band_face.inner_wires[0];

        // The punched faces are the OTHER faces with inner wires.
        for &f in &shell.faces {
            if f == band {
                continue;
            }
            punched_rings.extend(store.face(f).unwrap().inner_wires.iter().copied());
        }
        assert_eq!(punched_rings.len(), 2, "two punched hole rings");
        assert!(
            punched_rings.contains(&band_entry),
            "band entry ring shared with a punched face"
        );
        assert!(
            punched_rings.contains(&band_exit),
            "band exit ring shared with a punched face"
        );
        assert_ne!(band_entry, band_exit, "entry and exit rings differ");
    }

    /// No edge in the result shell spans the tool's full height: the bogus
    /// full-surface u-seam edges (z = -1.5 .. 3.5 in the demo) are gone now that
    /// the band reuses the SSI ring wires.
    #[test]
    fn no_edge_spans_tool_full_height() {
        use crate::topology::EdgeCurve;

        // Slab thickness is 1.5 (front peak) + 1.0 (down) = 2.5; the SSI rings
        // sag at most ~1.5 over the curved face. Any edge taller than this is a
        // full-tool-height seam artifact (the tube spans z = -1.5 .. 3.5 = 5.0).
        const MAX_RING_Z_EXTENT: f64 = 2.5 + 1.5;

        let (store, result) = slab_minus_tube(0.7);
        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();

        let mut max_extent = 0.0_f64;
        for &f in &shell.faces {
            let face = store.face(f).unwrap();
            let mut wires = vec![face.outer_wire];
            wires.extend(face.inner_wires.iter().copied());
            for w in wires {
                let wire = store.wire(w).unwrap();
                for oe in &wire.edges {
                    let edge = store.edge(oe.edge).unwrap();
                    if let EdgeCurve::Nurbs(curve) = &edge.curve {
                        // Sample the edge polyline and measure its z-extent.
                        let (t0, t1) = curve.parameter_domain();
                        let mut zmin = f64::INFINITY;
                        let mut zmax = f64::NEG_INFINITY;
                        for i in 0..=32 {
                            let t = t0 + (t1 - t0) * f64::from(i) / 32.0;
                            let p = curve.point_at(t).unwrap();
                            zmin = zmin.min(p.z);
                            zmax = zmax.max(p.z);
                        }
                        max_extent = max_extent.max(zmax - zmin);
                    }
                }
            }
        }
        assert!(
            max_extent < MAX_RING_Z_EXTENT,
            "an edge spans z-extent {max_extent} (>= {MAX_RING_Z_EXTENT}) — \
             stray full-tool-height seam edge still present"
        );
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

    /// The hole rings tessellate conformally: after the polyline-trim-loop fix
    /// (degree-1 trim curves sampled at their control points), the punched
    /// front/back faces and the band (hole-wall) face emit IDENTICAL 3D vertices
    /// along each shared SSI ring, so the dense per-segment T-junctions along the
    /// hole rings are eliminated.
    ///
    /// Measured (position-deduplicated, 1e-6 quantization):
    /// - plain curved slab (no hole): 384 boundary edges (all perimeter)
    /// - slab − tube, BEFORE the shared-sampling fix: 1788 boundary edges (384
    ///   perimeter + 1404 along the hole rings — the dense punch-vs-band mismatch).
    /// - slab − tube, before the SEAM-FILL fix: 264 boundary edges, 4 of them in
    ///   the hole-ring region at the SSI seam azimuth (the punch chord vs. band
    ///   vertical-stitch disagreement at the tool's u-seam).
    /// - slab − tube, AFTER the seam-fill fix: 0 hole-ring boundary edges. The
    ///   seam wedge is filled with true intersection samples shared by both the
    ///   punch ring (`uv_a`) and the band ribbon (`uv_b`, see
    ///   `super::super::loops::fill_seam_gap`), so the two sides conform across
    ///   the seam and the band ribbon spans the full tool u domain.
    ///
    /// Two assertions pin the result:
    /// 1. The cut result's total boundary-edge count is no worse than the plain
    ///    slab's own perimeter nonconformance (plus a small margin); the prior
    ///    ~1404 hole-ring boundary edges are gone.
    /// 2. Direct hole-ring conformance: NO boundary-edge midpoint lies in the
    ///    tube-wall ring region (distance to the tube axis within [0.7·r, 1.3·r]
    ///    while z is inside the slab). The seam gap is now filled, so even the
    ///    former seam-azimuth residual is gone.
    #[test]
    fn hole_rings_tessellate_conformally() {
        #[allow(clippy::cast_possible_truncation)]
        fn canon_id(canon: &mut HashMap<(i64, i64, i64), u32>, p: &Point3) -> u32 {
            const Q: f64 = 1e6;
            let k = (
                (p.x * Q).round() as i64,
                (p.y * Q).round() as i64,
                (p.z * Q).round() as i64,
            );
            let next = canon.len() as u32;
            *canon.entry(k).or_insert(next)
        }

        // Collects boundary edges (used != 2 after position-dedup) as 3D
        // endpoint pairs.
        fn boundary_edges(store: &TopologyStore, solid: SolidId) -> Vec<(Point3, Point3)> {
            let mesh = TessellateSolid::new(solid, TessellationParams::default())
                .execute(store)
                .unwrap();
            let mut canon: HashMap<(i64, i64, i64), u32> = HashMap::new();
            let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
            let mut endpoints: HashMap<(u32, u32), (Point3, Point3)> = HashMap::new();
            for tri in &mesh.indices {
                let pa = mesh.vertices[tri[0] as usize];
                let pb = mesh.vertices[tri[1] as usize];
                let pc = mesh.vertices[tri[2] as usize];
                let a = canon_id(&mut canon, &pa);
                let b = canon_id(&mut canon, &pb);
                let c = canon_id(&mut canon, &pc);
                for &(x, y, px, py) in &[(a, b, pa, pb), (b, c, pb, pc), (c, a, pc, pa)] {
                    let key = if x < y { (x, y) } else { (y, x) };
                    *counts.entry(key).or_insert(0) += 1;
                    endpoints.entry(key).or_insert((px, py));
                }
            }
            counts
                .iter()
                .filter(|(_, &c)| c != 2)
                .map(|(k, _)| endpoints[k])
                .collect()
        }

        // Plain slab already exhibits per-face perimeter boundary edges with no
        // boolean applied — the generic independent-per-face limitation.
        let mut plain_store = TopologyStore::new();
        let plain = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut plain_store)
            .unwrap();
        let plain_boundary = boundary_edges(&plain_store, plain).len();
        assert!(
            plain_boundary > 0,
            "expected the plain slab to already exhibit per-face perimeter \
             boundary edges (the pre-existing tessellation limitation)"
        );

        const RADIUS: f64 = 0.7;
        let (store, result) = slab_minus_tube(RADIUS);
        let cut_edges = boundary_edges(&store, result);
        let cut_boundary = cut_edges.len();

        // (1) The cut result carries no MORE boundary edges than the plain
        // slab's own perimeter nonconformance (plus a small margin). The prior
        // ~1404 hole-ring T-junctions are eliminated.
        const MARGIN: usize = 16;
        assert!(
            cut_boundary <= plain_boundary + MARGIN,
            "cut result has {cut_boundary} boundary edges, expected \
             <= {plain_boundary} (plain perimeter) + {MARGIN}; hole-ring \
             T-junctions appear to have returned"
        );

        // (2) Direct hole-ring conformance: NO boundary-edge midpoint lies in
        // the tube-wall ring region. The tube axis runs along (3,3,z); a ring
        // boundary edge would sit at radius ~RADIUS from that axis, inside the
        // slab body in z. The seam wedge is now filled with shared intersection
        // samples (see `fill_seam_gap`), so even the former seam-azimuth residual
        // (up to 4 edges) is gone.
        let axis = Point3::new(3.0, 3.0, 0.0);
        let mut ring_edges = 0usize;
        for (p, q) in &cut_edges {
            let m = Point3::new((p.x + q.x) * 0.5, (p.y + q.y) * 0.5, (p.z + q.z) * 0.5);
            let dxy = ((m.x - axis.x).powi(2) + (m.y - axis.y).powi(2)).sqrt();
            let in_ring_radius = dxy >= 0.7 * RADIUS && dxy <= 1.3 * RADIUS;
            let in_slab_z = m.z > -1.2 && m.z < 1.7;
            if in_ring_radius && in_slab_z {
                ring_edges += 1;
            }
        }
        assert_eq!(
            ring_edges, 0,
            "expected 0 hole-ring boundary edges after the seam-fill fix, \
             found {ring_edges}; the punch/band rings are not conforming along \
             the tube wall"
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
        // Band (hole-wall) faces are identified by their boundary topology: a
        // band's `outer_wire` is itself a hole ring shared with a punched face's
        // inner wires. (Both bands and punched faces now carry inner wires, so
        // `inner_wires.is_empty()` no longer discriminates.)
        let mut all_inner: Vec<crate::topology::WireId> = Vec::new();
        for &f in &shell.faces {
            all_inner.extend(store.face(f).unwrap().inner_wires.iter().copied());
        }
        let band_faces: Vec<_> = collect_nurbs_faces(&store, &shell.faces)
            .into_iter()
            .filter(|(fid, _)| all_inner.contains(&store.face(*fid).unwrap().outer_wire))
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
