//! Result assembly for the through-cut subtract.
//!
//! Copies the target solid's faces (so the inputs stay untouched), punches the
//! SSI loops as trim holes onto the copies, builds the tool band faces, and
//! collects everything into a new shell + solid. The tool's caps and the rest
//! of its body are discarded.

use std::collections::{HashMap, HashSet};

use crate::error::{OperationError, Result};
use crate::math::Point2;
use crate::topology::{
    EdgeId, FaceData, FaceId, FaceSurface, OrientedEdge, SolidId, TopologyStore, WireData,
};

use super::band::{build_band_face, BandRingWires};
use super::loops::{collect_nurbs_faces, extract_cut_loops_trimmed};
use super::punch::{punch_loop, ChainRing};
use super::split::{self, ChainTopology, Fragment, TraceRun};

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
    op_id: Option<&crate::topology::OpId>,
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

    // Extract + validate the through-cut loops on the ORIGINAL faces,
    // filtered per face by its kept trim region (split fragments share one
    // parent surface; only the containing fragment receives a branch).
    let target_trims = super::loops::collect_trim_regions(store, &target_faces)?;
    let cuts = extract_cut_loops_trimmed(&target_nurbs, &tool_nurbs, &target_trims)?;

    // ---- F3b: plan the target face splits. -------------------------------
    // Chains that cross target face boundaries do not punch interior holes;
    // their per-face trace runs split the affected faces instead. The trace
    // edges are built ONCE here and shared by the fragments and the band.
    let (topos, affected) = plan_face_splits(store, &cuts, &target_faces)?;

    let (id_map, outcome, mut result_faces) =
        prepare_result_faces(store, &target_faces, &affected, op_id)?;
    let lookup = ResultLookup {
        id_map: &id_map,
        fragments: &outcome.fragments,
    };

    // Cap-plane closure edges of every open (cap-touching) chain, shared
    // between the band fragments and the notched cap rebuild.
    let mut closure_edges: Vec<EdgeId> = Vec::new();

    // Punch each loop onto its RESULT face (copy or containing fragment),
    // then build the band face that shares those exact hole-ring wires.
    // Through loops are ordered [entry, exit] (loops.rs sorts by mean v), so
    // the two punch results map directly to the band's entry/exit rings.
    for (ci, cut) in cuts.iter().enumerate() {
        match cut {
            super::loops::ToolFaceCut::Through { tool_face, loops } => {
                assemble_through(
                    store,
                    *tool_face,
                    loops,
                    [topos.get(&(ci, 0)), topos.get(&(ci, 1))],
                    &lookup,
                    &mut result_faces,
                    op_id,
                )?;
            }
            super::loops::ToolFaceCut::Pocket { tool_face, entry } => {
                // Punch the entry hole, band down to the buried ring, and keep
                // the buried tool cap (sense-flipped) as the pocket floor.
                let buried = super::pocket::resolve_buried_end(store, entry, &tool_faces)?;
                let entry_face = lookup.resolve_loop(entry)?;
                let entry_ring = punch_onto(store, entry, entry_face)?;
                let buried_uv =
                    super::pocket::buried_ring_uv(store, buried.ring_wire, buried.v_boundary)?;
                let band = super::band::build_pocket_band_face(
                    store,
                    *tool_face,
                    entry,
                    &buried_uv,
                    entry_ring,
                    buried.ring_wire,
                )?;
                result_faces.push(band);
                let floor = super::pocket::pocket_floor(store, buried.cap_face)?;
                result_faces.push(floor);
                if let Some(op) = op_id {
                    name_band(store, op, *tool_face, band);
                    name_rim(store, op, entry_face, entry_ring, 0);
                    name_floor(store, op, buried.cap_face, floor);
                }
            }
            super::loops::ToolFaceCut::MultiFaceThrough { chains } => {
                let entry_ring = ring_for_chain(store, &chains[0], topos.get(&(ci, 0)), &lookup)?;
                let exit_ring = ring_for_chain(store, &chains[1], topos.get(&(ci, 1)), &lookup)?;
                assemble_multiface_through(
                    store,
                    chains,
                    &entry_ring,
                    &exit_ring,
                    &lookup,
                    &mut result_faces,
                    &mut closure_edges,
                    op_id,
                )?;
            }
            super::loops::ToolFaceCut::MultiFacePocket { entry } => {
                assemble_multiface_pocket(
                    store,
                    entry,
                    &tool_faces,
                    &lookup,
                    &mut result_faces,
                    op_id,
                )?;
            }
        }
    }

    // Rebuild the notched caps of every cap-touching cut: wire surgery on
    // the kept sub-edges plus the band's cap-plane closure edges (the SAME
    // EdgeIds — watertight by construction).
    rebuild_caps_into_result(store, &outcome, &closure_edges, op_id, &mut result_faces)?;

    name_fragment_rims(store, &affected, &outcome.fragments, op_id);

    Ok(finish_solid(store, result_faces))
}

/// Rebuilds the pending notched caps (if any) and appends their fragments
/// to the result face list. The kept-edge set is collected from the split
/// wall fragments' outer wires — a sub-edge absent from it was removed
/// doorway material and drops out of the rebuilt cap wires.
fn rebuild_caps_into_result(
    store: &mut TopologyStore,
    outcome: &split::SplitOutcome,
    closure_edges: &[EdgeId],
    op_id: Option<&crate::topology::OpId>,
    result_faces: &mut Vec<FaceId>,
) -> Result<()> {
    if outcome.planar_pending.is_empty() && closure_edges.is_empty() {
        return Ok(());
    }
    let mut kept_edges: HashSet<EdgeId> = HashSet::new();
    for frags in outcome.fragments.values() {
        for frag in frags {
            let wire = store.face(frag.face)?.outer_wire;
            for oe in &store.wire(wire)?.edges {
                kept_edges.insert(oe.edge);
            }
        }
    }
    let cap_fragments = super::caps::rebuild_notched_caps(
        store,
        &outcome.planar_pending,
        &outcome.sub_edges,
        &kept_edges,
        closure_edges,
        op_id,
    )?;
    result_faces.extend(cap_fragments);
    Ok(())
}

/// Splits the affected target faces, then copies every UNAFFECTED target
/// face (preserving the input solid and carrying names over unchanged) —
/// EXCEPT the planar caps a cap-touching cut notched, which are rebuilt
/// after the band assembly provides their closure edges. Returns the copy
/// map, the split outcome, and the result face list seeded with the copies
/// + fragments.
#[allow(clippy::type_complexity)]
fn prepare_result_faces(
    store: &mut TopologyStore,
    target_faces: &[FaceId],
    affected: &[(FaceId, Vec<TraceRun>)],
    op_id: Option<&crate::topology::OpId>,
) -> Result<(HashMap<FaceId, FaceId>, split::SplitOutcome, Vec<FaceId>)> {
    let outcome = if affected.is_empty() {
        split::SplitOutcome::default()
    } else {
        split::split_target_faces(store, affected, target_faces, op_id)?
    };
    let mut id_map: HashMap<FaceId, FaceId> = HashMap::new();
    let mut result_faces: Vec<FaceId> = Vec::with_capacity(target_faces.len());
    let affected_ids: Vec<FaceId> = affected.iter().map(|(f, _)| *f).collect();
    for &fid in target_faces {
        if affected_ids.contains(&fid) || outcome.planar_pending.contains(&fid) {
            continue;
        }
        let copy = copy_face(store, fid)?;
        // Persistent names carry over UNCHANGED to the result copies (the
        // newest result owns the name; the input face drops out of the
        // registry). Independent of the boolean's own op id.
        store.names_mut().transfer_face(fid, copy);
        id_map.insert(fid, copy);
        result_faces.push(copy);
    }
    for (fid, _) in affected {
        if let Some(frags) = outcome.fragments.get(fid) {
            result_faces.extend(frags.iter().map(|f| f.face));
        }
    }
    Ok((id_map, outcome, result_faces))
}

/// Plans the F3b face splits: chain trace topologies per (cut, chain) and
/// the per-face trace runs, ordered by the target solid's face order.
#[allow(clippy::type_complexity)]
fn plan_face_splits(
    store: &mut TopologyStore,
    cuts: &[super::loops::ToolFaceCut],
    target_faces: &[FaceId],
) -> Result<(
    HashMap<(usize, usize), ChainTopology>,
    Vec<(FaceId, Vec<TraceRun>)>,
)> {
    let seam_crossing = |store: &TopologyStore, loop_: &super::loops::CutLoop| -> Result<bool> {
        let FaceSurface::Nurbs(surf) = &store.face(loop_.target_face)?.surface else {
            return Ok(false);
        };
        Ok(super::loops::crosses_target_seam(&loop_.branch, surf))
    };

    let mut topos: HashMap<(usize, usize), ChainTopology> = HashMap::new();
    let mut runs_by_face: HashMap<FaceId, Vec<TraceRun>> = HashMap::new();
    for (ci, cut) in cuts.iter().enumerate() {
        match cut {
            super::loops::ToolFaceCut::Through { loops, .. } => {
                // A closed loop straddling the target's parametric seam is
                // split at its exact seam samples; the halves become
                // boundary-notch trace runs on the (single) target face.
                for (li, loop_) in loops.iter().enumerate() {
                    if seam_crossing(store, loop_)? {
                        let FaceSurface::Nurbs(surf) =
                            store.face(loop_.target_face)?.surface.clone()
                        else {
                            return Err(OperationError::Failed(
                                "seam-straddling loop on a non-NURBS target face".into(),
                            )
                            .into());
                        };
                        #[allow(clippy::cast_possible_truncation)]
                        let (topo, runs) = split::split_seam_loop(store, loop_, &surf, li as u32)?;
                        for run in runs {
                            runs_by_face.entry(run.target_face).or_default().push(run);
                        }
                        topos.insert((ci, li), topo);
                    }
                }
            }
            super::loops::ToolFaceCut::MultiFaceThrough { chains } => {
                for (li, chain) in chains.iter().enumerate() {
                    for seg in &chain.segments {
                        if seam_crossing(store, seg)? {
                            return Err(OperationError::Failed(
                                "chained cut loop crossing the target face's \
                                 parametric seam is unsupported"
                                    .into(),
                            )
                            .into());
                        }
                    }
                    // Target-crossing chains split the faces they cross;
                    // OPEN (cap-touching) chains always split — even on a
                    // single target face their trace is a boundary notch,
                    // never an interior hole.
                    if chain.crosses_target_faces() || !chain.closed {
                        let topo = split::build_chain_topology(store, chain)?;
                        #[allow(clippy::cast_possible_truncation)]
                        let runs = split::trace_runs(chain, &topo, li as u32)?;
                        for run in runs {
                            runs_by_face.entry(run.target_face).or_default().push(run);
                        }
                        topos.insert((ci, li), topo);
                    }
                }
            }
            super::loops::ToolFaceCut::Pocket { entry, .. } => {
                if seam_crossing(store, entry)? {
                    return Err(OperationError::Failed(
                        "pocket entry loop crossing the target face's \
                         parametric seam is unsupported"
                            .into(),
                    )
                    .into());
                }
            }
            super::loops::ToolFaceCut::MultiFacePocket { entry } => {
                if entry.crosses_target_faces() {
                    return Err(OperationError::Failed(
                        "pocket cut crossing target face boundaries is unsupported \
                         (a blind cut must enter through a single face)"
                            .into(),
                    )
                    .into());
                }
                for seg in &entry.segments {
                    if seam_crossing(store, seg)? {
                        return Err(OperationError::Failed(
                            "pocket entry loop crossing the target face's \
                             parametric seam is unsupported"
                                .into(),
                        )
                        .into());
                    }
                }
            }
        }
    }
    // Deterministic order: the target solid's face order.
    let affected: Vec<(FaceId, Vec<TraceRun>)> = target_faces
        .iter()
        .filter_map(|f| runs_by_face.remove(f).map(|runs| (*f, runs)))
        .collect();
    Ok((topos, affected))
}

/// Rim names for split fragments: one `CutRim` per fragment carrying a
/// trace, composed from the fragment's (post-split) persistent name.
fn name_fragment_rims(
    store: &mut TopologyStore,
    affected: &[(FaceId, Vec<TraceRun>)],
    fragments: &HashMap<FaceId, Vec<Fragment>>,
    op_id: Option<&crate::topology::OpId>,
) {
    let Some(op) = op_id else {
        return;
    };
    for (fid, _) in affected {
        let Some(frags) = fragments.get(fid) else {
            continue;
        };
        for frag in frags {
            let Some(name) = store.names().name_of_face(frag.face).cloned() else {
                continue;
            };
            store.names_mut().bind_edge(
                frag.first_trace_edge,
                crate::topology::EdgeName::CutRim {
                    op: op.clone(),
                    target: Box::new(name),
                    loop_index: frag.first_loop_index,
                },
            );
        }
    }
}

/// Assembles one single-tool-face through cut: a seam-straddling loop's ring
/// is the wire of its shared trace edges (the fragment already carries the
/// notch in its outer trim); an interior loop is punched as a trim hole. The
/// band face joins the two rings.
fn assemble_through(
    store: &mut TopologyStore,
    tool_face: FaceId,
    loops: &[super::loops::CutLoop; 2],
    topos: [Option<&ChainTopology>; 2],
    lookup: &ResultLookup<'_>,
    result_faces: &mut Vec<FaceId>,
    op_id: Option<&crate::topology::OpId>,
) -> Result<()> {
    let mut rings: [Option<crate::topology::WireId>; 2] = [None, None];
    for (li, loop_) in loops.iter().enumerate() {
        if let Some(topo) = topos[li] {
            rings[li] = Some(
                store.add_wire(WireData {
                    edges: topo
                        .edges
                        .iter()
                        .map(|&e| OrientedEdge::new(e, true))
                        .collect(),
                    is_closed: true,
                }),
            );
        } else {
            let face = lookup.resolve_loop(loop_)?;
            let ring = punch_onto(store, loop_, face)?;
            if let Some(op) = op_id {
                #[allow(clippy::cast_possible_truncation)]
                name_rim(store, op, face, ring, li as u32);
            }
            rings[li] = Some(ring);
        }
    }
    let (Some(entry), Some(exit)) = (rings[0], rings[1]) else {
        return Err(OperationError::Failed("through cut without two rings".into()).into());
    };
    let band = build_band_face(store, tool_face, loops, BandRingWires { entry, exit })?;
    result_faces.push(band);
    if let Some(op) = op_id {
        name_band(store, op, tool_face, band);
    }
    Ok(())
}

/// Resolves original target faces to result faces: unaffected faces map to
/// their copies; split faces map to the fragment containing a UV sample.
struct ResultLookup<'a> {
    id_map: &'a HashMap<FaceId, FaceId>,
    fragments: &'a HashMap<FaceId, Vec<Fragment>>,
}

impl ResultLookup<'_> {
    /// The result face for an interior loop: the copy, or the fragment whose
    /// kept region contains the loop's UV centroid.
    fn resolve_loop(&self, loop_: &super::loops::CutLoop) -> Result<FaceId> {
        let centroid = uv_centroid(&loop_.branch.uv_a);
        self.resolve(loop_.target_face, centroid)
    }

    fn resolve(&self, original: FaceId, uv: Point2) -> Result<FaceId> {
        if let Some(copy) = self.id_map.get(&original) {
            return Ok(*copy);
        }
        if let Some(frags) = self.fragments.get(&original) {
            for frag in frags {
                if split::polygon_contains(&frag.polygon, uv) {
                    return Ok(frag.face);
                }
            }
            return Err(OperationError::Failed(
                "interior cut loop lies on a removed fragment of a split \
                 target face (inconsistent cut)"
                    .into(),
            )
            .into());
        }
        Err(OperationError::Failed("cut loop references an unknown target face".into()).into())
    }
}

/// Mean of a UV trace.
fn uv_centroid(uv: &[Point2]) -> Point2 {
    if uv.is_empty() {
        return Point2::new(0.0, 0.0);
    }
    #[allow(clippy::cast_precision_loss)]
    let inv = 1.0 / uv.len() as f64;
    let (mut su, mut sv) = (0.0, 0.0);
    for p in uv {
        su += p.x;
        sv += p.y;
    }
    Point2::new(su * inv, sv * inv)
}

/// The hole ring of one chained loop: interior chains punch their single
/// result face (trim hole + ring wire, as in Phase B); target-crossing
/// chains already carry their trace edges in the fragment wires, so the
/// ring is assembled from the shared trace topology without punching.
fn ring_for_chain(
    store: &mut TopologyStore,
    chain: &super::stitch::CutChain,
    topo: Option<&ChainTopology>,
    lookup: &ResultLookup<'_>,
) -> Result<ChainRing> {
    if let Some(topo) = topo {
        let wire = store.add_wire(WireData {
            edges: topo
                .edges
                .iter()
                .map(|&e| OrientedEdge::new(e, true))
                .collect(),
            // An open (cap-touching) chain's trace is terminal-to-terminal;
            // the cap-plane closure edges complete the circuit in the band
            // fragments and the notched caps.
            is_closed: chain.closed,
        });
        return Ok(ChainRing {
            wire,
            edges: topo.edges.clone(),
            junctions: topo.junctions.clone(),
            tool_pcurves: topo.tool_pcurves.clone(),
        });
    }
    let mut remapped = chain.clone();
    for seg in &mut remapped.segments {
        let centroid = uv_centroid(&seg.branch.uv_a);
        seg.target_face = lookup.resolve(seg.target_face, centroid)?;
    }
    super::punch::punch_chain(store, &remapped)
}

/// Assembles one multi-face through cut: the two hole rings (punched or
/// split-shared) are joined by one band fragment per tool side face, sharing
/// the ring edges and the new kink edges. Open (cap-touching) chains also
/// produce two cap-plane closure edges, appended to `closure_edges` for the
/// notched cap rebuild.
#[allow(clippy::too_many_arguments)]
fn assemble_multiface_through(
    store: &mut TopologyStore,
    chains: &[super::stitch::CutChain; 2],
    entry_ring: &ChainRing,
    exit_ring: &ChainRing,
    lookup: &ResultLookup<'_>,
    result_faces: &mut Vec<FaceId>,
    closure_edges: &mut Vec<EdgeId>,
    op_id: Option<&crate::topology::OpId>,
) -> Result<()> {
    let fragments = if chains[0].closed {
        super::band::build_band_fragments(store, &chains[0], &chains[1], entry_ring, exit_ring)?
    } else {
        let (fragments, closures) = super::band::build_open_band_fragments(
            store, &chains[0], &chains[1], entry_ring, exit_ring,
        )?;
        closure_edges.extend(closures);
        fragments
    };
    for fragment in &fragments {
        result_faces.push(fragment.face);
    }
    if let Some(op) = op_id {
        for fragment in &fragments {
            name_band(store, op, fragment.tool_face, fragment.face);
        }
        // Rims for interior (single-target) CLOSED chains bind here; rims
        // of split faces (target-crossing or open chains) bind per fragment
        // after the cut loop.
        for (li, chain) in chains.iter().enumerate() {
            if !chain.closed {
                continue;
            }
            if let Some(face) = chain.single_target_face() {
                let ring = if li == 0 { entry_ring } else { exit_ring };
                let resolved = lookup.resolve(face, uv_centroid(&chain.segments[0].branch.uv_a))?;
                #[allow(clippy::cast_possible_truncation)]
                name_rim(store, op, resolved, ring.wire, li as u32);
            }
        }
    }
    Ok(())
}

/// Assembles one multi-face pocket cut: the shared bottom ring is resolved
/// across all crossed side faces, one band fragment per face runs down to its
/// buried ring edge, and the flipped buried cap becomes the floor.
fn assemble_multiface_pocket(
    store: &mut TopologyStore,
    entry: &super::stitch::CutChain,
    tool_faces: &[FaceId],
    lookup: &ResultLookup<'_>,
    result_faces: &mut Vec<FaceId>,
    op_id: Option<&crate::topology::OpId>,
) -> Result<()> {
    let buried = super::pocket::resolve_buried_chain_end(store, entry, tool_faces)?;
    let mut remapped = entry.clone();
    for seg in &mut remapped.segments {
        let centroid = uv_centroid(&seg.branch.uv_a);
        seg.target_face = lookup.resolve(seg.target_face, centroid)?;
    }
    let entry_ring = super::punch::punch_chain(store, &remapped)?;
    let fragments =
        super::band::build_pocket_band_fragments(store, &remapped, &entry_ring, &buried)?;
    for fragment in &fragments {
        result_faces.push(fragment.face);
    }
    let floor = super::pocket::pocket_floor(store, buried.cap_face)?;
    result_faces.push(floor);
    if let Some(op) = op_id {
        for fragment in &fragments {
            name_band(store, op, fragment.tool_face, fragment.face);
        }
        if let Some(entry_face) = remapped.single_target_face() {
            name_rim(store, op, entry_face, entry_ring.wire, 0);
        }
        name_floor(store, op, buried.cap_face, floor);
    }
    Ok(())
}

/// Binds the pocket floor's `Floor { op, cap name }` when the buried cap is
/// named (unnamed tools propagate unnamed floors).
fn name_floor(
    store: &mut TopologyStore,
    op: &crate::topology::OpId,
    cap_face: FaceId,
    floor: FaceId,
) {
    if let Some(cap_name) = store.names().name_of_face(cap_face).cloned() {
        store.names_mut().bind_face(
            floor,
            crate::topology::FaceName::Floor {
                op: op.clone(),
                cap: Box::new(cap_name),
            },
        );
    }
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

/// Punches one cut loop onto the given RESULT face (a copy or a split
/// fragment) and returns the hole-ring [`WireId`] created on it (so the band
/// face can share it). Punching the preserved input face would attach the
/// ring to the wrong face.
fn punch_onto(
    store: &mut TopologyStore,
    loop_: &super::loops::CutLoop,
    result_face: FaceId,
) -> Result<crate::topology::WireId> {
    let mut remapped = loop_.clone();
    remapped.target_face = result_face;
    punch_loop(store, &remapped)
}

/// Names the band face `Band { op, tool_face name, 0 }` when the tool side
/// face is itself named (unnamed tools propagate unnamed bands).
fn name_band(
    store: &mut TopologyStore,
    op: &crate::topology::OpId,
    tool_face: FaceId,
    band: FaceId,
) {
    if let Some(tool_name) = store.names().name_of_face(tool_face).cloned() {
        store.names_mut().bind_face(
            band,
            crate::topology::FaceName::Band {
                op: op.clone(),
                tool_face: Box::new(tool_name),
                loop_index: 0,
            },
        );
    }
}

/// Names a punched hole-rim ring edge `CutRim { op, punched face name, loop }`
/// when the punched target face is named. The name binds to the ring wire's
/// first edge (for a chained ring, its first chain segment's edge).
fn name_rim(
    store: &mut TopologyStore,
    op: &crate::topology::OpId,
    punched_face: FaceId,
    ring_wire: crate::topology::WireId,
    loop_index: u32,
) {
    let Some(target_name) = store.names().name_of_face(punched_face).cloned() else {
        return;
    };
    let Ok(wire) = store.wire(ring_wire) else {
        return;
    };
    let Some(rim_edge) = wire.edges.first().map(|oe| oe.edge) else {
        return;
    };
    store.names_mut().bind_edge(
        rim_edge,
        crate::topology::EdgeName::CutRim {
            op: op.clone(),
            target: Box::new(target_name),
            loop_index,
        },
    );
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
        // The copy references the same shared boundary edges, so the per-edge
        // UV images remain valid on the copy.
        pcurves: src.pcurves.clone(),
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
#[allow(clippy::unwrap_used, clippy::expect_used)]
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
        let result = subtract_through_cut(&mut store, slab, tube, None).unwrap();
        (store, result)
    }

    /// F4 acceptance 1 (stability): rebuilding the same named model into a
    /// fresh store resolves every persistent name to geometrically identical
    /// entities.
    #[test]
    fn boolean_names_are_rebuild_stable() {
        use crate::topology::{EdgeName, FaceName, FaceRole, OpId};

        let build = |center_x: f64| {
            let mut store = TopologyStore::new();
            let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
                .with_op_id(OpId::new("slab1"))
                .execute(&mut store)
                .unwrap();
            let tube = MakeNurbsTube::new(Point3::new(center_x, 3.0, -1.5), 0.7, 5.0)
                .with_op_id(OpId::new("win1"))
                .execute(&mut store)
                .unwrap();
            let result =
                subtract_through_cut(&mut store, slab, tube, Some(&OpId::new("cut1"))).unwrap();
            (store, result)
        };

        let slab_top = FaceName::Created {
            op: OpId::new("slab1"),
            role: FaceRole::Top,
        };
        let tool_side = FaceName::Created {
            op: OpId::new("win1"),
            role: FaceRole::Side(0),
        };
        let band = FaceName::Band {
            op: OpId::new("cut1"),
            tool_face: Box::new(tool_side.clone()),
            loop_index: 0,
        };
        let rim_entry = EdgeName::CutRim {
            op: OpId::new("cut1"),
            target: Box::new(FaceName::Created {
                op: OpId::new("slab1"),
                role: FaceRole::Bottom,
            }),
            loop_index: 0,
        };
        let rim_exit = EdgeName::CutRim {
            op: OpId::new("cut1"),
            target: Box::new(slab_top.clone()),
            loop_index: 1,
        };

        let (store_a, result_a) = build(3.0);
        let (store_b, _) = build(3.0);

        for name in [&slab_top, &band] {
            let fa = store_a.names().face(name).expect("A resolves");
            let fb = store_b.names().face(name).expect("B resolves");
            let sample = |store: &TopologyStore, f| match &store.face(f).unwrap().surface {
                FaceSurface::Nurbs(s) => s.point_at(0.31, 0.62).unwrap(),
                FaceSurface::Plane(p) => *p.origin(),
                _ => panic!("unexpected surface"),
            };
            assert!(
                (sample(&store_a, fa) - sample(&store_b, fb)).norm() < 1e-9,
                "{name:?} moved across rebuilds"
            );
        }
        assert!(
            store_a.names().edge(&rim_entry).is_some(),
            "entry rim named"
        );
        assert!(store_a.names().edge(&rim_exit).is_some(), "exit rim named");

        // The resolved slab top face belongs to the RESULT solid (name moved
        // off the input): it carries the punched hole.
        let top_face = store_a.names().face(&slab_top).unwrap();
        let shell = store_a
            .shell(store_a.solid(result_a).unwrap().outer_shell)
            .unwrap();
        assert!(shell.faces.contains(&top_face), "name resolves into result");
        assert!(
            !store_a.face(top_face).unwrap().inner_wires.is_empty(),
            "punched top face carries the hole"
        );
    }

    /// F4 acceptance 2 (parameter change): moving the window keeps the same
    /// names resolving to the same ROLES (the punched top face still resolves,
    /// with its hole at the new location).
    #[test]
    fn boolean_names_survive_parameter_changes() {
        use crate::topology::{FaceName, FaceRole, OpId};

        let build = |center_x: f64| {
            let mut store = TopologyStore::new();
            let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
                .with_op_id(OpId::new("slab1"))
                .execute(&mut store)
                .unwrap();
            let tube = MakeNurbsTube::new(Point3::new(center_x, 3.0, -1.5), 0.7, 5.0)
                .with_op_id(OpId::new("win1"))
                .execute(&mut store)
                .unwrap();
            subtract_through_cut(&mut store, slab, tube, Some(&OpId::new("cut1"))).unwrap();
            store
        };

        let slab_top = FaceName::Created {
            op: OpId::new("slab1"),
            role: FaceRole::Top,
        };
        let band = FaceName::Band {
            op: OpId::new("cut1"),
            tool_face: Box::new(FaceName::Created {
                op: OpId::new("win1"),
                role: FaceRole::Side(0),
            }),
            loop_index: 0,
        };

        let moved = build(3.5);
        let top = moved.names().face(&slab_top).expect("top still resolves");
        assert!(
            !moved.face(top).unwrap().inner_wires.is_empty(),
            "moved window still punches the top face"
        );
        let band_face = moved.names().face(&band).expect("band still resolves");
        // The band lies on the moved tool: its surface contains the new axis.
        let FaceSurface::Nurbs(surf) = &moved.face(band_face).unwrap().surface else {
            panic!("band must be NURBS");
        };
        let p = surf.point_at(0.0, 0.5).unwrap();
        let r = ((p.x - 3.5).powi(2) + (p.y - 3.0).powi(2)).sqrt();
        assert!((r - 0.7).abs() < 1e-6, "band follows the moved tool");
    }

    /// F4 acceptance 3 (topology change): shortening the tool turns the
    /// through cut into a pocket — the band keeps resolving, the exit rim
    /// stops resolving, and the floor appears.
    #[test]
    fn boolean_names_survive_through_to_pocket_transition() {
        use crate::topology::{EdgeName, FaceName, FaceRole, OpId};

        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .with_op_id(OpId::new("slab1"))
            .execute(&mut store)
            .unwrap();
        // Short tube: ends inside the slab (pocket).
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -3.0), 0.7, 3.5)
            .with_op_id(OpId::new("win1"))
            .execute(&mut store)
            .unwrap();
        subtract_through_cut(&mut store, slab, tube, Some(&OpId::new("cut1"))).unwrap();

        let tool_side = FaceName::Created {
            op: OpId::new("win1"),
            role: FaceRole::Side(0),
        };
        let band = FaceName::Band {
            op: OpId::new("cut1"),
            tool_face: Box::new(tool_side),
            loop_index: 0,
        };
        assert!(store.names().face(&band).is_some(), "pocket band resolves");

        let floor = FaceName::Floor {
            op: OpId::new("cut1"),
            cap: Box::new(FaceName::Created {
                op: OpId::new("win1"),
                role: FaceRole::CapEnd,
            }),
        };
        assert!(
            store.names().face(&floor).is_some(),
            "pocket floor resolves"
        );

        let rim_exit = EdgeName::CutRim {
            op: OpId::new("cut1"),
            target: Box::new(FaceName::Created {
                op: OpId::new("slab1"),
                role: FaceRole::Top,
            }),
            loop_index: 1,
        };
        assert!(
            store.names().edge(&rim_exit).is_none(),
            "the exit rim genuinely no longer exists"
        );
    }

    /// F3a pocket subtract: a tube entering the slab from below and ending
    /// INSIDE it cuts a blind pocket — entry hole on the back face, a band
    /// down the tube wall, and the buried tool cap (sense-flipped) as the
    /// pocket floor.
    #[test]
    fn half_buried_tube_cuts_a_pocket() {
        use crate::topology::FaceSurface;

        let mut store = TopologyStore::new();
        // Slab: front z in [0, 1.5] (peak at center), back = front - 1.0.
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        // Tube from z = -3 up to z = 0.5: crosses the back face (z ≈ 0 at the
        // tube footprint) once and ends inside the slab (front ≈ 0.98 there).
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -3.0), 0.7, 3.5)
            .execute(&mut store)
            .unwrap();

        let result = subtract_through_cut(&mut store, slab, tube, None).unwrap();
        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();

        // Exactly one punched face (the entry face carries the hole).
        let punched = shell
            .faces
            .iter()
            .filter(|&&f| !store.face(f).unwrap().inner_wires.is_empty())
            .count();
        assert!(punched >= 2, "entry face + band expected, got {punched}");

        // The pocket floor exists: a planar face at z ≈ 0.5 whose effective
        // outward normal points DOWN into the cavity (sense-flipped tool cap).
        let mut floor_found = false;
        for &f in &shell.faces {
            let face = store.face(f).unwrap();
            if let FaceSurface::Plane(plane) = &face.surface {
                let origin_z = plane.origin().z;
                if (origin_z - 0.5).abs() < 1e-6 {
                    let n = plane.plane_normal();
                    let effective_z = if face.same_sense { n.z } else { -n.z };
                    assert!(
                        effective_z < 0.0,
                        "pocket floor must face down into the cavity"
                    );
                    floor_found = true;
                }
            }
        }
        assert!(floor_found, "expected a planar pocket floor at z = 0.5");

        // The result tessellates edge-manifold and stays inside the slab's
        // z-range (the tube below the slab is discarded).
        let mesh = TessellateSolid::new(result, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
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
        let zmin = mesh
            .vertices
            .iter()
            .map(|p| p.z)
            .fold(f64::INFINITY, f64::min);
        assert!(
            zmin > -1.0 - 1e-6,
            "mesh reaches below the slab (zmin = {zmin}); the tube stub \
             outside the target must be discarded"
        );
    }

    /// The deferred F1 acceptance case: a revolved solid (closed wall, u/v seam
    /// on the wall surface) cut by a HORIZONTAL tube. Both the entry and exit
    /// holes land on the SAME closed wall face, and the tool's own periodic
    /// direction wraps during SSI. The tube runs along +Y so its holes sit at
    /// wall azimuths ±π/2, safely away from the wall's parametric seam at +X.
    #[test]
    fn revolved_solid_minus_horizontal_tube_is_manifold() {
        use crate::geometry::nurbs::NurbsCurve3D;
        use crate::math::Vector3;
        use crate::operations::creation::{MakeNurbsPrism, MakeRevolvedSolid};

        let mut store = TopologyStore::new();
        // Vase-like profile: wall radius 2.0-2.6 over height 0-3.6.
        let vase = MakeRevolvedSolid::new(vec![(2.0, 0.0), (2.4, 1.2), (2.1, 2.4), (2.6, 3.6)])
            .execute(&mut store)
            .unwrap();
        // Horizontal tube along +Y through both walls at mid-height.
        let circle =
            NurbsCurve3D::circle(Point3::new(0.0, -4.0, 1.8), 0.5, Vector3::y(), Vector3::x())
                .unwrap();
        let tube = MakeNurbsPrism::new(circle, Vector3::new(0.0, 8.0, 0.0))
            .execute(&mut store)
            .unwrap();

        let result = subtract_through_cut(&mut store, vase, tube, None).unwrap();

        // Entry and exit holes both land on the single closed wall face: one
        // result face carries exactly 2 hole inner wires (the band face carries
        // 1 — its exit ring).
        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();
        let two_hole_faces = shell
            .faces
            .iter()
            .filter(|&&f| store.face(f).unwrap().inner_wires.len() == 2)
            .count();
        assert_eq!(
            two_hole_faces, 1,
            "exactly one face (the revolved wall) carries both holes"
        );

        // The whole result tessellates edge-manifold: every undirected edge is
        // used by 1 or 2 triangles.
        let mesh = TessellateSolid::new(result, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty(), "empty result mesh");
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
    /// - slab − tube, before seam conformance: 264 boundary edges, 4 of them in
    ///   the hole-ring region at the SSI seam azimuth (the punch chord vs. band
    ///   vertical-stitch disagreement at the tool's u-seam).
    /// - slab − tube, WITH seam conformance: 0 hole-ring boundary edges. The
    ///   SSI marcher wraps the tool's periodic u direction and emits exact seam
    ///   samples shared by both the punch ring (`uv_a`) and the band ribbon
    ///   (`uv_b`), so the two sides conform across the seam and the band ribbon
    ///   spans the full tool u domain.
    ///
    /// Two assertions pin the result:
    /// 1. The cut result's total boundary-edge count is no worse than the plain
    ///    slab's own perimeter nonconformance (plus a small margin); the prior
    ///    ~1404 hole-ring boundary edges are gone.
    /// 2. Direct hole-ring conformance: NO boundary-edge midpoint lies in the
    ///    tube-wall ring region (distance to the tube axis within [0.7·r, 1.3·r]
    ///    while z is inside the slab). The marcher's exact seam samples close
    ///    the seam, so even the former seam-azimuth residual is gone.
    #[test]
    fn hole_rings_tessellate_conformally() {
        const RADIUS: f64 = 0.7;
        const MARGIN: usize = 16;

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

        // Since shared-edge topology (F2), the plain slab position-welds fully
        // watertight: per-face perimeters conform exactly, so position-dedup
        // leaves ZERO boundary edges.
        let mut plain_store = TopologyStore::new();
        let plain = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut plain_store)
            .unwrap();
        let plain_boundary = boundary_edges(&plain_store, plain).len();
        assert_eq!(
            plain_boundary, 0,
            "plain slab must position-weld watertight (found {plain_boundary} \
             boundary edges)"
        );

        let (store, result) = slab_minus_tube(RADIUS);
        let cut_edges = boundary_edges(&store, result);
        let cut_boundary = cut_edges.len();

        // (1) The cut result carries no more boundary edges than a small
        // margin above the (now zero) plain-slab baseline. The prior ~1404
        // hole-ring T-junctions are eliminated.
        assert!(
            cut_boundary <= plain_boundary + MARGIN,
            "cut result has {cut_boundary} boundary edges, expected \
             <= {plain_boundary} (plain baseline) + {MARGIN}; hole-ring \
             T-junctions appear to have returned"
        );

        // (2) Direct hole-ring conformance: NO boundary-edge midpoint lies in
        // the tube-wall ring region. The tube axis runs along (3,3,z); a ring
        // boundary edge would sit at radius ~RADIUS from that axis, inside the
        // slab body in z. The marcher's exact seam samples are shared by punch
        // and band, so even the former seam-azimuth residual (up to 4 edges)
        // is gone.
        let axis = Point3::new(3.0, 3.0, 0.0);
        let mut ring_edges = 0usize;
        for (p, q) in &cut_edges {
            let m = Point3::new((p.x + q.x) * 0.5, (p.y + q.y) * 0.5, (p.z + q.z) * 0.5);
            let dxy = ((m.x - axis.x).powi(2) + (m.y - axis.y).powi(2)).sqrt();
            let in_ring_radius = (0.7 * RADIUS..=1.3 * RADIUS).contains(&dxy);
            let in_slab_z = m.z > -1.2 && m.z < 1.7;
            if in_ring_radius && in_slab_z {
                ring_edges += 1;
            }
        }
        assert_eq!(
            ring_edges, 0,
            "expected 0 hole-ring boundary edges with marcher seam conformance, \
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

    // ---- F5 Phase B acceptance: multi-face (box) cutters ----

    /// Straight segmented-prism wall: 6 x 0.4 footprint extruded 3 up.
    fn straight_wall_profile() -> Vec<crate::operations::creation::ProfileSegment> {
        use crate::operations::creation::ProfileSegment;
        let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
        let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
        vec![
            line(p(0.0, 0.0), p(6.0, 0.0)),
            line(p(6.0, 0.0), p(6.0, 0.4)),
            line(p(6.0, 0.4), p(0.0, 0.4)),
            line(p(0.0, 0.4), p(0.0, 0.0)),
        ]
    }

    /// Box-cutter window profile in the XZ plane at `y = -1`, spanning
    /// `x in [x0, x0 + 1.5]`, `z in [1, 2]`.
    fn box_window_profile(x0: f64) -> Vec<crate::operations::creation::ProfileSegment> {
        use crate::operations::creation::ProfileSegment;
        let p = |x: f64, z: f64| Point3::new(x, -1.0, z);
        let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
        let x1 = x0 + 1.5;
        vec![
            line(p(x0, 1.0), p(x1, 1.0)), // sill
            line(p(x1, 1.0), p(x1, 2.0)), // jamb-right
            line(p(x1, 2.0), p(x0, 2.0)), // head
            line(p(x0, 2.0), p(x0, 1.0)), // jamb-left
        ]
    }

    fn wall_tags() -> Vec<crate::topology::SegmentTag> {
        ["outer", "end-east", "inner", "end-west"]
            .iter()
            .map(|t| crate::topology::SegmentTag::new(*t))
            .collect()
    }

    fn box_tags() -> Vec<crate::topology::SegmentTag> {
        ["sill", "jamb-right", "head", "jamb-left"]
            .iter()
            .map(|t| crate::topology::SegmentTag::new(*t))
            .collect()
    }

    /// Builds wall − box window with op ids; `x0` positions the window,
    /// `depth` the box extrusion length from `y = -1` (2.4 = through,
    /// 1.25 = buried mid-wall).
    fn named_wall_minus_box(x0: f64, depth: f64) -> (TopologyStore, SolidId) {
        use crate::math::Vector3;
        use crate::operations::creation::MakeSegmentedPrism;
        use crate::topology::OpId;

        let mut store = TopologyStore::new();
        let wall = MakeSegmentedPrism::new(straight_wall_profile(), Vector3::new(0.0, 0.0, 3.0))
            .with_op_id(OpId::new("wall1"))
            .with_segment_tags(wall_tags())
            .execute(&mut store)
            .unwrap();
        let cutter = MakeSegmentedPrism::new(box_window_profile(x0), Vector3::new(0.0, depth, 0.0))
            .with_op_id(OpId::new("win1"))
            .with_segment_tags(box_tags())
            .execute(&mut store)
            .unwrap();
        let result =
            subtract_through_cut(&mut store, wall, cutter, Some(&OpId::new("cut1"))).unwrap();
        (store, result)
    }

    use super::super::test_support::welded_boundary_edges;

    /// Acceptance B1: segmented-prism wall − genuine 4-face box cutter → a
    /// through window. The result is position-weld watertight, the hole is
    /// genuinely open, both punched wall faces carry a 4-edge chained hole
    /// ring, and the band consists of exactly 4 fragment faces.
    #[test]
    fn wall_minus_box_window_is_watertight_with_open_hole() {
        use crate::math::Vector3;
        use crate::operations::creation::MakeSegmentedPrism;

        let mut store = TopologyStore::new();
        let wall = MakeSegmentedPrism::new(straight_wall_profile(), Vector3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();
        let cutter = MakeSegmentedPrism::new(box_window_profile(2.0), Vector3::new(0.0, 2.4, 0.0))
            .execute(&mut store)
            .unwrap();
        let result = subtract_through_cut(&mut store, wall, cutter, None).unwrap();

        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();
        assert_eq!(
            shell.faces.len(),
            10,
            "6 wall faces + 4 band fragments, got {}",
            shell.faces.len()
        );

        // The two punched wall faces each carry ONE hole ring wire of 4 edges
        // (one per chain segment).
        let mut punched = 0usize;
        for &f in &shell.faces {
            let face = store.face(f).unwrap();
            if face.inner_wires.is_empty() {
                continue;
            }
            punched += 1;
            assert_eq!(face.inner_wires.len(), 1, "one hole per punched face");
            let ring = store.wire(face.inner_wires[0]).unwrap();
            assert_eq!(ring.edges.len(), 4, "chained hole ring has 4 edges");
        }
        assert_eq!(punched, 2, "entry + exit faces punched");

        // Position-weld watertight: no boundary edges anywhere — the chained
        // rings, the kink crossings, and the wall perimeter all conform.
        let boundary = welded_boundary_edges(&store, result);
        assert_eq!(
            boundary, 0,
            "wall − box window must position-weld watertight \
             (found {boundary} boundary edges)"
        );

        // The hole is open: no mesh vertex intrudes into the tunnel interior.
        let mesh = TessellateSolid::new(result, TessellationParams::default())
            .execute(&store)
            .unwrap();
        for v in &mesh.vertices {
            let inside = v.x > 2.05 && v.x < 3.45 && v.z > 1.05 && v.z < 1.95;
            assert!(
                !(inside && v.y > 0.05 && v.y < 0.35),
                "vertex ({:.3},{:.3},{:.3}) intrudes into the window tunnel",
                v.x,
                v.y,
                v.z
            );
        }
    }

    /// Adjacent band fragments share their kink-crossing edges, and each
    /// fragment shares its entry/exit ring edge with the punched wall faces'
    /// chained hole rings (F2 shared-edge topology).
    #[test]
    fn box_window_band_fragments_share_ring_and_kink_edges() {
        use crate::math::Vector3;
        use crate::operations::creation::MakeSegmentedPrism;
        use crate::topology::EdgeId;

        let mut store = TopologyStore::new();
        let wall = MakeSegmentedPrism::new(straight_wall_profile(), Vector3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();
        let cutter = MakeSegmentedPrism::new(box_window_profile(2.0), Vector3::new(0.0, 2.4, 0.0))
            .execute(&mut store)
            .unwrap();
        let result = subtract_through_cut(&mut store, wall, cutter, None).unwrap();
        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();

        // Ring edges of the two punched faces.
        let mut ring_edges: Vec<EdgeId> = Vec::new();
        for &f in &shell.faces {
            for &w in &store.face(f).unwrap().inner_wires {
                ring_edges.extend(store.wire(w).unwrap().edges.iter().map(|oe| oe.edge));
            }
        }
        assert_eq!(ring_edges.len(), 8, "2 chained rings x 4 edges");

        // Band fragments: faces with a 4-edge outer wire referencing ring
        // edges (the wall faces' outer wires reference wall boundary edges).
        let mut kink_edge_uses: HashMap<EdgeId, usize> = HashMap::new();
        let mut fragments = 0usize;
        for &f in &shell.faces {
            let face = store.face(f).unwrap();
            let wire = store.wire(face.outer_wire).unwrap();
            let ring_refs = wire
                .edges
                .iter()
                .filter(|oe| ring_edges.contains(&oe.edge))
                .count();
            if ring_refs == 0 {
                continue;
            }
            fragments += 1;
            assert_eq!(wire.edges.len(), 4, "fragment wire: 2 rings + 2 kinks");
            assert_eq!(
                ring_refs, 2,
                "fragment references one entry + one exit ring edge"
            );
            for oe in &wire.edges {
                if !ring_edges.contains(&oe.edge) {
                    *kink_edge_uses.entry(oe.edge).or_insert(0) += 1;
                }
            }
        }
        assert_eq!(fragments, 4, "one band fragment per box side face");
        assert_eq!(kink_edge_uses.len(), 4, "4 shared kink-crossing edges");
        for (&edge, &uses) in &kink_edge_uses {
            assert_eq!(
                uses, 2,
                "kink edge {edge:?} must be shared by exactly 2 fragments"
            );
        }
    }

    /// Acceptance B2 (F4 rebuild-stability pattern): every band fragment binds
    /// `Band {{ op, tool_face: that side face's Created{{Tagged}} name }}`, the
    /// chained rims bind `CutRim`, and rebuilding the same model into a fresh
    /// store resolves every name to identical geometry.
    #[test]
    fn box_window_band_names_are_rebuild_stable() {
        use crate::topology::{EdgeName, FaceName, FaceRole, OpId, SegmentTag};

        let (store_a, result_a) = named_wall_minus_box(2.0, 2.4);
        let (store_b, _) = named_wall_minus_box(2.0, 2.4);

        let band_name = |tag: &str| FaceName::Band {
            op: OpId::new("cut1"),
            tool_face: Box::new(FaceName::Created {
                op: OpId::new("win1"),
                role: FaceRole::Tagged(SegmentTag::new(tag)),
            }),
            loop_index: 0,
        };

        let shell_a = store_a
            .shell(store_a.solid(result_a).unwrap().outer_shell)
            .unwrap();
        for tag in ["sill", "jamb-right", "head", "jamb-left"] {
            let name = band_name(tag);
            let fa = store_a.names().face(&name).expect("band resolves in A");
            let fb = store_b.names().face(&name).expect("band resolves in B");
            assert!(
                shell_a.faces.contains(&fa),
                "named band fragment lives in the result shell"
            );
            let sample = |store: &TopologyStore, f| match &store.face(f).unwrap().surface {
                FaceSurface::Nurbs(s) => s.point_at(0.4, 0.3).unwrap(),
                other => panic!("band fragment must be NURBS, got {other:?}"),
            };
            assert!(
                (sample(&store_a, fa) - sample(&store_b, fb)).norm() < 1e-9,
                "band fragment {tag} moved across rebuilds"
            );
        }

        // Chained hole rims: entry on the outer wall face, exit on the inner.
        for (target_tag, loop_index) in [("outer", 0u32), ("inner", 1u32)] {
            let rim = EdgeName::CutRim {
                op: OpId::new("cut1"),
                target: Box::new(FaceName::Created {
                    op: OpId::new("wall1"),
                    role: FaceRole::Tagged(SegmentTag::new(target_tag)),
                }),
                loop_index,
            };
            assert!(
                store_a.names().edge(&rim).is_some(),
                "{target_tag} rim named in A"
            );
            assert!(
                store_b.names().edge(&rim).is_some(),
                "{target_tag} rim named in B"
            );
        }

        // The wall's tagged faces transferred onto the punched result copies.
        let outer = FaceName::Created {
            op: OpId::new("wall1"),
            role: FaceRole::Tagged(SegmentTag::new("outer")),
        };
        let outer_face = store_a.names().face(&outer).expect("outer resolves");
        assert!(shell_a.faces.contains(&outer_face));
        assert!(
            !store_a.face(outer_face).unwrap().inner_wires.is_empty(),
            "punched outer face carries the chained hole ring"
        );
    }

    /// Acceptance B2 (parameter change): moving the box window keeps every
    /// canonical name resolving to the same ROLE at the new location.
    #[test]
    fn box_window_names_survive_parameter_changes() {
        use crate::topology::{FaceName, FaceRole, OpId, SegmentTag};

        let (moved, _) = named_wall_minus_box(2.3, 2.4);

        let sill_band = FaceName::Band {
            op: OpId::new("cut1"),
            tool_face: Box::new(FaceName::Created {
                op: OpId::new("win1"),
                role: FaceRole::Tagged(SegmentTag::new("sill")),
            }),
            loop_index: 0,
        };
        let band_face = moved.names().face(&sill_band).expect("sill band resolves");
        let FaceSurface::Nurbs(surf) = &moved.face(band_face).unwrap().surface else {
            panic!("band fragment must be NURBS");
        };
        // The sill fragment follows the moved tool: it stays the z = 1 plane
        // strip, now spanning x in [2.3, 3.8].
        let p = surf.point_at(0.5, 0.5).unwrap();
        assert!((p.z - 1.0).abs() < 1e-9, "sill band stays at z = 1");
        assert!(
            p.x > 2.3 - 1e-9 && p.x < 3.8 + 1e-9,
            "sill band moved with the box"
        );

        let outer = FaceName::Created {
            op: OpId::new("wall1"),
            role: FaceRole::Tagged(SegmentTag::new("outer")),
        };
        let outer_face = moved.names().face(&outer).expect("outer still resolves");
        assert!(
            !moved.face(outer_face).unwrap().inner_wires.is_empty(),
            "moved window still punches the outer face"
        );
    }

    /// Acceptance B3 (pocket variant): a box buried mid-wall cuts a blind
    /// niche — 4 band fragments down to the shared buried ring, the flipped
    /// buried cap as the floor (named `Floor`), and NO exit rim.
    #[test]
    fn buried_box_cuts_a_pocket_with_fragment_band_and_floor() {
        use crate::topology::{EdgeName, FaceName, FaceRole, OpId, SegmentTag};

        // Box from y = -1 to y = 0.25: enters the wall (y in [0, 0.4]) and
        // ends inside it.
        let (store, result) = named_wall_minus_box(2.0, 1.25);
        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();
        assert_eq!(
            shell.faces.len(),
            11,
            "6 wall faces + 4 band fragments + floor, got {}",
            shell.faces.len()
        );

        // The floor: a planar face at y = 0.25 whose effective normal points
        // back INTO the cavity (-Y).
        let mut floor_found = false;
        for &f in &shell.faces {
            let face = store.face(f).unwrap();
            if let FaceSurface::Plane(plane) = &face.surface {
                if (plane.origin().y - 0.25).abs() < 1e-9 {
                    let n = plane.plane_normal();
                    let effective_y = if face.same_sense { n.y } else { -n.y };
                    assert!(effective_y < 0.0, "pocket floor must face into the cavity");
                    floor_found = true;
                }
            }
        }
        assert!(floor_found, "expected a planar pocket floor at y = 0.25");

        // Names: all 4 fragment bands + the floor resolve; the exit rim does
        // not exist.
        for tag in ["sill", "jamb-right", "head", "jamb-left"] {
            let band = FaceName::Band {
                op: OpId::new("cut1"),
                tool_face: Box::new(FaceName::Created {
                    op: OpId::new("win1"),
                    role: FaceRole::Tagged(SegmentTag::new(tag)),
                }),
                loop_index: 0,
            };
            assert!(
                store.names().face(&band).is_some(),
                "pocket band fragment {tag} resolves"
            );
        }
        let floor = FaceName::Floor {
            op: OpId::new("cut1"),
            cap: Box::new(FaceName::Created {
                op: OpId::new("win1"),
                role: FaceRole::CapEnd,
            }),
        };
        assert!(
            store.names().face(&floor).is_some(),
            "pocket floor resolves"
        );
        let exit_rim = EdgeName::CutRim {
            op: OpId::new("cut1"),
            target: Box::new(FaceName::Created {
                op: OpId::new("wall1"),
                role: FaceRole::Tagged(SegmentTag::new("inner")),
            }),
            loop_index: 1,
        };
        assert!(
            store.names().edge(&exit_rim).is_none(),
            "the exit rim genuinely no longer exists"
        );

        // The pocket result position-welds watertight.
        let boundary = welded_boundary_edges(&store, result);
        assert_eq!(
            boundary, 0,
            "buried-box pocket must position-weld watertight \
             (found {boundary} boundary edges)"
        );
    }

    /// Acceptance B4 (curved host): an annular segmented-prism wall (two arc
    /// side faces + two radial end faces) cut by a radial box — the chained
    /// entry/exit loops land on the CURVED arc faces, the punched cylindrical
    /// faces stay conformal, and all 4 fragment bands resolve by name.
    #[test]
    fn curved_wall_box_window_through_arc_faces() {
        use crate::math::Vector3;
        use crate::operations::creation::{MakeSegmentedPrism, ProfileSegment};
        use crate::topology::{FaceName, FaceRole, OpId, SegmentTag};
        use std::f64::consts::PI;

        let deg = |d: f64| d * PI / 180.0;
        let mut store = TopologyStore::new();

        // Annular wall strip: outer arc r = 8.4, inner arc r = 8.0, azimuth
        // 60..120 degrees, extruded 2.5 up. The inner arc is traversed
        // backwards via the -Z normal (the Phase A fillet convention).
        let outer_start = Point3::new(8.4 * deg(60.0).cos(), 8.4 * deg(60.0).sin(), 0.0);
        let outer_end = Point3::new(8.4 * deg(120.0).cos(), 8.4 * deg(120.0).sin(), 0.0);
        let inner_start = Point3::new(8.0 * deg(120.0).cos(), 8.0 * deg(120.0).sin(), 0.0);
        let inner_end = Point3::new(8.0 * deg(60.0).cos(), 8.0 * deg(60.0).sin(), 0.0);
        let profile = vec![
            ProfileSegment::Arc {
                center: Point3::origin(),
                radius: 8.4,
                normal: Vector3::z(),
                ref_dir: Vector3::x(),
                start_angle: deg(60.0),
                end_angle: deg(120.0),
            },
            ProfileSegment::Line {
                start: outer_end,
                end: inner_start,
            },
            ProfileSegment::Arc {
                center: Point3::origin(),
                radius: 8.0,
                normal: -Vector3::z(),
                ref_dir: Vector3::x(),
                start_angle: deg(-120.0),
                end_angle: deg(-60.0),
            },
            ProfileSegment::Line {
                start: inner_end,
                end: outer_start,
            },
        ];
        let wall_tags: Vec<SegmentTag> = ["convex", "end-west", "concave", "end-east"]
            .iter()
            .map(|t| SegmentTag::new(*t))
            .collect();
        let wall = MakeSegmentedPrism::new(profile, Vector3::new(0.0, 0.0, 2.5))
            .with_op_id(OpId::new("wall1"))
            .with_segment_tags(wall_tags)
            .execute(&mut store)
            .unwrap();

        // Radial box cutter through the wall at 90 degrees azimuth: window
        // rectangle in the XZ plane at y = 6.5, extruded 2.6 along +Y.
        let p = |x: f64, z: f64| Point3::new(x, 6.5, z);
        let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
        let cutter_profile = vec![
            line(p(-0.7, 0.9), p(0.7, 0.9)),
            line(p(0.7, 0.9), p(0.7, 1.7)),
            line(p(0.7, 1.7), p(-0.7, 1.7)),
            line(p(-0.7, 1.7), p(-0.7, 0.9)),
        ];
        let cutter = MakeSegmentedPrism::new(cutter_profile, Vector3::new(0.0, 2.6, 0.0))
            .with_op_id(OpId::new("win1"))
            .with_segment_tags(box_tags())
            .execute(&mut store)
            .unwrap();

        let result =
            subtract_through_cut(&mut store, wall, cutter, Some(&OpId::new("cut1"))).unwrap();

        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();
        assert_eq!(shell.faces.len(), 10, "6 wall faces + 4 band fragments");

        // Both punched faces are the CURVED arc faces (cylindrical patches).
        for &f in &shell.faces {
            let face = store.face(f).unwrap();
            if face.inner_wires.is_empty() {
                continue;
            }
            let FaceSurface::Nurbs(surf) = &face.surface else {
                panic!("punched face must be NURBS");
            };
            let sample = surf.point_at(0.31, 0.62).unwrap();
            let r = (sample.x * sample.x + sample.y * sample.y).sqrt();
            assert!(
                (r - 8.0).abs() < 1e-9 || (r - 8.4).abs() < 1e-9,
                "punched face is not one of the arc walls (r = {r})"
            );
        }

        // All 4 fragment bands resolve by name.
        for tag in ["sill", "jamb-right", "head", "jamb-left"] {
            let band = FaceName::Band {
                op: OpId::new("cut1"),
                tool_face: Box::new(FaceName::Created {
                    op: OpId::new("win1"),
                    role: FaceRole::Tagged(SegmentTag::new(tag)),
                }),
                loop_index: 0,
            };
            assert!(
                store.names().face(&band).is_some(),
                "curved-host band fragment {tag} resolves"
            );
        }

        // Manifold and watertight under position welding.
        let boundary = welded_boundary_edges(&store, result);
        assert_eq!(
            boundary, 0,
            "curved wall − box window must position-weld watertight \
             (found {boundary} boundary edges)"
        );
    }

    // ---- F5 Phase C acceptance: window across a target kink edge ----

    /// Straight wall whose OUTER side is segmented at x = 3 (two collinear
    /// pieces sharing a vertical target kink edge).
    fn segmented_outer_wall_profile() -> Vec<crate::operations::creation::ProfileSegment> {
        use crate::operations::creation::ProfileSegment;
        let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
        let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
        vec![
            line(p(0.0, 0.0), p(3.0, 0.0)), // outer-a
            line(p(3.0, 0.0), p(6.0, 0.0)), // outer-b
            line(p(6.0, 0.0), p(6.0, 0.4)), // end-east
            line(p(6.0, 0.4), p(0.0, 0.4)), // inner
            line(p(0.0, 0.4), p(0.0, 0.0)), // end-west
        ]
    }

    fn segmented_outer_wall_tags() -> Vec<crate::topology::SegmentTag> {
        ["outer-a", "outer-b", "end-east", "inner", "end-west"]
            .iter()
            .map(|t| crate::topology::SegmentTag::new(*t))
            .collect()
    }

    /// Builds (segmented-outer wall) − (box window at `x0`, straddling the
    /// x = 3 joint when `x0 = 2.0`) with op ids.
    fn named_kink_wall_minus_box(x0: f64) -> (TopologyStore, SolidId) {
        use crate::math::Vector3;
        use crate::operations::creation::MakeSegmentedPrism;
        use crate::topology::OpId;

        let mut store = TopologyStore::new();
        let wall =
            MakeSegmentedPrism::new(segmented_outer_wall_profile(), Vector3::new(0.0, 0.0, 3.0))
                .with_op_id(OpId::new("wall1"))
                .with_segment_tags(segmented_outer_wall_tags())
                .execute(&mut store)
                .unwrap();
        let cutter = MakeSegmentedPrism::new(box_window_profile(x0), Vector3::new(0.0, 2.4, 0.0))
            .with_op_id(OpId::new("win1"))
            .with_segment_tags(box_tags())
            .execute(&mut store)
            .unwrap();
        let result =
            subtract_through_cut(&mut store, wall, cutter, Some(&OpId::new("cut1"))).unwrap();
        (store, result)
    }

    /// Acceptance C1: a box window straddling a target kink edge splits both
    /// outer wall faces along the trace (boundary notches), shares the split
    /// kink sub-edges between the fragments, keeps the hole genuinely open,
    /// and position-welds watertight.
    #[test]
    fn window_across_target_kink_is_watertight_with_open_hole() {
        use crate::topology::EdgeId;

        let (store, result) = named_kink_wall_minus_box(2.0);
        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();
        // 5 unaffected faces (copies) + 2 notched fragments + 4 band
        // fragments.
        assert_eq!(shell.faces.len(), 11, "got {}", shell.faces.len());

        // The exit face (single inner face) carries the chained interior
        // hole ring; the two notched outer fragments carry NO inner wires
        // (their hole halves live in the outer trim notch).
        let punched: Vec<_> = shell
            .faces
            .iter()
            .filter(|&&f| !store.face(f).unwrap().inner_wires.is_empty())
            .collect();
        assert_eq!(punched.len(), 1, "only the inner face is punched");

        // The split kink edge: its two kept sub-edges (below / above the
        // window) are each shared by exactly the two outer fragments.
        let mut edge_face_uses: HashMap<EdgeId, usize> = HashMap::new();
        for &f in &shell.faces {
            let face = store.face(f).unwrap();
            for oe in &store.wire(face.outer_wire).unwrap().edges {
                *edge_face_uses.entry(oe.edge).or_insert(0) += 1;
            }
        }
        // Kink sub-edges: vertical segments at x = 3, y = 0, spanning
        // z in [0,1] (below) and [2,3] (above).
        let mut kink_subs = 0usize;
        for (&edge, &uses) in &edge_face_uses {
            let data = store.edge(edge).unwrap();
            let a = store.vertex(data.start).unwrap().point;
            let b = store.vertex(data.end).unwrap().point;
            let on_kink = (a.x - 3.0).abs() < 1e-9
                && (b.x - 3.0).abs() < 1e-9
                && a.y.abs() < 1e-9
                && b.y.abs() < 1e-9;
            if on_kink {
                kink_subs += 1;
                assert_eq!(uses, 2, "kink sub-edge must be shared by 2 fragments");
            }
        }
        assert_eq!(kink_subs, 2, "below + above kink sub-edges");

        // Position-weld watertight.
        let boundary = welded_boundary_edges(&store, result);
        assert_eq!(
            boundary, 0,
            "window across the kink must position-weld watertight \
             (found {boundary} boundary edges)"
        );

        // The hole is open: no vertex intrudes into the window tunnel.
        let mesh = TessellateSolid::new(result, TessellationParams::default())
            .execute(&store)
            .unwrap();
        for v in &mesh.vertices {
            let inside = v.x > 2.05 && v.x < 3.45 && v.z > 1.05 && v.z < 1.95;
            assert!(
                !(inside && v.y > 0.05 && v.y < 0.35),
                "vertex ({:.3},{:.3},{:.3}) intrudes into the window tunnel",
                v.x,
                v.y,
                v.z
            );
        }
    }

    /// Acceptance C1 (names): the notched outer faces KEEP their tagged
    /// names (one kept fragment = the parent face), the rims and band
    /// fragments resolve, and everything re-resolves identically across two
    /// from-scratch builds and after a parameter change.
    #[test]
    fn window_across_kink_names_are_rebuild_stable() {
        use crate::topology::{EdgeName, FaceName, FaceRole, OpId, SegmentTag};

        let (store_a, result_a) = named_kink_wall_minus_box(2.0);
        let (store_b, _) = named_kink_wall_minus_box(2.0);

        let outer_name = |tag: &str| FaceName::Created {
            op: OpId::new("wall1"),
            role: FaceRole::Tagged(SegmentTag::new(tag)),
        };
        let shell_a = store_a
            .shell(store_a.solid(result_a).unwrap().outer_shell)
            .unwrap();
        for tag in ["outer-a", "outer-b"] {
            let fa = store_a.names().face(&outer_name(tag)).expect("resolves A");
            let fb = store_b.names().face(&outer_name(tag)).expect("resolves B");
            assert!(
                shell_a.faces.contains(&fa),
                "notched {tag} is in the result"
            );
            let sample = |store: &TopologyStore, f| match &store.face(f).unwrap().surface {
                FaceSurface::Nurbs(s) => s.point_at(0.1, 0.9).unwrap(),
                other => panic!("outer face must be NURBS, got {other:?}"),
            };
            assert!(
                (sample(&store_a, fa) - sample(&store_b, fb)).norm() < 1e-9,
                "{tag} moved across rebuilds"
            );
        }

        // Rims: one CutRim per notched fragment (entry side, loop 0),
        // composed from the fragment's kept name.
        for tag in ["outer-a", "outer-b"] {
            let rim = EdgeName::CutRim {
                op: OpId::new("cut1"),
                target: Box::new(outer_name(tag)),
                loop_index: 0,
            };
            assert!(store_a.names().edge(&rim).is_some(), "{tag} rim resolves");
        }

        // All 4 band fragments resolve.
        for tag in ["sill", "jamb-right", "head", "jamb-left"] {
            let band = FaceName::Band {
                op: OpId::new("cut1"),
                tool_face: Box::new(FaceName::Created {
                    op: OpId::new("win1"),
                    role: FaceRole::Tagged(SegmentTag::new(tag)),
                }),
                loop_index: 0,
            };
            assert!(
                store_a.names().face(&band).is_some(),
                "band fragment {tag} resolves"
            );
        }

        // Parameter change: sliding the window keeps everything resolving.
        let (moved, _) = named_kink_wall_minus_box(2.2);
        for tag in ["outer-a", "outer-b"] {
            assert!(
                moved.names().face(&outer_name(tag)).is_some(),
                "moved window: {tag} still resolves"
            );
        }
    }

    /// Acceptance C3 (the flipped F1 seam guard): a revolved (closed) wall
    /// cut by a tube ACROSS its own parametric seam now succeeds — the
    /// seam-straddling hole is applied as two boundary notches on the
    /// unrolled wall face — and the result is manifold with the hole
    /// genuinely open.
    #[test]
    fn revolved_wall_cut_across_its_seam_is_manifold_with_hole() {
        use crate::geometry::nurbs::NurbsCurve3D;
        use crate::math::Vector3;
        use crate::operations::creation::{MakeNurbsPrism, MakeRevolvedSolid};
        use crate::topology::{FaceName, FaceRole, OpId};

        let build = || -> (TopologyStore, SolidId) {
            let mut store = TopologyStore::new();
            // Plain cylindrical revolved wall (closed in u, seam at +X).
            let vase = MakeRevolvedSolid::new(vec![(2.0, 0.0), (2.0, 3.0)])
                .with_op_id(OpId::new("vase1"))
                .execute(&mut store)
                .unwrap();
            // Tube along +X through both walls: the +X hole straddles the
            // wall's parametric seam; the -X hole is a plain interior hole.
            let circle =
                NurbsCurve3D::circle(Point3::new(-4.0, 0.0, 1.5), 0.4, Vector3::x(), Vector3::y())
                    .unwrap();
            let tube = MakeNurbsPrism::new(circle, Vector3::new(8.0, 0.0, 0.0))
                .with_op_id(OpId::new("win1"))
                .execute(&mut store)
                .unwrap();
            let result =
                subtract_through_cut(&mut store, vase, tube, Some(&OpId::new("cut1"))).unwrap();
            (store, result)
        };
        let (store, result) = build();

        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();
        // 2 caps + the notched wall fragment + 1 band = 4 faces.
        assert_eq!(shell.faces.len(), 4, "got {}", shell.faces.len());

        // The wall keeps its persistent name (one kept fragment = the wall
        // face itself, now notched at the seam) and carries the interior
        // hole as an inner wire while the seam hole lives in the outer trim.
        let wall_name = FaceName::Created {
            op: OpId::new("vase1"),
            role: FaceRole::Wall,
        };
        let wall = store.names().face(&wall_name).expect("wall resolves");
        assert!(shell.faces.contains(&wall), "named wall is in the result");
        let wall_face = store.face(wall).unwrap();
        assert_eq!(
            wall_face.inner_wires.len(),
            1,
            "interior (-X) hole is an inner wire; the seam (+X) hole is a \
             boundary notch"
        );
        assert_eq!(
            wall_face.trim.as_ref().unwrap().holes.len(),
            1,
            "one interior trim hole"
        );

        // Manifold: the result position-welds watertight (the unrolled
        // wall's seam sides weld against each other; the notches weld
        // against the band).
        let boundary = welded_boundary_edges(&store, result);
        assert_eq!(
            boundary, 0,
            "seam-straddling cut must position-weld watertight \
             (found {boundary} boundary edges)"
        );

        // The hole is genuinely open at the seam: no mesh vertex intrudes
        // into the tube interior around the +X wall crossing.
        let mesh = TessellateSolid::new(result, TessellationParams::default())
            .execute(&store)
            .unwrap();
        for v in &mesh.vertices {
            if v.x > 1.6 && v.x < 2.4 {
                let d = (v.y * v.y + (v.z - 1.5) * (v.z - 1.5)).sqrt();
                assert!(
                    d > 0.32,
                    "vertex ({:.3},{:.3},{:.3}) intrudes into the seam hole",
                    v.x,
                    v.y,
                    v.z
                );
            }
        }

        // Rebuild stability: the notched wall resolves to identical geometry
        // in a from-scratch rebuild.
        let (store_b, _) = build();
        let wall_b = store_b.names().face(&wall_name).expect("wall resolves B");
        let sample = |store: &TopologyStore, f: FaceId| match &store.face(f).unwrap().surface {
            FaceSurface::Nurbs(s) => s.point_at(0.31, 0.62).unwrap(),
            other => panic!("wall must be NURBS, got {other:?}"),
        };
        assert!(
            (sample(&store, wall) - sample(&store_b, wall_b)).norm() < 1e-9,
            "wall moved across rebuilds"
        );
    }

    /// Builds a 3-outer-segment wall (narrow tagged `outer-mid` piece,
    /// x in [2.5, 3.0]) minus a box window at `x0` (window spans
    /// `[x0, x0 + 1.5]`, covering `outer-mid` fully for `x0 in [1.5, 2.5]`).
    fn named_middle_segment_wall_minus_box(x0: f64) -> (TopologyStore, SolidId) {
        use crate::math::Vector3;
        use crate::operations::creation::{MakeSegmentedPrism, ProfileSegment};
        use crate::topology::{OpId, SegmentTag};

        let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
        let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
        let mut store = TopologyStore::new();
        let profile = vec![
            line(p(0.0, 0.0), p(2.5, 0.0)), // outer-a
            line(p(2.5, 0.0), p(3.0, 0.0)), // outer-mid (narrow)
            line(p(3.0, 0.0), p(6.0, 0.0)), // outer-c
            line(p(6.0, 0.0), p(6.0, 0.4)),
            line(p(6.0, 0.4), p(0.0, 0.4)),
            line(p(0.0, 0.4), p(0.0, 0.0)),
        ];
        let tags: Vec<SegmentTag> = [
            "outer-a",
            "outer-mid",
            "outer-c",
            "end-east",
            "inner",
            "end-west",
        ]
        .iter()
        .map(|t| SegmentTag::new(*t))
        .collect();
        let wall = MakeSegmentedPrism::new(profile, Vector3::new(0.0, 0.0, 3.0))
            .with_op_id(OpId::new("wall1"))
            .with_segment_tags(tags)
            .execute(&mut store)
            .unwrap();
        let cutter = MakeSegmentedPrism::new(box_window_profile(x0), Vector3::new(0.0, 2.4, 0.0))
            .with_op_id(OpId::new("win1"))
            .with_segment_tags(box_tags())
            .execute(&mut store)
            .unwrap();
        let result =
            subtract_through_cut(&mut store, wall, cutter, Some(&OpId::new("cut1"))).unwrap();
        (store, result)
    }

    /// Acceptance C2: a window spanning a full (narrow) middle segment
    /// severs that face into TWO kept fragments named `Split{{l|r}}`; the
    /// names re-resolve across rebuilds and after sliding the window.
    #[test]
    fn window_spanning_middle_segment_binds_split_names() {
        use crate::topology::{FaceName, FaceRole, OpId, SegmentTag, SplitSide};

        let build = named_middle_segment_wall_minus_box;

        // Window x in [2.0, 3.5] covers outer-mid (x in [2.5, 3.0]) fully.
        let (store_a, result_a) = build(2.0);
        let (store_b, _) = build(2.0);

        let parent = FaceName::Created {
            op: OpId::new("wall1"),
            role: FaceRole::Tagged(SegmentTag::new("outer-mid")),
        };
        assert!(
            store_a.names().face(&parent).is_none(),
            "the severed parent name retires"
        );

        let split = |side: SplitSide| FaceName::Split {
            op: OpId::new("cut1"),
            parent: Box::new(parent.clone()),
            side,
        };
        let shell_a = store_a
            .shell(store_a.solid(result_a).unwrap().outer_shell)
            .unwrap();
        for side in [SplitSide::Left, SplitSide::Right] {
            let fa = store_a.names().face(&split(side)).expect("fragment in A");
            let fb = store_b.names().face(&split(side)).expect("fragment in B");
            assert!(shell_a.faces.contains(&fa));
            let sample = |store: &TopologyStore, f| match &store.face(f).unwrap().surface {
                FaceSurface::Nurbs(s) => s.point_at(0.5, 0.5).unwrap(),
                other => panic!("fragment must be NURBS, got {other:?}"),
            };
            assert!(
                (sample(&store_a, fa) - sample(&store_b, fb)).norm() < 1e-9,
                "{side:?} fragment moved across rebuilds"
            );
        }

        // Deterministic sides: Left = above the window head (positive
        // cross-product side of the canonical sill chord), Right = below
        // the sill. Verify via the fragments' z-extents.
        let z_extent = |store: &TopologyStore, f: FaceId| -> (f64, f64) {
            use crate::tessellation::{TessellateFace, TessellationParams};
            let mesh = TessellateFace::new(f, TessellationParams::default())
                .execute(store)
                .unwrap();
            let zmin = mesh
                .vertices
                .iter()
                .map(|p| p.z)
                .fold(f64::INFINITY, f64::min);
            let zmax = mesh
                .vertices
                .iter()
                .map(|p| p.z)
                .fold(f64::NEG_INFINITY, f64::max);
            (zmin, zmax)
        };
        let left = store_a.names().face(&split(SplitSide::Left)).unwrap();
        let right = store_a.names().face(&split(SplitSide::Right)).unwrap();
        let (left_zmin, _) = z_extent(&store_a, left);
        let (_, right_zmax) = z_extent(&store_a, right);
        assert!(
            left_zmin > 2.0 - 1e-9,
            "Left fragment lies above the window head (zmin = {left_zmin})"
        );
        assert!(
            right_zmax < 1.0 + 1e-9,
            "Right fragment lies below the sill (zmax = {right_zmax})"
        );

        // The neighbors keep their tagged names (single kept fragments).
        for tag in ["outer-a", "outer-c"] {
            let name = FaceName::Created {
                op: OpId::new("wall1"),
                role: FaceRole::Tagged(SegmentTag::new(tag)),
            };
            assert!(store_a.names().face(&name).is_some(), "{tag} resolves");
        }

        // Watertight.
        let boundary = welded_boundary_edges(&store_a, result_a);
        assert_eq!(
            boundary, 0,
            "middle-segment split must position-weld watertight \
             (found {boundary} boundary edges)"
        );

        // Parameter change: sliding the window (still covering outer-mid)
        // keeps the split names resolving.
        let (moved, _) = build(2.1);
        for side in [SplitSide::Left, SplitSide::Right] {
            assert!(
                moved.names().face(&split(side)).is_some(),
                "moved window: {side:?} fragment still resolves"
            );
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

        let _ = subtract_through_cut(&mut store, slab, tube, None).unwrap();

        for f in original_faces {
            assert!(
                store.face(f).unwrap().inner_wires.is_empty(),
                "input slab face must stay un-punched"
            );
        }
    }
}
