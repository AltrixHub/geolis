//! Cap-notch rebuild for cap-touching cuts (F6 R2).
//!
//! A cap-touching (full-height door) cut splits the wall side faces' shared
//! ring edges and closes the doorway's circuit with cap-plane closure edges
//! ([`super::band::build_open_band_fragments`]). The affected planar cap is
//! rebuilt by WIRE SURGERY only: every split parent edge in its wire is
//! replaced by the sub-edges the kept wall fragments retained (the notched
//! span simply disappears), and the band's closure edges bridge the notch —
//! the SAME `EdgeId`s on both sides (F2 shared-edge convention), so the
//! result is watertight by construction. The kept edge pool is then chained
//! into connected cycles; each cycle becomes one planar cap fragment.
//!
//! Names follow the F5 split rule: one kept fragment transfers the cap's
//! name; two bind [`FaceName::Split`] `Left` / `Right` by the canonical-
//! chord rule projected into the cap plane's own 2D frame (the closure
//! edges are the cap's cut traces); three or more kept fragments are a
//! typed error. Unnamed caps propagate unnamed fragments.
//!
//! [`FaceName::Split`]: crate::topology::FaceName::Split

use std::collections::{HashMap, HashSet};

use crate::error::{OperationError, Result};
use crate::math::{Point2, Point3};
use crate::topology::{
    EdgeId, FaceData, FaceId, FaceSurface, OpId, OrientedEdge, TopologyStore, WireData,
};

/// Rebuilds every pending notched cap and returns the new cap fragment
/// faces (appended to the result shell by the caller).
///
/// `sub_edges` maps each split parent edge to its sub-edges (parent
/// parameter order); `kept_edges` is the set of edges retained by the kept
/// wall fragments' wires (a sub-edge absent from it was removed doorway
/// material); `closure_edges` are the cap-plane closure edges of every
/// open-chain cut in this boolean.
///
/// # Errors
///
/// Typed errors when a pending cap is not planar, carries inner wires,
/// its notched boundary does not chain into closed cycles, it yields more
/// than two fragments, or a closure edge remains unconsumed.
pub(crate) fn rebuild_notched_caps(
    store: &mut TopologyStore,
    pending: &[FaceId],
    sub_edges: &HashMap<EdgeId, Vec<EdgeId>>,
    kept_edges: &HashSet<EdgeId>,
    closure_edges: &[EdgeId],
    op_id: Option<&OpId>,
) -> Result<Vec<FaceId>> {
    let mut used_closures: HashSet<EdgeId> = HashSet::new();
    let mut out = Vec::new();
    for &cap in pending {
        let fragments = rebuild_cap(
            store,
            cap,
            sub_edges,
            kept_edges,
            closure_edges,
            &mut used_closures,
            op_id,
        )?;
        out.extend(fragments);
    }
    if used_closures.len() != closure_edges.len() {
        return Err(OperationError::Failed(
            "a cap-plane closure edge was not consumed by any notched cap \
             (inconsistent cap-touching cut)"
                .into(),
        )
        .into());
    }
    Ok(out)
}

/// Rebuilds one notched cap into its kept fragments.
#[allow(clippy::too_many_arguments)]
fn rebuild_cap(
    store: &mut TopologyStore,
    cap: FaceId,
    sub_edges: &HashMap<EdgeId, Vec<EdgeId>>,
    kept_edges: &HashSet<EdgeId>,
    closure_edges: &[EdgeId],
    used_closures: &mut HashSet<EdgeId>,
    op_id: Option<&OpId>,
) -> Result<Vec<FaceId>> {
    let face = store.face(cap)?.clone();
    let FaceSurface::Plane(plane) = &face.surface else {
        return Err(
            OperationError::Failed("cap-notch rebuild requires a planar cap face".into()).into(),
        );
    };
    if !face.inner_wires.is_empty() {
        return Err(OperationError::Failed(
            "cap-notch rebuild does not support caps with inner wires".into(),
        )
        .into());
    }

    // Kept boundary pool: the cap's wire with every split parent edge
    // replaced by its KEPT sub-edges, in traversal order and orientation.
    let wire = store.wire(face.outer_wire)?.clone();
    let mut pool: Vec<OrientedEdge> = Vec::new();
    for oe in &wire.edges {
        match sub_edges.get(&oe.edge) {
            None => pool.push(*oe),
            Some(subs) => {
                let ordered: Vec<EdgeId> = if oe.forward {
                    subs.clone()
                } else {
                    subs.iter().rev().copied().collect()
                };
                for sub in ordered {
                    if kept_edges.contains(&sub) {
                        pool.push(OrientedEdge::new(sub, oe.forward));
                    }
                }
            }
        }
    }

    // Chain the pool + closure edges into closed cycles by shared vertices.
    let cycles = chain_cycles(store, &pool, closure_edges, used_closures)?;

    // One planar fragment per cycle (same plane, same sense; the planar
    // tessellation is wire-driven, so no trim is carried).
    let mut fragment_faces = Vec::with_capacity(cycles.len());
    for cycle in &cycles {
        let new_wire = store.add_wire(WireData {
            edges: cycle.clone(),
            is_closed: true,
        });
        fragment_faces.push(store.add_face(FaceData {
            surface: face.surface.clone(),
            outer_wire: new_wire,
            inner_wires: Vec::new(),
            same_sense: face.same_sense,
            trim: None,
            pcurves: Vec::new(),
        }));
    }

    // Name evolution: the F5 split rule.
    match fragment_faces.as_slice() {
        [] => Err(
            OperationError::Failed("cap-notch rebuild produced no kept fragments".into()).into(),
        ),
        [single] => {
            store.names_mut().transfer_face(cap, *single);
            Ok(fragment_faces)
        }
        [a, b] => {
            if let Some(op) = op_id {
                let to2d = |p: Point3| -> Point2 {
                    let rel = p - *plane.origin();
                    Point2::new(rel.dot(plane.u_dir()), rel.dot(plane.v_dir()))
                };
                let (left, right) =
                    order_cap_fragments(store, *a, *b, &cycles, closure_edges, &to2d)?;
                store.names_mut().split_face(cap, op, left, right);
            }
            Ok(fragment_faces)
        }
        _ => Err(OperationError::Failed(
            "cap-notch rebuild produced more than two kept fragments \
             (unsupported)"
                .into(),
        )
        .into()),
    }
}

/// Chains the kept boundary pool into closed cycles, pulling in closure
/// edges (oriented as needed) to bridge the notches. Deterministic: cycles
/// start at the earliest unused pool edge in wire order.
fn chain_cycles(
    store: &TopologyStore,
    pool: &[OrientedEdge],
    closure_edges: &[EdgeId],
    used_closures: &mut HashSet<EdgeId>,
) -> Result<Vec<Vec<OrientedEdge>>> {
    let oriented_ends =
        |oe: &OrientedEdge| -> Result<(crate::topology::VertexId, crate::topology::VertexId)> {
            let edge = store.edge(oe.edge)?;
            Ok(if oe.forward {
                (edge.start, edge.end)
            } else {
                (edge.end, edge.start)
            })
        };

    let mut used_pool = vec![false; pool.len()];
    let mut cycles: Vec<Vec<OrientedEdge>> = Vec::new();
    for first in 0..pool.len() {
        if used_pool[first] {
            continue;
        }
        used_pool[first] = true;
        let mut cycle = vec![pool[first]];
        let (start_vertex, mut current) = oriented_ends(&pool[first])?;
        while current != start_vertex {
            // Prefer continuing along the kept boundary; bridge a notch
            // with a closure edge otherwise.
            let mut next: Option<OrientedEdge> = None;
            for (idx, oe) in pool.iter().enumerate() {
                if used_pool[idx] {
                    continue;
                }
                if oriented_ends(oe)?.0 == current {
                    if next.is_some() {
                        return Err(OperationError::Failed(
                            "notched cap boundary is ambiguous (two kept \
                             edges start at one vertex)"
                                .into(),
                        )
                        .into());
                    }
                    next = Some(*oe);
                    used_pool[idx] = true;
                }
            }
            if next.is_none() {
                for &closure in closure_edges {
                    if used_closures.contains(&closure) {
                        continue;
                    }
                    let edge = store.edge(closure)?;
                    let oriented = if edge.start == current {
                        Some(OrientedEdge::new(closure, true))
                    } else if edge.end == current {
                        Some(OrientedEdge::new(closure, false))
                    } else {
                        None
                    };
                    if let Some(oriented) = oriented {
                        if next.is_some() {
                            return Err(OperationError::Failed(
                                "notched cap boundary is ambiguous (two \
                                 closure edges meet one vertex)"
                                    .into(),
                            )
                            .into());
                        }
                        used_closures.insert(closure);
                        next = Some(oriented);
                    }
                }
            }
            let Some(next) = next else {
                return Err(OperationError::Failed(
                    "notched cap boundary does not close (missing kept \
                     sub-edge or closure edge)"
                        .into(),
                )
                .into());
            };
            current = oriented_ends(&next)?.1;
            cycle.push(next);
        }
        cycles.push(cycle);
    }
    Ok(cycles)
}

/// Orders two cap fragments into `(left, right)` by the canonical-chord
/// rule in the cap plane's 2D frame: the canonical trace is the consumed
/// closure edge whose canonically-oriented projected chord starts
/// lexicographically first; a fragment is Left when its vertex-average
/// centroid lies on the positive cross-product side of that chord.
fn order_cap_fragments(
    store: &TopologyStore,
    a: FaceId,
    b: FaceId,
    cycles: &[Vec<OrientedEdge>],
    closure_edges: &[EdgeId],
    to2d: &impl Fn(Point3) -> Point2,
) -> Result<(FaceId, FaceId)> {
    // The consumed closure edges of THIS cap: those appearing in a cycle.
    let mut chords: Vec<(Point2, Point2)> = Vec::new();
    for cycle in cycles {
        for oe in cycle {
            if closure_edges.contains(&oe.edge) {
                let edge = store.edge(oe.edge)?;
                let p0 = to2d(store.vertex(edge.start)?.point);
                let p1 = to2d(store.vertex(edge.end)?.point);
                chords.push(canonical_chord(p0, p1));
            }
        }
    }
    let chord = chords
        .into_iter()
        .min_by(|x, y| {
            let kx = (x.0.x, x.0.y, x.1.x, x.1.y);
            let ky = (y.0.x, y.0.y, y.1.x, y.1.y);
            kx.partial_cmp(&ky).unwrap_or(std::cmp::Ordering::Equal)
        })
        .ok_or_else(|| OperationError::Failed("split cap without any closure edge".into()))?;

    let centroid = |face: FaceId| -> Result<Point2> {
        let cycle = if face == a { &cycles[0] } else { &cycles[1] };
        let mut sum = Point2::new(0.0, 0.0);
        let mut count = 0usize;
        for oe in cycle {
            let edge = store.edge(oe.edge)?;
            let p = to2d(store.vertex(edge.start)?.point);
            sum = Point2::new(sum.x + p.x, sum.y + p.y);
            count += 1;
        }
        #[allow(clippy::cast_precision_loss)]
        let inv = 1.0 / count.max(1) as f64;
        Ok(Point2::new(sum.x * inv, sum.y * inv))
    };

    let side = |center: Point2| -> f64 {
        let dir = chord.1 - chord.0;
        let rel = center - chord.0;
        dir.x * rel.y - dir.y * rel.x
    };
    let side_a = side(centroid(a)?);
    let side_b = side(centroid(b)?);
    if side_a > 0.0 && side_b < 0.0 {
        Ok((a, b))
    } else if side_a < 0.0 && side_b > 0.0 {
        Ok((b, a))
    } else {
        Err(OperationError::Failed(
            "cap fragments do not lie on opposite sides of the canonical \
             closure chord (ambiguous SplitSide)"
                .into(),
        )
        .into())
    }
}

/// Orients a chord so `end - start` is lexicographically positive.
fn canonical_chord(a: Point2, b: Point2) -> (Point2, Point2) {
    let d = b - a;
    if d.x > 0.0 || (d.x == 0.0 && d.y > 0.0) {
        (a, b)
    } else {
        (b, a)
    }
}
