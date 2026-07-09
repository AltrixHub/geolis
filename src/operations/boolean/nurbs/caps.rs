//! Cap-notch rebuild for cap-touching cuts (F6 R2).
//!
//! A cap-touching (full-height door) cut splits the wall side faces' shared
//! ring edges and closes the doorway's circuit with cap-plane closure edges
//! ([`super::band::build_open_band_fragments`]). The affected planar cap is
//! rebuilt by WIRE SURGERY only: every split parent edge in its wire is
//! replaced by the sub-edges the kept wall fragments retained (the notched
//! span simply disappears), and the band's closure edges bridge the notch —
//! the SAME `EdgeId`s on both sides (F2 shared-edge convention), so the
//! result is watertight by construction. The kept edge pool (outer AND
//! inner wires — annulus caps carry one inner wire per footprint hole) is
//! then chained into connected cycles, classified by their winding in the
//! cap plane's 2D frame: cycles winding like the original outer wire become
//! planar cap fragments; opposite-winding cycles are hole loops (an
//! untouched courtyard wire riding along, or a hole ring notched within
//! itself) and are re-attached as inner wires of the fragment containing
//! them. A doorway through an annulus wall consumes BOTH wires' kept
//! sub-edges into one outer cycle — the courtyard hole merges into the
//! fragment boundary and the inner wire disappears.
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
use crate::geometry::curve::Curve;
use crate::math::{Point2, Point3};
use crate::topology::{
    EdgeCurve, EdgeId, FaceData, FaceId, FaceSurface, OpId, OrientedEdge, TopologyStore, WireData,
    WireId,
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
/// Typed errors when a pending cap is not planar, its notched boundary
/// does not chain into closed cycles, it yields more than two fragments,
/// a surviving hole loop lies in no kept fragment, or a closure edge
/// remains unconsumed.
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
    let to2d = |p: Point3| -> Point2 {
        let rel = p - *plane.origin();
        Point2::new(rel.dot(plane.u_dir()), rel.dot(plane.v_dir()))
    };

    // Chain the kept boundary pool + closure edges into closed cycles by
    // shared vertices.
    let pool = kept_boundary_pool(store, &face, sub_edges, kept_edges)?;
    let cycles = chain_cycles(store, &pool, closure_edges, used_closures)?;

    // Classify cycles by winding in the cap plane's 2D frame: cycles that
    // wind like the ORIGINAL outer wire are fragment outer boundaries;
    // opposite-winding cycles are hole loops that survived the notch (cap
    // inner wires wind opposite the outer by construction, and wire
    // surgery preserves traversal orientation).
    let outer_edges = store.wire(face.outer_wire)?.edges.clone();
    let outer_ccw = polygon_signed_area(&cycle_polygon(store, &outer_edges, &to2d)?) > 0.0;

    let mut fragment_cycles: Vec<Vec<OrientedEdge>> = Vec::new();
    let mut fragment_polygons: Vec<Vec<Point2>> = Vec::new();
    let mut hole_cycles: Vec<(Vec<OrientedEdge>, Vec<Point2>)> = Vec::new();
    for cycle in cycles {
        let polygon = cycle_polygon(store, &cycle, &to2d)?;
        if (polygon_signed_area(&polygon) > 0.0) == outer_ccw {
            fragment_cycles.push(cycle);
            fragment_polygons.push(polygon);
        } else {
            hole_cycles.push((cycle, polygon));
        }
    }

    let fragment_holes = assign_hole_cycles(store, hole_cycles, &fragment_polygons)?;

    // One planar fragment per outer cycle (same plane, same sense; the
    // planar tessellation is wire-driven, so no trim is carried).
    let mut fragment_faces = Vec::with_capacity(fragment_cycles.len());
    for (cycle, holes) in fragment_cycles.iter().zip(fragment_holes) {
        let new_wire = store.add_wire(WireData {
            edges: cycle.clone(),
            is_closed: true,
        });
        fragment_faces.push(store.add_face(FaceData {
            surface: face.surface.clone(),
            outer_wire: new_wire,
            inner_wires: holes,
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
                let (left, right) =
                    order_cap_fragments(store, *a, *b, &fragment_cycles, closure_edges, &to2d)?;
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

/// Builds the kept boundary pool: the cap's outer AND inner wires with
/// every split parent edge replaced by its KEPT sub-edges, in traversal
/// order and orientation.
fn kept_boundary_pool(
    store: &TopologyStore,
    face: &FaceData,
    sub_edges: &HashMap<EdgeId, Vec<EdgeId>>,
    kept_edges: &HashSet<EdgeId>,
) -> Result<Vec<OrientedEdge>> {
    let mut pool: Vec<OrientedEdge> = Vec::new();
    let wire_ids: Vec<WireId> = std::iter::once(face.outer_wire)
        .chain(face.inner_wires.iter().copied())
        .collect();
    for wire_id in wire_ids {
        let wire = store.wire(wire_id)?;
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
    }
    Ok(pool)
}

/// Re-attaches each surviving hole cycle to the fragment containing it,
/// returning one inner-wire list per fragment (parallel to
/// `fragment_polygons`).
fn assign_hole_cycles(
    store: &mut TopologyStore,
    hole_cycles: Vec<(Vec<OrientedEdge>, Vec<Point2>)>,
    fragment_polygons: &[Vec<Point2>],
) -> Result<Vec<Vec<WireId>>> {
    let mut fragment_holes: Vec<Vec<WireId>> = vec![Vec::new(); fragment_polygons.len()];
    for (cycle, polygon) in hole_cycles {
        let sample = polygon.first().copied().ok_or_else(|| {
            OperationError::Failed("notched cap hole loop has no boundary samples".into())
        })?;
        let containing = fragment_polygons
            .iter()
            .position(|frag| super::split::polygon_contains(frag, sample))
            .ok_or_else(|| {
                OperationError::Failed(
                    "notched cap hole loop lies in no kept fragment \
                     (inconsistent cap-touching cut)"
                        .into(),
                )
            })?;
        fragment_holes[containing].push(store.add_wire(WireData {
            edges: cycle,
            is_closed: true,
        }));
    }
    Ok(fragment_holes)
}

/// Samples a cycle's oriented edges into a polygon in the cap plane's 2D
/// frame (interior samples per edge in traversal order; each edge's tail
/// sample is dropped — it coincides with the next edge's head).
fn cycle_polygon(
    store: &TopologyStore,
    cycle: &[OrientedEdge],
    to2d: &impl Fn(Point3) -> Point2,
) -> Result<Vec<Point2>> {
    /// Interior samples per edge (winding / containment only — never used
    /// for boundary geometry).
    const SAMPLES: usize = 8;

    let mut poly = Vec::with_capacity(cycle.len() * SAMPLES);
    for oe in cycle {
        let edge = store.edge(oe.edge)?;
        let (t0, t1) = if oe.forward {
            (edge.t_start, edge.t_end)
        } else {
            (edge.t_end, edge.t_start)
        };
        for i in 0..SAMPLES {
            #[allow(clippy::cast_precision_loss)]
            let frac = i as f64 / SAMPLES as f64;
            let t = t0 + (t1 - t0) * frac;
            let p = match &edge.curve {
                EdgeCurve::Line(c) => c.evaluate(t)?,
                EdgeCurve::Arc(c) => c.evaluate(t)?,
                EdgeCurve::Circle(c) => c.evaluate(t)?,
                EdgeCurve::Ellipse(c) => c.evaluate(t)?,
                EdgeCurve::Nurbs(c) => c.point_at(t)?,
            };
            poly.push(to2d(p));
        }
    }
    Ok(poly)
}

/// Shoelace signed area of a 2D polygon (positive = counter-clockwise).
fn polygon_signed_area(poly: &[Point2]) -> f64 {
    let n = poly.len();
    let mut area2 = 0.0;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        area2 += a.x * b.y - b.x * a.y;
    }
    0.5 * area2
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
/// lexicographically first; a fragment is Left when its representative
/// interior point lies on the positive cross-product side of that chord.
///
/// The representative point is a LOCAL probe taken just inside each
/// fragment at the midpoint of its own cut-trace (closure) edge — NOT the
/// fragment's global centroid. A cap severed by a second doorway can leave
/// long fragments that wrap around the annulus (e.g. one arm runs back
/// along the far wall), and their vertex-average centroid can land on the
/// wrong side of — or exactly on — the extended chord line, which the
/// centroid test misreads as an ambiguous split. The doorway separates the
/// two fragments, so each one's boundary carries a cut-trace edge and the
/// interior immediately inside that edge is unambiguously on that
/// fragment's side of the cut.
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

    let side = |center: Point2| -> f64 {
        let dir = chord.1 - chord.0;
        let rel = center - chord.0;
        dir.x * rel.y - dir.y * rel.x
    };
    let cycle_of = |face: FaceId| if face == a { &cycles[0] } else { &cycles[1] };
    let side_a = side(fragment_interior_probe(
        store,
        cycle_of(a),
        closure_edges,
        to2d,
    )?);
    let side_b = side(fragment_interior_probe(
        store,
        cycle_of(b),
        closure_edges,
        to2d,
    )?);
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

/// A representative interior point of a cap fragment, taken just inside its
/// boundary at the midpoint of its own cut-trace (closure) edge.
///
/// Unlike the fragment's global centroid this stays LOCAL to the doorway,
/// so it lands on the fragment's true side of the cut even when the
/// fragment wraps far around an annulus cap. The step direction is the
/// inward normal of the trace edge, chosen from the fragment's winding
/// (interior is left of a CCW boundary, right of a CW one), and the step is
/// a small fraction of the edge length so the probe never leaves the local
/// neighbourhood of the trace.
fn fragment_interior_probe(
    store: &TopologyStore,
    cycle: &[OrientedEdge],
    closure_edges: &[EdgeId],
    to2d: &impl Fn(Point3) -> Point2,
) -> Result<Point2> {
    // A small, geometry-scaled step into the interior (fraction of the
    // trace edge length via the inward normal below).
    const STEP: f64 = 1e-3;

    let trace = cycle
        .iter()
        .find(|oe| closure_edges.contains(&oe.edge))
        .ok_or_else(|| {
            OperationError::Failed(
                "cap fragment carries no cut-trace edge (cannot classify split side)".into(),
            )
        })?;
    let edge = store.edge(trace.edge)?;
    let (start, end) = if trace.forward {
        (edge.start, edge.end)
    } else {
        (edge.end, edge.start)
    };
    let p0 = to2d(store.vertex(start)?.point);
    let p1 = to2d(store.vertex(end)?.point);
    let mid = Point2::new((p0.x + p1.x) * 0.5, (p0.y + p1.y) * 0.5);
    let dir = p1 - p0;
    // Interior is left of the directed boundary for a CCW cycle, right for
    // a CW one; the inward normal has magnitude = edge length, so a fixed
    // small fraction of it is a geometry-scaled step into the interior.
    let ccw = polygon_signed_area(&cycle_polygon(store, cycle, to2d)?) > 0.0;
    let inward = if ccw {
        Point2::new(-dir.y, dir.x)
    } else {
        Point2::new(dir.y, -dir.x)
    };
    Ok(Point2::new(
        mid.x + STEP * inward.x,
        mid.y + STEP * inward.y,
    ))
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
