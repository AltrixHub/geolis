//! Constrained-Delaunay tessellation of trimmed NURBS faces.
//!
//! The trim loops are sampled into UV polylines and inserted into a spade
//! constrained Delaunay triangulation as constraint edges. Grid points from the
//! shared adaptive refinement (see [`super::tessellate_nurbs`]) that fall
//! strictly inside the trim region are added as Steiner points. Triangles whose
//! centroid lies outside the trim region (outside the outer loop or inside a
//! hole) are discarded, and the surviving vertices are mapped through the
//! surface to 3D.
//!
//! ## Boundary conformance (two shipped layers)
//!
//! **Design (b) — shared-edge `BRep` adjacency (shipped for prism / tube /
//! revolved solids):** adjacent faces reference the same [`EdgeId`] for their
//! common boundary curve; the edge is sampled once per
//! [`super::edge_samples::EdgeSampleCache`] and every face consumes the same
//! polyline. NURBS faces carrying per-edge pcurves take the edge-driven outer
//! loop ([`edge_driven_outer_uv`]) and edge-driven hole loops
//! ([`face_hole_loops_uv`]); planar caps consume the cached 3D samples
//! directly. On the edge-driven paths every boundary vertex is **3D-pinned**
//! ([`UvPinMap`]): the emitted position is the canonical cached edge sample,
//! not `surface.point_at(uv)`, so two faces meeting at an edge emit
//! bit-identical boundary vertices even when their surfaces disagree by the
//! SSI marcher's acceptance residual (the F6 R3 rim weld). Closed side
//! surfaces (whose seam is not a topological edge) close their UV rectangle
//! with seam connectors.
//!
//! **Design (a) — solid-agnostic boundary-conforming sampling (fallback for
//! faces without pcurves, e.g. the curved slab / wall builders):** any face
//! whose outer boundary is the full parameter rectangle — an untrimmed face
//! (via [`tessellate_untrimmed_conforming`]) or a punched target face whose
//! outer loop is the domain rectangle (detected by
//! [`loop_is_domain_rectangle`]) — samples that boundary at the
//! *curve-intrinsic* parameters of each boundary isocurve
//! ([`super::tessellate_nurbs::conforming_boundary_uv`]). Because that
//! sampling is a function of the boundary curve's geometry alone, two faces
//! sharing a curve independently arrive at the identical parameter set; the
//! deviation collapses to floating-point noise. The per-edge cache uses the
//! same chord-adaptive algorithm, so the two layers agree exactly where they
//! meet.

use std::collections::HashMap;

use spade::{
    ConstrainedDelaunayTriangulation, InsertionError, Point2 as SpadePoint2, Triangulation,
};

use crate::error::{Result, TessellationError};
use crate::geometry::nurbs::{NurbsCurve2D, NurbsSurface};
use crate::geometry::surface::Surface;
use crate::math::{Point2, Point3, Vector3};
use crate::topology::{FaceData, FaceTrim, TopologyStore, TrimLoop, WireId};

use super::edge_samples::EdgeSampleCache;
use super::tessellate_nurbs::{
    adaptive_grid_parameters, conforming_boundary_uv, BOUNDARY_CHORD_TOLERANCE,
};
use super::{SurfaceTessellationOptions, TriangleMesh};

/// Number of samples per pcurve segment when converting a trim loop to a UV
/// polyline. Trim loops in P5 are low-degree (lines / rational arcs), so a fixed
/// per-segment sampling is sufficient and avoids the spade duplicate-vertex
/// hazards that adaptive sampling near-coincident points would create.
const PCURVE_SEGMENT_SAMPLES: usize = 32;

/// Minimum UV separation between consecutive polyline points. Points closer than
/// this are merged so spade never receives duplicate constraint vertices.
const MERGE_EPS: f64 = 1e-9;

/// Canonical 3D positions for boundary UV vertices lying on shared edges,
/// keyed by the exact UV bit pattern.
///
/// Edge-driven boundary loops record one entry per (pcurve-mapped) edge
/// sample; when the CDT emits a vertex at a pinned UV, its 3D position is the
/// cached edge sample instead of `surface.point_at(uv)`. Every face
/// referencing the edge therefore emits the bit-identical 3D vertex — the
/// same canonical polyline the planar-cap path consumes directly — so
/// cross-face rim coincidence is exact (no marcher-residual near-duplicates,
/// no tolerance-based post-weld).
pub(crate) type UvPinMap = HashMap<(u64, u64), Point3>;

/// The exact-bits pin key of a UV point.
fn pin_key(p: &Point2) -> (u64, u64) {
    (p.x.to_bits(), p.y.to_bits())
}

/// Tessellates a trimmed NURBS face into a triangle mesh.
///
/// # Errors
///
/// Returns [`TessellationError::InvalidParameters`] if `options` are invalid,
/// [`TessellationError::Failed`] if a trim loop degenerates (fewer than 3
/// distinct UV points) or the constrained triangulation cannot be built, and
/// propagates any surface evaluation error.
pub fn tessellate_trimmed_nurbs_face(
    surface: &NurbsSurface,
    trim: &FaceTrim,
    options: &SurfaceTessellationOptions,
) -> Result<TriangleMesh> {
    validate_options(options)?;

    // 1. Sample each loop into a UV polyline. When the outer loop is the full
    //    parameter rectangle (e.g. a punched target face whose silhouette is the
    //    original NURBS boundary), sample it at the boundary-curve-intrinsic
    //    parameters so it conforms with adjacent faces that share those curves.
    let outer = outer_loop_uv(surface, &trim.outer)?;
    let holes: Vec<Vec<Point2>> = trim.holes.iter().map(sample_loop).collect::<Result<_>>()?;

    tessellate_cdt(surface, &outer, &holes, &UvPinMap::new(), options)
}

/// Tessellates an untrimmed NURBS face whose outer boundary is the full
/// parameter rectangle, using the boundary-conforming CDT path.
///
/// Adjacent faces sharing a boundary curve sample it identically (see
/// [`conforming_boundary_uv`]), so their meshes meet with no silhouette sliver —
/// unlike the tensor-grid tessellator, whose per-face boundary chords disagree.
///
/// # Errors
///
/// Propagates option-validation, boundary-sampling, and CDT construction errors.
pub(crate) fn tessellate_untrimmed_conforming(
    surface: &NurbsSurface,
    options: &SurfaceTessellationOptions,
) -> Result<TriangleMesh> {
    validate_options(options)?;
    let outer = conforming_boundary_uv(surface, BOUNDARY_CHORD_TOLERANCE)?;
    tessellate_cdt(surface, &outer, &[], &UvPinMap::new(), options)
}

/// Chooses the outer-loop UV sampling: boundary-conforming when the loop is the
/// full parameter rectangle, otherwise the loop's own (e.g. SSI ring) sampling.
fn outer_loop_uv(surface: &NurbsSurface, outer: &TrimLoop) -> Result<Vec<Point2>> {
    if loop_is_domain_rectangle(surface, outer) {
        conforming_boundary_uv(surface, BOUNDARY_CHORD_TOLERANCE)
    } else {
        sample_loop(outer)
    }
}

/// Tessellates a NURBS face whose outer UV loop and hole loops were assembled
/// by the caller (edge-driven path), with the caller's pin map anchoring
/// boundary vertices to the canonical shared-edge samples.
///
/// # Errors
///
/// Propagates option-validation and CDT construction errors.
pub(crate) fn tessellate_with_outer_uv(
    surface: &NurbsSurface,
    outer: &[Point2],
    holes: &[Vec<Point2>],
    pins: &UvPinMap,
    options: &SurfaceTessellationOptions,
) -> Result<TriangleMesh> {
    validate_options(options)?;
    tessellate_cdt(surface, outer, holes, pins, options)
}

/// One UV chain per wire edge (traversal order), mapped from the shared edge
/// samples through the face's pcurves, with every sample's canonical 3D point
/// recorded in `pins`. Returns `None` when the face records no pcurve for
/// some wire edge (legacy face) — the caller falls back to geometric paths.
fn wire_uv_chains(
    store: &TopologyStore,
    cache: &mut EdgeSampleCache,
    face: &FaceData,
    wire_id: WireId,
    pins: &mut UvPinMap,
) -> Result<Option<Vec<Vec<Point2>>>> {
    let Ok(wire) = store.wire(wire_id) else {
        return Ok(None);
    };
    if wire.edges.is_empty() {
        return Ok(None);
    }
    let oriented: Vec<crate::topology::OrientedEdge> = wire.edges.clone();

    let mut chains: Vec<Vec<Point2>> = Vec::with_capacity(oriented.len());
    for oe in &oriented {
        let Some(pcurve) = face.pcurve_for(oe.edge) else {
            return Ok(None);
        };
        let samples = cache.get(store, oe.edge)?;
        let mut uv = Vec::with_capacity(samples.params.len());
        for (&t, &p3) in samples.params.iter().zip(&samples.points) {
            let q = pcurve.point_at(t)?;
            pins.insert(pin_key(&q), p3);
            uv.push(q);
        }
        if !oe.forward {
            uv.reverse();
        }
        chains.push(uv);
    }
    Ok(Some(chains))
}

/// Builds the outer UV polyline of a face from its outer wire's shared edge
/// samples mapped through the face's pcurves — the edge-driven boundary path
/// of shared-edge topology. Every face referencing an edge consumes the same
/// cached samples (and pins their canonical 3D positions), so adjacent faces
/// emit identical 3D boundary vertices by construction.
///
/// Consecutive wire edges whose UV images do not meet (a geometrically closed
/// direction whose seam is not a topological edge) are joined by an
/// axis-aligned seam connector subdivided at the adaptive grid parameters, so
/// the CDT boundary follows the surface curvature along the seam.
///
/// Returns `None` when the face records no pcurve for some outer wire edge
/// (legacy face) — the caller falls back to the geometric paths.
///
/// # Errors
///
/// Propagates store lookups, pcurve evaluation, and grid-parameter errors.
pub(crate) fn edge_driven_outer_uv(
    store: &TopologyStore,
    cache: &mut EdgeSampleCache,
    face: &FaceData,
    surface: &NurbsSurface,
    options: &SurfaceTessellationOptions,
    pins: &mut UvPinMap,
) -> Result<Option<Vec<Point2>>> {
    let Some(chains) = wire_uv_chains(store, cache, face, face.outer_wire, pins)? else {
        return Ok(None);
    };

    // Assemble the loop, inserting seam connectors where consecutive chains do
    // not meet in UV.
    let (u_grid, v_grid) = adaptive_grid_parameters(surface, options);
    let mut pts: Vec<Point2> = Vec::new();
    let n = chains.len();
    for i in 0..n {
        pts.extend_from_slice(&chains[i]);
        let Some(&end) = chains[i].last() else {
            return Ok(None);
        };
        let start = chains[(i + 1) % n][0];
        if (end - start).norm() > MERGE_EPS {
            append_seam_connector(&mut pts, end, start, &u_grid, &v_grid);
        }
    }
    dedup_closed(&mut pts);
    if pts.len() < 3 {
        return Ok(None);
    }
    Ok(Some(pts))
}

/// Builds one hole's UV polyline from its inner wire's shared edge samples
/// mapped through the face's pcurves (with canonical 3D pins) — the
/// edge-driven counterpart of [`edge_driven_outer_uv`] for hole rings.
///
/// Hole rings are interior loops, so consecutive chains must meet exactly
/// (no seam connectors); a gap or a missing pcurve falls back to the trim
/// hole sampling (`None`).
fn edge_driven_hole_uv(
    store: &TopologyStore,
    cache: &mut EdgeSampleCache,
    face: &FaceData,
    wire_id: WireId,
    pins: &mut UvPinMap,
) -> Result<Option<Vec<Point2>>> {
    let Some(chains) = wire_uv_chains(store, cache, face, wire_id, pins)? else {
        return Ok(None);
    };
    let mut pts: Vec<Point2> = Vec::new();
    let n = chains.len();
    for i in 0..n {
        let Some(&end) = chains[i].last() else {
            return Ok(None);
        };
        if (end - chains[(i + 1) % n][0]).norm() > MERGE_EPS {
            return Ok(None);
        }
        pts.extend_from_slice(&chains[i]);
    }
    dedup_closed(&mut pts);
    if pts.len() < 3 {
        return Ok(None);
    }
    Ok(Some(pts))
}

/// Assembles a face's hole UV loops for the edge-driven tessellation path:
/// inner wires whose edges all carry pcurves sample edge-driven (pinned to
/// the canonical shared-edge 3D polylines); the rest fall back to the
/// matching trim hole loop (punch lockstep order). An inner wire with
/// neither pcurves nor a trim hole is 3D-only bookkeeping (e.g. a band
/// face's exit ring) and contributes no CDT hole.
///
/// # Errors
///
/// Propagates store lookups, pcurve evaluation, and trim sampling errors.
pub(crate) fn face_hole_loops_uv(
    store: &TopologyStore,
    cache: &mut EdgeSampleCache,
    face: &FaceData,
    pins: &mut UvPinMap,
) -> Result<Vec<Vec<Point2>>> {
    let trim_holes: &[TrimLoop] = face.trim.as_ref().map_or(&[], |t| t.holes.as_slice());
    let mut holes: Vec<Vec<Point2>> = Vec::with_capacity(face.inner_wires.len());
    for (i, &wire_id) in face.inner_wires.iter().enumerate() {
        if let Some(loop_uv) = edge_driven_hole_uv(store, cache, face, wire_id, pins)? {
            holes.push(loop_uv);
        } else if let Some(hole) = trim_holes.get(i) {
            holes.push(sample_loop(hole)?);
        }
    }
    // Trim holes beyond the paired range (no matching inner wire) keep the
    // plain trim sampling.
    for hole in trim_holes.iter().skip(face.inner_wires.len()) {
        holes.push(sample_loop(hole)?);
    }
    Ok(holes)
}

/// Appends the interior points of a straight axis-aligned UV connector from
/// `from` to `to` (exclusive on both ends), subdivided at the adaptive grid
/// parameters of the varying direction. A non-axis-aligned connector (not
/// expected from rectangular-domain faces) gets no interior points — a plain
/// chord.
fn append_seam_connector(
    pts: &mut Vec<Point2>,
    from: Point2,
    to: Point2,
    u_grid: &[f64],
    v_grid: &[f64],
) {
    let vertical = (from.x - to.x).abs() < MERGE_EPS;
    let horizontal = (from.y - to.y).abs() < MERGE_EPS;
    if vertical {
        let (lo, hi) = (from.y.min(to.y), from.y.max(to.y));
        let mut inner: Vec<f64> = v_grid
            .iter()
            .copied()
            .filter(|&v| v > lo + MERGE_EPS && v < hi - MERGE_EPS)
            .collect();
        if from.y > to.y {
            inner.reverse();
        }
        pts.extend(inner.into_iter().map(|v| Point2::new(from.x, v)));
    } else if horizontal {
        let (lo, hi) = (from.x.min(to.x), from.x.max(to.x));
        let mut inner: Vec<f64> = u_grid
            .iter()
            .copied()
            .filter(|&u| u > lo + MERGE_EPS && u < hi - MERGE_EPS)
            .collect();
        if from.x > to.x {
            inner.reverse();
        }
        pts.extend(inner.into_iter().map(|u| Point2::new(u, from.y)));
    }
}

/// Reports whether every control point of `loop_` lies on the surface's
/// parameter-domain boundary — i.e. the loop is the full-domain rectangle (as
/// built by the punch pipeline), not an interior SSI trim ring.
fn loop_is_domain_rectangle(surface: &NurbsSurface, loop_: &TrimLoop) -> bool {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_domain();
    let eps = 1e-7;
    let on_boundary = |p: &Point2| {
        ((p.x - u_min).abs() < eps || (p.x - u_max).abs() < eps)
            || ((p.y - v_min).abs() < eps || (p.y - v_max).abs() < eps)
    };
    loop_
        .curves
        .iter()
        .all(|c| c.control_points().iter().all(on_boundary))
}

/// Constrained-Delaunay tessellation of a NURBS surface from pre-sampled UV
/// outer/hole loops. Interior grid samples (from the shared adaptive refinement)
/// that fall strictly inside the trim region are added as Steiner points.
///
/// Triangles are kept by TOPOLOGICAL parity classification (flood fill from
/// spade's outer face, flipping at constraint edges — the same classifier the
/// planar CDT path uses), not by a geometric centroid-in-polygon test. The
/// geometric test is a coin flip for the degenerate sliver triangles that
/// arise when a trim polyline is straight up to floating-point noise (e.g. a
/// planar tool's SSI rim, constant `v` up to ~1e-16): the zigzag spawns
/// zero-area triangles just outside the constraint chain whose centroids sit
/// ON the boundary. Parity classifies them robustly as exterior.
fn tessellate_cdt(
    surface: &NurbsSurface,
    outer: &[Point2],
    holes: &[Vec<Point2>],
    pins: &UvPinMap,
    options: &SurfaceTessellationOptions,
) -> Result<TriangleMesh> {
    // Adaptive UV grid (shared with the full-domain tessellator).
    let (u_params, v_params) = adaptive_grid_parameters(surface, options);

    // Constrained Delaunay triangulation.
    let mut cdt = ConstrainedDelaunayTriangulation::<SpadePoint2<f64>>::new();

    insert_constraint_loop(&mut cdt, outer)?;
    for hole in holes {
        insert_constraint_loop(&mut cdt, hole)?;
    }

    // Steiner points: grid samples strictly inside the trim region. Points on
    // (or within noise of) a constraint segment are skipped: inserting a point
    // that lies exactly on a constraint makes spade SPLIT the constraint,
    // adding a boundary vertex the adjacent face does not have — a conformance
    // crack of one chord sagitta. Domain-boundary grid rows/columns land
    // exactly on a full-domain outer loop, so this guard is load-bearing.
    for &u in &u_params {
        for &v in &v_params {
            let p = Point2::new(u, v);
            if point_in_region(&p, outer, holes) && !near_any_segment(&p, outer, holes) {
                // Ignore individual insertion failures (e.g. a point landing on
                // an existing constraint vertex); the constraint loops already
                // pin the boundary.
                let _ = cdt.insert(SpadePoint2::new(u, v));
            }
        }
    }

    // 4. Keep interior triangles (parity flood fill); map to 3D.
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_domain();
    let u_span = u_max - u_min;
    let v_span = v_max - v_min;

    let mut mesh = TriangleMesh::default();
    let mut vertex_map: HashMap<usize, u32> = HashMap::new();

    let interior_faces = super::tessellate_face::classify_interior_faces(&cdt);
    for face in cdt.inner_faces() {
        if !interior_faces.contains(&face.fix().index()) {
            continue;
        }
        let verts = face.vertices();

        let mut tri = [0u32; 3];
        for (slot, vh) in verts.iter().enumerate() {
            let idx = vh.fix().index();
            let mesh_idx = if let Some(&existing) = vertex_map.get(&idx) {
                existing
            } else {
                let pos = vh.position();
                let (u, v) = (pos.x, pos.y);
                // Boundary vertices on shared edges take the canonical
                // per-edge 3D sample (bit-identical across the faces that
                // reference the edge); everything else evaluates the surface.
                let p3 = match pins.get(&pin_key(&Point2::new(u, v))) {
                    Some(&pinned) => pinned,
                    None => surface.point_at(u, v)?,
                };
                // A collapsed pole yields a zero normal; fall back to +Z so the
                // mesh stays well-formed, matching the full-domain tessellator.
                let n = surface.normal(u, v).unwrap_or_else(|_| Vector3::z());
                #[allow(clippy::cast_possible_truncation)]
                let new_idx = mesh.vertices.len() as u32;
                mesh.vertices.push(p3);
                mesh.normals.push(n);
                let su = if u_span.abs() > f64::EPSILON {
                    (u - u_min) / u_span
                } else {
                    0.0
                };
                let sv = if v_span.abs() > f64::EPSILON {
                    (v - v_min) / v_span
                } else {
                    0.0
                };
                mesh.uvs.push(Point2::new(su, sv));
                vertex_map.insert(idx, new_idx);
                new_idx
            };
            tri[slot] = mesh_idx;
        }
        mesh.indices.push(tri);
    }

    if mesh.indices.is_empty() {
        return Err(TessellationError::Failed(
            "trimmed tessellation produced no triangles (trim region may be empty)".into(),
        )
        .into());
    }

    Ok(mesh)
}

fn validate_options(options: &SurfaceTessellationOptions) -> Result<()> {
    if options.normal_tolerance <= 0.0 {
        return Err(TessellationError::InvalidParameters(
            "normal_tolerance must be strictly positive".to_owned(),
        )
        .into());
    }
    if options.min_divisions == 0 || options.min_divisions > options.max_divisions {
        return Err(TessellationError::InvalidParameters(
            "require 1 <= min_divisions <= max_divisions".to_owned(),
        )
        .into());
    }
    Ok(())
}

/// Samples a trim loop into a closed UV polyline with no duplicate vertices.
///
/// Each pcurve is sampled at [`PCURVE_SEGMENT_SAMPLES`] interior+start points
/// (the per-segment tail is shared with the next segment's head), then the whole
/// polyline is deduplicated, including the wrap-around closing point.
fn sample_loop(loop_: &TrimLoop) -> Result<Vec<Point2>> {
    let mut pts: Vec<Point2> = Vec::new();
    for curve in &loop_.curves {
        append_curve_samples(curve, &mut pts)?;
    }
    dedup_closed(&mut pts);

    if pts.len() < 3 {
        return Err(TessellationError::Failed(
            "trim loop degenerated to fewer than 3 distinct UV points".into(),
        )
        .into());
    }
    Ok(pts)
}

/// Appends the start point and interior samples of `curve` (dropping the tail,
/// which is the next segment's start or the loop closure).
///
/// A **degree-1 polyline** is sampled at exactly its control points (no
/// chord-based resampling, no densification): for degree 1 the control points
/// lie on the curve, so this reproduces the curve losslessly. This determinism
/// is what makes two faces whose trim loops were built from the *same* trace
/// points emit identical 3D vertices along a shared ring (the punch/band cut
/// rings), eliminating the T-junction cracks at the hole rim. Higher-degree
/// trim curves (rational arcs such as `circle_uv`) keep the fixed per-segment
/// adaptive sampling, which they need to approximate their curvature.
fn append_curve_samples(curve: &NurbsCurve2D, out: &mut Vec<Point2>) -> Result<()> {
    if curve.degree() == 1 {
        // Control points of a degree-1 curve are exactly its polyline vertices.
        // Drop the tail control point: it is the next segment's start (or, for
        // the final segment, the loop closure handled by `dedup_closed`).
        let cps = curve.control_points();
        for p in &cps[..cps.len() - 1] {
            out.push(*p);
        }
        return Ok(());
    }

    let (t0, t1) = curve.parameter_domain();
    for i in 0..PCURVE_SEGMENT_SAMPLES {
        #[allow(clippy::cast_precision_loss)]
        let frac = i as f64 / PCURVE_SEGMENT_SAMPLES as f64;
        let t = t0 + frac * (t1 - t0);
        out.push(curve.point_at(t)?);
    }
    Ok(())
}

/// Removes consecutive near-duplicate points and the closing duplicate (last
/// point coincident with the first).
fn dedup_closed(pts: &mut Vec<Point2>) {
    pts.dedup_by(|a, b| (*a - *b).norm() < MERGE_EPS);
    while pts.len() >= 2 {
        let first = pts[0];
        let last = pts[pts.len() - 1];
        if (first - last).norm() < MERGE_EPS {
            pts.pop();
        } else {
            break;
        }
    }
}

/// Inserts a closed polygon as constraint edges into the CDT.
///
/// Vertices are pre-sanitized (deduplicated) by [`sample_loop`], so spade never
/// receives a duplicate constraint vertex; constraints between coincident
/// handles are skipped defensively.
fn insert_constraint_loop(
    cdt: &mut ConstrainedDelaunayTriangulation<SpadePoint2<f64>>,
    points: &[Point2],
) -> Result<()> {
    if points.len() < 3 {
        return Err(
            TessellationError::Failed("constraint loop needs at least 3 points".into()).into(),
        );
    }

    let mut handles = Vec::with_capacity(points.len());
    for p in points {
        let h = cdt
            .insert(SpadePoint2::new(p.x, p.y))
            .map_err(|e: InsertionError| TessellationError::Failed(format!("CDT insert: {e}")))?;
        handles.push(h);
    }

    for i in 0..handles.len() {
        let from = handles[i];
        let to = handles[(i + 1) % handles.len()];
        if from != to {
            // `try_add_constraint` resolves (rather than panics on) any
            // intersection with an existing constraint; the loops are
            // pre-sanitized so this is a defensive guard.
            let _ = cdt.try_add_constraint(from, to);
        }
    }
    Ok(())
}

/// UV distance below which a Steiner candidate counts as lying ON a
/// constraint segment (and is skipped — see the Steiner insertion loop).
const SEGMENT_SKIP_EPS: f64 = 1e-9;

/// Whether `p` lies within [`SEGMENT_SKIP_EPS`] of any constraint segment of
/// the outer loop or a hole loop.
fn near_any_segment(p: &Point2, outer: &[Point2], holes: &[Vec<Point2>]) -> bool {
    let near_loop = |poly: &[Point2]| -> bool {
        let count = poly.len();
        (0..count).any(|i| {
            let seg_a = poly[i];
            let seg_b = poly[(i + 1) % count];
            let dir = seg_b - seg_a;
            let len_sq = dir.norm_squared();
            let dist_sq = if len_sq < 1e-30 {
                (*p - seg_a).norm_squared()
            } else {
                let frac = ((*p - seg_a).dot(&dir) / len_sq).clamp(0.0, 1.0);
                (*p - (seg_a + dir * frac)).norm_squared()
            };
            dist_sq < SEGMENT_SKIP_EPS * SEGMENT_SKIP_EPS
        })
    };
    near_loop(outer) || holes.iter().any(|h| near_loop(h))
}

/// Tests whether `p` is inside the trim region: inside the outer loop and
/// outside every hole.
fn point_in_region(p: &Point2, outer: &[Point2], holes: &[Vec<Point2>]) -> bool {
    if !point_in_polygon(p, outer) {
        return false;
    }
    for hole in holes {
        if point_in_polygon(p, hole) {
            return false;
        }
    }
    true
}

/// Even-odd ray-cast point-in-polygon test in UV.
fn point_in_polygon(p: &Point2, poly: &[Point2]) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let pi = poly[i];
        let pj = poly[j];
        let intersects = (pi.y > p.y) != (pj.y > p.y)
            && p.x < (pj.x - pi.x) * (p.y - pi.y) / (pj.y - pi.y) + pi.x;
        if intersects {
            inside = !inside;
        }
        j = i;
    }
    inside
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::nurbs::{KnotVector, NurbsCurve2D};
    use crate::math::Point3;
    use crate::operations::creation::MakeNurbsFace;
    use crate::operations::query::ClosestPointOnSurface;
    use crate::tessellation::{TessellateFace, TessellationParams};
    use crate::topology::TopologyStore;

    /// Bilinear planar patch over [0,1]x[0,1] in the z=0 plane.
    fn unit_patch() -> NurbsSurface {
        NurbsSurface::from_unweighted(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
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

    fn line2d(a: (f64, f64), b: (f64, f64)) -> NurbsCurve2D {
        NurbsCurve2D::from_unweighted(
            vec![Point2::new(a.0, a.1), Point2::new(b.0, b.1)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap()
    }

    fn rect_loop(u0: f64, v0: f64, u1: f64, v1: f64, ccw: bool) -> TrimLoop {
        let c = if ccw {
            [(u0, v0), (u1, v0), (u1, v1), (u0, v1)]
        } else {
            [(u0, v0), (u0, v1), (u1, v1), (u1, v0)]
        };
        TrimLoop::new(vec![
            line2d(c[0], c[1]),
            line2d(c[1], c[2]),
            line2d(c[2], c[3]),
            line2d(c[3], c[0]),
        ])
    }

    fn circle_loop_cw(cx: f64, cy: f64, r: f64, n: usize) -> TrimLoop {
        // Clockwise polygonal hole approximating a circle.
        let mut curves = Vec::with_capacity(n);
        for i in 0..n {
            #[allow(clippy::cast_precision_loss)]
            let a0 = -std::f64::consts::TAU * (i as f64) / (n as f64);
            #[allow(clippy::cast_precision_loss)]
            let a1 = -std::f64::consts::TAU * ((i + 1) as f64) / (n as f64);
            let p0 = (cx + r * a0.cos(), cy + r * a0.sin());
            let p1 = (cx + r * a1.cos(), cy + r * a1.sin());
            curves.push(line2d(p0, p1));
        }
        TrimLoop::new(curves)
    }

    #[test]
    fn full_domain_trim_covers_whole_patch() {
        let surface = unit_patch();
        let trim = FaceTrim::new(rect_loop(0.0, 0.0, 1.0, 1.0, true), vec![]);
        let mesh =
            tessellate_trimmed_nurbs_face(&surface, &trim, &SurfaceTessellationOptions::default())
                .unwrap();
        assert!(!mesh.indices.is_empty());
        // Every vertex lies on the patch (z = 0).
        for vtx in &mesh.vertices {
            assert!(vtx.z.abs() < 1e-9, "vertex off patch: z = {}", vtx.z);
        }
        // The mesh covers roughly the full [0,1]^2 area.
        let area: f64 = mesh
            .indices
            .iter()
            .map(|t| {
                let a = mesh.vertices[t[0] as usize];
                let b = mesh.vertices[t[1] as usize];
                let c = mesh.vertices[t[2] as usize];
                (b - a).cross(&(c - a)).norm() * 0.5
            })
            .sum();
        assert!((area - 1.0).abs() < 1e-6, "covered area = {area}");
    }

    #[test]
    fn circular_hole_excludes_interior() {
        let surface = unit_patch();
        let outer = rect_loop(0.0, 0.0, 1.0, 1.0, true);
        let hole = circle_loop_cw(0.5, 0.5, 0.2, 32);
        let trim = FaceTrim::new(outer, vec![hole]);
        let mesh =
            tessellate_trimmed_nurbs_face(&surface, &trim, &SurfaceTessellationOptions::default())
                .unwrap();

        // No triangle centroid falls inside the hole circle.
        for t in &mesh.indices {
            let a = mesh.vertices[t[0] as usize];
            let b = mesh.vertices[t[1] as usize];
            let c = mesh.vertices[t[2] as usize];
            let cx = (a.x + b.x + c.x) / 3.0;
            let cy = (a.y + b.y + c.y) / 3.0;
            let d = ((cx - 0.5).powi(2) + (cy - 0.5).powi(2)).sqrt();
            assert!(d > 0.2 - 1e-3, "centroid inside hole: d = {d}");
        }

        // Covered area is about the patch minus the hole disc.
        let area: f64 = mesh
            .indices
            .iter()
            .map(|t| {
                let a = mesh.vertices[t[0] as usize];
                let b = mesh.vertices[t[1] as usize];
                let c = mesh.vertices[t[2] as usize];
                (b - a).cross(&(c - a)).norm() * 0.5
            })
            .sum();
        let expected = 1.0 - std::f64::consts::PI * 0.2 * 0.2;
        assert!(
            (area - expected).abs() < 0.02,
            "area {area} not near {expected}"
        );
    }

    #[test]
    fn trimmed_vertices_lie_on_surface() {
        let mut store = TopologyStore::new();
        let surface = unit_patch();
        let outer = rect_loop(0.0, 0.0, 1.0, 1.0, true);
        let hole = circle_loop_cw(0.5, 0.5, 0.2, 24);
        let trim = FaceTrim::new(outer, vec![hole]);
        let face = MakeNurbsFace::new(surface)
            .with_trim(trim)
            .execute(&mut store)
            .unwrap();
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        for vtx in &mesh.vertices {
            let cp = ClosestPointOnSurface::new(face, *vtx)
                .execute(&store)
                .unwrap();
            assert!(
                cp.distance < 1e-6,
                "vertex off surface: d = {}",
                cp.distance
            );
        }
    }

    #[test]
    fn trimmed_mesh_is_manifold_along_boundaries() {
        let surface = unit_patch();
        let outer = rect_loop(0.0, 0.0, 1.0, 1.0, true);
        let hole_segments = 32;
        let hole = circle_loop_cw(0.5, 0.5, 0.2, hole_segments);
        let trim = FaceTrim::new(outer, vec![hole]);
        let mesh =
            tessellate_trimmed_nurbs_face(&surface, &trim, &SurfaceTessellationOptions::default())
                .unwrap();

        // Undirected edge -> use-count over every triangle.
        let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
        for tri in &mesh.indices {
            for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }

        // Manifold: every edge is shared by exactly 1 (boundary) or 2 (interior)
        // triangles; nothing is referenced 3+ times.
        for (&(a, b), &c) in &counts {
            assert!(
                c == 1 || c == 2,
                "edge ({a},{b}) used {c} times (expected 1 or 2)"
            );
        }

        // Boundary edges (used exactly once) cover both the hole ring and the
        // outer rectangle ring, so there are at least as many as the hole's
        // polyline segment count.
        let boundary = counts.values().filter(|&&c| c == 1).count();
        assert!(
            boundary >= hole_segments,
            "boundary edge count {boundary} below hole segment count {hole_segments}"
        );
    }

    #[test]
    fn hole_larger_than_domain_is_rejected() {
        let surface = unit_patch();
        let outer = rect_loop(0.0, 0.0, 1.0, 1.0, true);
        // Hole covering the whole domain leaves no interior triangles.
        let hole = circle_loop_cw(0.5, 0.5, 2.0, 16);
        let trim = FaceTrim::new(outer, vec![hole]);
        let result =
            tessellate_trimmed_nurbs_face(&surface, &trim, &SurfaceTessellationOptions::default());
        assert!(result.is_err(), "oversized hole must yield an error");
    }
}
