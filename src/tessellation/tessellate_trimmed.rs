//! Constrained-Delaunay tessellation of trimmed NURBS faces.
//!
//! The trim loops are sampled into UV polylines and inserted into a spade
//! constrained Delaunay triangulation as constraint edges. Grid points from the
//! shared adaptive refinement (see [`super::tessellate_nurbs`]) that fall
//! strictly inside the trim region are added as Steiner points. Triangles whose
//! centroid lies outside the trim region (outside the outer loop or inside a
//! hole) are discarded, and the surviving vertices are mapped through the
//! surface to 3D.

use std::collections::HashMap;

use spade::{
    ConstrainedDelaunayTriangulation, InsertionError, Point2 as SpadePoint2, Triangulation,
};

use crate::error::{Result, TessellationError};
use crate::geometry::nurbs::{NurbsCurve2D, NurbsSurface};
use crate::geometry::surface::Surface;
use crate::math::{Point2, Vector3};
use crate::topology::{FaceTrim, TrimLoop};

use super::tessellate_nurbs::adaptive_grid_parameters;
use super::{SurfaceTessellationOptions, TriangleMesh};

/// Number of samples per pcurve segment when converting a trim loop to a UV
/// polyline. Trim loops in P5 are low-degree (lines / rational arcs), so a fixed
/// per-segment sampling is sufficient and avoids the spade duplicate-vertex
/// hazards that adaptive sampling near-coincident points would create.
const PCURVE_SEGMENT_SAMPLES: usize = 32;

/// Minimum UV separation between consecutive polyline points. Points closer than
/// this are merged so spade never receives duplicate constraint vertices.
const MERGE_EPS: f64 = 1e-9;

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

    // 1. Sample each loop into a UV polyline.
    let outer = sample_loop(&trim.outer)?;
    let holes: Vec<Vec<Point2>> = trim.holes.iter().map(sample_loop).collect::<Result<_>>()?;

    // 2. Adaptive UV grid (shared with the full-domain tessellator).
    let (u_params, v_params) = adaptive_grid_parameters(surface, options);

    // 3. Constrained Delaunay triangulation.
    let mut cdt = ConstrainedDelaunayTriangulation::<SpadePoint2<f64>>::new();

    insert_constraint_loop(&mut cdt, &outer)?;
    for hole in &holes {
        insert_constraint_loop(&mut cdt, hole)?;
    }

    // Steiner points: grid samples strictly inside the trim region.
    for &u in &u_params {
        for &v in &v_params {
            let p = Point2::new(u, v);
            if point_in_region(&p, &outer, &holes) {
                // Ignore individual insertion failures (e.g. a point landing on
                // an existing constraint vertex); the constraint loops already
                // pin the boundary.
                let _ = cdt.insert(SpadePoint2::new(u, v));
            }
        }
    }

    // 4. Keep triangles whose centroid is inside the trim region; map to 3D.
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_domain();
    let u_span = u_max - u_min;
    let v_span = v_max - v_min;

    let mut mesh = TriangleMesh::default();
    let mut vertex_map: HashMap<usize, u32> = HashMap::new();

    for face in cdt.inner_faces() {
        let verts = face.vertices();
        let centroid = {
            let p0 = verts[0].position();
            let p1 = verts[1].position();
            let p2 = verts[2].position();
            Point2::new((p0.x + p1.x + p2.x) / 3.0, (p0.y + p1.y + p2.y) / 3.0)
        };
        if !point_in_region(&centroid, &outer, &holes) {
            continue;
        }

        let mut tri = [0u32; 3];
        for (slot, vh) in verts.iter().enumerate() {
            let idx = vh.fix().index();
            let mesh_idx = if let Some(&existing) = vertex_map.get(&idx) {
                existing
            } else {
                let pos = vh.position();
                let (u, v) = (pos.x, pos.y);
                let p3 = surface.point_at(u, v)?;
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
fn append_curve_samples(curve: &NurbsCurve2D, out: &mut Vec<Point2>) -> Result<()> {
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
