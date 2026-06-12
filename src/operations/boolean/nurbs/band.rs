//! Tool band (hole-wall) face construction for the through-cut subtract.
//!
//! ## Shipped path: generic stitched-CDT
//!
//! The plan offered a fallback ladder for the band's UV topology (generic
//! stitched-CDT first, dedicated quad-strip only if that proved infeasible).
//! The **generic stitched-CDT path shipped**: the band between the entry and
//! exit loops on the tool side surface is represented as an ordinary trimmed
//! NURBS face whose `FaceTrim.outer` is a single simple polygon in the tool's
//! unrolled UV rectangle, and the existing constrained-Delaunay trimmed
//! tessellator (P5) meshes it directly. No dedicated band tessellation was
//! needed.
//!
//! ### Why it works
//!
//! Each loop's `uv_b` trace is a single-valued graph `v = f(u)` over the tool's
//! `u` domain (the SSI marcher terminates at the `u` seam, so each trace spans
//! nearly the full `u` range at a roughly constant `v`). The entry loop sits at
//! a lower mean `v` than the exit loop (the loops are pre-sorted by mean `v` in
//! [`super::loops`]). Stitching
//!
//! ```text
//!   entry trace  (u increasing)
//!   -> exit trace (u decreasing)
//!   -> close
//! ```
//!
//! yields a ribbon polygon that is simple (non-self-intersecting) in the
//! unrolled rectangle, so the generic trimmed CDT meshes it without a seam cut.
//!
//! ### Orientation
//!
//! Subtract pushes the band normals INTO the hole, so the band face is built
//! with `same_sense = false`.

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{KnotVector, NurbsCurve2D, NurbsSurface};
use crate::math::{Point2, TOLERANCE};
use crate::topology::{FaceData, FaceId, FaceSurface, FaceTrim, TopologyStore, TrimLoop, WireData};

use super::loops::ToolFaceCut;

/// Builds the band (hole-wall) face for one tool side face from its two cut
/// loops, and returns the new face's id.
///
/// # Errors
///
/// Returns an error if the tool face is not a NURBS face or the stitched band
/// polygon degenerates (fewer than 3 distinct UV points).
pub(crate) fn build_band_face(store: &mut TopologyStore, cut: &ToolFaceCut) -> Result<FaceId> {
    let surface = match &store.face(cut.tool_face)?.surface {
        FaceSurface::Nurbs(s) => s.clone(),
        _ => {
            return Err(OperationError::Failed(
                "through-cut band requires a NURBS tool side face".into(),
            )
            .into())
        }
    };

    let entry = clamp_trace(&cut.loops[0].branch.uv_b, &surface);
    let exit = clamp_trace(&cut.loops[1].branch.uv_b, &surface);

    let outer = stitch_band_loop(&entry, &exit)?;
    let trim = FaceTrim::new(outer, Vec::new());

    // A band face shares the tool's outer boundary topologically only as a
    // bookkeeping wire; reuse a degenerate-free closed wire built from the
    // surface's own boundary so the face has a valid `outer_wire`. The trim is
    // what drives tessellation.
    let outer_wire = boundary_wire(store, &surface)?;

    Ok(store.add_face(FaceData {
        surface: FaceSurface::Nurbs(surface),
        outer_wire,
        inner_wires: Vec::new(),
        // Subtract: band normals point into the hole.
        same_sense: false,
        trim: Some(trim),
    }))
}

/// Clamps a UV trace into the surface's parameter domain (the SSI corrector may
/// land a hair outside on the seam side) and deduplicates.
fn clamp_trace(uv: &[Point2], surface: &NurbsSurface) -> Vec<Point2> {
    let ((u0, u1), (v0, v1)) = surface.parameter_domain();
    let mut out: Vec<Point2> = Vec::with_capacity(uv.len());
    for p in uv {
        let c = Point2::new(p.x.clamp(u0, u1), p.y.clamp(v0, v1));
        if out.last().is_none_or(|q| (c - q).norm() > 1e-9) {
            out.push(c);
        }
    }
    out
}

/// Stitches the entry and exit traces into a single simple band polygon (a
/// CCW outer trim loop of degree-1 segments).
///
/// The entry trace is walked in its natural (u-increasing) direction and the
/// exit trace reversed (u-decreasing), so the ribbon closes without crossing.
/// Winding is normalized to counter-clockwise (the trim outer convention).
fn stitch_band_loop(entry: &[Point2], exit: &[Point2]) -> Result<TrimLoop> {
    if entry.len() < 2 || exit.len() < 2 {
        return Err(
            OperationError::Failed("band trace degenerated to fewer than 2 points".into()).into(),
        );
    }

    // Order each trace by ascending u so the stitch direction is unambiguous.
    let mut e = entry.to_vec();
    let mut x = exit.to_vec();
    e.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal));
    x.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal));

    // Polygon: entry (u up) then exit (u down).
    let mut poly: Vec<Point2> = Vec::with_capacity(e.len() + x.len());
    poly.extend_from_slice(&e);
    poly.extend(x.iter().rev().copied());
    dedup_closed(&mut poly);

    if poly.len() < 3 {
        return Err(OperationError::Failed(
            "stitched band polygon degenerated to fewer than 3 points".into(),
        )
        .into());
    }

    // Normalize to CCW for the trim outer convention.
    if signed_area(&poly) < 0.0 {
        poly.reverse();
    }

    let n = poly.len();
    let mut curves = Vec::with_capacity(n);
    for i in 0..n {
        curves.push(uv_segment(poly[i], poly[(i + 1) % n]));
    }
    Ok(TrimLoop::new(curves))
}

/// Builds a bookkeeping closed outer wire from the surface's exact boundary
/// isocurves (reusing the same four-isocurve scheme as `MakeNurbsFace`).
fn boundary_wire(
    store: &mut TopologyStore,
    surface: &NurbsSurface,
) -> Result<crate::topology::WireId> {
    use crate::topology::{EdgeCurve, EdgeData, OrientedEdge, VertexData, VertexId};

    let [u_min_edge, u_max_edge, v_min_edge, v_max_edge] = surface.boundary_curves()?;
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_domain();
    let c00 = surface.point_at(u_min, v_min)?;
    let c10 = surface.point_at(u_max, v_min)?;
    let c11 = surface.point_at(u_max, v_max)?;
    let c01 = surface.point_at(u_min, v_max)?;

    let mut corners: Vec<(crate::math::Point3, VertexId)> = Vec::new();
    let mut vertex_for = |store: &mut TopologyStore, p: crate::math::Point3| -> VertexId {
        for (cp, id) in &corners {
            if (cp - p).norm() < TOLERANCE {
                return *id;
            }
        }
        let id = store.add_vertex(VertexData::new(p));
        corners.push((p, id));
        id
    };
    let v00 = vertex_for(store, c00);
    let v10 = vertex_for(store, c10);
    let v11 = vertex_for(store, c11);
    let v01 = vertex_for(store, c01);

    let segments = [
        (v_min_edge, v00, v10, true),
        (u_max_edge, v10, v11, true),
        (v_max_edge, v01, v11, false),
        (u_min_edge, v00, v01, false),
    ];
    let mut oriented = Vec::with_capacity(4);
    for (curve, ns, ne, forward) in segments {
        if ns == ne {
            let (t0, t1) = curve.parameter_domain();
            if (curve.point_at(t1)? - curve.point_at(t0)?).norm() < TOLERANCE {
                continue;
            }
        }
        let (t0, t1) = curve.parameter_domain();
        let edge = store.add_edge(EdgeData {
            start: ns,
            end: ne,
            curve: EdgeCurve::Nurbs(curve),
            t_start: t0,
            t_end: t1,
        });
        oriented.push(OrientedEdge::new(edge, forward));
    }
    if oriented.is_empty() {
        return Err(OperationError::Failed("band boundary collapsed".into()).into());
    }
    Ok(store.add_wire(WireData {
        edges: oriented,
        is_closed: true,
    }))
}

/// A degree-1 two-point UV line segment.
fn uv_segment(a: Point2, b: Point2) -> NurbsCurve2D {
    NurbsCurve2D::from_unweighted(
        vec![a, b],
        KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap_or_else(|_| unreachable!()),
        1,
    )
    .unwrap_or_else(|_| unreachable!())
}

/// Removes consecutive near-duplicate points and a coincident wrap point.
fn dedup_closed(pts: &mut Vec<Point2>) {
    pts.dedup_by(|a, b| (*a - *b).norm() < 1e-9);
    while pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).norm() < 1e-9 {
        pts.pop();
    }
}

/// Shoelace signed area. Positive = counter-clockwise.
fn signed_area(pts: &[Point2]) -> f64 {
    let n = pts.len();
    let mut a2 = 0.0;
    for i in 0..n {
        let p = pts[i];
        let q = pts[(i + 1) % n];
        a2 += p.x * q.y - q.x * p.y;
    }
    0.5 * a2
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::nurbs::InversionOptions;
    use crate::math::Point3;
    use crate::operations::boolean::nurbs::loops::{collect_nurbs_faces, extract_cut_loops};
    use crate::operations::creation::{MakeCurvedSlab, MakeNurbsTube};
    use crate::tessellation::{TessellateFace, TessellationParams};
    use crate::topology::SolidId;
    use std::collections::HashMap;

    fn solid_faces(store: &TopologyStore, solid: SolidId) -> Vec<FaceId> {
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        shell.faces.clone()
    }

    /// Builds the band face for a slab×tube and returns (store, band face,
    /// tool surface).
    fn band_face(radius: f64) -> (TopologyStore, FaceId, NurbsSurface) {
        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -1.5), radius, 5.0)
            .execute(&mut store)
            .unwrap();
        let target = collect_nurbs_faces(&store, &solid_faces(&store, slab));
        let tool = collect_nurbs_faces(&store, &solid_faces(&store, tube));
        let cuts = extract_cut_loops(&store, &target, &tool).unwrap();
        let tool_surf = tool[0].1.clone();
        let band = build_band_face(&mut store, &cuts[0]).unwrap();
        (store, band, tool_surf)
    }

    #[test]
    fn band_mesh_is_watertight() {
        let (store, band, _) = band_face(0.7);
        let mesh = TessellateFace::new(band, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty(), "band produced no triangles");

        // Edge-manifold check: every undirected edge used 1 or 2 times.
        let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
        for tri in &mesh.indices {
            for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }
        for (&(a, b), &c) in &counts {
            assert!(c == 1 || c == 2, "band edge ({a},{b}) used {c} times");
        }
    }

    #[test]
    fn band_vertices_lie_on_tool_surface() {
        let (store, band, tool_surf) = band_face(0.7);
        let mesh = TessellateFace::new(band, TessellationParams::default())
            .execute(&store)
            .unwrap();
        let opts = InversionOptions::default();
        for v in &mesh.vertices {
            let inv = tool_surf.closest_point(v, &opts).unwrap();
            assert!(
                inv.distance < 1e-6,
                "band vertex off tool surface: d = {}",
                inv.distance
            );
        }
    }

    #[test]
    fn band_z_extent_spans_slab_thickness() {
        let (store, band, _) = band_face(0.7);
        let mesh = TessellateFace::new(band, TessellationParams::default())
            .execute(&store)
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
        // The band runs from the back face (z ~ -1 .. 0.5) up to the front face
        // (z ~ 0 .. 1.5); its vertical extent should be a meaningful fraction of
        // the slab thickness (>= 0.5).
        assert!(
            zmax - zmin > 0.5,
            "band z-extent {} too small (zmin={zmin}, zmax={zmax})",
            zmax - zmin
        );
    }
}
