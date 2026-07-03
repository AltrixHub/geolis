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
//! `u` domain. The SSI marcher wraps the tool's periodic `u` direction, so the
//! loop arrives genuinely `closed` with its tool `u` kept wrapped into
//! `[u0, u1]` and with **exact seam samples** (the crossing point emitted at
//! both `u0` and `u1`), so each trace spans the **full** `u` domain at a
//! roughly constant `v`. The entry loop sits at a lower mean `v` than the exit
//! loop (the loops are pre-sorted by mean `v` in [`super::loops`]). Stitching
//!
//! ```text
//!   entry trace  (u increasing)
//!   -> exit trace (u decreasing)
//!   -> close
//! ```
//!
//! yields a ribbon polygon that is simple (non-self-intersecting) in the
//! unrolled rectangle, so the generic trimmed CDT meshes it without a seam cut.
//! Because the traces reach `u0` and `u1` exactly, the ribbon's left (`u0`)
//! and right (`u1`) closing edges land on the same seam azimuth and coincide in
//! 3D, covering the seam wedge. The two rings (entry/exit) share the exact seam
//! samples with the punched target faces, so the hole rim tessellates
//! conformally with no slit at the seam. If a marcher seam sample did not
//! converge, the ribbon degrades to the marched span (a sub-step gap at the
//! seam) — the honest fallback.
//!
//! ### Orientation
//!
//! Subtract pushes the band normals INTO the hole, so the band face is built
//! with `same_sense = false`.

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{KnotVector, NurbsCurve2D, NurbsSurface};
use crate::math::Point2;
use crate::topology::{FaceData, FaceId, FaceSurface, FaceTrim, TopologyStore, TrimLoop, WireId};

use super::loops::ToolFaceCut;

/// The two hole-ring wires shared with the punched target faces for one tool
/// side face: the entry ring (lower mean v) and the exit ring (upper mean v).
///
/// These are the exact [`WireId`]s returned by [`super::punch::punch_loop`] for
/// the same tool face's two loops, so the band face shares its boundary edges
/// with the punched target faces (correct `BRep` adjacency).
#[derive(Debug, Clone, Copy)]
pub(crate) struct BandRingWires {
    /// Entry ring wire (matches `cut.loops[0]`, the lower-v loop).
    pub entry: WireId,
    /// Exit ring wire (matches `cut.loops[1]`, the upper-v loop).
    pub exit: WireId,
}

/// Builds the band (hole-wall) face for one tool side face from its two cut
/// loops, and returns the new face's id.
///
/// `rings` carries the two hole-ring wires already created by the punch step for
/// this tool face's loops; the band face reuses them as its boundary so it
/// shares edges/wires with the punched target faces instead of fabricating a new
/// full-surface boundary.
///
/// # Errors
///
/// Returns an error if the tool face is not a NURBS face or the stitched band
/// polygon degenerates (fewer than 3 distinct UV points).
pub(crate) fn build_band_face(
    store: &mut TopologyStore,
    cut: &ToolFaceCut,
    rings: BandRingWires,
) -> Result<FaceId> {
    // Subtract: band normals point INTO the hole (`same_sense = false`).
    build_band_face_oriented(store, cut, rings, false)
}

/// Builds the band (hole-wall / plug-wall) face for one tool side face with an
/// explicit `same_sense`.
///
/// The band region (the tool side surface strip between the entry and exit
/// loops) is identical for the subtract and intersect through-cuts; only the
/// normal orientation differs. Subtract points the band normals into the hole
/// (`same_sense = false`); intersect points them outward from the kept plug
/// (`same_sense = true`).
///
/// # Errors
///
/// Returns an error if the tool face is not a NURBS face or the stitched band
/// polygon degenerates (fewer than 3 distinct UV points).
pub(crate) fn build_band_face_oriented(
    store: &mut TopologyStore,
    cut: &ToolFaceCut,
    rings: BandRingWires,
    same_sense: bool,
) -> Result<FaceId> {
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

    // The band's real boundary is the two SSI rings (entry + exit). Share the
    // exact wires the punch step attached to the target faces so the band face
    // has correct BRep adjacency (no fabricated full-surface seam wire).
    Ok(store.add_face(FaceData {
        surface: FaceSurface::Nurbs(surface),
        outer_wire: rings.entry,
        inner_wires: vec![rings.exit],
        same_sense,
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
    use crate::operations::boolean::nurbs::punch::punch_loop;
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
    /// tool surface). The two hole loops are punched first (as the real pipeline
    /// does) so the band shares their ring wires.
    fn band_face(radius: f64) -> (TopologyStore, FaceId, NurbsSurface) {
        let (store, band, surf, _rings) = band_face_with_rings(radius);
        (store, band, surf)
    }

    /// Like [`band_face`] but also returns the entry/exit ring wires the band
    /// shares with the punched faces.
    fn band_face_with_rings(radius: f64) -> (TopologyStore, FaceId, NurbsSurface, BandRingWires) {
        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -1.5), radius, 5.0)
            .execute(&mut store)
            .unwrap();
        let target = collect_nurbs_faces(&store, &solid_faces(&store, slab));
        let tool = collect_nurbs_faces(&store, &solid_faces(&store, tube));
        let cuts = extract_cut_loops(&target, &tool).unwrap();
        let tool_surf = tool[0].1.clone();
        // Punch both loops (entry then exit) and share their ring wires.
        let entry = punch_loop(&mut store, &cuts[0].loops[0]).unwrap();
        let exit = punch_loop(&mut store, &cuts[0].loops[1]).unwrap();
        let rings = BandRingWires { entry, exit };
        let band = build_band_face(&mut store, &cuts[0], rings).unwrap();
        (store, band, tool_surf, rings)
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

    #[test]
    fn band_boundary_is_exactly_the_two_ring_wires() {
        let (store, band, _surf, rings) = band_face_with_rings(0.7);
        let face = store.face(band).unwrap();
        // The band's outer wire is the entry ring; its single inner wire is the
        // exit ring — the exact wires the punch step created.
        assert_eq!(face.outer_wire, rings.entry, "outer wire = entry ring");
        assert_eq!(face.inner_wires.len(), 1, "exactly one inner wire");
        assert_eq!(face.inner_wires[0], rings.exit, "inner wire = exit ring");
        // The two rings are distinct wires (entry != exit).
        assert_ne!(rings.entry, rings.exit, "entry and exit rings differ");
    }
}
