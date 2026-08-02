//! z-slab decomposition: the 2D half of [`super::UnionPrisms`].
//!
//! Every profile is flattened once into line-only rings, the z
//! breakpoints are collected, and each interval between two consecutive
//! breakpoints gets its cross-section computed by a union followed by a
//! subtract. Intervals with no material are kept in the returned vector
//! (with an empty face list) so the mesh stage can address a slab's
//! neighbours by index — [`super::UnionPrisms::execute`] drops them from
//! the public result.

use crate::error::Result;
use crate::geometry::pline::Pline;
use crate::operations::boolean_2d::{
    signed_area, subtract_all_with_holes, union_all_with_holes, Polygon, PolygonWithHoles,
    WALL_EPS, WALL_EPS_SQ,
};

use super::{PrismProfile, PrismRegion, PrismSlab, Z_EPS};

/// A profile with every ring flattened to line segments, ready for the
/// arrangement engine. Degenerate rings are already dropped.
struct FlatProfile {
    faces: Vec<PolygonWithHoles>,
    z_base: f64,
    z_top: f64,
    cuts: Vec<FlatCut>,
}

/// A cut with its region flattened.
struct FlatCut {
    region: PolygonWithHoles,
    z_base: f64,
    z_top: f64,
}

/// Computes the fused cross-section of every z interval, bottom to top.
///
/// The returned vector is contiguous in z: entry `i`'s `z_top` is entry
/// `i + 1`'s `z_base`, which is what lets the mesh stage read a slab's
/// neighbour coverage off `i ± 1`.
///
/// # Errors
///
/// Propagates arrangement-engine failures from the union / subtract
/// stages.
pub(super) fn slab_regions(
    profiles: &[PrismProfile],
    arc_tolerance: f64,
) -> Result<Vec<PrismSlab>> {
    let flats: Vec<FlatProfile> = profiles
        .iter()
        .map(|profile| flatten_profile(profile, arc_tolerance))
        .filter(|flat| !flat.faces.is_empty())
        .collect();

    let breakpoints = z_breakpoints(&flats);
    let mut slabs = Vec::with_capacity(breakpoints.len().saturating_sub(1));
    let mut fused_by_cover = FusedByCover::default();
    for pair in breakpoints.windows(2) {
        let (z_base, z_top) = (pair[0], pair[1]);
        let mid = 0.5 * (z_base + z_top);

        let mut cover: Vec<usize> = Vec::new();
        let mut cuts: Vec<PolygonWithHoles> = Vec::new();
        for (index, flat) in flats.iter().enumerate() {
            if !covers(flat.z_base, flat.z_top, mid) {
                continue;
            }
            cover.push(index);
            cuts.extend(
                flat.cuts
                    .iter()
                    .filter(|cut| covers(cut.z_base, cut.z_top, mid))
                    .map(|cut| cut.region.clone()),
            );
        }

        let faces = if cover.is_empty() {
            Vec::new()
        } else {
            carve(fused_by_cover.get(&cover, &flats)?, &cuts)?
        };
        slabs.push(PrismSlab::new(z_base, z_top, faces));
    }
    Ok(slabs)
}

/// Memo of `⋃ material` keyed by the set of profiles covering a slab.
///
/// A slab's material is decided entirely by *which* profiles reach it, so
/// two slabs with the same cover have the same union — and a floor's
/// slabs almost always do: every opening in the group opens two
/// breakpoints, and the walls that span them are the same walls. Without
/// this the identical union of every element on the floor was recomputed
/// once per slab, which is how one window turned the fusion into three
/// full-floor unions and ten windows into three unions of ten times the
/// material.
///
/// The cover is a sorted index list, so equality is exact — no geometry
/// is compared and nothing is approximated. Lookup is a linear scan
/// because the number of DISTINCT covers on one floor is small (it grows
/// only when elements have genuinely different z spans), and each entry
/// is one `Vec<usize>` comparison.
#[derive(Default)]
struct FusedByCover {
    entries: Vec<(Vec<usize>, Vec<PolygonWithHoles>)>,
}

impl FusedByCover {
    /// The fused material of `cover`, computing it on first request.
    ///
    /// `cover` is built by ascending index in [`slab_regions`], so it is
    /// already sorted and two equal covers compare equal.
    fn get(&mut self, cover: &[usize], flats: &[FlatProfile]) -> Result<&[PolygonWithHoles]> {
        let hit = self.entries.iter().position(|(key, _)| key == cover);
        let index = if let Some(index) = hit {
            index
        } else {
            let material: Vec<PolygonWithHoles> = cover
                .iter()
                .flat_map(|&i| flats[i].faces.iter().cloned())
                .collect();
            let fused = union_all_with_holes(&material)?.faces;
            self.entries.push((cover.to_vec(), fused));
            self.entries.len() - 1
        };
        Ok(&self.entries[index].1)
    }
}

/// `fused ∩ ¬⋃ cuts`, as typed face topology.
fn carve(fused: &[PolygonWithHoles], cuts: &[PolygonWithHoles]) -> Result<Vec<PolygonWithHoles>> {
    if cuts.is_empty() {
        return Ok(fused.to_vec());
    }
    // The fused faces are disjoint, so subtracting from each in turn is
    // the same region as one global subtraction — and lets a face no cut
    // reaches keep its vertices verbatim.
    let mut carved = Vec::with_capacity(fused.len());
    for face in fused {
        carved.extend(subtract_all_with_holes(face.clone(), cuts)?);
    }
    Ok(carved)
}

/// Sorted, [`Z_EPS`]-deduplicated elevations at which the cross-section
/// can change.
///
/// A cut's bounds are clamped into its owning profile's span: a cut
/// reaching past the prism it belongs to changes nothing there, so it
/// must not open an empty slab either.
fn z_breakpoints(flats: &[FlatProfile]) -> Vec<f64> {
    let mut zs: Vec<f64> = Vec::new();
    for flat in flats {
        zs.push(flat.z_base);
        zs.push(flat.z_top);
        for cut in &flat.cuts {
            zs.push(cut.z_base.clamp(flat.z_base, flat.z_top));
            zs.push(cut.z_top.clamp(flat.z_base, flat.z_top));
        }
    }
    // Every value is finite (validated on construction), so `total_cmp`
    // is a total order here and no NaN can reach the dedup.
    zs.sort_by(f64::total_cmp);
    zs.dedup_by(|later, kept| (*later - *kept).abs() <= Z_EPS);
    zs
}

/// Whether the half-open span `[lo, hi)` covers a slab whose midpoint is
/// `mid`. Slabs are longer than [`Z_EPS`] and their bounds are the
/// breakpoints themselves, so a strict test is unambiguous.
fn covers(lo: f64, hi: f64, mid: f64) -> bool {
    mid > lo && mid < hi
}

/// Flattens a profile's rings, dropping the degenerate ones.
fn flatten_profile(profile: &PrismProfile, arc_tolerance: f64) -> FlatProfile {
    FlatProfile {
        faces: profile
            .faces()
            .iter()
            .filter_map(|region| flatten_region(region, arc_tolerance))
            .collect(),
        z_base: profile.z_base(),
        z_top: profile.z_top(),
        cuts: profile
            .cuts()
            .iter()
            .filter_map(|cut| {
                flatten_region(cut.region(), arc_tolerance).map(|region| FlatCut {
                    region,
                    z_base: cut.z_base(),
                    z_top: cut.z_top(),
                })
            })
            .collect(),
    }
}

/// Flattens one region. A degenerate outer ring drops the whole region; a
/// degenerate hole drops only that hole.
fn flatten_region(region: &PrismRegion, arc_tolerance: f64) -> Option<PolygonWithHoles> {
    let outer = flatten_ring(region.outer(), arc_tolerance)?;
    let holes = region
        .holes()
        .iter()
        .filter_map(|hole| flatten_ring(hole, arc_tolerance))
        .collect();
    Some(PolygonWithHoles { outer, holes })
}

/// Flattens one ring into a [`Polygon`] — bulges tessellated, the closing
/// duplicate and any coincident neighbours removed.
///
/// Returns `None` for a ring that carries no area: fewer than three
/// distinct points, or a shoelace area within `WALL_EPS_SQ` of zero.
fn flatten_ring(ring: &Pline, arc_tolerance: f64) -> Option<Polygon> {
    let sampled = ring.to_points(arc_tolerance.max(WALL_EPS));
    let mut points: Vec<(f64, f64)> = Vec::with_capacity(sampled.len());
    for point in &sampled {
        if let Some(&(x, y)) = points.last() {
            if (point.x - x).hypot(point.y - y) < WALL_EPS {
                continue;
            }
        }
        points.push((point.x, point.y));
    }
    // `to_points` walks a closed ring back to its first vertex; the
    // `Polygon` contract closes the loop implicitly instead.
    if let (Some(&first), Some(&last)) = (points.first(), points.last()) {
        if points.len() >= 2 && (first.0 - last.0).hypot(first.1 - last.1) < WALL_EPS {
            points.pop();
        }
    }
    if points.len() < 3 || signed_area(&points).abs() <= WALL_EPS_SQ {
        return None;
    }
    Some(points)
}
