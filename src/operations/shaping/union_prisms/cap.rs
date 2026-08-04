//! Horizontal caps of the z-slab decomposition.
//!
//! A slab is capped wherever its neighbour's coverage stops: a
//! down-facing cap over `region(i) \ region(i - 1)` at its base and an
//! up-facing cap over `region(i) \ region(i + 1)` at its top, both region
//! differences taken with the crate's 2D boolean machinery. The slab
//! below the first and above the last is the empty region, so the body is
//! closed at the extremes.
//!
//! # One computation, two consumers
//!
//! These faces are the body's ONLY horizontal surfaces, so they are also
//! the only place its outline carries a horizontal edge. They are
//! therefore computed here once and consumed twice: [`super::mesh`]
//! triangulates them, and [`super::FusedPrisms::caps`] hands their rings
//! to a caller drawing the body. Splitting the two would let the drawn
//! outline and the meshed surface describe different bodies — which is
//! exactly what a caller drawing each slab's ring at both its ends gets:
//! one opening turns the whole group into three slabs, so the ring of
//! every element in it would be drawn again at that opening's sill and
//! head, where the material is continuous and the mesh has no face at
//! all.

use crate::error::Result;
use crate::operations::boolean_2d::{
    intersect_all_with_holes, subtract_all_with_holes, union_all_with_holes, PolygonWithHoles,
    WALL_EPS,
};

use super::{CapFace, CapFacing, PrismSlab};

/// Collects every horizontal cap of the contiguous slab list produced by
/// [`super::slab::slab_regions`].
///
/// Ordered bottom to top, a slab's base caps before its top caps. Slabs
/// with no material contribute nothing and are still counted as
/// neighbours — their empty region caps whatever adjoins them.
///
/// # Errors
///
/// Propagates arrangement-engine failures from the region differences.
pub(super) fn cap_faces(slabs: &[PrismSlab]) -> Result<Vec<CapFace>> {
    let mut caps = Vec::new();
    for (index, slab) in slabs.iter().enumerate() {
        if slab.faces().is_empty() {
            continue;
        }
        // Coverage below / above. Out-of-range neighbours are the empty
        // region, which closes the body at the extremes.
        let below = index.checked_sub(1).and_then(|i| slabs.get(i));
        let above = slabs.get(index + 1);

        for region in uncovered_by(slab, below)? {
            caps.push(CapFace::new(slab.z_base(), CapFacing::Down, region));
        }
        for region in uncovered_by(slab, above)? {
            caps.push(CapFace::new(slab.z_top(), CapFacing::Up, region));
        }
    }
    Ok(caps)
}

/// The part of `slab`'s cross-section its `neighbour` does not cover —
/// i.e. the area needing a horizontal cap. `None` is the empty region
/// outside the stack, which caps the whole cross-section.
///
/// # Shared-material shortcut
///
/// Two neighbouring slabs carved out of the SAME material (equal cover —
/// the identity [`super::slab::slab_regions`]'s per-cover memo is built
/// on) differ ONLY where their cuts differ:
///
/// ```text
/// slab_i \ slab_j = (F \ Ci) \ (F \ Cj) = (F \ Ci) ∩ (Cj \ Ci)
/// ```
///
/// so the cap is an arrangement over the CUTS — window-sized rectangles —
/// clipped to the slab, instead of one over the whole floor's plan
/// boundary against itself. When the neighbour's cuts are a subset of
/// this slab's, the difference is empty and no arrangement runs at all.
/// This is what keeps the fused body's cap stage off the project's total
/// outline complexity: a floor of walls has one cover, so every interior
/// cap between two of its slabs takes this path.
///
/// Covers that differ (a genuinely different material set above / below)
/// fall back to the direct region difference.
fn uncovered_by(slab: &PrismSlab, neighbour: Option<&PrismSlab>) -> Result<Vec<PolygonWithHoles>> {
    let Some(neighbour) = neighbour else {
        return Ok(slab.faces().to_vec());
    };
    if neighbour.faces().is_empty() {
        return Ok(slab.faces().to_vec());
    }
    if slab.cover() != neighbour.cover() {
        return region_difference(slab.faces(), neighbour.faces());
    }
    // Same material: every cut the neighbour carries that this slab also
    // carries removes the same area from both, so only the neighbour's
    // EXTRA cuts can expose a cap.
    if neighbour
        .cut_keys()
        .iter()
        .all(|key| slab.cut_keys().contains(key))
    {
        return Ok(Vec::new());
    }
    let exposed = region_difference(neighbour.cuts(), slab.cuts())?;
    // `exposed` is disjoint from this slab's own cuts by construction,
    // so clipping it to the slab is the same as clipping it to the
    // UNCARVED material — and that material is a union of per-profile
    // parts, which lets each clip run against the two or three parts the
    // exposed region actually touches instead of the whole floor.
    let mut pieces = Vec::new();
    for region in &exposed {
        let Some(region_bounds) = Bounds::of(region) else {
            // No usable extent: clip against the fused faces rather than
            // guess which parts are near.
            pieces.extend(intersect_all_with_holes(region, slab.faces())?);
            continue;
        };
        for part in slab.parts() {
            if Bounds::of(part).is_some_and(|part_bounds| !region_bounds.overlaps(&part_bounds)) {
                continue;
            }
            pieces.extend(intersect_all_with_holes(
                region,
                std::slice::from_ref(part),
            )?);
        }
    }
    // Parts overlap wherever elements meet, so the clipped pieces can
    // too; one union folds them back into disjoint cap faces.
    if pieces.len() < 2 {
        return Ok(pieces);
    }
    Ok(union_all_with_holes(&pieces)?.faces)
}

/// Axis-aligned extent of a region's outer ring, grown by [`WALL_EPS`]
/// so a shared boundary counts as an overlap.
#[derive(Clone, Copy)]
struct Bounds {
    min_x: f64,
    min_y: f64,
    max_x: f64,
    max_y: f64,
}

impl Bounds {
    /// `None` when the ring is empty or carries a non-finite coordinate —
    /// the caller must then fall back to the exact path rather than
    /// reject anything.
    fn of(region: &PolygonWithHoles) -> Option<Self> {
        let mut bounds = Self {
            min_x: f64::INFINITY,
            min_y: f64::INFINITY,
            max_x: f64::NEG_INFINITY,
            max_y: f64::NEG_INFINITY,
        };
        for &(x, y) in &region.outer {
            if !x.is_finite() || !y.is_finite() {
                return None;
            }
            bounds.min_x = bounds.min_x.min(x);
            bounds.min_y = bounds.min_y.min(y);
            bounds.max_x = bounds.max_x.max(x);
            bounds.max_y = bounds.max_y.max(y);
        }
        (bounds.min_x <= bounds.max_x).then_some(bounds)
    }

    fn overlaps(&self, other: &Self) -> bool {
        self.min_x - WALL_EPS <= other.max_x
            && other.min_x - WALL_EPS <= self.max_x
            && self.min_y - WALL_EPS <= other.max_y
            && other.min_y - WALL_EPS <= self.max_y
    }
}

/// `⋃ region ∩ ¬⋃ other` — the part of one region set the other does not
/// cover.
fn region_difference(
    region: &[PolygonWithHoles],
    other: &[PolygonWithHoles],
) -> Result<Vec<PolygonWithHoles>> {
    if other.is_empty() {
        return Ok(region.to_vec());
    }
    let mut difference = Vec::with_capacity(region.len());
    for face in region {
        difference.extend(subtract_all_with_holes(face.clone(), other)?);
    }
    Ok(difference)
}
