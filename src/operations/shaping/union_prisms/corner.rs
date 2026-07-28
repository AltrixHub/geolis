//! Vertical corner edges of the fused body.
//!
//! The horizontal slab rings alone leave a wall without side definition:
//! the fused outline says where the body starts and stops in plan, but
//! nothing marks the vertical arris at a wall end or a junction. This
//! module walks each slab's boundary rings and emits one vertical segment
//! per genuine corner.
//!
//! # Corner test
//!
//! At every ring vertex the incoming and outgoing edge directions are
//! compared; the vertex is a corner when the boundary turns by more than
//! the angle tolerance. That single test separates the two kinds of
//! vertex the arrangement engine produces:
//!
//! | Vertex kind | Turn |
//! |---|---|
//! | Wall end, L / T / X junction corner | ~90° — kept |
//! | Chord joint inside a flattened arc | bounded by the arc tolerance — dropped |
//! | Collinear seam where two inputs' edges met | ~0° — dropped |
//!
//! The last row is a free win: the union of two collinear input edges
//! keeps the vertex where they met, and the angle test removes it from
//! the drawn outline without touching the region topology.
//!
//! # No cross-slab merging
//!
//! A corner shared by two stacked slabs is emitted twice, once per slab.
//! The two positions come from separate arrangement runs, so identifying
//! them would need a tolerance-fuzzy match — and a caller drawing lines
//! cannot tell a stacked pair from a merged one anyway.

use crate::math::Point3;
use crate::operations::boolean_2d::{Polygon, WALL_EPS};

use super::{CornerEdge, PrismSlab};

/// Collects the vertical corner edges of every slab, bottom to top.
///
/// Within one slab the order is: face order, then outer ring before its
/// holes, then ring vertex order — the same walk
/// [`super::mesh::build_mesh`] uses for the wall quads, so the two agree
/// on which boundary a corner belongs to.
pub(super) fn corner_edges(slabs: &[PrismSlab], angle_tolerance: f64) -> Vec<CornerEdge> {
    let mut edges = Vec::new();
    for slab in slabs {
        for face in slab.faces() {
            push_ring_corners(
                &mut edges,
                &face.outer,
                slab.z_base(),
                slab.z_top(),
                angle_tolerance,
            );
            for hole in &face.holes {
                push_ring_corners(
                    &mut edges,
                    hole,
                    slab.z_base(),
                    slab.z_top(),
                    angle_tolerance,
                );
            }
        }
    }
    edges
}

/// Emits one vertical segment per corner of a single boundary ring.
fn push_ring_corners(
    edges: &mut Vec<CornerEdge>,
    ring: &Polygon,
    z_base: f64,
    z_top: f64,
    angle_tolerance: f64,
) {
    let count = ring.len();
    if count < 3 {
        return;
    }
    for index in 0..count {
        let previous = ring[(index + count - 1) % count];
        let here = ring[index];
        let next = ring[(index + 1) % count];

        let incoming = (here.0 - previous.0, here.1 - previous.1);
        let outgoing = (next.0 - here.0, next.1 - here.1);
        // A zero-length edge carries no direction — and produces no wall
        // quad either, so the vertex has no arris to draw.
        if incoming.0.hypot(incoming.1) < WALL_EPS || outgoing.0.hypot(outgoing.1) < WALL_EPS {
            continue;
        }

        if turn_angle(incoming, outgoing) > angle_tolerance {
            edges.push(CornerEdge {
                base: Point3::new(here.0, here.1, z_base),
                top: Point3::new(here.0, here.1, z_top),
            });
        }
    }
}

/// Unsigned angle in `0..=π` between two direction vectors.
///
/// `atan2(|cross|, dot)` rather than `acos(dot / |a| / |b|)`: it stays
/// accurate for the near-collinear directions that dominate a flattened
/// arc, where `acos` loses its significant digits.
fn turn_angle(a: (f64, f64), b: (f64, f64)) -> f64 {
    let cross = a.0 * b.1 - a.1 * b.0;
    let dot = a.0 * b.0 + a.1 * b.1;
    cross.abs().atan2(dot)
}
