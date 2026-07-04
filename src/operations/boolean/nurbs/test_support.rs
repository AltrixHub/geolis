//! Shared test helpers for the NURBS boolean acceptance tests.

use std::collections::HashMap;

use crate::math::Point3;
use crate::tessellation::{TessellateSolid, TessellationParams};
use crate::topology::{SolidId, TopologyStore};

/// Position-welded boundary-edge count of a solid's tessellation (edges
/// used != 2 after vertex welding).
///
/// Welding is proximity-based (each vertex joins an existing
/// representative within `1e-6`, probing the neighboring quantization
/// cells) rather than raw grid bucketing: adjacent faces emit rim
/// vertices agreeing to ~1e-9, and axis-aligned fixtures land those
/// exactly on `1e-6` grid-cell boundaries where single-cell bucketing
/// splits coincident points spuriously.
#[allow(clippy::unwrap_used)]
pub(crate) fn welded_boundary_edges(store: &TopologyStore, solid: SolidId) -> usize {
    const WELD: f64 = 1e-6;
    #[allow(clippy::cast_possible_truncation)]
    fn cell(p: &Point3) -> (i64, i64, i64) {
        (
            (p.x / WELD).round() as i64,
            (p.y / WELD).round() as i64,
            (p.z / WELD).round() as i64,
        )
    }
    #[allow(clippy::cast_possible_truncation)]
    fn canon_id(
        cells: &mut HashMap<(i64, i64, i64), Vec<u32>>,
        reps: &mut Vec<Point3>,
        p: &Point3,
    ) -> u32 {
        let (cx, cy, cz) = cell(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(ids) = cells.get(&(cx + dx, cy + dy, cz + dz)) {
                        for &id in ids {
                            if (reps[id as usize] - p).norm() <= WELD {
                                return id;
                            }
                        }
                    }
                }
            }
        }
        let id = reps.len() as u32;
        reps.push(*p);
        cells.entry((cx, cy, cz)).or_default().push(id);
        id
    }
    let mesh = TessellateSolid::new(solid, TessellationParams::default())
        .execute(store)
        .unwrap();
    assert!(!mesh.indices.is_empty(), "empty mesh");
    let mut cells: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
    let mut reps: Vec<Point3> = Vec::new();
    let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
    for tri in &mesh.indices {
        let a = canon_id(&mut cells, &mut reps, &mesh.vertices[tri[0] as usize]);
        let b = canon_id(&mut cells, &mut reps, &mesh.vertices[tri[1] as usize]);
        let c = canon_id(&mut cells, &mut reps, &mesh.vertices[tri[2] as usize]);
        for &(x, y) in &[(a, b), (b, c), (c, a)] {
            let key = if x < y { (x, y) } else { (y, x) };
            *counts.entry(key).or_insert(0) += 1;
        }
    }
    counts.values().filter(|&&c| c != 2).count()
}
