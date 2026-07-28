//! Triangle-mesh assembly for the z-slab decomposition.
//!
//! Each slab contributes a vertical wall per boundary segment of its
//! region, plus the horizontal caps that close it where the neighbouring
//! slab's coverage stops. Caps are triangulated by the same constrained
//! Delaunay pass the planar-face tessellator uses, so a cap with holes is
//! handled by the shared odd-depth interior classification rather than a
//! second triangulator.
//!
//! # Orientation
//!
//! Slab regions arrive with the [`PolygonWithHoles`] winding contract —
//! CCW outer, CW holes — so for a boundary edge `(dx, dy)` the outward
//! normal (away from material) is `(dy, -dx)` on **both** ring kinds: on
//! an outer ring it points out of the body, on a hole ring it points into
//! the void. Triangle winding follows the same rule, and the caps flip
//! with their `+z` / `-z` normal.

use std::collections::HashMap;

use spade::{ConstrainedDelaunayTriangulation, Point2 as SpadePoint2, Triangulation};

use crate::error::Result;
use crate::math::{Point2, Point3, Vector3};
use crate::operations::boolean_2d::{subtract_all_with_holes, Polygon, PolygonWithHoles, WALL_EPS};
use crate::tessellation::tessellate_face::{classify_interior_faces, insert_constraint_loop};
use crate::tessellation::TriangleMesh;

use super::PrismSlab;

/// Builds the fused mesh from the contiguous slab list produced by
/// [`super::slab::slab_regions`].
///
/// # Errors
///
/// Propagates arrangement-engine failures from the cap region
/// differences, and tessellation failures from the cap triangulation.
pub(super) fn build_mesh(slabs: &[PrismSlab]) -> Result<TriangleMesh> {
    let mut mesh = TriangleMesh::default();
    for (index, slab) in slabs.iter().enumerate() {
        if slab.faces().is_empty() {
            continue;
        }

        for face in slab.faces() {
            push_wall(&mut mesh, &face.outer, slab.z_base(), slab.z_top());
            for hole in &face.holes {
                push_wall(&mut mesh, hole, slab.z_base(), slab.z_top());
            }
        }

        // Coverage below / above. Out-of-range neighbours are the empty
        // region, which closes the body at the extremes.
        let below = index
            .checked_sub(1)
            .and_then(|i| slabs.get(i))
            .map_or(&[][..], PrismSlab::faces);
        let above = slabs.get(index + 1).map_or(&[][..], PrismSlab::faces);

        for face in region_difference(slab.faces(), below)? {
            push_cap(&mut mesh, &face, slab.z_base(), false)?;
        }
        for face in region_difference(slab.faces(), above)? {
            push_cap(&mut mesh, &face, slab.z_top(), true)?;
        }
    }
    Ok(mesh)
}

/// `⋃ region ∩ ¬⋃ other` — the part of one slab's cross-section that the
/// neighbouring slab does not cover, i.e. the area needing a cap.
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

/// Extrudes one boundary ring into a band of vertical quads.
#[allow(clippy::cast_possible_truncation)]
fn push_wall(mesh: &mut TriangleMesh, ring: &Polygon, z_base: f64, z_top: f64) {
    let count = ring.len();
    if count < 3 {
        return;
    }
    let mut along = 0.0;
    for index in 0..count {
        let start = ring[index];
        let end = ring[(index + 1) % count];
        let (dx, dy) = (end.0 - start.0, end.1 - start.1);
        let length = dx.hypot(dy);
        if length < WALL_EPS {
            continue;
        }
        let normal = Vector3::new(dy / length, -dx / length, 0.0);

        let base = mesh.vertices.len() as u32;
        for (point, u) in [(start, along), (end, along + length)] {
            for z in [z_base, z_top] {
                mesh.vertices.push(Point3::new(point.0, point.1, z));
                mesh.normals.push(normal);
                mesh.uvs.push(Point2::new(u, z));
            }
        }
        // 0 = start/base, 1 = start/top, 2 = end/base, 3 = end/top.
        mesh.indices.push([base, base + 2, base + 3]);
        mesh.indices.push([base, base + 3, base + 1]);
        along += length;
    }
}

/// Triangulates one cap face at elevation `z`, facing up or down.
#[allow(clippy::cast_possible_truncation)]
fn push_cap(mesh: &mut TriangleMesh, face: &PolygonWithHoles, z: f64, up: bool) -> Result<()> {
    let mut cdt = ConstrainedDelaunayTriangulation::<SpadePoint2<f64>>::new();
    insert_constraint_loop(&mut cdt, &spade_ring(&face.outer))?;
    for hole in &face.holes {
        insert_constraint_loop(&mut cdt, &spade_ring(hole))?;
    }
    let interior = classify_interior_faces(&cdt);

    let normal = Vector3::new(0.0, 0.0, if up { 1.0 } else { -1.0 });
    let mut emitted: HashMap<usize, u32> = HashMap::new();
    for handle in cdt.inner_faces() {
        if !interior.contains(&handle.fix().index()) {
            continue;
        }
        let mut triangle = [0u32; 3];
        for (slot, vertex) in handle.vertices().iter().enumerate() {
            let key = vertex.fix().index();
            let index = if let Some(&existing) = emitted.get(&key) {
                existing
            } else {
                let position = vertex.position();
                let fresh = mesh.vertices.len() as u32;
                mesh.vertices.push(Point3::new(position.x, position.y, z));
                mesh.normals.push(normal);
                mesh.uvs.push(Point2::new(position.x, position.y));
                emitted.insert(key, fresh);
                fresh
            };
            triangle[slot] = index;
        }
        // CDT faces wind CCW in the xy projection, which is CCW seen from
        // +z; a down-facing cap needs the opposite winding.
        if !up {
            triangle.swap(1, 2);
        }
        mesh.indices.push(triangle);
    }
    Ok(())
}

fn spade_ring(ring: &Polygon) -> Vec<SpadePoint2<f64>> {
    ring.iter().map(|&(x, y)| SpadePoint2::new(x, y)).collect()
}
