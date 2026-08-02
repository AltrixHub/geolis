//! Fixtures for [`super::UnionPrisms`].
//!
//! Every fixture runs [`assert_output_sane`], which pins the invariants the
//! module docs promise: finite geometry, in-range indices, elevations
//! drawn only from the breakpoint set, a positive divergence-theorem
//! volume, and corner edges that are vertical and slab-aligned.
//! Single-slab fixtures additionally pin combinatorial watertightness —
//! with one cross-section there are no T-junctions, so every welded edge
//! must be shared by exactly two triangles.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use std::collections::HashMap;

use super::{
    FusedPrisms, PrismCut, PrismProfile, PrismRegion, PrismSlab, UnionPrisms, DEFAULT_ARC_TOLERANCE,
};
use crate::geometry::pline::{Pline, PlineVertex};
use crate::operations::boolean_2d::{signed_area, PolygonWithHoles};
use crate::tessellation::TriangleMesh;

/// Absolute tolerance for area / volume assertions on line-only fixtures.
const EXACT_EPS: f64 = 1e-6;

// ===== Input builders =====

fn ring(points: &[(f64, f64)]) -> Pline {
    Pline {
        vertices: points
            .iter()
            .map(|&(x, y)| PlineVertex::line(x, y))
            .collect(),
        closed: true,
    }
}

fn region(points: &[(f64, f64)]) -> PrismRegion {
    PrismRegion::try_from_parts(ring(points), vec![]).expect("ring must be well-formed")
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> PrismRegion {
    region(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
}

/// A full circle as a two-vertex closed polyline (two semicircular arcs,
/// `bulge = tan(π / 4) = 1`), wound CCW.
fn circle(cx: f64, cy: f64, radius: f64) -> PrismRegion {
    let pline = Pline {
        vertices: vec![
            PlineVertex::new(cx - radius, cy, 1.0),
            PlineVertex::new(cx + radius, cy, 1.0),
        ],
        closed: true,
    };
    PrismRegion::try_from_parts(pline, vec![]).expect("circle must be well-formed")
}

fn prism(faces: Vec<PrismRegion>, z_base: f64, z_top: f64) -> PrismProfile {
    PrismProfile::try_new(faces, z_base, z_top).expect("z span must be valid")
}

fn cut(region: PrismRegion, z_base: f64, z_top: f64) -> PrismCut {
    PrismCut::try_new(region, z_base, z_top).expect("cut z span must be valid")
}

fn run(profiles: Vec<PrismProfile>) -> FusedPrisms {
    UnionPrisms::new(profiles)
        .execute()
        .expect("union_prisms must succeed")
}

// ===== Region readers =====

fn face_area(face: &PolygonWithHoles) -> f64 {
    signed_area(&face.outer) + face.holes.iter().map(signed_area).sum::<f64>()
}

fn slab_area(faces: &[PolygonWithHoles]) -> f64 {
    faces.iter().map(face_area).sum()
}

/// Total number of boundary segments across every ring of a slab — the
/// number of vertical wall quads the slab must produce.
fn slab_segments(faces: &[PolygonWithHoles]) -> usize {
    faces
        .iter()
        .map(|f| f.outer.len() + f.holes.iter().map(Vec::len).sum::<usize>())
        .sum()
}

/// Number of genuine corners of a ring: the arrangement engine keeps the
/// vertex where two inputs' collinear edges met, so a straight run may
/// carry interior vertices.
fn corner_count(ring: &[(f64, f64)]) -> usize {
    let n = ring.len();
    (0..n)
        .filter(|&i| {
            let (prev, here, next) = (ring[(i + n - 1) % n], ring[i], ring[(i + 1) % n]);
            let (ax, ay) = (here.0 - prev.0, here.1 - prev.1);
            let (bx, by) = (next.0 - here.0, next.1 - here.1);
            (ax * by - ay * bx).abs() > EXACT_EPS
        })
        .count()
}

// ===== Mesh readers =====

fn tri_points(mesh: &TriangleMesh, tri: [u32; 3]) -> [crate::math::Point3; 3] {
    [
        mesh.vertices[tri[0] as usize],
        mesh.vertices[tri[1] as usize],
        mesh.vertices[tri[2] as usize],
    ]
}

/// Divergence-theorem volume of the (geometrically closed) surface.
///
/// Exact even where the mesh T-junctions, because a crack between a cap
/// and a longer wall edge has zero area.
fn signed_volume(mesh: &TriangleMesh) -> f64 {
    mesh.indices
        .iter()
        .map(|&tri| {
            let [a, b, c] = tri_points(mesh, tri);
            a.coords.dot(&b.coords.cross(&c.coords)) / 6.0
        })
        .sum()
}

/// Number of vertical wall quads: side triangles span two elevations and
/// come in pairs.
fn wall_quad_count(mesh: &TriangleMesh) -> usize {
    let side_triangles = mesh
        .indices
        .iter()
        .filter(|&&tri| {
            let pts = tri_points(mesh, tri);
            let min = pts.iter().map(|p| p.z).fold(f64::INFINITY, f64::min);
            let max = pts.iter().map(|p| p.z).fold(f64::NEG_INFINITY, f64::max);
            max - min > EXACT_EPS
        })
        .count();
    assert_eq!(side_triangles % 2, 0, "wall triangles must come in pairs");
    side_triangles / 2
}

/// Area of the horizontal caps at elevation `z` facing up (`up = true`)
/// or down.
fn cap_area(mesh: &TriangleMesh, z: f64, up: bool) -> f64 {
    mesh.indices
        .iter()
        .filter_map(|&tri| {
            let pts = tri_points(mesh, tri);
            if pts.iter().any(|p| (p.z - z).abs() > EXACT_EPS) {
                return None;
            }
            let cross = (pts[1] - pts[0]).cross(&(pts[2] - pts[0]));
            if (cross.z > 0.0) != up {
                return None;
            }
            Some(cross.norm() / 2.0)
        })
        .sum()
}

/// Corner edges as `(x, y, z_base, z_top)`, in emission order.
fn corner_spans(out: &FusedPrisms) -> Vec<(f64, f64, f64, f64)> {
    out.corner_edges()
        .iter()
        .map(|edge| (edge.base().x, edge.base().y, edge.base().z, edge.top().z))
        .collect()
}

/// The corner edges standing at plan position `(x, y)`, as `(z_base, z_top)`.
fn corners_at(out: &FusedPrisms, x: f64, y: f64) -> Vec<(f64, f64)> {
    corners_near(out, x, y, EXACT_EPS)
}

/// [`corners_at`] with an explicit plan tolerance, for positions that a
/// flattened arc only reaches within its chord deviation.
fn corners_near(out: &FusedPrisms, x: f64, y: f64, tolerance: f64) -> Vec<(f64, f64)> {
    corner_spans(out)
        .into_iter()
        .filter(|&(ex, ey, _, _)| (ex - x).abs() < tolerance && (ey - y).abs() < tolerance)
        .map(|(_, _, z0, z1)| (z0, z1))
        .collect()
}

/// Asserts each corner edge is vertical and spans exactly one slab.
fn assert_corner_edges_sane(out: &FusedPrisms) {
    for edge in out.corner_edges() {
        let (base, top) = (edge.base(), edge.top());
        assert!(
            (base.x - top.x).abs() < EXACT_EPS && (base.y - top.y).abs() < EXACT_EPS,
            "corner edge {base:?} → {top:?} is not vertical"
        );
        assert!(base.x.is_finite() && base.y.is_finite());
        assert!(top.z > base.z, "corner edge must rise");
        assert!(
            out.slabs().iter().any(|slab| {
                (slab.z_base() - base.z).abs() < EXACT_EPS
                    && (slab.z_top() - top.z).abs() < EXACT_EPS
            }),
            "corner edge {base:?} → {top:?} spans no slab"
        );
    }
}

/// Every invariant the module docs promise for an output.
fn assert_output_sane(out: &FusedPrisms) {
    assert_corner_edges_sane(out);
    let mesh = out.mesh();
    assert!(!mesh.indices.is_empty(), "mesh must carry triangles");
    assert_eq!(mesh.vertices.len(), mesh.normals.len());
    assert_eq!(mesh.vertices.len(), mesh.uvs.len());

    let mut breakpoints: Vec<f64> = out
        .slabs()
        .iter()
        .flat_map(|slab| [slab.z_base(), slab.z_top()])
        .collect();
    breakpoints.sort_by(f64::total_cmp);
    breakpoints.dedup_by(|a, b| (*a - *b).abs() <= EXACT_EPS);

    for vertex in &mesh.vertices {
        assert!(
            vertex.x.is_finite() && vertex.y.is_finite() && vertex.z.is_finite(),
            "non-finite vertex {vertex:?}"
        );
        assert!(
            breakpoints
                .iter()
                .any(|&z| (vertex.z - z).abs() <= EXACT_EPS),
            "vertex z {} is not a breakpoint of {breakpoints:?}",
            vertex.z
        );
    }
    for normal in &mesh.normals {
        assert!(
            (normal.norm() - 1.0).abs() < 1e-9,
            "normal {normal:?} is not unit-length"
        );
    }
    for tri in &mesh.indices {
        for &index in tri {
            assert!(
                (index as usize) < mesh.vertices.len(),
                "index {index} out of range"
            );
        }
    }
    assert!(
        signed_volume(mesh) > 0.0,
        "outward orientation must give a positive volume"
    );
}

/// Quantizes positions and asserts every undirected triangle edge is
/// shared by exactly two triangles. Valid only for single-slab bodies —
/// see the module docs on T-junctions between slabs.
#[allow(clippy::cast_possible_truncation)]
fn assert_position_weld_watertight(mesh: &TriangleMesh) {
    let quantize = |value: f64| (value * 1e6).round() as i64;
    let mut ids: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut vertex_ids = Vec::with_capacity(mesh.vertices.len());
    for p in &mesh.vertices {
        let next = ids.len();
        let id = *ids
            .entry((quantize(p.x), quantize(p.y), quantize(p.z)))
            .or_insert(next);
        vertex_ids.push(id);
    }
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in &mesh.indices {
        for k in 0..3 {
            let a = vertex_ids[tri[k] as usize];
            let b = vertex_ids[tri[(k + 1) % 3] as usize];
            if a == b {
                continue;
            }
            *edge_count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    let boundary = edge_count.values().filter(|&&count| count != 2).count();
    assert_eq!(
        boundary, 0,
        "mesh has {boundary} non-manifold/boundary edges"
    );
}

// ===== Fixtures =====

/// The horizontal arm of the crossing fixtures: 10 × 2, centred.
fn bar_x() -> PrismRegion {
    rect(-5.0, -1.0, 5.0, 1.0)
}

/// The vertical arm of the crossing fixtures: 2 × 10, centred.
fn bar_y() -> PrismRegion {
    rect(-1.0, -5.0, 1.0, 5.0)
}

/// Area of `bar_x ∪ bar_y`: 20 + 20 − 4.
const PLUS_AREA: f64 = 36.0;

#[test]
fn crossing_prisms_fuse_without_internal_side_faces() {
    let out = run(vec![
        prism(vec![bar_x()], 0.0, 3.0),
        prism(vec![bar_y()], 0.0, 3.0),
    ]);
    assert_output_sane(&out);

    assert_eq!(out.slabs().len(), 1, "equal spans give a single slab");
    let faces = out.slabs()[0].faces();
    assert_eq!(faces.len(), 1, "the two prisms fuse into one body");
    assert!(faces[0].holes.is_empty());
    assert_eq!(faces[0].outer.len(), 12, "the plus outline has 12 corners");
    assert!((slab_area(faces) - PLUS_AREA).abs() < EXACT_EPS);

    // The junction interior carries no side face: the wall count follows
    // the FUSED outline (12), not the sum of the inputs (4 + 4 = 8) and
    // not the sum plus interior seams.
    assert_eq!(wall_quad_count(out.mesh()), slab_segments(faces));
    assert_eq!(wall_quad_count(out.mesh()), 12);

    assert!((cap_area(out.mesh(), 0.0, false) - PLUS_AREA).abs() < EXACT_EPS);
    assert!((cap_area(out.mesh(), 3.0, true) - PLUS_AREA).abs() < EXACT_EPS);
    assert!((signed_volume(out.mesh()) - PLUS_AREA * 3.0).abs() < EXACT_EPS);
    assert_position_weld_watertight(out.mesh());
}

#[test]
fn different_heights_step_caps_the_lower_slab() {
    let out = run(vec![
        prism(vec![bar_x()], 0.0, 3.0),
        prism(vec![bar_y()], 0.0, 6.0),
    ]);
    assert_output_sane(&out);

    assert_eq!(out.slabs().len(), 2, "the step splits the body in two");
    let lower = &out.slabs()[0];
    let upper = &out.slabs()[1];
    assert!((lower.z_base() - 0.0).abs() < EXACT_EPS);
    assert!((lower.z_top() - 3.0).abs() < EXACT_EPS);
    assert!((upper.z_base() - 3.0).abs() < EXACT_EPS);
    assert!((upper.z_top() - 6.0).abs() < EXACT_EPS);

    // Lower interval fused, upper interval = the tall prism alone.
    assert_eq!(lower.faces().len(), 1);
    assert_eq!(lower.faces()[0].outer.len(), 12);
    assert!((slab_area(lower.faces()) - PLUS_AREA).abs() < EXACT_EPS);
    assert_eq!(upper.faces().len(), 1);
    assert_eq!(upper.faces()[0].outer.len(), 4);
    assert!((slab_area(upper.faces()) - 20.0).abs() < EXACT_EPS);

    // The step is a top cap on the lower slab covering exactly the area
    // the tall prism does not continue over. Nothing faces down there:
    // the upper cross-section is a subset of the lower one.
    assert!((cap_area(out.mesh(), 3.0, true) - (PLUS_AREA - 20.0)).abs() < EXACT_EPS);
    assert!(cap_area(out.mesh(), 3.0, false) < EXACT_EPS);
    assert!((cap_area(out.mesh(), 0.0, false) - PLUS_AREA).abs() < EXACT_EPS);
    assert!((cap_area(out.mesh(), 6.0, true) - 20.0).abs() < EXACT_EPS);

    let expected = PLUS_AREA * 3.0 + 20.0 * 3.0;
    assert!((signed_volume(out.mesh()) - expected).abs() < EXACT_EPS);
}

#[test]
fn mid_height_cut_pierces_only_its_own_interval() {
    let out = run(vec![prism(vec![rect(0.0, 0.0, 10.0, 10.0)], 0.0, 3.0)
        .with_cut(cut(rect(4.0, 4.0, 6.0, 6.0), 1.0, 2.0))]);
    assert_output_sane(&out);

    assert_eq!(out.slabs().len(), 3, "the cut opens two extra breakpoints");
    let [below, pierced, above] = [&out.slabs()[0], &out.slabs()[1], &out.slabs()[2]];

    // Only the interval containing the cut carries the hole.
    assert_eq!(below.faces().len(), 1);
    assert!(below.faces()[0].holes.is_empty(), "no full-height hole");
    assert!((slab_area(below.faces()) - 100.0).abs() < EXACT_EPS);

    assert_eq!(pierced.faces().len(), 1);
    assert_eq!(pierced.faces()[0].holes.len(), 1, "the window is a hole");
    assert!(
        signed_area(&pierced.faces()[0].holes[0]) < 0.0,
        "hole is CW"
    );
    assert!((slab_area(pierced.faces()) - 96.0).abs() < EXACT_EPS);

    assert_eq!(above.faces().len(), 1);
    assert!(above.faces()[0].holes.is_empty(), "no full-height hole");
    assert!((slab_area(above.faces()) - 100.0).abs() < EXACT_EPS);

    // The opening is closed by a soffit above and a sill below.
    assert!((cap_area(out.mesh(), 1.0, true) - 4.0).abs() < EXACT_EPS);
    assert!((cap_area(out.mesh(), 2.0, false) - 4.0).abs() < EXACT_EPS);
    assert!(cap_area(out.mesh(), 1.0, false) < EXACT_EPS);
    assert!(cap_area(out.mesh(), 2.0, true) < EXACT_EPS);

    assert!((signed_volume(out.mesh()) - (100.0 * 3.0 - 4.0)).abs() < EXACT_EPS);
}

#[test]
fn cut_spanning_a_fused_junction_removes_the_neighbour_material() {
    // The cut belongs to the horizontal arm but reaches past it in y, so
    // the material it removes is contributed by the vertical arm. Only
    // the fused region makes it an enclosed hole.
    let out = run(vec![
        prism(vec![bar_x()], 0.0, 3.0).with_cut(cut(rect(-0.5, -2.0, 0.5, 2.0), 1.0, 2.0)),
        prism(vec![bar_y()], 0.0, 3.0),
    ]);
    assert_output_sane(&out);

    assert_eq!(out.slabs().len(), 3);
    let pierced = &out.slabs()[1];
    assert!((pierced.z_base() - 1.0).abs() < EXACT_EPS);
    assert!((pierced.z_top() - 2.0).abs() < EXACT_EPS);
    assert_eq!(pierced.faces().len(), 1);
    assert_eq!(pierced.faces()[0].holes.len(), 1);
    assert_eq!(pierced.faces()[0].outer.len(), 12, "outline still the plus");
    assert!((slab_area(pierced.faces()) - (PLUS_AREA - 4.0)).abs() < EXACT_EPS);

    for slab in [&out.slabs()[0], &out.slabs()[2]] {
        assert_eq!(slab.faces().len(), 1);
        assert!(slab.faces()[0].holes.is_empty());
        assert!((slab_area(slab.faces()) - PLUS_AREA).abs() < EXACT_EPS);
    }

    assert!((signed_volume(out.mesh()) - (PLUS_AREA * 3.0 - 4.0)).abs() < EXACT_EPS);
}

#[test]
fn curved_profile_fuses_with_a_straight_one() {
    let out = UnionPrisms::new(vec![
        prism(vec![circle(0.0, 0.0, 2.0)], 0.0, 3.0),
        prism(vec![rect(-6.0, -0.5, 6.0, 0.5)], 0.0, 3.0),
    ])
    .with_arc_tolerance(0.001)
    .execute()
    .expect("curved union must succeed");
    assert_output_sane(&out);

    assert_eq!(out.slabs().len(), 1);
    let faces = out.slabs()[0].faces();
    assert_eq!(faces.len(), 1, "the bar pierces the disc: one fused face");
    assert!(faces[0].holes.is_empty());
    assert!(
        faces[0].outer.len() > 40,
        "the arc must survive as many chords, got {}",
        faces[0].outer.len()
    );

    // π·2² + 12·1 − ∫_{-0.5}^{0.5} 2·√(4 − y²) dy.
    let overlap = 2.0 * (0.5 * 3.75_f64.sqrt() + 4.0 * 0.25_f64.asin());
    let expected = std::f64::consts::PI * 4.0 + 12.0 - overlap;
    assert!(
        (slab_area(faces) - expected).abs() < 0.02,
        "fused area {} vs expected {expected}",
        slab_area(faces)
    );
    assert!((signed_volume(out.mesh()) - slab_area(faces) * 3.0).abs() < 1e-6);
    assert_eq!(wall_quad_count(out.mesh()), slab_segments(faces));
    assert_position_weld_watertight(out.mesh());
}

#[test]
fn two_l_shapes_enclose_a_courtyard() {
    // Bottom-left L and its point reflection through (5, 5): together a
    // 10 × 10 square ring around a 6 × 6 courtyard.
    let lower = region(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 2.0),
        (2.0, 2.0),
        (2.0, 10.0),
        (0.0, 10.0),
    ]);
    let upper = region(&[
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 8.0),
        (8.0, 8.0),
        (8.0, 0.0),
        (10.0, 0.0),
    ]);
    let out = run(vec![
        prism(vec![lower], 0.0, 3.0),
        prism(vec![upper], 0.0, 3.0),
    ]);
    assert_output_sane(&out);

    assert_eq!(out.slabs().len(), 1);
    let faces = out.slabs()[0].faces();
    assert_eq!(faces.len(), 1);
    assert_eq!(
        corner_count(&faces[0].outer),
        4,
        "the outer boundary is the square"
    );
    assert_eq!(faces[0].holes.len(), 1, "the courtyard is a hole ring");
    assert_eq!(corner_count(&faces[0].holes[0]), 4);
    assert!((signed_area(&faces[0].holes[0]) + 36.0).abs() < EXACT_EPS);
    assert!((slab_area(faces) - 64.0).abs() < EXACT_EPS);

    // Outer walls plus the inner walls that face the courtyard — one quad
    // per boundary segment, nothing more.
    assert_eq!(wall_quad_count(out.mesh()), slab_segments(faces));

    // The inner wall faces INTO the void: its normals point away from the
    // material, i.e. toward the courtyard centre.
    let centre = crate::math::Point3::new(5.0, 5.0, 1.5);
    let inward = out
        .mesh()
        .indices
        .iter()
        .filter(|&&tri| {
            let pts = tri_points(out.mesh(), tri);
            let mid = (pts[0].coords + pts[1].coords + pts[2].coords) / 3.0;
            if (pts[0].z - pts[1].z).abs() < EXACT_EPS && (pts[1].z - pts[2].z).abs() < EXACT_EPS {
                return false;
            }
            out.mesh().normals[tri[0] as usize].dot(&(centre.coords - mid)) > 0.0
        })
        .count();
    assert_eq!(
        inward / 2,
        corner_count(&faces[0].holes[0]),
        "the courtyard is walled on all four sides"
    );

    assert!((signed_volume(out.mesh()) - 64.0 * 3.0).abs() < EXACT_EPS);
    assert_position_weld_watertight(out.mesh());
}

// ===== Degenerate input is skipped =====

#[test]
fn degenerate_rings_are_skipped_not_rejected() {
    // A zero-area sliver ring and a profile with no faces at all.
    let sliver = region(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]);
    let out = run(vec![
        prism(vec![sliver], 0.0, 3.0),
        prism(vec![], 0.0, 9.0),
        prism(vec![rect(0.0, 0.0, 2.0, 2.0)], 0.0, 3.0),
    ]);
    assert_output_sane(&out);

    assert_eq!(out.slabs().len(), 1, "only the real prism opens a slab");
    assert!((slab_area(out.slabs()[0].faces()) - 4.0).abs() < EXACT_EPS);
    assert!((signed_volume(out.mesh()) - 12.0).abs() < EXACT_EPS);
    assert_eq!(
        out.corner_edges().len(),
        4,
        "the skipped rings contribute no arris"
    );
}

#[test]
fn a_degenerate_hole_drops_only_that_hole() {
    let outer = ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
    let real_hole = ring(&[(2.0, 2.0), (2.0, 4.0), (4.0, 4.0), (4.0, 2.0)]);
    let flat_hole = ring(&[(6.0, 6.0), (7.0, 6.0), (8.0, 6.0)]);
    let region =
        PrismRegion::try_from_parts(outer, vec![real_hole, flat_hole]).expect("rings are legal");

    let out = run(vec![prism(vec![region], 0.0, 2.0)]);
    assert_output_sane(&out);
    assert_eq!(out.slabs()[0].faces()[0].holes.len(), 1);
    assert!((slab_area(out.slabs()[0].faces()) - 96.0).abs() < EXACT_EPS);
}

#[test]
fn no_profiles_yields_an_empty_body() {
    let out = UnionPrisms::new(vec![])
        .execute()
        .expect("empty input is Ok");
    assert!(out.slabs().is_empty());
    assert!(out.mesh().indices.is_empty());
    assert!(out.mesh().vertices.is_empty());
}

#[test]
fn a_cut_outside_its_profile_changes_nothing() {
    let plain = run(vec![prism(vec![rect(0.0, 0.0, 4.0, 4.0)], 0.0, 3.0)]);
    let with_cut = run(vec![prism(vec![rect(0.0, 0.0, 4.0, 4.0)], 0.0, 3.0)
        .with_cut(cut(rect(1.0, 1.0, 2.0, 2.0), 5.0, 7.0))]);
    assert_output_sane(&with_cut);

    assert_eq!(with_cut.slabs().len(), plain.slabs().len());
    assert!((signed_volume(with_cut.mesh()) - signed_volume(plain.mesh())).abs() < EXACT_EPS);
}

// ===== Contract violations are errors =====

#[test]
fn non_positive_profile_span_is_rejected() {
    let err = PrismProfile::try_new(vec![rect(0.0, 0.0, 1.0, 1.0)], 2.0, 2.0)
        .expect_err("z_top == z_base must be rejected");
    let message = format!("{err}");
    assert!(message.contains("z_top must exceed z_base"), "{message}");

    let err = PrismProfile::try_new(vec![rect(0.0, 0.0, 1.0, 1.0)], 2.0, 1.0)
        .expect_err("inverted span must be rejected");
    assert!(format!("{err}").contains("z_top must exceed z_base"));
}

#[test]
fn non_positive_cut_span_is_rejected() {
    let err = PrismCut::try_new(rect(0.0, 0.0, 1.0, 1.0), 1.0, 1.0)
        .expect_err("zero-height cut must be rejected");
    assert!(format!("{err}").contains("PrismCut::try_new"));
}

#[test]
fn non_finite_z_is_rejected() {
    let err = PrismProfile::try_new(vec![], 0.0, f64::NAN).expect_err("NaN must be rejected");
    assert!(format!("{err}").contains("finite"));
}

#[test]
fn an_open_ring_is_rejected() {
    let open = Pline {
        vertices: vec![
            PlineVertex::line(0.0, 0.0),
            PlineVertex::line(1.0, 0.0),
            PlineVertex::line(1.0, 1.0),
        ],
        closed: false,
    };
    let err = PrismRegion::try_from_parts(open, vec![]).expect_err("open rings must be rejected");
    assert!(format!("{err}").contains("closed = true"));
}

#[test]
fn a_non_finite_vertex_is_rejected() {
    let broken = Pline {
        vertices: vec![
            PlineVertex::line(0.0, 0.0),
            PlineVertex::line(f64::INFINITY, 0.0),
            PlineVertex::line(1.0, 1.0),
        ],
        closed: true,
    };
    let err = PrismRegion::try_from_parts(broken, vec![]).expect_err("must be rejected");
    assert!(format!("{err}").contains("not finite"));
}

#[test]
fn a_non_positive_arc_tolerance_is_rejected() {
    let err = UnionPrisms::new(vec![prism(vec![rect(0.0, 0.0, 1.0, 1.0)], 0.0, 1.0)])
        .with_arc_tolerance(0.0)
        .execute()
        .expect_err("zero tolerance must be rejected");
    assert!(format!("{err}").contains("arc tolerance"));
}

// ===== Vertical corner edges =====

#[test]
fn x_cross_emits_one_corner_edge_per_fused_outline_corner() {
    let out = run(vec![
        prism(vec![bar_x()], 0.0, 3.0),
        prism(vec![bar_y()], 0.0, 3.0),
    ]);
    assert_output_sane(&out);

    // Every vertex of the plus outline is a 90° / 270° turn, so all 12
    // become arrises — and nothing else does.
    let outline = &out.slabs()[0].faces()[0].outer;
    assert_eq!(outline.len(), 12);
    assert_eq!(out.corner_edges().len(), 12);

    let mut expected: Vec<(f64, f64)> = outline.clone();
    let mut actual: Vec<(f64, f64)> = corner_spans(&out)
        .into_iter()
        .map(|(x, y, z0, z1)| {
            assert!((z0 - 0.0).abs() < EXACT_EPS && (z1 - 3.0).abs() < EXACT_EPS);
            (x, y)
        })
        .collect();
    let sort = |v: &mut Vec<(f64, f64)>| {
        v.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    };
    sort(&mut expected);
    sort(&mut actual);
    assert_eq!(actual, expected, "corner edges sit on the outline corners");
}

#[test]
fn a_step_stacks_the_tall_prisms_corner_edges() {
    let out = run(vec![
        prism(vec![bar_x()], 0.0, 3.0),
        prism(vec![bar_y()], 0.0, 6.0),
    ]);
    assert_output_sane(&out);

    // 12 corners on the fused lower slab + 4 on the tall prism alone.
    assert_eq!(out.corner_edges().len(), 16);

    // The tall prism's own corners survive into both slabs, stacked at
    // the step elevation rather than merged.
    for (x, y) in [(1.0, 5.0), (-1.0, 5.0), (1.0, -5.0), (-1.0, -5.0)] {
        let mut spans = corners_at(&out, x, y);
        spans.sort_by(|a, b| a.0.total_cmp(&b.0));
        assert_eq!(
            spans,
            vec![(0.0, 3.0), (3.0, 6.0)],
            "corner ({x}, {y}) must stack across the step"
        );
    }

    // A corner that only the short prism contributes stops at the step.
    assert_eq!(corners_at(&out, 5.0, 1.0), vec![(0.0, 3.0)]);
}

#[test]
fn a_default_tessellated_circle_has_no_corner_edges() {
    let out = run(vec![prism(vec![circle(0.0, 0.0, 2.0)], 0.0, 3.0)]);
    assert_output_sane(&out);

    assert_eq!(out.slabs().len(), 1);
    assert!(
        out.slabs()[0].faces()[0].outer.len() > 8,
        "the circle must actually be faceted"
    );
    assert!(
        out.corner_edges().is_empty(),
        "arc facets are not corners: got {:?}",
        corner_spans(&out)
    );
}

/// Pins the [`super::DEFAULT_CORNER_ANGLE_TOLERANCE`] derivation at its
/// worst case: the smallest radius the default claims to handle,
/// `10 × DEFAULT_ARC_TOLERANCE`, where the sagitta criterion allows the
/// widest chord angle. Still zero corner edges.
#[test]
fn the_worst_case_arc_facet_clears_the_default_corner_tolerance() {
    let radius = 10.0 * DEFAULT_ARC_TOLERANCE;
    let out = run(vec![prism(vec![circle(0.0, 0.0, radius)], 0.0, 1.0)]);
    assert_output_sane(&out);
    assert!(
        out.corner_edges().is_empty(),
        "r = 10 × arc tolerance must still tessellate below the threshold: {:?}",
        corner_spans(&out)
    );
}

#[test]
fn circle_and_bar_emit_corners_only_where_the_boundary_really_turns() {
    let out = UnionPrisms::new(vec![
        prism(vec![circle(0.0, 0.0, 2.0)], 0.0, 3.0),
        prism(vec![rect(-6.0, -0.5, 6.0, 0.5)], 0.0, 3.0),
    ])
    .with_arc_tolerance(0.001)
    .execute()
    .expect("curved union must succeed");
    assert_output_sane(&out);

    // 4 bar ends + the 4 points where the bar's flanks cross the circle.
    // Nothing on the arc itself.
    assert_eq!(out.corner_edges().len(), 8);

    for (x, y) in [(-6.0, -0.5), (-6.0, 0.5), (6.0, -0.5), (6.0, 0.5)] {
        assert_eq!(
            corners_at(&out, x, y).len(),
            1,
            "expected exactly one corner edge at the bar end ({x}, {y}); got {:?}",
            corner_spans(&out)
        );
    }
    // The crossings sit on a CHORD of the flattened circle, so they land
    // within the arc tolerance of the analytic x = √(2² − 0.5²), not on it.
    let crossing = 3.75_f64.sqrt();
    for (x, y) in [
        (-crossing, -0.5),
        (-crossing, 0.5),
        (crossing, -0.5),
        (crossing, 0.5),
    ] {
        assert_eq!(
            corners_near(&out, x, y, 0.01).len(),
            1,
            "expected exactly one corner edge at the crossing ({x}, {y}); got {:?}",
            corner_spans(&out)
        );
    }
}

#[test]
fn courtyard_corners_emit_and_collinear_seams_do_not() {
    let lower = region(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 2.0),
        (2.0, 2.0),
        (2.0, 10.0),
        (0.0, 10.0),
    ]);
    let upper = region(&[
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 8.0),
        (8.0, 8.0),
        (8.0, 0.0),
        (10.0, 0.0),
    ]);
    let out = run(vec![
        prism(vec![lower], 0.0, 3.0),
        prism(vec![upper], 0.0, 3.0),
    ]);
    assert_output_sane(&out);

    // The outer ring carries 8 vertices but only 4 turns: the two extra
    // are the collinear seams where the L-shapes' edges met.
    let face = &out.slabs()[0].faces()[0];
    assert_eq!(face.outer.len(), 8);
    assert_eq!(out.corner_edges().len(), 8, "4 outer + 4 courtyard corners");

    for (x, y) in [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)] {
        assert_eq!(corners_at(&out, x, y).len(), 1, "outer corner ({x}, {y})");
    }
    for (x, y) in [(2.0, 2.0), (8.0, 2.0), (8.0, 8.0), (2.0, 8.0)] {
        assert_eq!(corners_at(&out, x, y).len(), 1, "hole corner ({x}, {y})");
    }
    // The seam vertices themselves — mid-edge, zero turn — emit nothing.
    for (x, y) in [(8.0, 0.0), (0.0, 8.0), (10.0, 2.0), (2.0, 10.0)] {
        assert!(
            corners_at(&out, x, y).is_empty(),
            "collinear seam ({x}, {y}) must not emit an edge"
        );
    }
}

#[test]
fn an_opening_emits_corner_edges_only_in_its_own_slab() {
    let out = run(vec![prism(vec![rect(0.0, 0.0, 10.0, 10.0)], 0.0, 3.0)
        .with_cut(cut(rect(4.0, 4.0, 6.0, 6.0), 1.0, 2.0))]);
    assert_output_sane(&out);

    // 4 outer corners per slab (3 slabs) + the opening's 4, in the
    // pierced slab only.
    assert_eq!(out.corner_edges().len(), 16);
    for (x, y) in [(4.0, 4.0), (6.0, 4.0), (6.0, 6.0), (4.0, 6.0)] {
        assert_eq!(
            corners_at(&out, x, y),
            vec![(1.0, 2.0)],
            "the reveal at ({x}, {y}) spans only the opening"
        );
    }
    assert_eq!(
        corners_at(&out, 0.0, 0.0),
        vec![(0.0, 1.0), (1.0, 2.0), (2.0, 3.0)],
        "the outer corner is stacked once per slab"
    );
}

#[test]
fn an_empty_body_has_no_corner_edges() {
    let out = UnionPrisms::new(vec![])
        .execute()
        .expect("empty input is Ok");
    assert!(out.corner_edges().is_empty());
}

#[test]
fn a_shallower_corner_tolerance_keeps_shallower_turns() {
    // A chevron band whose apexes turn by 0.583 rad (~33°): below the
    // default threshold, above an explicitly lowered one.
    let chevron = region(&[
        (0.0, 0.0),
        (4.0, 1.2),
        (8.0, 0.0),
        (8.0, 1.0),
        (4.0, 2.2),
        (0.0, 1.0),
    ]);
    let profiles = || vec![prism(vec![chevron.clone()], 0.0, 2.0)];

    let coarse = UnionPrisms::new(profiles())
        .execute()
        .expect("must succeed");
    assert_output_sane(&coarse);
    assert_eq!(coarse.corner_edges().len(), 4, "only the four band ends");
    assert!(corners_at(&coarse, 4.0, 1.2).is_empty());

    let fine = UnionPrisms::new(profiles())
        .with_corner_angle_tolerance(0.5)
        .execute()
        .expect("must succeed");
    assert_output_sane(&fine);
    assert_eq!(fine.corner_edges().len(), 6, "the two apexes join in");
    assert_eq!(corners_at(&fine, 4.0, 1.2).len(), 1);
    assert_eq!(corners_at(&fine, 4.0, 2.2).len(), 1);
}

#[test]
fn a_non_positive_corner_angle_tolerance_is_rejected() {
    for bad in [0.0, -0.1, f64::NAN, std::f64::consts::PI, 4.0] {
        let err = UnionPrisms::new(vec![prism(vec![rect(0.0, 0.0, 1.0, 1.0)], 0.0, 1.0)])
            .with_corner_angle_tolerance(bad)
            .execute()
            .expect_err("must be rejected");
        assert!(
            format!("{err}").contains("corner angle tolerance"),
            "{bad}: {err}"
        );
    }
}

#[test]
fn union_prisms_is_deterministic() {
    let build = || {
        vec![
            prism(vec![bar_x()], 0.0, 3.0).with_cut(cut(rect(-0.5, -2.0, 0.5, 2.0), 1.0, 2.0)),
            prism(vec![bar_y()], 0.0, 6.0),
        ]
    };
    let first = run(build());
    let second = run(build());

    assert_eq!(first.slabs().len(), second.slabs().len());
    for (a, b) in first.slabs().iter().zip(second.slabs()) {
        assert_eq!(a.faces().len(), b.faces().len());
        for (fa, fb) in a.faces().iter().zip(b.faces()) {
            assert_eq!(fa, fb, "slab topology must be stable");
        }
    }
    assert_eq!(first.mesh().indices, second.mesh().indices);
    assert_eq!(first.mesh().vertices.len(), second.mesh().vertices.len());
    assert_eq!(first.corner_edges(), second.corner_edges());
}

/// The material a slab fuses is decided by its cover — which profiles
/// reach it — so [`super::slab::slab_regions`] computes each distinct
/// cover's union once and reuses it across every slab that shares it.
///
/// This fixture is what stops that reuse from being keyed wrongly:
/// staggered z spans make EVERY slab a different cover, so a memo that
/// confused two covers would hand a slab the wrong material. Each slab's
/// area is checked against a union of exactly the profiles that reach it,
/// computed independently of the memo.
#[test]
fn each_slab_fuses_exactly_the_profiles_that_reach_it() {
    // Three overlapping bands with staggered spans: 0–3, 1–4, 2–5. The
    // breakpoints 0,1,2,3,4,5 open five slabs whose covers are
    // {0}, {0,1}, {0,1,2}, {1,2}, {2} — every one distinct.
    let spans = [(0.0, 3.0), (1.0, 4.0), (2.0, 5.0)];
    let footprints = [
        rect(0.0, 0.0, 6.0, 1.0),
        rect(2.0, -2.0, 3.0, 4.0),
        rect(4.0, 0.5, 9.0, 1.5),
    ];
    let profiles: Vec<PrismProfile> = spans
        .iter()
        .zip(footprints.iter())
        .map(|(&(z0, z1), f)| prism(vec![f.clone()], z0, z1))
        .collect();

    let out = run(profiles);
    assert_eq!(out.slabs().len(), 5, "staggered spans open five slabs");

    for slab in out.slabs() {
        let mid = 0.5 * (slab.z_base() + slab.z_top());
        // The cover, recomputed here rather than read out of the memo.
        let cover: Vec<PrismProfile> = spans
            .iter()
            .zip(footprints.iter())
            .filter(|(&(z0, z1), _)| mid > z0 && mid < z1)
            .map(|(_, f)| prism(vec![f.clone()], 0.0, 1.0))
            .collect();
        assert!(!cover.is_empty(), "every reported slab has material");
        let expected: f64 = run(cover)
            .slabs()
            .iter()
            .flat_map(PrismSlab::faces)
            .map(face_area)
            .sum();
        let got: f64 = slab.faces().iter().map(face_area).sum();
        assert!(
            (got - expected).abs() < EXACT_EPS,
            "slab [{}, {}] fused {got} of material; its cover unions to {expected}",
            slab.z_base(),
            slab.z_top(),
        );
    }
}

/// The same reuse must survive cuts: two slabs can share a cover and
/// still differ, because a cut carves only the slab its own z span
/// reaches. A memo that cached the CARVED region instead of the fused
/// material would leak one slab's opening into the other.
#[test]
fn a_shared_cover_still_carves_each_slab_separately() {
    let footprint = rect(0.0, 0.0, 8.0, 1.0);
    let opening = rect(2.0, -1.0, 4.0, 2.0);
    let profile = prism(vec![footprint], 0.0, 3.0).with_cuts(vec![cut(opening, 1.0, 2.0)]);
    let out = run(vec![profile]);

    assert_eq!(out.slabs().len(), 3, "the cut opens sill and head");
    let area = |i: usize| -> f64 { out.slabs()[i].faces().iter().map(face_area).sum() };
    let whole = 8.0;
    assert!((area(0) - whole).abs() < EXACT_EPS, "below the sill: whole");
    assert!(
        (area(1) - (whole - 2.0)).abs() < EXACT_EPS,
        "across the opening the cut is removed, got {}",
        area(1),
    );
    assert!(
        (area(2) - whole).abs() < EXACT_EPS,
        "above the head the band is whole again, got {}",
        area(2),
    );
}
