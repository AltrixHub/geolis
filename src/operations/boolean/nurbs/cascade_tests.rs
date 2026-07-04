//! F6 Phase R1 acceptance: cascaded booleans on one wall.
//!
//! The exact scenario revion's real wall pipeline produces: one
//! `MakeSegmentedPrism` wall (straight rectangular footprint, tagged, with
//! an op id) cut SEQUENTIALLY by three 4-side-face box cutters (a door-sized
//! box plus two window-sized boxes, all tagged, all with op ids, at distinct
//! non-overlapping positions). The second and third `Subtract` operate on an
//! ALREADY-PUNCHED result solid: their target face lists contain punched
//! copies carrying earlier trim holes AND the earlier cuts' band fragments.
//!
//! Pinned contract:
//! - the final solid is position-weld watertight (0 boundary edges) with
//!   3 genuinely open holes;
//! - ALL names resolve: the wall's tagged faces via the transfer chain
//!   through every copy generation, each cut's band fragments under its own
//!   op id (earlier cuts' bands SURVIVE later cuts — they are target faces
//!   then), and each cut's entry/exit rims;
//! - cutting in a different order yields geometry-equal results (sampled
//!   per resolved name) and the same resolving name set.

use std::collections::HashMap;

use crate::math::{Point3, Vector3};
use crate::operations::boolean::Subtract;
use crate::operations::creation::{MakeSegmentedPrism, ProfileSegment};
use crate::tessellation::{TessellateSolid, TessellationParams};
use crate::topology::{
    EdgeName, FaceName, FaceRole, FaceSurface, OpId, SegmentTag, SolidId, TopologyStore,
};

use super::test_support::welded_boundary_edges;

/// Wall footprint: 6.0 x 0.4 at y in [0, 0.4], extruded 3.0 up.
fn wall_profile() -> Vec<ProfileSegment> {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
    vec![
        line(p(0.0, 0.0), p(6.0, 0.0)),
        line(p(6.0, 0.0), p(6.0, 0.4)),
        line(p(6.0, 0.4), p(0.0, 0.4)),
        line(p(0.0, 0.4), p(0.0, 0.0)),
    ]
}

fn wall_tags() -> Vec<SegmentTag> {
    ["outer", "end-east", "inner", "end-west"]
        .iter()
        .map(|t| SegmentTag::new(*t))
        .collect()
}

fn box_tags() -> Vec<SegmentTag> {
    ["sill", "jamb-right", "head", "jamb-left"]
        .iter()
        .map(|t| SegmentTag::new(*t))
        .collect()
}

/// One opening cutter: a box profile in the XZ plane at `y = -1`, spanning
/// `x in [x0, x1]`, `z in [z0, z1]`, extruded 2.4 along `+Y` (through the
/// wall thickness 0.4 with margin on both sides).
#[derive(Clone, Copy)]
struct Opening {
    op: &'static str,
    x0: f64,
    x1: f64,
    z0: f64,
    z1: f64,
}

/// Door-sized cut plus two window-sized cuts, all interior (no cap contact —
/// cap-touching doors are Phase R2) and mutually non-overlapping.
const DOOR: Opening = Opening {
    op: "cut-door",
    x0: 0.6,
    x1: 1.5,
    z0: 0.15,
    z1: 2.25,
};
const WIN_A: Opening = Opening {
    op: "cut-winA",
    x0: 2.2,
    x1: 3.4,
    z0: 1.0,
    z1: 2.0,
};
const WIN_B: Opening = Opening {
    op: "cut-winB",
    x0: 4.2,
    x1: 5.4,
    z0: 1.0,
    z1: 2.0,
};

/// The tool op id of an opening cut (`cut-door` -> `door1`, etc.).
fn tool_op(cut: &Opening) -> String {
    format!("tool-{}", cut.op)
}

fn build_cutter(store: &mut TopologyStore, cut: &Opening) -> SolidId {
    let q = |x: f64, z: f64| Point3::new(x, -1.0, z);
    let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
    let profile = vec![
        line(q(cut.x0, cut.z0), q(cut.x1, cut.z0)), // sill
        line(q(cut.x1, cut.z0), q(cut.x1, cut.z1)), // jamb-right
        line(q(cut.x1, cut.z1), q(cut.x0, cut.z1)), // head
        line(q(cut.x0, cut.z1), q(cut.x0, cut.z0)), // jamb-left
    ];
    MakeSegmentedPrism::new(profile, Vector3::new(0.0, 2.4, 0.0))
        .with_op_id(OpId::new(tool_op(cut)))
        .with_segment_tags(box_tags())
        .execute(store)
        .unwrap()
}

/// Builds wall − openings applied SEQUENTIALLY in the given order; each
/// `Subtract` result feeds the next as its target.
fn cascade(order: &[Opening]) -> (TopologyStore, SolidId) {
    let mut store = TopologyStore::new();
    let wall = MakeSegmentedPrism::new(wall_profile(), Vector3::new(0.0, 0.0, 3.0))
        .with_op_id(OpId::new("wall1"))
        .with_segment_tags(wall_tags())
        .execute(&mut store)
        .unwrap();
    let mut current = wall;
    for cut in order {
        let cutter = build_cutter(&mut store, cut);
        current = Subtract::new(current, cutter)
            .with_op_id(OpId::new(cut.op))
            .execute(&mut store)
            .unwrap_or_else(|e| panic!("{} failed: {e:?}", cut.op));
    }
    (store, current)
}

/// Every persistent face name the cascade must resolve: the wall's 4 tagged
/// faces plus one band fragment per (cut, tool side face).
fn expected_face_names(cuts: &[Opening]) -> Vec<FaceName> {
    let mut names: Vec<FaceName> = ["outer", "end-east", "inner", "end-west"]
        .iter()
        .map(|tag| FaceName::Created {
            op: OpId::new("wall1"),
            role: FaceRole::Tagged(SegmentTag::new(*tag)),
        })
        .collect();
    for cut in cuts {
        for tag in ["sill", "jamb-right", "head", "jamb-left"] {
            names.push(FaceName::Band {
                op: OpId::new(cut.op),
                tool_face: Box::new(FaceName::Created {
                    op: OpId::new(tool_op(cut)),
                    role: FaceRole::Tagged(SegmentTag::new(tag)),
                }),
                loop_index: 0,
            });
        }
    }
    names
}

/// Every persistent edge name the cascade must resolve: entry (outer, loop 0)
/// and exit (inner, loop 1) rims per cut.
fn expected_edge_names(cuts: &[Opening]) -> Vec<EdgeName> {
    let mut names = Vec::new();
    for cut in cuts {
        for (target_tag, loop_index) in [("outer", 0u32), ("inner", 1u32)] {
            names.push(EdgeName::CutRim {
                op: OpId::new(cut.op),
                target: Box::new(FaceName::Created {
                    op: OpId::new("wall1"),
                    role: FaceRole::Tagged(SegmentTag::new(target_tag)),
                }),
                loop_index,
            });
        }
    }
    names
}

/// Samples a resolved named face's surface at fixed UV probes (full
/// precision, used for order-independence geometry comparison).
fn face_samples(store: &TopologyStore, name: &FaceName) -> Vec<Point3> {
    let face_id = store
        .names()
        .face(name)
        .unwrap_or_else(|| panic!("{name:?} does not resolve"));
    let face = store.face(face_id).unwrap();
    let FaceSurface::Nurbs(surf) = &face.surface else {
        panic!(
            "{name:?} must resolve to a NURBS face, got {:?}",
            face.surface
        );
    };
    [(0.31, 0.62), (0.5, 0.5), (0.87, 0.13)]
        .iter()
        .map(|&(u, v)| surf.point_at(u, v).unwrap())
        .collect()
}

/// R1 acceptance 1: wall − door − window − window (sequential) is
/// position-weld watertight with 3 genuinely open holes, and every name
/// resolves into the final result shell.
#[test]
fn cascaded_wall_three_openings_is_watertight_with_open_holes() {
    let cuts = [DOOR, WIN_A, WIN_B];
    let (store, result) = cascade(&cuts);

    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();
    // 6 wall faces (4 sides + 2 caps) + 3 cuts x 4 band fragments.
    assert_eq!(
        shell.faces.len(),
        18,
        "6 wall faces + 12 band fragments, got {}",
        shell.faces.len()
    );

    // The outer and inner wall faces each accumulated 3 chained hole rings.
    let mut hole_counts: Vec<usize> = shell
        .faces
        .iter()
        .map(|&f| store.face(f).unwrap().inner_wires.len())
        .filter(|&n| n > 0)
        .collect();
    hole_counts.sort_unstable();
    assert_eq!(
        hole_counts,
        vec![3, 3],
        "outer + inner faces carry 3 hole rings each"
    );

    // Watertight: zero position-welded boundary edges.
    let boundary = welded_boundary_edges(&store, result);
    assert_eq!(
        boundary, 0,
        "cascaded wall must position-weld watertight (found {boundary} \
         boundary edges)"
    );

    // Each hole is genuinely open: no mesh vertex intrudes into any tunnel.
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    for cut in &cuts {
        let (x0, x1) = (cut.x0 + 0.05, cut.x1 - 0.05);
        let (z0, z1) = (cut.z0 + 0.05, cut.z1 - 0.05);
        for v in &mesh.vertices {
            let inside = v.x > x0 && v.x < x1 && v.z > z0 && v.z < z1;
            assert!(
                !(inside && v.y > 0.05 && v.y < 0.35),
                "vertex ({:.3},{:.3},{:.3}) intrudes into the {} tunnel",
                v.x,
                v.y,
                v.z,
                cut.op
            );
        }
    }

    // All face names resolve INTO the final result shell: the wall tags via
    // the transfer chain, and every cut's band fragments — the door's and
    // window A's bands survived the later cuts (they were target faces).
    for name in expected_face_names(&cuts) {
        let f = store
            .names()
            .face(&name)
            .unwrap_or_else(|| panic!("{name:?} does not resolve after the cascade"));
        assert!(
            shell.faces.contains(&f),
            "{name:?} resolves outside the final result shell"
        );
    }
    // All rim edge names resolve.
    for name in expected_edge_names(&cuts) {
        assert!(
            store.names().edge(&name).is_some(),
            "{name:?} does not resolve after the cascade"
        );
    }
}

/// R1 acceptance 2 (order independence): cutting in a different order yields
/// geometry-equal results (full-precision sampled comparison per resolved
/// name) and the same resolving name set.
#[test]
fn cascade_order_does_not_change_geometry_or_names() {
    let orders: [[Opening; 3]; 3] = [
        [DOOR, WIN_A, WIN_B],
        [WIN_B, WIN_A, DOOR],
        [WIN_A, DOOR, WIN_B],
    ];
    let reference_cuts = orders[0];
    let names = expected_face_names(&reference_cuts);
    let edge_names = expected_edge_names(&reference_cuts);

    let (ref_store, ref_result) = cascade(&orders[0]);
    let ref_samples: HashMap<String, Vec<Point3>> = names
        .iter()
        .map(|n| (format!("{n:?}"), face_samples(&ref_store, n)))
        .collect();
    let ref_shell_len = ref_store
        .shell(ref_store.solid(ref_result).unwrap().outer_shell)
        .unwrap()
        .faces
        .len();

    for order in &orders[1..] {
        let (store, result) = cascade(order);
        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();
        assert_eq!(shell.faces.len(), ref_shell_len, "face count differs");
        assert_eq!(
            welded_boundary_edges(&store, result),
            0,
            "reordered cascade must stay watertight"
        );
        for name in &names {
            let samples = face_samples(&store, name);
            let reference = &ref_samples[&format!("{name:?}")];
            for (s, r) in samples.iter().zip(reference) {
                assert!(
                    (*s - *r).norm() < 1e-9,
                    "{name:?} geometry differs across cut orders: {s:?} vs {r:?}"
                );
            }
        }
        for name in &edge_names {
            assert!(
                store.names().edge(name).is_some(),
                "{name:?} does not resolve in reordered cascade"
            );
        }
    }
}

/// Pinned characterization of the F5 lead "near-duplicate unwelded vertices
/// on the sill line" (verified real, 2026-07-04, F6 R1):
///
/// The wall − box-window mesh carries cross-face vertex pairs along the
/// shared hole-rim edges (sill, head, jambs) separated by up to ~7.9e-10.
/// Root cause: the trimmed CDT maps every vertex through its OWN face's
/// surface — the punched wall face evaluates `S_target(uv_a)` from the trim
/// hole polyline while the band fragment evaluates `S_tool(uv_b)` via its
/// pcurve — so the SSI marcher's point-acceptance residual (bounded by its
/// 1e-7 acceptance, `stitch::JUNCTION_TOLERANCE`) surfaces as a cross-face
/// mismatch. The true weld belongs in the tessellation layer (shared-edge 3D
/// pinning so both faces consume one 3D sample per shared-edge parameter) —
/// scheduled as part of F6 R3.
///
/// Until then this test pins the bound: every near-duplicate pair stays
/// within the marcher's acceptance, so the 1e-6 position weld (the shipped
/// watertightness contract) always merges them. A regression pushing rim
/// disagreement past the marcher bound fails here.
#[test]
fn rim_cross_face_mismatch_stays_within_marcher_acceptance() {
    /// The SSI marcher's point-acceptance bound (`stitch::JUNCTION_TOLERANCE`).
    const MARCHER_ACCEPTANCE: f64 = 1e-7;

    let mut store = TopologyStore::new();
    let wall = MakeSegmentedPrism::new(wall_profile(), Vector3::new(0.0, 0.0, 3.0))
        .with_op_id(OpId::new("wall1"))
        .with_segment_tags(wall_tags())
        .execute(&mut store)
        .unwrap();
    let cutter = build_cutter(&mut store, &WIN_A);
    let result = Subtract::new(wall, cutter)
        .with_op_id(OpId::new("cut1"))
        .execute(&mut store)
        .unwrap();
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();

    // Any two distinct-but-nearby vertices are the two faces' images of one
    // shared rim sample; their separation must stay within the marcher's
    // acceptance so the 1e-6 position weld closes them.
    let mut max_near_dup = 0.0_f64;
    for i in 0..mesh.vertices.len() {
        for j in (i + 1)..mesh.vertices.len() {
            let d = (mesh.vertices[i] - mesh.vertices[j]).norm();
            if d > 0.0 && d < 1e-4 {
                max_near_dup = max_near_dup.max(d);
            }
        }
    }
    assert!(
        max_near_dup <= MARCHER_ACCEPTANCE,
        "cross-face rim vertex mismatch {max_near_dup:.3e} exceeds the SSI \
         marcher's acceptance {MARCHER_ACCEPTANCE:.0e}; the punch/band rim \
         images have drifted apart"
    );
}
