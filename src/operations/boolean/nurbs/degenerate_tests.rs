//! F6 Phase R3 acceptance: degenerate contacts.
//!
//! Each degenerate cutter placement is either a CLEAN cut or a PINNED typed
//! error with a full-precision regression test — no silent wrong geometry.
//!
//! ## Verdicts
//!
//! | Contact | Verdict |
//! |---|---|
//! | Flush sill (`z0 == 0.0`, sill coplanar with the bottom cap) | **CLEAN** — degenerates to the R2 full-height-door notch: the sill's tangential along-boundary branches are dropped by the loop extraction, and the jamb corner terminals (tool kink on the target boundary) are accepted as target-boundary terminals |
//! | Flush sill AND head (`z0 == 0.0`, `z1 == 3.0`, coplanar with both caps) | **CLEAN** — degenerates to the R2 floor-to-ceiling door (wall faces sever, both caps notch) |
//! | Wall-end-flush opening (cutter jamb coplanar with the wall's end face) | **PINNED typed error** (`open branch`) — the coplanar jamb yields no transversal SSI branches that close the circuit around the wall end, and the end face is a NURBS side face, not a notchable planar cap |
//! | Boundary-vertex hit (cutter jamb through a collinear profile joint, cut trace on the shared joint edge) | **PINNED typed error** (`splits a boundary edge shared with a face the cut does not cross`) — the jamb trace collapses onto the existing joint edge (dropped as an along-boundary tangential contact), the sill / head traces terminate on the joint, and the splitter refuses to split the joint edge shared with the untouched neighbour face |
//! | Near-degenerate band (margin > 0 but smaller than one SSI marcher step, e.g. sill 0.01 above the cap) | **PINNED typed error** (`open branch`) — pre-existing marcher behaviour: the branch terminates one step short of the tool kink with an interior endpoint. Margins of one marcher step or more cut cleanly (see the sanity test) |

use crate::math::{Point3, Vector3};
use crate::operations::boolean::Subtract;
use crate::operations::creation::{MakeSegmentedPrism, ProfileSegment};
use crate::tessellation::{TessellateSolid, TessellationParams};
use crate::topology::{
    EdgeName, FaceName, FaceRole, OpId, SegmentTag, SolidId, SplitSide, TopologyStore,
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

/// Builds one box cutter spanning `x in [x0, x1]`, `z in [z0, z1]`,
/// extruded 2.4 along `+Y` through the wall thickness.
fn build_cutter(store: &mut TopologyStore, x0: f64, x1: f64, z0: f64, z1: f64) -> SolidId {
    let q = |x: f64, z: f64| Point3::new(x, -1.0, z);
    let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
    let profile = vec![
        line(q(x0, z0), q(x1, z0)), // sill
        line(q(x1, z0), q(x1, z1)), // jamb-right
        line(q(x1, z1), q(x0, z1)), // head
        line(q(x0, z1), q(x0, z0)), // jamb-left
    ];
    MakeSegmentedPrism::new(profile, Vector3::new(0.0, 2.4, 0.0))
        .with_op_id(OpId::new("tool1"))
        .with_segment_tags(box_tags())
        .execute(store)
        .unwrap()
}

/// Builds the straight wall and one box cutter.
fn wall_and_cutter(x0: f64, x1: f64, z0: f64, z1: f64) -> (TopologyStore, SolidId, SolidId) {
    let mut store = TopologyStore::new();
    let wall = MakeSegmentedPrism::new(wall_profile(), Vector3::new(0.0, 0.0, 3.0))
        .with_op_id(OpId::new("wall1"))
        .with_segment_tags(wall_tags())
        .execute(&mut store)
        .unwrap();
    let cutter = build_cutter(&mut store, x0, x1, z0, z1);
    (store, wall, cutter)
}

fn wall_tag_name(tag: &str) -> FaceName {
    FaceName::Created {
        op: OpId::new("wall1"),
        role: FaceRole::Tagged(SegmentTag::new(tag)),
    }
}

fn cap_name(role: FaceRole) -> FaceName {
    FaceName::Created {
        op: OpId::new("wall1"),
        role,
    }
}

fn band_name(tool_tag: &str) -> FaceName {
    FaceName::Band {
        op: OpId::new("cut1"),
        tool_face: Box::new(FaceName::Created {
            op: OpId::new("tool1"),
            role: FaceRole::Tagged(SegmentTag::new(tool_tag)),
        }),
        loop_index: 0,
    }
}

fn split_name(parent: FaceName, side: SplitSide) -> FaceName {
    FaceName::Split {
        op: OpId::new("cut1"),
        parent: Box::new(parent),
        side,
    }
}

/// Asserts no mesh vertex intrudes into the opening tunnel (margins shrink
/// the probe box so rim vertices do not trip it).
fn assert_open(mesh: &crate::tessellation::TriangleMesh, x0: f64, x1: f64, z0: f64, z1: f64) {
    let (x0, x1) = (x0 + 0.05, x1 - 0.05);
    let (z0, z1) = (z0.max(0.0) + 0.05, z1.min(3.0) - 0.05);
    for v in &mesh.vertices {
        let inside = v.x > x0 && v.x < x1 && v.z > z0 && v.z < z1;
        assert!(
            !(inside && v.y > 0.05 && v.y < 0.35),
            "vertex ({:.3},{:.3},{:.3}) intrudes into the opening",
            v.x,
            v.y,
            v.z,
        );
    }
}

/// R3 flush sill (CLEAN): a door whose sill plane is exactly coplanar with
/// the wall's bottom cap (`z0 == 0.0`) degenerates to the R2 full-height
/// door notch — watertight, open doorway, the bottom cap split into two
/// named fragments, jamb / head bands named, NO sill band.
#[test]
fn flush_sill_door_degenerates_to_full_height_notch() {
    let (x0, x1, z0, z1) = (0.6, 1.5, 0.0, 2.25);
    let (mut store, wall, cutter) = wall_and_cutter(x0, x1, z0, z1);
    let result = Subtract::new(wall, cutter)
        .with_op_id(OpId::new("cut1"))
        .execute(&mut store)
        .unwrap_or_else(|e| panic!("flush sill must cut cleanly: {e:?}"));

    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();
    // outer + inner notched faces, 2 ends, top cap, 2 bottom cap fragments,
    // 3 band fragments (jamb-left / head / jamb-right) — exactly the R2
    // full-height-door topology.
    assert_eq!(shell.faces.len(), 10, "got {}", shell.faces.len());

    let boundary = welded_boundary_edges(&store, result);
    assert_eq!(boundary, 0, "flush sill must weld watertight ({boundary})");

    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, x0, x1, z0, z1);
    // No stray geometry below the wall (the flush sill contributes nothing).
    for v in &mesh.vertices {
        assert!(v.z >= -1e-12, "vertex below the wall base: {v:?}");
    }

    // Wall side / end faces keep their tagged names (notch → transfer).
    for tag in ["outer", "end-east", "inner", "end-west"] {
        let f = store
            .names()
            .face(&wall_tag_name(tag))
            .unwrap_or_else(|| panic!("{tag} does not resolve"));
        assert!(shell.faces.contains(&f), "{tag} outside the result shell");
    }
    // Top cap survives; the bottom cap split into Left / Right.
    assert!(store.names().face(&cap_name(FaceRole::CapEnd)).is_some());
    assert!(
        store.names().face(&cap_name(FaceRole::CapStart)).is_none(),
        "CapStart must retire when the cap splits"
    );
    for side in [SplitSide::Left, SplitSide::Right] {
        let name = split_name(cap_name(FaceRole::CapStart), side);
        let f = store
            .names()
            .face(&name)
            .unwrap_or_else(|| panic!("{name:?} does not resolve"));
        assert!(shell.faces.contains(&f), "{name:?} outside the shell");
    }
    // Jamb / head bands resolve; NO sill band (the sill is coplanar with
    // the cap and bounds no hole wall).
    for tag in ["jamb-left", "head", "jamb-right"] {
        assert!(
            store.names().face(&band_name(tag)).is_some(),
            "band {tag} does not resolve"
        );
    }
    assert!(
        store.names().face(&band_name("sill")).is_none(),
        "a sill band must not exist for a flush-sill door"
    );
    // Entry / exit rims resolve.
    for (target, loop_index) in [("outer", 0u32), ("inner", 1u32)] {
        let rim = EdgeName::CutRim {
            op: OpId::new("cut1"),
            target: Box::new(wall_tag_name(target)),
            loop_index,
        };
        assert!(store.names().edge(&rim).is_some(), "{rim:?} unresolved");
    }
}

/// R3 flush sill + flush head (CLEAN): a door coplanar with BOTH caps
/// (`z0 == 0.0`, `z1 == 3.0`) degenerates to the R2 floor-to-ceiling door:
/// the wall faces sever into Split L/R and both caps notch.
#[test]
fn flush_sill_and_head_door_severs_wall_faces() {
    let (x0, x1, z0, z1) = (0.6, 1.5, 0.0, 3.0);
    let (mut store, wall, cutter) = wall_and_cutter(x0, x1, z0, z1);
    let result = Subtract::new(wall, cutter)
        .with_op_id(OpId::new("cut1"))
        .execute(&mut store)
        .unwrap_or_else(|e| panic!("flush sill+head must cut cleanly: {e:?}"));

    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();
    // outer L/R + inner L/R + 2 ends + 2 fragments per cap + 2 jamb bands.
    assert_eq!(shell.faces.len(), 12, "got {}", shell.faces.len());
    assert_eq!(
        welded_boundary_edges(&store, result),
        0,
        "flush sill+head must weld watertight"
    );
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, x0, x1, z0, z1);

    // The severed outer / inner faces retire; Split L/R resolve.
    for tag in ["outer", "inner"] {
        assert!(
            store.names().face(&wall_tag_name(tag)).is_none(),
            "{tag} must retire when the door severs it"
        );
        for side in [SplitSide::Left, SplitSide::Right] {
            assert!(
                store
                    .names()
                    .face(&split_name(wall_tag_name(tag), side))
                    .is_some(),
                "{tag} {side:?} does not resolve"
            );
        }
    }
    // Both caps split; only the jamb bands exist.
    for role in [FaceRole::CapStart, FaceRole::CapEnd] {
        for side in [SplitSide::Left, SplitSide::Right] {
            assert!(
                store
                    .names()
                    .face(&split_name(cap_name(role.clone()), side))
                    .is_some(),
                "{role:?} {side:?} does not resolve"
            );
        }
    }
    for tag in ["jamb-left", "jamb-right"] {
        assert!(store.names().face(&band_name(tag)).is_some());
    }
    for tag in ["sill", "head"] {
        assert!(
            store.names().face(&band_name(tag)).is_none(),
            "band {tag} must not exist"
        );
    }
}

/// R3 wall-end-flush (PINNED typed error): a cutter whose jamb is exactly
/// coplanar with the wall's end face (`x1 == 6.0`). The coplanar jamb
/// yields no SSI branch on the end face, so the sill / head branches
/// dead-end at the wall-end corner — there is no cap to notch (the end is
/// a NURBS side face). Historically the marcher stopped one step short of
/// the corner (open-branch error); the analytic extrusion×plane path lands
/// the endpoints exactly on the corner, so chaining proceeds one stage
/// further and rejects the configuration as a chain crossing one tool face
/// twice. Still a typed error, never silent wrong geometry.
#[test]
fn wall_end_flush_opening_is_a_typed_error() {
    let (mut store, wall, cutter) = wall_and_cutter(5.0, 6.0, 0.5, 2.0);
    let result = Subtract::new(wall, cutter)
        .with_op_id(OpId::new("cut1"))
        .execute(&mut store);
    let err = result.expect_err("wall-end-flush opening must be a typed error");
    let message = format!("{err}");
    assert!(
        message.contains("crosses one tool side face twice"),
        "unexpected error for the wall-end-flush opening: {message}"
    );
}

/// R3 boundary-vertex hit (PINNED typed error): a two-segment collinear
/// wall (joint at x = 3) cut by a window whose jamb passes exactly through
/// the profile joint — the cut trace lands on the shared joint edge. The
/// jamb trace collapses onto the joint edge (an along-boundary tangential
/// contact), the sill / head traces terminate on the joint, and the
/// splitter refuses to split the joint edge shared with the untouched
/// neighbour face. The typed error is pinned; no silent wrong geometry.
#[test]
fn cut_trace_on_profile_joint_vertex_is_a_typed_error() {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
    // Straight wall split into two collinear outer / inner segments with a
    // joint at x = 3.
    let profile = vec![
        line(p(0.0, 0.0), p(3.0, 0.0)),
        line(p(3.0, 0.0), p(6.0, 0.0)),
        line(p(6.0, 0.0), p(6.0, 0.4)),
        line(p(6.0, 0.4), p(3.0, 0.4)),
        line(p(3.0, 0.4), p(0.0, 0.4)),
        line(p(0.0, 0.4), p(0.0, 0.0)),
    ];
    let mut store = TopologyStore::new();
    let wall = MakeSegmentedPrism::new(profile, Vector3::new(0.0, 0.0, 3.0))
        .with_op_id(OpId::new("wall1"))
        .execute(&mut store)
        .unwrap();
    // Window jamb-right exactly on the joint plane x = 3.
    let cutter = build_cutter(&mut store, 2.2, 3.0, 1.0, 2.0);
    let result = Subtract::new(wall, cutter)
        .with_op_id(OpId::new("cut1"))
        .execute(&mut store);
    let err = result.expect_err("cut trace on a profile joint must be a typed error");
    let message = format!("{err}");
    assert!(
        message.contains("splits a boundary edge shared with a face the cut does not cross"),
        "unexpected error for the joint-vertex cut: {message}"
    );
}

/// R3 sanity: the same geometry moved OFF the degenerate contacts by one
/// SSI marcher step or more keeps cutting cleanly (the degenerate handling
/// introduces no regression in the ordinary neighbourhood).
#[test]
fn near_degenerate_neighbours_still_cut_cleanly() {
    // Sill one marcher step above the cap: ordinary interior window cut.
    let (mut store, wall, cutter) = wall_and_cutter(0.6, 1.5, 0.15, 2.25);
    let result = Subtract::new(wall, cutter)
        .with_op_id(OpId::new("cut1"))
        .execute(&mut store)
        .unwrap_or_else(|e| panic!("near-flush sill must cut cleanly: {e:?}"));
    assert_eq!(welded_boundary_edges(&store, result), 0);

    // Jamb one marcher step short of the end face: ordinary interior window.
    let (mut store, wall, cutter) = wall_and_cutter(5.0, 5.9, 0.5, 2.0);
    let result = Subtract::new(wall, cutter)
        .with_op_id(OpId::new("cut1"))
        .execute(&mut store)
        .unwrap_or_else(|e| panic!("near-end jamb must cut cleanly: {e:?}"));
    assert_eq!(welded_boundary_edges(&store, result), 0);
}

/// R3 near-degenerate band: a margin greater than zero but smaller than
/// one SSI marcher step (sill 0.01 above the cap). Historically a PINNED
/// typed error — the marcher left the jamb branch terminating one step
/// short of the tool kink with an interior endpoint. The analytic
/// extrusion×plane fast path computes the exact boundary crossing, so this
/// margin now cuts CLEANLY like its neighbours (flush contact and
/// one-step-plus margins were already clean). Pinned as a clean, manifold
/// cut — never silent wrong geometry.
#[test]
fn sub_step_margin_cuts_cleanly() {
    let (mut store, wall, cutter) = wall_and_cutter(0.6, 1.5, 0.01, 2.25);
    let result = Subtract::new(wall, cutter)
        .with_op_id(OpId::new("cut1"))
        .execute(&mut store)
        .unwrap_or_else(|e| panic!("sub-step margin must cut cleanly: {e:?}"));
    assert_eq!(welded_boundary_edges(&store, result), 0);
}
