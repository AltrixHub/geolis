//! F6 Phase R2 acceptance: cap-touching (full-height door) openings.
//!
//! A door box reaching below the wall's bottom cap cuts OPEN chains: the
//! entry / exit traces are boundary-to-boundary notches on the wall side
//! faces, the band fragments close the doorway with cap-plane closure
//! edges, and the notched bottom cap is rebuilt by wire surgery — the SAME
//! `EdgeId`s on the cap, the wall fragments, and the band fragments, so the
//! doorway is watertight by construction.
//!
//! Pinned contract:
//! - a bottom-touching door (head inside the wall) keeps the wall faces'
//!   tagged names (boundary notch → transfer), splits the bottom cap into
//!   `Split { Left / Right }` of `CapStart`, names the jamb / head bands,
//!   and binds the entry / exit rims; the result is position-weld
//!   watertight with a genuinely open doorway, and rebuild-stable;
//! - a floor-to-ceiling door splits the wall front / back faces into
//!   `Split { Left / Right }` and notches BOTH caps;
//! - cascades keep working after a door: door + window in a split fragment
//!   is watertight and geometry-equal across cut orders (the window hole
//!   transfers onto the containing fragment when the door splits a punched
//!   face), and a second door notches the already-notched faces (the
//!   trimmed-outer perimeter generalization).

use std::collections::HashMap;

use crate::math::{Point3, Vector3};
use crate::operations::boolean::Subtract;
use crate::operations::creation::{MakeSegmentedPrism, ProfileSegment};
use crate::tessellation::{TessellateSolid, TessellationParams};
use crate::topology::{
    EdgeName, FaceName, FaceRole, FaceSurface, OpId, SegmentTag, SolidId, SplitSide, TopologyStore,
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

/// One opening cutter spanning `x in [x0, x1]`, `z in [z0, z1]`, extruded
/// 2.4 along `+Y` through the wall thickness. `z0 < 0` makes it
/// cap-touching (the sill face lies below the wall).
#[derive(Clone, Copy)]
struct Opening {
    op: &'static str,
    x0: f64,
    x1: f64,
    z0: f64,
    z1: f64,
}

/// Bottom-touching door: sill below the wall, head inside it.
const DOOR: Opening = Opening {
    op: "cut-door",
    x0: 0.6,
    x1: 1.5,
    z0: -0.5,
    z1: 2.25,
};

/// Floor-to-ceiling door: sill below the wall, head above it.
const TALL_DOOR: Opening = Opening {
    op: "cut-door",
    x0: 0.6,
    x1: 1.5,
    z0: -0.5,
    z1: 3.5,
};

/// Second bottom-touching door (for the door-then-door cascade).
const DOOR_B: Opening = Opening {
    op: "cut-doorB",
    x0: 3.0,
    x1: 4.0,
    z0: -0.5,
    z1: 2.25,
};

/// Interior window inside the tall door's high-x split fragment.
const WINDOW: Opening = Opening {
    op: "cut-win",
    x0: 2.2,
    x1: 3.4,
    z0: 1.0,
    z1: 2.0,
};

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

/// Builds wall − openings applied SEQUENTIALLY in the given order.
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

fn band_name(cut: &Opening, tool_tag: &str) -> FaceName {
    FaceName::Band {
        op: OpId::new(cut.op),
        tool_face: Box::new(FaceName::Created {
            op: OpId::new(tool_op(cut)),
            role: FaceRole::Tagged(SegmentTag::new(tool_tag)),
        }),
        loop_index: 0,
    }
}

fn split_name(op: &str, parent: FaceName, side: SplitSide) -> FaceName {
    FaceName::Split {
        op: OpId::new(op),
        parent: Box::new(parent),
        side,
    }
}

/// A geometric sample of a resolved named face: NURBS faces sample the
/// surface at fixed UVs; planar faces (cap fragments) average their outer
/// wire's vertices — which distinguishes the Left / Right fragments sharing
/// one plane.
fn face_sample(store: &TopologyStore, name: &FaceName) -> Vec<Point3> {
    let face_id = store
        .names()
        .face(name)
        .unwrap_or_else(|| panic!("{name:?} does not resolve"));
    let face = store.face(face_id).unwrap();
    match &face.surface {
        FaceSurface::Nurbs(surf) => [(0.31, 0.62), (0.5, 0.5), (0.87, 0.13)]
            .iter()
            .map(|&(u, v)| surf.point_at(u, v).unwrap())
            .collect(),
        FaceSurface::Plane(_) => {
            let wire = store.wire(face.outer_wire).unwrap();
            let mut sum = Point3::new(0.0, 0.0, 0.0);
            let mut count = 0.0;
            for oe in &wire.edges {
                let edge = store.edge(oe.edge).unwrap();
                let p = store.vertex(edge.start).unwrap().point;
                sum = Point3::new(sum.x + p.x, sum.y + p.y, sum.z + p.z);
                count += 1.0;
            }
            vec![Point3::new(sum.x / count, sum.y / count, sum.z / count)]
        }
        other => panic!("{name:?} resolved to unexpected surface {other:?}"),
    }
}

/// Asserts no mesh vertex intrudes into the opening tunnel (margins shrink
/// the probe box so rim vertices do not trip it).
fn assert_open(mesh: &crate::tessellation::TriangleMesh, cut: &Opening) {
    let (x0, x1) = (cut.x0 + 0.05, cut.x1 - 0.05);
    let (z0, z1) = (cut.z0.max(0.0) + 0.05, cut.z1.min(3.0) - 0.05);
    for v in &mesh.vertices {
        let inside = v.x > x0 && v.x < x1 && v.z > z0 && v.z < z1;
        assert!(
            !(inside && v.y > 0.05 && v.y < 0.35),
            "vertex ({:.3},{:.3},{:.3}) intrudes into the {} opening",
            v.x,
            v.y,
            v.z,
            cut.op
        );
    }
}

/// R2 acceptance 1: wall − bottom-touching door. Watertight, open doorway,
/// the bottom cap notched into two named fragments sharing the closure and
/// sub-edges with the bands and the notched wall faces, all names resolve.
#[test]
fn full_height_door_notches_bottom_cap() {
    let (store, result) = cascade(&[DOOR]);
    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();

    // outer + inner notched fragments, 2 ends, top cap, 2 bottom cap
    // fragments, 3 band fragments (jamb-left / head / jamb-right).
    assert_eq!(shell.faces.len(), 10, "got {}", shell.faces.len());

    let boundary = welded_boundary_edges(&store, result);
    assert_eq!(boundary, 0, "door wall must weld watertight ({boundary})");

    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &DOOR);

    // Wall side / end faces keep their tagged names (notch → transfer).
    for tag in ["outer", "end-east", "inner", "end-west"] {
        let f = store
            .names()
            .face(&wall_tag_name(tag))
            .unwrap_or_else(|| panic!("{tag} does not resolve"));
        assert!(shell.faces.contains(&f), "{tag} outside the result shell");
    }
    // Top cap survives unchanged; the bottom cap split into Left / Right.
    assert!(store.names().face(&cap_name(FaceRole::CapEnd)).is_some());
    assert!(
        store.names().face(&cap_name(FaceRole::CapStart)).is_none(),
        "CapStart must retire when the cap splits"
    );
    for side in [SplitSide::Left, SplitSide::Right] {
        let name = split_name(DOOR.op, cap_name(FaceRole::CapStart), side);
        let f = store
            .names()
            .face(&name)
            .unwrap_or_else(|| panic!("{name:?} does not resolve"));
        assert!(shell.faces.contains(&f), "{name:?} outside the shell");
        assert!(
            matches!(store.face(f).unwrap().surface, FaceSurface::Plane(_)),
            "cap fragment must stay planar"
        );
    }
    // Jamb / head bands resolve; there is no sill band (the sill face lies
    // below the wall).
    for tag in ["jamb-left", "head", "jamb-right"] {
        assert!(
            store.names().face(&band_name(&DOOR, tag)).is_some(),
            "band {tag} does not resolve"
        );
    }
    assert!(
        store.names().face(&band_name(&DOOR, "sill")).is_none(),
        "a sill band must not exist for a cap-touching door"
    );
    // Entry (outer, loop 0) and exit (inner, loop 1) rims resolve.
    for (target, loop_index) in [("outer", 0u32), ("inner", 1u32)] {
        let rim = EdgeName::CutRim {
            op: OpId::new(DOOR.op),
            target: Box::new(wall_tag_name(target)),
            loop_index,
        };
        assert!(store.names().edge(&rim).is_some(), "{rim:?} unresolved");
    }

    // Shared-edge topology: each bottom-cap fragment shares at least one
    // edge with the notched outer face (a kept sub-edge) and at least one
    // with a band fragment (a cap-plane closure edge).
    let wire_edges = |f| -> Vec<crate::topology::EdgeId> {
        store
            .wire(store.face(f).unwrap().outer_wire)
            .unwrap()
            .edges
            .iter()
            .map(|oe| oe.edge)
            .collect()
    };
    let outer_edges = wire_edges(store.names().face(&wall_tag_name("outer")).unwrap());
    let band_edges: Vec<crate::topology::EdgeId> = ["jamb-left", "jamb-right"]
        .iter()
        .flat_map(|tag| wire_edges(store.names().face(&band_name(&DOOR, tag)).unwrap()))
        .collect();
    for side in [SplitSide::Left, SplitSide::Right] {
        let name = split_name(DOOR.op, cap_name(FaceRole::CapStart), side);
        let cap_edges = wire_edges(store.names().face(&name).unwrap());
        assert!(
            cap_edges.iter().any(|e| outer_edges.contains(e)),
            "{name:?} shares no sub-edge with the notched outer face"
        );
        assert!(
            cap_edges.iter().any(|e| band_edges.contains(e)),
            "{name:?} shares no closure edge with the jamb bands"
        );
    }
}

/// R2 acceptance 1b (F4 pattern): rebuilding the same door wall resolves
/// every persistent name to geometrically identical entities.
#[test]
fn full_height_door_is_rebuild_stable() {
    let names: Vec<FaceName> = {
        let mut names = vec![
            wall_tag_name("outer"),
            wall_tag_name("inner"),
            cap_name(FaceRole::CapEnd),
            split_name(DOOR.op, cap_name(FaceRole::CapStart), SplitSide::Left),
            split_name(DOOR.op, cap_name(FaceRole::CapStart), SplitSide::Right),
        ];
        for tag in ["jamb-left", "head", "jamb-right"] {
            names.push(band_name(&DOOR, tag));
        }
        names
    };

    let (store_a, _) = cascade(&[DOOR]);
    let (store_b, _) = cascade(&[DOOR]);
    for name in &names {
        let a = face_sample(&store_a, name);
        let b = face_sample(&store_b, name);
        for (pa, pb) in a.iter().zip(&b) {
            assert!(
                (*pa - *pb).norm() < 1e-9,
                "{name:?} moved across rebuilds: {pa:?} vs {pb:?}"
            );
        }
    }
}

/// R2 acceptance 2: wall − floor-to-ceiling door. The wall front / back
/// faces split into Left / Right, BOTH caps are notched, and the result is
/// watertight with an open full-height doorway.
#[test]
fn floor_to_ceiling_door_splits_wall_faces() {
    let (store, result) = cascade(&[TALL_DOOR]);
    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();

    // outer L/R + inner L/R + 2 ends + 2 fragments per cap + 2 jamb bands.
    assert_eq!(shell.faces.len(), 12, "got {}", shell.faces.len());

    let boundary = welded_boundary_edges(&store, result);
    assert_eq!(boundary, 0, "tall door must weld watertight ({boundary})");

    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &TALL_DOOR);

    // The severed outer / inner faces retire; Split L/R resolve.
    for tag in ["outer", "inner"] {
        assert!(
            store.names().face(&wall_tag_name(tag)).is_none(),
            "{tag} must retire when the door severs it"
        );
        for side in [SplitSide::Left, SplitSide::Right] {
            let name = split_name(TALL_DOOR.op, wall_tag_name(tag), side);
            let f = store
                .names()
                .face(&name)
                .unwrap_or_else(|| panic!("{name:?} does not resolve"));
            assert!(shell.faces.contains(&f), "{name:?} outside the shell");
            // Split fragments bind their rims per fragment.
            let loop_index = u32::from(tag == "inner");
            let rim = EdgeName::CutRim {
                op: OpId::new(TALL_DOOR.op),
                target: Box::new(name.clone()),
                loop_index,
            };
            assert!(store.names().edge(&rim).is_some(), "{rim:?} unresolved");
        }
    }
    // Both caps split.
    for role in [FaceRole::CapStart, FaceRole::CapEnd] {
        assert!(
            store.names().face(&cap_name(role.clone())).is_none(),
            "{role:?} must retire when the cap splits"
        );
        for side in [SplitSide::Left, SplitSide::Right] {
            let name = split_name(TALL_DOOR.op, cap_name(role.clone()), side);
            assert!(
                store.names().face(&name).is_some(),
                "{name:?} does not resolve"
            );
        }
    }
    // Only the jamb bands exist (sill below, head above the wall).
    for tag in ["jamb-left", "jamb-right"] {
        assert!(
            store.names().face(&band_name(&TALL_DOOR, tag)).is_some(),
            "band {tag} does not resolve"
        );
    }
    for tag in ["sill", "head"] {
        assert!(
            store.names().face(&band_name(&TALL_DOOR, tag)).is_none(),
            "band {tag} must not exist"
        );
    }
}

/// R2 acceptance 3: cascade after the door. Door THEN window (the window
/// lies inside the door's high-x split fragment) is watertight, the hole
/// lands on the containing fragment, and the cut order does not change the
/// resolved face geometry (window-then-door transfers the punched hole onto
/// the split fragment).
#[test]
fn door_then_window_cascade_is_order_independent() {
    // Face names shared by both orders (rim EDGE names differ by design:
    // a rim snapshots the punched face's name at punch time).
    let mut names: Vec<FaceName> = vec![wall_tag_name("end-east"), wall_tag_name("end-west")];
    for tag in ["outer", "inner"] {
        for side in [SplitSide::Left, SplitSide::Right] {
            names.push(split_name(TALL_DOOR.op, wall_tag_name(tag), side));
        }
    }
    for role in [FaceRole::CapStart, FaceRole::CapEnd] {
        for side in [SplitSide::Left, SplitSide::Right] {
            names.push(split_name(TALL_DOOR.op, cap_name(role.clone()), side));
        }
    }
    for tag in ["jamb-left", "jamb-right"] {
        names.push(band_name(&TALL_DOOR, tag));
    }
    for tag in ["sill", "jamb-right", "head", "jamb-left"] {
        names.push(band_name(&WINDOW, tag));
    }

    let (store_a, result_a) = cascade(&[TALL_DOOR, WINDOW]);
    let (store_b, result_b) = cascade(&[WINDOW, TALL_DOOR]);
    for (store, result, label) in [
        (&store_a, result_a, "door→win"),
        (&store_b, result_b, "win→door"),
    ] {
        assert_eq!(
            welded_boundary_edges(store, result),
            0,
            "{label} cascade must weld watertight"
        );
        let mesh = TessellateSolid::new(result, TessellationParams::default())
            .execute(store)
            .unwrap();
        assert_open(&mesh, &TALL_DOOR);
        assert_open(&mesh, &WINDOW);
    }

    let reference_samples: HashMap<String, Vec<Point3>> = names
        .iter()
        .map(|n| (format!("{n:?}"), face_sample(&store_a, n)))
        .collect();
    for name in &names {
        let sample_b = face_sample(&store_b, name);
        let reference = &reference_samples[&format!("{name:?}")];
        for (pa, pb) in reference.iter().zip(&sample_b) {
            assert!(
                (*pa - *pb).norm() < 1e-9,
                "{name:?} geometry differs across cut orders: {pa:?} vs {pb:?}"
            );
        }
    }

    // The window hole rides on the high-x outer / inner fragments in BOTH
    // orders (door-then-window punches the fragment; window-then-door
    // transfers the hole when the door splits the punched face).
    for store in [&store_a, &store_b] {
        for tag in ["outer", "inner"] {
            let mut hole_fragments = 0;
            for side in [SplitSide::Left, SplitSide::Right] {
                let name = split_name(TALL_DOOR.op, wall_tag_name(tag), side);
                let f = store.names().face(&name).unwrap();
                let face = store.face(f).unwrap();
                if !face.inner_wires.is_empty() {
                    hole_fragments += 1;
                    assert_eq!(face.inner_wires.len(), 1, "{name:?} hole count");
                    // The window lies at x in [2.2, 3.4]: the punched
                    // fragment is the high-x one — every outer-wire vertex
                    // sits at or beyond the doorway's right jamb.
                    let wire = store.wire(face.outer_wire).unwrap();
                    for oe in &wire.edges {
                        let edge = store.edge(oe.edge).unwrap();
                        for v in [edge.start, edge.end] {
                            let p = store.vertex(v).unwrap().point;
                            assert!(
                                p.x >= TALL_DOOR.x1 - 1e-9,
                                "window hole landed on the wrong fragment \
                                 (wire vertex at x = {})",
                                p.x
                            );
                        }
                    }
                }
            }
            assert_eq!(hole_fragments, 1, "{tag}: exactly one fragment punched");
        }
    }
}

/// R2 acceptance 4 (trimmed-outer perimeter generalization): a SECOND
/// bottom-touching door notches the already-notched wall faces and
/// re-splits an already-split cap fragment — the F5 rectangular-outer
/// restriction is gone.
#[test]
fn door_then_door_notches_notched_faces() {
    let (store, result) = cascade(&[DOOR, DOOR_B]);
    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();

    assert_eq!(
        welded_boundary_edges(&store, result),
        0,
        "two-door wall must weld watertight"
    );
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &DOOR);
    assert_open(&mesh, &DOOR_B);

    // The wall side faces are notched twice — their tagged names transfer
    // through both cuts.
    for tag in ["outer", "inner"] {
        let f = store
            .names()
            .face(&wall_tag_name(tag))
            .unwrap_or_else(|| panic!("{tag} does not resolve after two doors"));
        assert!(shell.faces.contains(&f), "{tag} outside the shell");
    }

    // The bottom cap splits twice: door A leaves Split{A, CapStart, L/R};
    // door B re-splits the fragment containing x in [3, 4], nesting its
    // name. Exactly one first-generation fragment survives, and the other
    // one's two second-generation fragments resolve.
    let first_gen: Vec<FaceName> = [SplitSide::Left, SplitSide::Right]
        .iter()
        .map(|&side| split_name(DOOR.op, cap_name(FaceRole::CapStart), side))
        .collect();
    let surviving: Vec<&FaceName> = first_gen
        .iter()
        .filter(|n| store.names().face(n).is_some())
        .collect();
    assert_eq!(
        surviving.len(),
        1,
        "exactly one first-generation cap fragment survives door B"
    );
    let renotched = first_gen
        .iter()
        .find(|n| store.names().face(n).is_none())
        .unwrap();
    for side in [SplitSide::Left, SplitSide::Right] {
        let nested = split_name(DOOR_B.op, renotched.clone(), side);
        assert!(
            store.names().face(&nested).is_some(),
            "{nested:?} does not resolve"
        );
    }

    // Three planar bottom-cap fragments at z = 0 in the final shell.
    let bottom_caps = shell
        .faces
        .iter()
        .filter(|&&f| {
            let face = store.face(f).unwrap();
            match &face.surface {
                FaceSurface::Plane(plane) => plane.origin().z.abs() < 1e-9,
                _ => false,
            }
        })
        .count();
    assert_eq!(bottom_caps, 3, "two doors leave three bottom-cap fragments");

    // Both doors' jamb / head bands resolve.
    for cut in [&DOOR, &DOOR_B] {
        for tag in ["jamb-left", "head", "jamb-right"] {
            assert!(
                store.names().face(&band_name(cut, tag)).is_some(),
                "{} band {tag} does not resolve",
                cut.op
            );
        }
    }
}

/// R2 acceptance 3b (demo variant geometry): a bottom-touching door THEN an
/// interior window — the window punches the already-NOTCHED wall faces (a
/// trim hole on a single kept fragment, filtered by its trim region).
#[test]
fn door_then_window_punches_notched_faces() {
    let (store, result) = cascade(&[DOOR, WINDOW]);
    assert_eq!(
        welded_boundary_edges(&store, result),
        0,
        "door + window wall must weld watertight"
    );
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &DOOR);
    assert_open(&mesh, &WINDOW);

    // The notched wall faces keep their names and carry the window hole.
    for tag in ["outer", "inner"] {
        let f = store
            .names()
            .face(&wall_tag_name(tag))
            .unwrap_or_else(|| panic!("{tag} does not resolve"));
        assert_eq!(
            store.face(f).unwrap().inner_wires.len(),
            1,
            "{tag} carries the window hole"
        );
    }
    // Both cuts' bands resolve (no sill band for the door; all four for
    // the window).
    for tag in ["jamb-left", "head", "jamb-right"] {
        assert!(store.names().face(&band_name(&DOOR, tag)).is_some());
    }
    for tag in ["sill", "jamb-right", "head", "jamb-left"] {
        assert!(store.names().face(&band_name(&WINDOW, tag)).is_some());
    }
}

/// Open-chain extraction contract: a cap-touching door yields one
/// multi-face through cut whose two chains are OPEN, single-target,
/// direction-aligned, and terminal-pinned exactly on the bottom ring
/// (z = 0).
#[test]
fn door_branches_chain_into_open_through_loops() {
    use super::loops::{collect_nurbs_faces, extract_cut_loops, ToolFaceCut};

    let mut store = TopologyStore::new();
    let wall = MakeSegmentedPrism::new(wall_profile(), Vector3::new(0.0, 0.0, 3.0))
        .execute(&mut store)
        .unwrap();
    let cutter = build_cutter(&mut store, &DOOR);
    let faces = |solid: SolidId| {
        store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap()
            .faces
            .clone()
    };
    let target = collect_nurbs_faces(&store, &faces(wall));
    let tool = collect_nurbs_faces(&store, &faces(cutter));
    let cuts = extract_cut_loops(&target, &tool).unwrap();
    assert_eq!(cuts.len(), 1, "one multi-face through cut");
    let ToolFaceCut::MultiFaceThrough { chains } = &cuts[0] else {
        panic!("expected a multi-face through cut, got {:?}", cuts[0]);
    };
    assert!(chains[0].mean_v() < chains[1].mean_v());
    for chain in chains {
        assert!(!chain.closed, "cap-touching chains must be open");
        assert_eq!(chain.segments.len(), 3, "jamb-left + head + jamb-right");
        assert!(
            chain.single_target_face().is_some(),
            "each open chain stays on one wall face"
        );
        // Terminals pinned exactly on the bottom ring.
        let head = chain.segments[0].branch.points[0];
        let tail = *chain.segments.last().unwrap().branch.points.last().unwrap();
        assert!(head.z == 0.0, "head terminal not pinned to z = 0");
        assert!(tail.z == 0.0, "tail terminal not pinned to z = 0");
        // Deterministic orientation: lexicographically smaller terminal
        // first.
        assert!(head.x < tail.x, "open chain not normalized");
    }
}
