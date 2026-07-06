//! F6 Phase R3b acceptance: curved-wall opening coverage.
//!
//! Every fixture is a plan-arc `MakeSegmentedPrism` wall whose profile
//! includes exact rational `Arc` segments (an annular strip: outer arc +
//! end lines + inner arc, tagged, with op ids), cut by 4-side-face box
//! cutters extruded through the curved faces — the exact scenarios
//! revion's curved walls with windows and doors produce (item 4).
//!
//! ## Verdicts (R3b)
//!
//! | Case | Verdict |
//! |---|---|
//! | Window fully inside the arc face | **CLEAN** — the F5 B4 result carried through the R1-style named cascade fixture; bands, rims and tags resolve, rebuild-stable |
//! | Window straddling the arc↔line tangent joint | **CLEAN** — the R2/C split machinery is geometry-agnostic in UV: both joint-adjacent faces take boundary-notch hole halves, the split joint sub-edges are shared, tags transfer |
//! | Full-height (cap-touching) door through the arc face | **CLEAN** — the R2 open-chain + cap-notch path works on curved faces; the annular bottom cap splits into two fragments whose notch rides on ARC sub-span edges shared with the curved wall faces |
//! | Door + two windows cascade on the arc wall | **CLEAN** — order-independent (1e-9 sampled geometry) with a stable resolving name set |

use std::collections::HashMap;
use std::f64::consts::PI;

use crate::math::{Point3, Vector3};
use crate::operations::boolean::Subtract;
use crate::operations::creation::{MakeSegmentedPrism, ProfileSegment};
use crate::tessellation::{TessellateSolid, TessellationParams};
use crate::topology::{
    EdgeId, EdgeName, FaceName, FaceRole, FaceSurface, OpId, SegmentTag, SolidId, SplitSide,
    TopologyStore,
};

use super::test_support::welded_boundary_edges;

const R_OUTER: f64 = 8.4;
const R_INNER: f64 = 8.0;
const HEIGHT: f64 = 3.0;

fn deg(d: f64) -> f64 {
    d * PI / 180.0
}

/// Annular plan-arc wall strip: outer arc r = 8.4, inner arc r = 8.0,
/// azimuth 60..120 degrees, extruded 3.0 up. The inner arc is traversed
/// backwards via the -Z normal (the Phase A fillet convention).
fn arc_wall_profile() -> Vec<ProfileSegment> {
    let outer_start = Point3::new(R_OUTER * deg(60.0).cos(), R_OUTER * deg(60.0).sin(), 0.0);
    let outer_end = Point3::new(R_OUTER * deg(120.0).cos(), R_OUTER * deg(120.0).sin(), 0.0);
    let inner_start = Point3::new(R_INNER * deg(120.0).cos(), R_INNER * deg(120.0).sin(), 0.0);
    let inner_end = Point3::new(R_INNER * deg(60.0).cos(), R_INNER * deg(60.0).sin(), 0.0);
    vec![
        ProfileSegment::Arc {
            center: Point3::origin(),
            radius: R_OUTER,
            normal: Vector3::z(),
            ref_dir: Vector3::x(),
            start_angle: deg(60.0),
            end_angle: deg(120.0),
        },
        ProfileSegment::Line {
            start: outer_end,
            end: inner_start,
        },
        ProfileSegment::Arc {
            center: Point3::origin(),
            radius: R_INNER,
            normal: -Vector3::z(),
            ref_dir: Vector3::x(),
            start_angle: deg(-120.0),
            end_angle: deg(-60.0),
        },
        ProfileSegment::Line {
            start: inner_end,
            end: outer_start,
        },
    ]
}

fn arc_wall_tags() -> Vec<SegmentTag> {
    ["convex", "end-west", "concave", "end-east"]
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

/// One opening cutter: a box profile in the XZ plane at `y = 6`, spanning
/// `x in [x0, x1]`, `z in [z0, z1]`, extruded 3.5 along `+Y` — through the
/// annular wall (crossing the concave face first, then the convex face)
/// with margin on both sides. `z0 < 0` makes it cap-touching.
#[derive(Clone, Copy)]
struct Opening {
    op: &'static str,
    x0: f64,
    x1: f64,
    z0: f64,
    z1: f64,
}

/// Interior window centered on the arc apex (fully inside the arc faces).
const WINDOW: Opening = Opening {
    op: "cut-win",
    x0: -0.7,
    x1: 0.7,
    z0: 0.9,
    z1: 1.7,
};

/// Full-height door: sill below the wall, head inside it.
const DOOR: Opening = Opening {
    op: "cut-door",
    x0: -0.45,
    x1: 0.45,
    z0: -0.5,
    z1: 2.25,
};

/// Interior window west of the door.
const WIN_A: Opening = Opening {
    op: "cut-winA",
    x0: -2.6,
    x1: -1.4,
    z0: 1.0,
    z1: 2.0,
};

/// Interior window east of the door.
const WIN_B: Opening = Opening {
    op: "cut-winB",
    x0: 1.4,
    x1: 2.6,
    z0: 1.0,
    z1: 2.0,
};

/// The tool op id of an opening cut (`cut-door` -> `tool-cut-door`).
fn tool_op(cut: &Opening) -> String {
    format!("tool-{}", cut.op)
}

fn build_cutter(store: &mut TopologyStore, cut: &Opening) -> SolidId {
    let q = |x: f64, z: f64| Point3::new(x, 6.0, z);
    let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
    let profile = vec![
        line(q(cut.x0, cut.z0), q(cut.x1, cut.z0)), // sill
        line(q(cut.x1, cut.z0), q(cut.x1, cut.z1)), // jamb-right
        line(q(cut.x1, cut.z1), q(cut.x0, cut.z1)), // head
        line(q(cut.x0, cut.z1), q(cut.x0, cut.z0)), // jamb-left
    ];
    MakeSegmentedPrism::new(profile, Vector3::new(0.0, 3.5, 0.0))
        .with_op_id(OpId::new(tool_op(cut)))
        .with_segment_tags(box_tags())
        .execute(store)
        .unwrap()
}

/// Builds (arc wall) − openings applied SEQUENTIALLY in the given order;
/// each `Subtract` result feeds the next as its target.
fn cascade(order: &[Opening]) -> (TopologyStore, SolidId) {
    let mut store = TopologyStore::new();
    let wall = MakeSegmentedPrism::new(arc_wall_profile(), Vector3::new(0.0, 0.0, HEIGHT))
        .with_op_id(OpId::new("wall1"))
        .with_segment_tags(arc_wall_tags())
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

/// Entry rim (concave face, crossed first along `+Y`) and exit rim
/// (convex face) of an interior opening cut.
fn rim_names(cut: &Opening) -> [EdgeName; 2] {
    [("concave", 0u32), ("convex", 1u32)].map(|(target, loop_index)| EdgeName::CutRim {
        op: OpId::new(cut.op),
        target: Box::new(wall_tag_name(target)),
        loop_index,
    })
}

/// A geometric sample of a resolved named face: NURBS faces sample the
/// surface at fixed UVs; planar faces (cap fragments) average their outer
/// wire's vertices — which distinguishes the Left / Right fragments
/// sharing one plane.
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

/// Asserts no mesh vertex intrudes into the opening tunnel through the
/// annular wall (radial probe between the arc faces; margins shrink the
/// box so rim vertices do not trip it).
fn assert_open(mesh: &crate::tessellation::TriangleMesh, cut: &Opening) {
    let (x0, x1) = (cut.x0 + 0.05, cut.x1 - 0.05);
    let (z0, z1) = (cut.z0.max(0.0) + 0.05, cut.z1.min(HEIGHT) - 0.05);
    for v in &mesh.vertices {
        let r = (v.x * v.x + v.y * v.y).sqrt();
        let inside = v.x > x0 && v.x < x1 && v.z > z0 && v.z < z1;
        assert!(
            !(inside && r > R_INNER + 0.05 && r < R_OUTER - 0.05),
            "vertex ({:.3},{:.3},{:.3}) intrudes into the {} tunnel",
            v.x,
            v.y,
            v.z,
            cut.op
        );
    }
}

/// The outer-wire edge ids of a face.
fn wire_edges(store: &TopologyStore, f: crate::topology::FaceId) -> Vec<EdgeId> {
    store
        .wire(store.face(f).unwrap().outer_wire)
        .unwrap()
        .edges
        .iter()
        .map(|oe| oe.edge)
        .collect()
}

/// R3b case 1: a window fully inside the curved (arc) faces, through the
/// R1-style named cascade fixture — watertight, all band / tag / rim names
/// resolve into the result shell, and a from-scratch rebuild resolves every
/// name to geometrically identical entities.
#[test]
fn arc_window_names_resolve_and_are_rebuild_stable() {
    let (store, result) = cascade(&[WINDOW]);
    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();
    // 6 wall faces (4 sides + 2 caps) + 4 band fragments.
    assert_eq!(shell.faces.len(), 10, "got {}", shell.faces.len());

    let boundary = welded_boundary_edges(&store, result);
    assert_eq!(
        boundary, 0,
        "arc-window wall must position-weld watertight ({boundary})"
    );
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &WINDOW);

    // Both punched faces are the CURVED arc faces, carrying one hole each.
    for tag in ["convex", "concave"] {
        let f = store
            .names()
            .face(&wall_tag_name(tag))
            .unwrap_or_else(|| panic!("{tag} does not resolve"));
        assert!(shell.faces.contains(&f), "{tag} outside the result shell");
        assert_eq!(
            store.face(f).unwrap().inner_wires.len(),
            1,
            "{tag} carries the window hole"
        );
    }
    for tag in ["end-west", "end-east"] {
        assert!(store.names().face(&wall_tag_name(tag)).is_some());
    }
    for tag in ["sill", "jamb-right", "head", "jamb-left"] {
        let f = store
            .names()
            .face(&band_name(&WINDOW, tag))
            .unwrap_or_else(|| panic!("band {tag} does not resolve"));
        assert!(shell.faces.contains(&f), "band {tag} outside the shell");
    }
    for rim in rim_names(&WINDOW) {
        assert!(store.names().edge(&rim).is_some(), "{rim:?} unresolved");
    }

    // Rebuild stability: every persistent face name resolves to identical
    // geometry in a from-scratch rebuild.
    let (store_b, _) = cascade(&[WINDOW]);
    let mut names: Vec<FaceName> = ["convex", "end-west", "concave", "end-east"]
        .iter()
        .map(|tag| wall_tag_name(tag))
        .collect();
    for tag in ["sill", "jamb-right", "head", "jamb-left"] {
        names.push(band_name(&WINDOW, tag));
    }
    for name in &names {
        let a = face_sample(&store, name);
        let b = face_sample(&store_b, name);
        for (pa, pb) in a.iter().zip(&b) {
            assert!(
                (*pa - *pb).norm() < 1e-9,
                "{name:?} moved across rebuilds: {pa:?} vs {pb:?}"
            );
        }
    }
}

// ---- Case 2: window straddling the arc <-> line tangent joint ----

/// Plan wall whose outer / inner boundaries are an arc (azimuth 90..150,
/// r = 8.4 / 8.0) continuing TANGENTIALLY into a straight run (x in [0, 3]
/// at y = 8.4 / 8.0). The arc and line side faces share a vertical joint
/// edge at x = 0 — the curved analogue of the F5 Phase C collinear kink.
fn tangent_wall_profile() -> Vec<ProfileSegment> {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let outer_150 = Point3::new(R_OUTER * deg(150.0).cos(), R_OUTER * deg(150.0).sin(), 0.0);
    let inner_150 = Point3::new(R_INNER * deg(150.0).cos(), R_INNER * deg(150.0).sin(), 0.0);
    vec![
        ProfileSegment::Line {
            start: p(3.0, R_OUTER),
            end: p(0.0, R_OUTER),
        },
        ProfileSegment::Arc {
            center: Point3::origin(),
            radius: R_OUTER,
            normal: Vector3::z(),
            ref_dir: Vector3::x(),
            start_angle: deg(90.0),
            end_angle: deg(150.0),
        },
        ProfileSegment::Line {
            start: outer_150,
            end: inner_150,
        },
        ProfileSegment::Arc {
            center: Point3::origin(),
            radius: R_INNER,
            normal: -Vector3::z(),
            ref_dir: Vector3::x(),
            start_angle: deg(-150.0),
            end_angle: deg(-90.0),
        },
        ProfileSegment::Line {
            start: p(0.0, R_INNER),
            end: p(3.0, R_INNER),
        },
        ProfileSegment::Line {
            start: p(3.0, R_INNER),
            end: p(3.0, R_OUTER),
        },
    ]
}

fn tangent_wall_tags() -> Vec<SegmentTag> {
    [
        "outer-line",
        "outer-arc",
        "end-west",
        "inner-arc",
        "inner-line",
        "end-east",
    ]
    .iter()
    .map(|t| SegmentTag::new(*t))
    .collect()
}

/// Window straddling the arc↔line joint at x = 0.
const STRADDLE: Opening = Opening {
    op: "cut-straddle",
    x0: -0.5,
    x1: 0.5,
    z0: 1.0,
    z1: 2.0,
};

/// Builds (tangent arc+line wall) − (window straddling the joint).
fn tangent_wall_minus_straddle() -> (TopologyStore, SolidId) {
    let mut store = TopologyStore::new();
    let wall = MakeSegmentedPrism::new(tangent_wall_profile(), Vector3::new(0.0, 0.0, HEIGHT))
        .with_op_id(OpId::new("wall1"))
        .with_segment_tags(tangent_wall_tags())
        .execute(&mut store)
        .unwrap();
    let cutter = build_cutter(&mut store, &STRADDLE);
    let result = Subtract::new(wall, cutter)
        .with_op_id(OpId::new(STRADDLE.op))
        .execute(&mut store)
        .unwrap_or_else(|e| panic!("straddle cut failed: {e:?}"));
    (store, result)
}

/// R3b case 2: a window STRADDLING the arc↔line tangent joint — both
/// joint-adjacent faces take their hole half as a boundary notch, keep
/// their tagged names, share the split joint sub-edges, and the result is
/// watertight with a genuinely open window.
#[test]
fn window_straddling_arc_line_joint_is_watertight() {
    let (store, result) = tangent_wall_minus_straddle();
    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();
    // 2 caps + 2 ends + 4 notched side faces + 4 band fragments.
    assert_eq!(shell.faces.len(), 12, "got {}", shell.faces.len());

    // The hole halves live in boundary notches: no interior hole wires.
    for &f in &shell.faces {
        assert!(
            store.face(f).unwrap().inner_wires.is_empty(),
            "straddling window must notch, not punch"
        );
    }

    let boundary = welded_boundary_edges(&store, result);
    assert_eq!(
        boundary, 0,
        "straddling window must position-weld watertight ({boundary})"
    );

    // The window is genuinely open (the wall here spans y in [~8.0, 8.4]
    // for |x| <= 0.5 on both the arc and line sides).
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    for v in &mesh.vertices {
        let inside = v.x > STRADDLE.x0 + 0.05
            && v.x < STRADDLE.x1 - 0.05
            && v.z > STRADDLE.z0 + 0.05
            && v.z < STRADDLE.z1 - 0.05;
        assert!(
            !(inside && v.y > R_INNER + 0.05 && v.y < R_OUTER - 0.05),
            "vertex ({:.3},{:.3},{:.3}) intrudes into the straddle tunnel",
            v.x,
            v.y,
            v.z
        );
    }

    // All 4 notched side faces keep their tagged names (boundary notch →
    // transfer) and resolve into the result shell.
    for tag in ["outer-line", "outer-arc", "inner-arc", "inner-line"] {
        let f = store
            .names()
            .face(&wall_tag_name(tag))
            .unwrap_or_else(|| panic!("{tag} does not resolve"));
        assert!(shell.faces.contains(&f), "{tag} outside the result shell");
    }
    // All 4 band fragments resolve.
    for tag in ["sill", "jamb-right", "head", "jamb-left"] {
        assert!(
            store.names().face(&band_name(&STRADDLE, tag)).is_some(),
            "band {tag} does not resolve"
        );
    }

    // The split joint sub-edges (below / above the window, on the outer
    // joint at (0, 8.4) and the inner joint at (0, 8.0)) are each shared
    // by exactly the two adjacent fragments.
    let mut edge_face_uses: HashMap<EdgeId, usize> = HashMap::new();
    for &f in &shell.faces {
        for e in wire_edges(&store, f) {
            *edge_face_uses.entry(e).or_insert(0) += 1;
        }
    }
    let mut joint_subs = 0usize;
    for (&edge, &uses) in &edge_face_uses {
        let data = store.edge(edge).unwrap();
        let a = store.vertex(data.start).unwrap().point;
        let b = store.vertex(data.end).unwrap().point;
        let on_joint = a.x.abs() < 1e-9
            && b.x.abs() < 1e-9
            && ((a.y - R_OUTER).abs() < 1e-9 && (b.y - R_OUTER).abs() < 1e-9
                || (a.y - R_INNER).abs() < 1e-9 && (b.y - R_INNER).abs() < 1e-9)
            && (a.z - b.z).abs() > 1e-9;
        if on_joint {
            joint_subs += 1;
            assert_eq!(uses, 2, "joint sub-edge must be shared by 2 fragments");
        }
    }
    assert_eq!(
        joint_subs, 4,
        "below + above sub-edges on the outer and inner joints"
    );
}

/// R3b case 2b: the straddling cut re-resolves identically across two
/// from-scratch builds.
#[test]
fn window_straddling_arc_line_joint_is_rebuild_stable() {
    let (store_a, _) = tangent_wall_minus_straddle();
    let (store_b, _) = tangent_wall_minus_straddle();
    let mut names: Vec<FaceName> = ["outer-line", "outer-arc", "inner-arc", "inner-line"]
        .iter()
        .map(|tag| wall_tag_name(tag))
        .collect();
    for tag in ["sill", "jamb-right", "head", "jamb-left"] {
        names.push(band_name(&STRADDLE, tag));
    }
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

/// R3b case 3: a full-height (cap-touching) door through the ARC face —
/// the R2 open-chain + cap-notch path on a curved face. The annular bottom
/// cap splits into two named fragments whose notch edges include ARC
/// sub-spans shared with the curved wall faces; watertight, open doorway,
/// all names resolve.
#[test]
fn full_height_door_through_arc_face_notches_bottom_cap() {
    let (store, result) = cascade(&[DOOR]);
    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();
    // convex + concave notched faces, 2 ends, top cap, 2 bottom cap
    // fragments, 3 band fragments (jamb-left / head / jamb-right).
    assert_eq!(shell.faces.len(), 10, "got {}", shell.faces.len());

    let boundary = welded_boundary_edges(&store, result);
    assert_eq!(
        boundary, 0,
        "curved door wall must weld watertight ({boundary})"
    );
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &DOOR);

    // Wall side / end faces keep their tagged names (notch → transfer).
    for tag in ["convex", "end-west", "concave", "end-east"] {
        let f = store
            .names()
            .face(&wall_tag_name(tag))
            .unwrap_or_else(|| panic!("{tag} does not resolve"));
        assert!(shell.faces.contains(&f), "{tag} outside the result shell");
    }
    // Top cap survives; the bottom cap split into Left / Right planar
    // fragments.
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
    // Jamb / head bands resolve; NO sill band (the sill lies below the
    // wall).
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
    // Entry / exit rims resolve.
    for rim in rim_names(&DOOR) {
        assert!(store.names().edge(&rim).is_some(), "{rim:?} unresolved");
    }

    // Shared-edge topology on the CURVED boundary: each bottom-cap
    // fragment's outer wire contains an ARC sub-span — an edge whose
    // endpoints both lie on one of the wall arcs (r = 8.0 or 8.4) at
    // z = 0 and which is SHARED (same EdgeId) with the notched curved
    // wall face — plus a cap-plane closure edge shared with a jamb band.
    let convex_edges = wire_edges(
        &store,
        store.names().face(&wall_tag_name("convex")).unwrap(),
    );
    let concave_edges = wire_edges(
        &store,
        store.names().face(&wall_tag_name("concave")).unwrap(),
    );
    let band_edges: Vec<EdgeId> = ["jamb-left", "jamb-right"]
        .iter()
        .flat_map(|tag| wire_edges(&store, store.names().face(&band_name(&DOOR, tag)).unwrap()))
        .collect();
    for side in [SplitSide::Left, SplitSide::Right] {
        let name = split_name(DOOR.op, cap_name(FaceRole::CapStart), side);
        let cap_edges = wire_edges(&store, store.names().face(&name).unwrap());
        let arc_subs = cap_edges
            .iter()
            .filter(|&&e| {
                let data = store.edge(e).unwrap();
                let a = store.vertex(data.start).unwrap().point;
                let b = store.vertex(data.end).unwrap().point;
                let ra = (a.x * a.x + a.y * a.y).sqrt();
                let rb = (b.x * b.x + b.y * b.y).sqrt();
                let on_arc = ((ra - R_OUTER).abs() < 1e-9 && (rb - R_OUTER).abs() < 1e-9)
                    || ((ra - R_INNER).abs() < 1e-9 && (rb - R_INNER).abs() < 1e-9);
                on_arc
                    && a.z.abs() < 1e-9
                    && b.z.abs() < 1e-9
                    && (convex_edges.contains(&e) || concave_edges.contains(&e))
            })
            .count();
        assert!(
            arc_subs >= 2,
            "{name:?} must ride on arc sub-spans shared with both curved \
             wall faces (found {arc_subs})"
        );
        assert!(
            cap_edges.iter().any(|e| band_edges.contains(e)),
            "{name:?} shares no closure edge with the jamb bands"
        );
    }
}

/// R3b case 4: cascade on the curved wall — full-height door + two
/// windows subtracted sequentially. Watertight with genuinely open
/// openings in every cut order, geometry-equal per resolved name (1e-9
/// sampled), and the resolving name set is stable across orders.
#[test]
fn curved_cascade_is_order_independent() {
    let orders: [[Opening; 3]; 3] = [
        [DOOR, WIN_A, WIN_B],
        [WIN_B, WIN_A, DOOR],
        [WIN_A, DOOR, WIN_B],
    ];

    // Face names shared by all orders (rim EDGE names are cascade-order-
    // dependent by design: a rim snapshots the punched face's name).
    let mut names: Vec<FaceName> = ["convex", "end-west", "concave", "end-east"]
        .iter()
        .map(|tag| wall_tag_name(tag))
        .collect();
    names.push(cap_name(FaceRole::CapEnd));
    for side in [SplitSide::Left, SplitSide::Right] {
        names.push(split_name(DOOR.op, cap_name(FaceRole::CapStart), side));
    }
    for tag in ["jamb-left", "head", "jamb-right"] {
        names.push(band_name(&DOOR, tag));
    }
    for win in [&WIN_A, &WIN_B] {
        for tag in ["sill", "jamb-right", "head", "jamb-left"] {
            names.push(band_name(win, tag));
        }
    }

    let (ref_store, ref_result) = cascade(&orders[0]);
    let ref_shell_len = ref_store
        .shell(ref_store.solid(ref_result).unwrap().outer_shell)
        .unwrap()
        .faces
        .len();
    // 4 sides + top cap + 2 bottom fragments + 3 door bands + 8 window
    // bands.
    assert_eq!(ref_shell_len, 18, "got {ref_shell_len}");
    let ref_samples: HashMap<String, Vec<Point3>> = names
        .iter()
        .map(|n| (format!("{n:?}"), face_sample(&ref_store, n)))
        .collect();

    for order in &orders {
        let (store, result) = cascade(order);
        let shell = store
            .shell(store.solid(result).unwrap().outer_shell)
            .unwrap();
        assert_eq!(shell.faces.len(), ref_shell_len, "face count differs");
        assert_eq!(
            welded_boundary_edges(&store, result),
            0,
            "curved cascade must weld watertight"
        );
        let mesh = TessellateSolid::new(result, TessellationParams::default())
            .execute(&store)
            .unwrap();
        for cut in order {
            assert_open(&mesh, cut);
        }
        // Name-set stability: every shared name resolves into the shell,
        // the retired CapStart resolves in none.
        assert!(store.names().face(&cap_name(FaceRole::CapStart)).is_none());
        for name in &names {
            let f = store.names().face(name).unwrap_or_else(|| {
                panic!(
                    "{name:?} does not resolve in {order:?}",
                    order = order.map(|c| c.op)
                )
            });
            assert!(
                shell.faces.contains(&f),
                "{name:?} resolves outside the result shell"
            );
            let samples = face_sample(&store, name);
            let reference = &ref_samples[&format!("{name:?}")];
            for (s, r) in samples.iter().zip(reference) {
                assert!(
                    (*s - *r).norm() < 1e-9,
                    "{name:?} geometry differs across cut orders: {s:?} vs {r:?}"
                );
            }
        }
    }
}
