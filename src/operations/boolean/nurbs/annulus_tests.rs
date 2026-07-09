//! Annulus (hole-capable segmented prism) boolean acceptance.
//!
//! A closed room wall is ONE segmented prism whose footprint is an outer
//! ring plus a courtyard hole. Openings cut through the wall thickness
//! enter an OUTER ring face and exit the INNER (hole ring) face — two
//! NURBS target faces, exactly the multi-face through-cut class — and the
//! caps now carry inner wires, exercising the annulus-aware cap-notch
//! rebuild:
//!
//! - a window punches both ring faces and leaves the caps' courtyard
//!   wires untouched;
//! - a bottom-touching door notches BOTH cap wires — the courtyard hole
//!   merges into the bottom cap's outer boundary (single fragment, no
//!   inner wire, name transfers);
//! - a floor-to-ceiling door severs one wall segment's faces into
//!   `Split { Left / Right }` while BOTH caps stay single C-shaped
//!   fragments (an annulus stays connected after one radial notch);
//! - with TWO courtyards, a door into courtyard A merges A into the cap
//!   boundary while courtyard B's wire rides along untouched as the
//!   fragment's inner wire.

use crate::math::{Point3, Vector3};
use crate::operations::boolean::Subtract;
use crate::operations::creation::{MakeSegmentedPrism, ProfileSegment};
use crate::tessellation::{TessellateSolid, TessellationParams};
use crate::topology::{FaceName, FaceRole, OpId, SegmentTag, SolidId, SplitSide, TopologyStore};

use super::test_support::welded_boundary_edges;

const HEIGHT: f64 = 3.0;

/// Axis-aligned rectangle chain in the z = 0 plane. `ccw = true` winds
/// counter-clockwise about `+Z`, `false` clockwise.
fn rect_ring(x0: f64, y0: f64, x1: f64, y1: f64, ccw: bool) -> Vec<ProfileSegment> {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let corners = if ccw {
        [p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)]
    } else {
        [p(x0, y0), p(x0, y1), p(x1, y1), p(x1, y0)]
    };
    (0..4)
        .map(|i| ProfileSegment::Line {
            start: corners[i],
            end: corners[(i + 1) % 4],
        })
        .collect()
}

fn tags(names: [&str; 4]) -> Vec<SegmentTag> {
    names.iter().map(|t| SegmentTag::new(*t)).collect()
}

/// A closed square room: outer 6 × 6, courtyard 5 × 5 (wall thickness
/// 0.5), one annulus segmented prism. The hole winds clockwise (the
/// natural annulus footprint convention), and both rings are tagged.
fn room_wall(store: &mut TopologyStore) -> SolidId {
    MakeSegmentedPrism::new(
        rect_ring(0.0, 0.0, 6.0, 6.0, true),
        Vector3::new(0.0, 0.0, HEIGHT),
    )
    .with_holes(vec![rect_ring(0.5, 0.5, 5.5, 5.5, false)])
    .with_op_id(OpId::new("room1"))
    .with_segment_tags(tags(["outer-s", "outer-e", "outer-n", "outer-w"]))
    .with_hole_tags(vec![tags(["inner-w", "inner-n", "inner-e", "inner-s"])])
    .execute(store)
    .unwrap()
}

/// A duplex: outer 9 × 6 footprint with TWO courtyards (A west, B east)
/// separated by a partition wall.
fn duplex_wall(store: &mut TopologyStore) -> SolidId {
    MakeSegmentedPrism::new(
        rect_ring(0.0, 0.0, 9.0, 6.0, true),
        Vector3::new(0.0, 0.0, HEIGHT),
    )
    .with_holes(vec![
        rect_ring(0.5, 0.5, 4.0, 5.5, false),
        rect_ring(5.0, 0.5, 8.5, 5.5, false),
    ])
    .with_op_id(OpId::new("duplex1"))
    .with_segment_tags(tags(["outer-s", "outer-e", "outer-n", "outer-w"]))
    .with_hole_tags(vec![
        tags(["a-w", "a-n", "a-e", "a-s"]),
        tags(["b-w", "b-n", "b-e", "b-s"]),
    ])
    .execute(store)
    .unwrap()
}

/// One opening cutter through the SOUTH wall (`y in [0, 0.5]`): profile in
/// the vertical XZ plane at `y = -1`, extruded 2.5 along `+Y` so the far
/// cap ends inside the courtyard air. `z0 < 0` makes it cap-touching.
#[derive(Clone, Copy)]
struct Opening {
    op: &'static str,
    x0: f64,
    x1: f64,
    z0: f64,
    z1: f64,
}

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
    MakeSegmentedPrism::new(profile, Vector3::new(0.0, 2.5, 0.0))
        .with_op_id(OpId::new(tool_op(cut)))
        .with_segment_tags(tags(["sill", "jamb-right", "head", "jamb-left"]))
        .execute(store)
        .unwrap()
}

fn subtract(store: &mut TopologyStore, wall: SolidId, cut: &Opening) -> SolidId {
    let cutter = build_cutter(store, cut);
    Subtract::new(wall, cutter)
        .with_op_id(OpId::new(cut.op))
        .execute(store)
        .unwrap_or_else(|e| panic!("{} failed: {e:?}", cut.op))
}

fn created(op: &str, role: FaceRole) -> FaceName {
    FaceName::Created {
        op: OpId::new(op),
        role,
    }
}

fn tag_name(op: &str, tag: &str) -> FaceName {
    created(op, FaceRole::Tagged(SegmentTag::new(tag)))
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

fn resolve(store: &TopologyStore, name: &FaceName) -> crate::topology::FaceId {
    store
        .names()
        .face(name)
        .unwrap_or_else(|| panic!("{name:?} does not resolve"))
}

/// Asserts no mesh vertex intrudes into the south-wall opening tunnel
/// (margins shrink the probe box so rim vertices do not trip it).
fn assert_open(mesh: &crate::tessellation::TriangleMesh, cut: &Opening) {
    let (x0, x1) = (cut.x0 + 0.05, cut.x1 - 0.05);
    let (z0, z1) = (cut.z0.max(0.0) + 0.05, cut.z1.min(HEIGHT) - 0.05);
    for v in &mesh.vertices {
        let inside = v.x > x0 && v.x < x1 && v.z > z0 && v.z < z1;
        assert!(
            !(inside && v.y > 0.05 && v.y < 0.45),
            "vertex ({:.3},{:.3},{:.3}) intrudes into the {} opening",
            v.x,
            v.y,
            v.z,
            cut.op
        );
    }
}

const WINDOW: Opening = Opening {
    op: "cut-win",
    x0: 2.0,
    x1: 3.5,
    z0: 1.0,
    z1: 2.0,
};

const DOOR: Opening = Opening {
    op: "cut-door",
    x0: 0.8,
    x1: 1.7,
    z0: -0.5,
    z1: 2.25,
};

const TALL_DOOR: Opening = Opening {
    op: "cut-door",
    x0: 0.8,
    x1: 1.7,
    z0: -0.5,
    z1: 3.5,
};

/// A window through-cuts the annulus wall: it enters the outer ring face
/// and exits the INNER ring face — two NURBS target faces — and the caps'
/// courtyard wires ride through untouched.
#[test]
fn window_through_annulus_wall_punches_outer_and_inner_ring_faces() {
    let mut store = TopologyStore::new();
    let wall = room_wall(&mut store);
    let result = subtract(&mut store, wall, &WINDOW);
    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();

    // 8 side faces + 2 caps + 4 band fragments.
    assert_eq!(shell.faces.len(), 14, "got {}", shell.faces.len());
    assert_eq!(
        welded_boundary_edges(&store, result),
        0,
        "windowed room must weld watertight"
    );
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &WINDOW);

    // Both punched ring faces keep their tagged names and carry the hole.
    for tag in ["outer-s", "inner-s"] {
        let face = resolve(&store, &tag_name("room1", tag));
        assert_eq!(
            store.face(face).unwrap().inner_wires.len(),
            1,
            "{tag} carries the window hole ring"
        );
    }
    // All four reveal bands resolve.
    for tag in ["sill", "jamb-right", "head", "jamb-left"] {
        assert!(
            store.names().face(&band_name(&WINDOW, tag)).is_some(),
            "band {tag} does not resolve"
        );
    }
    // The caps are unaffected copies: names carry over, courtyard intact.
    for role in [FaceRole::CapStart, FaceRole::CapEnd] {
        let cap = resolve(&store, &created("room1", role.clone()));
        assert_eq!(
            store.face(cap).unwrap().inner_wires.len(),
            1,
            "{role:?} must keep its courtyard wire"
        );
    }
}

/// A bottom-touching full-height door notches BOTH of the bottom cap's
/// wires: the courtyard hole merges into the cap's outer boundary — one
/// kept fragment, no inner wire, `CapStart` transfers. The top cap stays
/// an untouched annulus.
#[test]
fn full_height_door_merges_courtyard_into_bottom_cap_boundary() {
    let mut store = TopologyStore::new();
    let wall = room_wall(&mut store);
    let result = subtract(&mut store, wall, &DOOR);
    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();

    // 8 notched/copied side faces + top cap + 1 bottom cap fragment +
    // 3 bands (jamb-left / head / jamb-right; the sill lies below).
    assert_eq!(shell.faces.len(), 13, "got {}", shell.faces.len());
    assert_eq!(
        welded_boundary_edges(&store, result),
        0,
        "door room must weld watertight"
    );
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &DOOR);

    // The notched ring faces keep their tagged names (boundary notch).
    for tag in ["outer-s", "inner-s"] {
        assert!(
            store.names().face(&tag_name("room1", tag)).is_some(),
            "{tag} must keep its name after the notch"
        );
    }
    // Bottom cap: single kept fragment, name transfers, and the courtyard
    // wire has merged into the outer boundary (no inner wire left).
    let bottom = resolve(&store, &created("room1", FaceRole::CapStart));
    assert!(shell.faces.contains(&bottom), "CapStart outside the shell");
    assert!(
        store.face(bottom).unwrap().inner_wires.is_empty(),
        "the courtyard must merge into the notched bottom cap's boundary"
    );
    // Top cap: untouched annulus copy.
    let top = resolve(&store, &created("room1", FaceRole::CapEnd));
    assert_eq!(
        store.face(top).unwrap().inner_wires.len(),
        1,
        "CapEnd must keep its courtyard wire"
    );
    // Jamb / head bands resolve; no sill band for a cap-touching door.
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

    // Shared-edge topology: the rebuilt bottom cap shares kept sub-edges
    // with the notched outer face and closure edges with the jamb bands.
    let wire_edges = |f: crate::topology::FaceId| -> Vec<crate::topology::EdgeId> {
        store
            .wire(store.face(f).unwrap().outer_wire)
            .unwrap()
            .edges
            .iter()
            .map(|oe| oe.edge)
            .collect()
    };
    let cap_edges = wire_edges(bottom);
    let outer_edges = wire_edges(resolve(&store, &tag_name("room1", "outer-s")));
    let band_edges: Vec<crate::topology::EdgeId> = ["jamb-left", "jamb-right"]
        .iter()
        .flat_map(|tag| wire_edges(resolve(&store, &band_name(&DOOR, tag))))
        .collect();
    assert!(
        cap_edges.iter().any(|e| outer_edges.contains(e)),
        "bottom cap shares no kept sub-edge with the notched outer face"
    );
    assert!(
        cap_edges.iter().any(|e| band_edges.contains(e)),
        "bottom cap shares no closure edge with the jamb bands"
    );
}

/// A floor-to-ceiling door severs one wall segment: the outer and inner
/// ring faces split into `Left / Right`, and BOTH caps are notched — but
/// an annulus stays connected after one radial notch, so each cap stays a
/// single C-shaped fragment whose name TRANSFERS (no cap split).
#[test]
fn floor_to_ceiling_door_severs_annulus_wall_segment() {
    let mut store = TopologyStore::new();
    let wall = room_wall(&mut store);
    let result = subtract(&mut store, wall, &TALL_DOOR);
    let shell = store
        .shell(store.solid(result).unwrap().outer_shell)
        .unwrap();

    // outer-s L/R + inner-s L/R + 6 copied sides + 2 cap fragments +
    // 2 jamb bands.
    assert_eq!(shell.faces.len(), 14, "got {}", shell.faces.len());
    assert_eq!(
        welded_boundary_edges(&store, result),
        0,
        "severed room must weld watertight"
    );
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &TALL_DOOR);

    // The severed ring faces retire; Split L/R resolve for both rings.
    for tag in ["outer-s", "inner-s"] {
        assert!(
            store.names().face(&tag_name("room1", tag)).is_none(),
            "{tag} must retire when the door severs it"
        );
        for side in [SplitSide::Left, SplitSide::Right] {
            let name = FaceName::Split {
                op: OpId::new(TALL_DOOR.op),
                parent: Box::new(tag_name("room1", tag)),
                side,
            };
            assert!(
                store.names().face(&name).is_some(),
                "{name:?} does not resolve"
            );
        }
    }
    // Both caps stay single fragments: names transfer, courtyards merged.
    for role in [FaceRole::CapStart, FaceRole::CapEnd] {
        let cap = resolve(&store, &created("room1", role.clone()));
        assert!(shell.faces.contains(&cap), "{role:?} outside the shell");
        assert!(
            store.face(cap).unwrap().inner_wires.is_empty(),
            "{role:?} courtyard must merge into the notched boundary"
        );
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

/// With TWO courtyards, a bottom-touching door into courtyard A merges A
/// into the bottom cap's boundary while courtyard B's wire RIDES ALONG
/// untouched as the kept fragment's inner wire.
#[test]
fn door_preserves_untouched_courtyard_wire() {
    const DUPLEX_DOOR: Opening = Opening {
        op: "cut-door",
        x0: 1.5,
        x1: 2.5,
        z0: -0.5,
        z1: 2.25,
    };
    let mut store = TopologyStore::new();
    let wall = duplex_wall(&mut store);
    let result = subtract(&mut store, wall, &DUPLEX_DOOR);

    assert_eq!(
        welded_boundary_edges(&store, result),
        0,
        "duplex door must weld watertight"
    );
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &DUPLEX_DOOR);

    // Bottom cap: courtyard A merged into the boundary, courtyard B kept.
    let bottom = resolve(&store, &created("duplex1", FaceRole::CapStart));
    let bottom_face = store.face(bottom).unwrap();
    assert_eq!(
        bottom_face.inner_wires.len(),
        1,
        "exactly courtyard B must survive as an inner wire"
    );
    let survivor = store.wire(bottom_face.inner_wires[0]).unwrap();
    for oe in &survivor.edges {
        let edge = store.edge(oe.edge).unwrap();
        for v in [edge.start, edge.end] {
            let p = store.vertex(v).unwrap().point;
            assert!(
                p.x >= 5.0 - 1e-9,
                "surviving inner wire must be courtyard B (vertex at x = {})",
                p.x
            );
        }
    }
    // Top cap: untouched copy with both courtyards.
    let top = resolve(&store, &created("duplex1", FaceRole::CapEnd));
    assert_eq!(
        store.face(top).unwrap().inner_wires.len(),
        2,
        "CapEnd must keep both courtyard wires"
    );
}

/// Second bottom-touching door through the NORTH wall (`y in [5.5, 6.0]`):
/// profile in the XZ plane at `y = 7.0`, extruded 2.5 along `-Y` so the far
/// cap ends inside the courtyard air. Reversed winding keeps the prism
/// outward-facing under the flipped extrusion direction.
fn build_cutter_north(store: &mut TopologyStore, cut: &Opening) -> SolidId {
    let q = |x: f64, z: f64| Point3::new(x, 7.0, z);
    let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
    let profile = vec![
        line(q(cut.x0, cut.z0), q(cut.x0, cut.z1)), // jamb-left
        line(q(cut.x0, cut.z1), q(cut.x1, cut.z1)), // head
        line(q(cut.x1, cut.z1), q(cut.x1, cut.z0)), // jamb-right
        line(q(cut.x1, cut.z0), q(cut.x0, cut.z0)), // sill
    ];
    MakeSegmentedPrism::new(profile, Vector3::new(0.0, -2.5, 0.0))
        .with_op_id(OpId::new(tool_op(cut)))
        .with_segment_tags(tags(["jamb-left", "head", "jamb-right", "sill"]))
        .execute(store)
        .unwrap()
}

/// Bottom-touching door through the WEST wall (`x in [0, 0.5]`): profile in
/// the YZ plane at `x = -1.0`, extruded 2.5 along `+X`. `y0..y1` spans the
/// wall run. Mirrors the south cutter under an X↔Y swap (winding reversed
/// to keep the prism outward-facing after the reflection).
fn build_cutter_west(
    store: &mut TopologyStore,
    op: &'static str,
    y0: f64,
    y1: f64,
    z0: f64,
    z1: f64,
) -> SolidId {
    let q = |y: f64, z: f64| Point3::new(-1.0, y, z);
    let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };
    let profile = vec![
        line(q(y0, z0), q(y0, z1)), // jamb-left
        line(q(y0, z1), q(y1, z1)), // head
        line(q(y1, z1), q(y1, z0)), // jamb-right
        line(q(y1, z0), q(y0, z0)), // sill
    ];
    MakeSegmentedPrism::new(profile, Vector3::new(2.5, 0.0, 0.0))
        .with_op_id(OpId::new(format!("tool-{op}")))
        .with_segment_tags(tags(["jamb-left", "head", "jamb-right", "sill"]))
        .execute(store)
        .unwrap()
}

/// REPRO (bug): a second bottom-touching door on the perpendicular (west)
/// wall. Door 1 (south) opens the annular bottom cap into a single C
/// fragment; door 2 (west) then severs that C into two fragments whose
/// global centroids no longer straddle the canonical closure chord — the
/// ambiguous-SplitSide path, now classified by a local interior probe.
#[test]
fn two_doors_on_perpendicular_walls_cut_both() {
    const DOOR_W: Opening = Opening {
        op: "cut-doorW",
        // Reused as the west wall's y-run span.
        x0: 2.5,
        x1: 3.5,
        z0: -0.5,
        z1: 2.25,
    };
    let mut store = TopologyStore::new();
    let wall = room_wall(&mut store);
    let after_south = subtract(&mut store, wall, &DOOR);
    let cutter = build_cutter_west(
        &mut store, DOOR_W.op, DOOR_W.x0, DOOR_W.x1, DOOR_W.z0, DOOR_W.z1,
    );
    let result = Subtract::new(after_south, cutter)
        .with_op_id(OpId::new(DOOR_W.op))
        .execute(&mut store)
        .unwrap_or_else(|e| panic!("second (west) door failed: {e:?}"));

    assert_eq!(
        welded_boundary_edges(&store, result),
        0,
        "perpendicular two-door room must weld watertight"
    );
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &DOOR);
}

/// REPRO (bug): THREE bottom-touching doors around the room (south, west,
/// north). Every cut after the first severs an already-notched cap, so all
/// three must succeed and the room stay watertight.
#[test]
fn three_doors_around_room_all_cut() {
    const DOOR_W: Opening = Opening {
        op: "cut-doorW",
        x0: 2.5,
        x1: 3.5,
        z0: -0.5,
        z1: 2.25,
    };
    const DOOR_N: Opening = Opening {
        op: "cut-doorN",
        x0: 3.5,
        x1: 4.4,
        z0: -0.5,
        z1: 2.25,
    };
    let mut store = TopologyStore::new();
    let mut current = room_wall(&mut store);
    current = subtract(&mut store, current, &DOOR);

    let west = build_cutter_west(
        &mut store, DOOR_W.op, DOOR_W.x0, DOOR_W.x1, DOOR_W.z0, DOOR_W.z1,
    );
    current = Subtract::new(current, west)
        .with_op_id(OpId::new(DOOR_W.op))
        .execute(&mut store)
        .unwrap_or_else(|e| panic!("west door failed: {e:?}"));

    let north = build_cutter_north(&mut store, &DOOR_N);
    current = Subtract::new(current, north)
        .with_op_id(OpId::new(DOOR_N.op))
        .execute(&mut store)
        .unwrap_or_else(|e| panic!("north door failed: {e:?}"));

    assert_eq!(
        welded_boundary_edges(&store, current),
        0,
        "three-door room must weld watertight"
    );
    let mesh = TessellateSolid::new(current, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &DOOR);
    assert_open(&mesh, &DOOR_N);
}

/// REPRO (bug): a second bottom-touching door on the OPPOSITE (north) wall.
/// Door 1 (south) notches the bottom-cap annulus into a single C fragment;
/// door 2 (north) then severs that C into two fragments, driving the
/// `order_cap_fragments` classification — the ambiguous-SplitSide path.
#[test]
fn two_doors_on_opposite_walls_cut_both() {
    const DOOR_N: Opening = Opening {
        op: "cut-doorN",
        x0: 3.5,
        x1: 4.4,
        z0: -0.5,
        z1: 2.25,
    };
    let mut store = TopologyStore::new();
    let wall = room_wall(&mut store);
    let after_south = subtract(&mut store, wall, &DOOR);
    // Second door on the north wall — this is where the bug fires.
    let cutter = build_cutter_north(&mut store, &DOOR_N);
    let result = Subtract::new(after_south, cutter)
        .with_op_id(OpId::new(DOOR_N.op))
        .execute(&mut store)
        .unwrap_or_else(|e| panic!("second (north) door failed: {e:?}"));

    assert_eq!(
        welded_boundary_edges(&store, result),
        0,
        "two-door room must weld watertight"
    );
    let mesh = TessellateSolid::new(result, TessellationParams::default())
        .execute(&store)
        .unwrap();
    assert_open(&mesh, &DOOR);
    assert_open(&mesh, &DOOR_N);
}

/// R2-style cascade on the annulus in BOTH cut orders: door → window
/// punches the already-notched ring faces; window → door notches the
/// already-punched faces (the hole transfers onto the kept fragment).
#[test]
fn door_and_window_cascade_on_annulus_in_both_orders() {
    for order in [[&DOOR, &WINDOW], [&WINDOW, &DOOR]] {
        let mut store = TopologyStore::new();
        let mut current = room_wall(&mut store);
        for cut in order {
            current = subtract(&mut store, current, cut);
        }
        let label = format!("{} then {}", order[0].op, order[1].op);

        assert_eq!(
            welded_boundary_edges(&store, current),
            0,
            "{label}: door + window room must weld watertight"
        );
        let mesh = TessellateSolid::new(current, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert_open(&mesh, &DOOR);
        assert_open(&mesh, &WINDOW);

        // The notched ring faces keep their names and carry the window
        // hole.
        for tag in ["outer-s", "inner-s"] {
            let face = resolve(&store, &tag_name("room1", tag));
            assert_eq!(
                store.face(face).unwrap().inner_wires.len(),
                1,
                "{label}: {tag} carries the window hole"
            );
        }
        for tag in ["jamb-left", "head", "jamb-right"] {
            assert!(store.names().face(&band_name(&DOOR, tag)).is_some());
        }
        for tag in ["sill", "jamb-right", "head", "jamb-left"] {
            assert!(store.names().face(&band_name(&WINDOW, tag)).is_some());
        }
    }
}
