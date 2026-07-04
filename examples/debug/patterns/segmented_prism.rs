//! Segmented-prism wall showcase (F5 Phases A + B).
//!
//! Variant 1 builds an L-shaped wall footprint as a `MakeSegmentedPrism`:
//! six line segments plus one exact rational arc rounding the outer corner,
//! extruded 3 m along `+Z`. Every plan segment becomes its own named side
//! face (`FaceRole::Tagged`), adjacent side faces share vertical kink edges,
//! and the caps share the per-segment ring edges — so the tessellated solid
//! is watertight by construction.
//!
//! Variant 2 (Phase B) cuts a genuine 4-side-face box cutter through a
//! straight segmented-prism wall: every (wall face × box face) SSI branch is
//! open and ends on the box's kink edges, the stitcher chains them into the
//! entry/exit window loops, and the subtract emits one named band fragment
//! per box side face. The face mesh and the shared `BRep` edges are both
//! registered so kink and rim edges are visible in the viewport.
//!
//! Variant 3 (Phase C) slides the window ACROSS a wall joint: the outer and
//! inner sides are each segmented into two collinear tagged pieces sharing a
//! vertical kink edge at mid-span, and the box cutter straddles that joint.
//! The SSI branches end on the TARGET faces' shared boundary, the stitcher
//! chains them across the target kink, and the F3b splitter applies each
//! hole half as a boundary notch — the notched faces keep their tagged
//! names, and the split kink sub-edges (above and below the window) are
//! shared by the neighboring fragments.
//!
//! Variant 4 (F6 Phase R1) cascades three box cutters through ONE wall —
//! door + window + window subtracted SEQUENTIALLY, each cut operating on
//! the already-punched result of the previous one. The earlier cuts' band
//! fragments ride through the later cuts as target faces, so the finished
//! wall carries three named openings whose reveals all resolve by name.
//!
//! Variant 5 (F6 Phase R2) cuts a FULL-HEIGHT door: the cutter reaches
//! below the wall's bottom cap, so the entry / exit traces are open
//! boundary notches on the wall faces, the doorway walls close onto
//! cap-plane closure edges, and the bottom cap is rebuilt as two named
//! `Split` fragments sharing those exact edges — a genuinely open,
//! watertight doorway. A window cut follows on the same wall (the R2
//! cascade: a later cut on an already-notched face).
//!
//! Variant 6 (F6 Phase R3) cuts FLUSH doors: the first cutter's sill plane
//! is EXACTLY coplanar with the wall's bottom cap (`z0 == 0.0`) and the
//! second is flush with BOTH caps (`z0 == 0.0`, `z1 == 3.0`). The sill /
//! head faces' tangential along-boundary SSI branches are dropped and the
//! jamb corner terminals are accepted as target-boundary terminals, so the
//! degenerate contacts cleanly reduce to the R2 notch / sever cuts — no
//! sliver geometry, watertight by construction.
//!
//! Variant 7 (F6 Phase R3b) is a plan-ARC curved wall: an annular
//! segmented-prism strip (outer arc + end lines + inner arc, exact
//! rational `Arc` profile segments) minus a FULL-HEIGHT door and a window,
//! both cut through the curved faces by radial 4-side-face box cutters.
//! The door notches the annular bottom cap along ARC sub-span edges shared
//! with the curved wall faces; the window punches a chained hole through
//! both cylindrical faces. Both cuts are clean and watertight.

use std::f64::consts::{FRAC_PI_2, PI};

use geolis::math::{Point3, Vector3};
use geolis::operations::boolean::Subtract;
use geolis::operations::creation::{MakeSegmentedPrism, ProfileSegment};
use geolis::tessellation::{TessellateSolid, TessellationParams};
use geolis::topology::{OpId, SegmentTag, SolidId, TopologyStore};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label, SceneBounds};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const WALL_COLOR: Color = Color::rgb(200, 170, 130);
const CUT_WALL_COLOR: Color = Color::rgb(150, 180, 210);
const KINK_WALL_COLOR: Color = Color::rgb(170, 210, 160);
const CASCADE_WALL_COLOR: Color = Color::rgb(210, 160, 190);
const DOOR_WALL_COLOR: Color = Color::rgb(160, 200, 200);
const FLUSH_WALL_COLOR: Color = Color::rgb(220, 200, 150);
const CURVED_WALL_COLOR: Color = Color::rgb(180, 170, 220);
const EDGE_COLOR: Color = Color::rgb(255, 255, 255);

/// Wall thickness in plan.
const THICKNESS: f64 = 0.4;
/// Wall height along `+Z`.
const HEIGHT: f64 = 3.0;
/// Radius of the rounded outer corner arc.
const CORNER_RADIUS: f64 = 0.4;

/// One demo variant builder: builds the solid into a fresh store.
type WallBuilder = fn(&mut TopologyStore) -> geolis::Result<SolidId>;

pub fn register(storage: &MeshStorage, bounds: &mut SceneBounds) {
    let variants: [(f64, &str, Color, WallBuilder); 7] = [
        (8.0, "1", WALL_COLOR, build_l_wall),
        (13.0, "2", CUT_WALL_COLOR, build_window_wall),
        (18.0, "3", KINK_WALL_COLOR, build_kink_window_wall),
        (23.0, "4", CASCADE_WALL_COLOR, build_cascade_wall),
        (28.0, "5", DOOR_WALL_COLOR, build_full_height_door_wall),
        (33.0, "6", FLUSH_WALL_COLOR, build_flush_door_wall),
        (38.0, "7", CURVED_WALL_COLOR, build_curved_opening_wall),
    ];
    for (label_y, label, color, build) in variants {
        register_label(storage, -1.5, label_y, label, LABEL_SIZE, LABEL_COLOR);
        if !register_variant(storage, bounds, color, build) {
            return;
        }
    }
}

/// Builds one variant into a fresh store, tessellates it, and registers
/// its face mesh and shared `BRep` edges. Returns `false` when the build
/// or tessellation fails (the caller stops registering further variants).
fn register_variant(
    storage: &MeshStorage,
    bounds: &mut SceneBounds,
    color: Color,
    build: WallBuilder,
) -> bool {
    let mut store = TopologyStore::new();
    let Ok(solid) = build(&mut store) else {
        return false;
    };
    let Ok(mesh) = TessellateSolid::new(solid, TessellationParams::default()).execute(&store)
    else {
        return false;
    };
    register_face(storage, bounds, mesh, color);
    if let Ok(solid_data) = store.solid(solid) {
        register_edges(storage, bounds, &store, solid_data.outer_shell, EDGE_COLOR);
    }
    true
}

/// Builds the L-shaped segmented-prism wall: an L footprint (legs 6 m and
/// 4 m, thickness 0.4 m) whose outer corner at the leg junction is rounded
/// by an exact rational arc. One tagged side face per plan segment.
fn build_l_wall(store: &mut TopologyStore) -> geolis::Result<SolidId> {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);

    // CCW footprint: outer boundary of an L-shaped wall strip.
    let profile = vec![
        // Outer run along the long leg, stopping short of the corner arc.
        ProfileSegment::Line {
            start: p(0.0, 0.0),
            end: p(6.0 - CORNER_RADIUS, 0.0),
        },
        // Rounded outer corner: exact quarter arc about the corner center.
        ProfileSegment::Arc {
            center: Point3::new(6.0 - CORNER_RADIUS, CORNER_RADIUS, 0.0),
            radius: CORNER_RADIUS,
            normal: Vector3::z(),
            ref_dir: Vector3::x(),
            start_angle: -FRAC_PI_2,
            end_angle: 0.0,
        },
        // Outer run up the short leg.
        ProfileSegment::Line {
            start: p(6.0, CORNER_RADIUS),
            end: p(6.0, 4.0),
        },
        // Short-leg end cap.
        ProfileSegment::Line {
            start: p(6.0, 4.0),
            end: p(6.0 - THICKNESS, 4.0),
        },
        // Inner run down the short leg.
        ProfileSegment::Line {
            start: p(6.0 - THICKNESS, 4.0),
            end: p(6.0 - THICKNESS, THICKNESS),
        },
        // Inner run along the long leg.
        ProfileSegment::Line {
            start: p(6.0 - THICKNESS, THICKNESS),
            end: p(0.0, THICKNESS),
        },
        // Long-leg end cap.
        ProfileSegment::Line {
            start: p(0.0, THICKNESS),
            end: p(0.0, 0.0),
        },
    ];
    let tags = [
        "outer-long",
        "outer-corner",
        "outer-short",
        "end-short",
        "inner-short",
        "inner-long",
        "end-long",
    ]
    .iter()
    .map(|t| SegmentTag::new(*t))
    .collect();

    MakeSegmentedPrism::new(profile, Vector3::new(0.0, 0.0, HEIGHT))
        .with_op_id(OpId::new("demo-l-wall"))
        .with_segment_tags(tags)
        .execute(store)
}

/// Builds the Phase B variant: a straight segmented-prism wall (6 m long,
/// 0.4 m thick, at plan `y` 8..8.4) minus a genuine 4-side-face box cutter
/// extruded horizontally through it — a through window whose reveal is four
/// named band fragments (`Band { op, tool_face: Tagged(sill/head/jambs) }`).
fn build_window_wall(store: &mut TopologyStore) -> geolis::Result<SolidId> {
    const Y0: f64 = 8.0;
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);

    let wall_profile = vec![
        ProfileSegment::Line {
            start: p(0.0, Y0),
            end: p(6.0, Y0),
        },
        ProfileSegment::Line {
            start: p(6.0, Y0),
            end: p(6.0, Y0 + THICKNESS),
        },
        ProfileSegment::Line {
            start: p(6.0, Y0 + THICKNESS),
            end: p(0.0, Y0 + THICKNESS),
        },
        ProfileSegment::Line {
            start: p(0.0, Y0 + THICKNESS),
            end: p(0.0, Y0),
        },
    ];
    let wall_tags = ["outer", "end-east", "inner", "end-west"]
        .iter()
        .map(|t| SegmentTag::new(*t))
        .collect();
    let wall = MakeSegmentedPrism::new(wall_profile, Vector3::new(0.0, 0.0, HEIGHT))
        .with_op_id(OpId::new("demo-window-wall"))
        .with_segment_tags(wall_tags)
        .execute(store)?;

    // Window rectangle in the vertical XZ plane in front of the wall,
    // extruded through it along +Y.
    let q = |x: f64, z: f64| Point3::new(x, Y0 - 1.0, z);
    let cutter_profile = vec![
        ProfileSegment::Line {
            start: q(2.0, 1.0),
            end: q(3.6, 1.0),
        },
        ProfileSegment::Line {
            start: q(3.6, 1.0),
            end: q(3.6, 2.2),
        },
        ProfileSegment::Line {
            start: q(3.6, 2.2),
            end: q(2.0, 2.2),
        },
        ProfileSegment::Line {
            start: q(2.0, 2.2),
            end: q(2.0, 1.0),
        },
    ];
    let cutter_tags = ["sill", "jamb-right", "head", "jamb-left"]
        .iter()
        .map(|t| SegmentTag::new(*t))
        .collect();
    let cutter = MakeSegmentedPrism::new(cutter_profile, Vector3::new(0.0, 2.4, 0.0))
        .with_op_id(OpId::new("demo-window-box"))
        .with_segment_tags(cutter_tags)
        .execute(store)?;

    Subtract::new(wall, cutter)
        .with_op_id(OpId::new("demo-window-cut"))
        .execute(store)
}

/// Builds the Phase C variant: a straight wall whose outer and inner sides
/// are each segmented into two collinear tagged pieces joined at mid-span
/// (x = 3), minus a box window straddling that joint. The window's hole
/// halves become boundary notches on the four adjacent wall faces, which
/// keep their tagged names; the split joint edges above and below the
/// window are shared by the neighboring fragments.
fn build_kink_window_wall(store: &mut TopologyStore) -> geolis::Result<SolidId> {
    const Y0: f64 = 13.0;
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);

    let wall_profile = vec![
        ProfileSegment::Line {
            start: p(0.0, Y0),
            end: p(3.0, Y0),
        },
        ProfileSegment::Line {
            start: p(3.0, Y0),
            end: p(6.0, Y0),
        },
        ProfileSegment::Line {
            start: p(6.0, Y0),
            end: p(6.0, Y0 + THICKNESS),
        },
        ProfileSegment::Line {
            start: p(6.0, Y0 + THICKNESS),
            end: p(3.0, Y0 + THICKNESS),
        },
        ProfileSegment::Line {
            start: p(3.0, Y0 + THICKNESS),
            end: p(0.0, Y0 + THICKNESS),
        },
        ProfileSegment::Line {
            start: p(0.0, Y0 + THICKNESS),
            end: p(0.0, Y0),
        },
    ];
    let wall_tags = [
        "outer-west",
        "outer-east",
        "end-east",
        "inner-east",
        "inner-west",
        "end-west",
    ]
    .iter()
    .map(|t| SegmentTag::new(*t))
    .collect();
    let wall = MakeSegmentedPrism::new(wall_profile, Vector3::new(0.0, 0.0, HEIGHT))
        .with_op_id(OpId::new("demo-kink-wall"))
        .with_segment_tags(wall_tags)
        .execute(store)?;

    // Window straddling the x = 3 joint, extruded through the wall.
    let q = |x: f64, z: f64| Point3::new(x, Y0 - 1.0, z);
    let cutter_profile = vec![
        ProfileSegment::Line {
            start: q(2.2, 1.0),
            end: q(3.8, 1.0),
        },
        ProfileSegment::Line {
            start: q(3.8, 1.0),
            end: q(3.8, 2.2),
        },
        ProfileSegment::Line {
            start: q(3.8, 2.2),
            end: q(2.2, 2.2),
        },
        ProfileSegment::Line {
            start: q(2.2, 2.2),
            end: q(2.2, 1.0),
        },
    ];
    let cutter_tags = ["sill", "jamb-right", "head", "jamb-left"]
        .iter()
        .map(|t| SegmentTag::new(*t))
        .collect();
    let cutter = MakeSegmentedPrism::new(cutter_profile, Vector3::new(0.0, 2.4, 0.0))
        .with_op_id(OpId::new("demo-kink-window-box"))
        .with_segment_tags(cutter_tags)
        .execute(store)?;

    Subtract::new(wall, cutter)
        .with_op_id(OpId::new("demo-kink-window-cut"))
        .execute(store)
}

/// Builds the F6 Phase R1 variant: one straight tagged wall minus a door
/// and two windows, subtracted SEQUENTIALLY — each cut punches the result
/// of the previous one, and every cut carries its own op id so all three
/// openings' band fragments resolve by name on the final solid.
fn build_cascade_wall(store: &mut TopologyStore) -> geolis::Result<SolidId> {
    const Y0: f64 = 18.0;
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);

    let wall_profile = vec![
        ProfileSegment::Line {
            start: p(0.0, Y0),
            end: p(6.0, Y0),
        },
        ProfileSegment::Line {
            start: p(6.0, Y0),
            end: p(6.0, Y0 + THICKNESS),
        },
        ProfileSegment::Line {
            start: p(6.0, Y0 + THICKNESS),
            end: p(0.0, Y0 + THICKNESS),
        },
        ProfileSegment::Line {
            start: p(0.0, Y0 + THICKNESS),
            end: p(0.0, Y0),
        },
    ];
    let wall_tags = ["outer", "end-east", "inner", "end-west"]
        .iter()
        .map(|t| SegmentTag::new(*t))
        .collect();
    let wall = MakeSegmentedPrism::new(wall_profile, Vector3::new(0.0, 0.0, HEIGHT))
        .with_op_id(OpId::new("demo-cascade-wall"))
        .with_segment_tags(wall_tags)
        .execute(store)?;

    // Door + two windows: interior, non-overlapping openings, each a
    // genuine 4-side-face box extruded through the wall along +Y.
    let openings: [(&str, f64, f64, f64, f64); 3] = [
        ("demo-cascade-door", 0.6, 1.5, 0.15, 2.25),
        ("demo-cascade-win-a", 2.2, 3.4, 1.0, 2.0),
        ("demo-cascade-win-b", 4.2, 5.4, 1.0, 2.0),
    ];
    let mut current = wall;
    for (op, x0, x1, z0, z1) in openings {
        let q = |x: f64, z: f64| Point3::new(x, Y0 - 1.0, z);
        let cutter_profile = vec![
            ProfileSegment::Line {
                start: q(x0, z0),
                end: q(x1, z0),
            },
            ProfileSegment::Line {
                start: q(x1, z0),
                end: q(x1, z1),
            },
            ProfileSegment::Line {
                start: q(x1, z1),
                end: q(x0, z1),
            },
            ProfileSegment::Line {
                start: q(x0, z1),
                end: q(x0, z0),
            },
        ];
        let cutter_tags = ["sill", "jamb-right", "head", "jamb-left"]
            .iter()
            .map(|t| SegmentTag::new(*t))
            .collect();
        let cutter = MakeSegmentedPrism::new(cutter_profile, Vector3::new(0.0, 2.4, 0.0))
            .with_op_id(OpId::new(format!("{op}-box")))
            .with_segment_tags(cutter_tags)
            .execute(store)?;
        current = Subtract::new(current, cutter)
            .with_op_id(OpId::new(op))
            .execute(store)?;
    }
    Ok(current)
}

/// Builds the F6 Phase R3 variant: one straight tagged wall minus a
/// FLUSH-SILL door (`z0 == 0.0`, sill exactly coplanar with the bottom
/// cap) and a door flush with BOTH caps (`z0 == 0.0`, `z1 == 3.0`) —
/// exact degenerate contacts that cleanly reduce to the R2 notch / sever
/// cuts.
fn build_flush_door_wall(store: &mut TopologyStore) -> geolis::Result<SolidId> {
    const Y0: f64 = 28.0;
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);

    let wall_profile = vec![
        ProfileSegment::Line {
            start: p(0.0, Y0),
            end: p(6.0, Y0),
        },
        ProfileSegment::Line {
            start: p(6.0, Y0),
            end: p(6.0, Y0 + THICKNESS),
        },
        ProfileSegment::Line {
            start: p(6.0, Y0 + THICKNESS),
            end: p(0.0, Y0 + THICKNESS),
        },
        ProfileSegment::Line {
            start: p(0.0, Y0 + THICKNESS),
            end: p(0.0, Y0),
        },
    ];
    let wall_tags = ["outer", "end-east", "inner", "end-west"]
        .iter()
        .map(|t| SegmentTag::new(*t))
        .collect();
    let wall = MakeSegmentedPrism::new(wall_profile, Vector3::new(0.0, 0.0, HEIGHT))
        .with_op_id(OpId::new("demo-flush-wall"))
        .with_segment_tags(wall_tags)
        .execute(store)?;

    // Flush-sill door (coplanar with the bottom cap) then a door flush
    // with both caps, sequentially.
    let openings: [(&str, f64, f64, f64, f64); 2] = [
        ("demo-flush-door", 0.8, 1.7, 0.0, 2.25),
        ("demo-flush-tall-door", 3.2, 4.1, 0.0, HEIGHT),
    ];
    let mut current = wall;
    for (op, x0, x1, z0, z1) in openings {
        let q = |x: f64, z: f64| Point3::new(x, Y0 - 1.0, z);
        let cutter_profile = vec![
            ProfileSegment::Line {
                start: q(x0, z0),
                end: q(x1, z0),
            },
            ProfileSegment::Line {
                start: q(x1, z0),
                end: q(x1, z1),
            },
            ProfileSegment::Line {
                start: q(x1, z1),
                end: q(x0, z1),
            },
            ProfileSegment::Line {
                start: q(x0, z1),
                end: q(x0, z0),
            },
        ];
        let cutter_tags = ["sill", "jamb-right", "head", "jamb-left"]
            .iter()
            .map(|t| SegmentTag::new(*t))
            .collect();
        let cutter = MakeSegmentedPrism::new(cutter_profile, Vector3::new(0.0, 2.4, 0.0))
            .with_op_id(OpId::new(format!("{op}-box")))
            .with_segment_tags(cutter_tags)
            .execute(store)?;
        current = Subtract::new(current, cutter)
            .with_op_id(OpId::new(op))
            .execute(store)?;
    }
    Ok(current)
}

/// Builds the F6 Phase R3b variant: a plan-ARC curved wall — an annular
/// segmented-prism strip (outer arc r = 8.4, inner arc r = 8.0, azimuth
/// 60..120 degrees about a center at plan `y = 30`) minus a FULL-HEIGHT
/// door through the arc faces (the bottom cap notched along arc sub-span
/// edges) and a window punched through both cylindrical faces.
fn build_curved_opening_wall(store: &mut TopologyStore) -> geolis::Result<SolidId> {
    const CENTER_Y: f64 = 30.0;
    const R_OUTER: f64 = 8.4;
    const R_INNER: f64 = 8.0;
    let deg = |d: f64| d * PI / 180.0;
    let at = |r: f64, a: f64| Point3::new(r * a.cos(), CENTER_Y + r * a.sin(), 0.0);
    let center = Point3::new(0.0, CENTER_Y, 0.0);

    // Annular wall strip; the inner arc is traversed backwards via the -Z
    // normal (the Phase A fillet convention).
    let wall_profile = vec![
        ProfileSegment::Arc {
            center,
            radius: R_OUTER,
            normal: Vector3::z(),
            ref_dir: Vector3::x(),
            start_angle: deg(60.0),
            end_angle: deg(120.0),
        },
        ProfileSegment::Line {
            start: at(R_OUTER, deg(120.0)),
            end: at(R_INNER, deg(120.0)),
        },
        ProfileSegment::Arc {
            center,
            radius: R_INNER,
            normal: -Vector3::z(),
            ref_dir: Vector3::x(),
            start_angle: deg(-120.0),
            end_angle: deg(-60.0),
        },
        ProfileSegment::Line {
            start: at(R_INNER, deg(60.0)),
            end: at(R_OUTER, deg(60.0)),
        },
    ];
    let wall_tags = ["convex", "end-west", "concave", "end-east"]
        .iter()
        .map(|t| SegmentTag::new(*t))
        .collect();
    let wall = MakeSegmentedPrism::new(wall_profile, Vector3::new(0.0, 0.0, HEIGHT))
        .with_op_id(OpId::new("demo-curved-wall"))
        .with_segment_tags(wall_tags)
        .execute(store)?;

    // Full-height door through the arc apex, then a window east of it.
    // Both cutters extrude along +Y from inside the annulus, crossing the
    // concave face first and exiting through the convex face.
    let openings: [(&str, f64, f64, f64, f64); 2] = [
        ("demo-curved-door", -0.45, 0.45, -0.5, 2.25),
        ("demo-curved-win", 1.4, 2.6, 1.0, 2.0),
    ];
    let mut current = wall;
    for (op, x0, x1, z0, z1) in openings {
        let q = |x: f64, z: f64| Point3::new(x, CENTER_Y + 6.0, z);
        let cutter_profile = vec![
            ProfileSegment::Line {
                start: q(x0, z0),
                end: q(x1, z0),
            },
            ProfileSegment::Line {
                start: q(x1, z0),
                end: q(x1, z1),
            },
            ProfileSegment::Line {
                start: q(x1, z1),
                end: q(x0, z1),
            },
            ProfileSegment::Line {
                start: q(x0, z1),
                end: q(x0, z0),
            },
        ];
        let cutter_tags = ["sill", "jamb-right", "head", "jamb-left"]
            .iter()
            .map(|t| SegmentTag::new(*t))
            .collect();
        let cutter = MakeSegmentedPrism::new(cutter_profile, Vector3::new(0.0, 3.5, 0.0))
            .with_op_id(OpId::new(format!("{op}-box")))
            .with_segment_tags(cutter_tags)
            .execute(store)?;
        current = Subtract::new(current, cutter)
            .with_op_id(OpId::new(op))
            .execute(store)?;
    }
    Ok(current)
}

/// Builds the F6 Phase R2 variant: one straight tagged wall minus a
/// FULL-HEIGHT door (the cutter's sill lies below the wall, its head inside
/// it) and a window. The door's cut chains are OPEN — they terminate on the
/// wall's bottom ring — so the bottom cap is notched into two `Split`
/// fragments sharing the doorway's closure edges, and the window then cuts
/// the already-notched wall faces.
fn build_full_height_door_wall(store: &mut TopologyStore) -> geolis::Result<SolidId> {
    const Y0: f64 = 23.0;
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);

    let wall_profile = vec![
        ProfileSegment::Line {
            start: p(0.0, Y0),
            end: p(6.0, Y0),
        },
        ProfileSegment::Line {
            start: p(6.0, Y0),
            end: p(6.0, Y0 + THICKNESS),
        },
        ProfileSegment::Line {
            start: p(6.0, Y0 + THICKNESS),
            end: p(0.0, Y0 + THICKNESS),
        },
        ProfileSegment::Line {
            start: p(0.0, Y0 + THICKNESS),
            end: p(0.0, Y0),
        },
    ];
    let wall_tags = ["outer", "end-east", "inner", "end-west"]
        .iter()
        .map(|t| SegmentTag::new(*t))
        .collect();
    let wall = MakeSegmentedPrism::new(wall_profile, Vector3::new(0.0, 0.0, HEIGHT))
        .with_op_id(OpId::new("demo-door-wall"))
        .with_segment_tags(wall_tags)
        .execute(store)?;

    // Full-height door (sill below the wall) then a window, sequentially.
    let openings: [(&str, f64, f64, f64, f64); 2] = [
        ("demo-full-door", 0.8, 1.7, -0.5, 2.25),
        ("demo-door-win", 3.2, 4.4, 1.0, 2.0),
    ];
    let mut current = wall;
    for (op, x0, x1, z0, z1) in openings {
        let q = |x: f64, z: f64| Point3::new(x, Y0 - 1.0, z);
        let cutter_profile = vec![
            ProfileSegment::Line {
                start: q(x0, z0),
                end: q(x1, z0),
            },
            ProfileSegment::Line {
                start: q(x1, z0),
                end: q(x1, z1),
            },
            ProfileSegment::Line {
                start: q(x1, z1),
                end: q(x0, z1),
            },
            ProfileSegment::Line {
                start: q(x0, z1),
                end: q(x0, z0),
            },
        ];
        let cutter_tags = ["sill", "jamb-right", "head", "jamb-left"]
            .iter()
            .map(|t| SegmentTag::new(*t))
            .collect();
        let cutter = MakeSegmentedPrism::new(cutter_profile, Vector3::new(0.0, 2.4, 0.0))
            .with_op_id(OpId::new(format!("{op}-box")))
            .with_segment_tags(cutter_tags)
            .execute(store)?;
        current = Subtract::new(current, cutter)
            .with_op_id(OpId::new(op))
            .execute(store)?;
    }
    Ok(current)
}
