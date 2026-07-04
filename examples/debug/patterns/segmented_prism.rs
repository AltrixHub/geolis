//! Segmented-prism wall showcase (F5 Phase A).
//!
//! Builds an L-shaped wall footprint as a `MakeSegmentedPrism`: six line
//! segments plus one exact rational arc rounding the outer corner, extruded
//! 3 m along `+Z`. Every plan segment becomes its own named side face
//! (`FaceRole::Tagged`), adjacent side faces share vertical kink edges, and
//! the caps share the per-segment ring edges — so the tessellated solid is
//! watertight by construction. The face mesh and the shared `BRep` edges are
//! both registered so kink edges are visible in the viewport.

use std::f64::consts::FRAC_PI_2;

use geolis::math::{Point3, Vector3};
use geolis::operations::creation::{MakeSegmentedPrism, ProfileSegment};
use geolis::tessellation::{TessellateSolid, TessellationParams};
use geolis::topology::{OpId, SegmentTag, SolidId, TopologyStore};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label, SceneBounds};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const WALL_COLOR: Color = Color::rgb(200, 170, 130);
const EDGE_COLOR: Color = Color::rgb(255, 255, 255);

/// Wall thickness in plan.
const THICKNESS: f64 = 0.4;
/// Wall height along `+Z`.
const HEIGHT: f64 = 3.0;
/// Radius of the rounded outer corner arc.
const CORNER_RADIUS: f64 = 0.4;

pub fn register(storage: &MeshStorage, bounds: &mut SceneBounds) {
    register_label(storage, -1.5, 8.0, "1", LABEL_SIZE, LABEL_COLOR);

    let mut store = TopologyStore::new();
    let Ok(solid) = build_l_wall(&mut store) else {
        return;
    };

    let Ok(mesh) = TessellateSolid::new(solid, TessellationParams::default()).execute(&store)
    else {
        return;
    };
    register_face(storage, bounds, mesh, WALL_COLOR);

    if let Ok(solid_data) = store.solid(solid) {
        register_edges(storage, bounds, &store, solid_data.outer_shell, EDGE_COLOR);
    }
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
