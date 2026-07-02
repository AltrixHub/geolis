//! Curved-wall window showcase (P8) — the capstone of the NURBS boolean program.
//!
//! Builds a plan-arc curved wall with a rounded-rectangle window cut through it
//! (`Subtract(wall, outer_prism)`), a curved frame ring following the wall
//! curvature (`Subtract(Intersect(blank, outer_prism), inner_prism)` — a chained
//! keep-inside then keep-outside boolean), and a semi-transparent curved glass
//! pane tucked behind the frame. Every result is tessellated through the standard
//! pipeline, exercising the keep-inside intersect and the chained boolean
//! end-to-end.

use std::f64::consts::{FRAC_PI_2, PI};

use geolis::error::Result;
use geolis::geometry::nurbs::NurbsCurve3D;
use geolis::math::{Point3, Vector3};
use geolis::operations::boolean::{Intersect, Subtract};
use geolis::operations::creation::{MakeCurvedWall, MakeNurbsPrism};
use geolis::tessellation::{TessellateSolid, TessellationParams};
use geolis::topology::{SolidId, TopologyStore};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label, SceneBounds};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const WALL_COLOR: Color = Color::rgb(150, 180, 205);
const FRAME_COLOR: Color = Color::rgb(90, 105, 125);
const GLASS_COLOR: Color = Color::rgba(140, 200, 230, 90);
const EDGE_COLOR: Color = Color::rgb(255, 255, 255);

// Wall geometry (arc about the origin in plan).
const RADIUS: f64 = 8.0;
const WALL_THICKNESS: f64 = 0.4;
const HEIGHT: f64 = 6.0;
const WINDOW_Z: f64 = 3.0;

/// Unwraps a construction step, logging a warning naming the failed step so
/// pattern failures are visible in the run log instead of silently rendering
/// nothing.
fn ok_or_warn<T>(result: Result<T>, step: &str) -> Option<T> {
    match result {
        Ok(value) => Some(value),
        Err(err) => {
            tracing::warn!(target: "bim", "nurbs_window: {step} failed: {err}");
            None
        }
    }
}

/// Builds the radial window prism of the given rounded-rectangle size at the
/// mid-arc (90 deg). It starts at `RADIUS - 0.6` and is 1.2 deep radially, so it
/// pierces both the wall (thickness 0.4) and the thicker frame blank (0.52).
fn window_prism(
    store: &mut TopologyStore,
    width: f64,
    height: f64,
    corner: f64,
) -> Result<SolidId> {
    // Window center in the profile plane, seated 0.6 inside the mid radius.
    let center = Point3::new(0.0, RADIUS - 0.6, WINDOW_Z);
    let tangential = Vector3::x();
    let vertical = Vector3::z();
    let radial = Vector3::new(0.0, 1.2, 0.0);
    let profile =
        NurbsCurve3D::rounded_rectangle(center, tangential, vertical, width, height, corner)?;
    MakeNurbsPrism::new(profile, radial).execute(store)
}

/// Builds the curved glass pane: a thin curved solid at the mid radius spanning
/// the inner opening plus a small margin, so its edges tuck behind the frame.
///
/// The pane is a real (thin) solid rather than a single open surface: revion's
/// transparent pipeline culls back faces, so an open sheet would vanish when
/// viewed from behind — a closed thin solid presents a front face to either
/// side and blends exactly once per view direction.
fn glass_pane(store: &mut TopologyStore) -> Result<SolidId> {
    // Tangential half-extent 1.2 at RADIUS → angular half-sweep.
    let half_sweep = 1.2 / RADIUS;
    let base_z = 2.1;
    MakeCurvedWall::new(
        Point3::new(0.0, 0.0, base_z),
        RADIUS,
        FRAC_PI_2 - half_sweep,
        FRAC_PI_2 + half_sweep,
        1.8,
        0.02,
    )
    .execute(store)
}

pub fn register(storage: &MeshStorage, bounds: &mut SceneBounds) {
    register_label(storage, -5.0, 9.5, "1", LABEL_SIZE, LABEL_COLOR);

    let deg = |d: f64| d * PI / 180.0;
    let mut store = TopologyStore::new();

    // Curved wall, 55°..125° sweep.
    let Some(wall) = ok_or_warn(
        MakeCurvedWall::new(
            Point3::origin(),
            RADIUS,
            deg(55.0),
            deg(125.0),
            HEIGHT,
            WALL_THICKNESS,
        )
        .execute(&mut store),
        "curved wall",
    ) else {
        return;
    };

    // Outer window box (shared by the wall cut and the frame plate).
    let Some(outer_prism) = ok_or_warn(
        window_prism(&mut store, 2.6, 2.0, 0.35),
        "outer window prism",
    ) else {
        return;
    };
    let Some(wall_open) = ok_or_warn(
        Subtract::new(wall, outer_prism).execute(&mut store),
        "wall window subtract",
    ) else {
        return;
    };

    // Frame ring: keep-inside a thicker blank, then subtract the inner opening.
    let Some(blank) = ok_or_warn(
        MakeCurvedWall::new(
            Point3::origin(),
            RADIUS,
            deg(78.0),
            deg(102.0),
            HEIGHT,
            WALL_THICKNESS + 0.12,
        )
        .execute(&mut store),
        "frame blank",
    ) else {
        return;
    };
    let Some(plate) = ok_or_warn(
        Intersect::new(blank, outer_prism).execute(&mut store),
        "frame plate intersect",
    ) else {
        return;
    };
    let Some(inner_prism) = ok_or_warn(
        window_prism(&mut store, 2.1, 1.5, 0.25),
        "inner window prism",
    ) else {
        return;
    };
    let Some(frame) = ok_or_warn(
        Subtract::new(plate, inner_prism).execute(&mut store),
        "frame opening subtract",
    ) else {
        return;
    };

    // Register the wall and frame solids + their edges.
    if let Some(mesh) = ok_or_warn(
        TessellateSolid::new(wall_open, TessellationParams::default()).execute(&store),
        "wall tessellation",
    ) {
        register_face(storage, bounds, mesh, WALL_COLOR);
    }
    if let Some(mesh) = ok_or_warn(
        TessellateSolid::new(frame, TessellationParams::default()).execute(&store),
        "frame tessellation",
    ) {
        register_face(storage, bounds, mesh, FRAME_COLOR);
    }
    if let Ok(solid) = store.solid(wall_open) {
        register_edges(storage, bounds, &store, solid.outer_shell, EDGE_COLOR);
    }
    if let Ok(solid) = store.solid(frame) {
        register_edges(storage, bounds, &store, solid.outer_shell, EDGE_COLOR);
    }

    // Semi-transparent curved glass pane (alpha < 255 → transparent pipeline).
    if let Some(glass) = ok_or_warn(glass_pane(&mut store), "glass pane") {
        if let Some(mesh) = ok_or_warn(
            TessellateSolid::new(glass, TessellationParams::default()).execute(&store),
            "glass tessellation",
        ) {
            register_face(storage, bounds, mesh, GLASS_COLOR);
        }
    }
}
