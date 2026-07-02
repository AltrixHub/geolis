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

use geolis::geometry::nurbs::{NurbsCurve3D, NurbsSurface};
use geolis::math::{Point3, Vector3};
use geolis::operations::boolean::{Intersect, Subtract};
use geolis::operations::creation::{MakeCurvedWall, MakeNurbsPrism};
use geolis::tessellation::{
    tessellate_nurbs_surface, SurfaceTessellationOptions, TessellateSolid, TessellationParams,
};
use geolis::topology::TopologyStore;
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

/// Builds the radial window prism of the given rounded-rectangle size at the
/// mid-arc (90 deg). It starts at `RADIUS - 0.6` and is 1.2 deep radially, so it
/// pierces both the wall (thickness 0.4) and the thicker frame blank (0.52).
fn window_prism(
    store: &mut TopologyStore,
    width: f64,
    height: f64,
    corner: f64,
) -> Option<geolis::topology::SolidId> {
    // Window center in the profile plane, seated 0.6 inside the mid radius.
    let center = Point3::new(0.0, RADIUS - 0.6, WINDOW_Z);
    let tangential = Vector3::x();
    let vertical = Vector3::z();
    let radial = Vector3::new(0.0, 1.2, 0.0);
    let profile =
        NurbsCurve3D::rounded_rectangle(center, tangential, vertical, width, height, corner)
            .ok()?;
    MakeNurbsPrism::new(profile, radial).execute(store).ok()
}

/// Builds the curved glass pane: a mid-radius sub-arc (spanning the inner opening
/// plus a small margin) extruded vertically, so its edges tuck behind the frame.
fn glass_surface() -> Option<NurbsSurface> {
    // Tangential half-extent 1.2 at RADIUS → angular half-sweep.
    let half_sweep = 1.2 / RADIUS;
    let base_z = 2.1;
    let arc = NurbsCurve3D::arc(
        Point3::new(0.0, 0.0, base_z),
        RADIUS,
        Vector3::z(),
        Vector3::x(),
        FRAC_PI_2 - half_sweep,
        FRAC_PI_2 + half_sweep,
    )
    .ok()?;
    NurbsSurface::extrude(&arc, Vector3::new(0.0, 0.0, 1.8)).ok()
}

pub fn register(storage: &MeshStorage, bounds: &mut SceneBounds) {
    register_label(storage, -5.0, 9.5, "1", LABEL_SIZE, LABEL_COLOR);

    let deg = |d: f64| d * PI / 180.0;
    let mut store = TopologyStore::new();

    // Curved wall, 55°..125° sweep.
    let Ok(wall) = MakeCurvedWall::new(
        Point3::origin(),
        RADIUS,
        deg(55.0),
        deg(125.0),
        HEIGHT,
        WALL_THICKNESS,
    )
    .execute(&mut store) else {
        return;
    };

    // Outer window box (shared by the wall cut and the frame plate).
    let Some(outer_prism) = window_prism(&mut store, 2.6, 2.0, 0.35) else {
        return;
    };
    let Ok(wall_open) = Subtract::new(wall, outer_prism).execute(&mut store) else {
        return;
    };

    // Frame ring: keep-inside a thicker blank, then subtract the inner opening.
    let Ok(blank) = MakeCurvedWall::new(
        Point3::origin(),
        RADIUS,
        deg(78.0),
        deg(102.0),
        HEIGHT,
        WALL_THICKNESS + 0.12,
    )
    .execute(&mut store) else {
        return;
    };
    let Ok(plate) = Intersect::new(blank, outer_prism).execute(&mut store) else {
        return;
    };
    // Inner opening: corner 0.3 (the SSI marcher fragments a tighter 0.25 loop).
    let Some(inner_prism) = window_prism(&mut store, 2.1, 1.5, 0.3) else {
        return;
    };
    let Ok(frame) = Subtract::new(plate, inner_prism).execute(&mut store) else {
        return;
    };

    // Register the wall and frame solids + their edges.
    if let Ok(mesh) = TessellateSolid::new(wall_open, TessellationParams::default()).execute(&store)
    {
        register_face(storage, bounds, mesh, WALL_COLOR);
    }
    if let Ok(mesh) = TessellateSolid::new(frame, TessellationParams::default()).execute(&store) {
        register_face(storage, bounds, mesh, FRAME_COLOR);
    }
    if let Ok(solid) = store.solid(wall_open) {
        register_edges(storage, bounds, &store, solid.outer_shell, EDGE_COLOR);
    }
    if let Ok(solid) = store.solid(frame) {
        register_edges(storage, bounds, &store, solid.outer_shell, EDGE_COLOR);
    }

    // Semi-transparent curved glass pane (alpha < 255 → transparent pipeline).
    if let Some(glass) = glass_surface() {
        let options = SurfaceTessellationOptions::default();
        if let Ok(mesh) = tessellate_nurbs_surface(&glass, &options) {
            register_face(storage, bounds, mesh, GLASS_COLOR);
        }
    }
}
