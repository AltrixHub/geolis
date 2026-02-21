use geolis::geometry::pline::{Pline, PlineVertex};
use geolis::math::{Point3, Vector3};
use geolis::operations::boolean::Subtract;
use geolis::operations::creation::{MakeBox, MakeFace, MakeWire};
use geolis::operations::offset::WallOutline2D;
use geolis::operations::shaping::Extrude;
use geolis::tessellation::{StrokeStyle, TessellateSolid, TessellationParams};
use geolis::topology::TopologyStore;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label, register_stroke};

fn stroke_style() -> StrokeStyle {
    StrokeStyle::new(0.05).unwrap_or_else(|_| unreachable!())
}

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GREEN: Color = Color::rgb(100, 200, 100);
const GRAY: Color = Color::rgb(160, 160, 160);

const EDGE_COLOR: Color = Color::rgb(60, 60, 60);
const WALL_HEIGHT: f64 = 3.0;
const WALL_HALF_WIDTH: f64 = 0.15;

fn render_solid(
    storage: &MeshStorage,
    store: &TopologyStore,
    solid: geolis::topology::SolidId,
    mesh_color: Color,
    edge_color: Color,
) {
    if let Ok(mesh) =
        TessellateSolid::new(solid, TessellationParams::default()).execute(store)
    {
        register_face(storage, mesh, mesh_color);
    }
    if let Ok(solid_data) = store.solid(solid) {
        register_edges(storage, store, solid_data.outer_shell, edge_color);
    }
}

/// Draws a centerline and returns the wall outline points.
fn draw_centerline_and_offset(
    storage: &MeshStorage,
    centerline: &Pline,
    bx: f64,
    by: f64,
) -> Option<Vec<Pline>> {
    let center_pts: Vec<Point3> = centerline
        .vertices
        .iter()
        .map(|v| Point3::new(v.x + bx, v.y + by, 0.0))
        .collect();
    register_stroke(storage, &center_pts, stroke_style(), centerline.closed, GRAY);

    let wall = WallOutline2D::new(vec![centerline.clone()], WALL_HALF_WIDTH);
    wall.execute().ok()
}

pub fn register(storage: &MeshStorage) {
    // Case 1: Simple straight wall with one window
    case_straight_wall(storage, -12.0, 0.0);

    // Case 2: L-shaped wall with window on the longer segment
    case_l_wall(storage, 0.0, 0.0);

    // Case 3: Closed rectangular room with a window
    case_room_with_window(storage, -12.0, -14.0);
}

/// Case 1: Straight wall segment with a window opening.
fn case_straight_wall(storage: &MeshStorage, bx: f64, by: f64) {
    register_label(storage, bx - 1.5, by + 8.0, "1", LABEL_SIZE, LABEL_COLOR);

    let mut store = TopologyStore::new();

    let centerline = Pline {
        vertices: vec![
            PlineVertex::line(0.0, 0.0),
            PlineVertex::line(8.0, 0.0),
        ],
        closed: false,
    };

    let Some(outlines) = draw_centerline_and_offset(storage, &centerline, bx, by) else {
        return;
    };
    let Some(outline) = outlines.into_iter().next() else { return };

    let wall_pts: Vec<Point3> = outline
        .vertices
        .iter()
        .map(|v| Point3::new(v.x + bx, v.y + by, 0.0))
        .collect();

    let Ok(wire) = MakeWire::new(wall_pts, true).execute(&mut store) else { return };
    let Ok(face) = MakeFace::new(wire, vec![]).execute(&mut store) else { return };
    let Ok(wall_solid) = Extrude::new(face, Vector3::new(0.0, 0.0, WALL_HEIGHT))
        .execute(&mut store)
    else {
        return;
    };

    let Ok(window) = MakeBox::new(
        Point3::new(bx + 2.5, by - 1.0, 0.8),
        Point3::new(bx + 5.5, by + 1.0, 2.4),
    )
    .execute(&mut store)
    else {
        return;
    };

    let Ok(result) = Subtract::new(wall_solid, window).execute(&mut store) else {
        return;
    };

    render_solid(storage, &store, result, GREEN, EDGE_COLOR);
}

/// Case 2: L-shaped wall with a window on the longer segment.
fn case_l_wall(storage: &MeshStorage, bx: f64, by: f64) {
    register_label(storage, bx - 1.5, by + 10.0, "2", LABEL_SIZE, LABEL_COLOR);

    let mut store = TopologyStore::new();

    let centerline = Pline {
        vertices: vec![
            PlineVertex::line(0.0, 0.0),
            PlineVertex::line(8.0, 0.0),
            PlineVertex::line(8.0, 6.0),
        ],
        closed: false,
    };

    let Some(outlines) = draw_centerline_and_offset(storage, &centerline, bx, by) else {
        return;
    };
    let Some(outline) = outlines.into_iter().next() else { return };

    let wall_pts: Vec<Point3> = outline
        .vertices
        .iter()
        .map(|v| Point3::new(v.x + bx, v.y + by, 0.0))
        .collect();

    let Ok(wire) = MakeWire::new(wall_pts, true).execute(&mut store) else { return };
    let Ok(face) = MakeFace::new(wire, vec![]).execute(&mut store) else { return };
    let Ok(wall_solid) = Extrude::new(face, Vector3::new(0.0, 0.0, WALL_HEIGHT))
        .execute(&mut store)
    else {
        return;
    };

    let Ok(window) = MakeBox::new(
        Point3::new(bx + 2.0, by - 1.0, 0.8),
        Point3::new(bx + 5.0, by + 1.0, 2.4),
    )
    .execute(&mut store)
    else {
        return;
    };

    let Ok(result) = Subtract::new(wall_solid, window).execute(&mut store) else {
        return;
    };

    render_solid(storage, &store, result, GREEN, EDGE_COLOR);
}

/// Case 3: Closed rectangular room with a window on one wall.
fn case_room_with_window(storage: &MeshStorage, bx: f64, by: f64) {
    register_label(storage, bx - 1.5, by + 12.0, "3", LABEL_SIZE, LABEL_COLOR);

    let mut store = TopologyStore::new();

    let centerline = Pline {
        vertices: vec![
            PlineVertex::line(0.0, 0.0),
            PlineVertex::line(10.0, 0.0),
            PlineVertex::line(10.0, 8.0),
            PlineVertex::line(0.0, 8.0),
        ],
        closed: true,
    };

    let Some(outlines) = draw_centerline_and_offset(storage, &centerline, bx, by) else {
        return;
    };
    let Some(outer) = outlines.first() else { return };

    let wall_pts: Vec<Point3> = outer
        .vertices
        .iter()
        .map(|v| Point3::new(v.x + bx, v.y + by, 0.0))
        .collect();

    let Ok(outer_wire) = MakeWire::new(wall_pts, true).execute(&mut store) else {
        return;
    };

    // If there's an inner outline, use it as a hole in the face
    let inner_wires = if outlines.len() > 1 {
        let inner_pts: Vec<Point3> = outlines[1]
            .vertices
            .iter()
            .map(|v| Point3::new(v.x + bx, v.y + by, 0.0))
            .collect();
        match MakeWire::new(inner_pts, true).execute(&mut store) {
            Ok(w) => vec![w],
            Err(_) => vec![],
        }
    } else {
        vec![]
    };

    let Ok(face) = MakeFace::new(outer_wire, inner_wires).execute(&mut store) else {
        return;
    };
    let Ok(wall_solid) = Extrude::new(face, Vector3::new(0.0, 0.0, WALL_HEIGHT))
        .execute(&mut store)
    else {
        return;
    };

    let Ok(window) = MakeBox::new(
        Point3::new(bx + 3.0, by - 1.0, 0.8),
        Point3::new(bx + 7.0, by + 1.0, 2.4),
    )
    .execute(&mut store)
    else {
        return;
    };

    let Ok(result) = Subtract::new(wall_solid, window).execute(&mut store) else {
        return;
    };

    render_solid(storage, &store, result, GREEN, EDGE_COLOR);
}
