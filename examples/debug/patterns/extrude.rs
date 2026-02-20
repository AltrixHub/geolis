use geolis::math::{Point3, Vector3};
use geolis::operations::creation::{MakeFace, MakeWire};
use geolis::operations::shaping::Extrude;
use geolis::tessellation::{StrokeStyle, TessellateSolid, TessellationParams};
use geolis::topology::TopologyStore;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_face, register_label, register_stroke};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GRAY: Color = Color::rgb(180, 180, 180);
const GREEN: Color = Color::rgb(100, 200, 100);
const BLUE: Color = Color::rgb(100, 150, 255);

const HEIGHT: f64 = 3.0;

/// Runs `MakeWire` -> `MakeFace` -> `Extrude` -> `TessellateSolid` and renders the result.
fn render_extrude(
    storage: &MeshStorage,
    points: &[Point3],
    height: f64,
    outline_color: Color,
    mesh_color: Color,
) {
    if let Ok(style) = StrokeStyle::new(0.05) {
        register_stroke(storage, points, style, true, outline_color);
    }

    let mut topo = TopologyStore::new();
    let Ok(wire) = MakeWire::new(points.to_vec(), true).execute(&mut topo) else {
        return;
    };
    let Ok(face) = MakeFace::new(wire, vec![]).execute(&mut topo) else {
        return;
    };
    let Ok(solid) = Extrude::new(face, Vector3::new(0.0, 0.0, height)).execute(&mut topo) else {
        return;
    };
    if let Ok(mesh) = TessellateSolid::new(solid, TessellationParams::default()).execute(&topo) {
        register_face(storage, mesh, mesh_color);
    }
}

pub fn register(storage: &MeshStorage) {
    // Case 1: Triangle
    let bx = -12.0;
    let by = 2.0;
    register_label(storage, bx - 2.0, by + 5.5, "1", LABEL_SIZE, LABEL_COLOR);
    let tri = [
        Point3::new(bx, by, 0.0),
        Point3::new(bx + 4.0, by, 0.0),
        Point3::new(bx + 2.0, by + 4.0, 0.0),
    ];
    render_extrude(storage, &tri, HEIGHT, GRAY, GREEN);

    // Case 2: Square
    let bx = -4.0;
    let by = 2.0;
    register_label(storage, bx - 2.0, by + 5.5, "2", LABEL_SIZE, LABEL_COLOR);
    let sq = [
        Point3::new(bx, by, 0.0),
        Point3::new(bx + 4.0, by, 0.0),
        Point3::new(bx + 4.0, by + 4.0, 0.0),
        Point3::new(bx, by + 4.0, 0.0),
    ];
    render_extrude(storage, &sq, HEIGHT, GRAY, BLUE);

    // Case 3: L-shape (concave)
    let bx = 4.0;
    let by = 2.0;
    register_label(storage, bx - 2.0, by + 5.5, "3", LABEL_SIZE, LABEL_COLOR);
    let l_shape = [
        Point3::new(bx, by, 0.0),
        Point3::new(bx + 4.0, by, 0.0),
        Point3::new(bx + 4.0, by + 2.0, 0.0),
        Point3::new(bx + 2.0, by + 2.0, 0.0),
        Point3::new(bx + 2.0, by + 4.0, 0.0),
        Point3::new(bx, by + 4.0, 0.0),
    ];
    render_extrude(storage, &l_shape, HEIGHT, GRAY, GREEN);
}
