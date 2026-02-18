use geolis::math::Point3;
use geolis::operations::creation::{MakeFace, MakeWire};
use geolis::tessellation::{StrokeStyle, TessellateFace, TessellationParams};
use geolis::topology::TopologyStore;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_face, register_label, register_stroke};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GRAY: Color = Color::rgb(180, 180, 180);
const GREEN: Color = Color::rgb(100, 200, 100);
const BLUE: Color = Color::rgb(100, 150, 255);

/// Runs `MakeWire` -> `MakeFace` -> `TessellateFace` pipeline and renders the result.
fn render_face(storage: &MeshStorage, points: &[Point3], outline_color: Color, mesh_color: Color) {
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
    if let Ok(mesh) = TessellateFace::new(face, TessellationParams::default()).execute(&topo) {
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
    render_face(storage, &tri, GRAY, GREEN);

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
    render_face(storage, &sq, GRAY, BLUE);

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
    render_face(storage, &l_shape, GRAY, GREEN);

    // Case 4: Square with hole
    let bx = 14.0;
    let by = 2.0;
    register_label(storage, bx - 2.0, by + 5.5, "4", LABEL_SIZE, LABEL_COLOR);
    let outer = [
        Point3::new(bx, by, 0.0),
        Point3::new(bx + 6.0, by, 0.0),
        Point3::new(bx + 6.0, by + 6.0, 0.0),
        Point3::new(bx, by + 6.0, 0.0),
    ];
    let inner = [
        Point3::new(bx + 1.5, by + 1.5, 0.0),
        Point3::new(bx + 4.5, by + 1.5, 0.0),
        Point3::new(bx + 4.5, by + 4.5, 0.0),
        Point3::new(bx + 1.5, by + 4.5, 0.0),
    ];
    if let Ok(style) = StrokeStyle::new(0.05) {
        register_stroke(storage, &outer, style, true, GRAY);
        register_stroke(storage, &inner, style, true, GRAY);
    }

    let mut topo = TopologyStore::new();
    let Ok(outer_wire) = MakeWire::new(outer.to_vec(), true).execute(&mut topo) else {
        return;
    };
    let Ok(inner_wire) = MakeWire::new(inner.to_vec(), true).execute(&mut topo) else {
        return;
    };
    let Ok(face) = MakeFace::new(outer_wire, vec![inner_wire]).execute(&mut topo) else {
        return;
    };
    if let Ok(mesh) = TessellateFace::new(face, TessellationParams::default()).execute(&topo) {
        register_face(storage, mesh, BLUE);
    }
}
