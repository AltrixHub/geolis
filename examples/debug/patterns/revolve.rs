use geolis::math::{Point3, Vector3};
use geolis::operations::creation::{MakeFace, MakeWire};
use geolis::operations::shaping::Revolve;
use geolis::tessellation::{StrokeStyle, TessellateSolid, TessellationParams};
use geolis::topology::TopologyStore;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label, register_stroke};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GRAY: Color = Color::rgb(180, 180, 180);
const GREEN: Color = Color::rgb(100, 200, 100);
const BLUE: Color = Color::rgb(100, 150, 255);
const RED: Color = Color::rgb(220, 80, 80);

/// Runs `MakeWire` -> `MakeFace` -> `Revolve` -> `TessellateSolid` and renders.
fn render_revolve(
    storage: &MeshStorage,
    points: &[Point3],
    axis_origin: Point3,
    axis_dir: Vector3,
    outline_color: Color,
    mesh_color: Color,
) {
    // Draw profile outline
    if let Ok(style) = StrokeStyle::new(0.05) {
        register_stroke(storage, points, style, true, outline_color);
    }

    // Draw axis as a thin line
    if let Ok(style) = StrokeStyle::new(0.02) {
        let axis_len = 8.0;
        let axis_line = [
            axis_origin - axis_dir * (axis_len / 2.0),
            axis_origin + axis_dir * (axis_len / 2.0),
        ];
        register_stroke(storage, &axis_line, style, false, RED);
    }

    let mut topo = TopologyStore::new();
    let Ok(wire) = MakeWire::new(points.to_vec(), true).execute(&mut topo) else {
        return;
    };
    let Ok(face) = MakeFace::new(wire, vec![]).execute(&mut topo) else {
        return;
    };
    let Ok(solid) = Revolve::new(face, axis_origin, axis_dir).execute(&mut topo) else {
        return;
    };
    if let Ok(mesh) = TessellateSolid::new(solid, TessellationParams::default()).execute(&topo) {
        register_face(storage, mesh, mesh_color);
    }

    if let Ok(solid_data) = topo.solid(solid) {
        register_edges(storage, &topo, solid_data.outer_shell, outline_color);
    }
}

pub fn register(storage: &MeshStorage) {
    // Case 1: Square profile → hollow cylinder (annulus)
    //   Profile: rectangle at x=2..4, z=0..3
    //   Axis: Z-axis
    let bx = -14.0;
    let by = -4.0;
    register_label(storage, bx - 2.0, by + 10.0, "1", LABEL_SIZE, LABEL_COLOR);
    let square_profile = [
        Point3::new(bx + 16.0 + 2.0, 0.0, by),
        Point3::new(bx + 16.0 + 4.0, 0.0, by),
        Point3::new(bx + 16.0 + 4.0, 0.0, by + 6.0),
        Point3::new(bx + 16.0 + 2.0, 0.0, by + 6.0),
    ];
    render_revolve(
        storage,
        &square_profile,
        Point3::new(bx + 16.0, 0.0, 0.0),
        Vector3::z(),
        GRAY,
        GREEN,
    );

    // Case 2: Triangle with vertex on axis → cone
    //   Profile: triangle with apex at (0, 0, 5) on axis
    let offset_x = 12.0;
    register_label(storage, offset_x - 2.0, by + 10.0, "2", LABEL_SIZE, LABEL_COLOR);
    let tri_profile = [
        Point3::new(offset_x, 0.0, by + 6.0),        // on axis
        Point3::new(offset_x + 3.0, 0.0, by),         // off axis
        Point3::new(offset_x + 3.0, 0.0, by + 6.0),   // off axis
    ];
    render_revolve(
        storage,
        &tri_profile,
        Point3::new(offset_x, 0.0, 0.0),
        Vector3::z(),
        GRAY,
        BLUE,
    );

    // Case 3: Trapezoid → truncated cone (frustum)
    //   Profile: trapezoid with narrower top
    let offset_x = 24.0;
    register_label(storage, offset_x - 2.0, by + 10.0, "3", LABEL_SIZE, LABEL_COLOR);
    let trap_profile = [
        Point3::new(offset_x + 1.5, 0.0, by),
        Point3::new(offset_x + 4.0, 0.0, by),
        Point3::new(offset_x + 3.0, 0.0, by + 6.0),
        Point3::new(offset_x + 2.0, 0.0, by + 6.0),
    ];
    render_revolve(
        storage,
        &trap_profile,
        Point3::new(offset_x, 0.0, 0.0),
        Vector3::z(),
        GRAY,
        GREEN,
    );
}
