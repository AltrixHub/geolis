use geolis::math::{Point3, Vector3};
use geolis::operations::creation::{MakeFace, MakeWire};
use geolis::operations::shaping::Revolve;
use geolis::tessellation::{StrokeStyle, TessellateSolid, TessellationParams};
use geolis::topology::TopologyStore;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label, register_stroke, SceneBounds};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GRAY: Color = Color::rgb(180, 180, 180);
const GREEN: Color = Color::rgb(100, 200, 100);
const BLUE: Color = Color::rgb(100, 150, 255);
const RED: Color = Color::rgb(220, 80, 80);

/// Revolution axis (origin + direction) for [`render_revolve`].
struct RevolveAxis {
    origin: Point3,
    dir: Vector3,
}

/// Runs `MakeWire` -> `MakeFace` -> `Revolve` -> `TessellateSolid` and renders.
///
/// If `angle` is `Some`, creates a partial revolve; otherwise full 360°.
fn render_revolve(
    storage: &MeshStorage,
    bounds: &mut SceneBounds,
    points: &[Point3],
    axis: &RevolveAxis,
    angle: Option<f64>,
    outline_color: Color,
    mesh_color: Color,
) {
    // Draw profile outline
    if let Ok(style) = StrokeStyle::new(0.05) {
        register_stroke(storage, bounds, points, style, true, outline_color);
    }

    // Draw axis as a thin line
    if let Ok(style) = StrokeStyle::new(0.02) {
        let axis_len = 8.0;
        let axis_line = [
            axis.origin - axis.dir * (axis_len / 2.0),
            axis.origin + axis.dir * (axis_len / 2.0),
        ];
        register_stroke(storage, bounds, &axis_line, style, false, RED);
    }

    let mut topo = TopologyStore::new();
    let Ok(wire) = MakeWire::new(points.to_vec(), true).execute(&mut topo) else {
        return;
    };
    let Ok(face) = MakeFace::new(wire, vec![]).execute(&mut topo) else {
        return;
    };
    let mut revolve = Revolve::new(face, axis.origin, axis.dir);
    if let Some(a) = angle {
        revolve = revolve.with_angle(a);
    }
    let Ok(solid) = revolve.execute(&mut topo) else {
        return;
    };
    if let Ok(mesh) = TessellateSolid::new(solid, TessellationParams::default()).execute(&topo) {
        register_face(storage, bounds, mesh, mesh_color);
    }

    if let Ok(solid_data) = topo.solid(solid) {
        register_edges(
            storage,
            bounds,
            &topo,
            solid_data.outer_shell,
            outline_color,
        );
    }
}

/// Registers a revolve case at the given axis position.
///
/// - Profile is in the plane `y = ay`, spanning `z = 0..H`.
/// - Axis of revolution at `(ax, ay, 0)` along Z.
#[allow(clippy::too_many_arguments)]
fn register_case(
    storage: &MeshStorage,
    bounds: &mut SceneBounds,
    label: &str,
    ax: f64,
    ay: f64,
    profile: &[Point3],
    angle: Option<f64>,
    mesh_color: Color,
) {
    // Label positioned above-left of the revolve in the 2D (XY) projection.
    register_label(
        storage,
        bounds,
        ax - 5.0,
        ay + 6.0,
        label,
        LABEL_SIZE,
        LABEL_COLOR,
    );
    render_revolve(
        storage,
        bounds,
        profile,
        &RevolveAxis {
            origin: Point3::new(ax, ay, 0.0),
            dir: Vector3::z(),
        },
        angle,
        GRAY,
        mesh_color,
    );
}

#[allow(clippy::too_many_lines)]
pub fn register(storage: &MeshStorage, bounds: &mut SceneBounds) {
    // Profile height (z = 0..H) shared by all cases.
    let h = 6.0;

    // Column X positions (spacing > diameter of largest shape).
    let col = [0.0, 14.0, 28.0];

    // ── Row 1 (y = 0): Full revolve (360°) ─────────────────────

    let y1 = 0.0;

    // Case 1: Square profile → hollow cylinder
    let sq = [
        Point3::new(col[0] + 2.0, y1, 0.0),
        Point3::new(col[0] + 4.0, y1, 0.0),
        Point3::new(col[0] + 4.0, y1, h),
        Point3::new(col[0] + 2.0, y1, h),
    ];
    register_case(storage, bounds, "1", col[0], y1, &sq, None, GREEN);

    // Case 2: Triangle with vertex on axis → cone
    let tri = [
        Point3::new(col[1], y1, h),
        Point3::new(col[1] + 3.0, y1, 0.0),
        Point3::new(col[1] + 3.0, y1, h),
    ];
    register_case(storage, bounds, "2", col[1], y1, &tri, None, BLUE);

    // Case 3: Trapezoid → truncated cone (frustum)
    let trap = [
        Point3::new(col[2] + 1.5, y1, 0.0),
        Point3::new(col[2] + 4.0, y1, 0.0),
        Point3::new(col[2] + 3.0, y1, h),
        Point3::new(col[2] + 2.0, y1, h),
    ];
    register_case(storage, bounds, "3", col[2], y1, &trap, None, GREEN);

    // ── Row 2 (y = -14): Partial revolve ────────────────────────

    let y2 = -14.0;

    // Case 4: Square 90°
    let sq4 = [
        Point3::new(col[0] + 2.0, y2, 0.0),
        Point3::new(col[0] + 4.0, y2, 0.0),
        Point3::new(col[0] + 4.0, y2, h),
        Point3::new(col[0] + 2.0, y2, h),
    ];
    register_case(
        storage,
        bounds,
        "4",
        col[0],
        y2,
        &sq4,
        Some(std::f64::consts::FRAC_PI_2),
        GREEN,
    );

    // Case 5: Triangle on-axis 180°
    let tri5 = [
        Point3::new(col[1], y2, h),
        Point3::new(col[1] + 3.0, y2, 0.0),
        Point3::new(col[1] + 3.0, y2, h),
    ];
    register_case(
        storage,
        bounds,
        "5",
        col[1],
        y2,
        &tri5,
        Some(std::f64::consts::PI),
        BLUE,
    );

    // Case 6: Trapezoid 270°
    let trap6 = [
        Point3::new(col[2] + 1.5, y2, 0.0),
        Point3::new(col[2] + 4.0, y2, 0.0),
        Point3::new(col[2] + 3.0, y2, h),
        Point3::new(col[2] + 2.0, y2, h),
    ];
    register_case(
        storage,
        bounds,
        "6",
        col[2],
        y2,
        &trap6,
        Some(3.0 * std::f64::consts::FRAC_PI_2),
        GREEN,
    );
}
