//! NURBS curve tessellation showcase.
//!
//! Renders adaptively tessellated NURBS curves as polylines:
//! 1. an interpolated free-form 3D curve through scattered points,
//! 2. an exact rational unit circle, and
//! 3. a rational quarter arc.

use geolis::geometry::nurbs::NurbsCurve3D;
use geolis::math::{Point3, Vector3};
use geolis::tessellation::{tessellate_nurbs_curve, CurveTessellationOptions, StrokeStyle};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_label, register_stroke};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GREEN: Color = Color::rgb(100, 220, 100);
const BLUE: Color = Color::rgb(100, 150, 255);
const ORANGE: Color = Color::rgb(255, 170, 80);

fn stroke_style() -> StrokeStyle {
    StrokeStyle::new(0.05).unwrap_or_else(|_| unreachable!())
}

/// Tessellate a NURBS curve and register it as a translated stroke.
fn register_curve(storage: &MeshStorage, curve: &NurbsCurve3D, bx: f64, by: f64, color: Color) {
    let options = CurveTessellationOptions::default();
    let Ok(points) = tessellate_nurbs_curve(curve, &options) else {
        return;
    };
    let translated: Vec<Point3> = points
        .iter()
        .map(|p| Point3::new(p.x + bx, p.y + by, p.z))
        .collect();
    register_stroke(storage, &translated, stroke_style(), false, color);
}

pub fn register(storage: &MeshStorage) {
    // Case 1: interpolated free-form 3D curve.
    {
        let bx = 0.0;
        let by = 0.0;
        register_label(storage, bx - 1.5, by + 6.0, "1", LABEL_SIZE, LABEL_COLOR);

        let waypoints = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 3.0, 1.0),
            Point3::new(5.0, 1.0, 2.0),
            Point3::new(7.0, 4.0, 0.0),
            Point3::new(9.0, 0.0, 1.0),
        ];
        if let Ok((curve, _params)) = NurbsCurve3D::interpolate(&waypoints, 3) {
            register_curve(storage, &curve, bx, by, GREEN);
        }
    }

    // Case 2: exact rational unit circle (radius 2 for visibility).
    {
        let bx = 0.0;
        let by = -10.0;
        register_label(storage, bx - 1.5, by + 4.0, "2", LABEL_SIZE, LABEL_COLOR);

        if let Ok(curve) = NurbsCurve3D::circle(
            Point3::new(bx, by, 0.0),
            2.0,
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
        ) {
            // The circle is already centered at (bx, by); render without an
            // extra offset.
            register_curve(storage, &curve, 0.0, 0.0, BLUE);
        }
    }

    // Case 3: rational quarter arc.
    {
        let bx = 8.0;
        let by = -10.0;
        register_label(storage, bx - 1.5, by + 4.0, "3", LABEL_SIZE, LABEL_COLOR);

        if let Ok(curve) = NurbsCurve3D::arc(
            Point3::new(bx, by, 0.0),
            3.0,
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
            0.0,
            std::f64::consts::FRAC_PI_2,
        ) {
            register_curve(storage, &curve, 0.0, 0.0, ORANGE);
        }
    }
}
