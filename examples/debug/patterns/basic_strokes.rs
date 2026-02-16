//! Basic stroke shapes — simple polylines for quick visual testing.

use geolis::math::Point3;
use geolis::tessellation::StrokeStyle;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::register_stroke;

/// Register basic stroke meshes.
pub fn register(storage: &MeshStorage) {
    let strokes: &[(&[Point3], f64, bool, Color)] = &[
        // Straight line
        (
            &[
                Point3::new(-4.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
            ],
            0.3,
            false,
            Color::rgb(255, 100, 100),
        ),
        // L-shape (90°)
        (
            &[
                Point3::new(-4.0, 2.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
                Point3::new(0.0, 5.0, 0.0),
            ],
            0.25,
            false,
            Color::rgb(100, 200, 255),
        ),
        // Closed triangle
        (
            &[
                Point3::new(1.0, 2.0, 0.0),
                Point3::new(4.0, 2.0, 0.0),
                Point3::new(2.5, 5.0, 0.0),
            ],
            0.2,
            true,
            Color::rgb(100, 255, 150),
        ),
        // Smooth curve approximation
        (
            &[
                Point3::new(-4.0, -3.0, 0.0),
                Point3::new(-2.0, -1.5, 0.0),
                Point3::new(0.0, -2.5, 0.0),
                Point3::new(2.0, -1.5, 0.0),
                Point3::new(4.0, -3.0, 0.0),
            ],
            0.2,
            false,
            Color::rgb(255, 200, 80),
        ),
        // Closed square
        (
            &[
                Point3::new(-2.0, -6.0, 0.0),
                Point3::new(2.0, -6.0, 0.0),
                Point3::new(2.0, -4.0, 0.0),
                Point3::new(-2.0, -4.0, 0.0),
            ],
            0.15,
            true,
            Color::rgb(200, 130, 255),
        ),
    ];

    for &(points, width, closed, color) in strokes {
        if let Ok(style) = StrokeStyle::new(width) {
            register_stroke(storage, points, style, closed, color);
        }
    }
}
