//! Basic stroke shapes — simple polylines for quick visual testing.

use geolis::math::Point3;
use geolis::tessellation::StrokeStyle;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_label, register_stroke};

/// Label size and color.
const LABEL_SIZE: f64 = 0.8;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);

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

    // Case labels — positioned to the left or above each shape
    let label_pos: &[(f64, f64)] = &[
        (-6.0, 0.5),   // 1: Straight line
        (-6.0, 5.5),   // 2: L-shape
        (5.0, 5.5),    // 3: Closed triangle
        (-6.0, -1.0),  // 4: Smooth curve
        (-4.5, -3.5),  // 5: Closed square
    ];
    for (i, &(lx, ly)) in label_pos.iter().enumerate() {
        register_label(storage, lx, ly, &format!("{}", i + 1), LABEL_SIZE, LABEL_COLOR);
    }

    for &(points, width, closed, color) in strokes {
        if let Ok(style) = StrokeStyle::new(width) {
            register_stroke(storage, points, style, closed, color);
        }
    }
}
