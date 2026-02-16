//! `LineJoin` comparison — same polylines rendered with Miter / Auto / Bevel.

use geolis::math::Point3;
use geolis::tessellation::{LineJoin, StrokeStyle};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_label, register_stroke};

/// Label size and color.
const LABEL_SIZE: f64 = 0.8;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);

/// Register stroke-join comparison meshes.
#[allow(clippy::too_many_lines)]
pub fn register(storage: &MeshStorage) {
    // Polyline definitions: (base points, width, closed)
    let shapes: &[(&[Point3], f64, bool)] = &[
        // Zigzag — sharp (~30°) angles
        (
            &[
                Point3::new(-2.0, -3.0, 0.0),
                Point3::new(-1.0, 0.0, 0.0),
                Point3::new(0.0, -3.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, -3.0, 0.0),
            ],
            0.25,
            false,
        ),
        // Hairpin (~20°) — near reversal
        (
            &[
                Point3::new(-1.0, 1.0, 0.0),
                Point3::new(0.0, 4.0, 0.0),
                Point3::new(-0.8, 1.2, 0.0),
            ],
            0.25,
            false,
        ),
        // Right-angle staircase (90° turns)
        (
            &[
                Point3::new(-1.5, 5.0, 0.0),
                Point3::new(0.0, 5.0, 0.0),
                Point3::new(0.0, 6.5, 0.0),
                Point3::new(1.5, 6.5, 0.0),
                Point3::new(1.5, 8.0, 0.0),
            ],
            0.2,
            false,
        ),
        // Star (closed) — very sharp tips (~36°)
        (
            &[
                Point3::new(0.0, 13.5, 0.0),
                Point3::new(0.6, 11.5, 0.0),
                Point3::new(2.5, 11.5, 0.0),
                Point3::new(1.0, 10.5, 0.0),
                Point3::new(1.6, 8.5, 0.0),
                Point3::new(0.0, 9.5, 0.0),
                Point3::new(-1.6, 8.5, 0.0),
                Point3::new(-1.0, 10.5, 0.0),
                Point3::new(-2.5, 11.5, 0.0),
                Point3::new(-0.6, 11.5, 0.0),
            ],
            0.12,
            true,
        ),
    ];

    // Three columns: (x-offset, LineJoin, color palette)
    let columns: &[(f64, LineJoin, [Color; 4])] = &[
        (
            -6.0,
            LineJoin::Miter,
            [
                Color::rgb(255, 100, 100),
                Color::rgb(255, 140, 80),
                Color::rgb(230, 80, 80),
                Color::rgb(255, 120, 120),
            ],
        ),
        (
            0.0,
            LineJoin::Auto,
            [
                Color::rgb(100, 220, 130),
                Color::rgb(140, 255, 100),
                Color::rgb(80, 200, 100),
                Color::rgb(120, 240, 150),
            ],
        ),
        (
            6.0,
            LineJoin::Bevel,
            [
                Color::rgb(100, 150, 255),
                Color::rgb(80, 180, 255),
                Color::rgb(100, 100, 230),
                Color::rgb(130, 170, 255),
            ],
        ),
    ];

    // Case labels (one per shape row, to the left of the Miter column)
    let label_y: &[f64] = &[0.5, 4.5, 8.5, 14.0];
    for (i, &ly) in label_y.iter().enumerate() {
        register_label(storage, -10.0, ly, &format!("{}", i + 1), LABEL_SIZE, LABEL_COLOR);
    }

    for &(x_off, join, ref colors) in columns {
        for (idx, &(points, width, closed)) in shapes.iter().enumerate() {
            let shifted: Vec<Point3> = points
                .iter()
                .map(|p| Point3::new(p.x + x_off, p.y, p.z))
                .collect();

            if let Ok(style) = StrokeStyle::new(width) {
                let style = style.with_line_join(join);
                register_stroke(storage, &shifted, style, closed, colors[idx]);
            }
        }
    }
}
