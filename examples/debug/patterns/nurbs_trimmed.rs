//! Trimmed NURBS face showcase.
//!
//! Builds a wavy interpolated NURBS patch, trims it with a full-domain outer
//! rectangle and a circular hole, and renders the result through the full
//! `MakeNurbsFace` + face-tessellation pipeline (constrained Delaunay), so the
//! trimming path is exercised end-to-end rather than via a tessellation
//! shortcut.

use geolis::geometry::nurbs::{KnotVector, NurbsCurve2D, NurbsSurface};
use geolis::math::{Point2, Point3};
use geolis::operations::creation::MakeNurbsFace;
use geolis::tessellation::{TessellateFace, TessellationParams};
use geolis::topology::{FaceTrim, TopologyStore, TrimLoop};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_face, register_label};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const TEAL: Color = Color::rgb(90, 200, 200);

/// A wavy bicubic patch over the unit parameter square, spanning [0,6]x[0,6] in
/// XY with a gentle central rise.
fn wavy_patch() -> Option<NurbsSurface> {
    let heights = [
        [0.0, 0.5, 0.5, 0.0],
        [0.5, 2.0, 2.0, 0.5],
        [0.5, 2.0, 2.0, 0.5],
        [0.0, 0.5, 0.5, 0.0],
    ];
    let mut control = Vec::with_capacity(16);
    for (i, row) in heights.iter().enumerate() {
        for (j, &h) in row.iter().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            control.push(Point3::new(i as f64 * 2.0, j as f64 * 2.0, h));
        }
    }
    let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]).ok()?;
    NurbsSurface::from_unweighted(control, 4, 4, knots.clone(), knots, 3, 3).ok()
}

/// Counter-clockwise outer rectangle covering the full [0,1]^2 parameter domain.
fn full_domain_outer() -> TrimLoop {
    let line = |a: (f64, f64), b: (f64, f64)| {
        NurbsCurve2D::from_unweighted(
            vec![Point2::new(a.0, a.1), Point2::new(b.0, b.1)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).expect("valid degree-1 knots"),
            1,
        )
        .expect("valid degree-1 line")
    };
    TrimLoop::new(vec![
        line((0.0, 0.0), (1.0, 0.0)),
        line((1.0, 0.0), (1.0, 1.0)),
        line((1.0, 1.0), (0.0, 1.0)),
        line((0.0, 1.0), (0.0, 0.0)),
    ])
}

/// Clockwise circular hole in UV (rational circle reversed to wind clockwise).
fn circular_hole() -> Option<TrimLoop> {
    let circle = NurbsCurve2D::circle_uv(Point2::new(0.5, 0.5), 0.22).ok()?;
    // `circle_uv` winds counter-clockwise; reverse it so the hole winds CW.
    let cw = circle.reverse().ok()?;
    Some(TrimLoop::new(vec![cw]))
}

pub fn register(storage: &MeshStorage) {
    register_label(storage, -1.5, 8.0, "1", LABEL_SIZE, LABEL_COLOR);

    let Some(surface) = wavy_patch() else {
        return;
    };
    let Some(hole) = circular_hole() else {
        return;
    };
    let trim = FaceTrim::new(full_domain_outer(), vec![hole]);

    let mut store = TopologyStore::new();
    let Ok(face) = MakeNurbsFace::new(surface)
        .with_trim(trim)
        .execute(&mut store)
    else {
        return;
    };
    let Ok(mesh) = TessellateFace::new(face, TessellationParams::default()).execute(&store) else {
        return;
    };
    register_face(storage, mesh, TEAL);
}
