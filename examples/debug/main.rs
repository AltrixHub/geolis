//! Geolis Debug Viewer — minimal mesh viewer for `RawMesh2D` / `RawMesh3D`.
//!
//! ```text
//! main.rs          — entry point (this file)
//! viewer.rs        — pure viewer UI (AppState + components)
//! test_meshes.rs   — mesh dispatcher (patterns / test ground truth)
//! ```
//!
//! Usage:
//! ```text
//! cargo run --example debug                                # default (stroke_joins)
//! cargo run --example debug -- stroke_joins                # algorithm output
//! cargo run --example debug -- --test offset_intersection  # ground truth
//! ```
//!
//! Controls:
//! - 2D: Space+drag to pan, Cmd+scroll to zoom
//! - 3D: Right-click+drag to orbit, Middle-click+drag to pan, Scroll to zoom

mod test_meshes;
mod viewer;

use revion_app::App;
use revion_ui::{MeshStorage, RevionError};
use test_meshes::SceneBounds;

/// Derives an initial 3D camera pose framing the registered scene bounds.
///
/// Looks along revion's default diagonal (offset direction `(0, -2, 1)`,
/// Z-up) at the bounds center, from far enough away that the whole bounding
/// box fits comfortably in a 45-degree field of view.
fn camera_from_bounds(bounds: &SceneBounds) -> Option<([f32; 3], [f32; 3])> {
    if bounds.is_empty() {
        return None;
    }
    let target = bounds.center();
    let distance = (bounds.diagonal() * 1.2).max(5.0);
    let dir_len = (0.0_f64 * 0.0 + 2.0 * 2.0 + 1.0).sqrt();
    let dir = [0.0, -2.0 / dir_len, 1.0 / dir_len];
    #[allow(clippy::cast_possible_truncation)]
    let eye = [
        (target[0] + dir[0] * distance) as f32,
        (target[1] + dir[1] * distance) as f32,
        (target[2] + dir[2] * distance) as f32,
    ];
    #[allow(clippy::cast_possible_truncation)]
    let target = [target[0] as f32, target[1] as f32, target[2] as f32];
    Some((eye, target))
}

fn main() -> Result<(), RevionError> {
    // Default: WARN for everything, INFO for geolis.
    // Override with RUST_LOG env var (e.g. RUST_LOG=revion_renderer=debug).
    let env_filter = tracing_subscriber::EnvFilter::from_default_env()
        .add_directive(tracing_subscriber::filter::LevelFilter::WARN.into())
        .add_directive("debug=info".parse().unwrap_or_default())
        .add_directive("geolis=info".parse().unwrap_or_default());
    tracing_subscriber::fmt().with_env_filter(env_filter).init();

    // Build mesh storage with test data. `MeshStorage` is internally
    // reference-counted, so the clone captured by the component closure
    // shares the same backing storage.
    let storage = MeshStorage::new();
    let camera = test_meshes::register_test_meshes(&storage)
        .as_ref()
        .and_then(camera_from_bounds);

    let mut app = App::new("Geolis Debug Viewer")?;
    app.build_with_component(move |ctx| viewer::app_component(ctx, &storage, camera))?;
    app.run()
}
