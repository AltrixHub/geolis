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
use viewer::AppState;

fn main() -> Result<(), RevionError> {
    // Default: WARN for everything, INFO for geolis.
    // Override with RUST_LOG env var (e.g. RUST_LOG=revion_renderer=debug).
    let env_filter = tracing_subscriber::EnvFilter::from_default_env()
        .add_directive(tracing_subscriber::filter::LevelFilter::WARN.into())
        .add_directive("debug=info".parse().unwrap_or_default())
        .add_directive("geolis=info".parse().unwrap_or_default());
    tracing_subscriber::fmt().with_env_filter(env_filter).init();

    // Build mesh storage with test data
    let storage = MeshStorage::new();
    test_meshes::register_test_meshes(&storage);

    let state = AppState {
        mesh_storage: storage,
    };

    let mut app = App::new("Geolis Debug Viewer")?;
    app.context_mut().provide_store(state)?;
    app.build_with_component(viewer::app_component)?;
    app.run()
}
