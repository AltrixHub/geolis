//! Ground truth test patterns — hardcoded expected results for visual verification.
//!
//! Each module here draws hand-computed correct geometry (NOT algorithm output).
//! Compare these visually against algorithm-generated patterns in `patterns/`.

pub mod wall_offset;

use revion_ui::MeshStorage;

/// All available test pattern names.
pub const PATTERNS: &[&str] = &["wall_offset"];

/// Register meshes for the named test pattern. Returns `true` if found.
pub fn register(storage: &MeshStorage, name: &str) -> bool {
    match name {
        "wall_offset" => {
            wall_offset::register(storage);
            true
        }
        _ => false,
    }
}

// Re-export shared utilities from patterns for child modules.
#[allow(unused_imports)]
pub use super::patterns::{register_label, register_stroke};
