//! Test mesh dispatcher for the viewer example.
//!
//! Selects a pattern based on the first CLI argument.
//!
//! ```text
//! cargo run --example viewer                    # default (stroke_joins)
//! cargo run --example viewer -- stroke_joins    # LineJoin comparison
//! cargo run --example viewer -- basic_strokes   # basic shapes
//! ```

#[path = "patterns/mod.rs"]
mod patterns;

use revion_ui::MeshStorage;

/// Register test meshes, selecting the pattern from the CLI argument.
pub fn register_test_meshes(storage: &MeshStorage) {
    let name = std::env::args().nth(1);
    let name = name.as_deref().unwrap_or("stroke_joins");

    if !patterns::register(storage, name) {
        eprintln!("[viewer] unknown pattern: {name}");
        eprintln!("[viewer] available: {}", patterns::PATTERNS.join(", "));
        eprintln!("[viewer] falling back to: stroke_joins");
        patterns::register(storage, "stroke_joins");
    }
}
