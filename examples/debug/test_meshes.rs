//! Mesh dispatcher for the debug viewer.
//!
//! Selects a pattern based on CLI arguments:
//!
//! - `patterns/` — algorithm output visualization (default)
//! - `test/` — hardcoded ground truth (selected with `--test`)
//!
//! ```text
//! cargo run --example debug                                # default (stroke_joins)
//! cargo run --example debug -- stroke_joins                # algorithm output
//! cargo run --example debug -- --test offset_intersection  # ground truth
//! ```

#[path = "patterns/mod.rs"]
mod patterns;

#[path = "test/mod.rs"]
mod test;

use revion_ui::MeshStorage;

/// Parsed CLI arguments.
struct CliArgs {
    /// Use test (ground truth) patterns instead of algorithm output.
    test_mode: bool,
    /// Pattern name to display.
    pattern: String,
}

/// Parse CLI arguments, extracting `--test` flag and pattern name.
fn parse_args() -> CliArgs {
    let args: Vec<String> = std::env::args().skip(1).collect();

    let test_mode = args.iter().any(|a| a == "--test");
    let pattern = args
        .iter()
        .find(|a| !a.starts_with('-'))
        .cloned()
        .unwrap_or_else(|| "stroke_joins".to_string());

    CliArgs { test_mode, pattern }
}

/// Register test meshes, selecting the pattern from the CLI arguments.
pub fn register_test_meshes(storage: &MeshStorage) {
    let args = parse_args();
    let name = &args.pattern;

    if args.test_mode {
        if test::register(storage, name) {
            return;
        }
        eprintln!("[debug] unknown test pattern: {name}");
        eprintln!("[debug] available (--test): {}", test::PATTERNS.join(", "));
    } else {
        if patterns::register(storage, name) {
            return;
        }
        eprintln!("[debug] unknown pattern: {name}");
        eprintln!("[debug] available: {}", patterns::PATTERNS.join(", "));
        eprintln!("[debug] available (--test): {}", test::PATTERNS.join(", "));
    }

    eprintln!("[debug] falling back to: stroke_joins");
    patterns::register(storage, "stroke_joins");
}
