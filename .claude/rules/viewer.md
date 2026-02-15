# Viewer (Debug Tool)

## Overview

`examples/viewer/` is a debug viewer for visualising Geolis meshes.
Every `/dev-viewer` invocation **rewrites** a pattern file under `patterns/`, then runs the viewer.

## File Structure

```
examples/viewer/
├── main.rs              — Entry point (DO NOT EDIT)
├── viewer.rs            — Pure viewer UI (DO NOT EDIT)
├── test_meshes.rs       — CLI dispatcher (DO NOT EDIT)
└── patterns/
    ├── mod.rs           — Pattern registry + shared utilities (selectively editable)
    ├── stroke_joins.rs  — LineJoin comparison pattern
    └── basic_strokes.rs — Basic stroke shapes pattern
```

| File | Role | Editable? |
|------|------|-----------|
| `main.rs` | Logging init, app startup | No |
| `viewer.rs` | `AppState`, viewport layout, components | No |
| `test_meshes.rs` | CLI arg → pattern dispatch | No |
| `patterns/mod.rs` | Pattern registry + shared conversion utilities | Add new patterns here |
| `patterns/*.rs` | Individual pattern files | **Yes** (rewrite per `/dev-viewer`) |

## How It Works

1. `test_meshes.rs` reads the first CLI argument to select a pattern
2. `patterns::register(storage, name)` dispatches to the selected pattern file
3. Each pattern file calls shared utilities to register meshes into `MeshStorage`
4. The viewer displays them in 2D (left) and 3D (right) viewports

## Pattern System

### Running a Specific Pattern

```bash
# Default (stroke_joins)
cargo run --example viewer

# Specify pattern
cargo run --example viewer -- stroke_joins
cargo run --example viewer -- basic_strokes

# Unknown pattern → prints available list, falls back to stroke_joins
cargo run --example viewer -- foo
```

### Available Patterns

| Pattern | Description |
|---------|-------------|
| `stroke_joins` | `LineJoin` comparison — Miter / Auto / Bevel columns |
| `basic_strokes` | Simple stroke shapes (line, L-shape, triangle, curve, square) |

### Adding a New Pattern

1. Create `examples/viewer/patterns/new_pattern.rs`:

```rust
use revion_ui::MeshStorage;

use super::register_stroke; // or other shared utilities

pub fn register(storage: &MeshStorage) {
    // Register meshes here
}
```

2. Update `examples/viewer/patterns/mod.rs`:

```rust
// Add module declaration
pub mod new_pattern;

// Add to PATTERNS list
pub const PATTERNS: &[&str] = &["stroke_joins", "basic_strokes", "new_pattern"];

// Add match arm in register()
pub fn register(storage: &MeshStorage, name: &str) -> bool {
    match name {
        // ...existing...
        "new_pattern" => { new_pattern::register(storage); true }
        _ => false,
    }
}
```

### Shared Utilities (in `patterns/mod.rs`)

| Function | Description |
|----------|-------------|
| `into_raw_mesh_2d(mesh, color)` | Convert Geolis `TriangleMesh` → Revion `RawMesh2D` (XY projection) |
| `into_raw_mesh_3d(mesh, color)` | Convert Geolis `TriangleMesh` → Revion `RawMesh3D` |
| `register_stroke(storage, points, style, closed, color)` | Tessellate + register both 2D and 3D |

### `/dev-viewer` Behavior

When `/dev-viewer` is invoked:

- **Rewrite** the relevant pattern file under `patterns/` (or create a new one)
- **Do NOT** modify `main.rs`, `viewer.rs`, or `test_meshes.rs`
- **Update** `patterns/mod.rs` if adding a new pattern (add `pub mod`, `PATTERNS` entry, match arm)
- Run the viewer with `cargo run --example viewer -- <pattern_name>`

## Running

```bash
# Default pattern
cargo run --example viewer

# Specific pattern
cargo run --example viewer -- basic_strokes

# With Revion renderer debug logs
RUST_LOG=revion_renderer=debug cargo run --example viewer -- stroke_joins

# Verbose — all crates
RUST_LOG=debug cargo run --example viewer
```

## Logging

The viewer uses `tracing` with `EnvFilter`. Defaults:

| Target | Level |
|--------|-------|
| `geolis` | INFO |
| `viewer` | INFO |
| Everything else (revion, wgpu, ...) | WARN |

Override via `RUST_LOG` environment variable.
