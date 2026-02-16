# Debug Viewer

## Overview

`examples/debug/` is a debug viewer for visualising Geolis meshes.

## Development Workflow

The debug viewer is the core tool for a **visual TDD** workflow. Follow these 3 steps.

### Step 1: Build ground truth with `--test`

Hardcode hand-computed correct geometry from documentation or manual calculation, then visualize with `--test` mode to verify.

```bash
# Visualize ground truth
cargo run --example debug -- --test offset_intersection
```

- Hardcode hand-computed vertex coordinates in `test/`
- Visually confirm correctness in the viewer
- Color convention: Gray = original shape, Green = expected inward, Blue = expected outward

### Step 2: Create test cases from ground truth and iterate on the algorithm

Translate the verified ground truth from Step 1 into `#[test]` cases and fix the algorithm until tests pass.

```rust
#[test]
fn test_t_shape_inward_offset() {
    let result = PolylineOffset2D::new(t_shape(), 0.3, true).execute().unwrap();
    // Compare against ground truth coordinates from Step 1
    assert_approx_eq(&result, &expected_inward);
}
```

- Verify algorithm correctness automatically with `cargo test`
- Iterate on the algorithm until all tests pass

### Step 3: Visualize algorithm output in default mode

Once tests pass, run the viewer without `--test` to visualize the actual algorithm output and confirm overall appearance.

```bash
# Visualize actual algorithm output
cargo run --example debug -- offset_intersection
```

- Code in `patterns/` calls the algorithm and renders results
- Visually confirm that algorithm output matches expected results for the same input shapes
- Catch edge cases and subtle numerical precision issues here

### Workflow Summary

```
 +----------------------------------------------------------+
 | Step 1: Visualize ground truth with --test                |
 |   Hardcode in test/ -> visual confirmation in viewer      |
 +----------------------------+-----------------------------+
                              |
                              v
 +----------------------------------------------------------+
 | Step 2: Ground truth -> #[test] -> iterate on algorithm   |
 |   Loop with cargo test                                    |
 +----------------------------+-----------------------------+
                              |
                              v
 +----------------------------------------------------------+
 | Step 3: Visualize algorithm output in default mode        |
 |   Render actual output in patterns/ -> final visual check |
 +----------------------------------------------------------+
```

## File Structure

```
examples/debug/
+-- main.rs              -- Entry point (DO NOT EDIT)
+-- viewer.rs            -- Pure viewer UI (DO NOT EDIT)
+-- test_meshes.rs       -- CLI dispatcher (DO NOT EDIT)
+-- patterns/            -- Algorithm output visualization (Step 3)
|   +-- mod.rs           -- Pattern registry + shared utilities (selectively editable)
|   +-- stroke_joins.rs  -- LineJoin comparison pattern
|   +-- basic_strokes.rs -- Basic stroke shapes pattern
+-- test/                -- Ground truth visualization (Step 1)
    +-- mod.rs           -- Test pattern registry (selectively editable)
    +-- *.rs             -- Ground truth pattern files
```

| File | Role | Editable? |
|------|------|-----------|
| `main.rs` | Logging init, app startup | No |
| `viewer.rs` | `AppState`, viewport layout, components | No |
| `test_meshes.rs` | CLI arg -> pattern dispatch | No |
| `patterns/mod.rs` | Pattern registry + shared conversion utilities | Add new patterns here |
| `patterns/*.rs` | Algorithm output patterns | **Yes** |
| `test/mod.rs` | Test pattern registry + re-exports | Add new test patterns here |
| `test/*.rs` | Ground truth patterns | **Yes** |

## How It Works

1. `test_meshes.rs` parses CLI arguments (`--test` flag + pattern name)
2. Without `--test`: `patterns::register(storage, name)` dispatches to algorithm output
3. With `--test`: `test::register(storage, name)` dispatches to ground truth data
4. Each pattern file calls shared utilities to register meshes into `MeshStorage`
5. The viewer displays them in 2D (left) and 3D (right) viewports

## Pattern System

### Running Patterns

```bash
# Default (stroke_joins)
cargo run --example debug

# Algorithm output pattern (Step 3)
cargo run --example debug -- stroke_joins
cargo run --example debug -- offset_intersection

# Ground truth (Step 1)
cargo run --example debug -- --test offset_intersection

# Unknown pattern -> prints available list, falls back to stroke_joins
cargo run --example debug -- foo
```

### Available Patterns

| Pattern | Mode | Description |
|---------|------|-------------|
| `stroke_joins` | default | `LineJoin` comparison -- Miter / Auto / Bevel columns |
| `basic_strokes` | default | Simple stroke shapes (line, L-shape, triangle, curve, square) |
| `polyline_offset` | default | Polyline offset algorithm results |
| `offset_intersection` | default | Offset self-intersection algorithm results |
| `offset_intersection` | `--test` | Hand-computed ground truth for offset self-intersection |

### Adding a New Pattern

1. Create `examples/debug/patterns/new_pattern.rs`:

```rust
use revion_ui::MeshStorage;

use super::register_stroke; // or other shared utilities

pub fn register(storage: &MeshStorage) {
    // Register meshes here
}
```

2. Update `examples/debug/patterns/mod.rs`:

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

### Adding a New Test (Ground Truth) Pattern

1. Create `examples/debug/test/new_pattern.rs`:

```rust
use revion_ui::MeshStorage;

use super::register_stroke; // re-exported from patterns

pub fn register(storage: &MeshStorage) {
    // Register hardcoded expected meshes here
}
```

2. Update `examples/debug/test/mod.rs`:

```rust
pub mod new_pattern;

pub const PATTERNS: &[&str] = &["offset_intersection", "new_pattern"];

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
| `into_raw_mesh_2d(mesh, color)` | Convert Geolis `TriangleMesh` -> Revion `RawMesh2D` (XY projection) |
| `into_raw_mesh_3d(mesh, color)` | Convert Geolis `TriangleMesh` -> Revion `RawMesh3D` |
| `register_stroke(storage, points, style, closed, color)` | Tessellate + register both 2D and 3D |

Test patterns access these via re-export in `test/mod.rs`.

## Running

```bash
# Default pattern
cargo run --example debug

# Specific pattern
cargo run --example debug -- basic_strokes

# Ground truth
cargo run --example debug -- --test offset_intersection

# With Revion renderer debug logs
RUST_LOG=revion_renderer=debug cargo run --example debug -- stroke_joins

# Verbose -- all crates
RUST_LOG=debug cargo run --example debug
```

## Logging

The viewer uses `tracing` with `EnvFilter`. Defaults:

| Target | Level |
|--------|-------|
| `geolis` | INFO |
| `debug` | INFO |
| Everything else (revion, wgpu, ...) | WARN |

Override via `RUST_LOG` environment variable.
