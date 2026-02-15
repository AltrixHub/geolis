# Dev Viewer

Visualise meshes in the debug viewer. Each invocation rewrites a pattern file under `patterns/` to display the requested content.

## Input

$ARGUMENTS

**Arguments (required):**

| Argument | Description |
|----------|-------------|
| `{description}` | What to display — Geolis API calls, shapes, tessellation output, etc. |

**Flags (optional):**

| Flag | Description |
|------|-------------|
| `--verbose` | Enable Revion renderer debug logs |
| `--pattern <name>` | Target a specific pattern file (default: creates/rewrites as needed) |

## Instructions

### 1. Understand the Request

Parse the description to determine what meshes to generate. Common patterns:

| Pattern | Approach |
|---------|----------|
| Geolis API output | Call Geolis tessellation/geometry APIs, convert with shared utilities |
| Raw shapes | Use `revion_core::tessellation::{Circle, Sphere, Rectangle}` helpers |
| Specific data | Build `RawMesh2D`/`RawMesh3D` directly with vertices and indices |

### 2. Choose or Create a Pattern File

1. If `--pattern <name>` is given, target that pattern file
2. Otherwise, decide whether to rewrite an existing pattern or create a new one
3. For quick one-off tests, rewrite an existing pattern (e.g., `basic_strokes.rs`)
4. For new test categories, create a new pattern file

### 3. Write the Pattern File

Rewrite (or create) the pattern file under `examples/viewer/patterns/`:

```rust
//! Description of what this pattern visualises.

use geolis::math::Point3;
use geolis::tessellation::{StrokeStyle, TessellateStroke};
use revion_ui::MeshStorage;

use super::register_stroke; // shared utility

pub fn register(storage: &MeshStorage) {
    // Generate and register meshes
}
```

**Available shared utilities** (from `patterns/mod.rs`):

| Function | Description |
|----------|-------------|
| `into_raw_mesh_2d(mesh, color)` | Convert Geolis `TriangleMesh` → Revion `RawMesh2D` |
| `into_raw_mesh_3d(mesh, color)` | Convert Geolis `TriangleMesh` → Revion `RawMesh3D` |
| `register_stroke(storage, points, style, closed, color)` | Tessellate + register both 2D and 3D |

### 4. Register the Pattern (if new)

If creating a new pattern, update `examples/viewer/patterns/mod.rs`:

1. Add `pub mod new_pattern;`
2. Add name to `PATTERNS` array
3. Add match arm in `register()` function

### 5. Build Check

Verify the example compiles:

```bash
cargo check --example viewer
```

If it fails, fix the errors in the pattern file before proceeding.

### 6. Run the Viewer

```bash
# Run with the target pattern
cargo run --example viewer -- <pattern_name>

# Verbose — include Revion renderer logs
RUST_LOG=revion_renderer=debug cargo run --example viewer -- <pattern_name>
```

Use `--verbose` flag to decide which command to run.

### 7. Report

After launching, display:

```
Viewer launched: {brief description of what is displayed}
Pattern: <pattern_name>
- 2D viewport (left): Space+drag pan, Cmd+scroll zoom
- 3D viewport (right): Right-drag orbit, Middle-drag pan, Scroll zoom
```

**Rules:**
- Do NOT modify `main.rs`, `viewer.rs`, or `test_meshes.rs`
- Only modify files under `patterns/`
- Use shared utilities from `patterns/mod.rs` for mesh conversion
- Each `/dev-viewer` invocation replaces the target pattern file content

## File Structure

```
examples/viewer/
├── main.rs              — Entry point (DO NOT EDIT)
├── viewer.rs            — Pure viewer UI (DO NOT EDIT)
├── test_meshes.rs       — CLI dispatcher (DO NOT EDIT)
└── patterns/
    ├── mod.rs           — Pattern registry + shared utilities
    ├── stroke_joins.rs  — LineJoin comparison
    └── basic_strokes.rs — Basic stroke shapes
```

## Examples

```bash
# Visualise tessellation output from Geolis
/dev-viewer tessellate_stroke で L字ポリラインをテッセレートした結果

# Visualise raw geometry
/dev-viewer 半径1.0の円と半径0.5の円を重ねて表示

# Target a specific pattern file
/dev-viewer --pattern stroke_joins ジグザグとヘアピンのLineJoin比較

# With verbose logs
/dev-viewer --verbose Arc曲線のテッセレーション結果
```

## Notes

- The viewer window blocks the terminal until closed
- Logging defaults: `geolis=info`, `viewer=info`, others=`warn`
- Override with `RUST_LOG` env var for fine-grained control
- See `.claude/rules/viewer.md` for full reference
