# Dev Viewer

Visualise meshes in the debug viewer. Every invocation rewrites `test_meshes.rs` to display the requested content.

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

## Instructions

### 1. Understand the Request

Parse the description to determine what meshes to generate. Common patterns:

| Pattern | Approach |
|---------|----------|
| Geolis API output | Call Geolis tessellation/geometry APIs, convert result to `RawMesh2D`/`RawMesh3D` |
| Raw shapes | Use `revion_core::tessellation::{Circle, Sphere, Rectangle}` helpers |
| Specific data | Build `RawMesh2D`/`RawMesh3D` directly with vertices and indices |

### 2. Rewrite test_meshes.rs

1. Read `examples/viewer/test_meshes.rs`
2. **Rewrite** `register_test_meshes()` to generate the requested meshes
3. Add necessary `use` imports at the top
4. Keep the function signature: `pub fn register_test_meshes(storage: &MeshStorage)`
5. Register meshes via `storage.upsert_2d()` / `storage.upsert_3d()`

**Rules:**
- Do NOT modify `main.rs` or `viewer.rs`
- Each invocation fully replaces the previous test mesh code
- If Geolis types need conversion to Revion mesh types, do the conversion in this file

### 3. Build Check

Verify the example compiles:

```bash
cargo check --example viewer
```

If it fails, fix the errors in `test_meshes.rs` before proceeding.

### 4. Run the Viewer

```bash
# Default — geolis/viewer logs only
cargo run --example viewer

# Verbose — include Revion renderer logs
RUST_LOG=revion_renderer=debug cargo run --example viewer
```

Use `--verbose` flag to decide which command to run.

### 5. Report

After launching, display:

```
Viewer launched: {brief description of what is displayed}
- 2D viewport (left): Space+drag pan, Cmd+scroll zoom
- 3D viewport (right): Right-drag orbit, Middle-drag pan, Scroll zoom
```

## File Structure

```
examples/viewer/
├── main.rs          — Entry point (DO NOT EDIT)
├── viewer.rs        — Pure viewer UI (DO NOT EDIT)
└── test_meshes.rs   — Rewrite this every invocation
```

## Examples

```bash
# Visualise tessellation output from Geolis
/dev-viewer tessellate_stroke で L字ポリラインをテッセレートした結果

# Visualise raw geometry
/dev-viewer 半径1.0の円と半径0.5の円を重ねて表示

# Visualise a grid pattern
/dev-viewer 5x5グリッドに球体を並べて表示

# With verbose logs
/dev-viewer --verbose Arc曲線のテッセレーション結果
```

## Notes

- The viewer window blocks the terminal until closed
- Logging defaults: `geolis=info`, `viewer=info`, others=`warn`
- Override with `RUST_LOG` env var for fine-grained control
- See `.claude/rules/viewer.md` for full reference
