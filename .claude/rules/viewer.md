# Viewer (Debug Tool)

## Overview

`examples/viewer/` is a debug viewer for visualising Geolis meshes.
Every `/dev-viewer` invocation **rewrites** `test_meshes.rs` with the requested content, then runs the viewer.

## File Structure

```
examples/viewer/
├── main.rs          — Entry point (DO NOT EDIT)
├── viewer.rs        — Pure viewer UI (DO NOT EDIT)
└── test_meshes.rs   — Rewritten every /dev-viewer invocation
```

| File | Role | Editable? |
|------|------|-----------|
| `main.rs` | Logging init, app startup | No |
| `viewer.rs` | `AppState`, viewport layout, components | No |
| `test_meshes.rs` | Mesh generation — `register_test_meshes()` | **Yes** (rewrite each time) |

## How It Works

1. `/dev-viewer` rewrites `test_meshes.rs` with code that generates the requested meshes
2. Meshes are registered into `MeshStorage` via `upsert_2d()` / `upsert_3d()`
3. The viewer displays them in 2D (left) and 3D (right) viewports

## test_meshes.rs Contract

The function signature must always be:

```rust
pub fn register_test_meshes(storage: &MeshStorage) { ... }
```

Inside, register meshes using:

```rust
// 2D mesh
storage.upsert_2d(RawMesh2DId::new(), Arc::new(mesh_2d));

// 3D mesh
storage.upsert_3d(RawMesh3DId::new(), Arc::new(mesh_3d));
```

### Mesh Sources

| Source | How to Use |
|--------|------------|
| Geolis tessellation APIs | Call `tessellate_stroke()` etc., convert result to `RawMesh2D`/`RawMesh3D` |
| Revion helpers | `Circle`, `Sphere`, `Rectangle` from `revion_core::tessellation` |
| Raw data | Construct `RawMesh2D`/`RawMesh3D` directly with vertices + indices |

## Running

```bash
# Default (geolis logs only)
cargo run --example viewer

# With Revion renderer debug logs
RUST_LOG=revion_renderer=debug cargo run --example viewer

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
