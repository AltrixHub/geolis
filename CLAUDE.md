# Geolis - CAD Kernel for Architectural Modeling

Open-source CAD kernel written in Rust for architectural modeling and general geometric computation.

## Project Structure

```
src/
├── lib.rs              # Library entry point
├── error.rs            # GeolisError (thiserror)
├── math/               # Foundation (nalgebra)
├── geometry/           # Curves (Line, Arc) / Surfaces (Plane)
├── topology/           # BRep: Vertex, Edge, Wire, Face, Shell, Solid
├── operations/         # Creation, Shaping, Boolean, Offset, Transform, Query
└── tessellation/       # Mesh generation (curves, faces, solids)
```

## Architecture

Layered architecture with one-way dependency (top to bottom).

```
Math (nalgebra)
  -> Geometry (curves, surfaces) + Topology (BRep)
    -> Operations (extrude, boolean, offset, etc.)
      -> Tessellation (mesh generation)
```

## Build & Verify

```bash
cargo check --workspace --all-targets        # Quick check
cargo build --workspace --all-targets        # Full build
cargo clippy --workspace --all-targets       # Lint
cargo test --workspace                       # Tests
```

`--all-targets` is required. Without it, examples/tests/benches are not checked.

## Key Rules

- **No `unwrap()` / `expect()`** in production code (tests only)
- **`?` operator** to propagate errors; define custom error types with `thiserror`
- **Minimize `unsafe`** - `// SAFETY:` comment required when used
- **Clippy**: warn on `clippy::all`, `clippy::pedantic`; deny `unwrap_used` / `expect_used`
- Prefer borrowing; avoid unnecessary clone/allocation
- Use iterators; avoid manual indexing

## Commit Convention

```
<type>(<scope>): <summary>
```

- **type**: `feat`, `fix`, `refactor`, `perf`, `docs`, `test`, `chore`
- **scope**: `geometry`, `topology`, `operations`, `tessellation`, `math`
- Summary under 50 chars, imperative mood, no trailing period
- Do NOT add `Co-Authored-By` or `Generated with Claude Code` footer

## Custom Commands

| Command | Description |
|---------|-------------|
| `/dev-implement` | Implement features (includes quality checks) |
| `/dev-viewer` | Rewrite test meshes and run viewer for visual debugging |
| `/impl-verify` | Verify staged changes against implementation plan |
| `/git-commit` | Generate commit message -> `.ai/outputs/git.md` |
| `/git-pr` | Generate PR description -> `.ai/outputs/git.md` |
| `/git-issue` | Create GitHub issue |
| `/doc-config` | Update `.claude/` configuration files |
| `/doc-update` | Update documentation |

## Agents

| Agent | Role |
|-------|------|
| `quality-checker` | Verify build + test + clippy |
| `refactor` | Code refactoring |

## AI Output

Generated artifacts go to `.ai/outputs/`:

- `.ai/outputs/git.md` - Commit/PR messages
- `.ai/outputs/plans/` - Implementation plans
- `.ai/outputs/summaries/` - Task summaries

## Current Phase

**Phase 1: Windows on Planar Walls** - Line/Arc curves, Plane surfaces, BRep topology, basic operations.
