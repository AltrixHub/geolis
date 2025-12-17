# Geolis

An open CAD kernel for everyone.

## Overview

Geolis is an open-source CAD kernel written in Rust that combines accurate geometric representation with robust topology management. The name comes from "Geo" (geometry) + "Polis" (public space) — our mission is to open up the traditionally closed world of CAD kernels.

Starting with features essential for architectural modeling, we aim to progressively achieve commercial-grade quality.

## Features

- **NURBS-based**: Mathematically accurate representation of curves and surfaces
- **BRep structure**: Rigorous topology management enabling area/volume calculations and edge extraction
- **Incremental API**: Use only the features you need
- **Reliability-focused**: Limited features that work correctly

## Architecture

```
┌─────────────────────────────────────────┐
│              Application                │
└─────────────────────────────────────────┘
                    │
                    ▼
┌─────────────────────────────────────────┐
│           Operations Layer              │
│    Extrude / Revolve / Trim / Compute   │
└─────────────────────────────────────────┘
                    │
        ┌───────────┴───────────┐
        ▼                       ▼
┌───────────────────┐   ┌───────────────────┐
│  Topology Layer   │   │  Geometry Layer   │
│                   │   │                   │
│ Solid             │   │ Surface (NURBS)   │
│  └─ Shell         │──▶│ Curve (NURBS)     │
│      └─ Face      │   │ Point             │
│          └─ Loop  │   │                   │
│              └─Edge   │                   │
│                  └─Vertex                 │
└───────────────────┘   └───────────────────┘
                    │
                    ▼
┌─────────────────────────────────────────┐
│           Math Foundation               │
│       nalgebra (vectors & matrices)     │
└─────────────────────────────────────────┘
```

## Feature List

### Phase 1: Foundation (Architectural Modeling Level)

| Category | Feature | Description |
|----------|---------|-------------|
| **Curves** | NURBS curves | Basic curve representation |
| | Bézier curves | Special case of NURBS |
| | Arcs & circles | Accurate representation via NURBS |
| | Point projection | Find closest point on curve |
| | Offset | Used for wall thickness |
| | Curve intersection | Intersection points between curves |
| **Surfaces** | NURBS surfaces | Basic surface representation |
| | Extruded surfaces | Generate walls by extruding curves |
| | Revolved surfaces | Generate columns by revolving curves |
| | Trimmed surfaces | Create holes in surfaces (windows, openings) |
| | Point projection | Find closest point on surface |
| **Topology** | BRep structure | Vertex/Edge/Loop/Face management |
| | Half-edge | Efficient connectivity management |
| | Consistency validation | Ensure topological correctness |
| **Computation** | Area calculation | Accurate surface area computation |
| | Volume calculation | Accurate solid volume computation |
| **Output** | Tessellation | Mesh generation for display |
| | Trimmed surface meshing | Meshing surfaces with holes |

### Phase 2: Feature Expansion

| Category | Feature | Description |
|----------|---------|-------------|
| **Surfaces** | Surface offset | Move surface along normal direction |
| | Loft surfaces | Connect multiple curves with a surface |
| | Sweep surfaces | Move curve along a path |
| **Intersection** | Curve-surface intersection | Intersection points between curve and surface |
| | Plane-surface intersection | Cross-section generation |
| | Planar boolean | Union/difference/intersection in 2D |
| **Computation** | Centroid calculation | Compute center of mass |
| | Interference check | Detect solid intersections |
| **Output** | Adaptive meshing | Subdivision based on curvature |
| **Processing** | Shelling | Add thickness |

### Out of Scope (For Now)

The following features are highly complex and are currently out of scope:

- General 3D boolean operations
- Surface-surface intersection (SSI)
- Fillets and chamfers
- STEP/IGES import/export

## Design Principles

### 1. Separation of Geometry and Topology

We clearly separate the mathematical definition of shapes (geometry) from their connectivity (topology). This allows adding new surface types without affecting the topology layer.

### 2. Incremental API

```rust
// Use curves only
let curve = NurbsCurve::arc(center, radius, start, end);
let point = curve.evaluate(0.5);

// Use surfaces
let wall = ExtrudedSurface::new(curve, height);
let mesh = wall.tessellate();

// Use BRep
let solid = Solid::from_faces(faces);
let volume = solid.volume();
```

Designed so you can use only what you need, reducing the learning curve.

### 3. Focus on Reliability

- Implementation with numerical stability in mind
- API design that prevents creating invalid models
- Continuous topology consistency validation

## Tech Stack

- **Language**: Rust
- **Math foundation**: nalgebra
- **CAD core**: Custom implementation

## Roadmap

```
Phase 1: Foundation (Architectural Modeling Level)
├─ NURBS curve/surface implementation
├─ BRep structure implementation
├─ Extrude/revolve/trim operations
├─ Area/volume calculations
└─ Tessellation

Phase 2: Feature Expansion
├─ Surface offset
├─ Loft/sweep
├─ Planar boolean
└─ Interference check

Phase 3: Stabilization & Optimization
├─ Edge case handling
├─ Performance optimization
└─ Documentation
```

## License

MIT License

## Contributing

As an open-source project, contributions are welcome.

- Issue reports
- Pull requests
- Documentation improvements
- Test case additions
