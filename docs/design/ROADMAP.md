# Geolib - CAD Kernel Specification

Functional specification for "Geolib", a CAD kernel for architectural modeling.
Ultimate goal: Architectural modeling with windows on curved surfaces.

---

## 1. Geometry

### 1.1 Curves

| Type | Description | Parameters | Use Cases | Phase |
|------|-------------|------------|-----------|-------|
| Line | Straight line | Start point, end point | Wall edges, window frames | 1 |
| Arc | Circular arc | Center, radius, start angle, end angle | Arched windows, curved wall edges | 1 |
| Circle | Circle | Center, radius, normal | Circular windows, round column cross-sections | 2 |
| Ellipse | Ellipse | Center, major axis, minor axis | Elliptical windows | 2 |
| NurbsCurve | Free-form curve | Control points, knots, degree | Organic shapes | 3 |

### 1.2 Surfaces

| Type | Description | Parameters | Use Cases | Phase |
|------|-------------|------------|-----------|-------|
| Plane | Flat plane | Origin, normal | General walls, floors, ceilings | 1 |
| Cylinder | Cylindrical surface | Axis, radius | Curved walls, round columns | 2 |
| Cone | Conical surface | Axis, apex, angle | Tapered towers | 2 |
| Sphere | Spherical surface | Center, radius | Domes | 2 |
| Torus | Torus | Center, major radius, minor radius | Arch intersections | 2 |
| RuledSurface | Ruled surface | Two guide curves | Inclined curved walls | 3 |
| NurbsSurface | Free-form surface | Control points, knots, degree | Organic architecture | 3 |

---

## 2. Topology

| Element | Description | Stored Data | Role | Phase |
|---------|-------------|-------------|------|-------|
| Vertex | Vertex | Point3 | Edge endpoints | 1 |
| Edge | Edge | Curve, start Vertex, end Vertex | Constitutes face boundaries | 1 |
| Wire | Wire | List of Edges | Face contour (loop) | 1 |
| Face | Face | Surface, outer Wire, inner Wire[] | Constitutes solid boundaries | 1 |
| Shell | Shell | List of Faces | Closed set of faces | 1 |
| Solid | Solid | Shell | 3D solid body | 1 |

---

## 3. Operations

### 3.1 Creation Operations

| Operation | Input | Output | Description | Phase |
|-----------|-------|--------|-------------|-------|
| MakeBox | Width, height, depth | Solid | Generate a rectangular box | 1 |
| MakeCylinder | Radius, height | Solid | Generate a cylinder | 2 |
| MakeSphere | Radius | Solid | Generate a sphere | 2 |
| MakeCone | Radius, height, angle | Solid | Generate a cone | 2 |
| MakeFace | Wire | Face | Create a face from a wire | 1 |
| MakeWire | Edge[] | Wire | Create a wire from edges | 1 |
| MakeSolid | Shell | Solid | Create a solid from a shell | 1 |

### 3.2 Shaping Operations

| Operation | Input | Output | Description | Phase |
|-----------|-------|--------|-------------|-------|
| Extrude | Face, direction vector | Solid | Extrusion | 1 |
| Revolve | Face, axis, angle | Solid | Revolution | 2 |
| Sweep | Profile Face, Path Curve | Solid | Sweep | 3 |
| Loft | Face[] | Solid | Loft (cross-section interpolation) | 3 |

### 3.3 Boolean Operations

| Operation | Input | Output | Description | Phase |
|-----------|-------|--------|-------------|-------|
| Union | Solid, Solid | Solid | Union (merge) | 1 |
| Subtract | Solid, Solid | Solid | Subtraction (cut) | 1 |
| Intersect | Solid, Solid | Solid | Intersection (common part) | 2 |

### 3.4 Offset Operations

| Operation | Input | Output | Description | Phase |
|-----------|-------|--------|-------------|-------|
| CurveOffset2D | Curve, distance, corner style | Curve[] | 2D curve offset | 1 |
| WireOffset2D | Wire, distance, corner style | Wire[] | 2D wire offset | 1 |
| ThickenFace | Face, distance | Solid | Thicken a face | 1 |
| FaceOffset | Face, distance | Face | Face offset | 2 |
| SolidOffset | Solid, distance | Solid | Offset entire solid | 2 |
| NurbsCurveOffset | NurbsCurve, distance | NurbsCurve | NURBS curve offset | 3 |
| NurbsSurfaceOffset | NurbsSurface, distance | NurbsSurface | NURBS surface offset | 3 |

#### Offset Corner Styles

| Style | Description | Shape | Use Case |
|-------|-------------|-------|----------|
| Miter | Extend to intersection | `┐` | General walls |
| Round | Connect with arc | `╮` | Rounded corners |
| Bevel | Straight-line chamfer | `┒` | Chamfering |
| Square | Extend at right angle | `▐` | End treatment |

### 3.5 Modification Operations

| Operation | Input | Output | Description | Phase |
|-----------|-------|--------|-------------|-------|
| Fillet | Edge/Wire, radius | Solid/Wire | Round edges | 2 |
| Chamfer | Edge/Wire, distance | Solid/Wire | Chamfer edges | 2 |
| Trim | Face, Wire | Face | Trim face (create opening) | 1 |
| Split | Solid, Plane/Surface | Solid[] | Split | 2 |
| Shell | Solid, thickness, removed Face[] | Solid | Hollow out | 2 |

### 3.6 Transformation Operations

| Operation | Input | Output | Description | Phase |
|-----------|-------|--------|-------------|-------|
| Translate | Shape, vector | Shape | Translation | 1 |
| Rotate | Shape, axis, angle | Shape | Rotation | 1 |
| Scale | Shape, factor | Shape | Scaling | 1 |
| Mirror | Shape, plane | Shape | Mirroring | 1 |
| Transform | Shape, Matrix4x4 | Shape | Arbitrary transformation | 1 |

### 3.7 Query Operations

| Operation | Input | Output | Description | Phase |
|-----------|-------|--------|-------------|-------|
| PointOnCurve | Curve, t | Point3, Tangent | Curve evaluation | 1 |
| PointOnSurface | Surface, u, v | Point3, Normal | Surface evaluation | 1 |
| ClosestPointOnCurve | Curve, Point | Point3, t | Closest point on curve | 1 |
| ClosestPointOnSurface | Surface, Point | Point3, (u,v) | Closest point on surface | 2 |
| CurveCurveIntersect | Curve, Curve | Point3[] | Curve-curve intersection | 1 |
| CurveSurfaceIntersect | Curve, Surface | Point3[] | Curve-surface intersection | 2 |
| SurfaceSurfaceIntersect | Surface, Surface | Curve[] | Surface-surface intersection | 3 |
| BoundingBox | Shape | Box3 | Bounding box | 1 |
| Area | Face | f64 | Area | 2 |
| Volume | Solid | f64 | Volume | 2 |
| Length | Curve | f64 | Curve length | 1 |
| IsValid | Shape | bool | Validity check | 1 |

---

## 4. Tessellation (Meshing)

| Function | Input | Output | Description | Phase |
|----------|-------|--------|-------------|-------|
| TessellateFace | Face | TriangleMesh | Triangulate a face | 1 |
| TessellateSolid | Solid | TriangleMesh | Mesh an entire solid | 1 |
| TessellateCurve | Curve, tolerance | Polyline | Approximate curve with line segments | 1 |
| AdaptiveTessellation | Face, tolerance | TriangleMesh | Curvature-adaptive subdivision | 2 |
| TessellateWithHoles | Face (with holes) | TriangleMesh | Mesh a face with holes | 1 |

---

## 5. Phase Summary

### Phase 1: Windows on Planar Walls

| Category | Features |
|----------|----------|
| Curve | Line, Arc |
| Surface | Plane |
| Topology | Vertex, Edge, Wire, Face, Shell, Solid |
| Creation | MakeBox, MakeFace, MakeWire, MakeSolid |
| Shaping | Extrude |
| Boolean | Union, Subtract |
| Offset | CurveOffset2D, WireOffset2D, ThickenFace |
| Modification | Trim |
| Transformation | Translate, Rotate, Scale, Mirror, Transform |
| Query | PointOnCurve, PointOnSurface, ClosestPointOnCurve, CurveCurveIntersect, BoundingBox, Length, IsValid |
| Tessellation | TessellateFace, TessellateSolid, TessellateCurve, TessellateWithHoles |

### Phase 2: Windows on Analytic Surfaces

| Category | Additional Features |
|----------|---------------------|
| Curve | Circle, Ellipse |
| Surface | Cylinder, Cone, Sphere, Torus |
| Creation | MakeCylinder, MakeSphere, MakeCone |
| Shaping | Revolve |
| Boolean | Intersect |
| Offset | FaceOffset, SolidOffset |
| Modification | Fillet, Chamfer, Split, Shell |
| Query | ClosestPointOnSurface, CurveSurfaceIntersect, Area, Volume |
| Tessellation | AdaptiveTessellation |

### Phase 3: Windows on NURBS Surfaces

| Category | Additional Features |
|----------|---------------------|
| Curve | NurbsCurve |
| Surface | NurbsSurface, RuledSurface |
| Shaping | Sweep, Loft |
| Offset | NurbsCurveOffset, NurbsSurfaceOffset |
| Query | SurfaceSurfaceIntersect |

---

## 6. Mapping to Architectural Elements

| Architectural Element | Surface | Curve | Primary Operations | Phase |
|-----------------------|---------|-------|--------------------|-------|
| Straight wall | Plane | Line | Extrude, WireOffset2D | 1 |
| Curved wall | Cylinder | Arc | Extrude, WireOffset2D | 2 |
| Free-form wall | NurbsSurface | NurbsCurve | Extrude, NurbsSurfaceOffset | 3 |
| Floor / Ceiling | Plane | Line/Arc | Extrude | 1 |
| Rectangular window | - | Line | Trim, Subtract | 1 |
| Circular window | - | Circle | Trim, Subtract | 2 |
| Arched window | - | Line + Arc | Trim, Subtract | 1 |
| Window on curved surface | - | NurbsCurve | Trim | 3 |
| Column (rectangular) | Plane | Line | Extrude | 1 |
| Column (round) | Cylinder | Circle | Extrude | 2 |
| Dome | Sphere | - | Revolve | 2 |
| Stairs | Plane | Line | Extrude, Boolean | 1 |
| Ramp | RuledSurface | Line | Loft | 3 |

---

## 7. Precision and Computational Characteristics

| Type | Offset Precision | Intersection | Speed |
|------|-----------------|--------------|-------|
| Line | Exact | Exact | Fastest |
| Arc/Circle | Exact | Exact | Fast |
| Plane | Exact | Exact | Fast |
| Cylinder/Sphere/Cone | Exact | Exact | Fast |
| Torus | Exact | Analytical | Medium |
| NurbsCurve | Approximate | Numerical | Slow |
| NurbsSurface | Approximate | Numerical | Slowest |

---

## 8. Recommended Libraries

| Functionality | Recommended Library | Reason |
|---------------|---------------------|--------|
| Linear algebra | nalgebra | Rust standard, high performance |
| Robust geometric predicates | robust | Floating-point error handling |
| Triangulation | spade or earcutr | Delaunay / Earcut |
| NURBS | **Custom** | Full control, optimized for architectural use |
| 2D Geometry | **Custom** | Core of the CAD kernel |

---

## 9. Wall Generation Flow Example (Phase 1)

```
Input: Wall centerline (2D Wire), wall thickness, wall height

    Wall centerline
    ────────────
         │
         ↓ WireOffset2D(+thickness/2)
    ────────────  Outer wall line
         │
         ↓ WireOffset2D(-thickness/2)
    ────────────  Inner wall line
         │
         ↓ MakeFace(outer wall line - inner wall line)
    ┌──────────┐  Wall cross-section (2D Face)
    │          │
    └──────────┘
         │
         ↓ Extrude(height)
    ╔══════════╗
    ║          ║  Wall (3D Solid)
    ╚══════════╝
         │
         ↓ Subtract(window opening)
    ╔════╔══╗══╗
    ║    ║  ║  ║  Wall with window
    ╚════╚══╝══╝
```
