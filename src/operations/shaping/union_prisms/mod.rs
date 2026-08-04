//! Fused body of many vertical prisms, built by z-slab decomposition
//! instead of a 3D boolean.
//!
//! [`UnionPrisms`] takes a list of [`PrismProfile`]s — a 2D footprint
//! extruded between two z values, plus the rectangular-in-z cuts
//! (openings) that belong to it — and returns one triangle mesh for the
//! whole group together with the fused outline at every z interval.
//!
//! The typical caller is a BIM modeller fusing the walls of one floor
//! that share a kind: their prisms overlap at every junction, and
//! extruding them individually leaves internal faces inside the joint.
//! Fusing in 2D removes those faces before any 3D geometry exists.
//!
//! # Algorithm
//!
//! 1. **Breakpoints.** Every profile's `z_base` / `z_top` and every cut's
//!    `z_base` / `z_top` (clamped into its profile's span) are collected,
//!    sorted, and deduplicated with [`Z_EPS`]. Consecutive breakpoints
//!    delimit the z slabs.
//! 2. **Per-slab region.** Inside one slab the body's cross-section is
//!    constant, so it is a pure 2D problem: union the footprints of every
//!    profile whose span covers the slab
//!    ([`crate::operations::boolean_2d::union_all_with_holes`]), then
//!    subtract every cut active in the slab
//!    ([`crate::operations::boolean_2d::subtract_all_with_holes`]). A cut
//!    is active when its own z span covers the slab **and** its owning
//!    profile is active there — but it is subtracted from the *fused*
//!    region, so an opening carves the neighbouring profiles' material at
//!    a junction too.
//!
//!    The union depends only on the slab's **cover** — the set of
//!    profiles reaching it — so it is computed once per distinct cover
//!    and reused (`slab::FusedByCover`); only the carve is per slab.
//!    That matters because openings multiply slabs without changing who
//!    covers them: one window height turns a floor into three slabs that
//!    all fuse the same walls. Measured on a 31-band floor with ten
//!    windows (release, M2 Air), fusing three slabs from one cover
//!    instead of three: 1.43 ms → 0.79 ms per call.
//! 3. **Caps.** Horizontal faces appear wherever coverage changes: a
//!    down-facing cap over `region(i) \ region(i - 1)` at the slab's base
//!    and an up-facing cap over `region(i) \ region(i + 1)` at its top,
//!    both region differences taken with the same 2D boolean machinery.
//!    The slab below the first and above the last is the empty region, so
//!    the body is closed at the extremes. See [`cap`].
//! 4. **Mesh.** Each slab contributes a vertical wall per boundary
//!    segment of its region; each cap is triangulated at its elevation.
//! 5. **Corner edges.** Each slab's rings are walked once more and every
//!    vertex where the boundary genuinely turns yields a vertical segment
//!    spanning the slab — the arris a caller draws so a wall keeps its
//!    side definition. See [`corner`] for the angle test.
//!
//! # Drawing the result
//!
//! The body's edges are the cap rings ([`FusedPrisms::caps`], horizontal)
//! plus the corner edges ([`FusedPrisms::corner_edges`], vertical) — and
//! nothing else. A slab's own rings are NOT edges of the body: two
//! stacked slabs share the material along their whole common
//! cross-section, so the surface runs straight through the interface and
//! only the part a cap covers is a real horizontal arris. That is why the
//! caps are the drawn horizontal outline: one opening splits the whole
//! group into three slabs, and drawing each slab's ring at both its ends
//! would paint a line across every element of the group at that opening's
//! sill and head, where the material is in fact continuous.
//!
//! Because every z coordinate in the mesh is a breakpoint and every
//! horizontal boundary is a boolean result, the surface is geometrically
//! closed: it bounds exactly the union solid, and its divergence-theorem
//! volume is exact. It is not combinatorially conforming — a cap edge may
//! T-junction against a longer wall edge of the neighbouring slab — so it
//! is a display / volume mesh, not a `BRep`.
//!
//! # Degenerate input
//!
//! Rings that collapse (fewer than three distinct points, or zero area)
//! are **skipped**, as are profiles left with no footprint and slabs with
//! no material. Only contract violations — an unclosed ring, a
//! non-positive z interval, a non-finite coordinate — are errors.
//!
//! # Where this lives
//!
//! `shaping` is the crate's extrusion family (`Extrude`, `MakeLoft`,
//! `MakeHipRoof`): operations that turn planar profiles into a body. This
//! is that, for many profiles at once, and it consumes `boolean_2d`
//! exactly the way `offset::curve_band::carve_band_faces` does — one
//! level up from the arrangement engine, never inside it.

mod cap;
mod corner;
mod mesh;
mod slab;
#[cfg(test)]
mod tests;

use crate::error::{OperationError, Result};
use crate::geometry::pline::Pline;
use crate::math::Point3;
use crate::operations::boolean_2d::{PolygonWithHoles, WALL_EPS};
use crate::tessellation::TriangleMesh;

/// Chord tolerance used to flatten bulge-encoded arcs when the caller
/// does not set one. Matches
/// [`crate::tessellation::TessellationParams::default`]'s `tolerance`,
/// the crate-wide default deviation budget for curved boundaries.
pub const DEFAULT_ARC_TOLERANCE: f64 = 0.01;

/// Default angular threshold above which a boundary vertex counts as a
/// corner, in radians.
///
/// # Derivation
///
/// `Pline::to_points` bounds a chord's turn angle by the sagitta
/// criterion `θ_max = 2·acos(1 − tol / r)`, which for `tol / r ≤ 0.1`
/// equals its series form `2·√(2·tol / r)`. `θ_max` **decreases** as the
/// radius grows, so the worst-case facet sits at the smallest radius the
/// default promises to handle — `r = 10 × arc_tolerance`, a 0.1 fillet at
/// [`DEFAULT_ARC_TOLERANCE`] and the tightest an architectural footprint
/// realistically carries:
///
/// ```text
/// θ_max = 2·√(2 / 10) = 0.894427… rad ≈ 51.2°
/// ```
///
/// Rounded up to `0.9` rad (≈ 51.6°) for float margin. That leaves the
/// threshold above every arc facet a default tessellation can produce and
/// far below the 90° turn of a wall end or an L / T / X junction, so a
/// default-tessellated circle yields **zero** corner edges while every
/// real corner keeps its arris.
///
/// # Outside the derivation
///
/// An arc with `r < 10 × arc_tolerance` is flattened too coarsely for
/// *any* angle test to tell its facets from real corners — an `r = 0.05`
/// circle at the default tolerance turns 60° per chord. Tighten the arc
/// tolerance (which shrinks `θ_max` with `√tol`) rather than raising this
/// threshold, or set both explicitly.
pub const DEFAULT_CORNER_ANGLE_TOLERANCE: f64 = 0.9;

/// Tolerance for every z comparison: breakpoint deduplication, slab
/// height, and z-interval validation.
///
/// Deliberately the 2D pipeline's [`WALL_EPS`]. The z axis carries the
/// same length unit as x and y, so two elevations that the plan-view
/// arrangement would not be able to keep apart must not be able to open a
/// slab of their own either.
pub const Z_EPS: f64 = WALL_EPS;

/// A closed planar region: one outer ring plus zero or more hole rings.
///
/// Rings are [`Pline`]s, so a segment may be a bulge-encoded circular arc
/// — arcs are flattened once, with the operation's arc tolerance, on the
/// way into the boolean pipeline.
///
/// # Contract
///
/// - Every ring is closed and has at least two vertices (a two-vertex
///   ring is legal only when its bulges make it a lens or a circle).
/// - Holes lie inside the outer ring and do not overlap each other.
/// - Ring **winding is irrelevant**: containment is decided by the
///   arrangement engine's winding-number classifier, which is
///   orientation-blind. Outputs are re-wound to the
///   [`PolygonWithHoles`] contract regardless.
///
/// The first bullet is validated by [`PrismRegion::try_from_parts`]; the
/// second is the caller's, exactly as for
/// [`crate::operations::boolean_2d::subtract_all_with_holes`] inputs.
#[derive(Debug, Clone)]
pub struct PrismRegion {
    outer: Pline,
    holes: Vec<Pline>,
}

impl PrismRegion {
    /// Builds a region from an outer ring and its holes.
    ///
    /// # Errors
    ///
    /// [`OperationError::InvalidInput`] when a ring is not closed, has
    /// fewer than two vertices, or carries a non-finite coordinate or
    /// bulge. Rings that are well-formed but degenerate (zero area, all
    /// vertices coincident) are accepted here and skipped at execution —
    /// see the module docs.
    pub fn try_from_parts(outer: Pline, holes: Vec<Pline>) -> Result<Self> {
        validate_ring(&outer, "outer")?;
        for (index, hole) in holes.iter().enumerate() {
            validate_ring(hole, &format!("hole[{index}]"))?;
        }
        Ok(Self { outer, holes })
    }

    /// The outer ring.
    #[must_use]
    pub fn outer(&self) -> &Pline {
        &self.outer
    }

    /// The hole rings.
    #[must_use]
    pub fn holes(&self) -> &[Pline] {
        &self.holes
    }
}

/// A rectangular-in-z opening: `region` is removed from the fused body
/// between `z_base` and `z_top`.
///
/// A cut belongs to the profile that owns it, which decides *where* it is
/// active in z; what it removes is the fused material of the whole group,
/// so a window in one wall keeps cutting through the neighbour it is
/// joined to.
#[derive(Debug, Clone)]
pub struct PrismCut {
    region: PrismRegion,
    z_base: f64,
    z_top: f64,
}

impl PrismCut {
    /// Builds a cut spanning `z_base..z_top`.
    ///
    /// # Errors
    ///
    /// [`OperationError::InvalidInput`] when a bound is non-finite or
    /// `z_top` does not exceed `z_base` by more than [`Z_EPS`].
    pub fn try_new(region: PrismRegion, z_base: f64, z_top: f64) -> Result<Self> {
        validate_z_span("PrismCut::try_new", z_base, z_top)?;
        Ok(Self {
            region,
            z_base,
            z_top,
        })
    }

    /// The plan-view region removed inside the cut's z span.
    #[must_use]
    pub fn region(&self) -> &PrismRegion {
        &self.region
    }

    /// Lower bound of the cut's z span.
    #[must_use]
    pub fn z_base(&self) -> f64 {
        self.z_base
    }

    /// Upper bound of the cut's z span.
    #[must_use]
    pub fn z_top(&self) -> f64 {
        self.z_top
    }
}

/// One prism of the group: a footprint (possibly several disjoint faces)
/// extruded from `z_base` to `z_top`, with the openings that belong to it.
#[derive(Debug, Clone)]
pub struct PrismProfile {
    faces: Vec<PrismRegion>,
    z_base: f64,
    z_top: f64,
    cuts: Vec<PrismCut>,
}

impl PrismProfile {
    /// Builds a profile from its footprint faces and z span.
    ///
    /// `faces` may hold several disjoint regions — that is the shape a
    /// band footprint already has after
    /// [`crate::operations::offset::curve_band::carve_band_faces`].
    ///
    /// # Errors
    ///
    /// [`OperationError::InvalidInput`] when a bound is non-finite or
    /// `z_top` does not exceed `z_base` by more than [`Z_EPS`]. An empty
    /// `faces` list is legal and contributes nothing.
    pub fn try_new(faces: Vec<PrismRegion>, z_base: f64, z_top: f64) -> Result<Self> {
        validate_z_span("PrismProfile::try_new", z_base, z_top)?;
        Ok(Self {
            faces,
            z_base,
            z_top,
            cuts: Vec::new(),
        })
    }

    /// Adds one opening.
    #[must_use]
    pub fn with_cut(mut self, cut: PrismCut) -> Self {
        self.cuts.push(cut);
        self
    }

    /// Adds a batch of openings.
    #[must_use]
    pub fn with_cuts(mut self, cuts: impl IntoIterator<Item = PrismCut>) -> Self {
        self.cuts.extend(cuts);
        self
    }

    /// The footprint faces.
    #[must_use]
    pub fn faces(&self) -> &[PrismRegion] {
        &self.faces
    }

    /// Lower bound of the prism's z span.
    #[must_use]
    pub fn z_base(&self) -> f64 {
        self.z_base
    }

    /// Upper bound of the prism's z span.
    #[must_use]
    pub fn z_top(&self) -> f64 {
        self.z_top
    }

    /// The openings that belong to this profile.
    #[must_use]
    pub fn cuts(&self) -> &[PrismCut] {
        &self.cuts
    }
}

/// The fused cross-section over one z interval.
///
/// `faces` is the interval's boundary as typed face topology — CCW outer,
/// CW holes, every hole inside its outer — i.e. the plan shape the
/// interval's vertical surface is extruded from.
///
/// It is NOT the horizontal outline at this elevation: a boundary shared
/// with the neighbouring slab carries no horizontal edge at all. The
/// drawn horizontal edges are [`FusedPrisms::caps`].
#[derive(Debug, Clone)]
pub struct PrismSlab {
    z_base: f64,
    z_top: f64,
    faces: Vec<PolygonWithHoles>,
    cover: Vec<usize>,
    parts: Vec<PolygonWithHoles>,
    cut_keys: Vec<CutKey>,
    cuts: Vec<PolygonWithHoles>,
}

/// Identity of one active cut: `(profile index, cut index within that
/// profile)`. Compared as ids, never as geometry, so two slabs agree on
/// "the same cut" exactly.
pub(super) type CutKey = (usize, usize);

impl PrismSlab {
    pub(super) fn new(
        z_base: f64,
        z_top: f64,
        faces: Vec<PolygonWithHoles>,
        cover: Vec<usize>,
        parts: Vec<PolygonWithHoles>,
        cut_keys: Vec<CutKey>,
        cuts: Vec<PolygonWithHoles>,
    ) -> Self {
        Self {
            z_base,
            z_top,
            faces,
            cover,
            parts,
            cut_keys,
            cuts,
        }
    }

    /// The profiles whose z span reaches this interval, by index,
    /// ascending. Two slabs with equal covers were carved out of the
    /// SAME fused material — the identity the cap stage exploits.
    pub(super) fn cover(&self) -> &[usize] {
        &self.cover
    }

    /// The cover's material BEFORE fusing — one entry per covering
    /// profile face. `⋃ parts` is the same region [`Self::faces`] was
    /// carved out of, but keeping the parts separate lets a consumer
    /// clip a small region against only the few parts that reach it,
    /// instead of against the whole floor's fused boundary.
    pub(super) fn parts(&self) -> &[PolygonWithHoles] {
        &self.parts
    }

    /// The cuts active in this interval, by identity.
    pub(super) fn cut_keys(&self) -> &[CutKey] {
        &self.cut_keys
    }

    /// The flattened regions of [`Self::cut_keys`], in the same order.
    pub(super) fn cuts(&self) -> &[PolygonWithHoles] {
        &self.cuts
    }

    /// Lower bound of the interval.
    #[must_use]
    pub fn z_base(&self) -> f64 {
        self.z_base
    }

    /// Upper bound of the interval.
    #[must_use]
    pub fn z_top(&self) -> f64 {
        self.z_top
    }

    /// The fused boundary rings at this interval.
    #[must_use]
    pub fn faces(&self) -> &[PolygonWithHoles] {
        &self.faces
    }
}

/// Which way a horizontal cap faces, i.e. which side of the interface
/// carries the material.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CapFacing {
    /// Outward normal `-z`, at a slab's `z_base`: the body's underside
    /// there — the soffit over an opening, or the bottom of the stack.
    Down,
    /// Outward normal `+z`, at a slab's `z_top`: the body's upper
    /// surface there — the sill under an opening, the exposed step of a
    /// lower run, or the top of the stack.
    Up,
}

/// One horizontal face of the fused body: the part of a slab's
/// cross-section its neighbour across the interface does NOT cover, at
/// the elevation where the two meet.
///
/// # Contract
///
/// These faces are the body's only horizontal surfaces, so their rings —
/// outer and holes alike — are the body's only horizontal edges. A
/// caller drawing the body draws every ring of every cap at its own
/// [`Self::z`], and takes the vertical half from
/// [`FusedPrisms::corner_edges`]; see the module docs for why a slab's
/// own rings are not edges.
///
/// The same faces are what the mesh stage triangulates, so the drawn
/// outline and the meshed surface can never describe different bodies.
#[derive(Debug, Clone)]
pub struct CapFace {
    z: f64,
    facing: CapFacing,
    region: PolygonWithHoles,
}

impl CapFace {
    pub(super) fn new(z: f64, facing: CapFacing, region: PolygonWithHoles) -> Self {
        Self { z, facing, region }
    }

    /// The elevation the cap lies at — always a slab bound, hence a z
    /// breakpoint.
    #[must_use]
    pub fn z(&self) -> f64 {
        self.z
    }

    /// Which way the cap's outward normal points.
    #[must_use]
    pub fn facing(&self) -> CapFacing {
        self.facing
    }

    /// The capped area, as typed face topology — CCW outer, CW holes.
    #[must_use]
    pub fn region(&self) -> &PolygonWithHoles {
        &self.region
    }
}

/// A vertical arris of the fused body: the segment a caller draws at one
/// corner of one slab.
///
/// Both endpoints share the corner's `(x, y)`; only `z` differs, spanning
/// the slab the corner belongs to. A corner shared by two stacked slabs
/// produces two segments that meet at the shared elevation — see
/// [`corner`] for why they are not merged.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CornerEdge {
    base: Point3,
    top: Point3,
}

impl CornerEdge {
    /// Lower endpoint, at the slab's `z_base`.
    #[must_use]
    pub fn base(&self) -> Point3 {
        self.base
    }

    /// Upper endpoint, at the slab's `z_top`.
    #[must_use]
    pub fn top(&self) -> Point3 {
        self.top
    }
}

/// Result of [`UnionPrisms::execute`]: the fused body, the cross-sections
/// it was built from, and the horizontal + vertical edges that draw it.
#[derive(Debug, Clone)]
pub struct FusedPrisms {
    mesh: TriangleMesh,
    slabs: Vec<PrismSlab>,
    caps: Vec<CapFace>,
    corner_edges: Vec<CornerEdge>,
}

impl FusedPrisms {
    /// The fused triangle mesh. Empty when no profile carried material.
    #[must_use]
    pub fn mesh(&self) -> &TriangleMesh {
        &self.mesh
    }

    /// The non-empty slabs, ordered bottom to top. Intervals with no
    /// material are omitted.
    #[must_use]
    pub fn slabs(&self) -> &[PrismSlab] {
        &self.slabs
    }

    /// The horizontal caps of the body, bottom to top: a slab's base caps
    /// before its top caps.
    ///
    /// These are the body's horizontal edges — see [`CapFace`] and the
    /// module docs. Empty exactly when the body is.
    #[must_use]
    pub fn caps(&self) -> &[CapFace] {
        &self.caps
    }

    /// The vertical corner edges of every slab, bottom to top.
    ///
    /// Together with the cap rings — the body's horizontal edges — these
    /// are the edges that give the drawn body its shape: the caps say
    /// where it steps, the corner edges say where it turns. Facets of a
    /// flattened arc are excluded; see
    /// [`UnionPrisms::with_corner_angle_tolerance`].
    #[must_use]
    pub fn corner_edges(&self) -> &[CornerEdge] {
        &self.corner_edges
    }

    /// Consumes the result into its four parts.
    #[must_use]
    pub fn into_parts(self) -> (TriangleMesh, Vec<PrismSlab>, Vec<CapFace>, Vec<CornerEdge>) {
        (self.mesh, self.slabs, self.caps, self.corner_edges)
    }
}

/// Fuses a group of vertical prisms into one body without a 3D boolean.
///
/// See the module docs for the z-slab algorithm and the guarantees on the
/// resulting mesh.
#[derive(Debug, Clone)]
pub struct UnionPrisms {
    profiles: Vec<PrismProfile>,
    arc_tolerance: f64,
    corner_angle_tolerance: f64,
}

impl UnionPrisms {
    /// Creates the operation with [`DEFAULT_ARC_TOLERANCE`] and
    /// [`DEFAULT_CORNER_ANGLE_TOLERANCE`].
    #[must_use]
    pub fn new(profiles: Vec<PrismProfile>) -> Self {
        Self {
            profiles,
            arc_tolerance: DEFAULT_ARC_TOLERANCE,
            corner_angle_tolerance: DEFAULT_CORNER_ANGLE_TOLERANCE,
        }
    }

    /// Sets the maximum chord deviation used to flatten bulge-encoded
    /// arcs. Validated by [`UnionPrisms::execute`].
    #[must_use]
    pub fn with_arc_tolerance(mut self, tolerance: f64) -> Self {
        self.arc_tolerance = tolerance;
        self
    }

    /// Sets the turn angle, in radians, above which a boundary vertex
    /// gets a vertical corner edge. Validated by
    /// [`UnionPrisms::execute`].
    ///
    /// The default is derived from [`DEFAULT_ARC_TOLERANCE`] — see
    /// [`DEFAULT_CORNER_ANGLE_TOLERANCE`]. A caller that tessellates arcs
    /// finer than the default can lower this in step to keep shallower
    /// corners: the arc facet bound falls as `2·√(2·tol / r)`.
    #[must_use]
    pub fn with_corner_angle_tolerance(mut self, tolerance: f64) -> Self {
        self.corner_angle_tolerance = tolerance;
        self
    }

    /// Runs the z-slab decomposition and builds the fused mesh, outline
    /// and corner edges.
    ///
    /// # Errors
    ///
    /// - [`OperationError::InvalidInput`] — the arc tolerance is not
    ///   finite and strictly positive, or the corner angle tolerance is
    ///   not finite and inside `0 < t < π` (no turn can exceed π, so a
    ///   threshold at or above it could never fire).
    /// - [`OperationError::Failed`] — propagated from the 2D arrangement
    ///   engine on degenerate input (bilateral classification still
    ///   ambiguous after ε exhaustion, broken parent topology,
    ///   orientation/depth parity violation), exactly as for
    ///   [`crate::operations::boolean_2d::union_all_with_holes`].
    /// - [`crate::error::TessellationError::Failed`] — a cap region could
    ///   not be triangulated.
    pub fn execute(&self) -> Result<FusedPrisms> {
        if !self.arc_tolerance.is_finite() || self.arc_tolerance <= 0.0 {
            return Err(OperationError::InvalidInput(format!(
                "UnionPrisms::execute: arc tolerance must be finite and > 0; got {}",
                self.arc_tolerance
            ))
            .into());
        }
        if !self.corner_angle_tolerance.is_finite()
            || self.corner_angle_tolerance <= 0.0
            || self.corner_angle_tolerance >= std::f64::consts::PI
        {
            return Err(OperationError::InvalidInput(format!(
                "UnionPrisms::execute: corner angle tolerance must be finite and in \
                 (0, π) radians; got {}",
                self.corner_angle_tolerance
            ))
            .into());
        }

        let intervals = slab::slab_regions(&self.profiles, self.arc_tolerance)?;
        // Neighbour coverage is read off `i ± 1`, so the caps are taken
        // from the CONTIGUOUS interval list, before the empty ones are
        // dropped from the public result.
        let caps = cap::cap_faces(&intervals)?;
        let mesh = mesh::build_mesh(&intervals, &caps)?;
        let slabs: Vec<PrismSlab> = intervals
            .into_iter()
            .filter(|slab| !slab.faces().is_empty())
            .collect();
        let corner_edges = corner::corner_edges(&slabs, self.corner_angle_tolerance);
        Ok(FusedPrisms {
            mesh,
            slabs,
            caps,
            corner_edges,
        })
    }
}

/// Validates a ring's intrinsics: closed, at least two vertices, finite
/// coordinates and bulges.
fn validate_ring(ring: &Pline, label: &str) -> Result<()> {
    if !ring.closed {
        return Err(OperationError::InvalidInput(format!(
            "PrismRegion::try_from_parts: {label} must have closed = true"
        ))
        .into());
    }
    if ring.vertices.len() < 2 {
        return Err(OperationError::InvalidInput(format!(
            "PrismRegion::try_from_parts: {label} must have at least 2 vertices; got {}",
            ring.vertices.len()
        ))
        .into());
    }
    for (index, vertex) in ring.vertices.iter().enumerate() {
        if !vertex.x.is_finite() || !vertex.y.is_finite() || !vertex.bulge.is_finite() {
            return Err(OperationError::InvalidInput(format!(
                "PrismRegion::try_from_parts: {label} vertex {index} is not finite \
                 (x = {}, y = {}, bulge = {})",
                vertex.x, vertex.y, vertex.bulge
            ))
            .into());
        }
    }
    Ok(())
}

/// Validates a z interval: finite bounds, strictly positive height.
fn validate_z_span(op: &str, z_base: f64, z_top: f64) -> Result<()> {
    if !z_base.is_finite() || !z_top.is_finite() {
        return Err(OperationError::InvalidInput(format!(
            "{op}: z bounds must be finite; got z_base = {z_base}, z_top = {z_top}"
        ))
        .into());
    }
    if z_top - z_base <= Z_EPS {
        return Err(OperationError::InvalidInput(format!(
            "{op}: z_top must exceed z_base by more than {Z_EPS}; \
             got z_base = {z_base}, z_top = {z_top}"
        ))
        .into());
    }
    Ok(())
}
