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
//! 3. **Mesh.** Each slab contributes a vertical wall per boundary
//!    segment of its region, plus horizontal caps wherever coverage
//!    changes: a bottom cap over `region(i) \ region(i - 1)` and a top cap
//!    over `region(i) \ region(i + 1)`, both region differences taken with
//!    the same 2D boolean machinery. The slab below the first and above
//!    the last is the empty region, so the body is closed at the extremes.
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

mod mesh;
mod slab;
#[cfg(test)]
mod tests;

use crate::error::{OperationError, Result};
use crate::geometry::pline::Pline;
use crate::operations::boolean_2d::{PolygonWithHoles, WALL_EPS};
use crate::tessellation::TriangleMesh;

/// Chord tolerance used to flatten bulge-encoded arcs when the caller
/// does not set one. Matches
/// [`crate::tessellation::TessellationParams::default`]'s `tolerance`,
/// the crate-wide default deviation budget for curved boundaries.
pub const DEFAULT_ARC_TOLERANCE: f64 = 0.01;

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
/// CW holes, every hole inside its outer — i.e. exactly the outline a
/// caller draws as edges at this elevation.
#[derive(Debug, Clone)]
pub struct PrismSlab {
    z_base: f64,
    z_top: f64,
    faces: Vec<PolygonWithHoles>,
}

impl PrismSlab {
    pub(super) fn new(z_base: f64, z_top: f64, faces: Vec<PolygonWithHoles>) -> Self {
        Self {
            z_base,
            z_top,
            faces,
        }
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

/// Result of [`UnionPrisms::execute`]: the fused body plus the outline it
/// was built from.
#[derive(Debug, Clone)]
pub struct FusedPrisms {
    mesh: TriangleMesh,
    slabs: Vec<PrismSlab>,
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

    /// Consumes the result into its two halves.
    #[must_use]
    pub fn into_parts(self) -> (TriangleMesh, Vec<PrismSlab>) {
        (self.mesh, self.slabs)
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
}

impl UnionPrisms {
    /// Creates the operation with [`DEFAULT_ARC_TOLERANCE`].
    #[must_use]
    pub fn new(profiles: Vec<PrismProfile>) -> Self {
        Self {
            profiles,
            arc_tolerance: DEFAULT_ARC_TOLERANCE,
        }
    }

    /// Sets the maximum chord deviation used to flatten bulge-encoded
    /// arcs. Validated by [`UnionPrisms::execute`].
    #[must_use]
    pub fn with_arc_tolerance(mut self, tolerance: f64) -> Self {
        self.arc_tolerance = tolerance;
        self
    }

    /// Runs the z-slab decomposition and builds the fused mesh.
    ///
    /// # Errors
    ///
    /// - [`OperationError::InvalidInput`] — the arc tolerance is not
    ///   finite and strictly positive.
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

        let intervals = slab::slab_regions(&self.profiles, self.arc_tolerance)?;
        let mesh = mesh::build_mesh(&intervals)?;
        let slabs = intervals
            .into_iter()
            .filter(|slab| !slab.faces().is_empty())
            .collect();
        Ok(FusedPrisms { mesh, slabs })
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
