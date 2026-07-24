//! 2D polygon boolean-subtract via the shared arrangement engine.
//!
//! [`subtract_all_with_holes`] computes `base ∩ (¬⋃subtracts)` — the
//! parts of `base` that are not covered by any of `subtracts[i]`.
//! Output is typed face topology (zero, one, or many
//! [`PolygonWithHoles`]).
//!
//! The engine (split → snap → bilateral classify → face-walk →
//! containment-matrix assemble) is identical to the one used by
//! [`crate::operations::boolean_2d::union_all_with_holes`]; only the
//! fill oracle differs (see [`super::engine::SubtractOracle`]).

use crate::error::Result;

use super::diagnose::assess;
use super::engine::{run_arrangement, SubtractOracle};
use super::types::{signed_area, PolygonWithHoles};
use crate::diagnostics::{InputSnapshot, OpDiagnostic, OpHealth};

/// Subtract a list of regions from a base region. Returns the
/// remaining filled regions as typed face topology.
///
/// Semantics: `result = base ∩ (¬⋃subtracts)` — the parts of `base`
/// not covered by any `subtracts[i]`.
///
/// Special cases:
/// - `subtracts.is_empty()` returns `vec![base]` unchanged (no
///   subtraction performed).
/// - `subtracts` fully covering `base` returns an empty `Vec`.
/// - `subtracts` outside `base` returns `vec![base]` (subtracts
///   outside the base do not affect the result).
/// - A subtract that exactly matches an existing hole of `base` leaves
///   `base` unchanged (subtracting empty space is a no-op).
/// - Overlapping subtracts are de-overlapped by the arrangement engine
///   — the union of the subtract regions is what gets removed.
///
/// # Errors
///
/// Propagates [`crate::error::OperationError::Failed`] from the
/// arrangement engine on the same degenerate-input cases as
/// [`super::union_all_with_holes`] (ambiguous bilateral classification
/// after ε exhaustion, broken parent topology, orientation/depth
/// parity violation).
pub fn subtract_all_with_holes(
    base: PolygonWithHoles,
    subtracts: &[PolygonWithHoles],
) -> Result<Vec<PolygonWithHoles>> {
    if subtracts.is_empty() {
        return Ok(vec![base]);
    }

    // Feed every ring of base + subtracts into the arrangement so their
    // boundaries split each other where they cross.
    let mut segment_inputs: Vec<PolygonWithHoles> = Vec::with_capacity(1 + subtracts.len());
    segment_inputs.push(base.clone());
    segment_inputs.extend(subtracts.iter().cloned());

    let oracle = SubtractOracle {
        base: &base,
        subtracts,
    };
    run_arrangement(&segment_inputs, &oracle)
}

/// [`subtract_all_with_holes`] with a health verdict attached for diagnostics.
///
/// The result is returned verbatim (`Ok`/`Err` unchanged). The verdict is
/// non-`Ok` when the op errored (`Failed`), when a real cut consumed the whole
/// base (`Degenerate` / `EmptyResult`), or when a sliver face survived
/// (`Suspicious`); in those cases a readable [`InputSnapshot`] is attached for
/// the app to log or dump. Cheap input facts are gathered before `base` moves,
/// and the snapshot is built only when the verdict is non-`Ok` — so the clean
/// path costs nothing beyond the assessment.
#[must_use]
pub fn subtract_all_with_holes_diagnosed(
    base: PolygonWithHoles,
    subtracts: &[PolygonWithHoles],
) -> OpDiagnostic<Result<Vec<PolygonWithHoles>>> {
    let base_outer_verts = base.outer.len();
    let base_holes = base.holes.len();
    let base_area = signed_area(&base.outer).abs();
    let cutters = subtracts.len();
    // A real cut against a real base should leave geometry behind; a fully
    // empty result then means the cut swallowed the whole base.
    let expect_nonempty = !subtracts.is_empty() && base_area > super::diagnose::MIN_FACE_AREA;

    let result = subtract_all_with_holes(base, subtracts);

    let health = match &result {
        Err(e) => OpHealth::Failed(e.to_string()),
        Ok(faces) => assess(faces, expect_nonempty),
    };

    if health.is_ok() {
        OpDiagnostic::ok(result)
    } else {
        let snapshot = InputSnapshot::new("boolean_2d::subtract")
            .with("base_outer_verts", base_outer_verts)
            .with("base_holes", base_holes)
            .with("base_area", base_area)
            .with("cutters", cutters);
        OpDiagnostic::flagged(result, health, snapshot)
    }
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::cast_precision_loss,
    reason = "test code: panics are the failure signal; arc facet counts are tiny"
)]
mod tests {
    use std::panic::{catch_unwind, AssertUnwindSafe};

    use super::super::intersect_all_with_holes;
    use super::super::types::{signed_area, Polygon};
    use super::*;

    fn rect(x: f64, y: f64, w: f64, h: f64) -> Polygon {
        vec![(x, y), (x + w, y), (x + w, y + h), (x, y + h)]
    }

    /// CW (hole-winding) rect for use as a hole in a `PolygonWithHoles`.
    fn cw_rect(x: f64, y: f64, w: f64, h: f64) -> Polygon {
        vec![(x, y), (x, y + h), (x + w, y + h), (x + w, y)]
    }

    fn pwh_no_holes(outer: Polygon) -> PolygonWithHoles {
        PolygonWithHoles {
            outer,
            holes: Vec::new(),
        }
    }

    #[test]
    fn diagnosed_clean_subtract_is_ok_without_snapshot() {
        let base = pwh_no_holes(rect(0.0, 0.0, 4.0, 4.0));
        let hole = pwh_no_holes(rect(1.0, 1.0, 1.0, 1.0));
        let d = subtract_all_with_holes_diagnosed(base, &[hole]);
        assert!(d.is_ok());
        assert!(d.inputs.is_none());
        assert!(d.value.is_ok());
    }

    #[test]
    fn diagnosed_full_consumption_is_flagged_with_snapshot() {
        // A cut covering the whole base leaves nothing -> non-Ok verdict.
        let base = pwh_no_holes(rect(0.0, 0.0, 2.0, 2.0));
        let cover = pwh_no_holes(rect(-1.0, -1.0, 4.0, 4.0));
        let d = subtract_all_with_holes_diagnosed(base, &[cover]);
        assert!(!d.is_ok(), "full consumption should be flagged");
        let snap = d.inputs.expect("snapshot attached on non-Ok");
        assert_eq!(snap.op, "boolean_2d::subtract");
        assert!(snap.summary.iter().any(|(k, _)| k == "base_area"));
    }

    #[test]
    fn subtract_empty_list_returns_base_unchanged() {
        let base = pwh_no_holes(rect(0.0, 0.0, 10.0, 10.0));
        let result = subtract_all_with_holes(base.clone(), &[]).expect("subtract must succeed");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], base);
    }

    #[test]
    fn subtract_disjoint_small_rect_punches_hole() {
        // Large outer minus a small disjoint inner rect → outer with one hole.
        let base = pwh_no_holes(rect(0.0, 0.0, 10.0, 10.0));
        let cut = pwh_no_holes(rect(3.0, 3.0, 2.0, 2.0));
        let result = subtract_all_with_holes(base, &[cut]).expect("subtract must succeed");
        assert_eq!(result.len(), 1, "expected one face");
        assert_eq!(result[0].holes.len(), 1, "expected one hole");
        // Hole area = -4 (CW), outer area = 100.
        assert!((signed_area(&result[0].outer) - 100.0).abs() < 0.1);
        assert!((signed_area(&result[0].holes[0]) + 4.0).abs() < 0.1);
    }

    #[test]
    fn subtract_rect_fully_containing_base_returns_empty() {
        let base = pwh_no_holes(rect(2.0, 2.0, 4.0, 4.0));
        let cut = pwh_no_holes(rect(0.0, 0.0, 10.0, 10.0));
        let result = subtract_all_with_holes(base, &[cut]).expect("subtract must succeed");
        assert!(result.is_empty(), "everything removed → empty result");
    }

    #[test]
    fn subtract_rect_covers_all_but_corner_returns_l_shape() {
        // Base 10x10 minus a 9x10 column on the right → 1x10 column on the left.
        // The remaining region is a simple rectangle (no holes).
        let base = pwh_no_holes(rect(0.0, 0.0, 10.0, 10.0));
        let cut = pwh_no_holes(rect(1.0, 0.0, 9.0, 10.0));
        let result = subtract_all_with_holes(base, &[cut]).expect("subtract must succeed");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].holes.len(), 0);
        let area = signed_area(&result[0].outer);
        assert!((area - 10.0).abs() < 0.1, "area={area}, expected 10");
    }

    #[test]
    fn subtract_corner_creates_l_shape() {
        // Base 10x10 minus a 9x9 corner block → L-shaped region (no holes,
        // single outer with 6 vertices).
        let base = pwh_no_holes(rect(0.0, 0.0, 10.0, 10.0));
        let cut = pwh_no_holes(rect(1.0, 1.0, 9.0, 9.0));
        let result = subtract_all_with_holes(base, &[cut]).expect("subtract must succeed");
        assert_eq!(result.len(), 1, "L-shape is a single face");
        assert_eq!(result[0].holes.len(), 0, "L-shape has no holes");
        let area = signed_area(&result[0].outer);
        // Area = 100 - 81 = 19.
        assert!((area - 19.0).abs() < 0.1, "area={area}, expected 19");
    }

    #[test]
    fn subtract_multiple_overlapping_subtracts_handled_correctly() {
        // Base 10x10 minus two overlapping 3x3 rects whose union forms an
        // L-shaped cut. The arrangement engine de-overlaps them
        // automatically.
        let base = pwh_no_holes(rect(0.0, 0.0, 10.0, 10.0));
        let cut_a = pwh_no_holes(rect(3.0, 3.0, 3.0, 3.0));
        let cut_b = pwh_no_holes(rect(5.0, 5.0, 3.0, 3.0));
        let result = subtract_all_with_holes(base, &[cut_a, cut_b]).expect("subtract must succeed");
        assert_eq!(result.len(), 1, "single connected outer remains");
        assert_eq!(result[0].holes.len(), 1, "single merged hole");
        // Outer area = 100; hole area = -(3*3 + 3*3 - 1*1) = -17.
        assert!((signed_area(&result[0].outer) - 100.0).abs() < 0.1);
        let hole_area = signed_area(&result[0].holes[0]);
        assert!(
            (hole_area + 17.0).abs() < 0.1,
            "hole_area={hole_area}, expected -17"
        );
    }

    #[test]
    fn subtract_matching_existing_hole_is_noop() {
        // Base = donut (10x10 with a 4x4 hole at (3,3)). Subtracting a rect
        // exactly matching that hole is a no-op — the subtract sits entirely
        // in empty space.
        let base = PolygonWithHoles {
            outer: rect(0.0, 0.0, 10.0, 10.0),
            holes: vec![cw_rect(3.0, 3.0, 4.0, 4.0)],
        };
        let cut = pwh_no_holes(rect(3.0, 3.0, 4.0, 4.0));
        let result = subtract_all_with_holes(base.clone(), &[cut]).expect("subtract must succeed");
        assert_eq!(result.len(), 1, "donut unchanged");
        assert_eq!(result[0].holes.len(), 1, "still one hole");
        // The hole boundary should be unchanged in area.
        let original_hole_area = signed_area(&base.holes[0]);
        let result_hole_area = signed_area(&result[0].holes[0]);
        assert!(
            (original_hole_area - result_hole_area).abs() < 0.1,
            "hole area changed: orig={original_hole_area}, got={result_hole_area}"
        );
    }

    #[test]
    fn subtract_grows_existing_hole() {
        // Base = donut (10x10 with a 2x2 hole at (4,4)). Subtract a 4x4 rect
        // at (3,3) that fully contains the existing hole and grows it.
        let base = PolygonWithHoles {
            outer: rect(0.0, 0.0, 10.0, 10.0),
            holes: vec![cw_rect(4.0, 4.0, 2.0, 2.0)],
        };
        let cut = pwh_no_holes(rect(3.0, 3.0, 4.0, 4.0));
        let result = subtract_all_with_holes(base, &[cut]).expect("subtract must succeed");
        assert_eq!(result.len(), 1, "single face");
        assert_eq!(result[0].holes.len(), 1, "single merged hole");
        // Merged hole area = -16 (the 4x4 subtract dominates the original 2x2).
        let hole_area = signed_area(&result[0].holes[0]);
        assert!(
            (hole_area + 16.0).abs() < 0.1,
            "hole_area={hole_area}, expected -16"
        );
    }

    #[test]
    fn subtract_outside_base_is_noop() {
        // Subtract a rect that does not overlap the base at all.
        let base = pwh_no_holes(rect(0.0, 0.0, 5.0, 5.0));
        let cut = pwh_no_holes(rect(20.0, 20.0, 3.0, 3.0));
        let result = subtract_all_with_holes(base.clone(), &[cut]).expect("subtract must succeed");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].holes.len(), 0);
        let area = signed_area(&result[0].outer);
        assert!((area - 25.0).abs() < 0.1, "area unchanged");
    }

    #[test]
    fn subtract_splits_base_into_two_pieces() {
        // Subtract a vertical strip that splits the base into two
        // disjoint outer faces. Base 10x4, strip at x∈[4,6] full height.
        let base = pwh_no_holes(rect(0.0, 0.0, 10.0, 4.0));
        let cut = pwh_no_holes(rect(4.0, 0.0, 2.0, 4.0));
        let result = subtract_all_with_holes(base, &[cut]).expect("subtract must succeed");
        assert_eq!(result.len(), 2, "expected two disjoint faces");
        for f in &result {
            assert_eq!(f.holes.len(), 0);
            let area = signed_area(&f.outer);
            assert!(
                (area - 16.0).abs() < 0.1,
                "each piece area=4*4=16, got {area}"
            );
        }
    }

    // ===== Kernel-robustness regression tests (grazing / degenerate input) ==

    /// Net area of a face set: outer areas minus hole areas.
    fn total_area(faces: &[PolygonWithHoles]) -> f64 {
        faces
            .iter()
            .map(|f| {
                signed_area(&f.outer).abs()
                    - f.holes.iter().map(|h| signed_area(h).abs()).sum::<f64>()
            })
            .sum()
    }

    /// Annular band (a curved wall) between inner radius `ri` and outer radius
    /// `ro`, swept over `[a0, a1]`, tessellated into `n` facets per arc and
    /// returned as a CCW polygon.
    fn arc_band(ri: f64, ro: f64, a0: f64, a1: f64, n: usize) -> Polygon {
        let mut pts = Vec::new();
        for i in 0..=n {
            let t = a0 + (a1 - a0) * (i as f64) / (n as f64);
            pts.push((ro * t.cos(), ro * t.sin()));
        }
        for i in 0..=n {
            let t = a1 + (a0 - a1) * (i as f64) / (n as f64);
            pts.push((ri * t.cos(), ri * t.sin()));
        }
        if signed_area(&pts) < 0.0 {
            pts.reverse();
        }
        pts
    }

    /// Every face a boolean op emits must satisfy the winding contract
    /// (CCW outer, CW holes). A wrong-area/topology corruption would break it.
    fn assert_valid_windings(faces: &[PolygonWithHoles]) {
        for f in faces {
            assert!(
                signed_area(&f.outer) > 0.0,
                "outer must be CCW, area={}",
                signed_area(&f.outer)
            );
            for h in &f.holes {
                assert!(
                    signed_area(h) < 0.0,
                    "hole must be CW, area={}",
                    signed_area(h)
                );
            }
        }
    }

    /// A near-degenerate "sliver" subtract must return a typed `Err`, never
    /// panic. Covers the converted CDT post-condition (defect: a spade
    /// `TooSmall` vertex rejection used to `panic!` at the post-condition)
    /// **and** the sub-tolerance bilateral-classification `Err` path. The
    /// `catch_unwind` proves no unwind escapes on reachable input.
    #[test]
    fn cdt_sliver_returns_err_not_panic() {
        // Silence the default panic hook so a (regressed) panic does not spam
        // the test log; `catch_unwind` still reports it as an Err below.
        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));

        let cases: Vec<(&str, PolygonWithHoles, PolygonWithHoles)> = vec![
            // (a) Collinear near-zero-width cut: a full-height sliver only
            //     `WALL_EPS` wide. Its two walls are closer than bilateral
            //     classification can resolve -> typed Err.
            (
                "near-zero-width sliver cut",
                pwh_no_holes(rect(0.0, 0.0, 10.0, 10.0)),
                pwh_no_holes(vec![
                    (5.0, -1.0),
                    (5.0 + 1e-6, -1.0),
                    (5.0 + 1e-6, 11.0),
                    (5.0, 11.0),
                ]),
            ),
            // (b) An input coordinate below spade's MIN_ALLOWED_VALUE
            //     (~1.79e-43) flows to an output vertex; the CDT post-condition
            //     rejects it as `TooSmall`. This is the exact production panic,
            //     now a typed Err.
            (
                "sub-representable output coordinate",
                PolygonWithHoles {
                    outer: vec![(0.0, 1e-50), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
                    holes: Vec::new(),
                },
                pwh_no_holes(rect(3.0, 3.0, 2.0, 2.0)),
            ),
        ];

        let mut verdicts = Vec::new();
        for (label, base, cut) in &cases {
            let outcome = catch_unwind(AssertUnwindSafe(|| {
                subtract_all_with_holes(base.clone(), std::slice::from_ref(cut))
            }));
            verdicts.push((label, outcome));
        }
        std::panic::set_hook(prev_hook);

        for (label, outcome) in verdicts {
            let inner = outcome
                .unwrap_or_else(|_| panic!("'{label}' unwound: degenerate input must not panic"));
            assert!(
                inner.is_err(),
                "'{label}' must return a typed Err on degenerate input, got Ok"
            );
        }
    }

    /// A wide opening rectangle subtracted from a curved wall band, where one
    /// rectangle edge grazes an outer-arc facet near-tangentially. Sweeping the
    /// grazing offset across `1e-7..1e-3`, every result must be either `Ok` with
    /// the analytically-consistent area or a typed `Err` — never a panic, never
    /// a wrong area. The reproduced (moderate-offset) configuration must
    /// succeed. Regression for the "ambiguous half-edge / ε exhausted" failure
    /// and the world-space split-tolerance topology corruption.
    #[test]
    fn grazing_arc_chord_subtract_succeeds() {
        // Curved wall band, mean radius ~5, width 0.6, swept ±0.5 rad.
        let band = pwh_no_holes(arc_band(4.7, 5.3, -0.5, 0.5, 32));
        let band_area = signed_area(&band.outer).abs();

        // The rect's bottom edge is placed a grazing offset `d` below the
        // topmost outer-arc facet vertex (angle a1 = 0.5), so the edge is
        // near-tangent to that facet.
        let apex_y = 5.3 * 0.5_f64.sin();
        let cut = |d: f64| -> PolygonWithHoles {
            let y_bot = apex_y - d;
            pwh_no_holes(vec![(3.0, y_bot), (6.0, y_bot), (6.0, 4.0), (3.0, 4.0)])
        };

        // Analytic anchor via the partition invariant: for the same base and
        // cut, `subtract` and `intersect` split the base, so
        //   area(base ∩ ¬cut) + area(base ∩ cut) == area(base).
        // Evaluated at a robustly-resolvable offset.
        let d_ref = 1e-3;
        let sub_ref = subtract_all_with_holes(band.clone(), &[cut(d_ref)])
            .expect("reference subtract must succeed");
        assert_valid_windings(&sub_ref);
        let area_ref = total_area(&sub_ref);
        let int_ref = intersect_all_with_holes(&band, std::slice::from_ref(&cut(d_ref)))
            .expect("reference intersect must succeed");
        assert!(
            (area_ref + total_area(&int_ref) - band_area).abs() < 1e-6,
            "partition invariant violated: sub={area_ref}, int={}, band={band_area}",
            total_area(&int_ref)
        );

        // The reproduced grazing configuration (moderate offset) must succeed
        // with the analytic area.
        let d_primary = 1e-5;
        let sub_primary = subtract_all_with_holes(band.clone(), &[cut(d_primary)])
            .expect("grazing subtract at 1e-5 must succeed (not Err)");
        assert_valid_windings(&sub_primary);
        assert!(
            (total_area(&sub_primary) - area_ref).abs() < 1e-3,
            "grazing area {} deviates from analytic {area_ref}",
            total_area(&sub_primary)
        );

        // Full ε sweep: never panic, never wrong area; Ok-with-correct-area or
        // typed Err only.
        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let offsets = [1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3];
        let mut results = Vec::new();
        for &d in &offsets {
            let outcome = catch_unwind(AssertUnwindSafe(|| {
                subtract_all_with_holes(band.clone(), &[cut(d)])
            }));
            results.push((d, outcome));
        }
        std::panic::set_hook(prev_hook);

        let mut ok_count = 0;
        for (d, outcome) in results {
            let inner =
                outcome.unwrap_or_else(|_| panic!("grazing subtract panicked at offset {d:e}"));
            // A typed Err on a sub-tolerance offset is allowed; only Ok results
            // are area-checked (they must never be wrong).
            if let Ok(faces) = inner {
                assert_valid_windings(&faces);
                assert!(
                    (total_area(&faces) - area_ref).abs() < 1e-3,
                    "offset {d:e}: area {} is wrong (analytic {area_ref})",
                    total_area(&faces)
                );
                ok_count += 1;
            }
        }
        assert!(
            ok_count >= 5,
            "expected most grazing offsets to resolve, only {ok_count}/9 succeeded"
        );
    }

    /// Oblique float-coordinate split regression, pinning the world-space
    /// `arrangement_split` tolerance (K-2).
    ///
    /// A `~8e-6`-wide full-height slot cut through a long (length-20) wall
    /// places the slot's two walls a fixed WORLD distance apart. On the
    /// length-20 top/bottom host edges that world gap maps to a parameter gap
    /// of `~4e-7` — *below* the old fixed `WALL_EPS` parameter split tolerance
    /// but *above* the length-scaled `WALL_EPS / seg_len` world tolerance. With
    /// a parameter-space tolerance the two crossings collapse into one, a split
    /// point is dropped, a vertex dangles and the face topology corrupts; the
    /// world-space tolerance keeps the crossings distinct.
    ///
    /// The whole configuration is rotated `~17°`, so every coordinate is an
    /// irrational-scale float and no edge is axis-aligned. This is load-bearing:
    /// the integer / axis-aligned fixtures are invariant to the parameter↔world
    /// switch by construction (their crossings never land in the sensitive
    /// tolerance band) and cannot catch this regression. Rotation is
    /// area-preserving, so the expected area is analytic and exact — both the
    /// subtract area and the subtract+intersect partition invariant are pinned
    /// to `< 1e-6`, mirroring `grazing_arc_chord_subtract_succeeds`.
    #[test]
    fn oblique_float_split_partition_is_exact() {
        // Long wall (length 20, height 4). The long top/bottom edges are the
        // hosts whose length makes the parameter↔world tolerance distinction
        // bite.
        const WALL_LEN: f64 = 20.0;
        const WALL_H: f64 = 4.0;
        // Slot width a few WALL_EPS wide: robustly above the bilateral / snap
        // thresholds (~3·WALL_EPS) so the op resolves, yet its two walls map to
        // a parameter gap (~SLOT_W / WALL_LEN ≈ 4e-7) that sits inside the old
        // fixed-WALL_EPS merge band and outside the world-scaled band.
        const SLOT_W: f64 = 8e-6;
        const SLOT_X: f64 = 8.0; // interior: both slot walls cross the long edges

        // ~17°, well off any axis; sin/cos are both irrational-scale floats.
        let theta = 17.0_f64.to_radians();
        let (s, c) = theta.sin_cos();
        let rot = |poly: &Polygon| -> Polygon {
            poly.iter()
                .map(|&(x, y)| (x * c - y * s, x * s + y * c))
                .collect()
        };
        let rot_pwh = |p: &PolygonWithHoles| -> PolygonWithHoles {
            PolygonWithHoles {
                outer: rot(&p.outer),
                holes: p.holes.iter().map(&rot).collect(),
            }
        };

        let wall = pwh_no_holes(rect(0.0, 0.0, WALL_LEN, WALL_H));
        // Full-height through-cut: extends past both long edges so it splits the
        // wall cleanly into two pieces (no partial-height corner cases).
        let slot = pwh_no_holes(vec![
            (SLOT_X, -1.0),
            (SLOT_X + SLOT_W, -1.0),
            (SLOT_X + SLOT_W, WALL_H + 1.0),
            (SLOT_X, WALL_H + 1.0),
        ]);

        let wall_r = rot_pwh(&wall);
        let slot_r = rot_pwh(&slot);
        let wall_area = signed_area(&wall_r.outer).abs();

        // Subtract: the wall splits into two rectangles either side of the slot.
        let sub = subtract_all_with_holes(wall_r.clone(), std::slice::from_ref(&slot_r))
            .expect("oblique split subtract must succeed");
        assert_valid_windings(&sub);
        assert_eq!(
            sub.len(),
            2,
            "full-height slot splits the wall into two faces"
        );

        // Analytic (rotation-preserved) area: the removed slot is SLOT_W × WALL_H.
        let expected_sub = WALL_LEN * WALL_H - SLOT_W * WALL_H;
        let sub_area = total_area(&sub);
        assert!(
            (sub_area - expected_sub).abs() < 1e-6,
            "oblique subtract area {sub_area} deviates from analytic {expected_sub}"
        );

        // Partition invariant: subtract + intersect exactly repartition the wall.
        let int = intersect_all_with_holes(&wall_r, std::slice::from_ref(&slot_r))
            .expect("oblique split intersect must succeed");
        assert_valid_windings(&int);
        let int_area = total_area(&int);
        assert!(
            (sub_area + int_area - wall_area).abs() < 1e-6,
            "partition invariant violated: sub={sub_area}, int={int_area}, wall={wall_area}"
        );
    }
}
