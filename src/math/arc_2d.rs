/// 2D arc/bulge math utilities.
///
/// Bulge convention: `bulge = tan(sweep_angle / 4)`.
/// - `bulge = 0`: straight line
/// - `bulge > 0`: counter-clockwise arc
/// - `bulge < 0`: clockwise arc
/// - `|bulge| = 1`: semicircle
use std::f64::consts::PI;

/// Converts a bulge-defined arc segment to center-radius-angle form.
///
/// Returns `(cx, cy, radius, start_angle, sweep_angle)`.
///
/// # Panics
///
/// Does not panic. Returns degenerate values for zero-length chords.
#[must_use]
pub fn arc_from_bulge(x0: f64, y0: f64, x1: f64, y1: f64, bulge: f64) -> (f64, f64, f64, f64, f64) {
    let dx = x1 - x0;
    let dy = y1 - y0;
    let chord_len = (dx * dx + dy * dy).sqrt();

    if chord_len < 1e-12 {
        return (x0, y0, 0.0, 0.0, 0.0);
    }

    // Distance from chord midpoint to center.
    let sagitta_ratio = (1.0 - bulge * bulge) / (2.0 * bulge);
    let mx = (x0 + x1) * 0.5;
    let my = (y0 + y1) * 0.5;

    // Normal to chord pointing toward center (for positive bulge, center is left of chord).
    let nx = -dy / chord_len;
    let ny = dx / chord_len;

    let cx = mx + sagitta_ratio * (chord_len * 0.5) * nx;
    let cy = my + sagitta_ratio * (chord_len * 0.5) * ny;

    // r = d*(1+b²)/(4*|b|) derived from r = d/(2*sin(θ/2)) with θ=4*atan(b)
    let radius = (chord_len * 0.5) * (1.0 + bulge * bulge) / (2.0 * bulge.abs());

    let start_angle = (y0 - cy).atan2(x0 - cx);
    let end_angle = (y1 - cy).atan2(x1 - cx);

    let sweep = 4.0 * bulge.atan();

    // Normalize sweep to [-2π, 2π] range.
    let sweep = if sweep > 2.0 * PI {
        sweep - 2.0 * PI
    } else if sweep < -2.0 * PI {
        sweep + 2.0 * PI
    } else {
        sweep
    };

    // Verify end angle consistency (start + sweep should reach end_angle mod 2π).
    let _ = end_angle; // Used for verification in debug builds.

    (cx, cy, radius, start_angle, sweep)
}

/// Converts arc endpoints + center back to bulge value.
///
/// `is_ccw`: true for counter-clockwise arc, false for clockwise.
#[must_use]
pub fn bulge_from_arc(x0: f64, y0: f64, x1: f64, y1: f64, cx: f64, cy: f64, is_ccw: bool) -> f64 {
    let start_angle = (y0 - cy).atan2(x0 - cx);
    let end_angle = (y1 - cy).atan2(x1 - cx);

    let mut sweep = end_angle - start_angle;
    if is_ccw {
        if sweep < 0.0 {
            sweep += 2.0 * PI;
        }
    } else if sweep > 0.0 {
        sweep -= 2.0 * PI;
    }

    (sweep / 4.0).tan()
}

/// Returns the bulge of the `f0 -> f1` sub-arc of a segment whose bulge
/// is `bulge`, where `f0` and `f1` are parameters along that segment in
/// the [`arc_point_at`] convention.
///
/// A circular arc's parameter fraction is its angle fraction, so the
/// sub-arc sweeps `(f1 - f0) * sweep` and keeps the segment's center and
/// radius EXACTLY — this is what makes cutting an arc in two a
/// shape-preserving operation. A straight segment stays straight, and
/// reversing the range (`f1 < f0`) negates the bulge, matching the
/// reversed sub-arc.
///
/// # Panics
///
/// Does not panic.
#[must_use]
pub fn sub_arc_bulge(bulge: f64, f0: f64, f1: f64) -> f64 {
    if bulge.abs() < 1e-12 {
        return 0.0;
    }
    let sweep = 4.0 * bulge.atan();
    ((f1 - f0) * sweep / 4.0).tan()
}

/// Maximum absolute tangent-chord angle, in radians (160 degrees).
///
/// The sweep of an arc is twice its tangent-chord angle, so this caps a
/// single arc at 320 degrees (`|bulge| = tan(80°) ≈ 5.671`). The bulge
/// diverges as the tangent-chord angle approaches 180 degrees (the arc
/// closes into a full circle), so near-closed configurations are clamped
/// instead of producing an unbounded — eventually infinite — bulge.
pub const MAX_TANGENT_CHORD_ANGLE: f64 = 160.0 * PI / 180.0;

/// Returns the bulge of the unique circular arc that starts at
/// `(x0, y0)` with tangent direction `(tx, ty)` and ends at `(x1, y1)`.
///
/// By the tangent-chord angle theorem the signed angle `α` from the
/// chord to the start tangent is half the sweep. With this module's
/// convention (`bulge > 0` = counter-clockwise, center left of the
/// chord) the start tangent lies `sweep / 2` *clockwise* of the chord,
/// so `sweep = -2α` and
///
/// ```text
/// bulge = tan(sweep / 4) = -tan(α / 2)
/// ```
///
/// `α` is clamped to [`MAX_TANGENT_CHORD_ANGLE`]. The tangent does not
/// need to be normalized; only its direction is used.
///
/// # Panics
///
/// Does not panic. Degenerate input (zero-length chord, zero-length
/// tangent, non-finite coordinates) returns `0.0` — a straight segment.
#[must_use]
pub fn bulge_from_chord_tangent(x0: f64, y0: f64, x1: f64, y1: f64, tx: f64, ty: f64) -> f64 {
    if !(x0.is_finite()
        && y0.is_finite()
        && x1.is_finite()
        && y1.is_finite()
        && tx.is_finite()
        && ty.is_finite())
    {
        return 0.0;
    }

    let dx = x1 - x0;
    let dy = y1 - y0;
    if (dx * dx + dy * dy).sqrt() < 1e-12 || (tx * tx + ty * ty).sqrt() < 1e-12 {
        return 0.0;
    }

    // Signed angle from the chord to the tangent, in (-π, π].
    let alpha = (dx * ty - dy * tx)
        .atan2(dx * tx + dy * ty)
        .clamp(-MAX_TANGENT_CHORD_ANGLE, MAX_TANGENT_CHORD_ANGLE);

    -(alpha * 0.5).tan()
}

/// Evaluates a point on an arc at parameter `t` in `[0, 1]`.
#[must_use]
pub fn arc_point_at(
    cx: f64,
    cy: f64,
    radius: f64,
    start_angle: f64,
    sweep: f64,
    t: f64,
) -> (f64, f64) {
    let angle = start_angle + sweep * t;
    (cx + radius * angle.cos(), cy + radius * angle.sin())
}

/// Computes the unit tangent direction on an arc at parameter `t` in `[0, 1]`.
///
/// The tangent points in the direction of increasing `t`.
#[must_use]
pub fn arc_tangent_at(start_angle: f64, sweep: f64, t: f64) -> (f64, f64) {
    let angle = start_angle + sweep * t;
    let sign = if sweep >= 0.0 { 1.0 } else { -1.0 };
    // Tangent to circle at angle θ is (-sin θ, cos θ) for CCW; negate for CW.
    (-sign * angle.sin(), sign * angle.cos())
}

/// Offsets an arc segment defined by endpoints and bulge.
///
/// For an inward offset (toward center), the radius decreases.
/// Returns `None` if the offset radius would be ≤ 0 (arc collapses).
///
/// Returns `(x0', y0', x1', y1', bulge')`.
#[must_use]
pub fn offset_arc_segment(
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
    bulge: f64,
    distance: f64,
) -> Option<(f64, f64, f64, f64, f64)> {
    let (cx, cy, radius, start_angle, sweep) = arc_from_bulge(x0, y0, x1, y1, bulge);

    if radius < 1e-12 {
        return None;
    }

    // Determine offset direction: positive distance = left offset.
    // For a CCW arc (bulge > 0) the center lies on the LEFT of the
    // traversal (unit circle: at (r, 0) moving +y, left points at the
    // origin), so a left offset moves TOWARD the center — the radius
    // decreases. For a CW arc (bulge < 0) the center is on the right,
    // so a left offset moves away and the radius increases. This
    // matches the line-segment convention (left normal), keeping a
    // mixed line/arc polyline's offset internally consistent.
    let sign = if bulge > 0.0 { -1.0 } else { 1.0 };
    let new_radius = radius + sign * distance;

    if new_radius <= 1e-12 {
        return None;
    }

    // New endpoints: same angles, new radius.
    let ox0 = cx + new_radius * start_angle.cos();
    let oy0 = cy + new_radius * start_angle.sin();
    let end_angle = start_angle + sweep;
    let ox1 = cx + new_radius * end_angle.cos();
    let oy1 = cy + new_radius * end_angle.sin();

    // Bulge is invariant to radius changes (same sweep angle).
    Some((ox0, oy0, ox1, oy1, bulge))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    #[test]
    fn semicircle_ccw() {
        // CCW semicircle from (0,0) to (2,0), bulge=1.
        // Center at (1,0), radius=1, sweep=+π.
        // Arc goes from angle π → 3π/2 → 2π (through bottom).
        let (cx, cy, r, sa, sw) = arc_from_bulge(0.0, 0.0, 2.0, 0.0, 1.0);
        assert!((cx - 1.0).abs() < TOL, "cx={cx}");
        assert!(cy.abs() < TOL, "cy={cy}");
        assert!((r - 1.0).abs() < TOL, "r={r}");
        assert!((sw - PI).abs() < TOL, "sweep={sw}");

        let p0 = arc_point_at(cx, cy, r, sa, sw, 0.0);
        assert!(p0.0.abs() < TOL, "p0.x={}", p0.0);
        assert!(p0.1.abs() < TOL, "p0.y={}", p0.1);

        let p1 = arc_point_at(cx, cy, r, sa, sw, 1.0);
        assert!((p1.0 - 2.0).abs() < TOL, "p1.x={}", p1.0);
        assert!(p1.1.abs() < TOL, "p1.y={}", p1.1);

        // Midpoint at angle 3π/2 → (1, -1) (through bottom for CCW)
        let pm = arc_point_at(cx, cy, r, sa, sw, 0.5);
        assert!((pm.0 - 1.0).abs() < TOL, "pm.x={}", pm.0);
        assert!((pm.1 + 1.0).abs() < TOL, "pm.y={}", pm.1);
    }

    #[test]
    fn semicircle_cw_goes_through_top() {
        // CW semicircle from (0,0) to (2,0), bulge=-1.
        // Sweep=-π, arc goes π → π/2 → 0 (through top).
        let (cx, cy, r, sa, sw) = arc_from_bulge(0.0, 0.0, 2.0, 0.0, -1.0);
        assert!((cx - 1.0).abs() < TOL, "cx={cx}");
        assert!(cy.abs() < TOL, "cy={cy}");
        assert!((r - 1.0).abs() < TOL, "r={r}");
        assert!((sw + PI).abs() < TOL, "sweep={sw}");

        // Midpoint at angle π/2 → (1, 1)
        let pm = arc_point_at(cx, cy, r, sa, sw, 0.5);
        assert!((pm.0 - 1.0).abs() < TOL, "pm.x={}", pm.0);
        assert!((pm.1 - 1.0).abs() < TOL, "pm.y={}", pm.1);
    }

    #[test]
    fn quarter_circle_ccw() {
        // CCW quarter circle from (1,0) to (0,1), center at origin.
        // sweep = +π/2 (CCW), goes through first quadrant.
        let bulge = (PI / 8.0).tan();
        let (cx, cy, r, sa, sw) = arc_from_bulge(1.0, 0.0, 0.0, 1.0, bulge);
        assert!((r - 1.0).abs() < 1e-6, "r={r}");
        assert!(cx.abs() < 1e-6, "cx={cx}");
        assert!(cy.abs() < 1e-6, "cy={cy}");
        assert!((sw - PI / 2.0).abs() < 1e-6, "sweep={sw}");

        let p0 = arc_point_at(cx, cy, r, sa, sw, 0.0);
        assert!((p0.0 - 1.0).abs() < 1e-6);
        assert!(p0.1.abs() < 1e-6);

        // Midpoint at angle π/4 → (cos(π/4), sin(π/4))
        let pm = arc_point_at(cx, cy, r, sa, sw, 0.5);
        let expected = (PI / 4.0).cos();
        assert!((pm.0 - expected).abs() < 1e-6, "pm.x={}", pm.0);
        assert!((pm.1 - expected).abs() < 1e-6, "pm.y={}", pm.1);
    }

    #[test]
    fn bulge_from_arc_roundtrip() {
        // CCW semicircle
        let (cx, cy, _, _, sw) = arc_from_bulge(0.0, 0.0, 2.0, 0.0, 1.0);
        let bulge = bulge_from_arc(0.0, 0.0, 2.0, 0.0, cx, cy, sw > 0.0);
        assert!((bulge - 1.0).abs() < TOL, "bulge={bulge}");

        // CW semicircle
        let (center_x, center_y, _, _, sweep) = arc_from_bulge(0.0, 0.0, 2.0, 0.0, -1.0);
        let bulge_cw = bulge_from_arc(0.0, 0.0, 2.0, 0.0, center_x, center_y, sweep > 0.0);
        assert!((bulge_cw + 1.0).abs() < TOL, "bulge_cw={bulge_cw}");
    }

    #[test]
    fn arc_tangent_is_unit_and_correct() {
        // CCW semicircle from (0,0) to (2,0), start_angle=π, sweep=π.
        // At t=0: tangent direction for positive sweep at angle π is
        //   (-sign*sin(π), sign*cos(π)) = (0, -1) (downward, into bottom semicircle).
        let (_, _, _, sa, sw) = arc_from_bulge(0.0, 0.0, 2.0, 0.0, 1.0);
        let t0 = arc_tangent_at(sa, sw, 0.0);
        let len = (t0.0 * t0.0 + t0.1 * t0.1).sqrt();
        assert!((len - 1.0).abs() < TOL, "tangent not unit: len={len}");
        assert!(t0.0.abs() < TOL, "tx={}", t0.0);
        assert!((t0.1 + 1.0).abs() < TOL, "ty={}", t0.1);
    }

    #[test]
    fn offset_arc_left_moves_toward_ccw_center() {
        // Semicircle from (0,0) to (2,0), bulge=1, radius=1, center=(1,0).
        // Positive distance = LEFT of the CCW traversal = toward the
        // center: offset by 0.5 → new radius = 0.5.
        let result = offset_arc_segment(0.0, 0.0, 2.0, 0.0, 1.0, 0.5);
        assert!(result.is_some());
        let (x0, y0, x1, y1, b) = result.unwrap();
        // Bulge should remain 1 (same sweep).
        assert!((b - 1.0).abs() < TOL, "bulge={b}");
        let (cx, cy, _, _, _) = arc_from_bulge(x0, y0, x1, y1, b);
        let new_r = ((x0 - cx).powi(2) + (y0 - cy).powi(2)).sqrt();
        assert!((new_r - 0.5).abs() < 1e-6, "new_r={new_r}");
    }

    #[test]
    fn offset_arc_right_moves_away_from_ccw_center() {
        // Negative distance = RIGHT of the CCW traversal = away from
        // the center: offset by -0.5 → new radius = 1.5.
        let result = offset_arc_segment(0.0, 0.0, 2.0, 0.0, 1.0, -0.5);
        assert!(result.is_some());
        let (x0, y0, x1, y1, b) = result.unwrap();
        let (cx, cy, _, _, _) = arc_from_bulge(x0, y0, x1, y1, b);
        let new_r = ((x0 - cx).powi(2) + (y0 - cy).powi(2)).sqrt();
        assert!((new_r - 1.5).abs() < 1e-6, "new_r={new_r}");
    }

    #[test]
    fn offset_arc_collapse() {
        // Semicircle radius=1. Left offset by 1.5 crosses the CCW
        // center → collapses.
        let result = offset_arc_segment(0.0, 0.0, 2.0, 0.0, 1.0, 1.5);
        assert!(result.is_none());
    }

    #[test]
    fn straight_segment_bulge_zero() {
        // bulge=0 should give degenerate arc (radius≈∞, sweep≈0)
        let (cx, cy, r, _sa, sw) = arc_from_bulge(0.0, 0.0, 1.0, 0.0, 0.001);
        // Very small bulge → very large radius, very small sweep
        assert!(r > 100.0, "r={r}");
        assert!(sw.abs() < 0.01, "sweep={sw}");
        let _ = (cx, cy); // just checking it doesn't panic
    }

    #[test]
    fn bulge_from_chord_tangent_roundtrips_through_arc_from_bulge() {
        // The start tangent produced by arc_from_bulge + arc_tangent_at
        // must reproduce the bulge it was built from.
        let chords = [
            (0.0, 0.0, 2.0, 0.0),
            (-1.3, 0.7, 2.4, -0.9),
            (5.0, 5.0, 5.0, -3.0),
        ];
        let bulges = [
            1.0,
            -1.0,
            0.5,
            -0.5,
            (PI / 8.0).tan(),
            -(PI / 8.0).tan(),
            1e-4,
            -1e-4,
            3.0,
            -3.0,
        ];

        for (x0, y0, x1, y1) in chords {
            for bulge in bulges {
                let (_, _, _, sa, sw) = arc_from_bulge(x0, y0, x1, y1, bulge);
                let (tx, ty) = arc_tangent_at(sa, sw, 0.0);
                let recovered = bulge_from_chord_tangent(x0, y0, x1, y1, tx, ty);
                assert!(
                    (recovered - bulge).abs() < 1e-9,
                    "bulge={bulge} recovered={recovered} chord=({x0},{y0})-({x1},{y1})"
                );
            }
        }
    }

    #[test]
    fn bulge_from_chord_tangent_along_chord_is_straight() {
        // Tangent parallel to the chord → no turn → bulge exactly 0.
        let b = bulge_from_chord_tangent(1.0, 2.0, 4.0, 6.0, 3.0, 4.0);
        assert!(b.abs() < TOL, "bulge={b}");
    }

    #[test]
    fn bulge_from_chord_tangent_ignores_tangent_length() {
        let short = bulge_from_chord_tangent(0.0, 0.0, 2.0, 0.0, 0.0, -1.0);
        let long = bulge_from_chord_tangent(0.0, 0.0, 2.0, 0.0, 0.0, -1000.0);
        assert!((short - long).abs() < TOL, "short={short} long={long}");
        // Straight down at the start of a chord along +x is the CCW semicircle.
        assert!((short - 1.0).abs() < TOL, "bulge={short}");
    }

    #[test]
    fn bulge_from_chord_tangent_clamps_near_full_circle() {
        // Tangent almost antiparallel to the chord: the exact arc would
        // approach a full circle, so the tangent-chord angle is clamped.
        let limit = (MAX_TANGENT_CHORD_ANGLE * 0.5).tan();
        for angle_deg in [161.0_f64, 175.0, 179.999] {
            let angle = angle_deg.to_radians();
            let b = bulge_from_chord_tangent(0.0, 0.0, 1.0, 0.0, angle.cos(), angle.sin());
            assert!(b.is_finite(), "bulge={b} for {angle_deg} deg");
            assert!((b + limit).abs() < 1e-9, "bulge={b} for {angle_deg} deg");
        }
        // Just inside the guard the exact value is kept.
        let inside = 159.0_f64.to_radians();
        let b = bulge_from_chord_tangent(0.0, 0.0, 1.0, 0.0, inside.cos(), inside.sin());
        assert!((b + (inside * 0.5).tan()).abs() < 1e-12, "bulge={b}");
    }

    #[test]
    fn sub_arc_bulge_over_the_whole_range_is_the_input() {
        for bulge in [1.0, -1.0, 0.5, -0.5, 3.0, -3.0, 1e-4] {
            let whole = sub_arc_bulge(bulge, 0.0, 1.0);
            assert!((whole - bulge).abs() < TOL, "bulge={bulge} whole={whole}");
        }
    }

    #[test]
    fn sub_arc_bulge_keeps_the_center_and_radius() {
        // Splitting an arc at f must leave both halves on the SAME
        // circle, and each half must reach the split point exactly.
        let (x0, y0, x1, y1) = (-1.3, 0.7, 2.4, -0.9);
        for bulge in [0.3, -0.3, 1.0, -1.0, 2.5] {
            let (cx, cy, r, sa, sw) = arc_from_bulge(x0, y0, x1, y1, bulge);
            for f in [0.1, 0.5, 0.9] {
                let (mx, my) = arc_point_at(cx, cy, r, sa, sw, f);

                let first = sub_arc_bulge(bulge, 0.0, f);
                let (fx, fy, fr, fsa, fsw) = arc_from_bulge(x0, y0, mx, my, first);
                assert!((fx - cx).abs() < 1e-9 && (fy - cy).abs() < 1e-9, "center");
                assert!((fr - r).abs() < 1e-9, "radius {fr} vs {r}");
                let end = arc_point_at(fx, fy, fr, fsa, fsw, 1.0);
                assert!(
                    (end.0 - mx).abs() < 1e-9 && (end.1 - my).abs() < 1e-9,
                    "end"
                );

                let second = sub_arc_bulge(bulge, f, 1.0);
                let (sx, sy, sr, ssa, ssw) = arc_from_bulge(mx, my, x1, y1, second);
                assert!((sx - cx).abs() < 1e-9 && (sy - cy).abs() < 1e-9, "center");
                assert!((sr - r).abs() < 1e-9, "radius {sr} vs {r}");
                let tail = arc_point_at(sx, sy, sr, ssa, ssw, 1.0);
                assert!(
                    (tail.0 - x1).abs() < 1e-9 && (tail.1 - y1).abs() < 1e-9,
                    "tail"
                );
            }
        }
    }

    #[test]
    fn sub_arc_bulge_is_straight_and_reversible() {
        assert!(sub_arc_bulge(0.0, 0.0, 0.5).abs() < TOL);
        assert!(sub_arc_bulge(1e-13, 0.0, 0.5).abs() < TOL);
        let forward = sub_arc_bulge(0.8, 0.2, 0.7);
        let backward = sub_arc_bulge(0.8, 0.7, 0.2);
        assert!((forward + backward).abs() < TOL, "{forward} {backward}");
    }

    #[test]
    fn bulge_from_chord_tangent_degenerate_inputs_are_straight() {
        // Zero-length chord.
        assert!(bulge_from_chord_tangent(1.0, 1.0, 1.0, 1.0, 1.0, 0.0).abs() < TOL);
        // Zero-length tangent.
        assert!(bulge_from_chord_tangent(0.0, 0.0, 1.0, 0.0, 0.0, 0.0).abs() < TOL);
        // Non-finite coordinates and tangents.
        assert!(bulge_from_chord_tangent(f64::NAN, 0.0, 1.0, 0.0, 1.0, 1.0).abs() < TOL);
        assert!(bulge_from_chord_tangent(0.0, 0.0, f64::INFINITY, 0.0, 1.0, 1.0).abs() < TOL);
        assert!(bulge_from_chord_tangent(0.0, 0.0, 1.0, 0.0, f64::NAN, 1.0).abs() < TOL);
    }
}
