use crate::geometry::pline::{Pline, PlineVertex};

use super::slice::PlineSlice;

/// Stitches valid slices back into polylines by matching endpoints.
///
/// Uses a greedy approach: for each slice's end point, find the closest
/// unconnected slice whose start point matches, preferring forward index
/// distance (`cavalier_contours` style).
///
/// When `input_closed` is true, all results are marked closed (original
/// behavior for closed-polyline offsets).  When false, each chain is
/// checked: if the first and last vertices coincide the result is closed,
/// otherwise it stays open.
#[must_use]
pub fn connect(slices: &[&PlineSlice], input_closed: bool) -> Vec<Pline> {
    if slices.is_empty() {
        return Vec::new();
    }

    let n = slices.len();
    let mut used = vec![false; n];
    let mut results = Vec::new();
    let tol_sq = 1e-8;

    for start in 0..n {
        if used[start] {
            continue;
        }

        used[start] = true;
        let mut chain_verts: Vec<PlineVertex> = slices[start].vertices.clone();
        let mut current = start;

        // Try to extend the chain by finding a slice whose start matches our end.
        loop {
            let end_v = chain_verts.last().copied();
            let Some(end_pt) = end_v else { break };

            // Find the best unvisited slice whose start is close to our end.
            let mut best: Option<usize> = None;
            let mut best_dist_sq = tol_sq;

            for candidate in 0..n {
                if used[candidate] {
                    continue;
                }
                let cand_start = &slices[candidate].vertices[0];
                let dx = cand_start.x - end_pt.x;
                let dy = cand_start.y - end_pt.y;
                let dist_sq = dx * dx + dy * dy;
                if dist_sq < best_dist_sq {
                    best_dist_sq = dist_sq;
                    best = Some(candidate);
                }
            }

            if let Some(next) = best {
                used[next] = true;
                // Append the next slice's vertices (skip the first since it overlaps).
                chain_verts.extend_from_slice(&slices[next].vertices[1..]);
                current = next;
            } else {
                break;
            }
        }

        if chain_verts.len() < 2 {
            continue;
        }

        // Check if the chain forms a closed loop.
        let first = &chain_verts[0];
        let last = &chain_verts[chain_verts.len() - 1];
        let dx = last.x - first.x;
        let dy = last.y - first.y;
        let endpoints_coincide = dx * dx + dy * dy < tol_sq;

        if endpoints_coincide {
            // Remove the duplicate closing vertex.
            chain_verts.pop();
        }

        let is_closed = if input_closed {
            true
        } else {
            endpoints_coincide
        };

        if !is_closed || chain_verts.len() >= 3 {
            results.push(Pline {
                vertices: chain_verts,
                closed: is_closed,
            });
        }

        let _ = current; // suppress unused warning
    }

    results
}
