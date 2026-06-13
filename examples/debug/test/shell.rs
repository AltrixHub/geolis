use geolis::math::{Point3, Vector3};
use geolis::operations::creation::{MakeFace, MakeWire};
use geolis::operations::shaping::Extrude;
use geolis::tessellation::{TessellateSolid, TessellationParams};
use geolis::topology::TopologyStore;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label, SceneBounds};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const BLUE: Color = Color::rgb(100, 150, 255);
const RED: Color = Color::rgb(230, 100, 100);
const EDGE: Color = Color::rgb(40, 40, 40);

/// min/max で直方体を作る
fn make_box(
    store: &mut TopologyStore,
    min: Point3,
    max: Point3,
) -> Option<geolis::topology::SolidId> {
    let (x0, y0, z0) = (min.x, min.y, min.z);
    let (x1, y1, z1) = (max.x, max.y, max.z);
    let pts = vec![
        Point3::new(x0, y0, z0),
        Point3::new(x1, y0, z0),
        Point3::new(x1, y1, z0),
        Point3::new(x0, y1, z0),
    ];
    let wire = MakeWire::new(pts, true).execute(store).ok()?;
    let face = MakeFace::new(wire, vec![]).execute(store).ok()?;
    Extrude::new(face, Vector3::new(0.0, 0.0, z1 - z0))
        .execute(store)
        .ok()
}

fn render(
    storage: &MeshStorage,
    bounds: &mut SceneBounds,
    store: &TopologyStore,
    solid: geolis::topology::SolidId,
    color: Color,
) {
    if let Ok(mesh) = TessellateSolid::new(solid, TessellationParams::default()).execute(store) {
        register_face(storage, bounds, mesh, color);
    }
    if let Ok(s) = store.solid(solid) {
        register_edges(storage, bounds, store, s.outer_shell, EDGE);
    }
}

pub fn register(storage: &MeshStorage) {
    let mut bounds = SceneBounds::empty();
    let spacing = 14.0_f64;

    // ─────────────────────────────────────────────────────────────
    // Case 2 正解:
    //   Shell(4×4×4, bx=14, +x 除去, t=0.5) を y=2.5 で断面
    //   外箱 (14,0,0)-(18,4,4) / 内空洞 (14.5,0.5,0.5)-(18,3.5,3.5)
    //
    //   素材領域を 4 つの直方体に分解 (前半 y<2.5):
    //     A: 左壁  (14  , 0  , 0  )-(14.5, 2.5, 4  )
    //     B: 前板  (14.5, 0  , 0  )-(18  , 0.5, 4  )  ← y=0..0.5 は空洞なし
    //     C: 底板  (14.5, 0.5, 0  )-(18  , 2.5, 0.5)
    //     D: 天板  (14.5, 0.5, 3.5)-(18  , 2.5, 4  )
    // ─────────────────────────────────────────────────────────────
    {
        let bx = spacing;
        let by = 0.0_f64;
        register_label(
            storage,
            &mut bounds,
            bx - 2.0,
            by + 8.0,
            "2",
            LABEL_SIZE,
            LABEL_COLOR,
        );

        let mut store = TopologyStore::new();
        let pieces: &[(Point3, Point3)] = &[
            // A: 左壁
            (
                Point3::new(bx, by, 0.0),
                Point3::new(bx + 0.5, by + 2.5, 4.0),
            ),
            // B: 前板 (y=0..0.5 — 内空洞なし)
            (
                Point3::new(bx + 0.5, by, 0.0),
                Point3::new(bx + 4.0, by + 0.5, 4.0),
            ),
            // C: 底板
            (
                Point3::new(bx + 0.5, by + 0.5, 0.0),
                Point3::new(bx + 4.0, by + 2.5, 0.5),
            ),
            // D: 天板
            (
                Point3::new(bx + 0.5, by + 0.5, 3.5),
                Point3::new(bx + 4.0, by + 2.5, 4.0),
            ),
        ];
        for &(min, max) in pieces {
            if let Some(solid) = make_box(&mut store, min, max) {
                render(storage, &mut bounds, &store, solid, BLUE);
            }
        }
    }

    // ─────────────────────────────────────────────────────────────
    // Case 3 正解:
    //   Shell(6×6×6, bx=28, top 除去, t=1.5) を y=3.5 で断面
    //   外箱 (28,0,0)-(34,6,6) / 内空洞 (29.5,1.5,1.5)-(32.5,4.5,6)
    //
    //   素材領域を 4 つの直方体に分解 (前半 y<3.5):
    //     A: 左壁  (28  , 0  , 0  )-(29.5, 3.5, 6  )
    //     B: 右壁  (32.5, 0  , 0  )-(34  , 3.5, 6  )
    //     C: 前板  (29.5, 0  , 0  )-(32.5, 1.5, 6  )  ← y=0..1.5 は空洞なし
    //     D: 底板  (29.5, 1.5, 0  )-(32.5, 3.5, 1.5)
    // ─────────────────────────────────────────────────────────────
    {
        let bx = spacing * 2.0;
        let by = 0.0_f64;
        register_label(
            storage,
            &mut bounds,
            bx - 2.0,
            by + 8.0,
            "3",
            LABEL_SIZE,
            LABEL_COLOR,
        );

        let mut store = TopologyStore::new();
        let pieces: &[(Point3, Point3)] = &[
            // A: 左壁
            (
                Point3::new(bx, by, 0.0),
                Point3::new(bx + 1.5, by + 3.5, 6.0),
            ),
            // B: 右壁
            (
                Point3::new(bx + 4.5, by, 0.0),
                Point3::new(bx + 6.0, by + 3.5, 6.0),
            ),
            // C: 前板 (y=0..1.5 — 内空洞なし)
            (
                Point3::new(bx + 1.5, by, 0.0),
                Point3::new(bx + 4.5, by + 1.5, 6.0),
            ),
            // D: 底板
            (
                Point3::new(bx + 1.5, by + 1.5, 0.0),
                Point3::new(bx + 4.5, by + 3.5, 1.5),
            ),
        ];
        for &(min, max) in pieces {
            if let Some(solid) = make_box(&mut store, min, max) {
                render(storage, &mut bounds, &store, solid, RED);
            }
        }
    }
}
