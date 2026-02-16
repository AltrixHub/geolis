//! Pure viewer UI — state definition and component tree.
//!
//! This module is independent of what meshes are loaded.
//! Feed any `MeshStorage` via `AppState` and the viewer renders it.

#![allow(clippy::too_many_lines)]

use revion_macro::rsx;
use revion_ui::value_objects::{Dimension, Edges};
use revion_ui::{
    style::{AlignItems, FlexDirection, JustifyContent, LayoutStyle, VisualStyle},
    MeshStorage, RenderContext, Theme, VNode, ViewerType,
};

/// Minimal application state — holds only mesh storage.
#[derive(Debug, Clone, Default)]
pub struct AppState {
    pub mesh_storage: MeshStorage,
}

/// Root UI component — dual 2D / 3D viewports with a status bar.
pub fn app_component(ctx: &mut RenderContext) -> VNode {
    let theme = Theme::dark();

    let mesh_storage = ctx
        .use_store::<AppState>()
        .map(|s| s.with(|state| state.mesh_storage.clone()))
        .unwrap_or_default();

    let viewport_2d = build_viewport_2d(ctx, &theme, &mesh_storage);
    let viewport_3d = build_viewport_3d(ctx, &theme, &mesh_storage);

    rsx!(ctx,
        <div
            style={VisualStyle::new().background_color(theme.colors.background)}
            layout={
                LayoutStyle::new()
                    .width(Dimension::Percent(100.0))
                    .height(Dimension::Percent(100.0))
                    .flex_direction(FlexDirection::Column)
            }
        >
            // Viewports row
            <div
                layout={
                    LayoutStyle::new()
                        .width(Dimension::Percent(100.0))
                        .height(Dimension::Percent(100.0))
                        .flex_grow(1.0)
                        .flex_direction(FlexDirection::Row)
                        .gap(theme.spacing.xs)
                }
            >
                // 2D viewport with label
                <div
                    layout={
                        LayoutStyle::new()
                            .height(Dimension::Percent(100.0))
                            .flex_grow(1.0)
                            .flex_direction(FlexDirection::Column)
                            .padding(Edges::all(theme.spacing.xs.into()))
                    }
                >
                    <text
                        style={
                            VisualStyle::new()
                                .font_color(theme.colors.primary_hover)
                                .font_size(theme.font_size.sm)
                        }
                        layout={LayoutStyle::new().height(20.0).width(80.0).flex_shrink(0.0)}
                    >
                        "2D View"
                    </text>
                    {viewport_2d}
                </div>

                // 3D viewport with label
                <div
                    layout={
                        LayoutStyle::new()
                            .height(Dimension::Percent(100.0))
                            .flex_grow(1.0)
                            .flex_direction(FlexDirection::Column)
                            .padding(Edges::all(theme.spacing.xs.into()))
                    }
                >
                    <text
                        style={
                            VisualStyle::new()
                                .font_color(theme.colors.warning)
                                .font_size(theme.font_size.sm)
                        }
                        layout={LayoutStyle::new().height(20.0).width(80.0).flex_shrink(0.0)}
                    >
                        "3D View"
                    </text>
                    {viewport_3d}
                </div>
            </div>

            // Status bar
            <div
                style={VisualStyle::new().background_color(theme.colors.surface)}
                layout={
                    LayoutStyle::new()
                        .width(Dimension::Percent(100.0))
                        .height(30.0)
                        .flex_shrink(0.0)
                        .align_items(AlignItems::Center)
                        .justify_content(JustifyContent::Center)
                }
            >
                <text
                    style={
                        VisualStyle::new()
                            .font_color(theme.colors.text_secondary)
                            .font_size(theme.font_size.xs)
                    }
                    layout={LayoutStyle::new().width(600.0).height(18.0)}
                >
                    "Geolis Viewer | 2D: Space+drag pan, Cmd+scroll zoom | 3D: Right-drag orbit, Middle-drag pan"
                </text>
            </div>
        </div>
    )
}

fn build_viewport_2d(ctx: &mut RenderContext, theme: &Theme, mesh_storage: &MeshStorage) -> VNode {
    rsx!(ctx,
        <viewport
            viewer_type={ViewerType::Viewer2D}
            style={VisualStyle::new().border(1.0, theme.colors.primary)}
            layout={
                LayoutStyle::new()
                    .width(Dimension::Percent(100.0))
                    .height(Dimension::Percent(100.0))
                    .flex_grow(1.0)
            }
            mesh_storage={mesh_storage.clone()}
        />
    )
}

fn build_viewport_3d(ctx: &mut RenderContext, theme: &Theme, mesh_storage: &MeshStorage) -> VNode {
    rsx!(ctx,
        <viewport
            viewer_type={ViewerType::Viewer3D}
            style={VisualStyle::new().border(1.0, theme.colors.warning)}
            layout={
                LayoutStyle::new()
                    .width(Dimension::Percent(100.0))
                    .height(Dimension::Percent(100.0))
                    .flex_grow(1.0)
            }
            mesh_storage={mesh_storage.clone()}
        />
    )
}
