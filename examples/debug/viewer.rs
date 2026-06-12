//! Pure viewer UI — state definition and component tree.
//!
//! This module is independent of what meshes are loaded.
//! Feed any `MeshStorage` via `AppState` and the viewer renders it.

#![allow(clippy::too_many_lines)]

use revion_design_system::Theme;
use revion_ui::cx_builders::{div, text, viewport};
use revion_ui::value_objects::{Dimension, Edges};
use revion_ui::{
    style::{AlignItems, FlexDirection, JustifyContent, LayoutStyle, VisualStyle},
    MeshStorage, RenderContext, View, ViewerType,
};

/// Root UI component — dual 2D / 3D viewports with a status bar.
///
/// The mesh storage is captured at app startup and rendered as-is; the
/// viewer holds no reactive state of its own.
pub fn app_component(ctx: &mut RenderContext, mesh_storage: &MeshStorage) -> View {
    let theme = Theme::dark();

    let viewport_2d = build_viewport_2d(ctx, &theme, mesh_storage);
    let viewport_3d = build_viewport_3d(ctx, &theme, mesh_storage);

    div()
        .style(VisualStyle::new().background_color(theme.colors.background))
        .layout(
            LayoutStyle::new()
                .width(Dimension::Percent(100.0))
                .height(Dimension::Percent(100.0))
                .flex_direction(FlexDirection::Column),
        )
        .child(
            // Viewports row
            div()
                .layout(
                    LayoutStyle::new()
                        .width(Dimension::Percent(100.0))
                        .height(Dimension::Percent(100.0))
                        .flex_grow(1.0)
                        .flex_direction(FlexDirection::Row)
                        .gap(theme.spacing.xs),
                )
                .child(
                    // 2D viewport with label
                    div()
                        .layout(
                            LayoutStyle::new()
                                .height(Dimension::Percent(100.0))
                                .flex_grow(1.0)
                                .flex_direction(FlexDirection::Column)
                                .padding(Edges::all(theme.spacing.xs.into())),
                        )
                        .child(
                            text("2D View")
                                .style(
                                    VisualStyle::new()
                                        .font_color(theme.colors.primary_hover)
                                        .font_size(theme.font_size.sm),
                                )
                                .layout(
                                    LayoutStyle::new()
                                        .height(20.0)
                                        .width(80.0)
                                        .flex_shrink(0.0),
                                ),
                        )
                        .child(viewport_2d),
                )
                .child(
                    // 3D viewport with label
                    div()
                        .layout(
                            LayoutStyle::new()
                                .height(Dimension::Percent(100.0))
                                .flex_grow(1.0)
                                .flex_direction(FlexDirection::Column)
                                .padding(Edges::all(theme.spacing.xs.into())),
                        )
                        .child(
                            text("3D View")
                                .style(
                                    VisualStyle::new()
                                        .font_color(theme.colors.warning)
                                        .font_size(theme.font_size.sm),
                                )
                                .layout(
                                    LayoutStyle::new()
                                        .height(20.0)
                                        .width(80.0)
                                        .flex_shrink(0.0),
                                ),
                        )
                        .child(viewport_3d),
                ),
        )
        .child(
            // Status bar
            div()
                .style(VisualStyle::new().background_color(theme.colors.surface))
                .layout(
                    LayoutStyle::new()
                        .width(Dimension::Percent(100.0))
                        .height(30.0)
                        .flex_shrink(0.0)
                        .align_items(AlignItems::Center)
                        .justify_content(JustifyContent::Center),
                )
                .child(
                    text(
                        "Geolis Viewer | 2D: Space+drag pan, Cmd+scroll zoom | 3D: Right-drag orbit, Middle-drag pan",
                    )
                    .style(
                        VisualStyle::new()
                            .font_color(theme.colors.text_secondary)
                            .font_size(theme.font_size.xs),
                    )
                    .layout(LayoutStyle::new().width(600.0).height(18.0)),
                ),
        )
        .build_cx(ctx)
}

fn build_viewport_2d(ctx: &mut RenderContext, theme: &Theme, mesh_storage: &MeshStorage) -> View {
    viewport(ViewerType::Viewer2D)
        .style(VisualStyle::new().border(1.0, theme.colors.primary))
        .layout(
            LayoutStyle::new()
                .width(Dimension::Percent(100.0))
                .height(Dimension::Percent(100.0))
                .flex_grow(1.0),
        )
        .mesh_storage(mesh_storage.clone())
        .build_cx(ctx)
}

fn build_viewport_3d(ctx: &mut RenderContext, theme: &Theme, mesh_storage: &MeshStorage) -> View {
    viewport(ViewerType::Viewer3D)
        .style(VisualStyle::new().border(1.0, theme.colors.warning))
        .layout(
            LayoutStyle::new()
                .width(Dimension::Percent(100.0))
                .height(Dimension::Percent(100.0))
                .flex_grow(1.0),
        )
        .mesh_storage(mesh_storage.clone())
        .build_cx(ctx)
}
