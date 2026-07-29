# 06_plot_zoom_maps.R
# CWDmax and gap-filled months for the three zoom windows, combined into a
# single 3 x 2 panel figure with shared legends.

# ---- Setup ------------------------------------------------------------------
library(ggplot2)
library(tidyterra)
library(ggspatial)
library(patchwork)
library(here)

# ---- Zoom windows -----------------------------------------------------------
# Converted from WGS84 rectangles drawn in the GEE code editor.
# Table order is figure row order.
df_windows <- tibble::tribble(
  ~id,   ~label,    ~xmin,      ~ymin,      ~xmax,      ~ymax,
  "z05", "Zoom A",  -498940.6,  1597756.5,  -376923.9,  1708584.7,
  "z06", "Zoom B",  1134811.7,  1653744.2,  1265610.5,  1774527.1,
  "z04", "Zoom C", -2138564.8,  2060683.0, -2002554.8,  2188453.0
)

# ---- Load data --------------------------------------------------------------
r_conus <- terra::rast(here("data", "cwd_conus.tif"))

# ---- Colour scales ----------------------------------------------------------
# Fixed limits, so every panel shares the same colour scale. Defined as objects
# because patchwork's guides = "collect" only merges guides that match exactly.
scale_cwd <- scale_fill_viridis_c(
  name      = expression(CWD[max] ~ "[mm]"),
  direction = -1,
  limits    = c(0, 840),
  oob       = scales::squish,
  na.value  = NA
)

scale_gap <- scico::scale_fill_scico(
  name      = "Gap-filled months [count]",
  palette   = "lajolla",
  direction = -1,
  limits    = c(0, 36),
  na.value  = NA
)

# ---- Map annotations --------------------------------------------------------
# ggspatial sizes these in absolute units, so they do not shrink with the panel.
# Scaled down from the single-map version, where the defaults overflowed the
# axis expansion and covered the raster.
# which_north = "grid" keeps the arrow pointing straight up (grid north) instead
# of true north, which can look skewed under the Albers projection away from its
# central meridian.
layer_scalebar <- annotation_scale(
  location   = "bl",
  width_hint = 0.32,
  style      = "bar",
  bar_cols   = c("black", "white"),
  height     = grid::unit(0.1, "cm"),
  pad_x      = grid::unit(0.15, "cm"),
  pad_y      = grid::unit(0.12, "cm"),
  text_pad   = grid::unit(0.08, "cm"),
  text_cex   = 0.55,
  line_width = 0.4
)

layer_north_arrow <- annotation_north_arrow(
  location    = "br",
  which_north = "grid",
  height      = grid::unit(0.5, "cm"),
  width       = grid::unit(0.5, "cm"),
  pad_x       = grid::unit(0.15, "cm"),
  pad_y       = grid::unit(0.15, "cm"),
  style       = north_arrow_fancy_orienteering(
    fill       = c("black", "grey40"),
    line_col   = "grey20",
    line_width = 0.5,
    text_size  = 5
  )
)

# ---- Local helper -----------------------------------------------------------
# Used only in this script (called once per panel).
# Scale bar in the left column only: the windows differ in ground extent by up
# to 11 %, so every row carries its own map scale, but both panels of a row
# share one. North arrow once in the bottom-right panel, since all panels are
# EPSG:5070 and grid north is identical everywhere.
build_zoom_panel <- function(r_crop, fill_var, scale_fill,
                             show_scalebar = FALSE, show_arrow = FALSE,
                             title = NULL) {
  p <- ggplot() +
    geom_spatraster(
      data    = r_crop,
      aes(fill = .data[[fill_var]]),
      maxcell = 1e6
    ) +
    scale_fill +
    # Breaks under coord_sf are degrees whatever the plotting CRS. Fixed 0.5 deg
    # steps replace the automatic breaks, which collide in the narrower windows.
    # The wider expansion at the bottom leaves room for the scale bar.
    scale_x_continuous(breaks = scales::breaks_width(0.5)) +
    scale_y_continuous(
      breaks = scales::breaks_width(0.5),
      expand = expansion(mult = c(0.09, 0.05))
    ) +
    coord_sf(crs = terra::crs(r_crop)) +
    labs(x = "Longitude", y = "Latitude", title = title) +
    # base_size 9 rather than the default 11: panels are roughly a third of the
    # width they had as single maps.
    theme_classic(base_size = 9) +
    theme(
      plot.title        = element_text(size = 9, face = "plain", hjust = 0),
      plot.margin       = margin(2, 2, 2, 2),
      plot.tag          = element_text(size = 9, face = "bold"),
      legend.title      = element_text(size = 8),
      legend.text       = element_text(size = 7),
      legend.key.width  = grid::unit(1.4, "cm"),
      legend.key.height = grid::unit(0.3, "cm")
    )
  
  if (show_scalebar) {
    p <- p + layer_scalebar
  }
  
  if (show_arrow) {
    p <- p + layer_north_arrow
  }
  
  p
}

# ---- Build panels -----------------------------------------------------------
# Every panel keeps a fully labelled x axis. The windows sit at different
# eastings, so a single shared axis in the bottom row would imply an alignment
# that does not exist.
n_rows <- nrow(df_windows)

list_panels <- df_windows |>
  dplyr::mutate(row_index = dplyr::row_number()) |>
  purrr::pmap(function(id, label, xmin, ymin, xmax, ymax, row_index) {
    r_crop <- terra::crop(
      r_conus,
      terra::ext(xmin, xmax, ymin, ymax)
    )
    
    list(
      cwd = build_zoom_panel(
        r_crop, "cwd_max", scale_cwd,
        show_scalebar = TRUE,
        title         = label
      ),
      gap = build_zoom_panel(
        r_crop, "gap_filled_months", scale_gap,
        show_arrow = row_index == n_rows
      )
    )
  })

# ---- Assemble figure --------------------------------------------------------
# list_flatten yields cwd, gap, cwd, gap, ..., filling the two columns row-wise.
# Legend position is set on the patchwork, not on the panels: it applies to the
# collected guides, which belong to the figure rather than to any single panel.
p_zooms <- list_panels |>
  purrr::list_flatten() |>
  patchwork::wrap_plots(ncol = 2, byrow = TRUE) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    tag_levels = "a",
    tag_suffix = ")",
    theme      = theme(legend.position = "bottom")
  )

# Save 
ggsave(
  filename = here("fig", "cwd_gapfill_zooms.png"),
  plot     = p_zooms,
  device   = ragg::agg_png,
  width    = 16, height = 24, units = "cm",
  dpi      = 600
)