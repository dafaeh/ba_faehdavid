# 06_plot_zoom_maps.R
# CWDmax overview map with zoom-ins, each saved as a separate map from the
# gap_filled_months raster for the same extent.

# ---- Setup ------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyterra)
library(ggspatial)
library(scico)
library(scales)
library(here)
library(ragg)
# Zoom windows (EPSG:5070) 
# Converted from WGS84 rectangles drawn in the GEE code editor.
df_windows <- tibble::tribble(
  ~id,   ~xmin,      ~ymin,      ~xmax,      ~ymax,
  "z04", -2138564.8, 2060683.0, -2002554.8, 2188453.0,
  "z05", -498940.6,  1597756.5, -376923.9,  1708584.7,
  "z06", 1134811.7,  1653744.2, 1265610.5,  1774527.1
)
# Load data 
r_conus <- terra::rast(here("data", "cwd_conus.tif"))

# CWDmax colour scale 
# Fixed limits, so the overview map and every zoom panel share the same colour scale
cwd_limits    <- c(0, 840)
cwd_direction <- -1

# Gap-fill colour scale 
gap_limits <- c(0, 36)

# ---- Local helper -------------------------------------------------------
# Used only in this script (called once per window). Returns the two panels as a list,
# so each can be saved as its own map. Scale bar and north arrow appear
# on both panels, since each is now saved separately. which_north = "grid"
# keeps the arrow pointing straight up (grid north) instead of true north,
# which can look skewed under the Albers projection away from its central
# meridian.
plot_zoom_pair <- function(r, window, cwd_limits, gap_limits) {
  ext_win <- terra::ext(window$xmin, window$xmax, window$ymin, window$ymax)
  r_crop  <- terra::crop(r, ext_win)
  
  p_cwd <- ggplot() +
    geom_spatraster(data = r_crop, aes(fill = cwd_max)) +
    scale_fill_viridis_c(
      name      = expression(CWD[max]~"[mm]"),
      direction = cwd_direction,
      limits    = cwd_limits,
      oob       = scales::squish,
      na.value  = NA
    ) +
    annotation_scale(
      location   = "bl",
      width_hint = 0.15,
      style      = "bar",
      bar_cols   = c("black", "white"),
      height     = grid::unit(0.12, "cm"),
      text_cex   = 0.8
    ) +
    annotation_north_arrow(
      location    = "br",
      which_north = "grid",
      style       = north_arrow_fancy_orienteering(
        fill     = c("black", "grey40"),
        line_col = "grey20"
      ),
      height = grid::unit(0.8, "cm"),
      width  = grid::unit(0.8, "cm")
    ) +
    coord_sf(crs = terra::crs(r_crop)) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic()
  
  p_gap <- ggplot() +
    geom_spatraster(data = r_crop, aes(fill = gap_filled_months)) +
    scale_fill_scico(
      name      = "Gap-filled\nmonths",
      palette   = "lajolla",
      direction = -1,
      limits    = gap_limits,
      na.value  = NA
    ) +
    annotation_scale(
      location   = "bl",
      width_hint = 0.15,
      style      = "bar",
      bar_cols   = c("black", "white"),
      height     = grid::unit(0.12, "cm"),
      text_cex   = 0.8
    ) +
    annotation_north_arrow(
      location    = "br",
      which_north = "grid",
      style       = north_arrow_fancy_orienteering(
        fill     = c("black", "grey40"),
        line_col = "grey20"
      ),
      height = grid::unit(0.8, "cm"),
      width  = grid::unit(0.8, "cm")
    ) +
    coord_sf(crs = terra::crs(r_crop)) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic()
  
  list(cwd = p_cwd, gap = p_gap)
}

# Build plots 
list_pairs <- df_windows |>
  purrr::pmap(\(...) {
    window <- tibble::tibble(...)
    plot_zoom_pair(r_conus, window, cwd_limits, gap_limits)
  })
names(list_pairs) <- df_windows$id

# Save
# Two separate files per window
purrr::iwalk(list_pairs, \(pair, id) {
  ggsave(
    filename = here("fig", paste0("cwd_zoom_", id, ".png")),
    plot     = pair$cwd,
    device   = ragg::agg_png,
    width    = 22, height = 15, units = "cm",
    dpi      = 600
  )
  ggsave(
    filename = here("fig", paste0("gapfilled_zoom_", id, ".png")),
    plot     = pair$gap,
    device   = ragg::agg_png,
    width    = 22, height = 15, units = "cm",
    dpi      = 600
  )
})