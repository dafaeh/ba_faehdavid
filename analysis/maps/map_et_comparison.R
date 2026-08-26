# map_et_comparison.R
# 2 x 2 panel figure: DisALEXI and OpenET Ensemble mean monthly ET 2020-2022 in
# the top row, DisALEXI-based CWD_max 2020-2022 bottom left, bottom right empty.
# The ET window is a subset of the Northern California focus region, so CWD_max
# is cropped out of the focus region export.


# ---- Setup ------------------------------------------------------------------
library(terra)
library(tidyterra)
library(ggplot2)
library(ggspatial)
library(cowplot)
library(sf)
library(here)

# ---- Load data --------------------------------------------------------------
r_et         <- terra::rast(here("data-raw", "et_comparison_2020_2022.tif"))
r_cwd_norcal <- terra::rast(here("data", "cwd_northern_california.tif"))[["cwd_max"]]

# ---- Crop CWD_max to the ET window ------------------------------------------
# Both rasters come out of the same export settings (EPSG:5070, 30 m), so a
# plain crop is enough. The fallback reprojects onto the ET grid, bilinear
# because CWD_max is continuous.
if (terra::same.crs(r_cwd_norcal, r_et)) {
  r_cwd <- terra::crop(r_cwd_norcal, terra::ext(r_et))
} else {
  poly_et <- terra::project(
    terra::as.polygons(terra::ext(r_et), crs = terra::crs(r_et)),
    terra::crs(r_cwd_norcal)
  )
  r_cwd <- terra::project(
    terra::crop(r_cwd_norcal, poly_et),
    r_et,
    method = "bilinear"
  )
}

# ---- Colour scale limits ----------------------------------------------------
et_limits  <- c(70, 115)     # shared by both ET panels
cwd_limits <- c(610, 1080)   # 2nd/98th percentile of the cropped raster

# ---- Map annotations --------------------------------------------------------
# ggspatial sizes these in absolute units, so they do not shrink with the panel.
# which_north = "grid" keeps the arrow straight up; true north looks skewed
# under Albers away from the central meridian.
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

# ---- Map function -----------------------------------------------------------
make_map <- function(r, band_name, limits, legend_name) {
  
  ggplot() +
    tidyterra::geom_spatraster(data = r, aes(fill = .data[[band_name]])) +
    scale_fill_viridis_c(
      name      = legend_name,
      limits    = limits,
      direction = -1,            # reversed: dark = high values
      oob       = scales::squish,
      na.value  = NA
    ) +
    # Default axis expansion, so the scale bar sits in the margin below the
    # raster rather than on top of it.
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
    # Native CRS of the raster, no reprojection. datum = 4326 keeps the
    # graticule and axis labels in degrees.
    coord_sf(
      crs   = terra::crs(r),
      datum = sf::st_crs(4326)
    ) +
    
    layer_scalebar +
    layer_north_arrow +
    
    labs(x = "Longitude", y = "Latitude") +
    theme_classic(base_size = 13) +
    theme(
      axis.ticks        = element_line(colour = "grey30", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text         = element_text(size = 8),
      axis.title        = element_text(size = 11)
    )
}

# ---- Build panels -----------------------------------------------------------
map_disalexi <- make_map(r_et,  "et_disalexi", et_limits,
                         expression(ET~"[mm month"^-1*"]"))
map_ensemble <- make_map(r_et,  "et_ensemble", et_limits,
                         expression(ET~"[mm month"^-1*"]"))
map_cwd      <- make_map(r_cwd, "cwd_max",     cwd_limits,
                         expression(CWD[max]~"[mm]"))

# ---- Assemble figure --------------------------------------------------------
# align/axis keep all panels the same size despite the wider CWD_max legend.
# NULL leaves the bottom-right cell empty, the empty label suppresses its tag.
fig_et_comparison <- plot_grid(
  map_disalexi, map_ensemble,
  map_cwd,      NULL,
  ncol   = 2,
  labels = c("a", "b", "c", ""),
  align  = "hv",
  axis   = "tblr"
)

fig_et_comparison

# ---- Save -------------------------------------------------------------------
# PNG rather than PDF: the dense 30 m rasters bloat vector output.
ggsave(
  filename = here("fig", "map_et_comparison.png"),
  plot     = fig_et_comparison,
  device   = ragg::agg_png,
  width    = 24, height = 20, units = "cm",
  dpi      = 600
)
